"""
    get_rts_raw_sys(rts_src_dir, rts_siip_dir;)

Get the raw RTS system.
"""
function get_rts_raw_sys(rts_src_dir, rts_siip_dir;)
    return rawsys = PSY.PowerSystemTableData(
        rts_src_dir,
        100.0,
        joinpath(rts_siip_dir, "user_descriptors.yaml");
        timeseries_metadata_file=joinpath(rts_siip_dir, "timeseries_pointers.json"),
        generator_mapping_file=joinpath(rts_siip_dir, "generator_mapping.yaml"),
    )
end

"""
    get_rts_sys(rts_src_dir, 
        rts_siip_dir; 
        time_series_resolution = [Dates.Hour(1), Dates.Minute(5)],
        modifier_function::Function = rts_modifier_function!,
        kwargs...
    )

Prepare the RTS system for the Day Ahead (DA) and Real Time (RT) problems
"""
function get_rts_sys(
    rts_src_dir,
    rts_siip_dir;
    time_series_resolution=[Dates.Hour(1), Dates.Minute(5)],
    modifier_function::Function=rts_modifier_function!,
    kwargs...,
)
    rawsys = get_rts_raw_sys(rts_src_dir, rts_siip_dir)
    sys_DA = System(rawsys; time_series_resolution=time_series_resolution[1], kwargs...)
    sys_rt = System(rawsys; time_series_resolution=time_series_resolution[2], kwargs...)
    modifier_function(sys_DA, sys_rt;)
    return sys_DA, sys_rt
end

"""
    rts_modifier_function!(
        sys_DA=sys_DA, sys_rt=sys_rt; mult=1.5, DISPATCH_INCREASE=2.0, FIX_DECREASE=0.3
    )

Modify the RTS
"""
function rts_modifier_function!(
    sys_DA=sys_DA, sys_rt=sys_rt; mult=1.5, DISPATCH_INCREASE=2.0, FIX_DECREASE=0.3
)
    # Add area renewable energy forecasts for RT model
    area_mapping = get_aggregation_topology_mapping(Area, sys_rt)
    for (k, buses_in_area) in area_mapping
        k == "1" && continue
        remove_component!(sys_rt, get_component(Area, sys_rt, k))
        for b in buses_in_area
            PSY.set_area!(b, get_component(Area, sys_rt, "1"))
        end
    end

    # Make data more realistic
    for sys in [sys_DA, sys_rt]
        # Adjust Reserve Provisions
        # Remove Flex Reserves
        set_units_base_system!(sys, "SYSTEM_BASE")
        res_up = get_component(VariableReserve{ReserveUp}, sys, "Flex_Up")
        remove_component!(sys, res_up)
        res_dn = get_component(VariableReserve{ReserveDown}, sys, "Flex_Down")
        remove_component!(sys, res_dn)
        reg_reserve_up = get_component(VariableReserve, sys, "Reg_Up")
        mult = (sys == sys_DA) ? 1.5 : 1.0
        set_requirement!(reg_reserve_up, mult * get_requirement(reg_reserve_up))
        reg_reserve_dn = get_component(VariableReserve, sys, "Reg_Down")
        mult = (sys == sys_DA) ? 1.5 : 1.0
        set_requirement!(reg_reserve_dn, mult * get_requirement(reg_reserve_dn))
        spin_reserve_R1 = get_component(VariableReserve, sys, "Spin_Up_R1")
        spin_reserve_R2 = get_component(VariableReserve, sys, "Spin_Up_R2")
        spin_reserve_R3 = get_component(VariableReserve, sys, "Spin_Up_R3")

        for g in get_components(
            ThermalStandard,
            sys,
            x -> get_prime_mover(x) in [PrimeMovers.CT, PrimeMovers.CC],
        )
            if get_fuel(g) == ThermalFuels.DISTILLATE_FUEL_OIL
                remove_component!(sys, g)
                continue
            end
            g.operation_cost.shut_down = g.operation_cost.start_up / 2.0

            if PSY.get_base_power(g) > 3
                #remove_service!(g, reg_reserve_dn)
                #remove_service!(g, reg_reserve_up)
                continue
            end
            clear_services!(g)
            add_service!(g, reg_reserve_dn)
            add_service!(g, reg_reserve_up)

            if get_prime_mover(g) == PrimeMovers.CT
                set_status!(g, false)
                set_active_power!(g, 0.0)
                old_pwl_array = get_cost(get_variable(get_operation_cost(g)))
                new_pwl_array = similar(old_pwl_array)
                for (ix, tup) in enumerate(old_pwl_array)
                    if ix ∈ [1, length(old_pwl_array)]
                        cost_noise = 50.0 * rand()
                        new_pwl_array[ix] = ((tup[1] + cost_noise), tup[2])
                    else
                        try_again = true
                        while try_again
                            cost_noise = 50.0 * rand()
                            power_noise = 0.01 * rand()
                            slope_previous =
                                ((tup[1] + cost_noise) - old_pwl_array[ix - 1][1]) /
                                ((tup[2] - power_noise) - old_pwl_array[ix - 1][2])
                            slope_next =
                                (-(tup[1] + cost_noise) + old_pwl_array[ix + 1][1]) /
                                (-(tup[2] - power_noise) + old_pwl_array[ix + 1][2])
                            new_pwl_array[ix] = (
                                (tup[1] + cost_noise), (tup[2] - power_noise)
                            )
                            try_again = slope_previous > slope_next
                        end
                    end
                end
                get_variable(get_operation_cost(g)).cost = new_pwl_array
            end
        end
        #Remove units that make no sense to include
        names = [
            "114_SYNC_COND_1",
            "314_SYNC_COND_1",
            "313_STORAGE_1",
            "214_SYNC_COND_1",
            "212_CSP_1",
        ]
        for d in get_components(Generator, sys, x -> x.name ∈ names)
            remove_component!(sys, d)
        end
        for br in get_components(DCBranch, sys)
            remove_component!(sys, br)
        end
        for d in get_components(Storage, sys)
            remove_component!(sys, d)
        end

        # Remove large Coal and Nuclear from reserves
        for d in get_components(
            ThermalStandard, sys, x -> (occursin(r"STEAM|NUCLEAR", get_name(x)))
        )
            get_fuel(d) == ThermalFuels.COAL &&
                (get_ramp_limits(d) = (up=0.001, down=0.001))
            if get_fuel(d) == ThermalFuels.DISTILLATE_FUEL_OIL
                remove_component!(sys, d)
                continue
            end
            get_operation_cost(d).shut_down = get_operation_cost(d).start_up / 2.0
            if get_rating(d) < 3
                set_status!(d, false)
                #clear_services!(d)
                #reserve_hydro && add_service!(d, reg_reserve_up)
                #reserve_hydro && add_service!(d, reg_reserve_dn)
                #add_service!(d, spin_reserve_R1)
                set_status!(d, false)
                set_active_power!(d, 0.0)
                continue
            end
            clear_services!(d)
            get_operation_cost(d).shut_down = get_operation_cost(d).start_up / 2.0
            if get_fuel(d) == ThermalFuels.NUCLEAR
                get_ramp_limits(d) = (up=0.0, down=0.0)
                get_time_limits(d) = (up=4380.0, down=4380.0)
            end
        end

        for d in get_components(RenewableDispatch, sys)
            clear_services!(d)
        end

        # Add Hydro to regulation reserves
        for d in get_components(HydroEnergyReservoir, sys)
            remove_component!(sys, d)
        end

        for d in get_components(HydroDispatch, sys)
            clear_services!(d)
        end

        for g in get_components(
            RenewableDispatch, sys, x -> get_prime_mover(x) == PrimeMovers.PVe
        )
            rat_ = get_rating(g)
            set_rating!(g, DISPATCH_INCREASE * rat_)
        end

        for g in
            get_components(RenewableFix, sys, x -> get_prime_mover(x) == PrimeMovers.PVe)
            rat_ = get_rating(g)
            set_rating!(g, FIX_DECREASE * rat_)
        end
    end
end
