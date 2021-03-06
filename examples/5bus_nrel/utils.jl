using Core: GeneratedFunctionStub
"""
    build_5_bus_matpower_DA(
        data_dir::AbstractString;
        case_file::AbstractString="case5_re_uc.m",
        FORECASTS_DIR::AbstractString=joinpath(data_dir, "forecasts"),
        forecasts_pointers_file::AbstractString=joinpath(
            FORECASTS_DIR, "timeseries_pointers_da_7day.json"
        ),
        add_reserves::Bool=true,
    )

Builds base system for the 5bus NREL case (a.k.a NESTA case) from:
 - A matpower file containing grid information (case_file);
 - A file describing forecasts locations and details (forecasts_pointers_file);
 - A simple definition of reserves.
"""
function build_5_bus_matpower_DA(
    data_dir::AbstractString;
    case_file::AbstractString="case5_re_uc.m",
    FORECASTS_DIR::AbstractString=joinpath(data_dir, "forecasts"),
    forecasts_pointers_file::AbstractString=joinpath(
        FORECASTS_DIR, "timeseries_pointers_da_7day.json"
    ),
    add_reserves::Bool=true,
)
    case_file_path = joinpath(data_dir, case_file)
    sys = System(case_file_path)
    if add_reserves
        reserves = [
            VariableReserve{ReserveUp}("REG1", true, 5.0, 0.2),
            VariableReserve{ReserveUp}("REG2", true, 5.0, 0.4),
            VariableReserve{ReserveUp}("REG3", true, 5.0, 0.3),
            VariableReserve{ReserveUp}("REG4", true, 5.0, 0.4),
            VariableReserve{ReserveDown}("REG1d", true, 5.0, 0.2),
            VariableReserve{ReserveDown}("REG2d", true, 5.0, 0.4),
            VariableReserve{ReserveDown}("REG3d", true, 5.0, 0.3),
            VariableReserve{ReserveDown}("REG4d", true, 5.0, 0.4),
        ]
        contributing_devices = get_components(Generator, sys)
        for r in reserves
            add_service!(sys, r, contributing_devices)
        end
    end

    add_time_series!(sys, forecasts_pointers_file)

    return sys
end

"""
    build_5_bus_matpower_RT(
        data_dir::AbstractString;
        case_file::AbstractString="case5_re_uc.m",
        FORECASTS_DIR::AbstractString=joinpath(data_dir, "forecasts"),
        forecasts_pointers_file::AbstractString=joinpath(
            FORECASTS_DIR, "timeseries_pointers_rt_7day.json"
        ),
    )

Builds base system for the 5bus NREL case (a.k.a NESTA case) from:
 - A matpower file containing grid information (case_file);
 - A file describing forecasts locations and details (forecasts_pointers_file);
"""
function build_5_bus_matpower_RT(
    data_dir::AbstractString;
    case_file::AbstractString="case5_re_uc.m",
    FORECASTS_DIR::AbstractString=joinpath(data_dir, "forecasts"),
    forecasts_pointers_file::AbstractString=joinpath(
        FORECASTS_DIR, "timeseries_pointers_rt_7day.json"
    ),
)
    case_file_path = joinpath(data_dir, case_file)
    sys = System(case_file_path)

    add_time_series!(sys, forecasts_pointers_file)
    transform_single_time_series!(sys, 12, Minute(5))

    return sys
end

"""
    prep_systems_UCED(
        system::System;
        horizon_uc::Int=24,
        horizon_ed::Int=1,
        interval_uc::TimePeriod=Hour(24),
        interval_ed::TimePeriod=Hour(1),
    )

Duplicates the system to represent UC and ED for DA, transforming the time series
to the appropriate interval and horizon.
"""
function prep_systems_UCED(
    system::System;
    horizon_uc::Int=24,
    horizon_ed::Int=1,
    interval_uc::TimePeriod=Hour(24),
    interval_ed::TimePeriod=Hour(1),
)
    system_uc = system
    system_ed = deepcopy(system)

    transform_single_time_series!(system_uc, horizon_uc, interval_uc)
    transform_single_time_series!(system_ed, horizon_ed, interval_ed)

    return system_uc, system_ed
end

"""
    run_set_of_simulations(
        df::DataFrame, 
        data_dir::String, 
        example_dir::String, 
        solver_uc, 
        solver_ed, 
        solver_rt, 
        range_quota::Vector{Float64}, 
        initial_time::Date, 
        initial_bidding_time::DateTime,
        path::String
    )

Run a set of simulations to 5bus_nrel with descripted sistems as in `df`
"""
function run_set_of_simulations(df::DataFrame, 
    data_dir::String, 
    example_dir::String, 
    solver_uc, 
    solver_ed, 
    solver_rt, 
    range_quota::Vector{Float64}, 
    initial_time::Date, 
    initial_bidding_time::DateTime,
    path::String
)
    for l=1:size(df)[1]

        directory_name= "Net_"*string(df.Network[l]["DA"])[1:4]*"_"*string(df.Network[l]["RT"])[1:4]* "__Ramp_"*string(df.Ramp[l]["DA"])*
        "_"*string(df.Ramp[l]["RT"])*"__Min_gen_"*string(df.Minimal_generation[l]["DA"])*"_"*string(df.Minimal_generation[l]["RT"])* "__Reserve_" * string(df.Reserve[l]) *
        "__Offer_" * string(df.Offer_Bus[l]) * "__Bid_period_1-" * string(length(df.bidding_period[l]))

        if isdir(joinpath(example_dir, path, directory_name))==false

            # call our data preparation to build base system
            # the case was modified to not have hydros nor transformers
            sys_rt = build_5_bus_matpower_RT(data_dir;)

            base_da_system = build_5_bus_matpower_DA(
                data_dir;
                # using a modified (mod) file that reduced load for feasibility in DC-OPF
                forecasts_pointers_file=joinpath(
                    data_dir, "forecasts", "timeseries_pointers_da_7day_mod.json"
                ),
                add_reserves=df.Reserve[l],
            )

            generator_metadata = Dict()
            generator_metadata["RT"] = [gen for gen in get_components(Generator, sys_rt)]
            generator_metadata["DA"] = [gen for gen in get_components(Generator, base_da_system)]

            for i in keys(generator_metadata)
                #change minimal generation if it is false
                if df.Minimal_generation[l][i]==false
                    for generator in generator_metadata[i]
                        try
                            active_power_limits = (min=0.0, max=generator.active_power_limits[:max])              
                            PowerSystems.set_active_power_limits!(generator, active_power_limits)
                        catch
                        end
                    end
                end
                #change ramp if it is false
                if df.Ramp[l][i]==false
                    for generator in generator_metadata[i]
                        try
                            ramp_limits = (up=0.0, down=0.0)                
                            PowerSystems.set_ramp_limits!(generator, ramp_limits) 
                        catch
                        end
                    end
                end
            end

            # Add single generator at a defined bus
            gen = add_generator!(base_da_system, df.Offer_Bus[l], (min=0.0, max=0.0))
            @test gen in get_components(Generator, base_da_system)

            # create and set variable cost time-series for the generator
            ts_array = create_generator_bids(;
                initial_bidding_time=initial_bidding_time,
                bidding_periods=df.bidding_period[l],
                system=base_da_system,
                costs=zeros(length(df.bidding_period[l])),
            )
            set_variable_cost!(base_da_system, gen, ts_array)

            # duplicate system and prepare times series for the time varying parameters (loads, renewables, ...)
            sys_uc, sys_ed = prep_systems_UCED(base_da_system)

            # generic market formulation templates with defined network formulation
            # CopperPlate-OPF: network=CopperPlatePowerModel
            # DC-OPF: network=DCPPowerModel
            # NFA-OPF (only line limit constraints): network=NFAPowerModel
            # DC-PTDF-OPF (what ISOs do): network=StandardPTDFModel
            template_uc = template_unit_commitment(; network=df.Network[l]["DA"])
            template_rt = template_economic_dispatch(; network=df.Network[l]["RT"])
            template_ed = template_economic_dispatch(; network=df.Network[l]["DA"])

            # build a market clearing simulator
            market_simulator = UCEDRT(;
                system_uc=sys_uc,
                system_rt=sys_rt,
                system_ed=sys_ed,
                template_uc=template_uc,
                template_rt=template_rt,
                template_ed=template_ed,
                solver_uc=solver_uc,
                solver_rt=solver_rt,
                solver_ed=solver_ed,
            )

            #Calculates the dispatch result for a bid curve
            name_generator = get_name(gen);
            steps = 1;

            mkdir(joinpath(example_dir, path, directory_name))
            simulation_folder = joinpath(example_dir, path, directory_name)

            lmps_df, results_df = pq_curves_virtuals!(
                market_simulator, name_generator, range_quota, initial_time, steps, simulation_folder
            )

        end
    end
end

"""
    run_set_of_simulations_mix(
        df::DataFrame, 
        data_dir::String, 
        example_dir::String, 
        solver_uc, 
        solver_ed, 
        solver_rt, 
        range_quota_load::Vector{Float64}, 
        range_quota_gen::Vector{Float64}, 
        initial_time::Date, 
        initial_bidding_time::DateTime,
        path::String
    )

Run a set of simulations to 5bus_nrel with descripted sistems as in 'df'
"""
function run_set_of_simulations_mix(df::DataFrame, 
    data_dir::String, 
    example_dir::String, 
    solver_uc, 
    solver_ed, 
    solver_rt, 
    range_quota_load::Vector{Float64}, 
    range_quota_gen::Vector{Float64}, 
    initial_time::Date, 
    initial_bidding_time::DateTime,
    path::String
)
    for l=1:size(df)[1]

        directory_name= df.directory_name[l]

        if isdir(joinpath(example_dir, path, directory_name))==false

            # call our data preparation to build base system
            # the case was modified to not have hydros nor transformers
            sys_rt = build_5_bus_matpower_RT(data_dir;)

            base_da_system = build_5_bus_matpower_DA(
                data_dir;
                # using a modified (mod) file that reduced load for feasibility in DC-OPF
                forecasts_pointers_file=joinpath(
                    data_dir, "forecasts", "timeseries_pointers_da_7day_mod.json"
                ),
                add_reserves=df.Reserve[l],
            )

            generator_metadata = Dict()
            generator_metadata["RT"] = [gen for gen in get_components(Generator, sys_rt)]
            generator_metadata["DA"] = [gen for gen in get_components(Generator, base_da_system)]

            for i in keys(generator_metadata)
                #change minimal generation if it is false
                if df.Minimal_generation[l][i]==false
                    for generator in generator_metadata[i]
                        try
                            active_power_limits = (min=0.0, max=generator.active_power_limits[:max])              
                            PowerSystems.set_active_power_limits!(generator, active_power_limits)
                        catch
                        end
                    end
                end
                #change ramp if it is false
                if df.Ramp[l][i]==false
                    for generator in generator_metadata[i]
                        try
                            ramp_limits = (up=0.0, down=0.0)                
                            PowerSystems.set_ramp_limits!(generator, ramp_limits) 
                        catch
                        end
                    end
                end
            end

            # create and set variable time-series for the load
            ts_array = create_demand_series(;
                initial_bidding_time=initial_bidding_time,
                bidding_periods=df.bidding_period[l],
                system=base_da_system,
                demands=ones(length(df.bidding_period[l])),
            )

            # Add single load at a defined bus
            node_load = df.Load_Bus[l] # define bus
            load = add_load!(base_da_system, node_load, 1.0)

            add_time_series!(base_da_system, load, ts_array)

            # Add single generator at a defined bus
            gen = add_generator!(base_da_system, df.Offer_Bus[l], (min=0.0, max=0.0))

            # create and set variable cost time-series for the generator
            ts_array = create_generator_bids(;
                initial_bidding_time=initial_bidding_time,
                bidding_periods=df.bidding_period[l],
                system=base_da_system,
                costs=zeros(length(df.bidding_period[l])),
            )
            set_variable_cost!(base_da_system, gen, ts_array)

            # duplicate system and prepare times series for the time varying parameters (loads, renewables, ...)
            sys_uc, sys_ed = prep_systems_UCED(base_da_system)

            # generic market formulation templates with defined network formulation
            # CopperPlate-OPF: network=CopperPlatePowerModel
            # DC-OPF: network=DCPPowerModel
            # NFA-OPF (only line limit constraints): network=NFAPowerModel
            # DC-PTDF-OPF (what ISOs do): network=StandardPTDFModel
            template_uc = template_unit_commitment(; network=df.Network[l]["DA"])
            template_rt = template_economic_dispatch(; network=df.Network[l]["RT"])
            template_ed = template_economic_dispatch(; network=df.Network[l]["DA"])

            # build a market clearing simulator
            market_simulator = UCEDRT(;
                system_uc=sys_uc,
                system_rt=sys_rt,
                system_ed=sys_ed,
                template_uc=template_uc,
                template_rt=template_rt,
                template_ed=template_ed,
                solver_uc=solver_uc,
                solver_rt=solver_rt,
                solver_ed=solver_ed,
            )

            #Calculates the dispatch result for a bid curve
            name_load = get_name(load);
            name_gen = get_name(gen);
            steps = 1;

            mkdir(joinpath(example_dir, path, directory_name))
            simulation_folder = joinpath(example_dir, path, directory_name)

            lmps_df, results_df = pq_curves_load_gen_virtuals!(
                market_simulator, name_load, range_quota_load, initial_time, name_gen, range_quota_gen, steps, simulation_folder
            )

        end
    end
end

"""
    load_set_of_simulations(
        df::DataFrame, 
        data_dir::String, 
        example_dir::String, 
        solver_uc, 
        solver_ed, 
        solver_rt, 
        range_quota::Vector{Float64}, 
        lines::Vector{Int64},
        initial_bidding_time::DateTime,
        path::String,
    )

Load a set of simulations of 5bus_nrel with descripted sistems in `lines` from `df`
"""
function load_set_of_simulations(
    df::DataFrame, 
    data_dir::String, 
    example_dir::String, 
    solver_uc, 
    solver_ed, 
    solver_rt, 
    range_quota::Vector{Float64}, 
    lines::Vector{Int64},
    initial_bidding_time::DateTime,
    path::String,
)
    #=
    if graphic == "plot_price_curves" 
        global plt=Array{Any}(nothing, length(lines), length(period_analysed))
    elseif graphic == "plot_generation_stack_virtual"
        global plt=Array{Any}(nothing, length(lines), length(period_analysed),2)
    elseif graphic == "plot_revenue_curves_renewable_plus_virtual" || graphic == "plot_revenue_curves"
        global plt=Array{Any}(nothing, length(lines))
    end
    =#
    global lmps_df=Array{Any}(nothing, length(lines))
    global results_df=Array{Any}(nothing, length(lines))
    for (x,l) in enumerate(lines)

        # call our data preparation to build base system
        # the case was modified to not have hydros nor transformers
        sys_rt = build_5_bus_matpower_RT(data_dir;)

        base_da_system = build_5_bus_matpower_DA(
            data_dir;
            # using a modified (mod) file that reduced load for feasibility in DC-OPF
            forecasts_pointers_file=joinpath(
                data_dir, "forecasts", "timeseries_pointers_da_7day_mod.json"
            ),
            add_reserves=df.Reserve[l],
        )

        # Add single generator at a defined bus
        gen = add_generator!(base_da_system, df.Offer_Bus[l], (min=0.0, max=0.0))
        @test gen in get_components(Generator, base_da_system)

        # create and set variable cost time-series for the generator
        ts_array = create_generator_bids(;
            initial_bidding_time=initial_bidding_time,
            bidding_periods=df.bidding_period[l],
            system=base_da_system,
            costs=zeros(length(df.bidding_period[l])),
        )
        set_variable_cost!(base_da_system, gen, ts_array)

        # duplicate system and prepare times series for the time varying parameters (loads, renewables, ...)
        sys_uc, sys_ed = prep_systems_UCED(base_da_system)

        # generic market formulation templates with defined network formulation
        # CopperPlate-OPF: network=CopperPlatePowerModel
        # DC-OPF: network=DCPPowerModel
        # NFA-OPF (only line limit constraints): network=NFAPowerModel
        # DC-PTDF-OPF (what ISOs do): network=StandardPTDFModel
        template_uc = template_unit_commitment(; network=df.Network[l]["DA"])
        template_rt = template_economic_dispatch(; network=df.Network[l]["RT"])
        template_ed = template_economic_dispatch(; network=df.Network[l]["DA"])

        # build a market clearing simulator
        market_simulator = UCEDRT(;
            system_uc=sys_uc,
            system_rt=sys_rt,
            system_ed=sys_ed,
            template_uc=template_uc,
            template_rt=template_rt,
            template_ed=template_ed,
            solver_uc=solver_uc,
            solver_rt=solver_rt,
            solver_ed=solver_ed,
        )

        directory_name= "Net_"*string(df.Network[l]["DA"])[1:4]*"_"*string(df.Network[l]["RT"])[1:4]* "__Ramp_"*string(df.Ramp[l]["DA"])*
        "_"*string(df.Ramp[l]["RT"])*"__Min_gen_"*string(df.Minimal_generation[l]["DA"])*"_"*string(df.Minimal_generation[l]["RT"])* "__Reserve_" * string(df.Reserve[l]) *
        "__Offer_" * string(df.Offer_Bus[l]) * "__Bid_period_1-" * string(length(df.bidding_period[l]))

        simulation_folder = joinpath(example_dir, path, directory_name)

        lmps_df[l], results_df[l] = load_pq_curves(market_simulator, range_quota, simulation_folder)
        #=
        if graphic == "plot_price_curves" 
            for (y,t) in enumerate(period_analysed)
                global plt[x,y] = plot_price_curves(lmps_df, period_analysed[y], unique(df.Offer_Bus), df.Offer_Bus[l], initial_time)
            end
        elseif graphic == "plot_generation_stack_virtual"
            for (y,t) in enumerate(period_analysed)
                global plt[x,y,1] = plot_generation_stack_virtual(sys_uc, results_df; type="DA", period=period_analysed[y], initial_time, xtickfontsize=8, margin=8mm, size=(800, 600),)
                global plt[x,y,2] = plot_generation_stack_virtual(sys_rt, results_df; type="RT", period=period_analysed[y], initial_time, xtickfontsize=8, margin=8mm, size=(800, 600),)
            end
        elseif graphic == "plot_revenue_curves_renewable_plus_virtual"
            global plt[x] = plot_revenue_curves_renewable_plus_virtual(market_simulator, lmps_df, results_df, [0.0, 1.0, 2.0],"WindBusA", df.Offer_Bus[l]*"_virtual_supply",)
        elseif graphic == "plot_revenue_curves"
            period=[period_analysed[i][1] for i=1:length(period_analysed)]
            global plt[x] = plot_revenue_curves(
                market_simulator, lmps_df, results_df, period, df.Offer_Bus[l]*"_virtual_supply", initial_time
            )
        end
        =#
    end
    return lmps_df, results_df
end

"""
    load_plot_set_of_simulations(
        df::DataFrame, 
        data_dir::String, 
        example_dir::String, 
        solver_uc, 
        solver_ed, 
        solver_rt, 
        range_quota::Vector{Float64}, 
        initial_time::Date,
        lines::Vector{Int64},
        period_analysed:: Vector{Vector{Int64}},
        initial_bidding_time::DateTime,
        path::String,
        graphic::String,
    )

Load and plot a set of simulations of 5bus_nrel with descripted sistems in `lines` from `df`
"""
function load_plot_set_of_simulations(
    df::DataFrame, 
    data_dir::String, 
    example_dir::String, 
    solver_uc, 
    solver_ed, 
    solver_rt, 
    range_quota::Vector{Float64}, 
    initial_time::Date,
    lines::Vector{Int64},
    period_analysed:: Vector{Vector{Int64}},
    initial_bidding_time::DateTime,
    path::String,
    graphic::String,
    bool::Bool,
    virtual_type::String="",
)

    if graphic == "plot_price_curves" 
        global plt=Array{Any}(nothing, length(lines), length(period_analysed))
    elseif graphic == "plot_generation_stack_virtual"
        global plt=Array{Any}(nothing, length(lines), length(period_analysed),2)
    elseif graphic == "plot_revenue_curves_renewable_plus_virtual" || graphic == "plot_revenue_curves" || graphic =="plot_sum_revenue_curves"
        global plt=Array{Any}(nothing, length(lines))
    elseif graphic == "plot_thermal_commit_virtual"
        global plt=Array{Any}(nothing, length(lines), length(period_analysed))
    end
    
    for (x,l) in enumerate(lines)

        # call our data preparation to build base system
        # the case was modified to not have hydros nor transformers
        sys_rt = build_5_bus_matpower_RT(data_dir;)

        base_da_system = build_5_bus_matpower_DA(
            data_dir;
            # using a modified (mod) file that reduced load for feasibility in DC-OPF
            forecasts_pointers_file=joinpath(
                data_dir, "forecasts", "timeseries_pointers_da_7day_mod.json"
            ),
            add_reserves=df.Reserve[l],
        )

        # Add single generator at a defined bus
        gen = add_generator!(base_da_system, df.Offer_Bus[l], (min=0.0, max=0.0))

        # create and set variable cost time-series for the generator
        ts_array = create_generator_bids(;
            initial_bidding_time=initial_bidding_time,
            bidding_periods=df.bidding_period[l],
            system=base_da_system,
            costs=zeros(length(df.bidding_period[l])),
        )
        set_variable_cost!(base_da_system, gen, ts_array)

        # duplicate system and prepare times series for the time varying parameters (loads, renewables, ...)
        sys_uc, sys_ed = prep_systems_UCED(base_da_system)

        # generic market formulation templates with defined network formulation
        # CopperPlate-OPF: network=CopperPlatePowerModel
        # DC-OPF: network=DCPPowerModel
        # NFA-OPF (only line limit constraints): network=NFAPowerModel
        # DC-PTDF-OPF (what ISOs do): network=StandardPTDFModel
        template_uc = template_unit_commitment(; network=df.Network[l]["DA"])
        template_rt = template_economic_dispatch(; network=df.Network[l]["RT"])
        template_ed = template_economic_dispatch(; network=df.Network[l]["DA"])

        # build a market clearing simulator
        market_simulator = UCEDRT(;
            system_uc=sys_uc,
            system_rt=sys_rt,
            system_ed=sys_ed,
            template_uc=template_uc,
            template_rt=template_rt,
            template_ed=template_ed,
            solver_uc=solver_uc,
            solver_rt=solver_rt,
            solver_ed=solver_ed,
        )

        directory_name= "Net_"*string(df.Network[l]["DA"])[1:4]*"_"*string(df.Network[l]["RT"])[1:4]* "__Ramp_"*string(df.Ramp[l]["DA"])*
        "_"*string(df.Ramp[l]["RT"])*"__Min_gen_"*string(df.Minimal_generation[l]["DA"])*"_"*string(df.Minimal_generation[l]["RT"])* "__Reserve_" * string(df.Reserve[l]) *
        "__Offer_" * string(df.Offer_Bus[l]) * "__Bid_period_1-" * string(length(df.bidding_period[l]))

        simulation_folder = joinpath(example_dir, path, directory_name)

        lmps_df, results_df = load_pq_curves(market_simulator, range_quota, simulation_folder)
        
        if graphic == "plot_price_curves" 
            for (y,t) in enumerate(period_analysed)
                global plt[x,y] = plot_price_curves(lmps_df, period_analysed[y], unique(df.Offer_Bus), df.Offer_Bus[l], initial_time, sys_uc, bool, virtual_type)
            end
        elseif graphic == "plot_generation_stack_virtual"
            for (y,t) in enumerate(period_analysed)
                global plt[x,y,1] = plot_generation_stack_virtual(sys_uc, results_df; type="DA", period=period_analysed[y], initial_time, xtickfontsize=8, margin=8mm, size=(800, 600),)
                global plt[x,y,2] = plot_generation_stack_virtual(sys_rt, results_df; type="RT", period=period_analysed[y], initial_time, xtickfontsize=8, margin=8mm, size=(800, 600),)
            end
        elseif graphic == "plot_thermal_commit_virtual"
            for (y,t) in enumerate(period_analysed)
                global plt[x,y] = plot_thermal_commit_virtual(sys_uc, results_df; period=period_analysed[y], initial_time, xtickfontsize=8, margin=8mm, size=(800, 600),)
            end
        elseif graphic == "plot_revenue_curves_renewable_plus_virtual"
            global plt[x] = plot_revenue_curves_renewable_plus_virtual(market_simulator, lmps_df, results_df, [0.0, 1.0, 2.0],"WindBusA", df.Offer_Bus[l]*"_virtual_supply", bool)
        elseif graphic == "plot_revenue_curves"
            period=[period_analysed[i][1] for i=1:length(period_analysed)]
            global plt[x] = plot_revenue_curves(
                market_simulator, lmps_df, results_df, period, df.Offer_Bus[l]*"_virtual_supply", initial_time, bool
            )
        elseif graphic == "plot_sum_revenue_curves"
            period=[period_analysed[i][1] for i=1:length(period_analysed)]
            global plt[x] = plot_sum_revenue_curves(
                market_simulator, lmps_df, results_df, period, df.Offer_Bus[l]*"_virtual_supply", initial_time
            )
        end
        
    end
    return plt
end

"""
    load_plot_set_of_simulations_mix(
        df::DataFrame, 
        data_dir::String, 
        example_dir::String, 
        solver_uc, 
        solver_ed, 
        solver_rt, 
        range_quota_load::Vector{Float64},
        range_quota_gen::Vector{Float64},  
        initial_time::Date,
        lines::Vector{Int64},
        initial_bidding_time::DateTime,
        path::String,
        graphic::String,
    )

Load and plot a set of simulations of 5bus_nrel with descripted sistems in 'lines' from 'df'
"""
function load_plot_set_of_simulations_mix(
    df::DataFrame, 
    data_dir::String, 
    example_dir::String, 
    solver_uc, 
    solver_ed, 
    solver_rt, 
    range_quota_load::Vector{Float64}, 
    range_quota_gen::Vector{Float64},
    initial_time::Date,
    lines::Vector{Int64},
    initial_bidding_time::DateTime,
    path::String,
)
    global plt = Dict()
    global h = Dict()
   
    for l in lines

        # call our data preparation to build base system
        # the case was modified to not have hydros nor transformers
        sys_rt = build_5_bus_matpower_RT(data_dir;)

        base_da_system = build_5_bus_matpower_DA(
            data_dir;
            # using a modified (mod) file that reduced load for feasibility in DC-OPF
            forecasts_pointers_file=joinpath(
                data_dir, "forecasts", "timeseries_pointers_da_7day_mod.json"
            ),
            add_reserves=df.Reserve[l],
        )

        # create and set variable time-series for the load
        ts_array = create_demand_series(;
            initial_bidding_time=initial_bidding_time,
            bidding_periods=df.bidding_period[l],
            system=base_da_system,
            demands=ones(length(df.bidding_period[l])),
        )

        # Add single load at a defined bus
        node_load = df.Load_Bus[l] # define bus
        load = add_load!(base_da_system, node_load, 1.0)

        add_time_series!(base_da_system, load, ts_array)

        # Add single generator at a defined bus
        gen = add_generator!(base_da_system, df.Offer_Bus[l], (min=0.0, max=0.0))

        # create and set variable cost time-series for the generator
        ts_array = create_generator_bids(;
            initial_bidding_time=initial_bidding_time,
            bidding_periods=df.bidding_period[l],
            system=base_da_system,
            costs=zeros(length(df.bidding_period[l])),
        )
        set_variable_cost!(base_da_system, gen, ts_array)
        # duplicate system and prepare times series for the time varying parameters (loads, renewables, ...)
        sys_uc, sys_ed = prep_systems_UCED(base_da_system)

        # generic market formulation templates with defined network formulation
        # CopperPlate-OPF: network=CopperPlatePowerModel
        # DC-OPF: network=DCPPowerModel
        # NFA-OPF (only line limit constraints): network=NFAPowerModel
        # DC-PTDF-OPF (what ISOs do): network=StandardPTDFModel
        template_uc = template_unit_commitment(; network=df.Network[l]["DA"])
        template_rt = template_economic_dispatch(; network=df.Network[l]["RT"])
        template_ed = template_economic_dispatch(; network=df.Network[l]["DA"])

        # build a market clearing simulator
        market_simulator = UCEDRT(;
            system_uc=sys_uc,
            system_rt=sys_rt,
            system_ed=sys_ed,
            template_uc=template_uc,
            template_rt=template_rt,
            template_ed=template_ed,
            solver_uc=solver_uc,
            solver_rt=solver_rt,
            solver_ed=solver_ed,
        )

        directory_name= df.directory_name[l]

        simulation_folder = joinpath(example_dir, path, directory_name)

        lmps_df, results_df = load_mix_pq_curves(market_simulator, range_quota_load, range_quota_gen, simulation_folder)

        plt[l]=Dict()
        h[l]=Dict()

        plt[l]["deficit"],h[l]["deficit"],sum_deficit=heat_map_deficit(
        sys_rt,
        results_df,
        range_quota_load,
        range_quota_gen,
        initial_time,
        df.bidding_period[l],
        )

        plt[l]["coal"],h[l]["coal"]=heat_map_coal_generation(
        sys_uc,
        results_df,
        range_quota_load,
        range_quota_gen,
        initial_time,
        sum_deficit,
        df.bidding_period[l],
        [:P__ThermalStandard, :P__RenewableDispatch],
        )

        plt[l]["revenue - solar"],h[l]["revenue - solar"]=heat_map_revenue_curves_mix(
            market_simulator,
            lmps_df,
            results_df,
            df.bidding_period[l],
            range_quota_load,
            range_quota_gen,
            initial_time,
            load,
            get_name(gen),
            "SolarBusC",
            sys_uc,
        )

        plt[l]["revenue - wind"],h[l]["revenue - wind"]=heat_map_revenue_curves_mix(
            market_simulator,
            lmps_df,
            results_df,
            df.bidding_period[l],
            range_quota_load,
            range_quota_gen,
            initial_time,
            load,
            get_name(gen),
            "WindBusA",
            sys_uc,
        )
        
        
    end
    return plt,h
end

