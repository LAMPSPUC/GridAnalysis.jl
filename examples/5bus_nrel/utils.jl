"""
    build_5_bus_matpower_DA(data_dir::AbstractString, case_file::AbstractString)

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
    pm_data = PowerSystems.PowerModelsData(case_file_path)

    tsp = InfrastructureSystems.read_time_series_file_metadata(forecasts_pointers_file)

    sys = System(pm_data)
    if add_reserves
        reserves = [
            VariableReserve{ReserveUp}("REG1", true, 5.0, 0.1),
            VariableReserve{ReserveUp}("REG2", true, 5.0, 0.06),
            VariableReserve{ReserveUp}("REG3", true, 5.0, 0.03),
            VariableReserve{ReserveUp}("REG4", true, 5.0, 0.02),
        ]
        contributing_devices = get_components(Generator, sys)
        for r in reserves
            add_service!(sys, r, contributing_devices)
        end
    end

    add_time_series!(sys, tsp)

    return sys
end

"""
    prep_systems_UCED(system::System)

Duplicates the system to represent UC and ED for DA, transforming the time series 
to the appropriate interval and horizon.

PS.: Beacuse of a bug in PSI, we have to set the horizon of the ED problem, that 
would normally be of 1 hour, to 2 hours.
"""
function prep_systems_UCED(
    system::System;
    horizon_uc::Int=24,
    horizon_ed::Int=2,
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
    get_duals(<:AbstractPowerModel)

Return the duals for the specified formulation.
"""
get_duals(::Type{CopperPlatePowerModel}) = [:CopperPlateBalance]
get_duals(::Union{Type{NFAPowerModel}, Type{DCPPowerModel}}) = [:nodal_balance_active__Bus]
get_duals(::Type{StandardPTDFModel}) = [:CopperPlateBalance, :network_flow__Line]



"""
    fuel_type_mapping(system)

    return a mapping between the bus name and the fuel type for the given `system`.
"""
function fuel_type_mapping(system::System)
    generator_metadata =  [gen for gen in get_components(Generator, system)]

    fuel_type_map = Dict()
    for generator in generator_metadata 
        name = generator.name
        try 
            fuel_type_map[name] =  get_fuel(generator)
        catch
            # Assumes the bus name has format "<Fuel>Bus..."
            fuel_type_map[name] = first(split(name, "Bus"))
        end
    end

    return fuel_type_map
end

############# Plots #############
"""
    plot_generation_stack(
        system::System,
        results::SimulationProblemResults)

Plot the generation mix over the time period covered by the `results`. 
"""
@userplot plot_generation_stack
@recipe function f(p::plot_generation_stack)
    system, system_results, = p.args

    # get mapping from busname to fuel type
    fuel_type_dict = fuel_type_mapping(system)

    # get the output data for all fuel types
    variable_results = read_realized_variables(system_results)
    thermal_results = variable_results[:P__ThermalStandard]
    renewable_results = variable_results[:P__RenewableDispatch]

    all_generator_data = innerjoin(thermal_results, renewable_results, on = :DateTime)

    # stack the data and aggregate by fuel type
    stacked_data = stack(all_generator_data, variable_name = "bus_name", value_name = "output")
    stacked_data.fuel_type = get.(Ref(fuel_type_dict), stacked_data.bus_name, missing)

    aggregated_data = combine(
        groupby(stacked_data, [:DateTime, :fuel_type,]),
        :output => sum => :output
    )

    # unstack aggregated data and make area plot 
    unstacked_data = unstack(aggregated_data, :fuel_type, :output)
    times = unstacked_data.DateTime

    plot_data = select(unstacked_data, Not(:DateTime))

    label --> reduce(hcat, names(plot_data))
    yguide --> "Output (MWh)"
    legend --> :outertopright
    seriestype --> :line
    xrotation --> 45
    
    
    # now stack the matrix to get the cumulative values over all fuel types
    data = cumsum(Matrix(plot_data), dims = 2)
    for i in Base.axes(data, 2)
        @series begin
            fillrange := i > 1 ? data[:, i - 1] : 0
            times, data[:, i]
        end
    end
end
