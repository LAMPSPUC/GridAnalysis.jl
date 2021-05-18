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
