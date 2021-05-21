
# Function to plot the Renewable Dispatch

@userplot plot_renweable_generation_stack
@recipe function f(p::plot_renweable_generation_stack)
    system, system_results, = p.args

    # get mapping from busname to fuel type
    fuel_type_dict = fuel_type_mapping(system)

    # get the output data for all fuel types
    variable_results = read_realized_variables(system_results)
    renewable_results = variable_results[:P__RenewableDispatch]

    # stack the data and aggregate by fuel type
    stacked_data = stack(renewable_results, variable_name = "bus_name", value_name = "output")
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

plot_renweable_generation_stack(base_system, ed_results, xtickfontsize = 8, size = (800, 600))


# Function to plot the Thermal Standard

@userplot plot_thermal_generation_stack
@recipe function f(p::plot_thermal_generation_stack)
    system, system_results, = p.args

    # get mapping from busname to fuel type
    fuel_type_dict = fuel_type_mapping(system)

    # get the output data for all fuel types
    variable_results = read_realized_variables(system_results)
    thermal_results = variable_results[:P__ThermalStandard]

    # stack the data and aggregate by fuel type
    stacked_data = stack(thermal_results, variable_name = "bus_name", value_name = "output")
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

plot_thermal_generation_stack(base_system, ed_results, xtickfontsize = 8, size = (800, 600))


# Function to plot the Thermal Standard Commit
# It needs the uc_results.

@userplot plot_thermal_commit_stack
@recipe function f(p::plot_thermal_commit_stack)
    system, system_results, = p.args

    # get mapping from busname to fuel type
    fuel_type_dict = fuel_type_mapping(system)

    # get the output data for all fuel types
    variable_results = read_realized_variables(system_results)
    thermal_commit_results = variable_results[:On__ThermalStandard]

    # stack the data and aggregate by fuel type
    stacked_data = stack(thermal_commit_results, variable_name = "bus_name", value_name = "output")
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

plot_thermal_commit_stack(base_system, uc_results, xtickfontsize = 8, size = (800, 600))


# Function to plot the Prices

@userplot plot_prices_stack
@recipe function f(p::plot_prices_stack)

    prices = evaluate_prices(template_ed, ed_results)

    plot_data = select(prices, Not(:DateTime))

    times = prices[!,1]

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

plot_prices_stack(base_system, ed_results, xtickfontsize = 8, size = (800, 600))


# Function to plot the Demmand
# It needs the uc_results.

@userplot plot_demand_stack
@recipe function f(p::plot_demand_stack)
    system, system_results, = p.args

    load = collect(get_components(PowerLoad, system))

    ts_array = zeros(24,length(load))
    ts_names = get_time_series_names(Deterministic, load[1])
    for i in 1:length(load)
        ts_array[:,i] = get_time_series_values(Deterministic, load[i], ts_names[1])
    end

    times = get_time_series_timestamps(Deterministic, load[1], ts_names[1])

    plot_data = DataFrame(ts_array, [:Bus4, :Bus2, :Bus3])

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

plot_demand_stack(sys_uc, uc_results, xtickfontsize = 8, size = (800, 600))


# Function to plot the Net Demmand
# It needs the sys_uc and uc_results.

@userplot plot_net_demand_stack
@recipe function f(p::plot_net_demand_stack)
    system, system_results, = p.args

    load = collect(get_components(PowerLoad, system))

    ts_array = zeros(24,length(load))
    ts_names = get_time_series_names(Deterministic, load[1])
    for i in 1:length(load)
        ts_array[:,i] = get_time_series_values(Deterministic, load[i], ts_names[1])
    end

    # Renewable Data
    variable_results = read_realized_variables(system_results)
    renewable_results = variable_results[:P__RenewableDispatch]
    row, col = size(renewable_results)
    values = renewable_results[!,2:col]

    all_data = sum(ts_array, dims = 2) - sum.(eachrow(values))

    times = get_time_series_timestamps(Deterministic, load[1], ts_names[1])

    plot_data = DataFrame(all_data, [:net_demand])

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

plot_net_demand_stack(sys_uc, uc_results, xtickfontsize = 8, size = (800, 600))

