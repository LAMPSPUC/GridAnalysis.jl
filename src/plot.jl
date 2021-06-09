"""
    plot_generation_stack(
        system::System,
        results::SimulationProblemResults;
        generator_fields::AbstractArray=[:P__ThermalStandard, :P__RenewableDispatch],
        bus_names::AbstractArray=[])

Plot the generation mix over the time period covered by the `results`. The `bus_names`
and `generator_fields` control which buses, and generator types we want to include the plot.
"""
@userplot plot_generation_stack
@recipe function f(
    p::plot_generation_stack;
    generator_fields::AbstractArray=[:P__ThermalStandard, :P__RenewableDispatch],
    bus_names::AbstractArray=[],
)
    system, system_results, = p.args

    # get mapping from busname to fuel type
    fuel_type_dict = fuel_type_mapping(system)

    # get the output data for given fuel types
    variable_results = read_realized_variables(system_results; names=generator_fields)
    generator_data = getindex.(Ref(variable_results), generator_fields)
    if length(generator_data) > 1
        generator_data = innerjoin(generator_data...; on=:DateTime)
    else
        generator_data = first(generator_data)
    end

    # stack the data and aggregate by fuel type
    stacked_data = stack(generator_data; variable_name="gen_name", value_name="output")
    bus_map = bus_mapping(system)
    stacked_data.bus_name = [bus_map[gen] for gen in stacked_data.gen_name]

    # select rows for the given bus names, default to all buses.
    if !isempty(bus_names)
        bus_names = String.(bus_names)
        @assert all(bus_names .∈ [stacked_data.bus_name])
        filter!(:bus_name => in(bus_names), stacked_data)
    end

    stacked_data.fuel_type = get.(Ref(fuel_type_dict), stacked_data.gen_name, missing)

    aggregated_data = combine(
        groupby(stacked_data, [:DateTime, :fuel_type]), :output => sum => :output
    )

    # convert the output units into MWh.
    aggregated_data.output = aggregated_data.output .* get_base_power(system)

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
    data = cumsum(Matrix(plot_data); dims=2)
    for i in Base.axes(data, 2)
        @series begin
            fillrange := i > 1 ? data[:, i - 1] : 0
            times, data[:, i]
        end
    end
end

"""
    plot_prices(
        market_simulator::MarketSimulator,
        results::SimulationProblemResults;
        bus_names::AbstractArray=[])

Plot the simulation prices over the time period covered by the `results`. The `bus_names`
control which buses we want to include the plot.
"""
@userplot plot_prices
@recipe function f(p::plot_prices; bus_names::AbstractArray=[])
    market_simulator, system_results, = p.args

    prices = evaluate_prices(market_simulator, system_results)

    plot_data = select(prices, Not(:DateTime))

    # select rows for the given bus names, default to all buses.
    if !isempty(bus_names)
        bus_names = String.(bus_names)
        @assert issubset(bus_names, names(plot_data))
        select!(plot_data, bus_names)
    end

    times = prices[!, 1]

    if isa(market_simulator, UCRT) || isa(market_simulator, UCEDRT)
        yguide --> "Prices (\$/MW-5min)"
    else
        yguide --> "Prices (\$/MWh)"
    end

    label --> reduce(hcat, names(plot_data))
    legend --> :outertopright
    seriestype --> :line
    xrotation --> 45

    for i in Base.axes(plot_data, 2)
        @series begin
            times, plot_data[:, i]
        end
    end
end

"""
    plot_prices(
        market_simulator::MarketSimulator,
        results::Tuple{SimulationResults, SimulationResults};
        bus_names::AbstractArray=[])

Plot the simulation prices over the time period covered by the `results`. The `bus_names`
control which buses we want to include the plot.
"""
@userplot plot_prices
@recipe function f(p::plot_prices; bus_names::AbstractArray=[])
    market_simulator, system_results, = p.args

    prices = evaluate_prices_UCEDRT(market_simulator, system_results)

    plot_data_ed = select(prices["ED"], Not(:DateTime))
    prices_rt = prices["RT"]
    values = select(prices_rt, Not(:DateTime))
    n_prev, n_bus = size(values)

    intervals = get_time_series_params(market_simulator.system_rt).interval
    
    n_prev_hour = Int(60/intervals.value)
    n_days = Int(n_prev/n_prev_hour)
    
    names_bus = names(values)
    prices_rt_sum = zeros(n_days,n_bus)

    i = 1
    while i < n_days
        for j in 1:(n_prev_hour):n_prev
            prices_hour = prices_rt[prices_rt[j,:DateTime] .<= prices_rt.DateTime .< prices_rt[j,:DateTime]+Hour(1), :]
            prices_hour = select(prices_hour, Not(:DateTime))
            prices_rt_sum[i,:] = sum(Matrix(prices_hour), dims = 1)
            i = i + 1
        end
    end
    
    dictionary = Dict()
    n_row, n_col = size(values)

    for price in 1:n_col
        if !haskey(dictionary, names_bus[price])
            dictionary[names_bus[price]] = prices_rt_sum[:,price]
        else
            dictionary[names_bus[price]] = dictionary[names_bus[price]] + prices_rt_sum[:,price]
        end
    end

    plot_data = DataFrame(dictionary)

    # select rows for the given bus names, default to all buses.
    if !isempty(bus_names)
        bus_names = String.(bus_names)
        @assert issubset(bus_names, names(plot_data))
        select!(plot_data, bus_names)
    end

    times = prices["ED"][!, 1]

    label --> reduce(hcat, names(plot_data))
    yguide --> "Prices (\$/MWh)"
    legend --> :outertopright
    seriestype --> :line
    xrotation --> 45

    for i in Base.axes(plot_data, 2)
        @series begin
            times, plot_data[:, i]
        end
    end
end


"""
    plot_prices_RT(
        market_simulator::MarketSimulator,
        results::SimulationResults;
        bus_names::AbstractArray=[])

Plot the simulation prices over the time period covered by the `results`. The `bus_names`
control which buses we want to include the plot.
"""
@userplot plot_prices_RT
@recipe function f(p::plot_prices_RT; bus_names::AbstractArray=[])
    market_simulator, system_results, = p.args

    prices = evaluate_prices(market_simulator, system_results)

    values = select(prices, Not(:DateTime))
    n_prev, n_bus = size(values)

    intervals = get_time_series_params(market_simulator.system_rt).interval
    
    n_prev_hour = Int(60/intervals.value)
    n_days = Int(n_prev/n_prev_hour)
    
    names_bus = names(values)
    prices_rt = zeros(n_days,n_bus)

    i = 1
    while i < n_days
        #for j in 1:(n_prev_hour):n_prev
        for j = i+(12*(i-1))
            prices_hour = prices[prices[j,:DateTime] .<= prices.DateTime .< prices[j,:DateTime]+Hour(1), :]
            prices_hour = select(prices_hour, Not(:DateTime))
            prices_rt[i,:] = sum(Matrix(prices_hour), dims = 1)
            i = i + 1
        end
    end

    times = prices[1:n_prev_hour:n_prev, 1]
    
    dictionary = Dict()
    n_row, n_col = size(prices_rt)

    for price in 1:n_col
        if !haskey(dictionary, names_bus[price])
            dictionary[names_bus[price]] = prices_rt[:,price]
        else
            dictionary[names_bus[price]] = dictionary[names_bus[price]] + prices_rt[:,price]
        end
    end

    plot_data = DataFrame(dictionary)

    # select rows for the given bus names, default to all buses.
    if !isempty(bus_names)
        bus_names = String.(bus_names)
        @assert issubset(bus_names, names(plot_data))
        select!(plot_data, bus_names)
    end

    label --> reduce(hcat, names(plot_data))
    yguide --> "Prices (\$/MWh)"
    legend --> :outertopright
    seriestype --> :line
    xrotation --> 45

    for i in Base.axes(plot_data, 2)
        @series begin
            times, plot_data[:, i]
        end
    end
end


"""
    plot_thermal_commit(
        system::System,
        results::SimulationProblemResults)

Function to plot the Thermal Standard Commit variables over the time period covered by the `results`.
The `results` should be from the unit commitment problem.
1 is ON, 0 is OFF.
"""
@userplot plot_thermal_commit
@recipe function f(p::plot_thermal_commit;
    bus_names::AbstractArray=[],
)
    system, system_results, = p.args

    # get the output data for all fuel types
    variable_results = read_realized_variables(system_results)
    thermal_commit_results = variable_results[:On__ThermalStandard]

    plot_data = select(thermal_commit_results, Not(:DateTime))

    times = thermal_commit_results[!, 1]

    # select rows for the given bus names, default to all buses.
    if !isempty(bus_names)
        generator_names = names(plot_data)
        bus_map = bus_mapping(system)
        bus_names = String.(bus_names)
        @assert issubset(bus_names, [bus_map[gen] for gen in generator_names])
        generator_names = [gen for gen in generator_names if in(bus_map[gen], bus_names)]
        select!(plot_data, generator_names)
    end

    label --> reduce(hcat, names(plot_data))
    yguide --> "Commitment"
    legend --> :outertopright
    seriestype --> :line
    xrotation --> 45
    title --> "Thermal Standard Commit over the hours"

    for i in Base.axes(plot_data, 2)
        @series begin
            times, plot_data[:, i]
        end
    end
end


"""
    plot_demand_stack(
        system::System)

Function to plot the Demand over the time period covered by the `results`.
The `bus_names` controls which buses we want to include in the plot.
"""
@userplot plot_demand_stack
@recipe function f(
    p::plot_demand_stack;
    bus_names::AbstractArray=[],
)
    system, = p.args

    loads = collect(get_components(PowerLoad, system))

    ts_array = Dict()
    ts_names = get_time_series_names(SingleTimeSeries, loads[1])
    for load in loads
        if !haskey(ts_array, get_bus_name(load))
            ts_array[get_bus_name(load)] = get_time_series_values(SingleTimeSeries, load, ts_names[1])
        else
            ts_array[get_bus_name(load)] = ts_array[get_bus_name(load)] .+ get_time_series_values(SingleTimeSeries, load, ts_names[1])
        end
    end

    ts_array = DataFrame(ts_array)

    times = get_time_series_timestamps(SingleTimeSeries, loads[1], ts_names[1])
    
    plot_data = ts_array .* get_base_power(system)

    if !isempty(bus_names)
        bus_names = String.(bus_names)
        @assert issubset(bus_names, names(plot_data))
        select!(plot_data, bus_names)
    end

    label --> reduce(hcat, names(plot_data))
    yguide --> "Demand (MWh)"
    legend --> :outertopright
    seriestype --> :line
    xrotation --> 45
    title --> "Demand over the hours"

    # now stack the matrix to get the cumulative values over all fuel types
    data = cumsum(Matrix(plot_data); dims=2)
    for i in Base.axes(data, 2)
        @series begin
            fillrange := i > 1 ? data[:, i - 1] : 0
            times, data[:, i]
        end
    end
end


"""
    plot_net_demand_stack_prev(
        system::System,
        results::SimulationProblemResults)
Function to plot the Demand over the time period covered by the `results`.
The `bus_names` controls which buses we want to include in the plot.
It uses the Renewable Dispatch from the `results`.
"""
@userplot plot_net_demand_stack_prev
@recipe function f(p::plot_net_demand_stack_prev;
    x_ticks::StepRange{Int64, Int64}=nothing,
    bus_names::AbstractArray=[],
)
    
    system, system_results, = p.args
    
    loads = collect(get_components(PowerLoad, system))
    
    ts_array = Dict()
    ts_names = get_time_series_names(Deterministic, loads[1])
    
    for load in loads
        if !haskey(ts_array, get_bus_name(load))
            ts_array[get_bus_name(load)] = get_time_series_values(Deterministic, load, ts_names[1])
        else
            ts_array[get_bus_name(load)] = ts_array[get_bus_name(load)] .+ get_time_series_values(Deterministic, load, ts_names[1])
        end
    end
    
    ts_array = DataFrame(ts_array)
    
    # Renewable Data
    variable_results = read_realized_variables(system_results)
    renewable_results = variable_results[:P__RenewableDispatch]
    values = select(renewable_results, Not(:DateTime))
    
    bus_map = bus_mapping(system)
    rename!(values, [bus_map[gen] for gen in names(values)])
    
    # select rows for the given bus names, default to all buses.
    if !isempty(bus_names)
        bus_names = String.(bus_names)
        @assert issubset(bus_names, [names(values); names(ts_array)])
        select!(values, intersect(bus_names, names(values)))
        select!(ts_array, intersect(bus_names, names(ts_array)))
    end
    
    all_data = sum.(eachrow(ts_array)) - sum.(eachrow(values))
    times = get_time_series_timestamps(Deterministic, loads[1], ts_names[1])
    hours = zeros(length(times))
    for h in 1:length(times)
        hours[h] = Dates.hour(times[h])+1
    end
    plot_data = DataFrame(net_demand = all_data) .* get_base_power(system)
    x_ticks = isnothing(x_ticks) ? (0:1:length(times)) : x_ticks
    
    label --> reduce(hcat, names(plot_data))
    yguide --> "Net Demand (MWh)"
    legend --> :outertopright
    seriestype --> :line
    xrotation --> 0
    title --> "Net Demand over the hours"
    xticks --> x_ticks
    
    # now stack the matrix to get the cumulative values over all fuel types
    data = cumsum(Matrix(plot_data); dims=2)
    for i in Base.axes(data, 2)
        @series begin
            fillrange := i > 1 ? data[:, i - 1] : 0
            hours, data[:, i]
        end
    end
end


"""
    plot_net_demand_stack(
        system::System)
Function to plot the Demand over the time period covered by the `results`.
The `bus_names` controls which buses we want to include in the plot.
Renewable Dispatch data is the time series from the `system`.
"""
@userplot plot_net_demand_stack
@recipe function f(p::plot_net_demand_stack;
    bus_names::AbstractArray=[],
)
    
    system, = p.args
    
    loads = collect(get_components(PowerLoad, system))
    
    ts_array = Dict()
    ts_names = get_time_series_names(Deterministic, loads[1])
    
    for load in loads
        if !haskey(ts_array, get_bus_name(load))
            ts_array[get_bus_name(load)] = get_time_series_values(Deterministic, load, ts_names[1])
        else
            ts_array[get_bus_name(load)] = ts_array[get_bus_name(load)] .+ get_time_series_values(Deterministic, load, ts_names[1])
        end
    end
    
    ts_array = DataFrame(ts_array)
    
    # Renewable Data
    renewables = collect(get_components(RenewableDispatch, system))
    
    ts_renewable = Dict()
    ts_renewable_names = get_time_series_names(Deterministic, renewables[1])
    
    for renewable in renewables
        if !haskey(ts_renewable, get_bus_name(renewable))
            ts_renewable[get_bus_name(renewable)] = get_time_series_values(Deterministic, renewable, ts_renewable_names[1])
        else
            ts_renewable[get_bus_name(renewable)] = ts_renewable[get_bus_name(renewable)] .+ get_time_series_values(Deterministic, renewable, ts_renewable_names[1])
        end
    end
    
    ts_renewable = DataFrame(ts_renewable)
    
    bus_map = bus_mapping(system)
    
    # select rows for the given bus names, default to all buses.
    if !isempty(bus_names)
        bus_names = String.(bus_names)
        @assert issubset(bus_names, [names(ts_renewable); names(ts_array)])
        select!(ts_renewable, intersect(bus_names, names(ts_renewable)))
        select!(ts_array, intersect(bus_names, names(ts_array)))
    end
    
    all_data = sum.(eachrow(ts_array)) - sum.(eachrow(ts_renewable))
    times = get_time_series_timestamps(Deterministic, loads[1], ts_names[1])
    plot_data = DataFrame(net_demand = all_data) .* get_base_power(system)
    
    label --> reduce(hcat, names(plot_data))
    yguide --> "Net Demand (MWh)"
    legend --> :outertopright
    seriestype --> :line
    xrotation --> 45
    title --> "Net Demand over the hours"
    
    # now stack the matrix to get the cumulative values over all fuel types
    data = cumsum(Matrix(plot_data); dims=2)
    for i in Base.axes(data, 2)
        @series begin
            fillrange := i > 1 ? data[:, i - 1] : 0
            times, data[:, i]
        end
    end
end
