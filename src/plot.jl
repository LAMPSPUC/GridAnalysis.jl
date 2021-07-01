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

    fuel_names = unique(keys(fuel_type_dict))
    for (i, key) in enumerate(fuel_names)
        if occursin("HYDRO", fuel_names[i]) == true
            fuel_type_dict[key] = "HYDRO"
        elseif occursin("WIND", fuel_names[i]) == true
            fuel_type_dict[key] = "WIND"
        elseif occursin("PV", fuel_names[i]) == true || occursin("CSP", fuel_names[i]) == true
            fuel_type_dict[key] = "SOLAR"
        end
    end

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
        results::SimulationResults;
        bus_names::AbstractArray=[])

Plot the simulation prices over the time period covered by the `results`. The `bus_names`
control which buses we want to include the plot. If the `market_simulator` is the one
that evaluate the prices on Real Time (RT), it is on \$/MW-15min.
"""
@userplot plot_prices
@recipe function f(p::plot_prices; bus_names::AbstractArray=[], type::String="ED")
    market_simulator, system_results, = p.args

    if isa(market_simulator, UCEDRT)
        prices = evaluate_prices_UCEDRT(market_simulator, system_results)
        
        if type == "ED"
            # Selects the plot data if is desired to plot the ED prices
            plot_data = select(prices["ED"], Not(:DateTime))

            yguide --> "Prices (\$/MWh)"

            times = prices["ED"][!, 1]
        else
            # Selects the plot data if is desired to plot the ED prices
            plot_data = select(prices["RT"], Not(:DateTime))

            yguide --> "Prices (\$/MW-15min)"

            times = prices["RT"][!, 1]
        end
    else
        prices = evaluate_prices(market_simulator, system_results)

        prices_keys = collect(keys(prices))
    
        plot_data = select(prices[prices_keys[1]], Not(:DateTime))
    
        times = prices[prices_keys[1]][!, 1]
    end

    # select rows for the given bus names, default to all buses.
    if !isempty(bus_names)
        bus_names = String.(bus_names)
        @assert issubset(bus_names, names(plot_data))
        select!(plot_data, bus_names)
    end

    if isa(market_simulator, UCRT)
        yguide --> "Prices (\$/MW-15min)"
    elseif isa(market_simulator, UCED)
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
    plot_thermal_commit(
        system::System,
        results::SimulationProblemResults;
        bus_names::AbstractArray=[],
    )

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

# TO-DO: select the right number of days simulated.
# Currently is selecting all the time series and not the number of days simulated.
"""
    plot_demand_stack(
        system::System;
        bus_names::AbstractArray=[]
    )

Function to plot the Demand over the time period covered by the `results`.
The `bus_names` controls which buses we want to include in the plot.
The `type` controls if is wanted to plot the whole time series or just the first 24 observations.
The `start_time` controls in which day is going to be ploted for the Deterministic case.
"""
@userplot plot_demand_stack
@recipe function f(
    p::plot_demand_stack;
    bus_names::AbstractArray=[],
    type::String="SingleTimeSeries",
    start_time::Union{Nothing, Dates.DateTime} = nothing,
)
    system, = p.args

    # Getting the time series of the Demand
    loads = collect(get_components(PowerLoad, system))

    ts_array = Dict()
    if type == "SingleTimeSeries"
        ts_names = get_time_series_names(SingleTimeSeries, loads[1])
        
        for load in loads
            if !haskey(ts_array, get_bus_name(load))
                ts_array[get_bus_name(load)] = get_time_series_values(SingleTimeSeries, load, ts_names[1])
            else
                ts_array[get_bus_name(load)] = ts_array[get_bus_name(load)] .+ get_time_series_values(SingleTimeSeries, load, ts_names[1])
            end
        end
    elseif type == "Deterministic"
        ts_names = get_time_series_names(Deterministic, loads[1])
        
        for load in loads
            if !haskey(ts_array, get_bus_name(load))
                ts_array[get_bus_name(load)] = get_time_series_values(Deterministic, load, ts_names[1]; start_time)
            else
                ts_array[get_bus_name(load)] = ts_array[get_bus_name(load)] .+ get_time_series_values(Deterministic, load, ts_names[1]; start_time)
            end
        end
    end

    ts_array = DataFrame(ts_array)

    if type == "SingleTimeSeries"
        times = get_time_series_timestamps(SingleTimeSeries, loads[1], ts_names[1])
    elseif type == "Deterministic"
        times = get_time_series_timestamps(Deterministic, loads[1], ts_names[1]; start_time)
    end
    
    plot_data = ts_array .* get_base_power(system)

    # select rows for the given bus names, default to all buses.
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

# TO-DO: select the right number of days simulated.
# Currently is selecting all the time series and not the number of days simulated.
"""
    plot_net_demand_stack(
        system::System;
        bus_names::AbstractArray=[],
    )

Function to plot the Demand over the time period covered by the `results`.
The `bus_names` controls which buses we want to include in the plot.
Renewable Dispatch data is the time series from the `system`.
The `type` controls if is wanted to plot the whole time series or just the first 24 observations.
The `start_time` controls in which day is going to be ploted for the Deterministic case.
"""
@userplot plot_net_demand_stack
@recipe function f(p::plot_net_demand_stack;
    bus_names::AbstractArray=[],
    type::String="SingleTimeSeries",
    start_time::Union{Nothing, Dates.DateTime} = nothing,
)
    
    system, = p.args
    
    # Getting the time series of the Demand
    loads = collect(get_components(PowerLoad, system))
    
    ts_array = Dict()
    if type == "SingleTimeSeries"
        ts_names = get_time_series_names(SingleTimeSeries, loads[1])
        
        for load in loads
            if !haskey(ts_array, get_bus_name(load))
                ts_array[get_bus_name(load)] = get_time_series_values(SingleTimeSeries, load, ts_names[1])
            else
                ts_array[get_bus_name(load)] = ts_array[get_bus_name(load)] .+ get_time_series_values(SingleTimeSeries, load, ts_names[1])
            end
        end
    elseif type == "Deterministic"
        ts_names = get_time_series_names(Deterministic, loads[1])
        
        for load in loads
            if !haskey(ts_array, get_bus_name(load))
                ts_array[get_bus_name(load)] = get_time_series_values(Deterministic, load, ts_names[1]; start_time)
            else
                ts_array[get_bus_name(load)] = ts_array[get_bus_name(load)] .+ get_time_series_values(Deterministic, load, ts_names[1]; start_time)
            end
        end
    end
    
    ts_array = DataFrame(ts_array)
    
    # Getting the Renewable Data from the time series
    renewables = collect(get_components(RenewableDispatch, system))
    
    ts_renewable = Dict()
    if type == "SingleTimeSeries"
        ts_renewable_names = get_time_series_names(SingleTimeSeries, renewables[1])
        
        for renewable in renewables
            if !haskey(ts_renewable, get_bus_name(renewable))
                ts_renewable[get_bus_name(renewable)] = get_time_series_values(SingleTimeSeries, renewable, ts_renewable_names[1])
            else
                ts_renewable[get_bus_name(renewable)] = ts_renewable[get_bus_name(renewable)] .+ get_time_series_values(SingleTimeSeries, renewable, ts_renewable_names[1])
            end
        end
    elseif type == "Deterministic"
        ts_renewable_names = get_time_series_names(Deterministic, renewables[1])
        
        for renewable in renewables
            if !haskey(ts_renewable, get_bus_name(renewable))
                ts_renewable[get_bus_name(renewable)] = get_time_series_values(Deterministic, renewable, ts_renewable_names[1]; start_time)
            else
                ts_renewable[get_bus_name(renewable)] = ts_renewable[get_bus_name(renewable)] .+ get_time_series_values(Deterministic, renewable, ts_renewable_names[1]; start_time)
            end
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
    
    # Evaluating the Net Demand (Demand-Renewable)
    all_data = sum.(eachrow(ts_array)) - sum.(eachrow(ts_renewable))

    if type == "SingleTimeSeries"
        times = get_time_series_timestamps(SingleTimeSeries, loads[1], ts_names[1])
    elseif type == "Deterministic"
        times = get_time_series_timestamps(Deterministic, loads[1], ts_names[1]; start_time)
    end
    
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
