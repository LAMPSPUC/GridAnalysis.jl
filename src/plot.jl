"""
    plot_generation_stack(
        system::System,
        results::SimulationProblemResults;
        generator_fields::AbstractArray=[:P__ThermalStandard, :P__RenewableDispatch],
        bus_names::AbstractArray=[],
    )

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
        elseif occursin("PV", fuel_names[i]) == true ||
               occursin("CSP", fuel_names[i]) == true
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
    title --> "Generation stacked over the hours"

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
        bus_names::AbstractArray=[]
    )

Plot the simulation prices over the time period covered by the `results`. The `bus_names`
control which buses we want to include the plot. If the `market_simulator` is the one
that evaluate the prices on Real Time (RT), it is on \$/MW-15min.
"""
@userplot plot_prices
@recipe function f(p::plot_prices; bus_names::AbstractArray=[], type::String="ED")
    market_simulator, system_results, = p.args

    if isa(market_simulator, UCEDRT)
        prices = evaluate_prices_UCEDRT(market_simulator, system_results)

        if type == "DA"
            # Selects the plot data if is desired to plot the ED prices
            plot_data = select(prices["DA"], Not(:DateTime))

            yguide --> "Prices (\$/MWh)"

            times = prices["DA"][!, 1]
        else
            # Selects the plot data if is desired to plot the RT prices
            plot_data = select(prices["RT"], Not(:DateTime))

            sys_rt = market_simulator.system_rt
            resolution_mili_seg = get_time_series_resolution(sys_rt)
            resolution_seg = resolution_mili_seg/1000
            resolution_min = resolution_seg/60
            string_interval = string(resolution_min.value)
            yguide --> "Prices (\$/MW-" * string_interval * "min)"

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
        sys_rt = market_simulator.system_rt
        params = get_time_series_params(sys_rt)
        interval = params.interval
        string_interval = string(interval.value)
        yguide --> "Prices (\$/MW-" * string_interval * "min)"
    elseif isa(market_simulator, UCED)
        yguide --> "Prices (\$/MWh)"
    end

    label --> reduce(hcat, names(plot_data))
    legend --> :outertopright
    seriestype --> :line
    xrotation --> 45
    title --> "Prices over the hours"

    for i in Base.axes(plot_data, 2)
        @series begin
            times, plot_data[:, i]
        end
    end
end

"""
    plot_thermal_commit_generator_stack(
        system::System,
        results::SimulationProblemResults;
        bus_names::AbstractArray=[],
    )

Function to plot the Thermal Standard Commit variables over the time period covered by the `results`.
The `results` should be from the unit commitment problem.
It stacks the data so that it is possible to know how many generators are ON in each hour.
It groups by the generators.
"""
@userplot plot_thermal_commit_generator_stack
@recipe function f(p::plot_thermal_commit_generator_stack; bus_names::AbstractArray=[])
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
    title --> "Thermal Standard Commit stacked over the hours"

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
    plot_thermal_commit_type_stack(
        system::System,
        results::SimulationProblemResults;
        bus_names::AbstractArray=[],
    )

Function to plot the Thermal Standard Commit variables over the time period covered by the `results`.
The `results` should be from the unit commitment problem.
It stacks the data so that it is possible to know how many generators are ON in each hour.
It groups by the fuel type.
"""
@userplot plot_thermal_commit_type_stack
@recipe function f(p::plot_thermal_commit_type_stack; bus_names::AbstractArray=[])
    system, system_results, = p.args

    # get the output data for all fuel types
    variable_results = read_realized_variables(system_results)
    thermal_commit_results = variable_results[:On__ThermalStandard]

    plot_data = select(thermal_commit_results, Not(:DateTime))

    names_plot_data = names(plot_data)
    steam = zeros(24)
    NG = zeros(24)
    nuclear = zeros(24)
    for i in 1:length(names_plot_data)
        if occursin("STEAM", names_plot_data[i]) == true
            steam = hcat(steam, plot_data[!,i])
        elseif occursin("CT", names_plot_data[i]) == true || occursin("CC", names_plot_data[i]) == true
            NG = hcat(NG, plot_data[!,i])
        elseif occursin("NUCLEAR", names_plot_data[i]) == true
            nuclear = hcat(nuclear, plot_data[!,i])
        end
    end

    steam = vec(sum(steam, dims = 2))
    NG = vec(sum(NG, dims = 2))
    nuclear = vec(sum(nuclear, dims = 2))
    plot_data = DataFrame(Coal = steam, Natural_Gas = NG, Nuclear = nuclear)
    
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
    title --> "Thermal Standard Commit stacked over the hours"

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
    plot_demand_stack(
        system::System;
        bus_names::AbstractArray=[],
        type::String="SingleTimeSeries",
        start_time::Union{Nothing,Dates.DateTime}=nothing,
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
    start_time::Union{Nothing,Dates.DateTime}=nothing,
)
    system, = p.args

    # Getting the time series of the Demand
    loads = collect(get_components(PowerLoad, system))

    ts_array = Dict()
    if type == "SingleTimeSeries"
        ts_names = get_time_series_names(SingleTimeSeries, loads[1])

        for load in loads
            if !haskey(ts_array, get_bus_name(load))
                ts_array[get_bus_name(load)] = get_time_series_values(
                    SingleTimeSeries, load, ts_names[1]
                )
            else
                ts_array[get_bus_name(load)] =
                    ts_array[get_bus_name(load)] .+
                    get_time_series_values(SingleTimeSeries, load, ts_names[1])
            end
        end
    elseif type == "Deterministic"
        ts_names = get_time_series_names(Deterministic, loads[1])

        for load in loads
            if !haskey(ts_array, get_bus_name(load))
                ts_array[get_bus_name(load)] = get_time_series_values(
                    Deterministic, load, ts_names[1]; start_time
                )
            else
                ts_array[get_bus_name(load)] =
                    ts_array[get_bus_name(load)] .+
                    get_time_series_values(Deterministic, load, ts_names[1]; start_time)
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
        type::String="SingleTimeSeries",
        start_time::Union{Nothing,Dates.DateTime}=nothing,
    )

Function to plot the Demand over the time period covered by the `results`.
The `bus_names` controls which buses we want to include in the plot.
Renewable Dispatch data is the time series from the `system`.
The `type` controls if is wanted to plot the whole time series or just the first 24 observations.
The `start_time` controls in which day is going to be ploted for the Deterministic case.
"""
@userplot plot_net_demand_stack
@recipe function f(
    p::plot_net_demand_stack;
    bus_names::AbstractArray=[],
    type::String="SingleTimeSeries",
    start_time::Union{Nothing,Dates.DateTime}=nothing,
)
    system, = p.args

    # Getting the time series of the Demand
    loads = collect(get_components(PowerLoad, system))

    ts_array = Dict()
    if type == "SingleTimeSeries"
        ts_names = get_time_series_names(SingleTimeSeries, loads[1])

        for load in loads
            if !haskey(ts_array, get_bus_name(load))
                ts_array[get_bus_name(load)] = get_time_series_values(
                    SingleTimeSeries, load, ts_names[1]
                )
            else
                ts_array[get_bus_name(load)] =
                    ts_array[get_bus_name(load)] .+
                    get_time_series_values(SingleTimeSeries, load, ts_names[1])
            end
        end
    elseif type == "Deterministic"
        ts_names = get_time_series_names(Deterministic, loads[1])

        for load in loads
            if !haskey(ts_array, get_bus_name(load))
                ts_array[get_bus_name(load)] = get_time_series_values(
                    Deterministic, load, ts_names[1]; start_time
                )
            else
                ts_array[get_bus_name(load)] =
                    ts_array[get_bus_name(load)] .+
                    get_time_series_values(Deterministic, load, ts_names[1]; start_time)
            end
        end
    end

    ts_array = DataFrame(ts_array)

    # Getting the Renewable Data from the time series
    renewables = collect(get_components(RenewableGen, system))

    ts_renewable = Dict()
    if type == "SingleTimeSeries"
        ts_renewable_names = get_time_series_names(SingleTimeSeries, renewables[1])

        for renewable in renewables
            if !haskey(ts_renewable, get_bus_name(renewable))
                ts_renewable[get_bus_name(renewable)] = get_time_series_values(
                    SingleTimeSeries, renewable, ts_renewable_names[1]
                )
            else
                ts_renewable[get_bus_name(renewable)] =
                    ts_renewable[get_bus_name(renewable)] .+ get_time_series_values(
                        SingleTimeSeries, renewable, ts_renewable_names[1]
                    )
            end
        end
    elseif type == "Deterministic"
        ts_renewable_names = get_time_series_names(Deterministic, renewables[1])

        for renewable in renewables
            if !haskey(ts_renewable, get_bus_name(renewable))
                ts_renewable[get_bus_name(renewable)] = get_time_series_values(
                    Deterministic, renewable, ts_renewable_names[1]; start_time
                )
            else
                ts_renewable[get_bus_name(renewable)] =
                    ts_renewable[get_bus_name(renewable)] .+ get_time_series_values(
                        Deterministic, renewable, ts_renewable_names[1]; start_time
                    )
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

    plot_data = DataFrame(; net_demand=all_data) .* get_base_power(system)

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

"""
    plot_generation_stack_virtual(
        system::System,
        results::SimulationProblemResults;
        generator_fields::AbstractArray=[:P__ThermalStandard, :P__RenewableDispatch],
        period::Int=1,
        initial_time::Date,
    )

Plot the generation mix during the time `period` for the range of virtual bids in `results`. 
"""
@userplot plot_generation_stack_virtual
@recipe function f(
    p::plot_generation_stack_virtual;
    generator_fields::AbstractArray=[:P__ThermalStandard, :P__RenewableDispatch],
    type::AbstractString,
    period::Int=1,
    initial_time::Date,
    bus_names::AbstractArray=[],
)
    system, results_df, = p.args
    results_df = sort(results_df)
    aux_period = DateTime(initial_time) + Hour(period[1] - 1)

    for q in keys(results_df)
        if type == "RT"
            system_results = get_problem_results(results_df[q][type], "RT")
        elseif type == "ED" || type == "DA"
            system_results = get_problem_results(results_df[q][type], "UC")
        else
            system_results = get_problem_results(results_df[q], "UC")
        end

        # get mapping from busname to fuel type
        fuel_type_dict = fuel_type_mapping(system)

        # get the output data for given fuel types
        variable_results = read_realized_variables(system_results; names=generator_fields)
        generator_data = getindex.(Ref(variable_results), generator_fields)
        for i in 1:length(generator_data)
            generation = generator_data[i][aux_period .<= generator_data[i].DateTime .< aux_period + Hour(1), :]
            generation[!, "DateTime"] .= aux_period
            generator_data[i] = combine(
                groupby(generation, :DateTime),
                names(generation, Not(:DateTime)) .=> sum;
                renamecols=false,
            )
            lin, col = size(generator_data[i])
            for j in 2:col
                generator_data[i][1, j] =
                    generator_data[i][1, j] / length(generation[!, "DateTime"])
            end
        end
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
            elseif occursin("PV", fuel_names[i]) == true ||
                occursin("CSP", fuel_names[i]) == true
                fuel_type_dict[key] = "SOLAR"
            end
        end

        # select rows for the given bus names, default to all buses.
        if !isempty(bus_names)
            bus_names = String.(bus_names)
            #@assert all(bus_names .∈ [stacked_data.bus_name])
            @assert issubset(bus_names, stacked_data.bus_name)
            filter!(:bus_name => in(bus_names), stacked_data)
        end

        stacked_data.fuel_type = get.(Ref(fuel_type_dict), stacked_data.gen_name, missing)

        aggregated_data = combine(
            groupby(stacked_data, [:DateTime, :fuel_type]), :output => sum => :output
        )

        # convert the output units into MWh.
        aggregated_data.output = aggregated_data.output .* get_base_power(system)

        # unstack aggregated data and make area plot 
        global unstacked_data = unstack(aggregated_data, :fuel_type, :output)
        if q == minimum(keys(results_df))
            global data_frame = unstacked_data
        else
            global data_frame = append!(data_frame, unstacked_data)
        end
    end

    if (type == "ED" || type == "DA") && "OTHER" in names(data_frame)
        aux = data_frame[!, "OTHER"]
        global data_frame = select(data_frame, Not(:OTHER))
        global data_frame[!, "OTHER"] = aux
    end

    times = collect(keys(results_df)) * get_base_power(system)
    palette = :Dark2_8

    plot_data = select(data_frame, Not(:DateTime))

    label --> reduce(hcat, names(plot_data))
    yguide --> "Output (MWh)"
    xguide --> "Quantity (MW)"
    legend --> :outertopright
    seriestype --> :line
    xrotation --> 0
    color_palette --> palette
    title --> "Virtual Generation stacked over the hours"

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
    plot_thermal_commit_virtual(
        system::System,
        results::SimulationProblemResults;
        period::Int=1,
        initial_time::Date,
    )

Function to plot the Thermal Standard Commit variables curve depending on the range_quota.
The `results` should be from the unit commitment problem.
1 is ON, 0 is OFF.
"""
@userplot plot_thermal_commit_virtual
@recipe function f(
    p::plot_thermal_commit_virtual; 
    period::Int=1,
    initial_time::Date,
    bus_names::AbstractArray=[]
)

    system, results_df, = p.args
    results_df = sort(results_df)
    aux_period = DateTime(initial_time) + Hour(period[1] - 1)

    for v in keys(results_df)
        # get the output data for all fuel types
        system_results = get_problem_results(results_df[v]["DA"], "UC")
        variable_results = read_realized_variables(system_results)
        commit_results = variable_results[:On__ThermalStandard]

        commit_results = commit_results[aux_period .<= commit_results.DateTime .< aux_period + Hour(1), :]
        commit_results[!, "DateTime"] .= aux_period

        if v == minimum(keys(results_df))
            global thermal_commit_results = commit_results
        else
            global thermal_commit_results = append!(thermal_commit_results, commit_results)
        end
    end

    plot_data = select(thermal_commit_results, Not(:DateTime))
    
    #=
    names_plot_data = names(plot_data)
    steam = zeros(length(results_df))
    NG = zeros(length(results_df))
    nuclear = zeros(length(results_df))
    for i in 1:length(names_plot_data)
        if occursin("STEAM", names_plot_data[i]) == true
            steam = hcat(steam, plot_data[!,i])
        elseif occursin("CT", names_plot_data[i]) == true || occursin("CC", names_plot_data[i]) == true
            NG = hcat(NG, plot_data[!,i])
        elseif occursin("NUCLEAR", names_plot_data[i]) == true
            nuclear = hcat(nuclear, plot_data[!,i])
        end
    end

    steam = vec(sum(steam, dims = 2))
    NG = vec(sum(NG, dims = 2))
    nuclear = vec(sum(nuclear, dims = 2))
    plot_data = DataFrame(Coal = steam, Natural_Gas = NG, Nuclear = nuclear)
    =#

    # select rows for the given bus names, default to all buses.
    if !isempty(bus_names)
        generator_names = names(plot_data)
        bus_map = bus_mapping(system)
        bus_names = String.(bus_names)
        @assert issubset(bus_names, [bus_map[gen] for gen in generator_names])
        generator_names = [gen for gen in generator_names if in(bus_map[gen], bus_names)]
        select!(plot_data, generator_names)
    end

    times = collect(keys(results_df)) * get_base_power(system)

    label --> reduce(hcat, names(plot_data))
    yguide --> "Commitment"
    legend --> :outertopright
    seriestype --> :line
    xrotation --> 45
    title --> "Thermal Standard Commit over the hours"

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
    plot_flow_capacity_stack(
        system::System;
        bus_names::AbstractArray=[],
        type::String="SingleTimeSeries",
        start_time::Union{Nothing,Dates.DateTime}=nothing,
    )

Function to plot flow_capacity/demand over the time period covered by the `results`.
The `bus_names` controls which buses we want to include in the plot.
The `type` controls if is wanted to plot the whole time series or just the first 24 observations.
The `start_time` controls in which day is going to be ploted for the Deterministic case.
"""
@userplot plot_flow_capacity_stack
@recipe function f(
    p::plot_flow_capacity_stack;
    bus_names::AbstractArray=[],
    type::String="SingleTimeSeries",
    start_time::Union{Nothing,Dates.DateTime}=nothing,
)
    system, = p.args

    # Getting the time series of the Demand
    loads = collect(get_components(PowerLoad, system))

    ts_array = Dict()
    if type == "SingleTimeSeries"
        ts_names = get_time_series_names(SingleTimeSeries, loads[1])

        for load in loads
            if !haskey(ts_array, get_bus_name(load))
                ts_array[get_bus_name(load)] = get_time_series_values(
                    SingleTimeSeries, load, ts_names[1]
                )
            else
                ts_array[get_bus_name(load)] =
                    ts_array[get_bus_name(load)] .+
                    get_time_series_values(SingleTimeSeries, load, ts_names[1])
            end
        end
    elseif type == "Deterministic"
        ts_names = get_time_series_names(Deterministic, loads[1])

        for load in loads

            if !haskey(ts_array, get_bus_name(load))
                ts_array[get_bus_name(load)] = get_time_series_values(
                    Deterministic, load, ts_names[1]; start_time
                )
            else
                ts_array[get_bus_name(load)] =
                    ts_array[get_bus_name(load)] .+
                    get_time_series_values(Deterministic, load, ts_names[1]; start_time)
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

    for bus in names(plot_data)
        r=[get_rate(i) for i in get_components(Branch, system) if 
        (get_name(get_area(get_from(get_arc(i)))) == get_name(get_area(get_to(get_arc(i)))) &&
          (get_name(get_to(get_arc(i))) == bus || get_name(get_from(get_arc(i))) == bus))].* get_base_power(system)
    
        if r!=[]
            plot_data[!,bus]=sum(r) ./ plot_data[!,bus]
        else
            plot_data[!,bus].=0
        end
    end

    label --> reduce(hcat, names(plot_data))
    yguide --> "FlowCapacity/demand"
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
plot_prev_renew_stack(
        system::System;
        bus_names::AbstractArray=[],
        type::String="SingleTimeSeries",
        start_time::Union{Nothing,Dates.DateTime}=nothing,
    )

Function to plot the renewable generation forecast over the time period covered by the `results`.
The `bus_names` controls which buses we want to include in the plot.
Renewable Dispatch data is the time series from the `system`.
The `type` controls if is wanted to plot the whole time series or just the first 24 observations.
The `start_time` controls in which day is going to be ploted for the Deterministic case.
"""
@userplot plot_prev_renew_stack
@recipe function f(
    p::plot_prev_renew_stack;
    bus_names::AbstractArray=[],
    type::String="SingleTimeSeries",
    start_time::Union{Nothing,Dates.DateTime}=nothing,
)
    system, = p.args

    # Getting the Renewable Data from the time series
    renewables = collect(get_components(RenewableGen, system))
    loads = collect(get_components(PowerLoad, system))

    ts_renewable = Dict()
    if type == "SingleTimeSeries"
        ts_renewable_names = get_time_series_names(SingleTimeSeries, renewables[1])

        for renewable in renewables
            if !haskey(ts_renewable, get_bus_name(renewable))
                ts_renewable[get_bus_name(renewable)] = get_time_series_values(
                    SingleTimeSeries, renewable, ts_renewable_names[1]
                )
            else
                ts_renewable[get_bus_name(renewable)] =
                    ts_renewable[get_bus_name(renewable)] .+ get_time_series_values(
                        SingleTimeSeries, renewable, ts_renewable_names[1]
                    )
            end
        end
    elseif type == "Deterministic"
        ts_renewable_names = get_time_series_names(Deterministic, renewables[1])

        for renewable in renewables
            if !haskey(ts_renewable, get_bus_name(renewable))
                ts_renewable[get_bus_name(renewable)] = get_time_series_values(
                    Deterministic, renewable, ts_renewable_names[1]; start_time
                )
            else
                ts_renewable[get_bus_name(renewable)] =
                    ts_renewable[get_bus_name(renewable)] .+ get_time_series_values(
                        Deterministic, renewable, ts_renewable_names[1]; start_time
                    )
            end
        end
    end

    ts_renewable = DataFrame(ts_renewable)

    bus_map = bus_mapping(system)

    # select rows for the given bus names, default to all buses.
    if !isempty(bus_names)
        bus_names = String.(bus_names)
        @assert issubset(bus_names, names(ts_renewable))
        select!(ts_renewable, intersect(bus_names, names(ts_renewable)))
    end

    all_data = ts_renewable
    
    ts_names = get_time_series_names(SingleTimeSeries, loads[1])

    if type == "SingleTimeSeries"
        times = get_time_series_timestamps(SingleTimeSeries, loads[1], ts_names[1])
    elseif type == "Deterministic"
        times = get_time_series_timestamps(Deterministic, loads[1], ts_names[1]; start_time)
    end

    plot_data = all_data .* get_base_power(system)

    label --> reduce(hcat, names(plot_data))
    yguide --> "Generation (MWh)"
    legend --> :outertopright
    seriestype --> :line
    xrotation --> 45
    title --> "Generation over the hours"

    # now stack the matrix to get the cumulative values over all fuel types
    data = cumsum(Matrix(plot_data); dims=2)
    for i in Base.axes(data, 2)
        @series begin
            fillrange := i > 1 ? data[:, i - 1] : 0
            times, data[:, i]
        end
    end
end
 
#= TODO: curves plots with recipe function
@userplot plot_offert_price()
@recipe function f(p::plot_offert_price; period::Vector{Int64}, bus_name::AbstractArray=["bus5"])
    lmps_df, = p.args

    lmps_df = sort(lmps_df)

    index=[]
    data=Array{Any}(nothing, (length(period),length(bus_name)+1,length(lmps_df))) #length(period)-size(lmps_df[0])[1]
    for (i, v) in enumerate(keys(lmps_df))
        for t in period #arrumar
            data[t,1,i]=lmps_df[v][!,"DateTime"][t]
            c=2
            for j in bus_name
                data[t,c,i]=lmps_df[v][!,j][t]
                c=c+1
            end
        end
        index = vcat(index, v)
    end #TODO: select period better

    plot_data = data #data[t,bus,max_gen]: selecionado a linha que quero
    #label --> reduce(hcat, names(plot_data))
    yguide --> "Prices (\$/MWh)"
    legend --> :outertopright
    seriestype --> :line
    xrotation --> 45

    for i in bus_name
        for t in Base.axes(plot_data, 1) # verificar
            @series begin
                index, plot_data[t,:,i] #pegar a coluna do bus_name
            end
        end
    end
end

@userplot plot_revenue()
@recipe function f(p::plot_revenue; g::, period::Vector{Int64}, generator_name::AbstractArray)
    lmps_df,=p.args
    lmps_df = sort(lmps_df)
    results_df=g
    gen=get_component(ThermalStandard, market_simulator.system_uc, generator_name)
    bus_name=get_bus(gen)

    #ED price (bus_name)
    #RT price (bus_name)
    #UC generation (bus_name)

    #Revenue = (price["ED"] - price["RT"])*generation("UC")
    #Revenue = (price["RT"] - price["ED"])*generation("UC")

    #Revenue = price*generation("UC")

    index=[]
    data=Array{Any}(nothing, (length(period),2,length(lmps_df)))
    for (i, v) in enumerate(keys(lmps_df))
        for t in period
            data[t,1,i]=lmps_df[v][!,"DateTime"][t]
            variable_results = read_realized_variables(get_problem_results(results_df[v], "UC"), names=[:P__ThermalStandard])
            generator_data = getindex.(Ref(variable_results), [:P__ThermalStandard])
            virtual_gen=generator_data[1][!,7][t] #name_generator
            data[t,2,i]=(lmps_df[v][t,bus_name][1])*virtual_gen #arrumar 
        end
        index = vcat(index, v)
    end#TODO: select only period

    plot_data = data #data[t,bus,max_gen]: selecionado a linha que quero
    #label --> reduce(hcat, names(plot_data))
    yguide --> "Revenues (\$/MWh)"
    legend --> :outertopright
    seriestype --> :line
    xrotation --> 45

    for t in Base.axes(plot_data, 2) # verificar
        @series begin
            index, plot_data[t,2,:]
        end
    end

end
=#
