"""
    fuel_type_mapping(system)

Return a mapping between the bus name and the fuel type for the given `system`.
"""
function fuel_type_mapping(system::System)
    generator_metadata = [gen for gen in get_components(Generator, system)]

    fuel_type_map = Dict()
    for generator in generator_metadata
        name = generator.name
        try
            fuel_type_map[name] = get_fuel(generator)
        catch
            # Assumes the bus name has format "<Fuel>Bus..."
            fuel_type_map[name] = first(split(name, "Bus"))
        end
    end

    return fuel_type_map
end

"""
    get_bus_name(gen::Generator)

Returns the bus name of a generator or load.
"""
function get_bus_name(gen::Union{Generator,PowerLoad})
    return get_name(get_bus(gen))
end

"""
    bus_mapping(system::System)

Returns generator to bus mapping.
"""
function bus_mapping(system::System)
    generator_metadata = [gen for gen in get_components(Generator, system)]

    bus_map = Dict()
    for generator in generator_metadata
        name = generator.name
        bus_map[name] = get_bus_name(generator)
    end

    return bus_map
end

"""
    duals_constraint_names(<:AbstractPowerModel)

Return the constraints for which we care about the duals (because they form the energy prices) when using a specified network formulation.
"""
duals_constraint_names(::Type{CopperPlatePowerModel}) = [:CopperPlateBalance]
function duals_constraint_names(::Union{Type{NFAPowerModel},Type{DCPPowerModel}})
    return [:nodal_balance_active__Bus]
end
function duals_constraint_names(::Type{StandardPTDFModel})
    return [:CopperPlateBalance, :network_flow__Line]
end
function duals_constraint_names(market_simulator::UCED)
    return duals_constraint_names(market_simulator.template_ed.transmission)
end
function duals_constraint_names(market_simulator::UCRT)
    return duals_constraint_names(market_simulator.template_rt.transmission)
end
function duals_constraint_names(market_simulator::UCEDRT)
    d1 = duals_constraint_names(market_simulator.template_ed.transmission)
    d2 = duals_constraint_names(market_simulator.template_rt.transmission)
    return [d1, d2]
end

"""
    MEC(system::System, problem_results::PSI.SimulationProblemResults)

Returns the Marginal Energy Component (MEC) of OPF problem's solution. This component is defined as the dual value
(or lagrangian multiplier) of a grid-wide energy balance constraint (named in PSI `CopperPlateBalance`)
"""
function MEC(system::System, problem_results::PSI.SimulationProblemResults)
    duals = read_dual(problem_results, :CopperPlateBalance)
    return DataFrame(;
        DateTime=collect(keys(duals)),
        lmp=[duals[i][1, 2] / get_base_power(system) for i in keys(duals)],
    )
end

"""
    evaluate_prices(::Type{CopperPlatePowerModel}, system::System, problem_results::PSI.SimulationProblemResults, kwargs)

Returns a grid-unique series of energy prices for the simulation's data-range. 
In this formulation, for each period of the problem, the energy prices are constructed using the dual value
(or lagrangian multiplier) of a grid-wide energy balance constraint (named in PSI `CopperPlateBalance`). 
"""
function evaluate_prices(
    ::Type{CopperPlatePowerModel},
    system::System,
    problem_results::PSI.SimulationProblemResults,
    kwargs,
)
    return MEC(system, problem_results)
end

"""
    evaluate_prices(::Union{Type{NFAPowerModel},Type{DCPPowerModel}}, system::System, problem_results::PSI.SimulationProblemResults, kwargs)

Returns a nodal-wide series of energy prices (locational marginal prices - LMPS) for the simulation's data-range. 
In this formulation, for each period of the problem, the energy prices are constructed using the dual values
(or lagrangian multipliers) of each nodal energy balance constraint (named in PSI `nodal_balance_active__Bus`). 
"""
function evaluate_prices(
    ::Union{Type{NFAPowerModel},Type{DCPPowerModel}},
    system::System,
    problem_results::PSI.SimulationProblemResults,
    kwargs,
)
    duals = read_dual(problem_results, :nodal_balance_active__Bus)
    # get the interval
    params = get_time_series_params(system)
    interval = params.interval
    if typeof(interval) == Hour
        int_interval = interval.value*60
    else
        int_interval = Int(interval.value)
    end
    # get the interval resolution
    resolution_mili_seg = get_time_series_resolution(system)
    resolution_seg = resolution_mili_seg/1000
    resolution_min = resolution_seg/60
    int_resolution = Int(resolution_min.value)
    # makes interval/resolution to know how many prices are needed
    tot = Int(int_interval/int_resolution)
    # gets the needed number of prices hours, depending on the interval 
    # and on the resolution of the problem
    df = vcat([DataFrame(duals[i][1:tot, :]) for i in keys(duals)]...)
    # in this case, the lmps are the negative of the duals
    df[:, 2:end] = df[:, 2:end] ./ -get_base_power(system)
    return df
end

"""
    evaluate_prices(transmission::Type{StandardPTDFModel}, system::System, problem_results::PSI.SimulationProblemResults, kwargs)

Returns a nodal-wide series of energy prices (locational marginal prices - LMPS) for the simulation's data-range. 
In this formulation, for each period of the problem, the energy prices are constructed using the dual value
(or lagrangian multiplier) of a grid-wide energy balance constraint (named in PSI `CopperPlateBalance`) and 
the dual values of a branch-wide power-flow constraint (named in PSI `network_flow__Line`).
"""
function evaluate_prices(
    ::Type{StandardPTDFModel},
    system::System,
    problem_results::PSI.SimulationProblemResults,
    kwargs,
)
    # MEC
    mec = MEC(system, problem_results)

    # MCC
    PTDF_matrix = kwargs[:PTDF]
    # get line limit duals
    # ignoring transformer, because it seems that PSI is also ignoring
    flow_duals = read_dual(problem_results, :network_flow__Line)
    # get the interval
    params = get_time_series_params(system)
    interval = params.interval
    if typeof(interval) == Hour
        int_interval = interval.value*60
    else
        int_interval = Int(interval.value)
    end
    # get the interval resolution
    resolution_mili_seg = get_time_series_resolution(system)
    resolution_seg = resolution_mili_seg/1000
    resolution_min = resolution_seg/60
    int_resolution = Int(resolution_min.value)
    # makes interval/resolution to know how many prices are needed
    tot = int_interval/int_resolution
    # gets the needed number of prices hours, depending on the interval 
    # and on the resolution of the problem
    flow_duals = vcat([DataFrame(flow_duals[i][1:tot, :]) for i in keys(flow_duals)]...)
    # get duals in a matrix form ordered by the line names available in the PTDF
    line_names = intersect(PTDF_matrix.axes[1], names(flow_duals[:, 2:end]))
    μ = Matrix(flow_duals[:, line_names])
    # calculate the congestion factor of the prices by multiplying the line duals by the PTDF
    buses = get_components(Bus, sys_ed)
    congestion_lmp = Dict()
    for bus in buses
        congestion_lmp[get_name(bus)] =
            μ * [PTDF_matrix[line, get_number(bus)] for line in line_names]
    end
    congestion_lmp["DateTime"] = collect(keys(flow_duals))
    mcc = DataFrame(congestion_lmp)

    # Calculate the final locational marginal prices: LMP (MEC + MCC)
    LMP = deepcopy(mcc)
    for row in eachrow(LMP)
        for name in names(row[2:end])
            row[name] /= get_base_power(base_system)
            row[name] += .+mec[mec[:, "DateTime"] .== row["DateTime"], :mec][1]
        end
    end
    return LMP
end

"""
    evaluate_prices(market_simulator::UCED, problem_results::Dict{String, SimulationResults})

Returns energy prices for the simulation's data-range.  
"""
function evaluate_prices(
    market_simulator::UCED, problem_results::Dict{String,SimulationResults}
)
    ed_results = get_problem_results(problem_results["DA"], "ED")

    return Dict(
        "DA" => evaluate_prices(
            market_simulator.template_ed.transmission,
            market_simulator.system_ed,
            ed_results,
            market_simulator.kwargs,
        ),
    )
end

"""
    evaluate_prices(market_simulator::UCRT, problem_results::PSI.SimulationResults)

Returns energy prices for the simulation's data-range.  
"""
function evaluate_prices(
    market_simulator::UCRT, problem_results::Dict{String,SimulationResults}
)
    rt_results = get_problem_results(problem_results["RT"], "RT")

    return Dict(
        "RT" => evaluate_prices(
            market_simulator.template_rt.transmission,
            market_simulator.system_rt,
            rt_results,
            market_simulator.kwargs,
        ),
    )
end

"""
    evaluate_prices(market_simulator::UCEDRT, problem_results::Dict{String, SimulationResults})

Returns energy prices for the simulation's data-range.  
"""
function evaluate_prices_UCEDRT(
    market_simulator::UCEDRT, problem_results::Dict{String,SimulationResults}
)
    ed_results = get_problem_results(problem_results["DA"], "ED")
    rt_results = get_problem_results(problem_results["RT"], "RT")

    return Dict(
        "DA" => evaluate_prices(
            market_simulator.template_ed.transmission,
            market_simulator.system_ed,
            ed_results,
            market_simulator.kwargs,
        ),
        "RT" => evaluate_prices(
            market_simulator.template_rt.transmission,
            market_simulator.system_rt,
            rt_results,
            market_simulator.kwargs,
        ),
    )
end

"""
    get_time_series_params(system::System)

Returns the parameters associated with the time-series attached to the system.
"""
get_time_series_params(system::System) = system.data.time_series_params.forecast_params

"""
    load_pq_curves(market_simulator::UCEDRT, range_quota::Vector{Float64}, simulation_folder::String=pwd())

Returns the results of a simulation done previously.
"""
function load_pq_curves(
    market_simulator::UCEDRT, range_quota::Vector{Float64}, simulation_folder::String=pwd()
)
    lmps_df = Dict()
    results_df = Dict()
    for max_gen in range_quota
        results_df[max_gen] = Dict(
            "DA" =>
                SimulationResults(joinpath(simulation_folder, "da_quota_$max_gen", "1")),
            "RT" =>
                SimulationResults(joinpath(simulation_folder, "rt_quota_$max_gen", "1")),
        )
        lmps_df[max_gen] = evaluate_prices_UCEDRT(market_simulator, results_df[max_gen])
    end
    return lmps_df, results_df
end

"""
    load_mix_pq_curves(market_simulator::UCEDRT, range_quota_load::Vector{Float64}, range_quota_gen::Vector{Float64}, simulation_folder::String=pwd())

Returns the results of a simulation done previously.
"""
function load_mix_pq_curves(
    market_simulator::UCEDRT, range_quota_load::Vector{Float64}, range_quota_gen::Vector{Float64}, simulation_folder::String=pwd()
)
    lmps_df = Dict()
    results_df = Dict()
    for max_load in range_quota_load
        for max_gen in range_quota_gen
            results_df[[max_load max_gen]] = Dict(
                "DA" =>
                    SimulationResults(joinpath(simulation_folder, "da_quota_l_$max_load"*"_quota_g_$max_gen", "1")),
                "RT" =>
                    SimulationResults(joinpath(simulation_folder, "rt_quota_l_$max_load"*"_quota_g_$max_gen", "1")),
            )
            lmps_df[[max_load max_gen]] = evaluate_prices_UCEDRT(market_simulator, results_df[[max_load max_gen]])
        end
    end
    return lmps_df, results_df
end


"""
    load_pq_curves(market_simulator::UCED, range_quota::Vector{Float64}, simulation_folder::String=pwd())

Returns the results of a simulation done previously.
"""
function load_pq_curves(
    market_simulator::UCED, range_quota::Vector{Float64}, simulation_folder::String=pwd()
)
    lmps_df = Dict()
    results_df = Dict()
    for max_gen in range_quota
        results_df[max_gen] = Dict(
            "DA" => SimulationResults(joinpath(simulation_folder, "da_quota_$max_gen", "1"))
        )
        lmps_df[max_gen] = evaluate_prices(market_simulator, results_df[max_gen])
    end
    return lmps_df, results_df
end

"""
    load_pq_curves(market_simulator::UCRT, range_quota::Vector{Float64}, simulation_folder::String=pwd())

Returns the results of a simulation done previously.
"""
function load_pq_curves(
    market_simulator::UCRT, range_quota::Vector{Float64}, simulation_folder::String=pwd()
)
    lmps_df = Dict()
    results_df = Dict()
    for max_gen in range_quota
        results_df[max_gen] = Dict(
            "RT" => SimulationResults(
                joinpath(simulation_folder, "rt_quota_$max_gen", "1"), #It is assumed that the simulation is done just once
            )
        )
        lmps_df[max_gen] = evaluate_prices(market_simulator, results_df[max_gen])
    end
    return lmps_df, results_df
end

"""
    plot_price_curves(
        lmps_df::Dict{Any, Any}, period::Vector{Int64}, 
        bus_name::AbstractArray=["bus5"],
        node::String="bus5",
        initial_time::Date,
        system::System,
        ylimit::Bool,
    )

Function to plot the price curve for the virtual offer bids. 
The 'bus_names' and 'periods' controls which buses and periods we want to include
in the plot, respectively.
"""

function plot_price_curves(
    lmps_df::Dict{Any,Any},
    period::Vector{Int64},
    bus_name::AbstractArray,
    node::String,
    initial_time::Date,
    system::System,
    ylimit::Bool,
)
    lmps_df = sort(lmps_df)
    bus="lmp"
    aux_period = []
    for t in period
        aux_period = vcat(aux_period, DateTime(initial_time) + Hour(t - 1))
    end
    index = []
    max_element = 0
    min_element = 0
    data = Array{Any}(
        nothing,
        (
            length(period),
            length(bus_name) + 1,
            length(lmps_df),
            length(lmps_df[collect(keys(lmps_df))[1]]),
        ),
    )
    for (i, v) in enumerate(keys(lmps_df))
        for (l, k) in enumerate(keys(lmps_df[collect(keys(lmps_df))[1]]))
            for t in 1:length(period)
                data[t, 1, i, l] = aux_period[t]
                c = 2
                for j in bus_name
                    try
                        global prices_hour = lmps_df[v][k][
                            aux_period[t] .<= lmps_df[v][k].DateTime .< aux_period[t] + Hour(1),
                            j,
                        ]
                    catch
                        global prices_hour = lmps_df[v][k][
                            aux_period[t] .<= lmps_df[v][k].DateTime .< aux_period[t] + Hour(1),
                            bus,
                        ]
                    end

                    data[t, c, i, l] = sum(prices_hour; dims=1)[1]

                    if data[t, c, i, l] > max_element && data[t, c, i, l] < 1e3
                        max_element = data[t, c, i, l]
                    elseif data[t, c, i, l] < min_element && data[t, c, i, l] > -1e3
                        min_element = data[t, c, i, l]
                    end

                    c = c + 1
                end
            end
        end
        index = vcat(index, v)
    end
    index = index*get_base_power(system)
    c = 1
    for (l, k) in enumerate(keys(lmps_df[collect(keys(lmps_df))[1]]))
        palette = :Dark2_8
        for b in 1:length(bus_name)
            for t in 1:length(period)
                if length(lmps_df[collect(keys(lmps_df))[1]]) > 1 && k == "RT"
                    if c == 1
                        if bus_name[b] == node
                            plot(
                                index,
                                data[t, b + 1, :, l];
                                label="hour: " *
                                    string(period[t] - 1) *
                                    " - " *
                                    string(bus_name[b]) *
                                    " - " *
                                    k,
                                legend=:outertopright,
                                linestyle=:dash,
                                color = "black",
                                linewidth=3,
                            )
                        else
                            plot(
                                index,
                                data[t, b + 1, :, l];
                                label="hour: " *
                                    string(period[t] - 1) *
                                    " - " *
                                    string(bus_name[b]) *
                                    " - " *
                                    k,
                                legend=:outertopright,
                                linestyle=:dash,
                                palette=palette,
                            )
                        end
                    else
                        if bus_name[b] == node
                            plot!(
                                index,
                                data[t, b + 1, :, l];
                                label="hour: " *
                                    string(period[t] - 1) *
                                    " - " *
                                    string(bus_name[b]) *
                                    " - " *
                                    k,
                                legend=:outertopright,
                                linestyle=:dash,
                                color = "black",
                                linewidth=3,
                            )
                        else
                            plot!(
                                index,
                                data[t, b + 1, :, l];
                                label="hour: " *
                                    string(period[t] - 1) *
                                    " - " *
                                    string(bus_name[b]) *
                                    " - " *
                                    k,
                                legend=:outertopright,
                                linestyle=:dash,
                                palette=palette,
                            )
                        end
                    end
                else
                    if c == 1
                        if bus_name[b] == node
                            plot(
                                index,
                                data[t, b + 1, :, l];
                                label="hour: " *
                                    string(period[t] - 1) *
                                    " - " *
                                    string(bus_name[b]) *
                                    " - " *
                                    k,
                                legend=:outertopright,
                                color = "black",
                                linewidth=3,
                            )
                        else
                            plot(
                                index,
                                data[t, b + 1, :, l];
                                label="hour: " *
                                    string(period[t] - 1) *
                                    " - " *
                                    string(bus_name[b]) *
                                    " - " *
                                    k,
                                legend=:outertopright,
                                palette=palette,
                            )
                        end
                    else
                        if bus_name[b] == node
                            plot!(
                                index,
                                data[t, b + 1, :, l];
                                label="hour: " *
                                    string(period[t] - 1) *
                                    " - " *
                                    string(bus_name[b]) *
                                    " - " *
                                    k,
                                legend=:outertopright,
                                color = "black",
                                linewidth=3,
                            )
                        else
                            plot!(
                                index,
                                data[t, b + 1, :, l];
                                label="hour: " *
                                    string(period[t] - 1) *
                                    " - " *
                                    string(bus_name[b]) *
                                    " - " *
                                    k,
                                legend=:outertopright,
                                palette=palette,
                            )
                        end
                    end
                end
                c = c + 1
            end
        end
    end

    if ylimit == true
        return plot!(;
        title="Price per Virtual Bid on " * node,
        ylabel="Prices (\$/MWh)",
        xlabel="Bid offers (MW)",
        ylims=(min_element - 1, max_element * 1.1 + 1),
    )
    else
        return plot!(;
        title="Price per Virtual Bid on " * node,
        ylabel="Prices (\$/MWh)",
        xlabel="Bid offers (MW)",
    )
    end

    
end

"""
    plot_revenue_curves(
        market_simulator::UCEDRT,
        lmps_df::Dict{Any,Any},
        results_df::Dict{Any,Any},
        period::Vector{Int64},
        generator_name::String,
        initial_time::Date,
        system::System,
        ylimit::Bool,
    )

Function to plot the revenue curve for the the virtual offer bids. 
The 'generator_name' defines which is the virtual generator that we want to plot it's results
and 'periods' controls which periods we want to include in the plot.
"""
function plot_revenue_curves(
    market_simulator::UCEDRT,
    lmps_df::Dict{Any,Any},
    results_df::Dict{Any,Any},
    period::Vector{Int64},
    generator_name::String,
    initial_time::Date,
    system::System,
    ylimit::Bool,
)
    lmps_df = sort(lmps_df)
    gen = get_component(ThermalStandard, market_simulator.system_uc, generator_name)
    bus_name = get_name(get_bus(gen))

    index = []
    aux_period = []
    min_element = 0
    max_element = 0
    for t in period
        aux_period = vcat(aux_period, DateTime(initial_time) + Hour(t - 1))
    end
    data = Array{Any}(nothing, (length(period), 2, length(lmps_df)))
    price = zeros(length(keys(lmps_df[collect(keys(lmps_df))[1]])))
    price = Dict()
    for k in keys(lmps_df[collect(keys(lmps_df))[1]])
        price[k] = 0
    end

    for (i, v) in enumerate(keys(lmps_df))
        for t in 1:length(period)
            data[t, 1, i] = aux_period[t]
            variable_results = read_realized_variables(
                get_problem_results(results_df[v]["DA"], "UC"); names=[:P__ThermalStandard]
            )
            generator_data = getindex.(Ref(variable_results), [:P__ThermalStandard])
            virtual_gen = generator_data[1][!, generator_name][[period[t]]][1]

            for k in (keys(lmps_df[collect(keys(lmps_df))[1]]))
                try
                    global prices_hour = lmps_df[v][k][
                        aux_period[t] .<= lmps_df[v][k].DateTime .< aux_period[t] + Hour(1),
                        bus_name,
                    ]
                catch
                    global prices_hour = lmps_df[v][k][
                        aux_period[t] .<= lmps_df[v][k].DateTime .< aux_period[t] + Hour(1),
                        "lmp",
                    ]
                end
                price[k] = sum(prices_hour; dims=1)[1]
            end
            data[t, 2, i] = (price["DA"] - price["RT"]) * virtual_gen

            if data[t, 2, i] > max_element && data[t, 2, i] < 1e3
                max_element = data[t, 2, i]
            elseif data[t, 2, i] < min_element && data[t, 2, i] > -1e3
                min_element = data[t, 2, i]
            end
        end
        index = vcat(index, v)
    end
    palette = :Dark2_8
    index = index*get_base_power(system)
    c = 1
    for t in 1:length(period)
        if c == 1
            plot(
                index,
                data[t, 2, :];
                label="hour: " * string(period[t] - 1),
                legend=:outertopright,
                palette=palette,
            )
        else
            plot!(
                index,
                data[t, 2, :];
                label="hour: " * string(period[t] - 1),
                legend=:outertopright,
                palette=palette,
            )
        end
        c = c + 1
    end

    if ylimit == true
        return plot!(;
            title="Virtual Revenue per Offer on " * bus_name,
            ylabel="Revenue (\$)",
            xlabel="Bid offers (MW)",
            ylims=(min_element - 1, max_element * 1.1 + 1),
        )
    else
        return plot!(;
        title="Virtual Revenue per Offer on " * bus_name,
        ylabel="Revenue (\$)",
        xlabel="Bid offers (MW)",
    )
    end
end

"""
    plot_sum_revenue_curves(
        market_simulator::UCEDRT,
        lmps_df::Dict{Any,Any},
        results_df::Dict{Any,Any},
        period::Vector{Int64},
        generator_name::String,
        initial_time::Date,
    )

Function to plot the total revenue curve for the the virtual offer bids. 
The 'generator_name' defines which is the virtual generator that we want to plot it's results
and 'periods' controls which periods we want to include in the sum of the plot.
"""
function plot_sum_revenue_curves(
    market_simulator::UCEDRT,
    lmps_df::Dict{Any,Any},
    results_df::Dict{Any,Any},
    period::Vector{Int64},
    generator_name::String,
    initial_time::Date,
)
    lmps_df = sort(lmps_df)
    gen = get_component(ThermalStandard, market_simulator.system_uc, generator_name)
    bus_name = get_name(get_bus(gen))

    indices = []
    aux_period = []
    min_element = 0
    max_element = 0
    for t in period
        aux_period = vcat(aux_period, DateTime(initial_time) + Hour(t - 1))
    end
    data = Array{Any}(nothing, (length(period), 2, length(lmps_df)))
    price = zeros(length(keys(lmps_df[collect(keys(lmps_df))[1]])))
    price = Dict()
    for k in keys(lmps_df[collect(keys(lmps_df))[1]])
        price[k] = 0
    end

    for (i, v) in enumerate(keys(lmps_df))
        for t in 1:length(period)
            data[t, 1, i] = aux_period[t]
            variable_results = read_realized_variables(
                get_problem_results(results_df[v]["DA"], "UC"); names=[:P__ThermalStandard]
            )
            generator_data = getindex.(Ref(variable_results), [:P__ThermalStandard])
            virtual_gen = generator_data[1][!, generator_name][[period[t]]][1]

            for k in (keys(lmps_df[collect(keys(lmps_df))[1]]))
                try
                    global prices_hour = lmps_df[v][k][
                        aux_period[t] .<= lmps_df[v][k].DateTime .< aux_period[t] + Hour(1),
                        bus_name,
                    ]
                catch
                    global prices_hour = lmps_df[v][k][
                        aux_period[t] .<= lmps_df[v][k].DateTime .< aux_period[t] + Hour(1),
                        "lmp",
                    ]
                end
                price[k] = sum(prices_hour; dims=1)[1]
            end
            data[t, 2, i] = (price["DA"] - price["RT"]) * virtual_gen

            if data[t, 2, i] > max_element && data[t, 2, i] < 1e3
                max_element = data[t, 2, i]
            elseif data[t, 2, i] < min_element && data[t, 2, i] > -1e3
                min_element = data[t, 2, i]
            end
        end
        indices = vcat(indices, v)
    end
   
    palette = :Dark2_8
    plot(
        indices,
        sum(data[t, 2, :] for t=1:length(period));
        label=false,
        legend=:outertopright,
        palette=palette,
    )

    return plot!(;
        title="Total Virtual Revenue per Offer on " * bus_name,
        ylabel="Revenue (\$)",
        xlabel="Bid offers (p.u)",
        ylims=(min_element - 1, max_element * 1.1 + 1),
    )
end

"""
    plot_revenue_curves(
        market_simulator::UCED,
        lmps_df::Dict{Any,Any},
        results_df::Dict{Any,Any},
        period::Vector{Int64},
        generator_name::String,
        initial_time::Date,
        system::System,
        ylimit::Bool,
    )

Function to plot the revenue curve for the the virtual offer bids. 
The 'generator_name' defines which is the virtual generator that we want to plot it's results
and 'periods' controls which periods we want to include in the plot.
"""
function plot_revenue_curves(
    market_simulator::UCED,
    lmps_df::Dict{Any,Any},
    results_df::Dict{Any,Any},
    period::Vector{Int64},
    generator_name::String,
    initial_time::Date,
    system::System,
    ylimit::Bool,
)
    lmps_df = sort(lmps_df)
    gen = get_component(ThermalStandard, market_simulator.system_uc, generator_name)
    bus_name = get_name(get_bus(gen))

    index = []
    aux_period = []
    min_element = 0
    max_element = 0
    for t in period
        aux_period = vcat(aux_period, DateTime(initial_time) + Hour(t - 1))
    end
    data = Array{Any}(nothing, (length(period), 2, length(lmps_df)))
    price = zeros(length(keys(lmps_df[collect(keys(lmps_df))[1]])))
    price = Dict()
    for k in keys(lmps_df[collect(keys(lmps_df))[1]])
        price[k] = 0
    end

    for (i, v) in enumerate(keys(lmps_df))
        for t in 1:length(period)
            data[t, 1, i] = aux_period[t]
            variable_results = read_realized_variables(
                get_problem_results(results_df[v]["DA"], "UC"); names=[:P__ThermalStandard]
            )
            generator_data = getindex.(Ref(variable_results), [:P__ThermalStandard])
            virtual_gen = generator_data[1][!, generator_name][[period[t]]][1]

            for k in (keys(lmps_df[collect(keys(lmps_df))[1]]))
                prices_hour = lmps_df[v][k][
                    aux_period[t] .<= lmps_df[v][k].DateTime .< aux_period[t] + Hour(1),
                    bus_name,
                ]
                price[k] = sum(prices_hour; dims=1)[1]
            end
            data[t, 2, i] = price["DA"] * virtual_gen

            if data[t, 2, i] > max_element && data[t, 2, i] < 1e3
                max_element = data[t, 2, i]
            elseif data[t, 2, i] < min_element && data[t, 2, i] > -1e3
                min_element = data[t, 2, i]
            end
        end
        index = vcat(index, v)
    end
    palette = :Dark2_8
    index = index*get_base_power(system)
    c = 1
    for t in 1:length(period)
        if c == 1
            plot(
                index,
                data[t, 2, :];
                label="hour: " * string(period[t] - 1),
                legend=:outertopright,
                palette=palette,
            )
        else
            plot!(
                index,
                data[t, 2, :];
                label="hour: " * string(period[t] - 1),
                legend=:outertopright,
                palette=palette,
            )
        end
        c = c + 1
    end

    if ylimit == true
        return plot!(;
            title="Virtual Revenue per Offer on " * bus_name,
            ylabel="Revenue (\$)",
            xlabel="Bid offers (MW)",
            ylims=(min_element - 1, max_element * 1.1 + 1),
        )
    else
        return plot!(;
        title="Virtual Revenue per Offer on " * bus_name,
        ylabel="Revenue (\$)",
        xlabel="Bid offers (MW)",
    )
    end
end

"""
    plot_revenue_curves_load(
        market_simulator::UCED,
        lmps_df::Dict{Any,Any},
        period::Vector{Int64},
        range_quota::Vector{Float64},
        initial_time::Date,
        load::PowerLoad,
        system::System,
        ylimit::Bool,
    )

Function to plot the revenue curve for the the virtual demand bids. 
The 'load' defines which is the virtual load that we want to plot it's results
and 'period' controls which periods we want to include in the plot.
"""
function plot_revenue_curves_load(
    market_simulator::UCED,
    lmps_df::Dict{Any,Any},
    period::Vector{Int64},
    range_quota::Vector{Float64},
    initial_time::Date,
    load::PowerLoad,
    system::System,
    ylimit::Bool,
)

    lmps_df = sort(lmps_df)
    
    index = []
    aux_period = []
    min_element = 0
    max_element = 0
    start_time = DateTime(initial_time)
    for t in period
        aux_period = vcat(aux_period, start_time + Hour(t - 1))
    end

    data = zeros(length(keys(lmps_df)), length(period))
    price = zeros(length(keys(lmps_df[collect(keys(lmps_df))[1]])))
    price = Dict()
    for k in keys(lmps_df[collect(keys(lmps_df))[1]])
        price[k] = 0
    end

    loads = collect(get_components(PowerLoad, system))
    bus_name = get_name(get_bus(load))
    ts_names = get_time_series_names(SingleTimeSeries, loads[1])
    ts_values = get_time_series_values(Deterministic, load, ts_names[1]; start_time)
    period = period .- 1

    for (i, v) in enumerate(keys(lmps_df))
        for t in 1:length(period)
            for k in (keys(lmps_df[collect(keys(lmps_df))[1]]))
                prices_hour = lmps_df[v][k][
                    aux_period[t] .<= lmps_df[v][k].DateTime .< aux_period[t] + Hour(1),
                    bus_name,
                ]
                price[k] = sum(prices_hour; dims=1)[1]
            end
            data[i,t] = -price["DA"] * ts_values[t] * range_quota[i]
        end
        index = vcat(index, v)
    end

    palette = :Dark2_8
    index = index*get_base_power(system)
    
    c = 1
    for t in 1:length(period)
        if c == 1
            plot(
                index,
                data[:,t];
                label="hour: " * string(period[t]),
                legend=:outertopright,
                palette=palette,
            )
        else
            plot!(
                index,
                data[:,t];
                label="hour: " * string(period[t]),
                legend=:outertopright,
                palette=palette,
            )
        end
        c = c + 1
    end

    if ylimit == true
        return plot!(;
            title="Virtual Revenue per Offer on " * bus_name,
            ylabel="Revenue (\$)",
            xlabel="Bid offers (MW)",
            ylims=(min_element - 1, max_element * 1.1 + 1),
        )
    else
        return plot!(;
        title="Virtual Revenue per Offer on " * bus_name,
        ylabel="Revenue (\$)",
        xlabel="Bid offers (MW)",
    )
    
    end
end

"""
    plot_revenue_curves_load(
        market_simulator::UCEDRT,
        lmps_df::Dict{Any,Any},
        period::Vector{Int64},
        range_quota::Vector{Float64},
        initial_time::Date,
        load::PowerLoad,
        system::System,
        ylimit::Bool,
    )

Function to plot the revenue curve for the the virtual demand bids. 
The 'load' defines which is the virtual load that we want to plot it's results
and 'period' controls which periods we want to include in the plot.
"""
function plot_revenue_curves_load(
    market_simulator::UCEDRT,
    lmps_df::Dict{Any,Any},
    period::Vector{Int64},
    range_quota::Vector{Float64},
    initial_time::Date,
    load::PowerLoad,
    system::System,
    ylimit::Bool,
)

    lmps_df = sort(lmps_df)
    
    index = []
    aux_period = []
    min_element = 0
    max_element = 0
    start_time = DateTime(initial_time)
    for t in period
        aux_period = vcat(aux_period, start_time + Hour(t - 1))
    end

    data = zeros(length(keys(lmps_df)), length(period))
    price = zeros(length(keys(lmps_df[collect(keys(lmps_df))[1]])))
    price = Dict()
    for k in keys(lmps_df[collect(keys(lmps_df))[1]])
        price[k] = 0
    end

    loads = collect(get_components(PowerLoad, system))
    bus_name = get_name(get_bus(load))
    ts_names = get_time_series_names(SingleTimeSeries, loads[1])
    ts_values = get_time_series_values(Deterministic, load, ts_names[1]; start_time)
    period = period .- 1

    for (i, v) in enumerate(keys(lmps_df))
        for t in 1:length(period)
            for k in (keys(lmps_df[collect(keys(lmps_df))[1]]))
                prices_hour = lmps_df[v][k][
                    aux_period[t] .<= lmps_df[v][k].DateTime .< aux_period[t] + Hour(1),
                    bus_name,
                ]
                price[k] = sum(prices_hour; dims=1)[1]
            end
            data[i,t] = (price["RT"]-price["DA"]) * ts_values[t] * range_quota[i]
        end
        index = vcat(index, v)
    end

    palette = :Dark2_8
    index = index*get_base_power(system)
    
    c = 1
    for t in 1:length(period)
        if c == 1
            plot(
                index,
                data[:,t];
                label="hour: " * string(period[t]),
                legend=:outertopright,
                palette=palette,
            )
        else
            plot!(
                index,
                data[:,t];
                label="hour: " * string(period[t]),
                legend=:outertopright,
                palette=palette,
            )
        end
        c = c + 1
    end

    if ylimit == true
        return plot!(;
            title="Virtual Revenue per Offer on " * bus_name,
            ylabel="Revenue (\$)",
            xlabel="Bid offers (MW)",
            ylims=(min_element - 1, max_element * 1.1 + 1),
        )
    else
        return plot!(;
        title="Virtual Revenue per Offer on " * bus_name,
        ylabel="Revenue (\$)",
        xlabel="Bid offers (MW)",
    )
    
    end
end

"""
    plot_revenue_curves_renewable(
        market_simulator::UCEDRT,
        lmps_df::Dict{Any,Any},
        results_df::Dict{Any,Any},
        bids::Vector{Float64},
        generator_name::String,
        node::String,
        ylimit::Bool,
    )

Function to plot the revenue curve for the the renewable generators. 
The 'generator_name' defines which is the generator that we want to plot it's results
and 'bids' controls which bids we want to include in the plot.
"""
function plot_revenue_curves_renewable(
    market_simulator::UCEDRT,
    lmps_df::Dict{Any,Any},
    results_df::Dict{Any,Any},
    bids::Vector{Float64},
    generator_name::String,
    node::String,
    ylimit::Bool,
)
    gen = get_component(RenewableDispatch, market_simulator.system_uc, generator_name)
    bus_name = get_name(get_bus(gen))
    min_element = 0
    max_element = 0
    data = Array{Any}(nothing, (24, length(bids) + 1))
    data[:, 1] = lmps_df[first(keys(lmps_df))]["DA"][!, "DateTime"]
    for (j, q) in enumerate(bids)
        variable_results = read_realized_variables(
            get_problem_results(results_df[q]["DA"], "UC"); names=[:P__RenewableDispatch]
        )
        generator_data = getindex.(Ref(variable_results), [:P__RenewableDispatch])
        gen_da = generator_data[1][!, generator_name]

        variable_results = read_realized_variables(
            get_problem_results(results_df[q]["RT"], "RT"); names=[:P__RenewableDispatch]
        )
        generator_data = getindex.(Ref(variable_results), [:P__RenewableDispatch])
        gen_rt_aux = generator_data[1][!, generator_name]
        gen_rt = []

        for i in 1:(round(Int, length(gen_rt_aux) / 12)) #TODO: Change to horizon
            gen_rt = vcat(
                gen_rt, [sum(gen_rt_aux[(1 + 12 * (i - 1)):(12 + 12 * (i - 1))]) / 12]
            )
        end
        price = Dict()
        for k in (keys(lmps_df[collect(keys(lmps_df))[1]])) #Problems 
            for t in 1:24
                prices_hour = lmps_df[q][k][
                    lmps_df[first(keys(lmps_df))]["DA"][!, "DateTime"][t] .<= lmps_df[q][k].DateTime .< lmps_df[first(keys(lmps_df))]["DA"][!, "DateTime"][t] + Hour(
                        1
                    ),
                    bus_name,
                ]
                if t == 1
                    price[k] = sum(prices_hour)
                else
                    price[k] = vcat(price[k], sum(prices_hour))
                end
            end
        end

        data[:, j + 1] = gen_da .* price["DA"] + (gen_rt - gen_da) .* price["RT"]

        if maximum(data[:, j + 1]) > max_element && minimum(data[:, j + 1]) < 1e3
            max_element = maximum(data[:, j + 1])
        elseif minimum(data[:, j + 1]) < min_element && minimum(data[:, j + 1]) > -1e3
            min_element = minimum(data[:, j + 1])
        end
    end
    palette = :Dark2_8

    c = 1
    for (j, q) in enumerate(bids)
        if c == 1
            plot(
                0:(length(data[:, 1]) - 1),
                data[:, j + 1];
                label="bid: " * string(q),
                legend=:outertopright,
                palette=palette,
            )
        else
            plot!(
                0:(length(data[:, 1]) - 1),
                data[:, j + 1];
                label="bid: " * string(q),
                legend=:outertopright,
                palette=palette,
            )
        end
        c = c + 1
    end

    if ylimit == true
        return plot!(;
            title=generator_name * " Revenue per Virtual Offer on " * node,
            ylabel="Revenue (\$)",
            xlabel="Period",
            ylims=(min_element - 1, max_element * 1.1 + 1),
        )
    else    
        return plot!(;
        title=generator_name * " Revenue per Virtual Offer on " * node,
        ylabel="Revenue (\$)",
        xlabel="Period",
    )
    end
end

"""
    plot_revenue_curves_renewable_plus_virtual(
        market_simulator::UCEDRT,
        lmps_df::Dict{Any,Any},
        results_df::Dict{Any,Any},
        bids::Vector{Float64},
        renewable_gen::String,
        virtual_gen::String,
        ylimit::Bool,
    )

Function to plot the revenue curve for the the renewable and virtual generators. 
The 'renewable_gen' and 'virtual_gen' defines which are the generators that we want to plot it's results
and 'bids' controls which bids we want to include in the plot.
"""
function plot_revenue_curves_renewable_plus_virtual(
    market_simulator::UCEDRT,
    lmps_df::Dict{Any,Any},
    results_df::Dict{Any,Any},
    bids::Vector{Float64},
    renewable_gen::String,
    virtual_gen::String,
    ylimit::Bool,
)
    gen_r = get_component(RenewableDispatch, market_simulator.system_uc, renewable_gen)
    bus_r = get_name(get_bus(gen_r))
    gen_v = get_component(ThermalStandard, market_simulator.system_uc, virtual_gen)
    bus_v = get_name(get_bus(gen_v))
    min_element = 0
    max_element = 0
    data = Array{Any}(nothing, (24, length(bids) + 1, 3))
    data[:, 1, 1] = lmps_df[first(keys(lmps_df))]["DA"][!, "DateTime"]
    data[:, 1, 2] = data[:, 1, 1]
    data[:, 1, 3] = data[:, 1, 1]
    for (j, q) in enumerate(bids)
        variable_results = read_realized_variables(
            get_problem_results(results_df[q]["DA"], "UC"); names=[:P__RenewableDispatch]
        )
        generator_data = getindex.(Ref(variable_results), [:P__RenewableDispatch])
        gen_da_r = generator_data[1][!, renewable_gen]

        variable_results = read_realized_variables(
            get_problem_results(results_df[q]["RT"], "RT"); names=[:P__RenewableDispatch]
        )
        generator_data = getindex.(Ref(variable_results), [:P__RenewableDispatch])
        gen_rt_aux = generator_data[1][!, renewable_gen]
        gen_rt_r = []

        for i in 1:(round(Int, length(gen_rt_aux) / 12)) #TODO: Change "12" to get horizon
            gen_rt_r = vcat(
                gen_rt_r, [sum(gen_rt_aux[(1 + 12 * (i - 1)):(12 + 12 * (i - 1))]) / 12]
            )
        end

        variable_results = read_realized_variables(
            get_problem_results(results_df[q]["DA"], "UC"); names=[:P__ThermalStandard]
        )
        generator_data = getindex.(Ref(variable_results), [:P__ThermalStandard])
        gen_da_v = generator_data[1][!, virtual_gen]

        price = Dict()
        for k in (keys(lmps_df[collect(keys(lmps_df))[1]]))
            for t in 1:24
                prices_hour = lmps_df[q][k][
                    lmps_df[first(keys(lmps_df))]["DA"][!, "DateTime"][t] .<= lmps_df[q][k].DateTime .< lmps_df[first(keys(lmps_df))]["DA"][!, "DateTime"][t] + Hour(
                        1
                    ),
                    :,
                ]
                if t == 1
                    prices_hour[!, "DateTime"] .= lmps_df[first(keys(lmps_df))]["DA"][
                        !, "DateTime"
                    ][t]
                    price[k] = combine(
                        groupby(prices_hour, :DateTime),
                        names(prices_hour, Not(:DateTime)) .=> sum;
                        renamecols=false,
                    )
                else
                    prices_hour[!, "DateTime"] .= lmps_df[first(keys(lmps_df))]["DA"][
                        !, "DateTime"
                    ][t]
                    price[k] = vcat(
                        price[k],
                        combine(
                            groupby(prices_hour, :DateTime),
                            names(prices_hour, Not(:DateTime)) .=> sum;
                            renamecols=false,
                        ),
                    )
                end
            end
        end

        try
            data[:, j + 1, 1] =
                gen_da_r .* price["DA"][!, bus_r] +
                (gen_rt_r - gen_da_r) .* price["RT"][!, bus_r]
            data[:, j + 1, 2] = gen_da_v .* (price["DA"][!, bus_v] - price["RT"][!, bus_v])
        catch
            try
                bus="lmp"
                data[:, j + 1, 1] =
                    gen_da_r .* price["DA"][!, bus] +
                    (gen_rt_r - gen_da_r) .* price["RT"][!, bus]
                data[:, j + 1, 2] = gen_da_v .* (price["DA"][!, bus] - price["RT"][!, bus])
            catch
                bus="lmp"
                data[:, j + 1, 1] =
                gen_da_r .* price["DA"][!, bus] +
                (gen_rt_r - gen_da_r) .* price["RT"][!, bus_r]
                data[:, j + 1, 2] = gen_da_v .* (price["DA"][!, bus] - price["RT"][!, bus_v])
            end
        end
        
        for c in 1:size(data)[1]
            for i in 1:2
                if data[c, j + 1, i] > max_element && data[c, j + 1, i] < 1e3
                    max_element = data[c, j + 1, i]
                end
                if data[c, j + 1, i] < min_element && data[c, j + 1, i] > -1e3
                    min_element = data[c, j + 1, i]
                end
            end
        end
    end
    data[:, 2:(length(bids) + 1), 3] =
        data[:, 2:(length(bids) + 1), 1] + data[:, 2:(length(bids) + 1), 2]
    palette = :Dark2_8
    title = [
        renewable_gen * " Revenue - Virtual Offer on " * bus_v,
        virtual_gen * " Revenue - Virtual Offer on " * bus_v,
        "Joint Revenue - Virtual Offer on " * bus_v,
    ]
    plt = Array{Any}(nothing, (3)) #TODO: Change to typeof(plot): Plots.Plot{Plots.PlotlyBackend}
    for i in 1:3
        c = 1
        for (j, q) in enumerate(bids)
            if i == 1
                label = "bid: " * string(q)
            else
                label = false
            end
            if c == 1
                plt[i] = plot(
                    0:(length(data[:, 1, i]) - 1),
                    data[:, j + 1, i];
                    label=label,
                    legend=:outertopright,
                    palette=palette,
                    title=title[i],
                    titlefont=font(10, "Arial"),
                )
            else
                plot!(
                    plt[i],
                    0:(length(data[:, 1, i]) - 1),
                    data[:, j + 1, i];
                    label=label,
                    legend=:outertopright,
                    palette=palette,
                )
            end
            c = c + 1
        end
    end

    if ylimit == true
        return plot(
            plt...;
            layout=(3, 1),
            ylabel="Revenue (\$)",
            ylims=(min_element - 1, max_element * 1.1 + 1),
        )
    else
        return plot(
            plt...;
            layout=(3, 1),
            ylabel="Revenue (\$)",
        )
    end
end

"""
    plot_revenue_curves_renewable_plus_virtual_load(
        market_simulator::UCEDRT,
        lmps_df::Dict{Any,Any},
        results_df::Dict{Any,Any},
        bids::Vector{Float64},
        renewable_gen::String,
        virtual_load::String,
        ylimit::Bool,
        load::PowerLoad,
    )

Function to plot the revenue curve for the the renewable and virtual generators. 
The 'renewable_gen' and 'virtual_gen' defines which are the generators that we want to plot it's results
and 'bids' controls which bids we want to include in the plot.
"""
function plot_revenue_curves_renewable_plus_virtual_load(
    market_simulator::UCEDRT,
    lmps_df::Dict{Any,Any},
    results_df::Dict{Any,Any},
    bids::Vector{Float64},
    renewable_gen::String,
    virtual_load::String,
    ylimit::Bool,
    load::PowerLoad,
)
    gen_r = get_component(RenewableDispatch, market_simulator.system_uc, renewable_gen)
    bus_r = get_name(get_bus(gen_r))
    load_v = get_component(PowerLoad, market_simulator.system_uc, virtual_load)
    bus_v = get_name(get_bus(load_v))
    min_element = 0
    max_element = 0
    data = Array{Any}(nothing, (24, length(bids) + 1, 3))
    data[:, 1, 1] = lmps_df[first(keys(lmps_df))]["DA"][!, "DateTime"]
    data[:, 1, 2] = data[:, 1, 1]
    data[:, 1, 3] = data[:, 1, 1]
    for (j, q) in enumerate(bids)
        variable_results = read_realized_variables(
            get_problem_results(results_df[q]["DA"], "UC"); names=[:P__RenewableDispatch]
        )
        generator_data = getindex.(Ref(variable_results), [:P__RenewableDispatch])
        gen_da_r = generator_data[1][!, renewable_gen]

        variable_results = read_realized_variables(
            get_problem_results(results_df[q]["RT"], "RT"); names=[:P__RenewableDispatch]
        )
        generator_data = getindex.(Ref(variable_results), [:P__RenewableDispatch])
        gen_rt_aux = generator_data[1][!, renewable_gen]
        gen_rt_r = []

        for i in 1:(round(Int, length(gen_rt_aux) / 12)) #TODO: Change "12" to get horizon
            gen_rt_r = vcat(
                gen_rt_r, [sum(gen_rt_aux[(1 + 12 * (i - 1)):(12 + 12 * (i - 1))]) / 12]
            )
        end

        loads = collect(get_components(PowerLoad, market_simulator.system_uc))
        ts_names = get_time_series_names(SingleTimeSeries, loads[1])
        start_time = data[:,:,1][1,1]
        ts_values = get_time_series_values(Deterministic, load, ts_names[1]; start_time)
        gen_da_v = ts_values

        price = Dict()
        for k in (keys(lmps_df[collect(keys(lmps_df))[1]]))
            for t in 1:24
                prices_hour = lmps_df[q][k][
                    lmps_df[first(keys(lmps_df))]["DA"][!, "DateTime"][t] .<= lmps_df[q][k].DateTime .< lmps_df[first(keys(lmps_df))]["DA"][!, "DateTime"][t] + Hour(
                        1
                    ),
                    :,
                ]
                if t == 1
                    prices_hour[!, "DateTime"] .= lmps_df[first(keys(lmps_df))]["DA"][
                        !, "DateTime"
                    ][t]
                    price[k] = combine(
                        groupby(prices_hour, :DateTime),
                        names(prices_hour, Not(:DateTime)) .=> sum;
                        renamecols=false,
                    )
                else
                    prices_hour[!, "DateTime"] .= lmps_df[first(keys(lmps_df))]["DA"][
                        !, "DateTime"
                    ][t]
                    price[k] = vcat(
                        price[k],
                        combine(
                            groupby(prices_hour, :DateTime),
                            names(prices_hour, Not(:DateTime)) .=> sum;
                            renamecols=false,
                        ),
                    )
                end
            end
        end

        data[:, j + 1, 1] =
            gen_da_r .* price["DA"][!, bus_r] +
            (gen_rt_r - gen_da_r) .* price["RT"][!, bus_r]
        data[:, j + 1, 2] = gen_da_v .* (price["RT"][!, bus_v] - price["DA"][!, bus_v])
        for c in 1:size(data)[1]
            for i in 1:2
                if data[c, j + 1, i] > max_element && data[c, j + 1, i] < 1e3
                    max_element = data[c, j + 1, i]
                end
                if data[c, j + 1, i] < min_element && data[c, j + 1, i] > -1e3
                    min_element = data[c, j + 1, i]
                end
            end
        end
    end
    data[:, 2:(length(bids) + 1), 3] =
        data[:, 2:(length(bids) + 1), 1] + data[:, 2:(length(bids) + 1), 2]
    palette = :Dark2_8
    title = [
        renewable_gen * " Revenue - Virtual Offer on " * bus_v,
        virtual_load * " Revenue - Virtual Offer on " * bus_v,
        "Joint Revenue - Virtual Offer on " * bus_v,
    ]
    plt = Array{Any}(nothing, (3)) #TODO: Change to typeof(plot): Plots.Plot{Plots.PlotlyBackend}
    for i in 1:3
        c = 1
        for (j, q) in enumerate(bids)
            if i == 1
                label = "bid: " * string(q)
            else
                label = false
            end
            if c == 1
                plt[i] = plot(
                    0:(length(data[:, 1, i]) - 1),
                    data[:, j + 1, i];
                    label=label,
                    legend=:outertopright,
                    palette=palette,
                    title=title[i],
                    titlefont=font(10, "Arial"),
                )
            else
                plot!(
                    plt[i],
                    0:(length(data[:, 1, i]) - 1),
                    data[:, j + 1, i];
                    label=label,
                    legend=:outertopright,
                    palette=palette,
                )
            end
            c = c + 1
        end
    end

    if ylimit == true
        return plot(
            plt...;
            layout=(3, 1),
            ylabel="Revenue (\$)",
            ylims=(min_element - 1, max_element * 1.1 + 1),
        )
    else
        return plot(
            plt...;
            layout=(3, 1),
            ylabel="Revenue (\$)",
        )
    end
end

"""
    plot_revenue_curves(
        market_simulator::UCEDRT,
        lmps_df::Dict{Any,Any},
        results_df::Dict{Any,Any},
        period::Vector{Int64},
        generator_name::String,
        initial_time::Date,
        system::System,
    )

Function to plot the virtual generation curve for the virtual offer bids. 
The 'generator_name' defines which is the virtual generator that we want to plot it's results
and 'periods' controls which periods we want to include in the plot.
"""

function plot_generation_curves(
    market_simulator,
    lmps_df::Dict{Any,Any},
    results_df::Dict{Any,Any},
    period::Vector{Int64},
    generator_name::String,
    initial_time::Date,
    system::System,
)
    lmps_df = sort(lmps_df)
    gen = get_component(ThermalStandard, market_simulator.system_uc, generator_name)
    bus_name = get_name(get_bus(gen))
    aux_period = []
    for t in period
        aux_period = vcat(aux_period, DateTime(initial_time) + Hour(t - 1))
    end

    index = []
    data = Array{Any}(nothing, (length(period), 2, length(lmps_df)))
    for (i, v) in enumerate(keys(lmps_df))
        for t in 1:length(period)
            data[t, 1, i] = aux_period[t]
            variable_results = read_realized_variables(
                get_problem_results(results_df[v]["DA"], "UC"); names=[:P__ThermalStandard]
            )
            generator_data = getindex.(Ref(variable_results), [:P__ThermalStandard])
            virtual_gen = generator_data[1][!, generator_name][[period[t]]][1]
            data[t, 2, i] = virtual_gen
        end
        index = vcat(index, v)
    end
    palette = :Dark2_8
    index = index*get_base_power(system)
    c = 1
    for t in 1:length(period)
        if c == 1
            plot(
                index,
                data[t, 2, :];
                label="hour: " * string(period[t] - 1),
                legend=:outertopright,
                palette=palette,
            )
        else
            plot!(
                index,
                data[t, 2, :];
                label="hour: " * string(period[t] - 1),
                legend=:outertopright,
                palette=palette,
            )
        end
        c = c + 1
    end

    return plot!(;
        title=generator_name * " generation per Offer on " * bus_name,
        ylabel="Generation(MWh)",
        xlabel="Bid offers (MW)",
    )
end

"""
plot_revenue_curves_renewable(
    lmps_df::Dict{Any,Any},
    results_df::Dict{Any,Any},
    bids::Vector{Float64},
    generator_name::String,
    node::String,
)

Function to plot the renewable generation curve for the virtual offer bids. 
The 'generator_name' defines which is the renewable generator that we want to plot it's results
and 'bids' controls which possible virtual bids we want to include in the plot.
"""

function plot_generation_curves_renewable(
    lmps_df::Dict{Any,Any},
    results_df::Dict{Any,Any},
    bids::Vector{Float64},
    generator_name::String,
    node::String,
)
    data = Array{Any}(nothing, (24, length(bids) + 1, 2))
    data[:, 1, 1] = lmps_df[first(keys(lmps_df))]["DA"][!, "DateTime"]
    data[:, 1, 2] = lmps_df[first(keys(lmps_df))]["DA"][!, "DateTime"]
    for (j, q) in enumerate(bids)
        variable_results = read_realized_variables(
            get_problem_results(results_df[q]["DA"], "UC"); names=[:P__RenewableDispatch]
        )
        generator_data = getindex.(Ref(variable_results), [:P__RenewableDispatch])
        gen_da = generator_data[1][!, generator_name]

        variable_results = read_realized_variables(
            get_problem_results(results_df[q]["RT"], "RT"); names=[:P__RenewableDispatch]
        )
        generator_data = getindex.(Ref(variable_results), [:P__RenewableDispatch])
        gen_rt_aux = generator_data[1][!, generator_name]
        gen_rt = []

        for i in 1:(round(Int, length(gen_rt_aux) / 12)) #TODO: Change "12" to get the horizon
            gen_rt = vcat(
                gen_rt, [sum(gen_rt_aux[(1 + 12 * (i - 1)):(12 + 12 * (i - 1))]) / 12]
            )
        end

        data[:, j + 1, 1] = gen_da
        data[:, j + 1, 2] = gen_rt
    end

    c = 1
    for i in 1:2
        palette = :Dark2_8
        for (j, q) in enumerate(bids)
            if i == 1
                if c == 1
                    plot(
                        0:(length(data[:, 1, i]) - 1),
                        data[:, j + 1, i];
                        label="bid: " * string(q),
                        legend=:outertopright,
                        palette=palette,
                    )
                else
                    plot!(
                        0:(length(data[:, 1, i]) - 1),
                        data[:, j + 1, i];
                        label="bid: " * string(q),
                        legend=:outertopright,
                        palette=palette,
                    )
                end
                c = c + 1
            else
                if c == 1
                    plot(
                        0:(length(data[:, 1, i]) - 1),
                        data[:, j + 1, i];
                        label="bid: " * string(q),
                        linestyle=:dash,
                        legend=:outertopright,
                        palette=palette,
                    )
                else
                    plot!(
                        0:(length(data[:, 1, i]) - 1),
                        data[:, j + 1, i];
                        label="bid: " * string(q),
                        linestyle=:dash,
                        legend=:outertopright,
                        palette=palette,
                    )
                end
                c = c + 1
            end
        end
    end

    return plot!(;
        title=generator_name * " Generation per Virtual Offer on " * node,
        ylabel="Generation(MWh)",
        xlabel="Period",
    )
end

"""
    heat_map_revenue_curves_mix(
        market_simulator::UCEDRT,
        lmps_df::Dict{Any,Any},
        period::Vector{Int64},
        range_quota::Vector{Float64},
        initial_time::Date,
        load::PowerLoad,
        system::System,
        ylimit::Bool,
    )

    Function to plot the virtual revenue heat map for the mix of virtual INC and DEC bids.
"""

function heat_map_revenue_curves_mix(
    market_simulator::UCEDRT,
    lmps_df::Dict{Any,Any},
    results_df::Dict{Any,Any},
    period::Vector{Int64},
    range_quota_load::Vector{Float64},
    range_quota_gen::Vector{Float64},
    initial_time::Date,
    load::PowerLoad,
    generator_name::String,
    system::System,
)

    gen = get_component(ThermalStandard, market_simulator.system_uc, generator_name)
    bus_gen = get_name(get_bus(gen))

    loads = collect(get_components(PowerLoad, system))
    bus_load = get_name(get_bus(load))
    ts_names = get_time_series_names(SingleTimeSeries, loads[1])
    start_time = DateTime(initial_time)
    ts_values = get_time_series_values(Deterministic, load, ts_names[1]; start_time)#./range_quota_load[length(range_quota_load)]

    aux_period = []

    for t in period
        aux_period = vcat(aux_period, start_time + Hour(t - 1))
    end
    
    data = zeros(length(range_quota_load),length(range_quota_gen), length(period))
    price = zeros(length(keys(lmps_df[collect(keys(lmps_df))[1]])))
    price = Dict()
    for k in keys(lmps_df[collect(keys(lmps_df))[1]])
        price[k] = 0
    end

    for (i,l) in enumerate(range_quota_load)
        for (j,g) in enumerate(range_quota_gen)
            for t in 1:length(period)
                variable_results = read_realized_variables(
                    get_problem_results(results_df[[l g]]["DA"], "UC"); names=[:P__ThermalStandard]
                )
                generator_data = getindex.(Ref(variable_results), [:P__ThermalStandard])
                virtual_gen = generator_data[1][!, generator_name][[period[t]]][1]
                for k in (keys(lmps_df[collect(keys(lmps_df))[1]]))
                    prices_hour = lmps_df[[l g]][k][
                        aux_period[t] .<= lmps_df[[l g]][k].DateTime .< aux_period[t] + Hour(1),
                        :,
                    ]
                    prices_hour[!, "DateTime"] .= lmps_df[first(keys(lmps_df))]["DA"][
                        !, "DateTime"
                    ][t]
                    price[k] = combine(
                        groupby(prices_hour, :DateTime),
                        names(prices_hour, Not(:DateTime)) .=> sum;
                        renamecols=false,
                    )
                end
                revenue_load = ((price["RT"][!,bus_load]-price["DA"][!,bus_load]) * ts_values[t] * range_quota_load[i])[1]
                revenue_gen = ((price["DA"][!,bus_gen]-price["RT"][!,bus_gen]) * virtual_gen)[1] 
                data[i,j,t] = revenue_load + revenue_gen
            end
        end
    end
    data_sum = sum(data[:,:,t] for t =1:length(period)) 
    data_norm = data_sum
    gr()
    data_p = data_norm 
    h=heatmap(range_quota_gen.*100,
    range_quota_load.*100, data_p,
    c=cgrad([:blue, :white,:red, :yellow]),
    xlabel="INC Offer (MW/h)", ylabel="DEC bid (MW/h)",
    title="Virtual Revenue (\$)")
    
    return plot(h)
end

"""
    heat_map_coal_generation(
        system::System,
        results_df::Dict{Any,Any},
        range_quota_load::Vector{Float64},
        range_quota_gen::Vector{Float64},
        initial_time::Date,
        period::Vector{Int64}=[1],
        generator_fields::AbstractArray=[:P__ThermalStandard, :P__RenewableDispatch],
    )

Function to plot the coal generation heat map for the the mix of virtual INC and DEC bids.
"""
function heat_map_coal_generation(
    system::System,
    results_df::Dict{Any,Any},
    range_quota_load::Vector{Float64},
    range_quota_gen::Vector{Float64},
    initial_time::Date,
    period::Vector{Int64}=[1],
    generator_fields::AbstractArray=[:P__ThermalStandard, :P__RenewableDispatch],
)
    data=zeros(length(range_quota_load),length(range_quota_gen), length(period))
    start_time = DateTime(initial_time)
    aux_period = []

    for t in period
        aux_period = vcat(aux_period, start_time + Hour(t - 1))
    end

    for (z,t) in enumerate(aux_period)
        for (x,l) in enumerate(range_quota_load)
            for (y,g) in enumerate(range_quota_gen)
                system_results = get_problem_results(results_df[[l g]]["RT"], "RT")
                # get mapping from busname to fuel type
                fuel_type_dict = fuel_type_mapping(system)

                # get the output data for given fuel types
                variable_results = read_realized_variables(system_results; names=generator_fields)
                generator_data = getindex.(Ref(variable_results), generator_fields)
                for i in 1:length(generator_data)
                    generation = generator_data[i][t .<= generator_data[i].DateTime .< t + Hour(1), :]
                    generation[!, "DateTime"] .= t
                    generator_data[i] = combine(
                        groupby(generation, :DateTime),
                        names(generation, Not(:DateTime)) .=> sum;
                        renamecols=false,
                    )
                    col = size(generator_data[i])[2]
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
        
                #= select rows for the given bus names, default to all buses.
                if !isempty(bus_names)
                    bus_names = String.(bus_names)
                    @assert all(bus_names .∈ [stacked_data.bus_name])
                    filter!(:bus_name => in(bus_names), stacked_data)
                end
                =#
        
                stacked_data.fuel_type = get.(Ref(fuel_type_dict), stacked_data.gen_name, missing)
        
                aggregated_data = combine(
                    groupby(stacked_data, [:DateTime, :fuel_type]), :output => sum => :output
                )
        
                # convert the output units into MWh.
                aggregated_data.output = aggregated_data.output .* get_base_power(system)
        
                # unstack aggregated data and make area plot 
                global unstacked_data = unstack(aggregated_data, :fuel_type, :output)
                #=
                if i==1 && j==1
                    global data_frame = unstacked_data
                else
                    global data_frame = append!(data_frame, unstacked_data)
                end
                =#
                data[x,y,z]=unstacked_data[!,"COAL"][1]
            end
        end
    end
    data_sum = sum(data[:,:,t] for t =1:length(period)) 
    data_norm = data_sum./data_sum[1,1]
    gr()
    data_p = data_norm 
    h=heatmap(range_quota_gen.*100,
    range_quota_load.*100, data_p,
    c=cgrad([:blue, :white,:red, :yellow]),
    xlabel="INC Offer (MW/h)", ylabel="DEC bid (MW/h)",
    title="Coal Emission (%)")
    
    return plot(h)
end
