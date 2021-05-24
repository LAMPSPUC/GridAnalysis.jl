"""
    fuel_type_mapping(system)

    return a mapping between the bus name and the fuel type for the given `system`.
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
    # we only want the first hour prices (stored in the first row)
    # in case each problem has more than one hour
    df = vcat([DataFrame(duals[i][1, :]) for i in keys(duals)]...)
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
    # we only want the first hour of the ED prices (stored in the first row here)
    flow_duals = vcat([DataFrame(flow_duals[i][1, :]) for i in keys(flow_duals)]...)
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
    evaluate_prices(market_simulator::UCED, problem_results::PSI.SimulationProblemResults)

Returns energy prices for the simulation's data-range.  
"""
function evaluate_prices(
    market_simulator::UCED, problem_results::PSI.SimulationProblemResults
)
    return evaluate_prices(
        market_simulator.template_ed.transmission,
        market_simulator.system_ed,
        problem_results,
        market_simulator.kwargs,
    )
end

"""
    get_time_series_params(system::System)

Returns the parameters associated with the time-series attached to the system.
"""
get_time_series_params(system::System) = system.data.time_series_params.forecast_params
