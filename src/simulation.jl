"""
    Simulation(simulator::UCED, initial_time::Date, steps::Int) -> PSI.Simulation

Builds a multiday UC-ED PSI.Simulation from a UCED market simulator, an initial date and number of simulation steps.
"""
function PSI.Simulation(
    simulator::UCED,
    initial_time::Date,
    steps::Int=1;
    system_to_file::Bool=false,
    services_slack_variables::Bool=true,
    balance_slack_variables::Bool=true,
    constraint_duals::Array{Symbol,1}=[:CopperPlateBalance, :network_flow],
    name::String="test_case",
    simulation_folder=pwd(),
)
    problems = SimulationProblems(;
        UC=OperationsProblem(
            simulator.template_uc,
            simulator.system_uc;
            optimizer=simulator.solver_uc,
            system_to_file=system_to_file,
            simulator.kwargs...,
        ),
        ED=OperationsProblem(
            simulator.template_ed,
            simulator.system_ed;
            optimizer=simulator.solver_ed,
            system_to_file=system_to_file,
            services_slack_variables=services_slack_variables,
            balance_slack_variables=balance_slack_variables,
            constraint_duals=constraint_duals,
            simulator.kwargs...,
        ),
    )

    feedforward_chronologies = Dict(("UC" => "ED") => Synchronize(; periods=24))

    feedforward = Dict(
        ("ED", :devices, :ThermalStandard) =>
            SemiContinuousFF(; binary_source_problem=ON, affected_variables=[ACTIVE_POWER]),
    )

    intervals = Dict("UC" => (Hour(24), Consecutive()), "ED" => (Hour(1), Consecutive()))

    uc_ed_sequence = SimulationSequence(;
        problems=problems,
        intervals=intervals,
        ini_cond_chronology=InterProblemChronology(),
        feedforward_chronologies=feedforward_chronologies,
        feedforward=feedforward,
    )

    sim = PSI.Simulation(;
        name="ed_"*name,
        steps=steps,
        problems=problems,
        sequence=uc_ed_sequence,
        simulation_folder=simulation_folder,
        initial_time=initial_time,
    )

    return sim
end

"""
    Simulation(simulator::UCRT, initial_time::Date, steps::Int) -> PSI.Simulation

Builds a multiday UC-RT PSI.Simulation from a UCRT market simulator, an initial date and 
number of simulation steps.
"""
function PSI.Simulation(
    simulator::UCRT,
    initial_time::Date,
    steps::Int=1;
    system_to_file::Bool=false,
    services_slack_variables::Bool=true,
    balance_slack_variables::Bool=true,
    constraint_duals::Array{Symbol,1}=[:CopperPlateBalance, :network_flow],
    name::String="test_case",
    simulation_folder=pwd(),
)
    problems = SimulationProblems(;
        UC=OperationsProblem(
            simulator.template_uc,
            simulator.system_uc;
            optimizer=simulator.solver_uc,
            system_to_file=system_to_file,
            simulator.kwargs...,
        ),
        RT=OperationsProblem(
            simulator.template_rt,
            simulator.system_rt;
            optimizer=simulator.solver_rt,
            system_to_file=system_to_file,
            services_slack_variables=services_slack_variables,
            balance_slack_variables=balance_slack_variables,
            constraint_duals=constraint_duals,
            simulator.kwargs...,
        ),
    )

    feedforward_chronologies = Dict(("UC" => "RT") => Synchronize(; periods=24))

    feedforward = Dict(
        ("RT", :devices, :ThermalStandard) =>
            SemiContinuousFF(; binary_source_problem=ON, affected_variables=[ACTIVE_POWER]),
    )

    intervals = Dict("UC" => (Hour(24), Consecutive()), "RT" => (Minute(5), Consecutive()))

    uc_rt_sequence = SimulationSequence(;
        problems=problems,
        intervals=intervals,
        ini_cond_chronology=InterProblemChronology(),
        feedforward_chronologies=feedforward_chronologies,
        feedforward=feedforward,
    )

    sim = PSI.Simulation(;
        name="rt_"*name,
        steps=steps,
        problems=problems,
        sequence=uc_rt_sequence,
        simulation_folder=simulation_folder,
        initial_time=initial_time,
    )

    return sim
end

"""
    Simulation(simulator::UCEDRT, initial_time::Date, steps::Int) -> PSI.Simulation

Builds a multiday UC-ED and a UC-RT PSI.Simulation from a UCEDRT market simulator, an initial date and 
number of simulation steps. It returns a tuple with both simulations.
"""
function PSI.Simulation(
    simulator::UCEDRT,
    initial_time::Date,
    steps::Int=1;
    system_to_file::Bool=false,
    services_slack_variables::Bool=true,
    balance_slack_variables::Bool=true,
    constraint_duals::Vector{Vector{Symbol}},
    name::String="test_case",
    simulation_folder=pwd(),
)
    problem1 = SimulationProblems(;
        UC=OperationsProblem(
            simulator.template_uc,
            simulator.system_uc;
            optimizer=simulator.solver_uc,
            system_to_file=system_to_file,
            simulator.ext...,
        ),
        ED=OperationsProblem(
            simulator.template_ed,
            simulator.system_ed;
            optimizer=simulator.solver_ed,
            system_to_file=system_to_file,
            services_slack_variables=services_slack_variables,
            balance_slack_variables=balance_slack_variables,
            constraint_duals=constraint_duals[1],
            simulator.ext...,
        ),
    )

    feedforward_chronologies = Dict(("UC" => "ED") => Synchronize(; periods=24))

    feedforward = Dict(
        ("ED", :devices, :ThermalStandard) =>
            SemiContinuousFF(; binary_source_problem=ON, affected_variables=[ACTIVE_POWER]),
    )

    intervals = Dict("UC" => (Hour(24), Consecutive()), "ED" => (Hour(1), Consecutive()))

    uc_ed_sequence = SimulationSequence(;
        problems=problem1,
        intervals=intervals,
        ini_cond_chronology=InterProblemChronology(),
        feedforward_chronologies=feedforward_chronologies,
        feedforward=feedforward,
    )

    sim = PSI.Simulation(;
        name="uc_"*name,
        steps=steps,
        problems=problem1,
        sequence=uc_ed_sequence,
        simulation_folder=simulation_folder,
        initial_time=initial_time,
    )

    problem2 = SimulationProblems(;
        UC=OperationsProblem(
            simulator.template_uc,
            simulator.system_uc;
            optimizer=simulator.solver_uc,
            system_to_file=system_to_file,
            simulator.ext...,
        ),
        RT=OperationsProblem(
            simulator.template_rt,
            simulator.system_rt;
            optimizer=simulator.solver_rt,
            system_to_file=system_to_file,
            services_slack_variables=services_slack_variables,
            balance_slack_variables=balance_slack_variables,
            constraint_duals=constraint_duals[2],
            simulator.ext...,
        ),
    )

    feedforward_chronologies2 = Dict(("UC" => "RT") => Synchronize(; periods=24))

    feedforward2 = Dict(
        ("RT", :devices, :ThermalStandard) =>
            SemiContinuousFF(; binary_source_problem=ON, affected_variables=[ACTIVE_POWER]),
    )

    intervals2 = Dict("UC" => (Hour(24), Consecutive()), "RT" => (Minute(5), Consecutive()))

    uc_rt_sequence = SimulationSequence(;
        problems=problem2,
        intervals=intervals2,
        ini_cond_chronology=InterProblemChronology(),
        feedforward_chronologies=feedforward_chronologies2,
        feedforward=feedforward2,
    )

    sim2 = PSI.Simulation(;
        name="rt_"*name,
        steps=steps,
        problems=problem2,
        sequence=uc_rt_sequence,
        simulation_folder=simulation_folder,
        initial_time=initial_time,
    )

    return sim, sim2
end

"""
    run_multiday_simulation(simulator::MarketSimulator, initial_time::Date, steps::Int) -> SimulationResults

Runs a multiday PSI.Simulation from a MarketSimulator, an initial date and number of simulation steps.
"""
function run_multiday_simulation(
    simulator::MarketSimulator,
    initial_time::Date,
    steps::Int=1;
    system_to_file::Bool=false,
    services_slack_variables::Bool=false,
    balance_slack_variables::Bool=false,
    constraint_duals::Array{Symbol,1}=[:CopperPlateBalance, :network_flow],
    name::String="test_case",
    simulation_folder=pwd(),
    console_level=Logging.Warn,
    recorders=[:simulation],
)
    sim = Simulation(
        simulator,
        initial_time,
        steps;
        simulation_folder=simulation_folder,
        system_to_file=system_to_file,
        services_slack_variables=services_slack_variables,
        balance_slack_variables=balance_slack_variables,
        constraint_duals=constraint_duals,
        name=name,
    )

    build!(sim; console_level=console_level, recorders=recorders)

    execute!(sim)

    sim_results = SimulationResults(sim)

    return sim_results
end

"""
    run_multiday_simulation(simulator::MarketSimulator, initial_time::Date, steps::Int) -> SimulationResults

Runs two multiday PSI.Simulation from a MarketSimulator, an initial date and number of simulation steps,
one for the UC-ED simulation and other for the UC-RT simulation. It returns a dictionary with both 
simulation results.
"""
function run_multiday_simulation(
    simulator::UCEDRT,
    initial_time::Date,
    steps::Int=1;
    system_to_file::Bool=false,
    services_slack_variables::Bool=false,
    balance_slack_variables::Bool=false,
    constraint_duals::Vector{Vector{Symbol}},
    name::String="test_case",
    simulation_folder=pwd(),
    console_level=Logging.Warn,
    recorders=[:simulation],
)
    sim1, sim2 = Simulation(
        simulator,
        initial_time,
        steps;
        simulation_folder=simulation_folder,
        system_to_file=system_to_file,
        services_slack_variables=services_slack_variables,
        balance_slack_variables=balance_slack_variables,
        constraint_duals=constraint_duals,
        name=name,
    )

    build!(sim1; console_level=console_level, recorders=recorders)

    execute!(sim1)

    sim_results_1 = SimulationResults(sim1)

    build!(sim2; console_level=console_level, recorders=recorders)

    execute!(sim2)

    sim_results_2 = SimulationResults(sim2)

    return Dict("ED" => sim_results_1, "RT" => sim_results_2)
end
