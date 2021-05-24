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
        name=name,
        steps=steps,
        problems=problems,
        sequence=uc_ed_sequence,
        simulation_folder=simulation_folder,
        initial_time=initial_time,
    )

    return sim
end

"""
    run_multiday_simulation(simulator::UCED, initial_time::Date, steps::Int) -> SimulationResults

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
