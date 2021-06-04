

function set_active_power_limits(
    market_simulator::UCED, name_generator, active_power_limits
)

        generator_uc = get_component(name_generator) #em cada sistema 
        generator_ed = get_component(name_generator)

        set_active_power_limits!(generator_uc, active_power_limits)
        set_active_power_limits!(generator_da, active_power_limits)

    return generator_uc, generator_ed
end

"""
    pq_curves_da(generator::ThermalStandard; base_system::System,range_quota::Vector{Int64}, initial_time::Date, steps::Int = 1, simulation_folder::String = pwd(), solver_uc::MathOptInterface.OptimizerWithAttributes, solver_ed::MathOptInterface.OptimizerWithAttributes)

Creates a curve of generation and nodes prices for a vector of 'range_quotas' for max generation of a virtual 'generator' in a 'base_system'. 
"""
function pq_curves_virtuals(
    name_generator::AbstractString,
    market_simulator::,
    range_quota::Vector{Int64},
    initial_time::Date,
    steps::Int = 1,
    simulation_folder::String = pwd(),
    solver_uc,
    solver_ed
)
    lmps_df = Dict()
    results_df = Dict()

    for max_gen in range_quota

        set_active_power_limits!(market_simulator, name_generator, (min = 0.0, max = max_gen))


        # for each formulation you will need to save different dual variables:
        constraint_duals = duals_constraint_names(market_simulator)

        @test isa(constraint_duals, AbstractVector{Symbol})

        # Simulate market
        # build and run simulation
        results = run_multiday_simulation(
            market_simulator,
            initial_time, # initial time for simulation
            steps; # number of steps in simulation (normally number of days to simulate)
            services_slack_variables=false,
            balance_slack_variables=false,
            constraint_duals=constraint_duals,
            #=name="test_case_5bus",
            simulation_folder=mktempdir(), # Locally can use: joinpath(example_dir, "results"),
            =#
            name = "quota_$max_gen",
            simulation_folder = simulation_folder#,
            #kwargs...
        );

        @test isa(results, SimulationResults)

        # calculate prices
        ed_results = get_problem_results(results, "ED");
        prices = evaluate_prices(market_simulator, ed_results)

        @test isa(prices, DataFrame)
        
        # results
        results_df[max_gen] = get_problem_results(results, "UC");
        lmps_df[max_gen] = prices

    end
    return lmps_df, results_df

end
