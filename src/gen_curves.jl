"""
    pq_curves_da(generator::ThermalStandard; base_system::System,range_quota::Vector{Int64}, initial_time::Date, steps::Int = 1, simulation_folder::String = pwd(), solver_uc::MathOptInterface.OptimizerWithAttributes, solver_ed::MathOptInterface.OptimizerWithAttributes)

Creates a curve of generation and nodes prices for a vector of 'range_quotas' for max generation of a virtual 'generator' in a 'base_system'. 
"""
function pq_curves_da(
    generator::ThermalStandard;
    base_system::System,
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
        set_active_power_limits!(generator, (min = 0.0, max = max_gen))
       
        # duplicate system and prepare times series for the time varying parameters (loads, renewables, ...)
        sys_uc, sys_ed = prep_systems_UCED(base_system)

        # generic market formulation templates with defined network formulation
        # CopperPlate-OPF: network=CopperPlatePowerModel
        # DC-OPF: network=DCPPowerModel
        # NFA-OPF (only line limit constraints): network=NFAPowerModel
        # DC-PTDF-OPF (what ISOs do): network=StandardPTDFModel
        template_uc = template_unit_commitment(; network=DCPPowerModel)
        template_ed = template_economic_dispatch(; network=DCPPowerModel)

        # build a market clearing simulator (run `@doc UCED` for more information)
        market_simulator = UCED(;
            system_uc=sys_uc,
            system_ed=sys_ed,
            template_uc=template_uc,
            template_ed=template_ed,
            solver_uc=solver_uc,
            solver_ed=solver_ed,
        )

        @test isa(market_simulator, UCED)
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
