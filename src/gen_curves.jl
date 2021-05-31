function pq_curves_da(
    sys;
    generator,
    range_quota,
    initial_time,
    market_simulator,
    steps = 1,
    simulation_folder = pwd(),
    kwargs...
)
    lmps_df = Dict()
    results_df = Dict()
    for max_gen in range_quota
        set_active_power_limits!(generator, (min = 0.0, max = max_gen))

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
            simulation_folder = simulation_folder,
            name = "quota_$max_gen",
            kwargs...
        );

        @test isa(results, SimulationResults)

        # calculate prices
        ed_results = get_problem_results(results, "ED");
        prices = evaluate_prices(market_simulator, ed_results)

        @test isa(prices, DataFrame)
        
        # results
        results_df[max_gen] = get_problem_results(results, "UC");
        lmps_df[max_gen] = prices

        #= virtual generation
        bus=get_bus(generator)
        p=prices[!,:bus]

        variable_results = read_realized_variables(uc_results, names=[:P__ThermalStandard])
        generator_data = getindex.(Ref(variable_results), [:P__ThermalStandard])

        vitual_gen=generator_data[1][!,:7] #?
        revenue=p.*virtual_gen
        =#

        #=
        results_df[max_gen] = get_problem_results(results, "UC")
        ed_results = get_problem_results(results, "ED")
        ptdf = PTDF(sys_ed)
        lmps_df[max_gen] = psi_ptdf_lmps(ed_results, sys_ed, ptdf)
        =#
    end
    return lmps_df, results_df

end