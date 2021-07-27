
"""
    set_active_power_limits!(
        market_simulator::MarketSimulator,
        name_generator::AbstractString,
        active_power_limits::NamedTuple{(:min, :max),Tuple{Float64,Float64}},
    )

Set `active_power_limits` to a generator whose name is `name_generator` to any clearing market.
"""

function set_active_power_limits!(
    market_simulator::MarketSimulator,
    name_generator::AbstractString,
    active_power_limits::NamedTuple{(:min, :max),Tuple{Float64,Float64}},
)
    generator_uc = get_component(
        ThermalStandard, market_simulator.system_uc, name_generator
    )
    generator_ed = get_component(
        ThermalStandard, market_simulator.system_ed, name_generator
    )

    PowerSystems.set_active_power_limits!(generator_uc, active_power_limits)
    return PowerSystems.set_active_power_limits!(generator_ed, active_power_limits)
end

"""
    pq_curves_virtuals!(
        market_simulator::MarketSimulator,
        name_generator::AbstractString,
        range_quota::Vector{Float64},
        initial_time::Date,
        steps::Int=1,
        simulation_folder::String=pwd(),
    )

Creates a curve of generation and nodes prices for a vector `range_quotas` that represents the max 
generation of a virtual `generator` in a market clearing `market_simulator`. 
"""
function pq_curves_virtuals!(
    market_simulator::MarketSimulator,
    name_generator::AbstractString,
    range_quota::Vector{Float64},
    initial_time::Date,
    steps::Int=1,
    simulation_folder::String=pwd(),
)
    lmps_df = Dict()
    results_df = Dict()

    for max_gen in range_quota
        set_active_power_limits!(market_simulator, name_generator, (min=0.0, max=max_gen))

        # for each formulation you will need to save different dual variables:
        constraint_duals = duals_constraint_names(market_simulator)

        # Simulate market
        # build and run simulation
        results = run_multiday_simulation(
            market_simulator,
            initial_time, # initial time for simulation
            steps; # number of steps in simulation (normally number of days to simulate)
            services_slack_variables=true,
            balance_slack_variables=true,
            constraint_duals=constraint_duals,
            name="quota_$max_gen",
            simulation_folder=simulation_folder,
        )

        # results
        results_df[max_gen] = results
        if isa(market_simulator, UCEDRT) #TODO: evaluate_prices unified
            lmps_df[max_gen] = evaluate_prices_UCEDRT(market_simulator, results)
        else
            lmps_df[max_gen] = evaluate_prices(market_simulator, results)
        end
    end
    return lmps_df, results_df
end

"""
    set_active_power_max!(
        market_simulator::MarketSimulator,
        name_load::AbstractString,
        active_power_max::Float64,
    )

Set `active_power_max` to a load whose name is `name_load` for any clearing market.
"""

function set_active_power_max!(
    market_simulator::MarketSimulator,
    name_load::AbstractString,
    active_power_max::Float64,
)
    load_uc = get_component(
        PowerLoad, market_simulator.system_uc, name_load
    )
    load_ed = get_component(
        PowerLoad, market_simulator.system_ed, name_load
    )

    PowerSystems.set_max_active_power!(load_uc, active_power_max)
    return PowerSystems.set_max_active_power!(load_ed, active_power_max)
end

"""
    pq_curves_load_virtuals!(
        market_simulator::MarketSimulator,
        name_load::AbstractString,
        range_quota::Vector{Float64},
        initial_time::Date,
        steps::Int=1,
        simulation_folder::String=pwd(),
    )

Creates a curve of load and nodes prices for a vector `range_quota` that represents 
the max load of a virtual `load` in a market clearing `market_simulator`. 
"""
function pq_curves_load_virtuals!(
    market_simulator::MarketSimulator,
    name_load::AbstractString,
    range_quota::Vector{Float64},
    initial_time::Date,
    steps::Int=1,
    simulation_folder::String=pwd(),
)
    lmps_df = Dict()
    results_df = Dict()

    for max in range_quota
        set_active_power_max!(market_simulator, name_load, max)

        # for each formulation you will need to save different dual variables:
        constraint_duals = duals_constraint_names(market_simulator)

        # Simulate market
        # build and run simulation
        results = run_multiday_simulation(
            market_simulator,
            initial_time, # initial time for simulation
            steps; # number of steps in simulation (normally number of days to simulate)
            services_slack_variables=true,
            balance_slack_variables=true,
            constraint_duals=constraint_duals,
            name="quota_$max",
            simulation_folder=simulation_folder,
        )

        # results
        results_df[max] = results
        if isa(market_simulator, UCEDRT) #TODO: evaluate_prices unified
            lmps_df[max] = evaluate_prices_UCEDRT(market_simulator, results)
        else
            lmps_df[max] = evaluate_prices(market_simulator, results)
        end
    end
    return lmps_df, results_df
end


"""
    pq_curves_load_gen_virtuals!(
        market_simulator::MarketSimulator,
        name_load::AbstractString,
        range_quota_load::Vector{Float64},
        initial_time::Date,
        steps::Int=1,
        name_generator::AbstractString,
        range_quota_gen::Vector{Float64},
        simulation_folder::String=pwd(),
    )

Creates a curve of generation and nodes prices for two vectors of 'range_quota' that represents 
the max load of a virtual 'load' and the max generation of a virtual 'generator' in a market clearing 'market_simulator'. 
"""

function pq_curves_load_gen_virtuals!(
    market_simulator::MarketSimulator,
    name_load::AbstractString,
    range_quota_load::Vector{Float64},
    initial_time::Date,
    name_generator::AbstractString,
    range_quota_gen::Vector{Float64},
    steps::Int=1,
    simulation_folder::String=pwd(),
)
    lmps_df = Dict()
    results_df = Dict()

    for max_l in range_quota_load
        for max_gen in range_quota_gen
            set_active_power_max!(market_simulator, name_load, max_l)
            set_active_power_limits!(market_simulator, name_generator, (min=0.0, max=max_gen))

            # for each formulation you will need to save different dual variables:
            constraint_duals = duals_constraint_names(market_simulator)

            # Simulate market
            # build and run simulation
            results = run_multiday_simulation(
                market_simulator,
                initial_time, # initial time for simulation
                steps; # number of steps in simulation (normally number of days to simulate)
                services_slack_variables=true,
                balance_slack_variables=true,
                constraint_duals=constraint_duals,
                name="quota_l_$max_l"*"_quota_g_$max_gen",
                simulation_folder=simulation_folder,
            )

            # results
            results_df[[max_l max_gen]] = results
            if isa(market_simulator, UCEDRT) #TODO: evaluate_prices unified
                lmps_df[[max_l max_gen]] = evaluate_prices_UCEDRT(market_simulator, results)
            else
                lmps_df[[max_l max_gen]] = evaluate_prices(market_simulator, results)
            end
        end
    end
    return lmps_df, results_df
end
