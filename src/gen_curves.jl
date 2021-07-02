
"""
    set_active_power_limits(market_simulator::UCED, name_generator::AbstractString, active_power_limits::NamedTuple{(:min, :max), Tuple{Float64, Float64}})

Set 'active_power_limits' to a generator whose name is 'name_generator' if the clearing market is UCED.
"""

function set_active_power_limits!(
    market_simulator::UCED,
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
    set_active_power_limits(market_simulator::UCRT, name_generator::AbstractString, active_power_limits::NamedTuple{(:min, :max), Tuple{Float64, Float64}})

Set 'active_power_limits' to a generator whose name is 'name_generator' if the clearing market is UCRT.
"""

function set_active_power_limits!(
    market_simulator::UCRT,
    name_generator::AbstractString,
    active_power_limits::NamedTuple{(:min, :max),Tuple{Float64,Float64}},
)
    generator_uc = get_component(
        ThermalStandard, market_simulator.system_uc, name_generator
    )
    generator_rt = get_component(
        ThermalStandard, market_simulator.system_rt, name_generator
    )

    PowerSystems.set_active_power_limits!(generator_uc, active_power_limits)
    return PowerSystems.set_active_power_limits!(generator_rt, active_power_limits)
end

"""
    set_active_power_limits(market_simulator::UCEDRT, name_generator::AbstractString, active_power_limits::NamedTuple{(:min, :max), Tuple{Float64, Float64}})

Set 'active_power_limits' to a generator whose name is 'name_generator' if the clearing market is UCEDRT.
"""

function set_active_power_limits!(
    market_simulator::UCEDRT,
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
    pq_curves_virtuals(market_simulator::MarketSimulator, name_generator::AbstractString, range_quota::Vector{Int64}, initial_time::Date, steps::Int = 1, simulation_folder::String = pwd())

Creates a curve of generation and nodes prices for a vector 'range_quotas' that represents the max generation of a virtual 'generator' in a market clearing 'market_simulator'. 
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
