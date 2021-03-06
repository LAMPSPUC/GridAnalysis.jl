
"""
    create_generator_bids(;
    initial_bidding_time::DateTime,
    bidding_periods::Vector{Int},
    system::System,
    costs::Vector{T},
    ) where {T<:AbstractFloat}

Creates `SingleTimeSeries` for a generator `variable_cost`. Allows the user to define bidding periods 
(`bidding_periods`) relative to a `initial_bidding_time` where `costs` will be applied.
Every other period (not defined by the user) is set to `1e7 \$/MWh`.
"""
function create_generator_bids(;
    initial_bidding_time::DateTime,
    bidding_periods::Vector{Int},
    system::System,
    costs::Vector{T},
) where {T<:AbstractFloat}
    load = first(get_components(PowerLoad, system))
    timestamps = get_time_series_timestamps(SingleTimeSeries, load, "max_active_power")
    @assert initial_bidding_time in timestamps

    bidding_stamps =
        (bidding_periods .- 1) .+ findall(x -> x == initial_bidding_time, timestamps)[1]
    ts_costs = fill(1e7, length(timestamps))
    for i in 1:length(bidding_stamps)
        t = bidding_stamps[i]
        ts_costs[t] = costs[i]
    end
    bids = TimeArray(timestamps, ts_costs)
    ts_array = SingleTimeSeries(; name="variable_cost", data=bids)
    return ts_array
end

"""
    add_generator!(
        system::System, node::String, active_power_limits::NamedTuple{(:min, :max),Tuple{T,T}}
    ) where {T<:AbstractFloat}

Function to create and add generator to the system following an especified node with a defined 
active power limits.
"""
function add_generator!(
    system::System, node::String, active_power_limits::NamedTuple{(:min, :max),Tuple{T,T}}
) where {T<:AbstractFloat}
    bus = get_component(Bus, system, node)
    gen = ThermalStandard(;
        name=get_name(bus) * "_virtual_supply",
        available=true,
        status=true,
        bus=bus,
        active_power=0.0,
        reactive_power=0.0,
        rating=active_power_limits.max,
        prime_mover=PrimeMovers.ST,
        fuel=ThermalFuels.OTHER, #TODO: Creat virtual thermalfuel
        active_power_limits=active_power_limits,
        reactive_power_limits=(min=-0.0, max=0.0), # won't influence our simulations
        time_limits=(up=0.0, down=0.0),
        ramp_limits=(up=9999.0, down=9999.0),
        operation_cost=MarketBidCost(;
            no_load=0.0, start_up=(hot=0.0, warm=0.0, cold=0.0), shut_down=0.0
        ),
        base_power=get_base_power(system),
    )
    add_component!(system, gen)
    return gen
end


"""
    add_load!(
        system::System, 
        node::String, 
        max_active_power::Float64
    )

Function to create and add a load to the system following an especified node with a defined 
maximum active power.
"""
function add_load!(
    system::System,
    node::String, 
    max_active_power::Float64
)
    bus = get_component(Bus, system, node)
    load = PowerLoad(;
        name=get_name(bus) * "_virtual_load",
        available=true,
        bus=bus,
        model=LoadModels.ConstantPower,
        active_power=0.0,
        reactive_power=0.0,
        base_power=get_base_power(system),
        max_active_power=max_active_power,
        max_reactive_power=0.0, # won't influence our simulations
    )
    add_component!(system, load)
    return load
end

"""
    create_demand_series(;
        initial_bidding_time::DateTime,
        bidding_periods::Vector{Int},
        system::System,
        demands::Vector{T},
    ) where {T<:AbstractFloat}

Creates `SingleTimeSeries` for a load `max_active_power`. 
Allows the user to define bidding periods (`bidding_periods`) relative to a `initial_bidding_time` 
where `max_active_power` will be applied.
Every other period (not defined by the user) is set to `0.0 \$/MWh`.
"""
function create_demand_series(;
    initial_bidding_time::DateTime,
    bidding_periods::Vector{Int},
    system::System,
    demands::Vector{T},
) where {T<:AbstractFloat}

    load = first(get_components(PowerLoad, system))
    timestamps = get_time_series_timestamps(SingleTimeSeries, load, "max_active_power")
    @assert initial_bidding_time in timestamps

    bidding_stamps =
        (bidding_periods .- 1) .+ findall(x -> x == initial_bidding_time, timestamps)[1]
    ts_demand = fill(0.0, length(timestamps))
    for i in 1:length(bidding_stamps)
        t = bidding_stamps[i]
        ts_demand[t] = demands[i]
    end
    bids = TimeArray(timestamps, ts_demand)
    ts_array = SingleTimeSeries(; name="max_active_power", data=bids, scaling_factor_multiplier = get_max_active_power,)
    return ts_array
end

#  TODO: Finish implementation 
# function add_bidding_generators(system, nodal_bids)

#     generators = add_generator!.(system, nodes, active_power_limits)
#     add_timeSeries!.(generators, costs) 

# end
