
"""
    create_generator_bids(initial_bidding_time::DateTime, bidding_periods::Vector{Int}, system::System, costs::Vector{T}) where {T<:Real}

Function to help user to creat single generator bids from 'initial_bidding_time' for all periods defined in 'bidding_periods' with the apropriate costs defined in 'costs'.
"""
function create_generator_bids(
    initial_bidding_time::DateTime,
    bidding_periods::Vector{Int},
    system::System,
    costs::Vector{T},
) where {T<:Real}
    params = get_time_series_params(system)

    initial_time = params.initial_timestamp
    interval = params.interval
    count = params.count
    frequency = Minute(system.frequency)

    MultiDay = collect(
        initial_time:frequency:((initial_time + count * interval) - frequency)
    )

    bidding_periods = (bidding_periods .- 1) .* frequency + initial_bidding_time
    bids = Dict{DateTime,T}()
    for t in 1:length(bidding_periods)
        bids[bidding_periods[t]] = costs[t]
    end

    for t in setdiff(MultiDay, bidding_periods)
        bids[t] = 1e7
    end

    return SingleTimeSeries(; name="variable_cost", data=bids, resolution=frequency)
end

#= TODO: Finish implementation 
function add_bidding_generators(system, nodal_bids)

    generators = add_generator!.(system, nodes, active_power_limits)
    add_timeSeries!.(generators, costs) 
    

end

"""
    add_gerator!(system::System, node::String, active_power_limits::NamedTuple{(:min, :max), Tuple{Real, Real}})

Function to creat and add generator to the system following an especified node with a defined active power limits.
"""
function add_gerator!(system::System, node::String, active_power_limits::NamedTuple{(:min, :max), Tuple{T, T}}) where {T<:Real}
    bus = get_component(Bus, system, node)
    gen = ThermalStandard(;
        name=get_name(bus)*"_virtual_supply",
        available=true,
        status=true,
        bus=bus,
        active_power=0.0, 
        reactive_power=0.0,
        rating=active_power_limits.max,
        prime_mover=PrimeMovers.ST,
        fuel=ThermalFuels.COAL, #TO DO: Creat virtual thermalfuel
        active_power_limits=active_power_limits,
        reactive_power_limits=(min=-0.2, max=0.2), # won't influence our simulations
        time_limits=nothing,
        ramp_limits=nothing,
        operation_cost=ThreePartCost([(0.0, 0.0), (20.0, 1.0)], 0.0, 0.0, 0.0), 
        base_power=get_base_power(system),
    )
    add_component!(system, gen)
    return gen
end


function add_timeSeries!(generators, costs)


    remove_time_series!(sys, DeterministicSingleTimeSeries)

    for gen in generators #get_components(ThermalGen, sys)
        varcost = get_operation_cost(gen)
        data = TimeArray(MultiDay, repeat([get_cost(get_variable(varcost))], 8784))
        _time_series = SingleTimeSeries("variable_cost", data)
        add_time_series!(sys, gen, _time_series)
        #set_variable_cost!(sys, gen, _time_series)
    end
end

function add_bids_da!(sys, bids::T, bid_initial_time, gen_name) where T
    horizon = get_horizon(sys)
    initial_time = get_initial_timestamp(sys)
    interval = get_interval(sys)
    count = get_count(sys)
    generator = get_component(ThermalStandard, sys, gen_name)
    numB = length(bids)
    data = Dict{DateTime,T}()
    for t = 0:interval.value:numB-horizon
        data[bid_initial_time+t*Hour(1)] = bids[t+1:t+horizon]
    end
    kys = keys(data)
    for t in setdiff(initial_time:interval:initial_time+interval*(count-1), kys)
        data[t] = fill(1e7, horizon)
    end
    time_series_data = Deterministic(
        name = "variable_cost",
        data = data,
        resolution = Dates.Hour(1)
    )
    set_variable_cost!(sys, generator, time_series_data)
end
=#
