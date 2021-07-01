"""
    prep_systems_UCED(
        system::System;
        horizon_uc::Int=24,
        horizon_ed::Int=1,
        interval_uc::TimePeriod=Hour(24),
        interval_ed::TimePeriod=Hour(1),
    )

Duplicates the system to represent UC and ED for DA, transforming the time series
to the appropriate interval and horizon.
"""
function prep_systems_UCED(
    system::System;
    horizon_uc::Int=24,
    horizon_ed::Int=1,
    interval_uc::TimePeriod=Hour(24),
    interval_ed::TimePeriod=Hour(1),
)
    system_uc = system
    system_ed = deepcopy(system)

    transform_single_time_series!(system_uc, horizon_uc, interval_uc)
    transform_single_time_series!(system_ed, horizon_ed, interval_ed)

    return system_uc, system_ed
end

"""
    plot_prices_RT_hour(prices, ylim::Tuple)

Plot the RT prices with all simulator forecast
"""
function plot_prices_RT_hour(prices, ylim::Tuple)
    prices_keys = collect(keys(prices))
    prices_rt_df = prices[prices_keys[1]]

    values = select(prices_rt_df, Not(:DateTime))
    n_prev, n_bus = size(values)

    intervals = get_time_series_params(market_simulator.system_rt).interval

    n_prev_hour = Int(60 / intervals.value)
    n_days = Int(n_prev / n_prev_hour)

    names_bus = names(values)
    prices_rt = zeros(n_days, n_bus)

    i = 1
    while i < n_days
        for j in 1:(n_prev_hour):n_prev
            prices_hour = prices_rt_df[
                prices_rt_df[j, :DateTime] .<= prices_rt_df.DateTime .< prices_rt_df[j, :DateTime] + Hour(
                    1
                ),
                :,
            ]
            prices_hour = select(prices_hour, Not(:DateTime))
            prices_rt[i, :] = sum(Matrix(prices_hour); dims=1)
            i = i + 1
        end
    end

    times = prices_rt_df[1:n_prev_hour:n_prev, 1]

    labels = permutedims(names_bus)

    return plot(
        prices_rt;
        legend=:outertopright,
        xlab="Hours",
        ylab="Prices (\$/MWh)",
        label=labels,
        ylim=ylim,
    )
end

"""
    plot_prices_hour(prices, ylim::Tuple)

Plot Both ED and RT in the same plot
"""
function plot_DA_RT(prices, ylim::Tuple)
    prices_keys = collect(keys(prices))
    prices_rt_df = prices[prices_keys[1]]

    values_rt = select(prices_rt_df, Not(:DateTime))
    n_prev, n_bus = size(values_rt)

    intervals = get_time_series_params(market_simulator.system_rt).interval

    n_prev_hour = Int(60 / intervals.value)
    n_days = Int(n_prev / n_prev_hour)

    names_bus = names(values_rt)
    prices_rt = zeros(n_days, n_bus)

    i = 1
    while i < n_days
        for j in 1:(n_prev_hour):n_prev
            prices_hour = prices_rt_df[
                prices_rt_df[j, :DateTime] .<= prices_rt_df.DateTime .< prices_rt_df[j, :DateTime] + Hour(
                    1
                ),
                :,
            ]
            prices_hour = select(prices_hour, Not(:DateTime))
            prices_rt[i, :] = sum(Matrix(prices_hour); dims=1)
            i = i + 1
        end
    end

    times = prices_rt_df[1:n_prev_hour:n_prev, 1]

    labels = permutedims(names_bus)

    prices_ed_df = prices[prices_keys[2]]
    values_ed = select(prices_ed_df, Not(:DateTime))
    values_ed = Matrix(values_ed)

    palette = ["RoyalBlue", "Aquamarine", "DeepPink", "Coral", "Green"]

    plot(
        prices_rt;
        legend=:outertopright,
        xlab="Hours",
        ylab="Prices (\$/MWh)",
        label=labels,
        linestyle=:dash,
        palette=palette,
        ylim=ylim,
    )
    return plot!(
        values_ed;
        legend=:outertopright,
        xlab="Hours",
        ylab="Prices (\$/MWh)",
        label=labels,
        palette=palette,
        ylim=ylim,
    )
end
