"""
    build_5_bus_matpower_DA(
        data_dir::AbstractString;
        case_file::AbstractString="case5_re_uc.m",
        FORECASTS_DIR::AbstractString=joinpath(data_dir, "forecasts"),
        forecasts_pointers_file::AbstractString=joinpath(
            FORECASTS_DIR, "timeseries_pointers_da_7day.json"
        ),
        add_reserves::Bool=true,
    )

Builds base system for the 5bus NREL case (a.k.a NESTA case) from:
 - A matpower file containing grid information (case_file);
 - A file describing forecasts locations and details (forecasts_pointers_file);
 - A simple definition of reserves.
"""
function build_5_bus_matpower_DA(
    data_dir::AbstractString;
    case_file::AbstractString="case5_re_uc.m",
    FORECASTS_DIR::AbstractString=joinpath(data_dir, "forecasts"),
    forecasts_pointers_file::AbstractString=joinpath(
        FORECASTS_DIR, "timeseries_pointers_da_7day.json"
    ),
    add_reserves::Bool=true,
)
    case_file_path = joinpath(data_dir, case_file)
    sys = System(case_file_path)
    if add_reserves
        reserves = [
            VariableReserve{ReserveUp}("REG1", true, 5.0, 0.2),
            VariableReserve{ReserveUp}("REG2", true, 5.0, 0.4),
            VariableReserve{ReserveUp}("REG3", true, 5.0, 0.3),
            VariableReserve{ReserveUp}("REG4", true, 5.0, 0.4),
            VariableReserve{ReserveDown}("REG1d", true, 5.0, 0.2),
            VariableReserve{ReserveDown}("REG2d", true, 5.0, 0.4),
            VariableReserve{ReserveDown}("REG3d", true, 5.0, 0.3),
            VariableReserve{ReserveDown}("REG4d", true, 5.0, 0.4),
        ]
        contributing_devices = get_components(Generator, sys)
        for r in reserves
            add_service!(sys, r, contributing_devices)
        end
    end

    add_time_series!(sys, forecasts_pointers_file)

    return sys
end

"""
    build_5_bus_matpower_RT(
        data_dir::AbstractString;
        case_file::AbstractString="case5_re_uc.m",
        FORECASTS_DIR::AbstractString=joinpath(data_dir, "forecasts"),
        forecasts_pointers_file::AbstractString=joinpath(
            FORECASTS_DIR, "timeseries_pointers_rt_7day.json"
        ),
    )

Builds base system for the 5bus NREL case (a.k.a NESTA case) from:
 - A matpower file containing grid information (case_file);
 - A file describing forecasts locations and details (forecasts_pointers_file);
"""
function build_5_bus_matpower_RT(
    data_dir::AbstractString;
    case_file::AbstractString="case5_re_uc.m",
    FORECASTS_DIR::AbstractString=joinpath(data_dir, "forecasts"),
    forecasts_pointers_file::AbstractString=joinpath(
        FORECASTS_DIR, "timeseries_pointers_rt_7day.json"
    ),
)
    case_file_path = joinpath(data_dir, case_file)
    sys = System(case_file_path)

    add_time_series!(sys, forecasts_pointers_file)
    transform_single_time_series!(sys, 12, Minute(5))

    return sys
end

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
function plot_prices_RT_hour(prices)
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
        prices_rt; legend=:outertopright, xlab="Hours", ylab="Prices (\$/MWh)", label=labels
    )
end

"""
    plot_prices_hour(prices, ylim::Tuple)

Plot both ED and RT in the same plot
"""
function plot_DA_RT(prices)
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
    )
    return plot!(
        values_ed;
        legend=:outertopright,
        xlab="Hours",
        ylab="Prices (\$/MWh)",
        label=labels,
        palette=palette,
    )
end
