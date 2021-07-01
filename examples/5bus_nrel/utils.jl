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
    transform_single_time_series!(sys, 12, Minute(15))

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
