"""
    build_5_bus_matpower_DA(data_dir::AbstractString, case_file::AbstractString)

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
    pm_data = PowerSystems.PowerModelsData(case_file_path)

    tsp = InfrastructureSystems.read_time_series_file_metadata(forecasts_pointers_file)

    sys = System(pm_data)
    if add_reserves
        reserves = [
            VariableReserve{ReserveUp}("REG1", true, 5.0, 0.1),
            VariableReserve{ReserveUp}("REG2", true, 5.0, 0.06),
            VariableReserve{ReserveUp}("REG3", true, 5.0, 0.03),
            VariableReserve{ReserveUp}("REG4", true, 5.0, 0.02),
        ]
        contributing_devices = get_components(Generator, sys)
        for r in reserves
            add_service!(sys, r, contributing_devices)
        end
    end

    add_time_series!(sys, tsp)

    return sys
end

"""
    prep_systems_UCED(system::System)

Duplicates the system to represent UC and ED for DA, transforming the time series 
to the appropriate interval and horizon.

PS.: Beacuse of a bug in PSI, we have to set the horizon of the ED problem, that 
would normally be of 1 hour, to 2 hours.
"""
function prep_systems_UCED(
    system::System;
    horizon_uc::Int=24,
    horizon_ed::Int=2,
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
    get_duals(<:AbstractPowerModel)

Return the duals for the specified formulation.
"""
get_duals(::Type{CopperPlatePowerModel}) = [:CopperPlateBalance]
get_duals(::Union{Type{NFAPowerModel}, Type{DCPPowerModel}}) = [:nodal_balance_active__Bus]
get_duals(::Type{StandardPTDFModel}) = [:CopperPlateBalance, :network_flow__Line]
