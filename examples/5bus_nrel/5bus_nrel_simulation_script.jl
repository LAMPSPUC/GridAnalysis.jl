using Cbc
using Dates
using DataFrames
using GLPK
using GridAnalysis
using InfrastructureSystems
using PowerSystems
using PowerSimulations
using PowerSystemCaseBuilder

include("utils.jl") # case utilities

#' Data Prep and Build Market Simulator
# define solvers for UC and ED
solver_uc = optimizer_with_attributes(Cbc.Optimizer, "logLevel" => 1, "ratioGap" => 0.5)
solver_ed = optimizer_with_attributes(GLPK.Optimizer)

# call our data preparation to build base system
# corrently your REPL has to be on the folder for this example
base_system = build_5_bus_matpower_DA("./data";
    # using a modified (mod) file that reduced load for feasibility in DC-OPF
    forecasts_pointers_file=joinpath("./data", "forecasts", "timeseries_pointers_da_7day_mod.json"),
    add_reserves=false
)

# duplicate system and prepare times series for the time varying paramters (loads, renewables, ...)
sys_uc, sys_ed = prep_systems_UCED(base_system)

# generic market formulation templates with defined network formulation
# CopperPlate-OPF: network=CopperPlatePowerModel
# DC-OPF: network=DCPPowerModel
# NFA-OPF (only line limit constraints): network=NFAPowerModel
# DC-PTDF-OPF (what ISOs do): network=PTDFPowerModel
template_uc = template_unit_commitment(; network=DCPPowerModel)
template_ed = template_economic_dispatch(; network=DCPPowerModel)

# TODO: add the following to a utility functiuon in GridAnalysis:
# for each formulation you will need to save different dual variables:
constraint_duals = if template_ed.transmission == CopperPlatePowerModel
    [:CopperPlateBalance]
elseif template_ed.transmission in [NFAPowerModel; DCPPowerModel]
    [:nodal_balance_active__Bus]
elseif template_ed.transmission == StandardPTDFModel
    [:CopperPlateBalance, :network_flow__Line]
end

# build a market clearing simulator (run `@doc UCED` for more information)
market_simulator = UCED(
    system_uc=sys_uc,
    system_ed=sys_ed,
    template_uc=template_uc,
    template_ed=template_ed,
    solver_uc=solver_uc,
    solver_ed=solver_ed,
)

#' Simulate market
# build and run simulation
results = run_multiday_simulation(
    market_simulator,
    Date("2020-01-01"), # initial time for simulation
    1; # number of steps in simulation (normally number of days to simulate)
    services_slack_variables=false,
    balance_slack_variables=true, # true because of a bug in PSI, but it wont affect prices
    constraint_duals=constraint_duals,
    name= "test_case_5bus",
    simulation_folder="./results",
);

# separate results
uc_results = get_problem_results(results, "UC");
ed_results = get_problem_results(results, "ED");

# TODO: add the following to a utility functiuon in GridAnalysis (and make it better):
# calculate prices (it will be a bit more complicated for PTDFPowerModel)
prices = if template_ed.transmission == CopperPlatePowerModel
    duals = read_dual(ed_results, :CopperPlateBalance)
    DataFrame(DateTime=collect(keys(duals)), lmp=[duals[i][1,2] / get_base_power(base_system) for i in keys(duals)])
elseif template_ed.transmission in [NFAPowerModel; DCPPowerModel]
    duals = read_dual(ed_results, :nodal_balance_active__Bus)
    # we only want the first hour of the ED prices (stored in the first row here)
    # since the second hour is here given the bug in PSI (see `@doc prep_systems_UCED`)
    df = vcat([DataFrame(duals[i][1,:]) for i in keys(duals)]...)
    # in this case, the lmps are the negative of the duals
    df[:, 2:end] = df[:, 2:end] ./ -get_base_power(base_system)
    df
elseif template_ed.transmission == StandardPTDFModel
    # MEC
    duals_MEC = read_dual(ed_results, :CopperPlateBalance)
    MEC = DataFrame(DateTime=collect(keys(duals_MEC)), MEC=[duals_MEC[i][1,2] / get_base_power(base_system) for i in keys(duals_MEC)])
    
    # MCC
    PTDF_matrix = market_simulator.kwargs[:PTDF]
    # get line limit duals
    # ignoring transformer, because it seems that PSI is also ignoring
    flow_duals = read_dual(ed_results, :network_flow__Line)
    # we only want the first hour of the ED prices (stored in the first row here)
    # since the second hour is here given the bug in PSI (see `@doc prep_systems_UCED`)
    flow_duals = vcat([DataFrame(flow_duals[i][1,:]) for i in keys(flow_duals)]...)
    # get duals in a matrix form ordered by the line names available in the PTDF
    line_names = intersect(PTDF_matrix.axes[1], names(flow_duals[:,2:end]))
    μ = Matrix(flow_duals[:, line_names])
    # calculate the congestion factor of the prices by multiplying the line duals by the PTDF
    buses = get_components(Bus, sys_ed)
    congestion_lmp = Dict()
    for bus in buses
        congestion_lmp[get_name(bus)] = - μ * [PTDF_matrix[line, get_number(bus)] for line in line_names]
    end
    congestion_lmp["DateTime"] = collect(keys(duals_MEC))
    MCC = DataFrame(congestion_lmp)

    # Calculate the final locational marginal prices: LMP (MEC + MCC)
    LMP = deepcopy(MCC)
    for row in eachrow(LMP)
        for name in names(row[2:end])
            row[name] /= get_base_power(base_system)
            row[name] += .+ MEC[MEC[:,"DateTime"] .== row["DateTime"], :MEC][1]
        end
    end
    LMP
end
