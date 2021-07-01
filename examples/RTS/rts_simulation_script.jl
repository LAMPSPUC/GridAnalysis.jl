# Make sure to run this file while in the examples RTS enviroment.
# '] activate ./examples/RTS'
using Cbc
using Dates
using DataFrames
using GLPK
using GridAnalysis
using Gurobi
using PowerSystems
using PowerSimulations
using Test
using Measures
using Plots

const PSY = PowerSystems

# set directory
rts_dir = download("https://github.com/GridMod/RTS-GMLC", "master", mktempdir())
#rts_dir = "/home/rafaela/Documents/PUC/LAMPS/github/RTS-GMLC"
rts_src_dir = joinpath(rts_dir, "RTS_Data", "SourceData")
rts_siip_dir = joinpath(rts_dir, "RTS_Data", "FormattedData", "SIIP");

# might not work if running lines manually 
# (solution: edit to be the path for this examples directory 
# for example: 'example_dir = "./examples/RTS/"')
example_dir = dirname(@__FILE__)

include(joinpath(example_dir, "utils.jl")) # case utilities
include(joinpath(example_dir, "modify_RTS.jl")) # functions that modify the RTS problem

#' Data Prep and Build Market Simulator
# define solvers for Unit Commitment (UC), Real Time (RT) and Economic Dispatch (ED)
solver_uc = optimizer_with_attributes(Gurobi.Optimizer)
solver_rt = optimizer_with_attributes(Gurobi.Optimizer)
solver_ed = optimizer_with_attributes(Gurobi.Optimizer)

# define systems with resolutions
sys_DA, sys_rt = get_rts_sys(rts_src_dir, rts_siip_dir;)
transform_single_time_series!(sys_rt, 12, Minute(15))
sys_uc, sys_ed = prep_systems_UCED(sys_DA)

# generic market formulation templates with defined network formulation
# CopperPlate-OPF: network=CopperPlatePowerModel
# DC-OPF: network=DCPPowerModel
# NFA-OPF (only line limit constraints): network=NFAPowerModel
# DC-PTDF-OPF (what ISOs do): network=StandardPTDFModel
template_uc = template_unit_commitment(; network=DCPPowerModel)
template_rt = template_economic_dispatch(; network=DCPPowerModel)
template_ed = template_economic_dispatch(; network=DCPPowerModel)

# UC-ED

# build a market clearing simulator (run `@doc UCED` for more information)
market_simulator = UCED(;
    system_uc=sys_uc,
    system_ed=sys_ed,
    template_uc=template_uc,
    template_ed=template_ed,
    solver_uc=solver_uc,
    solver_ed=solver_ed,
)

@test isa(market_simulator, UCED)

# for each formulation you will need to save different dual variables:
constraint_duals = duals_constraint_names(market_simulator)

@test isa(constraint_duals, AbstractVector{Symbol})

# Simulate market
# build and run simulation
results = run_multiday_simulation(
    market_simulator,
    Date("2020-09-01"), # initial time for simulation
    1; # number of steps in simulation (normally number of days to simulate)
    services_slack_variables=false,
    balance_slack_variables=false,
    constraint_duals=constraint_duals,
    name="test_case_5bus",
    simulation_folder=mktempdir(), # Locally can use: joinpath(example_dir, "results"),
);

@test isa(results, SimulationResults)

# separate results
uc_results = get_problem_results(results, "UC");
ed_results = get_problem_results(results, "ED");

# calculate prices
prices = evaluate_prices(market_simulator, results)

@test isa(prices, Dict{String,DataFrame})

names(prices["ED"][!, :])

# get the generator data
generator_metadata = [gen for gen in get_components(Generator, sys_ed)]

# get the fuel type for all buses
fuel_type_dict = fuel_type_mapping(sys_ed)

# get all names of the fuel types and make sure its unique
fuel_names = unique(keys(fuel_type_dict))

# create the vectors, one for each type of fuel
hydro = Vector{String}()
wind = Vector{String}()
solar = Vector{String}()
coal = Vector{String}()
natural_gas = Vector{String}()
nuclear = Vector{String}()

# put the name of the bus in its respective vector
# these vectors are useful for selecting which bus is desired to plot
# for example: if is wanted to plot only the wind or only the solar generators
for i in 1:length(fuel_names)
    if occursin("HYDRO", fuel_names[i]) == true
        push!(hydro, get_bus_name(generator_metadata[i]))
    elseif occursin("WIND", fuel_names[i]) == true
        push!(wind, get_bus_name(generator_metadata[i]))
    elseif occursin("PV", fuel_names[i]) == true || occursin("CSP", fuel_names[i]) == true
        push!(solar, get_bus_name(generator_metadata[i]))
    elseif occursin("STEAM", fuel_names[i]) == true
        push!(coal, get_bus_name(generator_metadata[i]))
    elseif occursin("CT", fuel_names[i]) == true || occursin("CC", fuel_names[i]) == true
        push!(natural_gas, get_bus_name(generator_metadata[i]))
    else
        push!(nuclear, get_bus_name(generator_metadata[i]))
    end
end

# change the name of the fuel type so that it follows a pattern
for (i, key) in enumerate(fuel_names)
    if occursin("HYDRO", fuel_names[i]) == true
        fuel_type_dict[key] = "HYDRO"
    elseif occursin("WIND", fuel_names[i]) == true
        fuel_type_dict[key] = "WIND"
    elseif occursin("PV", fuel_names[i]) == true || occursin("CSP", fuel_names[i]) == true
        fuel_type_dict[key] = "SOLAR"
    elseif occursin("STEAM", fuel_names[i]) == true
        fuel_type_dict[key] = "COAL"
    elseif occursin("CT", fuel_names[i]) == true || occursin("CC", fuel_names[i]) == true
        fuel_type_dict[key] = "NATURAL GAS"
    else
        fuel_type_dict[key] = "NUCLEAR"
    end
end

# Plots
plot_generation_stack(sys_DA, ed_results; xtickfontsize=8, margin=8mm, size=(800, 600))
#plot_generation_stack(sys_DA, ed_results; bus_names=wind, xtickfontsize=8, margin=8mm, size=(800, 600))
#plot_generation_stack(sys_DA, ed_results; bus_names=coal, xtickfontsize=8, margin=8mm, size=(800, 600))
plot_generation_stack(
    sys_DA,
    ed_results;
    generator_fields=[:P__RenewableDispatch],
    xtickfontsize=8,
    margin=8mm,
    size=(800, 600),
)
plot_generation_stack(
    sys_DA,
    ed_results;
    generator_fields=[:P__ThermalStandard],
    xtickfontsize=8,
    margin=8mm,
    size=(800, 600),
)
plot_generation_stack(
    sys_DA,
    uc_results;
    generator_fields=[:P__RenewableDispatch],
    xtickfontsize=8,
    margin=8mm,
    size=(800, 600),
)

plot_prices(market_simulator, results; xtickfontsize=8, size=(800, 600))
plot_prices(
    market_simulator,
    results;
    bus_names=["Calvin", "Beethoven", "Anna", "Cole", "Curie"],
    xtickfontsize=8,
    size=(800, 600),
)
plot_prices(
    market_simulator, results; bus_names=unique(solar), xtickfontsize=8, size=(800, 600)
)
plot_prices(
    market_simulator, results; bus_names=unique(wind), xtickfontsize=8, size=(800, 600)
)
plot_prices(
    market_simulator, results; bus_names=unique(nuclear), xtickfontsize=8, size=(800, 600)
)

plot_thermal_commit(sys_DA, uc_results; xtickfontsize=8, size=(800, 600))

plot_demand_stack(sys_uc, uc_results; xtickfontsize=8, size=(800, 600))
plot_demand_stack(sys_ed, ed_results; xtickfontsize=8, size=(800, 600))
plot_demand_stack(
    sys_uc,
    uc_results;
    xtickfontsize=8,
    size=(800, 600),
    type="Deterministic",
    start_time=DateTime("2020-09-01"),
)

plot_net_demand_stack(sys_uc, uc_results; xtickfontsize=8, size=(800, 600))
plot_net_demand_stack(sys_ed, ed_results; xtickfontsize=8, size=(800, 600))
plot_net_demand_stack(
    sys_uc,
    uc_results;
    xtickfontsize=8,
    size=(800, 600),
    type="Deterministic",
    start_time=DateTime("2020-09-01"),
)
#plot_net_demand_stack(sys_uc, uc_results; bus_names = ["101_CT_1", "101_CT_2", "102_CT_1", "102_CT_2"], xtickfontsize=8, size=(800, 600))
#plot_net_demand_stack(sys_uc, uc_results; bus_names = ["Calvin", "Beethoven", "Anna", "Cole"], xtickfontsize=8, size=(800, 600))

# UCRT

# build a market clearing simulator (run `@doc UCED` for more information)
market_simulator = UCRT(;
    system_uc=sys_uc,
    system_rt=sys_rt,
    template_uc=template_uc,
    template_rt=template_rt,
    solver_uc=solver_uc,
    solver_rt=solver_rt,
)

@test isa(market_simulator, UCRT)

# for each formulation you will need to save different dual variables:
constraint_duals = duals_constraint_names(market_simulator)

@test isa(constraint_duals, AbstractVector{Symbol})

# Simulate market
# build and run simulation
results = run_multiday_simulation(
    market_simulator,
    Date("2020-09-01"), # initial time for simulation
    1; # number of steps in simulation (normally number of days to simulate)
    services_slack_variables=false,
    balance_slack_variables=true,
    constraint_duals=constraint_duals,
    name="test_case_5bus",
    simulation_folder=mktempdir(), # Locally can use: joinpath(example_dir, "results"),
);

@test isa(results, SimulationResults)

# separate results
uc_results = get_problem_results(results, "UC");
rt_results = get_problem_results(results, "RT");

# calculate prices
prices = evaluate_prices(market_simulator, results)

@test isa(prices, Dict{String,DataFrame})

names(prices["RT"][!, :])

# Plots
plot_generation_stack(sys_rt, rt_results; xtickfontsize=8, margin=8mm, size=(800, 600))
plot_generation_stack(
    sys_rt,
    rt_results;
    generator_fields=[:P__RenewableDispatch],
    xtickfontsize=8,
    margin=8mm,
    size=(800, 600),
)
plot_generation_stack(
    sys_rt,
    rt_results;
    generator_fields=[:P__ThermalStandard],
    xtickfontsize=8,
    margin=8mm,
    size=(800, 600),
)
plot_generation_stack(
    sys_DA,
    uc_results;
    generator_fields=[:P__RenewableDispatch],
    xtickfontsize=8,
    margin=8mm,
    size=(800, 600),
)

plot_prices(market_simulator, results; xtickfontsize=8, size=(800, 600))
plot_prices(
    market_simulator,
    results;
    bus_names=["Calvin", "Beethoven", "Anna", "Cole", "Curie"],
    xtickfontsize=8,
    size=(800, 600),
)

plot_thermal_commit(sys_DA, uc_results; xtickfontsize=8, size=(800, 600))

plot_demand_stack(sys_uc, uc_results; xtickfontsize=8, size=(800, 600))
# plot_demand_stack(sys_rt, rt_results; xtickfontsize=8, size=(800, 600)) # too much computer demmanding
plot_demand_stack(
    sys_uc,
    uc_results;
    xtickfontsize=8,
    size=(800, 600),
    type="Deterministic",
    start_time=DateTime("2020-09-01"),
)
plot_demand_stack(
    sys_rt,
    rt_results;
    xtickfontsize=8,
    size=(800, 600),
    type="Deterministic",
    start_time=DateTime("2020-09-01"),
)

plot_net_demand_stack(sys_uc, uc_results; xtickfontsize=8, size=(800, 600))
plot_net_demand_stack(sys_rt, rt_results; xtickfontsize=8, size=(800, 600))
plot_net_demand_stack(
    sys_uc,
    uc_results;
    xtickfontsize=8,
    size=(800, 600),
    type="Deterministic",
    start_time=DateTime("2020-09-01"),
)
plot_net_demand_stack(
    sys_rt,
    rt_results;
    xtickfontsize=8,
    size=(800, 600),
    type="Deterministic",
    start_time=DateTime("2020-09-01"),
)

plot_prices_RT_hour(prices, (-5, 20))

# UCEDRT

# build a market clearing simulator
market_simulator = UCEDRT(;
    system_uc=sys_uc,
    system_rt=sys_rt,
    system_ed=sys_ed,
    template_uc=template_uc,
    template_rt=template_rt,
    template_ed=template_ed,
    solver_uc=solver_uc,
    solver_rt=solver_rt,
    solver_ed=solver_ed,
)

@test isa(market_simulator, UCEDRT)

# for each formulation you will need to save different dual variables:
constraint_duals = duals_constraint_names(market_simulator)

@test isa(constraint_duals[1], AbstractVector{Symbol})
@test isa(constraint_duals[2], AbstractVector{Symbol})

# Simulate market
# build and run simulation
results = run_multiday_simulation(
    market_simulator,
    Date("2020-09-01"), # initial time for simulation
    1; # number of steps in simulation (normally number of days to simulate)
    services_slack_variables=false,
    balance_slack_variables=true,
    constraint_duals=constraint_duals,
    name="test_case_5bus",
    simulation_folder=mktempdir(), # Locally can use: joinpath(example_dir, "results"),
);

@test isa(results, Dict{String,SimulationResults})

# separate results
uc_results = get_problem_results(results["ED"], "UC");
rt_results = get_problem_results(results["RT"], "RT");
ed_results = get_problem_results(results["ED"], "ED");

@test isa(rt_results, PowerSimulations.SimulationProblemResults)

# calculate prices
prices = evaluate_prices_UCEDRT(market_simulator, results)

@test isa(prices, Dict{String,DataFrame})

# Plots
plot_generation_stack(sys_rt, rt_results; xtickfontsize=8, margin=8mm, size=(800, 600))
#plot_generation_stack(sys_rt, rt_results; bus_names=["101_CT_1", "101_CT_2", "102_CT_1", "102_CT_2"], xtickfontsize=8, margin=8mm, size=(800, 600))
#plot_generation_stack(sys_rt, rt_results; bus_names=["Calvin", "Beethoven", "Anna", "Cole", "Curie"], xtickfontsize=8, margin=8mm, size=(800, 600))
plot_generation_stack(
    sys_rt,
    rt_results;
    generator_fields=[:P__RenewableDispatch],
    xtickfontsize=8,
    margin=8mm,
    size=(800, 600),
)
plot_generation_stack(
    sys_rt,
    rt_results;
    generator_fields=[:P__ThermalStandard],
    xtickfontsize=8,
    margin=8mm,
    size=(800, 600),
)
plot_generation_stack(
    sys_DA,
    uc_results;
    generator_fields=[:P__RenewableDispatch],
    xtickfontsize=8,
    margin=8mm,
    size=(800, 600),
)
plot_generation_stack(
    sys_DA,
    ed_results;
    generator_fields=[:P__RenewableDispatch],
    xtickfontsize=8,
    margin=8mm,
    size=(800, 600),
)

plot_prices(market_simulator, results; xtickfontsize=8, size=(800, 600), type="ED")
plot_prices(market_simulator, results; xtickfontsize=8, size=(800, 600), type="RT")
plot_prices(
    market_simulator,
    results;
    bus_names=["Calvin", "Beethoven", "Anna", "Cole", "Curie"],
    xtickfontsize=8,
    size=(800, 600),
)

plot_thermal_commit(sys_DA, uc_results; xtickfontsize=8, size=(800, 600))
#plot_thermal_commit(sys_DA, uc_results; bus_names=["Calvin", "Beethoven", "Anna", "Cole"], xtickfontsize=8, size=(800, 600))

plot_demand_stack(sys_uc, uc_results; xtickfontsize=8, size=(800, 600))
# plot_demand_stack(sys_rt, rt_results; xtickfontsize=8, size=(800, 600)) # too much computer demmanding
plot_demand_stack(
    sys_uc, uc_results; xtickfontsize=8, size=(800, 600), type="Deterministic"
)
plot_demand_stack(
    sys_rt, rt_results; xtickfontsize=8, size=(800, 600), type="Deterministic"
)
#plot_demand_stack(sys_uc, uc_results; bus_names = ["101_CT_1", "101_CT_2", "102_CT_1", "102_CT_2"], xtickfontsize=8, size=(800, 600))
#plot_demand_stack(sys_uc, uc_results; bus_names = ["Calvin", "Beethoven", "Anna", "Cole"], xtickfontsize=8, size=(800, 600))

plot_net_demand_stack(sys_uc, uc_results; xtickfontsize=8, size=(800, 600))
plot_net_demand_stack(sys_rt, rt_results; xtickfontsize=8, size=(800, 600))
plot_net_demand_stack(
    sys_uc, uc_results; xtickfontsize=8, size=(800, 600), type="Deterministic"
)
plot_net_demand_stack(
    sys_rt, rt_results; xtickfontsize=8, size=(800, 600), type="Deterministic"
)
#plot_net_demand_stack(sys_uc, uc_results; bus_names = ["101_CT_1", "101_CT_2", "102_CT_1", "102_CT_2"], xtickfontsize=8, size=(800, 600))
#plot_net_demand_stack(sys_uc, uc_results; bus_names = ["Calvin", "Beethoven", "Anna", "Cole"], xtickfontsize=8, size=(800, 600))

plot_prices_RT_hour(prices, (-5, 20))
plot_DA_RT(prices, (-10, 50))
