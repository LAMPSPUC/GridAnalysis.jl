# Make sure to run this file while in the examples 5bus_nrel enviroment.
# '] activate ./examples/5bus_nrel'
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

# might not work if running lines manually 
# (solution: edit to be the path for this examples directory 
# for example: 'example_dir = "./examples/5bus_nrel/"')
example_dir = dirname(@__FILE__)

data_dir = joinpath(example_dir, "data")

include(joinpath(example_dir, "utils.jl")) # case utilities

#' Data Prep and Build Market Simulator
# define solvers for Unit Commitment (UC) and Economic Dispatch (ED)
solver_uc = optimizer_with_attributes(Gurobi.Optimizer)
solver_ed = optimizer_with_attributes(Gurobi.Optimizer)

# call our data preparation to build base system
# the case was modified to not have hydros nor transformers
base_system = build_5_bus_matpower_DA(
    data_dir;
    # using a modified (mod) file that reduced load for feasibility in DC-OPF
    forecasts_pointers_file=joinpath(
        data_dir, "forecasts", "timeseries_pointers_da_7day_mod.json"
    ),
    add_reserves=true,
)

# duplicate system and prepare times series for the time varying parameters (loads, renewables, ...)
sys_uc, sys_ed = prep_systems_UCED(base_system)

# generic market formulation templates with defined network formulation
# CopperPlate-OPF: network=CopperPlatePowerModel
# DC-OPF: network=DCPPowerModel
# NFA-OPF (only line limit constraints): network=NFAPowerModel
# DC-PTDF-OPF (what ISOs do): network=StandardPTDFModel
template_uc = template_unit_commitment(; network=DCPPowerModel)
template_ed = template_economic_dispatch(; network=DCPPowerModel)

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
    Date("2020-01-01"), # initial time for simulation
    1; # number of steps in simulation (normally number of days to simulate)
    services_slack_variables=false,
    balance_slack_variables=false,
    constraint_duals=constraint_duals,
    name="test_case_5bus",
    simulation_folder=mktempdir(), # Locally can use: joinpath(example_dir, "results"),
);

@test isa(results, Dict{String,SimulationResults})

# separate results
uc_results = get_problem_results(results["DA"], "UC");
ed_results = get_problem_results(results["DA"], "ED");

# calculate prices
prices = evaluate_prices(market_simulator, results)

@test isa(prices, Dict{String,DataFrame})

# Plots
plot_generation_stack(base_system, ed_results; xtickfontsize=8, margin=8mm, size=(800, 600))
plot_generation_stack(
    base_system,
    ed_results;
    bus_names=["bus1", "bus3"],
    xtickfontsize=8,
    margin=8mm,
    size=(800, 600),
)
plot_generation_stack(
    base_system,
    ed_results;
    generator_fields=[:P__RenewableDispatch],
    xtickfontsize=8,
    margin=8mm,
    size=(800, 600),
)
plot_generation_stack(
    base_system,
    ed_results;
    generator_fields=[:P__ThermalStandard],
    bus_names=["bus1", "bus3"],
    xtickfontsize=8,
    margin=8mm,
    size=(800, 600),
)
plot_generation_stack(
    base_system,
    uc_results;
    generator_fields=[:P__RenewableDispatch],
    bus_names=["bus3"],
    xtickfontsize=8,
    margin=8mm,
    size=(800, 600),
)

plot_prices(market_simulator, results; xtickfontsize=8, size=(800, 600))
plot_prices(
    market_simulator, results; bus_names=["bus1", "bus3"], xtickfontsize=8, size=(800, 600)
)

plot_thermal_commit_generator_stack(base_system, uc_results; xtickfontsize=8, size=(800, 600))
plot_thermal_commit_generator_stack(
    base_system, uc_results; bus_names=["bus1", "bus3"], xtickfontsize=8, size=(800, 600)
)

plot_demand_stack(sys_uc; xtickfontsize=8, size=(800, 600))
plot_demand_stack(sys_uc; bus_names=["bus2", "bus3"], xtickfontsize=8, size=(800, 600))

plot_net_demand_stack(sys_uc; xtickfontsize=8, size=(800, 600))
plot_net_demand_stack(sys_uc; bus_names=["bus2", "bus3"], xtickfontsize=8, size=(800, 600))
