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
# Or clone the directory and open as:
# for example: rts_dir = "/home/rafaela/Documents/PUC/LAMPS/github/RTS-GMLC"
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
sys_uc, sys_ed = prep_systems_UCED(sys_DA)

bus_names_thermal = [get_name(get_bus(i)) for i in get_components(ThermalStandard,sys_rt)]
bus_names_renewable = [get_name(get_bus(i)) for i in get_components(RenewableDispatch,sys_rt)]
bus_names = vcat(bus_names_thermal, bus_names_renewable)
bus_names = [sort(unique(bus_names))]

# generic market formulation templates with defined network formulation
# CopperPlate-OPF: network=CopperPlatePowerModel
# DC-OPF: network=DCPPowerModel
# NFA-OPF (only line limit constraints): network=NFAPowerModel
# DC-PTDF-OPF (what ISOs do): network=StandardPTDFModel
template_uc = template_unit_commitment(; network=DCPPowerModel)
template_rt = template_economic_dispatch(; network=DCPPowerModel)
template_ed = template_economic_dispatch(; network=DCPPowerModel)

reserves = get_components(VariableReserve{ReserveUp}, sys_uc)
list_reserve = [i for i in get_components(VariableReserve{ReserveUp}, sys_uc)]
list_reserve_original = deepcopy(list_reserve)

#= 

This code is the loop used to find the cases that make this problem feasible.
It was found that multiplying for two the ploblem remained infeasible but
when multiplying for three the ploblem would now be feasible.
Because of that, it was decided to find the smallest factor that belongs to that interval
moving with steps of size of 0.1.
It was found that, for this case, the problem is feasible when k is 2.5, 2.7, 2.9 or 3.0.

list = []
for k in 2:0.1:3
    for r in 1:length(reserves)
        list_reserve[r].requirement = list_reserve_original[r].requirement*k
    end
    try 
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
            balance_slack_variables=false,
            constraint_duals=constraint_duals,
            name="test_case_5bus",
            simulation_folder=mktempdir(), # Locally can use: joinpath(example_dir, "results"),
        );
        println(k)
        append!(list, k)
    catch
        println("Deficit")
    end
end

=#

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

@test isa(results, Dict{String,SimulationResults})

# separate results
uc_results = get_problem_results(results["DA"], "UC");
ed_results = get_problem_results(results["DA"], "ED");

# calculate prices
prices = evaluate_prices(market_simulator, results)

@test isa(prices, Dict{String,DataFrame})

names(prices["DA"][!, :])

# Plots
plot_generation_stack(sys_DA, ed_results; xtickfontsize=8, margin=8mm, size=(800, 600))
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
    market_simulator,
    results;
    bus_names=bus_names[1],
    xtickfontsize=8,
    size=(800, 600),
)

plot_thermal_commit_generator_stack(sys_DA, uc_results; xtickfontsize=8, size=(800, 600))
plot_thermal_commit_type_stack(sys_DA, uc_results; xtickfontsize=8, size=(800, 600))

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

# UCRT

#= 

To make the problem feasible, make sure to run the code commented.
It was tried to multiply the reserves for a factor between 2 and 3,
adding 0.1 each try.
It was found that, for this case, the problem is feasible when k is 2.5, 2.7, 2.9 or 3.0.

for r in 1:length(reserves)
    list_reserve[r].requirement = list_reserve_original[r].requirement*k
end

=#

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

@test isa(results, Dict{String,SimulationResults})

# separate results
uc_results = get_problem_results(results["RT"], "UC");
rt_results = get_problem_results(results["RT"], "RT");

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
    bus_names=bus_names[1],
    xtickfontsize=8,
    size=(800, 600),
)

plot_thermal_commit_generator_stack(sys_DA, uc_results; xtickfontsize=8, size=(800, 600))
plot_thermal_commit_type_stack(sys_DA, uc_results; xtickfontsize=8, size=(800, 600))

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
) # TO-DO: make it plot the whole day

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
) # TO-DO: make it plot the whole day

plot_prices_RT_hour(market_simulator, results, sys_rt, (-10, 50))

# UCEDRT

#= 

To make the problem feasible, make sure to run the code commented.
It was tried to multiply the reserves for a factor between 2 and 3,
adding 0.1 each try.
It was found that, for this case, the problem is feasible when k is 2.5, 2.7, 2.9 or 3.0.

for r in 1:length(reserves)
    list_reserve[r].requirement = list_reserve_original[r].requirement*k
end

=#

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
uc_results = get_problem_results(results["DA"], "UC");
rt_results = get_problem_results(results["RT"], "RT");
ed_results = get_problem_results(results["DA"], "ED");

@test isa(rt_results, PowerSimulations.SimulationProblemResults)

# calculate prices
prices = evaluate_prices_UCEDRT(market_simulator, results)

@test isa(prices, Dict{String,DataFrame})

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
plot_generation_stack(
    sys_DA,
    ed_results;
    generator_fields=[:P__RenewableDispatch],
    xtickfontsize=8,
    margin=8mm,
    size=(800, 600),
)

plot_prices(market_simulator, results; xtickfontsize=8, size=(800, 600), type="DA")
plot_prices(market_simulator, results; xtickfontsize=8, size=(800, 600), type="RT")
plot_prices(
    market_simulator,
    results;
    bus_names=["Calvin", "Beethoven", "Anna", "Cole", "Curie"],
    xtickfontsize=8,
    size=(800, 600),
)

plot_thermal_commit_generator_stack(sys_DA, uc_results; xtickfontsize=8, size=(800, 600))
plot_thermal_commit_type_stack(sys_DA, uc_results; xtickfontsize=8, size=(800, 600))

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
) # TO-DO: make it plot the whole day

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
) # TO-DO: make it plot the whole day

plot_prices_RT_hour(market_simulator, results, sys_rt, (-10, 50))
plot_DA_RT(market_simulator, results, sys_rt, (-10, 50))
