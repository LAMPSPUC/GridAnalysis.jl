# Make sure to run this file while in the examples 5bus_nrel enviroment.
# '] activate ./examples/5bus_nrel'
using Dates
using DataFrames
using GridAnalysis
using Gurobi
using Measures
using Plots
using PowerSystems
using PowerSimulations
using Test

# might not work if running lines manually 
# (solution: edit to be the path for this examples directory 
# for example: 'example_dir = "./examples/5bus_nrel/"')
example_dir = dirname(@__FILE__)

data_dir = joinpath(example_dir, "data")

include(joinpath(example_dir, "utils.jl")) # case utilities

#' Data Prep and Build Market Simulator
# define solvers for Unit Commitment (UC), Real Time (RT) and Economic Dispatch (ED)
solver_uc = optimizer_with_attributes(Gurobi.Optimizer)
solver_rt = optimizer_with_attributes(Gurobi.Optimizer)
solver_ed = optimizer_with_attributes(Gurobi.Optimizer)

# call our data preparation to build base system
# the case was modified to not have hydros nor transformers
sys_rt = build_5_bus_matpower_RT(data_dir;)

base_da_system = build_5_bus_matpower_DA(
    data_dir;
    # using a modified (mod) file that reduced load for feasibility in DC-OPF
    forecasts_pointers_file=joinpath(
        data_dir, "forecasts", "timeseries_pointers_da_7day_mod.json"
    ),
    add_reserves=false,
)

# Add single generator at a defined bus
node = "bus1" # define bus
gen = add_generator!(base_da_system, node, (min=0.0, max=0.0))
@test gen in get_components(Generator, base_da_system)

# create and set variable cost time-series for the generator
bidding_period = collect(1:24)
ts_array = create_generator_bids(;
    initial_bidding_time=DateTime("2020-01-01"),
    bidding_periods=bidding_period,
    system=base_da_system,
    costs=zeros(length(bidding_period)),
)
set_variable_cost!(base_da_system, gen, ts_array)

#Define range quota
range_quota = Float64.(collect(0:0.1:4));

# duplicate system and prepare times series for the time varying parameters (loads, renewables, ...)
sys_uc, sys_ed = prep_systems_UCED(base_da_system)

# generic market formulation templates with defined network formulation
# CopperPlate-OPF: network=CopperPlatePowerModel
# DC-OPF: network=DCPPowerModel
# NFA-OPF (only line limit constraints): network=NFAPowerModel
# DC-PTDF-OPF (what ISOs do): network=StandardPTDFModel
template_uc = template_unit_commitment(; network=DCPPowerModel)
template_rt = template_economic_dispatch(; network=DCPPowerModel)
template_ed = template_economic_dispatch(; network=DCPPowerModel)

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

#Calculates the dispatch result for a bid curve
name_generator = get_name(gen);
initial_time = Date("2020-01-01");
steps = 1;
#simulation_folder = joinpath(example_dir, "results", "...");
simulation_folder = mktempdir()

lmps_df, results_df = pq_curves_virtuals!(
    market_simulator, name_generator, range_quota, initial_time, steps, simulation_folder
)

@test isa(results_df[range_quota[1]], Dict{String,SimulationResults})
@test isa(lmps_df[range_quota[1]], Dict{String,DataFrame})

#Load the simulation done previously 
lmps_df, results_df = load_pq_curves(market_simulator, range_quota, simulation_folder)

@test isa(results_df[range_quota[1]], Dict{String,SimulationResults})
@test isa(lmps_df[range_quota[1]], Dict{String,DataFrame})

#Select data to plot
generator_name = "bus1_virtual_supply"
period = [19]
bus_name = ["bus1", "bus2", "bus3", "bus4", "bus5"]

# Plots
plot_price_curves(lmps_df, period, bus_name, node, initial_time, sys_rt, true)
plot_revenue_curves(
    market_simulator, lmps_df, results_df, period, generator_name, initial_time, sys_rt, true
)
plot_generation_curves(
    market_simulator, lmps_df, results_df, period, generator_name, initial_time, sys_rt
)

type = "DA";
plot_generation_stack_virtual(
    sys_uc,
    results_df;
    type,
    period=period,
    initial_time,
    xtickfontsize=8,
    margin=8mm,
    size=(800, 600),
)
type = "RT";
plot_generation_stack_virtual(
    sys_rt,
    results_df;
    type,
    period=period,
    initial_time,
    xtickfontsize=8,
    margin=8mm,
    size=(800, 600),
)

plot_revenue_curves_renewable(
    market_simulator, lmps_df, results_df, [0.0, 1.0], "SolarBusC", node, false
)

plot_revenue_curves_renewable(
    market_simulator, lmps_df, results_df, [0.0, 1.0, 2.0], "WindBusA", node, false
)

plot_revenue_curves_renewable_plus_virtual(
    market_simulator,
    lmps_df,
    results_df,
    [0.0, 1.0, 2.0],
    "WindBusA",
    name_generator,
    false
)

plot_generation_curves_renewable(
    lmps_df, results_df, [0.0, 1.0, 2.0], "WindBusA", node
)

plot_thermal_commit_virtual(sys_uc, results_df; period=period, initial_time, xtickfontsize=8, margin=8mm, size=(800, 600),)