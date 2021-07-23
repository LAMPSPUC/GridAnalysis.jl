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
    add_reserves=true,
)

# create and set variable cost time-series for the load
bidding_period = collect(1:24)
ts_array = create_demand_series(;
    initial_bidding_time=DateTime("2020-01-01"),
    bidding_periods=bidding_period,
    system=base_da_system,
    demands=ones(length(bidding_period)),
)

# Add single load at a defined bus
node_load = "bus3" # define bus
load = add_load!(base_da_system, node_load, 1.0)
@test load in get_components(PowerLoad, base_da_system)

add_time_series!(base_da_system, load, ts_array)

# Add single generator at a defined bus
node_gen = "bus5" # define bus
gen = add_generator!(base_da_system, node_gen, (min=0.0, max=0.0))
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

# Define range quota
range_quota = Float64.(collect(0:0.25:4));

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
name_load = get_name(load);
name_gen = get_name(gen);
initial_time = Date("2020-01-01");
steps = 1;
simulation_folder = joinpath(example_dir, "results", "Net_DCPP_Load_bus3_Offer_bus5_period_1-24")
#simulation_folder = mktempdir();
lmps_df, results_df = pq_curves_load_gen_virtuals!(
    market_simulator, name_load, range_quota, initial_time, name_gen, range_quota, steps, simulation_folder
)

@test isa(results_df[[range_quota[1] range_quota[1]]], Dict{String,SimulationResults})
@test isa(lmps_df[[range_quota[1] range_quota[1]]], Dict{String,DataFrame})

#Load the simulation done previously 
lmps_df, results_df = load_mix_pq_curves(market_simulator, range_quota, range_quota, simulation_folder)

@test isa(results_df[[range_quota[1] range_quota[1]]], Dict{String,SimulationResults})
@test isa(lmps_df[[range_quota[1] range_quota[1]]], Dict{String,DataFrame})

# Plots

plt,h,sum_deficit=heat_map_deficit(
    sys_rt,
    results_df,
    range_quota,
    range_quota,
    initial_time,
    bidding_period,
)

h,plt=heat_map_coal_generation(
    sys_uc,
    results_df,
    range_quota,
    range_quota,
    initial_time,
    sum_deficit,
    bidding_period,
    [:P__ThermalStandard, :P__RenewableDispatch],
)

h,plt=heat_map_revenue_curves_mix(
    market_simulator,
    lmps_df,
    results_df,
    bidding_period,
    range_quota,
    range_quota,
    initial_time,
    load,
    name_gen,
    "SolarBusC",
    sys_uc,
)




rt_results = get_problem_results(results_df[[4.0 4.0]]["RT"],"RT")
uc_results = get_problem_results(results_df[[4.0 4.0]]["DA"],"UC")
plt=plot_generation_stack(sys_rt, rt_results; xtickfontsize=8, margin=8mm, size=(800, 600))
plt=plot_thermal_commit_generator_stack(sys_uc, uc_results)

savefig(plt,"C:/Users/Daniele/Desktop/revenue-solar.png")
