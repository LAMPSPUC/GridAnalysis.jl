# Make sure to run this file while in the examples RTS enviroment.
# '] activate ./examples/RTS'
using Dates
using DataFrames
using GridAnalysis
using Gurobi
using Measures
using Plots
using PowerSystems
using PowerSimulations
using Test

const PSY = PowerSystems

# set directory
rts_dir = download("https://github.com/GridMod/RTS-GMLC", "master", mktempdir())
# Or clone the directory and open as:
# for example: rts_dir = "/home/rafaela/Documents/PUC/LAMPS/github/RTS-GMLC"
#rts_dir = "/home/rafaela/Documents/PUC/LAMPS/github/RTS-GMLC"
rts_src_dir = joinpath(rts_dir, "RTS_Data", "SourceData")
rts_siip_dir = joinpath(rts_dir, "RTS_Data", "FormattedData", "SIIP");

# might not work if running lines manually 
# (solution: edit to be the path for this examples directory 
# for example: 'example_dir = "./examples/RTS/"')
example_dir = dirname(@__FILE__)

data_dir = joinpath(example_dir, "data")

include(joinpath(example_dir, "utils.jl")) # case utilities
include(joinpath(example_dir, "modify_RTS.jl")) # functions that modify the RTS problem

#' Data Prep and Build Market Simulator
# define solvers for Unit Commitment (UC) and Economic Dispatch (ED)
solver_uc = optimizer_with_attributes(Gurobi.Optimizer)#(Cbc.Optimizer, "logLevel" => 1, "ratioGap" => 0.5)
solver_ed = optimizer_with_attributes(Gurobi.Optimizer)#(GLPK.Optimizer)

# define systems with resolutions
sys_DA, sys_rt = get_rts_sys(rts_src_dir, rts_siip_dir;)

# Add single generator at a defined bus
node = "Bach" # define bus
gen = add_generator!(sys_DA, node, (min=0.0, max=0.0))
@test gen in get_components(Generator, sys_DA)

# create and set variable cost time-series for the generator
bidding_period = collect(1:24)
ts_array = create_generator_bids(;
    initial_bidding_time=DateTime("2020-09-01"),
    bidding_periods=bidding_period,
    system=sys_DA,
    costs=zeros(length(bidding_period)),
)
set_variable_cost!(sys_DA, gen, ts_array)

#Define range quota
range_quota = Float64.(collect(0:1:10));

# duplicate system and prepare times series for the time varying parameters (loads, renewables, ...)
sys_uc, sys_ed = prep_systems_UCED(sys_DA)

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

#Calculates the dispatch result for a bid curve
name_generator = get_name(gen);
initial_time = Date("2020-09-01");
steps = 1;
simulation_folder = mktempdir();
#simulation_folder = joinpath(example_dir, "results", "virtual_1bus", "reserve_false");
lmps_df, results_df = pq_curves_virtuals!(
    market_simulator, name_generator, range_quota, initial_time, steps, simulation_folder
)

@test isa(results_df[range_quota[1]], Dict{String,SimulationResults})
@test isa(lmps_df[range_quota[1]], Dict{String,DataFrame})

#Select data to plot
generator_name = "Bach_virtual_supply"
period = [5] #bidding_period #[5,19]
bus_name = [get_name(i) for i in get_components(Bus,sys_ed)]

# Plots
plot_price_curves(lmps_df, period, bus_name, node, initial_time, sys_ed, false)
plot_revenue_curves(
    market_simulator, lmps_df, results_df, period, generator_name, initial_time, sys_ed, false
)
plot_generation_curves(
    market_simulator, lmps_df, results_df, period, generator_name, initial_time, sys_ed
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
