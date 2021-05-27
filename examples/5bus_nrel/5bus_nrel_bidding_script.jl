# Make sure to run this file while in the examples 5bus_nrel enviroment.
# '] activate ./examples/5bus_nrel'
using Cbc
using Dates
using DataFrames
using GLPK
using GridAnalysis
using PowerSystems
using PowerSimulations
using Test
using Measures
using Plots
using Gurobi

# might not work if running lines manually 
# (solution: edit to be the path for this examples directory 
# for example: 'example_dir = "./examples/5bus_nrel/"')
example_dir = dirname(@__FILE__)

data_dir = joinpath(example_dir, "data")

include(joinpath(example_dir, "utils.jl")) # case utilities

#' Data Prep and Build Market Simulator
# define solvers for Unit Commitment (UC) and Economic Dispatch (ED)
solver_uc = optimizer_with_attributes(Gurobi.Optimizer)#(Cbc.Optimizer, "logLevel" => 1, "ratioGap" => 0.5)
solver_ed = optimizer_with_attributes(Gurobi.Optimizer)#(GLPK.Optimizer)

# call our data preparation to build base system
# the case was modified to not have hydros nor transformers
base_system = build_5_bus_matpower_DA(
    data_dir;
    # using a modified (mod) file that reduced load for feasibility in DC-OPF
    forecasts_pointers_file=joinpath(
        data_dir, "forecasts", "timeseries_pointers_da_7day_mod.json"
    ),
    add_reserves=false,
)

# Add single generator
node = "bus1"
active_power_limits = (min=0.0, max=1.0)
ts_array = create_generator_bids(; initial_bidding_time=DateTime("2020-01-01"), bidding_periods=[1,5,16], system=base_system, costs=[10.0,15.0,16.0])

gen = add_gerator!(base_system, node, active_power_limits)

set_variable_cost!(base_system, gen, ts_array)

sys_uc, sys_ed = prep_systems_UCED(base_system)
