# Make sure to run this file while in the examples 5bus_nrel enviroment.
# '] activate ./examples/5bus_nrel'
using Dates
using DataFrames
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

# Add single generator at a defined bus
node = "bus5" # define bus
gen = add_gerator!(base_system, node, (min = 0.0, max = 0.0))

# create and set variable cost time-series for the generator
ts_array = create_generator_bids(; initial_bidding_time=DateTime("2020-01-01"), bidding_periods=collect(1:24), system=base_system, costs=ones(24).*0)
set_variable_cost!(base_system, gen, ts_array)

#Define range quota
range_quota=collect(0:1:5)

# Virtual Bids Simulation 
initial_time = Date("2020-01-01")
lmps_df, results_df = pq_curves_da(
    gen;
    base_system,
    range_quota,
    initial_time, #: TODO: The same as ts_array
    steps=1,
    simulation_folder = joinpath(example_dir, "results"),
    solver_uc,
    solver_ed
) #:TODO: Plots


max_gen=4
bus=get_name(get_bus(gen))
variable_results = read_realized_variables(results_df[max_gen], names=[:P__ThermalStandard])
generator_data = getindex.(Ref(variable_results), [:P__ThermalStandard])
lmps=lmps_df[max_gen]
virtual_gen=generator_data[1][!,:7]
p=lmps[!,:"bus5"] # TODO: Put bus here
revenue=p.*virtual_gen

#=
# Plots
plot_generation_stack(base_system, ed_results; xtickfontsize=8, margin=8mm, size=(800, 600))
plot_generation_stack(
    base_system,
    ed_results;
    bus_names=["bus5"],
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

plot_prices(market_simulator, ed_results; xtickfontsize=8, size=(800, 600))
plot_prices(
    market_simulator,
    ed_results;
    bus_names=["bus1", "bus3"],
    xtickfontsize=8,
    size=(800, 600),
)

plot_thermal_commit(base_system, uc_results; xtickfontsize=8, size=(800, 600))
plot_thermal_commit(
    base_system, uc_results; bus_names=["bus1", "bus2"], xtickfontsize=8, size=(800, 600)
)

plot_demand_stack(sys_uc, uc_results; xtickfontsize=8, size=(800, 600))
plot_demand_stack(
    sys_uc, uc_results; bus_names=["bus2", "bus3"], xtickfontsize=8, size=(800, 600)
)

plot_net_demand_stack_prev(sys_uc, uc_results; xtickfontsize=8, size=(800, 600))
plot_net_demand_stack_prev(
    sys_uc, uc_results; bus_names=["bus2", "bus3"], xtickfontsize=8, size=(800, 600)
)

plot_net_demand_stack(sys_uc, uc_results; xtickfontsize=8, size=(800, 600))
plot_net_demand_stack(
    sys_uc, uc_results; bus_names=["bus2", "bus3"], xtickfontsize=8, size=(800, 600)
)
=#