var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = GridAnalysis","category":"page"},{"location":"#GridAnalysis","page":"Home","title":"GridAnalysis","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for GridAnalysis.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [GridAnalysis]","category":"page"},{"location":"#GridAnalysis.DayAhead","page":"Home","title":"GridAnalysis.DayAhead","text":"DayAhead <: MarketSimulator\n\nA type that represents DayAhead (DA) market clearings. \n\n\n\n\n\n","category":"type"},{"location":"#GridAnalysis.MarketSimulator","page":"Home","title":"GridAnalysis.MarketSimulator","text":"MarketSimulator\n\nA type to represent different market clearing structures.\n\n\n\n\n\n","category":"type"},{"location":"#GridAnalysis.UCED","page":"Home","title":"GridAnalysis.UCED","text":"UCED <: DayAhead\n\nStructure with two DA systems stored (system_uc and system_ed). One represents an  Unit Commitment (UC) used to fix binary variables, the other an Economic Dispatch (ED) defined by template template_ed. Each model needs its own solver (solver_uc and solver_ed).  system_ed and the resulting problem are defined as a heuristic to calculate DA energy prices.\n\n\n\n\n\n","category":"type"},{"location":"#GridAnalysis.plot_generation_stack-Tuple","page":"Home","title":"GridAnalysis.plot_generation_stack","text":"plot_generation_stack(\n    system::System,\n    results::SimulationProblemResults;\n    generator_fields::AbstractArray=[:P__ThermalStandard, :P__RenewableDispatch],\n    bus_names::AbstractArray=[])\n\nPlot the generation mix over the time period covered by the results. The bus_names and generator_fields control which buses, and generator types we want to include the plot.\n\n\n\n\n\n","category":"method"},{"location":"#PowerSimulations.Simulation","page":"Home","title":"PowerSimulations.Simulation","text":"Simulation(simulator::UCED, initial_time::Date, steps::Int) -> PSI.Simulation\n\nBuilds a multiday UC-ED PSI.Simulation from a UCED market simulator, an initial date and number of simulation steps.\n\n\n\n\n\n","category":"type"},{"location":"#GridAnalysis.MEC-Tuple{PowerSystems.System, PowerSimulations.SimulationProblemResults}","page":"Home","title":"GridAnalysis.MEC","text":"MEC(system::System, problem_results::PSI.SimulationProblemResults)\n\nReturns the Marginal Energy Component (MEC) of OPF problem's solution. This component is defined as the dual value (or lagrangian multiplier) of a grid-wide energy balance constraint (named in PSI CopperPlateBalance)\n\n\n\n\n\n","category":"method"},{"location":"#GridAnalysis.duals_constraint_names-Tuple{Type{PowerSimulations.CopperPlatePowerModel}}","page":"Home","title":"GridAnalysis.duals_constraint_names","text":"duals_constraint_names(<:AbstractPowerModel)\n\nReturn the constraints for which we care about the duals (because they form the energy prices) when using a specified network formulation.\n\n\n\n\n\n","category":"method"},{"location":"#GridAnalysis.evaluate_prices-Tuple{Type{PowerSimulations.CopperPlatePowerModel}, PowerSystems.System, PowerSimulations.SimulationProblemResults, Any}","page":"Home","title":"GridAnalysis.evaluate_prices","text":"evaluate_prices(::Type{CopperPlatePowerModel}, system::System, problem_results::PSI.SimulationProblemResults, kwargs)\n\nReturns a grid-unique series of energy prices for the simulation's data-range.  In this formulation, for each period of the problem, the energy prices are constructed using the dual value (or lagrangian multiplier) of a grid-wide energy balance constraint (named in PSI CopperPlateBalance). \n\n\n\n\n\n","category":"method"},{"location":"#GridAnalysis.evaluate_prices-Tuple{Type{PowerSimulations.StandardPTDFModel}, PowerSystems.System, PowerSimulations.SimulationProblemResults, Any}","page":"Home","title":"GridAnalysis.evaluate_prices","text":"evaluate_prices(transmission::Type{StandardPTDFModel}, system::System, problem_results::PSI.SimulationProblemResults, kwargs)\n\nReturns a nodal-wide series of energy prices (locational marginal prices - LMPS) for the simulation's data-range.  In this formulation, for each period of the problem, the energy prices are constructed using the dual value (or lagrangian multiplier) of a grid-wide energy balance constraint (named in PSI CopperPlateBalance) and  the dual values of a branch-wide power-flow constraint (named in PSI network_flow__Line).\n\n\n\n\n\n","category":"method"},{"location":"#GridAnalysis.evaluate_prices-Tuple{UCED, PowerSimulations.SimulationProblemResults}","page":"Home","title":"GridAnalysis.evaluate_prices","text":"evaluate_prices(market_simulator::UCED, problem_results::PSI.SimulationProblemResults)\n\nReturns energy prices for the simulation's data-range.  \n\n\n\n\n\n","category":"method"},{"location":"#GridAnalysis.evaluate_prices-Tuple{Union{Type{PowerModels.DCPPowerModel}, Type{PowerModels.NFAPowerModel}}, PowerSystems.System, PowerSimulations.SimulationProblemResults, Any}","page":"Home","title":"GridAnalysis.evaluate_prices","text":"evaluate_prices(::Union{Type{NFAPowerModel},Type{DCPPowerModel}}, system::System, problem_results::PSI.SimulationProblemResults, kwargs)\n\nReturns a nodal-wide series of energy prices (locational marginal prices - LMPS) for the simulation's data-range.  In this formulation, for each period of the problem, the energy prices are constructed using the dual values (or lagrangian multipliers) of each nodal energy balance constraint (named in PSI nodal_balance_active__Bus). \n\n\n\n\n\n","category":"method"},{"location":"#GridAnalysis.fuel_type_mapping-Tuple{PowerSystems.System}","page":"Home","title":"GridAnalysis.fuel_type_mapping","text":"fuel_type_mapping(system)\n\nreturn a mapping between the bus name and the fuel type for the given `system`.\n\n\n\n\n\n","category":"method"},{"location":"#GridAnalysis.get_time_series_params-Tuple{PowerSystems.System}","page":"Home","title":"GridAnalysis.get_time_series_params","text":"get_time_series_params(system::System)\n\nReturns the parameters associated with the time-series attached to the system.\n\n\n\n\n\n","category":"method"},{"location":"#GridAnalysis.plot_generation_stack!-Tuple","page":"Home","title":"GridAnalysis.plot_generation_stack!","text":"plot_generation_stack(\n    system::System,\n    results::SimulationProblemResults;\n    generator_fields::AbstractArray=[:P__ThermalStandard, :P__RenewableDispatch],\n    bus_names::AbstractArray=[])\n\nPlot the generation mix over the time period covered by the results. The bus_names and generator_fields control which buses, and generator types we want to include the plot.\n\n\n\n\n\n","category":"method"},{"location":"#GridAnalysis.plot_generation_stack!-Tuple{RecipesBase.AbstractPlot, Vararg{Any, N} where N}","page":"Home","title":"GridAnalysis.plot_generation_stack!","text":"plot_generation_stack(\n    system::System,\n    results::SimulationProblemResults;\n    generator_fields::AbstractArray=[:P__ThermalStandard, :P__RenewableDispatch],\n    bus_names::AbstractArray=[])\n\nPlot the generation mix over the time period covered by the results. The bus_names and generator_fields control which buses, and generator types we want to include the plot.\n\n\n\n\n\n","category":"method"},{"location":"#GridAnalysis.run_multiday_simulation","page":"Home","title":"GridAnalysis.run_multiday_simulation","text":"run_multiday_simulation(simulator::UCED, initial_time::Date, steps::Int) -> SimulationResults\n\nRuns a multiday PSI.Simulation from a MarketSimulator, an initial date and number of simulation steps.\n\n\n\n\n\n","category":"function"}]
}
