"""
    MarketSimulator

A type to represent different market clearing structures.
"""
abstract type MarketSimulator end

"""
    DayAhead <: MarketSimulator

A type that represents DayAhead (DA) market clearings. 
"""
abstract type DayAhead <: MarketSimulator end

"""
    UCED <: DayAhead

Structure with two DA systems stored (`system_uc` and `system_ed`). One represents an 
Unit Commitment (UC) used to fix binary variables, the other an Economic Dispatch (ED)
defined by template `template_ed`. Each model needs its own solver (`solver_uc` and `solver_ed`). 
`system_ed` and the resulting problem are defined as an heuristic to retrieve DA energy prices.
"""
struct UCED <: DayAhead
    system_uc::System
    system_ed::System
    template_uc::OperationsProblemTemplate
    template_ed::OperationsProblemTemplate
    solver_uc::Any
    solver_ed::Any
    kwargs::Dict
end

function UCED(;
    system_uc::System,
    system_ed::System,
    template_uc::OperationsProblemTemplate,
    template_ed::OperationsProblemTemplate,
    solver_uc::Any,
    solver_ed::Any,
    kwargs=if template_uc.transmission == StandardPTDFModel
        Dict(:PTDF => PSY.PTDF(system_ed))
    else
        Dict()
    end,
)
    return UCED(
        system_uc, system_ed, template_uc, template_ed, solver_uc, solver_ed, kwargs
    )
end

# """
#     UC <: DayAhead

# Structure with one DA system is stored (`system_uc`) representing a Unit Commitment (UC) problem
# according to a defined template (`template_uc`), solvable by a passed solver (`solver_uc`). 
# In this case, no duals can be retrieved thus no energy prices.
# """
# struct UC <: DayAhead
#     system_uc::System
#     ptdf::Union{PSY.PTDF, Nothing}
#     template_uc::Dict
#     solver_uc
# end
