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
`system_ed` and the resulting problem are defined as a heuristic to calculate DA energy prices.
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

"""
    UCED(;
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

Function that returns Unit Commitment (UC) and Economic Dispatch (ED) problems.
"""
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

"""
    DART <: MarketSimulator

A type that represents DayAhead (DA) and Real Time (RT) market clearings. 
"""
abstract type DART <: MarketSimulator end

"""
    UCRT <: DART

Structure with one DA system stored (`system_uc`) and one RT system stored (`system_ed`). 
One represents an Unit Commitment (UC) used to fix binary variables, the other the Real Time (RT)
defined by template `template_rt`. Each model needs its own solver (`solver_uc` and `solver_rt`).
This one is used to make possible to evaluate the primal dispatch for RT prices. 
"""
struct UCRT <: DART
    system_uc::System
    system_rt::System
    template_uc::OperationsProblemTemplate
    template_rt::OperationsProblemTemplate
    solver_uc::Any
    solver_rt::Any
    kwargs::Dict
end

"""
    UCRT(;
        system_uc::System,
        system_rt::System,
        template_uc::OperationsProblemTemplate,
        template_rt::OperationsProblemTemplate,
        solver_uc::Any,
        solver_rt::Any,
        kwargs=if template_uc.transmission == StandardPTDFModel
            Dict(:PTDF => PSY.PTDF(system_ed))
        else
            Dict()
        end,
)

Function that returns Unit Commitment (UC) and Real Time (RT) problems.
"""

function UCRT(;
    system_uc::System,
    system_rt::System,
    template_uc::OperationsProblemTemplate,
    template_rt::OperationsProblemTemplate,
    solver_uc::Any,
    solver_rt::Any,
    kwargs=if template_uc.transmission == StandardPTDFModel
        Dict(:PTDF => PSY.PTDF(system_ed))
    else
        Dict()
    end,
)
    return UCRT(
        system_uc, system_rt, template_uc, template_rt, solver_uc, solver_rt, kwargs
    )
end

"""
    UCEDRT <: DART

Structure with two DA systems stored (`system_uc` and `system_ed`) and one RT system stored (`system_ed`). 
One represents an Unit Commitment (UC) used to fix binary variables, other an Economic Dispatch (ED)
defined by template `template_ed` and the other the Real Time (RT) defined by template `template_rt`. 
Each model needs its own solver (`solver_uc`, `solver_ed` and `solver_rt`).
This one is used to make possible to evaluate DA and RT prices. 
"""
struct UCEDRT <: DART
    system_uc::System
    system_ed::System
    system_rt::System
    template_uc::OperationsProblemTemplate
    template_ed::OperationsProblemTemplate
    template_rt::OperationsProblemTemplate
    solver_uc::Any
    solver_ed::Any
    solver_rt::Any
    kwargs::Dict
end

"""
    UCEDRT(;
        system_uc::System,
        system_ed::System,
        system_rt::System,
        template_uc::OperationsProblemTemplate,
        template_ed::OperationsProblemTemplate,
        template_rt::OperationsProblemTemplate,
        solver_uc::Any,
        solver_ed::Any,
        solver_rt::Any,
        kwargs=if template_uc.transmission == StandardPTDFModel
            Dict(:PTDF => PSY.PTDF(system_ed))
        else
            Dict()
        end,

)

Function that returns Unit Commitment (UC), Economic Dispatch (ED) and Real Time (RT) problems.
"""
function UCEDRT(;
    system_uc::System,
    system_ed::System,
    system_rt::System,
    template_uc::OperationsProblemTemplate,
    template_ed::OperationsProblemTemplate,
    template_rt::OperationsProblemTemplate,
    solver_uc::Any,
    solver_ed::Any,
    solver_rt::Any,
    kwargs=if template_uc.transmission == StandardPTDFModel
        Dict(:PTDF => PSY.PTDF(system_ed))
    else
        Dict()
    end,
)
    return UCEDRT(
        system_uc,
        system_ed,
        system_rt,
        template_uc,
        template_rt,
        template_ed,
        solver_uc,
        solver_rt,
        solver_ed,
        kwargs,
    )
end

# TODO: Create funtion to build UC simulator. Issue: https://github.com/LAMPSPUC/GridAnalysis.jl/issues/3
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
