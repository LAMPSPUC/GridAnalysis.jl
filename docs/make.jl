using GridAnalysis
using Documenter

DocMeta.setdocmeta!(GridAnalysis, :DocTestSetup, :(using GridAnalysis); recursive=true)

makedocs(;
    modules=[GridAnalysis],
    authors="LAMPS",
    repo="https://github.com/LAMPSPUC/GridAnalysis.jl/blob/{commit}{path}#{line}",
    sitename="GridAnalysis.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://github.com/LAMPSPUC/GridAnalysis.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
    checkdocs=:exports,
    strict=true,
)

deploydocs(;
    repo="github.com/LAMPSPUC/GridAnalysis.jl",
)
