using GridAnalysis
using Documenter

DocMeta.setdocmeta!(GridAnalysis, :DocTestSetup, :(using GridAnalysis); recursive=true)

makedocs(;
    modules=[GridAnalysis],
    authors="Invenia Technical Computing Corporation",
    repo="https://github.com/invenia/GridAnalysis.jl/blob/{commit}{path}#{line}",
    sitename="GridAnalysis.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://invenia.github.io/GridAnalysis.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
    checkdocs=:exports,
    strict=true,
)

deploydocs(;
    repo="github.com/invenia/GridAnalysis.jl",
)
