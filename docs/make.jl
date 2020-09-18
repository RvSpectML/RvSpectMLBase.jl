using RvSpectMLBase
using Documenter

makedocs(;
    modules=[RvSpectMLBase],
    authors="Eric Ford",
    repo="https://github.com/RvSpectML/RvSpectMLBase.jl/blob/{commit}{path}#L{line}",
    sitename="RvSpectMLBase.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://RvSpectML.github.io/RvSpectMLBase.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/RvSpectML/RvSpectMLBase.jl",
)
