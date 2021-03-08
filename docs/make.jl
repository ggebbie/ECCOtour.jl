using Reemergence
using Documenter

DocMeta.setdocmeta!(Reemergence, :DocTestSetup, :(using Reemergence); recursive=true)

makedocs(;
    modules=[Reemergence],
    authors="G Jake Gebbie",
    repo="https://github.com/ggebbie/Reemergence.jl/blob/{commit}{path}#{line}",
    sitename="Reemergence.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ggebbie.github.io/Reemergence.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ggebbie/Reemergence.jl",
)
