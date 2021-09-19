using ECCOtour
using Documenter

DocMeta.setdocmeta!(ECCOtour, :DocTestSetup, :(using ECCOtour); recursive=true)

makedocs(;
    modules=[ECCOtour],
    authors="G Jake Gebbie",
    repo="https://github.com/ggebbie/ECCOtour.jl/blob/{commit}{path}#{line}",
    sitename="ECCOtour.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ggebbie.github.io/ECCOtour.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ggebbie/ECCOtour.jl",
)
