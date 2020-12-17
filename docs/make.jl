using SeisPGVPGA
using Documenter

makedocs(;
    modules=[SeisPGVPGA],
    authors="kura-okubo",
    repo="https://github.com/kura-okubo/SeisPGVPGA.jl/blob/{commit}{path}#L{line}",
    sitename="SeisPGVPGA.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://kura-okubo.github.io/SeisPGVPGA.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/kura-okubo/SeisPGVPGA.jl",
)
