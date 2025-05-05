using CNNSeq2exp
using Documenter

DocMeta.setdocmeta!(CNNSeq2exp, :DocTestSetup, :(using CNNSeq2exp); recursive=true)

makedocs(;
    modules=[CNNSeq2exp],
    authors="anon",
    sitename="CNNSeq2exp.jl",
    format=Documenter.HTML(;
        canonical="https://anon-research-2025.github.io/CNNSeq2exp.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/anon-research-2025/CNNSeq2exp.jl",
    devbranch="main",
)
