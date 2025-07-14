using WaveFront
using Documenter

DocMeta.setdocmeta!(WaveFront, :DocTestSetup, :(using WaveFront); recursive=true)

makedocs(;
    modules=[WaveFront],
    authors="William Tegtow <w.tegtow@gmail.com> and contributors",
    sitename="WaveFront.jl",
    format=Documenter.HTML(;
        canonical="https://wtegtow.github.io/WaveFront.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/wtegtow/WaveFront.jl",
    devbranch="main",
)
