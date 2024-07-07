using TensorBranching
using Documenter

DocMeta.setdocmeta!(TensorBranching, :DocTestSetup, :(using TensorBranching); recursive=true)

makedocs(;
    modules=[TensorBranching],
    authors="Xuanzhao Gao <gaoxuanzhao@gmail.com> and contributors",
    sitename="TensorBranching.jl",
    format=Documenter.HTML(;
        canonical="https://ArrogantGao.github.io/TensorBranching.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "API" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/ArrogantGao/TensorBranching.jl",
    devbranch="main",
)
