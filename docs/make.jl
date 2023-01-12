using Deepburn
using Documenter

DocMeta.setdocmeta!(Deepburn, :DocTestSetup, :(using Deepburn); recursive=true)

makedocs(;
    modules=[Deepburn],
    authors="Van den Eynde Gert <gert.van.den.eynde@sckcen.be> and contributors",
    repo="https://github.com/gvdeynde/Deepburn.jl/blob/{commit}{path}#{line}",
    sitename="Deepburn.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://gvdeynde.github.io/Deepburn.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/gvdeynde/Deepburn.jl",
    devbranch="master",
)
