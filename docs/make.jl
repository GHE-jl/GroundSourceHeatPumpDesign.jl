using GroundSourceHeatPumpDesign
using Documenter

# Make `using GroundSourceHeatPumpDesign` available to every doctest in docstrings and pages.
DocMeta.setdocmeta!(
    GroundSourceHeatPumpDesign,
    :DocTestSetup,
    :(using GroundSourceHeatPumpDesign);
    recursive = true,
)

makedocs(;
    modules = [GroundSourceHeatPumpDesign],
    authors = "Gabriel-Dion <dion.gabriel100@gmail.com>",
    sitename = "GroundSourceHeatPumpDesign.jl",
    format = Documenter.HTML(;
        canonical = "https://GHE-jl.github.io/GroundSourceHeatPumpDesign.jl",
        edit_link = "main",
        assets = String[],
        mathengine = Documenter.KaTeX(),
        sidebar_sitename = false,
    ),
    pages = [
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Modeling theory" => "theory.md",
        "API reference" => "api.md",
        "References" => "references.md",
    ],
    # Keep the build strict so broken cross-references or missing docstrings fail CI.
    checkdocs = :exports,
)

deploydocs(;
    repo = "github.com/GHE-jl/GroundSourceHeatPumpDesign.jl",
    devbranch = "main",
)
