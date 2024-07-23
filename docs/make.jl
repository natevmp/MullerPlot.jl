using Documenter
using MullerPlot

makedocs(
    sitename = "MullerPlot.jl Documentation",
    format = Documenter.HTML(),
    modules = [MullerPlot],
    # doctest = false,
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/natevmp/MullerPlot.jl.git",
)