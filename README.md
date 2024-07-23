# MullerPlot

[![Build Status](https://github.com/natevmp/MullerPlot.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/natevmp/MullerPlot.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://natevmp.github.io/MullerPlot.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://natevmp.github.io/MullerPlot.jl/dev)

This package contains tools for designing a [Muller plot](https://en.wikipedia.org/wiki/Muller_plot), which is a visual representation of a time-varying population in which subclones expand and/or contract, that simultaneously demonstrates the clone's growth curves and their parent-child relationships. The core functionality takes subclone size information and returns 2D coordinate data which can be used in any plotting package.



## Example usage

First we load population size and parental relationship data.

```julia
using CSV, DataFrames

n_t_vid = CSV.File("data/dataMullerplot_sizes.csv", header=false) |> Tables.matrix
parentVid_vid = CSV.File("data/dataMullerplot_parentLineage.csv", header=false) |> Tables.matrix |> vec
tBirth_vid = CSV.File("data/dataMullerplot_arrivalTimes.csv", header=false) |> Tables.matrix |> vec
```

To use the package functionality, the input data has several requirements. `n_t_vid` must be a matrix of clone sizes, where the first dimension is indexed by the time points at which the sizes are measured, and the second by the different clones. Each clone in the system is designated an integer id that coincides with its index in any array indexing over the different clones:

|`n_t_vid`|  clone 1   |  clone 2   |  clone 3   | ... |
|:--- |:----------: |:----------: | :---------: |:-----:|
|$t_1$| $n_1(t_1)$ | $n_2(t_1)$ | $n_3(t_1)$ |...
|$t_2$| $n_1(t_2)$ | $n_2(t_2)$ | $n_3(t_2)$ |...
|$t_3$| $n_1(t_3)$ | $n_2(t_3)$ | $n_3(t_3)$ |...
|...| ... | ... | ... |...

Clone sizes $n_i(t)$ should be specified as a fraction of the total population size and exclude any subclones, so that $\forall j : \sum_i n_i(t_j) = 1$. 

`parentVid_vid` is a vector indexed by the different clone id's, containing the id (/index) of the clone's direct parent.  For example, `parentVid_vid[5]` returns the id/index of the parent of the clone at index 5. All clones that do not have direct parent in the system, or who's direct parent is the wild type, should be given the parent with id 0.


Next we create bounds for plotting each clone.
```julia
using MullerPlot

xL_t_vid, xU_t_vid = mullerBounds(n_t_vid, parentVid_vid)
```

The method `mullerBounds(...)` returns two matrices containing the coordinates of respectively the upper and lower bounds of each clone (newly arising population) in the population. The first dimension of these  matrices indexes the measurement times (typically the x-axis of a Muller plot), and the second the different clones. For example `xL_t_vid[:,5]` is a vector containing the size-axis coordinate (usually the y-axis) of the lower bound of the 5th clone at all timepoints, and `xL_t_vid[5,:]` is a vector containing the size-axis coordinates of all clones at the 5th timepoint.

It is also possible to select only a subset of clones for the plot, using a boolean vector as third argument. For example, if we wish to only show clones that consist of at least 1% of the population at the final measurement time:
```julia
visible_vid = n_t_vid[end,:] .> 0.01
xL_t_vid, xU_t_vid = MullerPlot.mullerBounds(n_t_vid, parentVid_vid, visible_vid)
```
Now we use the `Makie.jl` package to construct the plot.
```julia
_t = range(0,100,length=size(n_t_vid,1)) # set equidistant measurement times

using CairoMakie

fig = Figure(size=(800,600), fontsize=25)
Axis(
    fig[1,1],
    backgroundcolor=:grey75,
    xgridvisible=false,
    ygridvisible=false,
    yticksvisible=false,
    yticklabelsvisible=false,
    xticks=[0,25,50,75,100],
    xlabel="time",
)
for i in 2:size(xL_t_vid,2) # the first index is the wild type
    band!(
        _t, xL_t_vid[:,i], xU_t_vid[:,i],
    )
    scatter!(
        tBirth_vid[visible_vid][i],
        xU_t_vid[Int(round(tBirth_vid[visible_vid][i]))+1, i],
        markersize=15,
        marker=:star4,
    )
    println(tBirth_vid[i])
end
ylims!(0,1)
xlims!(0,100)
display(fig)```