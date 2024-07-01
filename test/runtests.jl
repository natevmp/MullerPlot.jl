using Revise
using MullerPlot
using Test
using CSV, DataFrames

# load the data
n_t_vid = CSV.File("data/dataMullerplot_sizes.csv", header=false) |> Tables.matrix
parentVid_vid = CSV.File("data/dataMullerplot_parentLineage.csv", header=false) |> Tables.matrix |> vec
tBirth_vid = CSV.File("data/dataMullerplot_arrivalTimes.csv", header=false) |> Tables.matrix |> vec
visible_vid = n_t_vid[end,:] .> 0.01

vid_child_Vid = MullerPlot.buildFamilyArray(parentVid_vid)

@test vid_child_Vid isa Vector{Vector{Int64}}

xL_t_vid, xU_t_vid = MullerPlot.mullerBounds(n_t_vid, parentVid_vid)
@test xL_t_vid isa Array{Float64,2}
@test xU_t_vid isa Array{Float64,2}

xL_t_vid, xU_t_vid = MullerPlot.mullerBounds(n_t_vid, parentVid_vid, visible_vid)
@test xL_t_vid isa Array{Float64,2}
@test xU_t_vid isa Array{Float64,2}
