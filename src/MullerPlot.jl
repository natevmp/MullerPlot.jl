module MullerPlot

export mullerBounds, buildFamilyArray, variantBirthTimes

"""
    buildFamilyArray(parentVid_vid)

Create a vector `vid_child_Vid`, indexed by the variant id's `vid`, that has as elements the children of each variant stored as a vector of id's. For example, `vid_child_Vid[5]` is a vector `vid_child` whose elements are the `vid`'s of the children of the variant with `vid=5`.

`parentVid_vid`:    Vector indexed by the variant id's `vid`, where each element refers to that variant's parent `vid`. For example `parentVid_vid[5]` is the `vid` of the parent of the variant with `vid=5`.
"""
function buildFamilyArray(parentVid_vid::Vector{T} where T<:Integer)
    vid_child_Vid = [Int64[] for _ in eachindex(parentVid_vid)]
    for (vid, parentVid) in enumerate(parentVid_vid)
        if parentVid==0 continue end
        push!(vid_child_Vid[parentVid], vid)
    end
    return vid_child_Vid
end

"""
    sizeVariantRec(vid, n_vid, vid_child_Vid)

Get combined size of all child clones of a variant.

`vid`:            variant id of parent \\
`n_vid`:          vector of clone sizes indexed by variant id \\
`vid_child_Vid`:  vector (indexed by variant id) where each element is a vector of that variant's children (given as their variant id's)
"""
function sizeVariantRec(vid::Integer, n_vid::AbstractArray{<:Float64,1}, vid_child_Vid::AbstractArray{<:AbstractArray{<:Integer}})
    childSizes = 0
    for child in vid_child_Vid[vid]
        childSizes += sizeVariantRec(child, n_vid, vid_child_Vid)
    end
    n_vid[vid] + childSizes
end

"""
    getLowerCoordinateFromParent(parent, xL_vid, xU_vid, n_vid, vid_child_Vid)

Get y-axis coordinate for the lower bound of all child clones of a parent clone.

`parent`:   variant id (`vid`) of parent clone \\
`xL_vid`:   vector of lower bounds in draw space indexed by variant id's (`vid`) \\
`xU_vid`:   vector of upper bounds in draw space indexed by variant id's (`vid`) \\
`n_vid`:    vector of clone sizes indexed by variant id \\
`vid_child_Vid`: vector (indexed by variant id) where each element is a vector of that variant's children (given as their variant id's)
"""
function getLowerCoordinateFromParent(parent, xL_vid, xU_vid, n_vid, vid_child_Vid)
    Δparent = xU_vid[parent]-xL_vid[parent]
    Δchildren = sizeVariantRec(parent, n_vid, vid_child_Vid) - n_vid[parent]
    xL_vid[parent] + (Δparent-Δchildren)/2
end

function setCoordinates!(parent, xL_vid, xU_vid, n_vid, vid_child_Vid)
    x0 = getLowerCoordinateFromParent(parent, xL_vid, xU_vid, n_vid, vid_child_Vid)
    for child in vid_child_Vid[parent]
        xL_vid[child] = x0
        xU_vid[child] = x0 + sizeVariantRec(child, n_vid, vid_child_Vid)
        x0 = xU_vid[child] # update x0 so next child has lower bound at current upper bound
    end
    for child in vid_child_Vid[parent]
        setCoordinates!(child, xL_vid, xU_vid, n_vid, vid_child_Vid)
    end
end

"""
    mullerBounds(n_t_vid, vid_child_Vid)

Return lower and upper bounds `xL_t_vid` and `xU_t_vid` of variant bands for creating a Muller plot.

`n_t_vid`: two dimensional Array whose elements are variant clone sizes. The first dimension `_t` is the timepoints at which the measurements occured; the second dimension `_vid` identifies the variant (i.e. the variant id). The first index refers to the wild type. \\
`vid_child_Vid`: A nested array indexed by the variant id's. Its elements are 1D Vectors which contain the variant id's of each variant's children.
"""
function mullerBounds(n_t_vid::AbstractArray{T,2} where T<:Real, vid_child_Vid::Vector{Vector{V}} where V<:Integer)
    _t = 1:size(n_t_vid,1)
    xL_t_vid = Array{Float64}(undef, length(_t), size(vid_child_Vid,1))
    xU_t_vid = Array{Float64}(undef, length(_t), size(vid_child_Vid,1))
    xL_t_vid[:,1] .= 0
    xU_t_vid[:,1] .= 1
    for t in eachindex(_t)
        # recursively set upper and lower bound coordinates for each clone and all of its children
        setCoordinates!(1, (@view xL_t_vid[t,:]), (@view xU_t_vid[t,:]), (@view n_t_vid[t,:]), vid_child_Vid)
    end
    return xL_t_vid, xU_t_vid
end

function mullerBounds(n_t_vid::AbstractArray{T,2} where T<:Real, vid_child_Vid::Vector{Vector{V}} where V<:Integer, visible_vid::AbstractVector{Bool})
    vidT_child_VidT = reduceCloneSpace(vid_child_Vid, visible_vid)
    mullerBounds(n_t_vid[:, visible_vid], vidT_child_VidT)
end

function mullerBounds(n_t_vid::AbstractArray{T,2} where T<:Real, parentVid_vid::Vector{V} where V<:Integer)
    childVid_child_Vid = buildFamilyArray(parentVid_vid)
    mullerBounds(n_t_vid, childVid_child_Vid)
end

function mullerBounds(n_t_vid::AbstractArray{<:Real,2} where T<:Real, parentVid_vid::Vector{V} where V<:Integer, visible_vid::AbstractVector{Bool})
    childVid_child_Vid = buildFamilyArray(parentVid_vid)
    mullerBounds(n_t_vid, childVid_child_Vid, visible_vid)
end

"""
    variantBirthTimes(n_t_vid)

Get the arrival (birth) times of each clone in the system
"""
function variantBirthTimes(n_t_vid::AbstractArray{<:Real,2})
    tBirth_vid = Vector{Int64}(undef, size(n_t_vid,2))
    for (vid, n_t) in enumerate(eachcol(n_t_vid))
        tBirth_vid[vid] = findfirst(n_t .> 0)
    end
    return tBirth_vid
end

function transformVidCoordinate(vidIn::Integer, vid1_vid2)
    return findfirst(vid1_vid2.==vidIn)
end

function reduceCloneSpace(vid_child_Vid::AbstractArray{<:AbstractArray{<:Integer}}, visible_vid::AbstractVector{Bool})
    visible_vid[1] = true # wild type must be visible
    nClones2 = sum(visible_vid)
    _vid = range(1,length(visible_vid))
    vid2_child_Vid2 = Vector{Vector{Int}}(undef, nClones2)
    for (vid2, vid_child) in enumerate(vid_child_Vid[visible_vid])
        vid2_child = Int[]
        for childVid in vid_child
            if !visible_vid[childVid] continue end
            childVid2 = transformVidCoordinate(childVid, _vid[visible_vid])
            push!(vid2_child, childVid2)
        end
        vid2_child_Vid2[vid2] = vid2_child
        if vid2==1 continue end
        # if visible clone has invisible parent it will be missed. Find cases and set parent to wild type
        visibleParent = false
        for i in 1:vid2-1
            visibleParent = visibleParent || in(vid2, vid2_child_Vid2[i])
        end
        if !visibleParent
            push!(vid2_child_Vid2[1], vid2)
        end
    end
    println(vid2_child_Vid2)
    return vid2_child_Vid2
end

end
