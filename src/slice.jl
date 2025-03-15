# generate slices from the kernelized graph

function slice(branch::SlicedBranch, slicer::AbstractSlicer)
    size_dict = uniformsize(branch.code, 2)
    slices = Vector{SlicedBranch{Int}}()
    _slice!(slices, branch, slicer, size_dict)
    return slices
end

function _slice!(slices::Vector{SlicedBranch{Int}}, branch::SlicedBranch{Int}, slicer::AbstractSlicer, size_dict::Dict{Int, Int})

    if (contraction_complexity(branch.code, size_dict).sc ≤ slicer.sc_target)
        push!(slices, branch)
        return nothing
    end

    # res is a vector of (mask, code), each corresponding to a slice
    region, loss = ob_region(branch.g, branch.code, slicer, slicer.region_selector, size_dict)
    branches = optimal_branches(branch.g, branch.code, slicer, region, size_dict)

    for newbranch in branches
        _slice!(slices, SlicedBranch(newbranch.g, newbranch.code, branch.r + newbranch.r), slicer, size_dict)
    end

    return nothing
end