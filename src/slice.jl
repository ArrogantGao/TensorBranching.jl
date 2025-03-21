# generate slices from the kernelized graph

function slice(branch::SlicedBranch, slicer::AbstractSlicer, reducer::AbstractReducer; verbose::Int = 0)
    size_dict = uniformsize(branch.code, 2)
    slices = Vector{SlicedBranch{Int}}()
    _slice!(slices, branch, slicer, reducer, size_dict, verbose)
    return slices
end

function _slice!(slices::Vector{SlicedBranch{Int}}, branch::SlicedBranch{Int}, slicer::AbstractSlicer, reducer::AbstractReducer, size_dict::Dict{Int, Int}, verbose::Int)

    if nv(branch.g) == 0
        push!(slices, branch)
        return nothing
    end

    if (contraction_complexity(branch.code, size_dict).sc ≤ slicer.sc_target)
        push!(slices, branch)
        (verbose ≥ 2) && (@info "current num of slices: $(length(slices))")
        return nothing
    end

    # res is a vector of (mask, code), each corresponding to a slice
    region, loss = ob_region(branch.g, branch.code, slicer, slicer.region_selector, size_dict, verbose)
    branches = optimal_branches(branch.g, branch.code, slicer, reducer, region, size_dict, verbose)

    for newbranch in branches
        _slice!(slices, SlicedBranch(newbranch.g, newbranch.code, branch.r + newbranch.r), slicer, reducer, size_dict, verbose)
    end

    return nothing
end