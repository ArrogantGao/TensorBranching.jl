# generate slices from the kernelized graph

function slice(g::SimpleGraph{Int}, code::DynamicNestedEinsum{Int}, slicer::AbstractSlicer)
    size_dict = uniformsize(code, 2)
    slices = Vector{SlicedBranch{Int}}()
    _slice!(slices, SlicedBranch(g, code, 0), slicer, size_dict)
    return slices
end

function _slice!(slices::Vector{SlicedBranch{Int}}, branch::SlicedBranch{Int}, slicer::AbstractSlicer, size_dict::Dict{Int, Int})

    if (contraction_complexity(branch.code, size_dict).sc ≤ slicer.sc_target)
        push!(slices, branch)
        return nothing
    end

    r0 = branch.r

    # res is a vector of (mask, code), each corresponding to a slice
    branches = dynamic_ob(branch.g, branch.code, slicer, size_dict)

    for newbranch in branches
        _slice!(slices, SlicedBranch(newbranch.g, newbranch.code, r0 + newbranch.r), slicer, size_dict)
    end

    return nothing
end

function dynamic_ob(g::SimpleGraph{Int}, code::DynamicNestedEinsum{Int}, slicer::ContractionTreeSlicer, size_dict::Dict{Int, Int})

    region, loss = ob_region(g, code, slicer, slicer.region_selector, size_dict)
    branches = optimal_branches(g, code, slicer, region, size_dict)

    return branches
end