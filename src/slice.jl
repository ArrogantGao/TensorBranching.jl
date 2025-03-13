# generate slices from the kernelized graph

function slice(g::SimpleGraph{Int}, code::DynamicNestedEinsum{Int}, slicer::ContractionTreeSlicer)
    size_dict = uniformsize(code, 2)
    mask = bmask(BitBasis.longinttype(nv(g), 2), 1:nv(g))
    slices = Vector{Tuple{typeof(mask), DynamicNestedEinsum{Int}}}()
    _slice!(slices, g, code, slicer, mask, size_dict)
    return slices
end

function _slice!(slices::Vector{Tuple{INT, DynamicNestedEinsum{Int}}}, g::SimpleGraph{Int}, code::DynamicNestedEinsum{Int}, slicer::ContractionTreeSlicer, mask::INT, size_dict::Dict{Int, Int}) where{INT}
    if (contraction_complexity(code, size_dict).sc ≤ slicer.sc_target)
        push!(slices, (mask, code))
        return nothing
    end

    # res is a vector of (mask, code), each corresponding to a slice
    branches = dynamic_ob_slicing(g, code, slicer, mask, size_dict)

    for (mask, code) in branches
        _slice!(slices, g, code, slicer, mask, size_dict)
    end
    return nothing
end

function dynamic_ob_slicing(g::SimpleGraph{Int}, code::DynamicNestedEinsum{Int}, slicer::ContractionTreeSlicer, mask::INT, size_dict::Dict{Int, Int}) where{INT}

    vs = ob_region()
    branches = optimal_branches(vs, g, code, slicer, mask, size_dict)

    return branches
end

function ob_region(g::SimpleGraph{Int}, code::DynamicNestedEinsum{Int}, slicer::ContractionTreeSlicer, mask::INT, size_dict::Dict{Int, Int}) where{INT}

    large_tensors = list_subtree(code, size_dict, slicer.threshold)
    large_tensors_iy = [t.eins.iy for t in large_tensors]

    # a heuristic way: try to find the 

end