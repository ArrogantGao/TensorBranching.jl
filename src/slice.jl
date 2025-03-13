# generate slices from the kernelized graph

function slice(g::SimpleGraph{Int}, code::DynamicNestedEinsum{Int}, slicer::AbstractSlicer)
    size_dict = uniformsize(code, 2)
    slices = Vector{Tuple{SimpleGraph{Int}, DynamicNestedEinsum{Int}}}()
    _slice!(slices, g, code, slicer, size_dict)
    return slices
end

function _slice!(slices::Vector{Tuple{SimpleGraph{Int}, DynamicNestedEinsum{Int}}}, g::SimpleGraph{Int}, code::DynamicNestedEinsum{Int}, slicer::AbstractSlicer, size_dict::Dict{Int, Int})

    if (contraction_complexity(code, size_dict).sc ≤ slicer.sc_target)
        push!(slices, (g, code))
        return nothing
    end

    # res is a vector of (mask, code), each corresponding to a slice
    branches = dynamic_ob_slicing(g, code, slicer, size_dict)

    for (g, code) in branches
        _slice!(slices, g, code, slicer, size_dict)
    end

    return nothing
end

function dynamic_ob_slicing(g::SimpleGraph{Int}, code::DynamicNestedEinsum{Int}, slicer::ContractionTreeSlicer, size_dict::Dict{Int, Int})

    region, loss = ob_region(g, code, slicer, slicer.region_selector, size_dict)
    branches = optimal_branches(g, code, slicer, region, size_dict)

    return branches
end