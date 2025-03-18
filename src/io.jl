function saveslices(slices::Vector{SlicedBranch{TS}}, tensors::Vector{TA}, r::Int64, filename::String) where {TS,TA}
    @save filename slices tensors r
    return nothing
end

function loadslices(filename::String)
    slices, tensors, r = load(filename, "slices", "tensors", "r")
    return slices, tensors, r
end