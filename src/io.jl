function saveslices(filename::String, slices::Vector{SlicedBranch{TS}}) where {TS}
    @save filename slices
    return nothing
end

function loadslices(filename::String)
    slices = load(filename, "slices")
    return slices
end