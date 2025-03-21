function saveslices(slices::Vector{SlicedBranch{TS}}, filename::String) where {TS}
    @save filename slices
    return nothing
end

function loadslices(filename::String)
    slices = load(filename, "slices")
    return slices
end