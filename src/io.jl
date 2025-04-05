function saveslices(filename::String, slices::Vector{SlicedBranch{TS}}) where {TS}
    @save filename slices
    return nothing
end

function threaded_saveslices(dirname::String, slices::Vector{SlicedBranch{TS}}, ids::Vector{Int}) where {TS}
    @assert length(slices) == length(ids)
    Threads.@threads for i in 1:length(slices)
        br = slices[i]
        id = ids[i]
        if !isnothing(br.code)
            graph_name = joinpath(dirname, "graph_$id.dot")
            code_name = joinpath(dirname, "eincode_$id.json")
            savegraph(graph_name, br.g)
            writejson(code_name, uncompress(br.code))
        end
    end
    return nothing
end

function loadslices(filename::String)
    slices = load(filename, "slices")
    return slices
end