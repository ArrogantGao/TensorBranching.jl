function saveslices(filename::String, slices::Union{Vector{SlicedBranch{TS}},Vector{SlicedWeightedBranch{TS}}}) where {TS}
    @save filename slices
    return nothing
end

function threaded_saveslices(dirname::String, slices::Union{Vector{SlicedBranch{TS}},Vector{SlicedWeightedBranch{TS}}}, ids::Vector{Int}) where {TS}
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

function save_finished(dirname::String, branch::Union{SlicedBranch, SlicedWeightedBranch}, id::Int)
    graph_name = joinpath(dirname, "graph_$(id).dot")
    code_name = joinpath(dirname, "eincode_$(id).json")
    savegraph(graph_name, branch.g)
    writejson(code_name, uncompress(branch.code))
    return nothing
end

function save_unfinished(dirname::String, branch::Union{SlicedBranch, SlicedWeightedBranch}, reducer::AbstractReducer, id::Int)
    filename = joinpath(dirname, "unfinished_$(id).jld")
    @save filename branch reducer
    return nothing
end

function load_unfinished(dirname::String, id::Int)
    filename = joinpath(dirname, "unfinished_$(id).jld")
    branch = load(filename, "branch")
    reducer = load(filename, "reducer")
    return branch, reducer
end