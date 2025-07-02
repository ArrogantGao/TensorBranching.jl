function save_weights(filename::String, weights::UnitWeight)
    f = open(filename, "w")
    write(f, "UnitWeight\n")
    write(f, string(weights.n) * "\n")
    close(f)
    return nothing
end

function save_weights(filename::String, weights::Vector{T}) where T
    f = open(filename, "w")
    write(f, "$(T)\n")
    for w in weights
        write(f, string(w) * "\n")
    end
    close(f)
    return nothing
end

function load_weights(filename::String)
    open(filename, "r") do f
        line1 = readline(f)
        if line1 == "UnitWeight"
            return UnitWeight(parse(Int, readline(f)))
        else
            if line1 == "Float64"
                T = Float64
            elseif line1 == "Int64"
                T = Int64
            else
                error("Unknown weight type: $line1")
            end
            return [parse(T, line) for line in eachline(f)]
        end
    end
    return nothing
end

function threaded_saveslices(dirname::String, slices::Vector{SlicedBranch}, ids::Vector{Int})
    @assert length(slices) == length(ids)
    Threads.@threads for i in 1:length(slices)
        br = slices[i]
        id = ids[i]
        if !isnothing(br.code)
            save_finished(dirname, br, id)
        end
    end
    return nothing
end

function save_finished(dirname::String, branch::SlicedBranch, id::Int)
    graph_name = joinpath(dirname, "graph_$(id).dot")
    code_name = joinpath(dirname, "eincode_$(id).json")
    weights_name = joinpath(dirname, "weights_$id.txt")
    savegraph(graph_name, branch.p.g)
    writejson(code_name, uncompress(branch.code))
    save_weights(weights_name, branch.p.weights)
    return nothing
end

function save_unfinished(dirname::String, branch::SlicedBranch, id::Int)
    filename = joinpath(dirname, "unfinished_$(id).jld")
    @save filename branch
    return nothing
end

function load_unfinished(dirname::String, id::Int)
    filename = joinpath(dirname, "unfinished_$(id).jld")
    branch = load(filename, "branch")
    return branch
end