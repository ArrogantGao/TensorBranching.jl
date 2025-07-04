function save_weights(filename::String, weights::UnitWeight)
    f = open(filename, "w")
    write(f, "UnitWeight\n")
    write(f, string(weights.n) * "\n")
    close(f)
    return nothing
end

function save_weights(filename::String, weights::Vector{T}) where T
    # if f ex
    if isfile(filename)
        rm(filename)
    end
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
            elseif line1 == "Float32"
                T = Float32
            elseif line1 == "Int64"
                T = Int64
            elseif line1 == "Int32"
                T = Int32
            else
                error("Unknown weight type: $line1")
            end
            weights = [parse(T, line) for line in eachline(f)]
            return weights
        end
    end
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

function save_code(filename::String, code::DynamicNestedEinsum)
    writejson(filename, code)
    return nothing
end

function save_code(filename::String, code::Nothing)
    open(filename, "w") do f
        write(f, "nothing")
    end
    return nothing
end

function load_code(filename::String)
    f = open(filename, "r")
    line1 = readline(f)
    if line1 == "nothing"
        return nothing
    else
        close(f)
        return readjson(filename)
    end
end

function save_finished(dirname::String, branch::SlicedBranch, id::Int)
    graph_name = joinpath(dirname, "graph_$(id).dot")
    code_name = joinpath(dirname, "eincode_$(id).json")
    weights_name = joinpath(dirname, "weights_$id.txt")
    savegraph(graph_name, branch.p.g)
    save_code(code_name, uncompress(branch.code))
    save_weights(weights_name, branch.p.weights)
    return nothing
end

function load_finished(dirname::String, id::Int)
    graph_name = joinpath(dirname, "graph_$(id).dot")
    code_name = joinpath(dirname, "eincode_$(id).json")
    weights_name = joinpath(dirname, "weights_$id.txt")
    g = loadgraph(graph_name)
    code = load_code(code_name)
    weights = load_weights(weights_name)
    return g, code, weights
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

function load_all_finished(dirname::String)
    df = CSV.read(joinpath(dirname, "slices.csv"), DataFrame)
    branches = SlicedBranch[]
    for id in df.id
        g, code, weights = load_finished(dirname, id)
        push!(branches, SlicedBranch(MISProblem(g, weights), code, df.r[df.id .== id][1]))
    end
    return branches
end