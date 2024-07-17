function removed_vertices(vertices::Vector{Int}, g::SimpleGraph, clause::Clause{N}) where N
    rvs = Int[]
    for (k, v) in enumerate(vertices)
        if readbit(clause.mask, k) == 1
            push!(rvs, v)
            if readbit(clause.val, k) == 1
                append!(rvs, neighbors(g, v))
            end
        end
    end
    return unique!(rvs)
end

function open_vertices(g::SimpleGraph, vertices::Vector{Int})
    return unique!([v for v in vertices if !all(x->x ∈ vertices, neighbors(g, v))])
end

function open_neighbors(g::SimpleGraph, vertices::Vector{Int})
    ov = Vector{Int}()
    for v in vertices
        for n in neighbors(g, v)
            push!(ov, n)
        end
    end
    return unique!(setdiff(ov, vertices))
end

function closed_neighbors(g::SimpleGraph, vertices::Vector{Int})
    return vertices ∪ open_neighbors(g, vertices)
end

function neighbor_cover(g::SimpleGraph, v::Int, k::Int)
    @assert k >= 0
    vertices = [v]
    for _ = 1:k
        vertices = union(vertices, (neighbors(g, w) for w in vertices)...)
    end
    openvertices = open_vertices(g, vertices)
    return vertices, openvertices
end

function neighbors_2nd(g::SimpleGraph, v::Int)
    return open_neighbors(g, v ∪ neighbors(g, v))
end

# vs a subgraph, return N(vs)
function Graphs.neighbors(g::SimpleGraph, vs::Vector{Int})
    set_vs = Set(vs)
    set_neighbors = Set{Int}()
    for v in vs
        neighbors_v = neighbors(g, v)
        for n in neighbors_v
            if n ∉ set_vs
                push!(set_neighbors, n)
            end
        end
    end
    return set_neighbors
end