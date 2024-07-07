function graph_from_tuples(n::Int, edgs)
    g = SimpleGraph(n)
    for (i, j) in edgs
        add_edge!(g, i, j)
    end
    g
end

# Let us create a function for finding reduced ``\alpha``-tensors."
function reduced_alpha(g::SimpleGraph, openvertices::Vector{Int})
	problem = GenericTensorNetwork(IndependentSet(g); openvertices, optimizer = GreedyMethod(nrepeat=1))
	alpha_tensor = solve(problem, SizeMax())
	return mis_compactify!(alpha_tensor)
end

function _reduced_alpha_configs(g::SimpleGraph, openvertices::Vector{Int})
	problem = GenericTensorNetwork(IndependentSet(g); openvertices, optimizer = GreedyMethod(nrepeat=1))
	alpha_tensor = solve(problem, SizeMax())
	alpha_configs = solve(problem, ConfigsMax(; bounded=false))
	reduced_alpha_tensor = mis_compactify!(alpha_tensor)
	# set the corresponding entries to 0.
	alpha_configs[map(iszero, reduced_alpha_tensor)] .= Ref(zero(eltype(alpha_configs)))
	# post processing
	configs = alpha_configs
	return configs
end

# Now we collect these configurations into a vector.
function collect_configs(cfg::CountingTropical{<:Real, <:ConfigEnumerator}, symbols::Union{Nothing, String}=nothing)
    cs = cfg.c.data
    symbols === nothing ? cs : [String([symbols[i] for (i, v) in enumerate(x) if v == 1]) for x in cs]
end

"""
    BranchingTable{N, C}

A table of branching configurations. The table is a vector of vectors of `StaticBitVector{N, C}`. Type parameters are:
- `N`: The number of bits in the bit vector.
- `C`: The number of integers as the storage.

### Examples
```jldoctest
julia> graph_sat = graph_from_tuples(3, [(1, 2), (2, 3), (1, 3)])
{3, 3} undirected simple Int64 graph

julia> tbl = reduced_alpha_configs(graph_sat, [1, 2])
BranchingTable{N}
001
```

To cover the branching table, at least one clause in each row must be satisfied.
"""
struct BranchingTable{N, C}
    table::Vector{Vector{StaticBitVector{N, C}}}
end
function BranchingTable(arr::AbstractArray{<:CountingTropical{<:Real, <:ConfigEnumerator}})
    return BranchingTable(filter(!isempty, vec(map(collect_configs, arr))))
end
function BranchingTable(n::Int, arr::AbstractVector{<:AbstractVector})
    return BranchingTable([StaticBitVector(n, x) for x in arr])
end
Base.:(==)(t1::BranchingTable, t2::BranchingTable) = all(x -> Set(x[1]) == Set(x[2]), zip(t1.table, t2.table))
function Base.show(io::IO, t::BranchingTable{N}) where N
    println(io, "BranchingTable{N}")
    for (i, row) in enumerate(t.table)
        print(io, join(["$r" for r in row], ", "))
        i < length(t.table) && println(io)
    end
end
Base.show(io::IO, ::MIME"text/plain", t::BranchingTable) = show(io, t)

# And a combination of the above two procedures.
function reduced_alpha_configs(graph::SimpleGraph, openvertices::Vector{Int})
	configs = _reduced_alpha_configs(graph, openvertices)
    return BranchingTable(configs)
end

struct DNF{N, T}
    clauses::Vector{Clause{N, T}}
end
DNF(c::Clause{N, T}, cs::Clause{N, T}...) where {N, T} = DNF([c, cs...])
Base.:(==)(x::DNF, y::DNF) = x.clauses == y.clauses
Base.length(x::DNF) = length(x.clauses)

function covered_by(t::BranchingTable, dnf::DNF)
    all(x->any(y->covered_by(y, dnf), x), t.table)
end
function covered_by(s::StaticBitVector, dnf::DNF)
    @assert length(s.data) == 1 "length of StaticBitVector too long, not yet supported."
    any(c->covered_by(s.data[1], c), dnf.clauses)
end

function Tbl2BitStrs(tbl::BranchingTable{N}) where N
    return [BitStr.(x) for x in tbl.table]
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


"""
    env_filter(tbl::BranchingTable{N, C}, g::SimpleGraph, vertices::Vector{Int}, openvertices::Vector{Int}) -> BranchingTable

Filter the given `tbl` based on the environment of the vertices.

# Arguments
- `tbl::BranchingTable{N, C}`: The branching table to filter.
- `g::SimpleGraph`: The graph representing the environment.
- `vertices::Vector{Int}`: The vertices to consider.
- `openvertices::Vector{Int}`: The open vertices.

# Returns
A new `BranchingTable` object containing the filtered rows.

"""
function env_filter(tbl::BranchingTable{N, C}, g::SimpleGraph, vertices::Vector{Int}, openvertices::Vector{Int}) where {N, C}
    tbl_list = vcat([BitStr.(x) for x in tbl.table]...)
    mis_S = maximum(count_ones.(tbl_list))

    ns = neighbors(g, vertices)
    so = Set(openvertices)

    new_table = Vector{Vector{StaticBitVector{N, C}}}()
    for row in tbl.table
        x = row[1]
        s1 = Int[]
        for i in 1:length(x)
            if (x[i] == 1) && (vertices[i] ∈ so)
                push!(s1, vertices[i])
            end
        end
        ns1 = neighbors(g, s1)
        ns0 = ns ∩ ns1

        if length(ns) - length(ns0) + count_ones(BitStr(x)) ≥ mis_S
            push!(new_table, row)
        end
    end

    return BranchingTable(new_table)
end