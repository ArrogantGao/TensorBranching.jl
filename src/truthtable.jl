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
    BranchingTable{INT}

A table of branching configurations. The table is a vector of vectors of `INT`. Type parameters are:
- `INT`: The number of integers as the storage.

# Fields
- `bit_length::Int`: The length of the bit string.
- `table::Vector{Vector{INT}}`: The table of bitstrings used for branching.

### Examples
```jldoctest
julia> graph_sat = graph_from_tuples(3, [(1, 2), (2, 3), (1, 3)])
{3, 3} undirected simple Int64 graph

julia> tbl = reduced_alpha_configs(TensorNetworkSolver(), graph_sat, [1, 2])
BranchingTable{LongLongUInt{1}}
100
```

To cover the branching table, at least one clause in each row must be satisfied.
"""
struct BranchingTable{INT <: Integer}
    bit_length::Int
    table::Vector{Vector{INT}}
end
function BranchingTable(arr::AbstractArray{<:CountingTropical{<:Real, <:ConfigEnumerator{N}}}) where N
    return BranchingTable(N, filter(!isempty, vec(map(collect_configs, arr))))
end
function BranchingTable(n::Int, arr::AbstractVector{<:AbstractVector})
    return BranchingTable(n, [_vec2int.(LongLongUInt, x) for x in arr])
end
nbits(t::BranchingTable) = t.bit_length
Base.:(==)(t1::BranchingTable, t2::BranchingTable) = all(x -> Set(x[1]) == Set(x[2]), zip(t1.table, t2.table))
function Base.show(io::IO, t::BranchingTable{INT}) where INT
    println(io, "BranchingTable{$INT}")
    for (i, row) in enumerate(t.table)
        print(io, join(["$(bitstring(r)[end-nbits(t)+1:end])" for r in row], ", "))
        i < length(t.table) && println(io)
    end
end
Base.show(io::IO, ::MIME"text/plain", t::BranchingTable) = show(io, t)

abstract type AbstractMISSolver end
struct TensorNetworkSolver <: AbstractMISSolver end

# And a combination of the above two procedures.
function reduced_alpha_configs(::TensorNetworkSolver, graph::SimpleGraph, openvertices::Vector{Int})
	configs = _reduced_alpha_configs(graph, openvertices)
    return BranchingTable(configs)
end

struct DNF{INT}
    clauses::Vector{Clause{INT}}
end
DNF(c::Clause{INT}, cs::Clause{INT}...) where {INT} = DNF([c, cs...])
Base.:(==)(x::DNF, y::DNF) = x.clauses == y.clauses
Base.length(x::DNF) = length(x.clauses)

function covered_by(t::BranchingTable, dnf::DNF)
    all(x->any(y->covered_by(y, dnf), x), t.table)
end
function covered_by(s::Integer, dnf::DNF)
    any(c->covered_by(s, c), dnf.clauses)
end