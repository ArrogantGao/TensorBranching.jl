abstract type BranchingStrategy end

"""
    struct NaiveBranching <: BranchingStrategy
A struct representing the NaiveBranching branching strategy.
"""
struct NaiveBranching <: BranchingStrategy end

"""
    struct SetCoverBranching <: BranchingStrategy

A struct representing a branching strategy for set cover problems.

# Fields
- `max_itr::Int`: The maximum number of iterations.

# Constructors
- `SetCoverBranching()`: Constructs a `SetCoverBranching` object with a default value of `2` for `max_itr`.
- `SetCoverBranching(max_itr::Int)`: Constructs a `SetCoverBranching` object with the specified `max_itr` value.

"""
struct SetCoverBranching <: BranchingStrategy 
    max_itr::Int
    SetCoverBranching() = new(2)
    SetCoverBranching(max_itr::Int) = new(max_itr)
end

"""
    missolve(g::SimpleGraph; show_count = false, strategy::BranchingStrategy = NaiveBranching(), kneighbor::Int = 2, use_rv::Bool = true)

Solves the maximum independent set (MIS) problem for a given graph `g` using a specified branching strategy.

## Arguments
- `g::SimpleGraph`: The input graph.
- `show_count::Bool = false`: Whether to return the count of iterations performed during the solving process.
- `strategy::BranchingStrategy = NaiveBranching()`: The branching strategy to use. Defaults to `NaiveBranching()`.
- `kneighbor::Int = 2`: The number of neighbors to consider during the branching process. Defaults to 2.
- `use_rv::Bool = true`: Whether to use number of removed vertices during the branching process. Defaults to true.
- `use_filter::Bool = true`: Whether to use the environment filter. Defaults to true.

## Returns
- If `show_count` is false, returns the maximum independent set (MIS) of the graph `g`.
- If `show_count` is true, returns a tuple `(mis, count)` where `mis` is the maximum independent set (MIS) of the graph `g` and `count` is the number of iterations performed during the solving process.
"""
function missolve(g::SimpleGraph; show_count = false, strategy::BranchingStrategy = NaiveBranching(), kneighbor::Int = 2, use_rv::Bool = true, use_filter::Bool = true)
    mis = mis_solver(g, strategy, kneighbor, use_rv, use_filter)
    return show_count ? (mis.mis, mis.count) : mis.mis
end

function rv(vertices::Vector{Int}, g::SimpleGraph, clause::Clause{N}) where N
    removed_vertices = Int[]
    for (k, v) in enumerate(vertices)
        if readbit(clause.mask, k) == 1
            push!(removed_vertices, v)
            if readbit(clause.val, k) == 1
                append!(removed_vertices, neighbors(g, v))
            end
        end
    end
    return unique!(removed_vertices)
end

function mis_solver(g::SimpleGraph, strategy::BranchingStrategy, kneighbor::Int, use_rv::Bool, use_filter::Bool)
    dg = degree(g)
    if nv(g) == 0 || nv(g) == 1
        return CountingMIS(nv(g))
    elseif (0 ∈ dg) || (1 ∈ dg)
        v = findfirst(x -> (x==0)||(x==1), dg)
        return 1 + mis_solver(induced_subgraph(g, setdiff(1:nv(g), v ∪ neighbors(g, v)))[1], strategy, kneighbor, use_rv, use_filter)
    elseif maximum(dg) ≥ 6
        v = argmax(dg)
        return max(1 + mis_solver(induced_subgraph(g, setdiff(1:nv(g), v ∪ neighbors(g, v)))[1], strategy, kneighbor, use_rv, use_filter), mis_solver(induced_subgraph(g, setdiff(1:nv(g), v))[1], strategy, kneighbor, use_rv, use_filter))
    else
        vertices, openvertices, dnf = optimal_branching_dnf(g, strategy, kneighbor, use_rv, use_filter)
        @assert !isempty(vertices)

        mis_count = Vector{CountingMIS}(undef, length(dnf.clauses))
        
        @threads for i in 1:length(dnf.clauses)
            clause = dnf.clauses[i]
            removed_vertices = rv(vertices, g, clause)
            gi = copy(g)
            rem_vertices!(gi, removed_vertices)
            @assert !isempty(removed_vertices)
            mis_count[i] = mis_solver(gi, strategy, kneighbor, use_rv, use_filter) + count_ones(clause.val)
        end
        max_mis = max(mis_count...)

        return max_mis
    end
end

function neighbor_cover(g::SimpleGraph, v::Int, k::Int)
    @assert k >= 0
    vertices = [v]
    for _ = 1:k
        vertices = union(vertices, (neighbors(g, w) for w in vertices)...)
    end
    openvertices = [v for v in vertices if !all(x->x ∈ vertices, neighbors(g, v))]
    return vertices, openvertices
end

"""
    optimal_branching_dnf(g::SimpleGraph, strategy::BranchingStrategy, kneighbor::Int, use_rv::Bool)

Compute the optimal branching in a directed acyclic graph (DAG) using a specified branching strategy.

# Arguments
- `g::SimpleGraph`: The input graph.
- `strategy::BranchingStrategy`: The branching strategy to use.
- `kneighbor::Int`: The number of neighbors to consider when selecting the branching vertex.
- `use_rv::Bool`: Whether to use the number of removed vertices during the branching process.
- `use_filter::Bool`: Whether to use the environment filter.

# Returns
- `vertices`: The set of vertices selected for the optimal branching.
- `openvertices`: The set of open vertices in the optimal branching.
- `impl_strategy`: The implementation strategy used for the optimal branching.

"""
function optimal_branching_dnf(g::SimpleGraph, strategy::BranchingStrategy, kneighbor::Int, use_rv::Bool, use_filter::Bool)
    # reference: Exaxt Exponential Algorithms by Fomin and Kratsch, chapter 2.3
    degree_g = degree(g)

    # use minimum boundary instead
    vs_min = Int[]
    ovs_min = [1:nv(g)...]
    for v in 1:nv(g)
        vs, ovs = neighbor_cover(g, v, kneighbor)
        if length(ovs) < length(ovs_min)
            vs_min = vs
            ovs_min = ovs
        end
    end
    vertices = vs_min
    openvertices = ovs_min

    subg, vmap = induced_subgraph(g, vertices)
    tbl = reduced_alpha_configs(subg, Int[findfirst(==(v), vertices) for v in openvertices])
    if use_filter
        filtered_tbl = env_filter(tbl, g, vertices, openvertices)
    else
        filtered_tbl = tbl
    end
    return vertices, openvertices, impl_strategy(filtered_tbl, strategy, vertices, g, use_rv)
end

function impl_strategy(tbl::BranchingTable{N}, strategy::NaiveBranching, vertices, g, use_rv) where N
    return DNF([Clause(bmask(BitStr{N, Int}, 1:N), BitStr(first(x))) for x in tbl.table])
end

function impl_strategy(tbl::BranchingTable{N}, strategy::SetCoverBranching, vertices::Vector{Int}, g::SimpleGraph, use_rv::Bool) where N
    return setcover_strategy(tbl, vertices, g, max_itr = strategy.max_itr, use_rv)
end

function BitBasis.BitStr(sv::StaticBitVector)
    @assert length(sv.data) == 1 "bit string too long!"
    return BitBasis.BitStr{length(sv), Int}(sv.data[1])
end