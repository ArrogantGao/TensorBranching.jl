abstract type AbsractBranching end
abstract type AbstractMeasurement end
abstract type AbstractVertexSelector end
abstract type AbstractTruthFilter end

"""
    struct NaiveBranching <: AbsractBranching
A struct representing the NaiveBranching branching strategy.
"""
struct NaiveBranching <: AbsractBranching end

function impl_strategy(g::SimpleGraph, vertices::Vector{Int}, tbl::BranchingTable{N}, strategy::NaiveBranching, measurement::AbstractMeasurement) where{N}
    return DNF([Clause(bmask(BitStr{N, Int}, 1:N), BitStr(first(x))) for x in tbl.table])
end

"""
    struct SetCoverBranching <: AbsractBranching

A struct representing a branching strategy for set cover problems.

# Fields
- `max_itr::Int`: The maximum number of iterations.

# Constructors
- `SetCoverBranching()`: Constructs a `SetCoverBranching` object with a default value of `2` for `max_itr`.
- `SetCoverBranching(max_itr::Int)`: Constructs a `SetCoverBranching` object with the specified `max_itr` value.

"""
struct SetCoverBranching <: AbsractBranching 
    max_itr::Int
    SetCoverBranching() = new(3)
    SetCoverBranching(max_itr::Int) = new(max_itr)
end

function impl_strategy(g::SimpleGraph, vertices::Vector{Int}, tbl::BranchingTable{N}, strategy::SetCoverBranching, measurement::AbstractMeasurement) where{N}
    return setcover_strategy(tbl, vertices, g, strategy.max_itr, measurement)
end

struct NaiveMeasure <: AbstractMeasurement end # each vertex is counted as 1

measure(g::SimpleGraph, ::NaiveMeasure) = nv(g)

function nv_removed(g::SimpleGraph, vertices::Vector{Int}, clause::Clause{N, T}, measurement::NaiveMeasure) where{N, T}
    return length(removed_vertices(vertices, g, clause))
end

struct D3Measure <: AbstractMeasurement end # n = sum max{d - 2, 0}

function measure(g::SimpleGraph, ::D3Measure)
    if nv(g) == 0
        return 0
    else
        dg = degree(g)
        return sum(max(d - 2, 0) for d in dg)
    end
end

function nv_removed(g::SimpleGraph, vertices::Vector{Int}, clause::Clause{N, T}, measurement::D3Measure) where{N, T}
    return measure(g, measurement) - measure(induced_subgraph(g, setdiff(1:nv(g), removed_vertices(vertices, g, clause)))[1], measurement)
end

struct MinBoundSelector <: AbstractVertexSelector
    k::Int # select the subgraph with minimum open vertices by k-layers of neighbors
end

function select_vertex(g::SimpleGraph{Int}, vertex_select::MinBoundSelector)

    kneighbor = vertex_select.k

    vs_min = Int[]
    ovs_min = Int[1:nv(g)...]
    for v in 1:nv(g)
        vs, ovs = neighbor_cover(g, v, kneighbor)
        if length(ovs) < length(ovs_min)
            vs_min = vs
            ovs_min = ovs
        end
    end
    vertices = vs_min
    openvertices = ovs_min

    return vertices, openvertices
end

struct ManulSelector <: AbstractVertexSelector
    vertices::Vector{Int}
end

function select_vertex(g::SimpleGraph, vertex_select::ManulSelector)
    vertices = vertex_select.vertex_select
    openvertices = open_vertices(g, vertices)
    return vertices, openvertices
end

struct NoFilter <: AbstractTruthFilter end

filt(g::SimpleGraph, vertices::Vector{Int}, openvertices::Vector{Int}, tbl::BranchingTable{N, C}, ::NoFilter) where{N, C} = tbl

struct EnvFilter <: AbstractTruthFilter end

function filt(g::SimpleGraph, vertices::Vector{Int}, openvertices::Vector{Int}, tbl::BranchingTable{N, C}, ::EnvFilter) where{N, C}
    ns = neighbors(g, vertices)
    so = Set(openvertices)

    new_table = Vector{Vector{StaticBitVector{N, C}}}()

    open_vertices_1 = [Int[] for i in 1:length(tbl.table)]
    neibs_0 = Set{Int}[]
    for i in 1:length(tbl.table)
        row = tbl.table[i]
        x = row[1]
        for n in 1:length(x)
            if (x[n] == 1) && (vertices[n] ∈ so)
                push!(open_vertices_1[i], vertices[n])
            end
        end
        push!(neibs_0, setdiff(ns, neighbors(g, open_vertices_1[i]) ∩ ns))
    end

    for i in 1:length(tbl.table)
        flag = true
        for j in 1:length(tbl.table)
            if i != j
                pink_block = setdiff(neibs_0[i], neibs_0[j])
                sg_pink, sg_vec = induced_subgraph(g, collect(pink_block))
                mis_pink = mis2(EliminateGraph(sg_pink))
                if (count_ones(BitStr(tbl.table[i][1])) + mis_pink ≤ count_ones(BitStr(tbl.table[j][1]))) && (!iszero(mis_pink))
                    flag = false
                    break
                end
            end
        end
        if flag
            push!(new_table, tbl.table[i])
        end
    end

    return BranchingTable(new_table)
end