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
    return Branches([Branch(Clause(bmask(BitStr{N, Int}, 1:N), BitStr(first(x))), vertices, g, measurement) for x in tbl.table])
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

struct NumOfVertices <: AbstractMeasurement end # each vertex is counted as 1

measure(g::SimpleGraph, ::NumOfVertices) = nv(g)

function size_reduced(g::SimpleGraph, vertices::Vector{Int}, clause::Clause{N, T}, measurement::NumOfVertices) where{N, T}
    return length(removed_vertices(vertices, g, clause))
end

struct NumOfDegree <: AbstractMeasurement end # n = sum max{d - 2, 0}

function measure(g::SimpleGraph, ::NumOfDegree)
    if nv(g) == 0
        return 0
    else
        dg = degree(g)
        return sum(max(d - 2, 0) for d in dg)
    end
end

function size_reduced(g::SimpleGraph, vertices::Vector{Int}, clause::Clause{N, T}, measurement::NumOfDegree) where{N, T}
    return measure(g, measurement) - measure(induced_subgraph(g, setdiff(1:nv(g), removed_vertices(vertices, g, clause)))[1], measurement)
end

struct MinBoundarySelector <: AbstractVertexSelector
    k::Int # select the subgraph with minimum open vertices by k-layers of neighbors
end

function select_vertex(g::SimpleGraph{Int}, vertex_select::MinBoundarySelector)

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

struct ManualSelector <: AbstractVertexSelector
    vertices::Vector{Int}
end

function select_vertex(g::SimpleGraph, vertex_select::ManualSelector)
    vertices = vertex_select.vertices
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

struct Branch{T}
    vertices_removed::Vector{Int}
    size_reduced::T
    mis::Int
end

function Branch(clause::Clause{N, C}, vertices::Vector{Int}, g::SimpleGraph, measurement::AbstractMeasurement) where {N, C}
    vertices_removed = removed_vertices(vertices, g, clause)
    sr = size_reduced(g, vertices, clause, measurement)
    return Branch(vertices_removed, sr, count_ones(clause.val))
end

struct Branches{T}
    branches::Vector{Branch{T}}
end

Base.:(==)(b1::Branch, b2::Branch) = (Set(b1.vertices_removed) == Set(b2.vertices_removed)) && (b1.size_reduced == b2.size_reduced) && (b1.mis == b2.mis)
Base.:(==)(b1::Branches, b2::Branches) = ((b1.branches) == (b2.branches))

function effective_γ(branches::Branches{T}) where{T}
    return complexity([branches.branches[i].size_reduced for i in 1:length(branches.branches)])
end