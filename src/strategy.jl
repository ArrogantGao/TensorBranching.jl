"""
    abstract type AbstractBranching

Abstract type representing a branching strategy. For each branching strategy, a corresponding `impl_strategy` function should be implemented.
"""
abstract type AbstractBranching end


"""
    abstract type AbstractMeasure

Abstract type representing a measurement of the size of given graphs. For each measurement, a corresponding `measure` function should be implemented.
"""
abstract type AbstractMeasure end

size_reduced(g::SimpleGraph, vertices::Vector{Int}, clause::Clause{T}, m::AbstractMeasure) where{T} = measure(g, m) - measure(induced_subgraph(g, setdiff(1:nv(g), removed_vertices(vertices, g, clause)))[1], m)
size_reduced(g::SimpleGraph, vertices_removed::Vector{Int}, m::AbstractMeasure) = measure(g, m) - measure(induced_subgraph(g, setdiff(1:nv(g), vertices_removed))[1], m)

"""
    abstract type AbstractVertexSelector

Abstract type representing a vertex selector, which selects a subset of vertices from a given graph for branching. A corresponding `select_vertex` function should be implemented.
"""
abstract type AbstractVertexSelector end


"""
    abstract type AbstractTruthFilter

Abstract type representing a truth filter.

This type serves as a base for defining different truth filter strategies.
"""
abstract type AbstractTruthFilter end

"""
    abstract type AbstractSetCoverSolver

Abstract type representing a solver for the set cover problem, including linear programming and integer programming solvers.
"""
abstract type  AbstractSetCoverSolver end

struct LPSetCoverSolver <: AbstractSetCoverSolver end
struct IPSetCoverSolver <: AbstractSetCoverSolver end

"""
    struct NaiveBranching <: AbstractBranching
A struct representing the NaiveBranching branching strategy.
"""
struct NaiveBranching <: AbstractBranching end

function impl_strategy(g::SimpleGraph, vertices::Vector{Int}, tbl::BranchingTable{INT}, ::NaiveBranching, measure::AbstractMeasure) where INT
    return Branches([Branch(Clause(bmask(INT, 1:nbits(tbl)), first(x)), vertices, g, measure) for x in tbl.table])
end

"""
    struct SetCoverBranching <: AbstractBranching

A struct representing a branching strategy for set cover problems.

# Fields
- `max_itr::Int`: The maximum number of iterations.

# Constructors
- `SetCoverBranching()`: Constructs a `SetCoverBranching` object with a default value of `2` for `max_itr`.
- `SetCoverBranching(max_itr::Int)`: Constructs a `SetCoverBranching` object with the specified `max_itr` value.
- `verbose::Bool`: A boolean indicating whether to print the output.

"""
Base.@kwdef struct SetCoverBranching <: AbstractBranching 
    max_itr::Int = 5
    solver::AbstractSetCoverSolver = IPSetCoverSolver()
    verbose::Bool = false
end

function impl_strategy(g::SimpleGraph, vertices::Vector{Int}, tbl::BranchingTable{INT}, strategy::SetCoverBranching, measure::AbstractMeasure) where INT
    return setcover_strategy(tbl, vertices, g, strategy.max_itr, measure, strategy.solver, strategy.verbose)
end

struct NumOfVertices <: AbstractMeasure end # each vertex is counted as 1

measure(g::SimpleGraph, ::NumOfVertices) = nv(g)

"""
    D3Measure <: AbstractMeasure

A measure of complexity by counting the number of vertices with degree at least 3.
"""
struct D3Measure <: AbstractMeasure end # n = sum max{d - 2, 0}

function measure(g::SimpleGraph, ::D3Measure)
    if nv(g) == 0
        return 0
    else
        dg = degree(g)
        return Int(sum(max(d - 2, 0) for d in dg))
    end
end

"""
    struct MinBoundarySelector <: AbstractVertexSelector

The `MinBoundarySelector` struct represents a strategy for selecting a subgraph with the minimum number of open vertices by k-layers of neighbors.

# Fields
- `k::Int`: The number of layers of neighbors to consider when selecting the subgraph.

"""
struct MinBoundarySelector <: AbstractVertexSelector
    k::Int # select the subgraph with minimum open vertices by k-layers of neighbors
end

function select_vertex(g::SimpleGraph{Int}, vertex_select::MinBoundarySelector)
    kneighbor = vertex_select.k

    local vs_min
    novs_min = nv(g)
    for v in 1:nv(g)
        vs, ovs = neighbor_cover(g, v, kneighbor)
        if length(ovs) < novs_min
            vs_min = vs
            novs_min = length(ovs)
        end
    end
    return vs_min
end

struct NoFilter <: AbstractTruthFilter end

filt(g::SimpleGraph, vertices::Vector{Int}, openvertices::Vector{Int}, tbl::BranchingTable{INT}, ::NoFilter) where{INT} = tbl

"""
    struct EnvFilter <: AbstractTruthFilter
A struct representing an environment filter, using the information of the environment to filter out some configurations which can not reach the optimal mis solution.
"""
struct EnvFilter <: AbstractTruthFilter end

# consider two different branching rule (A, and B) applied on the same set of vertices, with open vertices ovs.
# the neighbors of 1 vertices in A is label as NA1, and the neighbors of 1 vertices in B is label as NB1, and the pink_block is the set of vertices that are not in NB1 but in NA1.
# once mis(A) + mis(pink_block) ≤ mis(B), then A is not a good branching rule, and should be removed.
function filt(g::SimpleGraph, vertices::Vector{Int}, openvertices::Vector{Int}, tbl::BranchingTable{INT}, ::EnvFilter) where INT
    ns = neighbors(g, vertices)
    so = Set(openvertices)

    new_table = Vector{Vector{INT}}()

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
                if (count_ones(tbl.table[i][1]) + mis_pink ≤ count_ones(tbl.table[j][1])) && (!iszero(mis_pink))
                    flag = false
                    break
                end
            end
        end
        if flag
            push!(new_table, tbl.table[i])
        end
    end

    return BranchingTable(nbits(tbl), new_table)
end

"""
    struct Branch

A struct representing a branching strategy.

# Fields
- `vertices_removed::Vector{Int}`: A vector of integers representing the vertices removed in the branching strategy.
- `mis::Int`: An integer representing the maximum independent set (MIS) size of the branching strategy.

"""
struct Branch
    vertices_removed::Vector{Int}
    mis::Int
end

function Branch(clause::Clause{INT}, vertices::Vector{Int}, g::SimpleGraph, measure::AbstractMeasure) where {INT}
    vertices_removed = removed_vertices(vertices, g, clause)
    return Branch(vertices_removed, count_ones(clause.val))
end

"""
    struct Branches

A struct representing a collection of branches.

# Fields
- `branches::Vector{Branch}`: A vector of `Branch` objects.

"""
struct Branches
    branches::Vector{Branch}
end

Base.:(==)(b1::Branch, b2::Branch) = (Set(b1.vertices_removed) == Set(b2.vertices_removed)) && (b1.mis == b2.mis)
Base.:(==)(b1::Branches, b2::Branches) = ((b1.branches) == (b2.branches))

"""
    effective_γ(branches::Branches, g::SimpleGraph, measure::AbstractMeasure)

Compute the effective γ value for a set of branches in a graph under a given measure.

This function calculates the effective γ value for a given set of branches in a graph. The effective γ value is a measure of complexity that takes into account the size reduction achieved by removing vertices from the graph.

# Arguments
- `branches::Branches`: The set of branches for which to calculate the effective γ value.
- `g::SimpleGraph`: The graph for which to calculate the effective γ value.
- `measure::AbstractMeasure`: The measure used to calculate the size reduction.

# Returns
- The effective γ value for the set of branches.

"""
function effective_γ(branches::Branches, g::SimpleGraph, measure::AbstractMeasure)
    return complexity([size_reduced(g, b.vertices_removed, measure) for b in branches.branches])
end
