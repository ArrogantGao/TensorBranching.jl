# a tool to generate the branching process as a tree

using AbstractTrees

"""
    mutable struct BranchingNode

A mutable struct representing a node in a branching tree.

# Fields
- `children::Vector{BranchingNode}`: The children nodes of the current node.
- `graph::SimpleGraph{Int}`: The graph associated with the current node.
- `removed::Vector{Vector{Int}}`: The removed edges from the graph.

# Constructor
- `BranchingNode(graph::SimpleGraph{Int}; children::Vector{BranchingNode} = Vector{BranchingNode}(), removed::Vector{Vector{Int}} = Vector{Vector{Int}}())`: Constructs a new `BranchingNode` object.

"""
mutable struct BranchingNode
    children::Vector{BranchingNode}
    graph::SimpleGraph{Int}
    removed::Vector{Vector{Int}}
    BranchingNode(graph::SimpleGraph{Int}; children::Vector{BranchingNode} = Vector{BranchingNode}(), removed::Vector{Vector{Int}} = Vector{Vector{Int}}()) = new(children, graph, removed)
end

AbstractTrees.children(node::BranchingNode) = node.children
AbstractTrees.printnode(io::IO, node::BranchingNode) = print(io, "G{$(nv(node.graph))}")
AbstractTrees.nodevalue(node::BranchingNode) = node.graph

isleaf(node::BranchingNode) = isempty(node.children)
Base.show(io::IO, tree::BranchingNode) = print_tree(io, tree)

function add_child!(node::BranchingNode, child::BranchingNode)
    push!(node.children, child)
    return node
end

function add_removed!(node::BranchingNode, removed::Vector{Int})
    push!(node.removed, removed)
    return node
end

"""
    branching_tree(g::SimpleGraph, strategy::AbsractBranching, kneighbor::Int, use_rv::Bool)

Constructs a branching tree based on the given graph `g` using the specified `strategy`.
The `kneighbor` parameter determines the number of neighbors to consider during branching.
If `use_rv` is `true`, number of vertex removed is used; otherwise, only count the clause.

# Arguments
- `g::SimpleGraph`: The input graph.
- `strategy::AbsractBranching`: The branching strategy to use.
- `kneighbor::Int`: The number of neighbors to consider during branching.
- `use_rv::Bool`: Whether to use number of vertex removed.

# Returns
- `tree`: The constructed branching tree.
- `branch_num`: The number of branches in the tree.
"""
function branching_tree(g::SimpleGraph, strategy::AbsractBranching, kneighbor::Int, use_rv::Bool)
    tree = _branching_tree(g, strategy, kneighbor, use_rv)
    branch_num = length(collect(Leaves(tree)))
    return tree, branch_num
end

function _branching_tree(g::SimpleGraph, strategy::AbsractBranching, kneighbor::Int, use_rv::Bool)
    dg = degree(g)
    if nv(g) == 0 || nv(g) == 1
        return BranchingNode(g)
    elseif (0 ∈ dg) || (1 ∈ dg)
        v = findfirst(x -> (x==0)||(x==1), dg)
        return BranchingNode(g, children=[_branching_tree(induced_subgraph(g, setdiff(1:nv(g), v ∪ neighbors(g, v)))[1], strategy, kneighbor, use_rv)])
    elseif maximum(degree(g)) ≥ 6
        v = argmax(degree(g))
        return BranchingNode(g, children=[_branching_tree(induced_subgraph(g, setdiff(1:nv(g), v ∪ neighbors(g, v)))[1], strategy, kneighbor, use_rv), _branching_tree(induced_subgraph(g, setdiff(1:nv(g), v))[1], strategy, kneighbor, use_rv)], removed=[v ∪ neighbors(g, v), [v]])
    end
    
    vertices, openvertices, dnf = optimal_branches(g, strategy, kneighbor, use_rv)
    # @assert !isempty(vertices)
    
    root = BranchingNode(g)
    for i in 1:length(dnf.clauses)
        clause = dnf.clauses[i]
        removed_vertices = Int[]
        for (k, v) in enumerate(vertices)
            if readbit(clause.mask, k) == 1
                push!(removed_vertices, v)
                if readbit(clause.val, k) == 1
                    append!(removed_vertices, neighbors(g, v))
                end
            end
        end
        gi = copy(g)
        rem_vertices!(gi, removed_vertices)
        @assert !isempty(removed_vertices)
        add_removed!(root, unique(removed_vertices))
        add_child!(root, _branching_tree(gi, strategy, kneighbor, use_rv))
    end

    return root
end