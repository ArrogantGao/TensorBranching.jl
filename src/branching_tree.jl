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
    branching_tree(g::SimpleGraph, config::SolverConfig)

Constructs a branching tree based on the given graph `g` using the specified `config`.
"""
function branching_tree(g::SimpleGraph, config::SolverConfig)
    tree = _branching_tree(g, config)
    branch_num = length(collect(Leaves(tree)))
    return tree, branch_num
end

function _branching_tree(g::SimpleGraph, config::SolverConfig)
    dg = degree(g)
    if nv(g) == 0 || nv(g) == 1
        return BranchingNode(g)
    elseif (0 ∈ dg) || (1 ∈ dg)
        v = findfirst(x -> (x==0)||(x==1), dg)
        return BranchingNode(g, children=[_branching_tree(induced_subgraph(g, setdiff(1:nv(g), v ∪ neighbors(g, v)))[1], config)])
    elseif (2 ∈ dg)
        v = findfirst(x -> (x==2), dg)
        a, b = neighbors(g, v)
        if has_edge(g, a, b)
            return BranchingNode(g, children=[_branching_tree(induced_subgraph(g, setdiff(1:nv(g), [v, a, b]))[1], config)], removed=[[v, a, b]])
        else
            # apply the graph rewrite rule
            add_vertex!(g)
            nn = open_neighbors(g, [v, a, b])
            for n in nn
                add_edge!(g, nv(g), n)
            end
            return BranchingNode(g, children=[_branching_tree(induced_subgraph(g, setdiff(1:nv(g), [v, a, b]))[1], config)], removed=[[v, a, b]])
        end
    elseif maximum(dg) ≥ 6
        v = argmax(dg)
        return BranchingNode(g, children=[_branching_tree(induced_subgraph(g, setdiff(1:nv(g), v ∪ neighbors(g, v)))[1], config), _branching_tree(induced_subgraph(g, setdiff(1:nv(g), v))[1], config)], removed=[v ∪ neighbors(g, v), [v]])
    else
        root = BranchingNode(g)
        vertices = select_vertex(g, config.vertex_selector)
        branches = optimal_branches(g, vertices, config.branching_strategy; config.measure, config.table_filter, config.table_solver)
        
        for i in 1:length(branches.branches)
            rvs = branches.branches[i].vertices_removed
            gi = copy(g)
            rem_vertices!(gi, rvs)
            add_removed!(root, unique(rvs))
            add_child!(root, _branching_tree(gi, config))
        end

        return root
    end
end