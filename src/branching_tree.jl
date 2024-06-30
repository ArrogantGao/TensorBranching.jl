# a tool to generate the branching process as a tree

using AbstractTrees

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

function branching_tree(g::SimpleGraph, strategy::BranchingStrategy, kneighbor::Int)
    tree = _branching_tree(g, strategy, kneighbor)
    branch_num = length(collect(Leaves(tree)))
    return tree, branch_num
end

function _branching_tree(g::SimpleGraph, strategy::BranchingStrategy, kneighbor::Int)
    if nv(g) == 0 || nv(g) == 1
        return BranchingNode(g)
    elseif maximum(degree(g)) ≥ 6
        v = argmax(degree(g))
        return BranchingNode(g, children=[_branching_tree(induced_subgraph(g, setdiff(1:nv(g), v ∪ neighbors(g, v)))[1], strategy, kneighbor), _branching_tree(induced_subgraph(g, setdiff(1:nv(g), v))[1], strategy, kneighbor)], removed=[v ∪ neighbors(g, v), [v]])
    end
    
    vertices, openvertices, dnf = optimal_branching_dnf(g, strategy, kneighbor)
    @assert !isempty(vertices)
    
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
        add_child!(root, _branching_tree(gi, strategy, kneighbor))
    end

    return root
end