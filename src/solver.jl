function missolve(g::SimpleGraph; show_count = false, strategy::AbsractBranching = NaiveBranching(), measurement::AbstractMeasurement = NumOfVertices(), vertex_select::AbstractVertexSelector = MinBoundarySelector(2), filter::AbstractTruthFilter = EnvFilter())
    mis = mis_solver(g, strategy, measurement, vertex_select, filter)
    return show_count ? (mis.mis, mis.count) : mis.mis
end

function mis_solver(g::SimpleGraph, strategy::AbsractBranching, measurement::AbstractMeasurement, vertex_select::AbstractVertexSelector, filter::AbstractTruthFilter)
    dg = degree(g)
    if nv(g) == 0 || nv(g) == 1
        return CountingMIS(nv(g))
    elseif (0 ∈ dg) || (1 ∈ dg)
        v = findfirst(x -> (x==0)||(x==1), dg)
        return 1 + mis_solver(induced_subgraph(g, setdiff(1:nv(g), v ∪ neighbors(g, v)))[1], strategy,measurement, vertex_select, filter)
    elseif (2 ∈ dg)
        v = findfirst(x -> (x==2), dg)
        return folding(g, v, strategy, measurement, vertex_select, filter)
    elseif maximum(dg) ≥ 6
        v = argmax(dg)
        return max(1 + mis_solver(induced_subgraph(g, setdiff(1:nv(g), v ∪ neighbors(g, v)))[1], strategy, measurement, vertex_select, filter), mis_solver(induced_subgraph(g, setdiff(1:nv(g), v))[1], strategy, measurement, vertex_select, filter))
    else
        vertices = select_vertex(g, vertex_select)
        branches = optimal_branches(g, vertices, strategy; measurement, filter)

        mis_count = Vector{CountingMIS}(undef, length(branches.branches))
        
        @threads for i in 1:length(branches.branches)
            rvs = branches.branches[i].vertices_removed
            gi = copy(g)
            rem_vertices!(gi, rvs)
            # @assert !isempty(rvs)
            mis_count[i] = mis_solver(gi, strategy, measurement, vertex_select, filter) + branches.branches[i].mis
        end
        max_mis = maximum(mis_count)

        return max_mis
    end
end

"""
    optimal_branches(g::SimpleGraph, vertices::AbstractVector, strategy::AbsractBranching; measurement::AbstractMeasurement = NumOfVertices(), filter::AbstractTruthFilter = NoFilter())
    
Find the optimal branches for the given graph and vertices.

# Arguments
- `g::SimpleGraph`: the input graph
- `vertices::AbstractVector`: the vertices to be considered
- `strategy::AbsractBranching`: the branching strategy, e.g., `NaiveBranching()`, `SetCoverBranching()`

# Keyword Arguments
- `measurement::AbstractMeasurement`: the measurement for the branching strategy, e.g., `NumOfVertices()`, `NumOfDegree()`
- `filter::AbstractTruthFilter`: the filter for the branching strategy, e.g., `EnvFilter()`, `NoFilter()`
"""
function optimal_branches(g::SimpleGraph, vertices::AbstractVector, strategy::AbsractBranching;
            measurement::AbstractMeasurement = NumOfVertices(),
            filter::AbstractTruthFilter = NoFilter(),
        )
    # open vertices and induced subgraph
    openvertices = open_vertices(g, vertices)
    subg, vmap = induced_subgraph(g, vertices)

    # reduced alpha configurations
    tbl = reduced_alpha_configs(subg, Int[findfirst(==(v), vertices) for v in openvertices])
    # filter out some configurations using the environmental information
    filtered_tbl = filt(g, vertices, openvertices, tbl, filter)

    # implement the branching strategy
    branches = impl_strategy(g, vertices, filtered_tbl, strategy, measurement)
    return branches
end


function folding(g::SimpleGraph, v::Int, strategy::AbsractBranching, measurement::AbstractMeasurement, vertex_select::AbstractVertexSelector, filter::AbstractTruthFilter)
    @assert degree(g, v) == 2
    a, b = neighbors(g, v)
    if has_edge(g, a, b)
        return 1 + mis_solver(induced_subgraph(g, setdiff(1:nv(g), [v, a, b]))[1], strategy, measurement, vertex_select, filter)
    else
        # apply the graph rewrite rule
        add_vertex!(g)
        nn = open_neighbors(g, [v, a, b])
        for n in nn
            add_edge!(g, nv(g), n)
        end
        return 1 + mis_solver(induced_subgraph(g, setdiff(1:nv(g), [v, a, b]))[1], strategy, measurement, vertex_select, filter)
    end
end