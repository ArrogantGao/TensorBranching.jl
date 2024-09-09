"""
    SolverConfig{ST<:AbstractMISSolver, BT<:AbstractBranching, MT<:AbstractMeasure, VT<:AbstractVertexSelector, FT<:AbstractTruthFilter}

The configuration for the optimal branching MIS solver.

# Fields
- `table_solver::ST`: the MIS solver, e.g., `TensorNetworkSolver()`
- `branching_strategy::BT`: the branching strategy, e.g., `NaiveBranching()`, `SetCoverBranching()`
- `measure::MT`: the measure for the branching strategy, e.g., `NumOfVertices()`, `D3Measure()`
- `vertex_selector::VT`: the vertex selector, e.g., `MinBoundarySelector(2)`, `ManualSelector([1, 2, 3])`
- `table_filter::FT`: the truth filter, e.g., `EnvFilter()`, `NoFilter()`
"""
Base.@kwdef struct SolverConfig{ST<:AbstractMISSolver, BT<:AbstractBranching, MT<:AbstractMeasure, VT<:AbstractVertexSelector, FT<:AbstractTruthFilter}
    table_solver::ST = TensorNetworkSolver()
    branching_strategy::BT = NaiveBranching()
    measure::MT = NumOfVertices()
    vertex_selector::VT = MinBoundarySelector(2)
    table_filter::FT = EnvFilter()
end

"""
    solve_mis(g::SimpleGraph, config::SolverConfig)

Solve the maximum independent set problem for the given graph and configuration.

# Arguments
- `g::SimpleGraph`: the input graph
- `config::SolverConfig`: the solver configuration
"""
function solve_mis(g::SimpleGraph, config::SolverConfig)
    mis = mis_solver(g, config)
    return mis.mis
end

"""
    count_mis(g::SimpleGraph, config::SolverConfig)

Count the number of MIS for the given graph and configuration.

# Arguments
- `g::SimpleGraph`: the input graph
- `config::SolverConfig`: the solver configuration
"""
function count_mis(g::SimpleGraph, config::SolverConfig)
    mis = mis_solver(g, config)
    return mis.mis, mis.count
end

"""
    mis_solver(g::SimpleGraph, config::SolverConfig)

Solves the maximum independent set (MIS) problem for a given graph `g` using a recursive algorithm.

# Arguments
- `g::SimpleGraph`: The input graph.
- `config::SolverConfig`: Configuration options for the solver, see `SolverConfig`.

# Returns
- `max_mis`: The size of the maximum independent set.

"""
function mis_solver(g::SimpleGraph, config::SolverConfig)
    dg = degree(g)
    if nv(g) == 0 || nv(g) == 1
        return CountingMIS(nv(g))
    elseif (0 ∈ dg) || (1 ∈ dg)
        v = findfirst(x -> (x==0)||(x==1), dg)
        return 1 + mis_solver(induced_subgraph(g, setdiff(1:nv(g), v ∪ neighbors(g, v)))[1], config)
    elseif (2 ∈ dg)
        v = findfirst(x -> (x==2), dg)
        return folding(g, v, config)
    elseif maximum(dg) ≥ 6
        v = argmax(dg)
        return max(1 + mis_solver(induced_subgraph(g, setdiff(1:nv(g), v ∪ neighbors(g, v)))[1], config), mis_solver(induced_subgraph(g, setdiff(1:nv(g), v))[1], config))
    else
        vertices = select_vertex(g, config.vertex_selector)
        branches = optimal_branches(g, vertices, config.branching_strategy; config.measure, config.table_filter, config.table_solver)

        mis_count = Vector{CountingMIS}(undef, length(branches.branches))
        
        Threads.@threads for i in 1:length(branches.branches)
            rvs = branches.branches[i].vertices_removed
            gi = copy(g)
            rem_vertices!(gi, rvs)
            mis_count[i] = mis_solver(gi, config) + branches.branches[i].mis
        end
        max_mis = maximum(mis_count)

        return max_mis
    end
end

"""
    optimal_branches(g::SimpleGraph, vertices::AbstractVector, strategy::AbstractBranching; measure::AbstractMeasure = NumOfVertices(), filter::AbstractTruthFilter = NoFilter())
    
Find the optimal branches for the given graph and vertices.

# Arguments
- `g::SimpleGraph`: the input graph
- `vertices::AbstractVector`: the vertices to be considered
- `strategy::AbstractBranching`: the branching strategy, e.g., `NaiveBranching()`, `SetCoverBranching()`

# Keyword Arguments
- `measure::AbstractMeasure`: the measure for the branching strategy, e.g., `NumOfVertices()`, `D3Measure()`
- `filter::AbstractTruthFilter`: the filter for the branching strategy, e.g., `EnvFilter()`, `NoFilter()`
- `table_solver`: the solver configuration, e.g., `TensorNetworkSolver()`
"""
function optimal_branches(g::SimpleGraph, vertices::AbstractVector, strategy::AbstractBranching;
            measure::AbstractMeasure = NumOfVertices(),
            table_filter::AbstractTruthFilter = NoFilter(),
            table_solver = TensorNetworkSolver()
        )
    # open vertices and induced subgraph
    openvertices = open_vertices(g, vertices)
    subg, vmap = induced_subgraph(g, vertices)

    # reduced alpha configurations
    tbl = reduced_alpha_configs(table_solver, subg, Int[findfirst(==(v), vertices) for v in openvertices])
    # filter out some configurations using the environmental information
    filtered_tbl = filt(g, vertices, openvertices, tbl, table_filter)

    # implement the branching strategy
    branches = impl_strategy(g, vertices, filtered_tbl, strategy, measure)
    return branches
end


function folding(g::SimpleGraph, v::Int, config::SolverConfig)
    @assert degree(g, v) == 2
    a, b = neighbors(g, v)
    if has_edge(g, a, b)
        return 1 + mis_solver(induced_subgraph(g, setdiff(1:nv(g), [v, a, b]))[1], config)
    else
        # apply the graph rewrite rule
        add_vertex!(g)
        nn = open_neighbors(g, [v, a, b])
        for n in nn
            add_edge!(g, nv(g), n)
        end
        return 1 + mis_solver(induced_subgraph(g, setdiff(1:nv(g), [v, a, b]))[1], config)
    end
end
