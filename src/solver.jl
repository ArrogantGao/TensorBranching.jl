function missolve(g::SimpleGraph; show_count = false, strategy::AbsractBranching = NaiveBranching(), measurement::AbstractMeasurement = NaiveMeasure(), vertex_select::AbstractVertexSelector = MinBoundSelector(2), filter::AbstractTruthFilter = EnvFilter())
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
        vertices, openvertices, dnf, γ = optimal_branching_dnf(g, strategy, measurement, vertex_select, filter)
        # @assert !isempty(vertices)

        mis_count = Vector{CountingMIS}(undef, length(dnf.clauses))
        
        @threads for i in 1:length(dnf.clauses)
            clause = dnf.clauses[i]
            rvs = removed_vertices(vertices, g, clause)
            gi = copy(g)
            rem_vertices!(gi, rvs)
            # @assert !isempty(rvs)
            mis_count[i] = mis_solver(gi, strategy, measurement, vertex_select, filter) + count_ones(clause.val)
        end
        max_mis = maximum(mis_count)

        return max_mis
    end
end

function optimal_branching_dnf(g::SimpleGraph, strategy::AbsractBranching, measurement::AbstractMeasurement, vertex_select::AbstractVertexSelector, filter::AbstractTruthFilter)
    
    vertices, openvertices = select_vertex(g, vertex_select)

    subg, vmap = induced_subgraph(g, vertices)
    tbl = reduced_alpha_configs(subg, Int[findfirst(==(v), vertices) for v in openvertices])

    filtered_tbl = filt(g, vertices, openvertices, tbl, filter)
    opt_dnf, γ = impl_strategy(g, vertices, filtered_tbl, strategy, measurement)

    return vertices, openvertices, opt_dnf, γ
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