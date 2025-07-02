# kernlize the graph, via reduction
using OptimalBranching.OptimalBranchingMIS: reduce_graph, ReductionResult

function kernelize(g::SimpleGraph, weights::VT, reducer::AbstractReducer; verbose::Int = 0, vmap::Vector{Int} = collect(1:nv(g))) where VT
    (verbose ≥ 2) && (@info "kernelizing graph: $(nv(g)) vertices, $(ne(g)) edges")
    r = zero(eltype(weights))
    while true
        res = reduce_graph(g, weights, reducer) # res = (g_new, r_new, vmap_new)
        vmap = vmap[res.vmap]
        weights = res.weights
        r += res.r
        if g == res.g
            (verbose ≥ 2) && (@info "kernelized graph: $(nv(g)) vertices, $(ne(g)) edges")
            return ReductionResult(g, weights, r, vmap)
        end
        g = res.g
    end
end