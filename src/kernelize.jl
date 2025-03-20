# kernlize the graph, via reduction
using OptimalBranching.OptimalBranchingMIS: reduce_graph

function kernelize(g::SimpleGraph, reducer::AbstractReducer; verbose::Int = 0, vmap::Vector{Int} = collect(1:nv(g)))
    (verbose ≥ 2) && (@info "kernelizing graph: $(nv(g)) vertices, $(ne(g)) edges")
    r = 0
    while true
        res = reduce_graph(g, reducer) # res = (g_new, r_new, vmap_new)
        vmap = vmap[res[3]]
        r += res[2]
        if g == res[1]
            (verbose ≥ 2) && (@info "kernelized graph: $(nv(g)) vertices, $(ne(g)) edges")
            return (g, r, vmap)
        end
        g = res[1]
    end
end
