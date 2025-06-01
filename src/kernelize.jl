# kernlize the graph, via reduction
using OptimalBranching.OptimalBranchingMIS: reduce_graph, tn_reduce_graph

function kernelize(g::SimpleGraph, reducer::AbstractReducer; verbose::Int = 0, vmap::Vector{Int} = collect(1:nv(g)))
    (verbose ≥ 2) && (@info "kernelizing graph: $(nv(g)) vertices, $(ne(g)) edges")
    r = 0
    while true
        res = reduce_graph(g, reducer) # res = (g_new, r_new, vmap_new)
        vmap = vmap[res[3]]
        r += res[2]
        if g == res[1]
            (verbose ≥ 2) && (@info "kernelized graph: $(nv(g)) vertices, $(ne(g)) edges")
            return (g, r, vmap, reducer)
        end
        g = res[1]
    end
end

function kernelize(g::SimpleGraph, weights::Vector, reducer::AbstractReducer; verbose::Int = 0, vmap::Vector{Int} = collect(1:nv(g)))
    (verbose ≥ 2) && (@info "kernelizing graph: $(nv(g)) vertices, $(ne(g)) edges")
    r = 0
    while true
        res = reduce_graph(g, weights, reducer) # res = (g_new, weights_new, r_new, vmap_new)
        vmap = vmap[res[4]]
        r += res[3]
        if g == res[1]
            (verbose ≥ 2) && (@info "kernelized graph: $(nv(g)) vertices, $(ne(g)) edges")
            return (g, weights, r, vmap, reducer)
        end
        g = res[1]
        weights = res[2]
    end
end

# this will change the reducer
function kernelize(g::SimpleGraph, reducer::TensorNetworkReducer; verbose::Int = 0, vmap::Vector{Int} = collect(1:nv(g)))
    (verbose ≥ 2) && (@info "kernelizing graph: $(nv(g)) vertices, $(ne(g)) edges")
    r = 0

    vmap_0 = vmap
    reducer = deepcopy(reducer)

    while true
        res = reduce_graph(g, reducer, vmap_0 = vmap_0) # res = (g_new, r_new, vmap_new)
        vmap_0 = res[3]
        vmap = vmap[res[3]]
        r += res[2]
        if g == res[1]
            (verbose ≥ 2) && (@info "kernelized graph: $(nv(g)) vertices, $(ne(g)) edges")
            return (g, r, vmap, reducer)
        end
        g = res[1]
    end
end

# this will change the reducer
function kernelize(g::SimpleGraph, weights::Vector, reducer::TensorNetworkReducer; verbose::Int = 0, vmap::Vector{Int} = collect(1:nv(g)))
    (verbose ≥ 2) && (@info "kernelizing graph: $(nv(g)) vertices, $(ne(g)) edges")
    r = 0

    vmap_0 = vmap
    reducer = deepcopy(reducer)

    while true
        res = reduce_graph(g, weights, reducer, vmap_0 = vmap_0) # res = (g_new, r_new, vmap_new)
        vmap_0 = res[4]
        vmap = vmap[res[4]]
        r += res[3]
        if g == res[1]
            (verbose ≥ 2) && (@info "kernelized graph: $(nv(g)) vertices, $(ne(g)) edges")
            return (g, weights, r, vmap, reducer)
        end
        g = res[1]
        weights = res[2]
    end
end