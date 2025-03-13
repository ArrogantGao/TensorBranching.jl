# tools for graph, the main purpose is to use a mask to represent the vertices that are still in the graph

struct MaskedGraph{T, INT} <: AbstractGraph{T}
    g::SimpleGraph{T}
    mask::INT
    function MaskedGraph(g::SimpleGraph{T}) where T
        mask = bmask(BitBasis.longinttype(nv(g), 2), 1:nv(g))
        return new{T, typeof(mask)}(g, mask)
    end
    function MaskedGraph(g::SimpleGraph{T}, mask::INT) where {T, INT<:Integer}
        return new{T, INT}(g, mask)
    end
end
Base.show(io::IO, mg::MaskedGraph{T}) where T = print(io, "{$(nv(mg)), $(ne(mg))} masked $(T) graph")

function mask_map(mg::MaskedGraph)
    return [i for i in 1:nv(mg.g) if isone(readbit(mg.mask, i))]
end

function inverse_mask_map(mg::MaskedGraph)
    mmap = mask_map(mg)
    immap = zeros(Int, nv(mg.g))
    for (i, v) in enumerate(mmap)
        immap[v] = i
    end
    return immap
end

function Graphs.nv(mg::MaskedGraph)
    return count_ones(mg.mask)
end

function Graphs.ne(mg::MaskedGraph)
    num_edges = 0
    for e in edges(mg.g)
        num_edges += isone(readbit(mg.mask, e.src)) && isone(readbit(mg.mask, e.dst))
    end
    return num_edges
end

function Graphs.has_edge(mg::MaskedGraph, src::Int, dst::Int)
    return has_edge(mg.g, src, dst) ? isone(readbit(mg.mask, src)) && isone(readbit(mg.mask, dst)) : false
end

# the resulting simple graph's vertices are in the same order as the masked graph
function Graphs.SimpleGraph(mg::MaskedGraph)
    sg = SimpleGraph(nv(mg))
    immap = inverse_mask_map(mg)
    for e in edges(mg.g)
        has_edge(mg.g, e.src, e.dst) && add_edge!(sg, immap[e.src], immap[e.dst])
    end
    return sg
end

function Graphs.neighbors(mg::MaskedGraph, v::Int)
    return [n for n in neighbors(mg.g, v) if isone(readbit(mg.mask, n))]
end

