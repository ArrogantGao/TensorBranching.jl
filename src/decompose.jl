# transform optimized eincode to elimination order
function eincode2order(code::NestedEinsum{L}) where {L}
    elimination_order = Vector{L}()
    OMEinsum.isleaf(code) && return elimination_order
    for node in PostOrderDFS(code)
        (node isa LeafString) && continue
        for id in setdiff(vcat(getixsv(node.eins)...), getiyv(node.eins))
            push!(elimination_order, id)
        end
    end
    return reverse!(elimination_order)
end

function eincode2graph(code::Union{NestedEinsum, EinCode})
    fcode = code isa NestedEinsum ? flatten(code) : code
    indices = uniquelabels(fcode)
    (indices isa Vector{Int}) && sort!(indices)
    g = SimpleGraph(length(indices))
    id_dict = Dict(id => i for (i, id) in enumerate(indices))
    for xs in [getixsv(fcode); getiyv(fcode)]
        for i in 1:length(xs)-1
            for j in i+1:length(xs)
                add_edge!(g, id_dict[xs[i]], id_dict[xs[j]])
            end
        end
    end
    return g, id_dict
end

function decompose(code::NestedEinsum{L}) where {L}
    g, id_dict = eincode2graph(code)
    labels = collect(keys(id_dict))[sortperm(collect(values(id_dict)))]
    return decomposition_tree(g, eincode2order(code), labels = labels)
end

function max_bag(tree::DecompositionTreeNode)
    max_bag = tree.bag
    max_size = length(max_bag)
    for node in PostOrderDFS(tree)
        if length(node.bag) > max_size
            max_bag = node.bag
            max_size = length(node.bag)
        end
    end
    return max_bag
end

# this function maps an elimination order on a old graph to a new graph with some vertices removed or reordered
function update_order(eo_old::Vector{Int}, vmap::Vector{Int})
    ivmap = inverse_vmap_dict(vmap)
    eo_new = Vector{Int}()
    for v in eo_old
        haskey(ivmap, v) && push!(eo_new, ivmap[v])
    end
    return eo_new
end
function update_tree(g_new::SimpleGraph{Int}, eo_old::Vector{Int}, vmap::Vector{Int})
    eo_new = update_order(eo_old, vmap)
    return decomposition_tree(g_new, eo_new)
end

# reconstruct the contraction order from the grouped elimination order
# if set use_tree to true, the decomposition tree will be constructed to get a better elimination order
function order2eincode(g::SimpleGraph{Int}, eo::Vector{Int}; use_tree::Bool = true)
    rcode = rawcode(IndependentSet(g))
    ixs = getixsv(rcode)
    incidence_list = IncidenceList(Dict([i=>ix for (i, ix) in enumerate(ixs)]))
    trees = Vector{Union{ContractionTree, Int}}()
    for sub_vs in connected_components(g)
        if length(sub_vs) == 1
            # a special corner case for the mis problem: the connected component is a single vertex has no edges
            push!(trees, incidence_list.e2v[sub_vs[1]][1])
        else
            if use_tree
                sub_g, sub_vmap = induced_subgraph(g, sub_vs)
                sub_ivmap = inverse_vmap_dict(sub_vmap)
                sub_eo = [sub_ivmap[i] for i in eo if i in sub_vs]
                tree = decomposition_tree(sub_g, sub_eo)
                grouped_eo = EliminationOrder(tree).order

                map!(x -> map!(y -> sub_vmap[y], x, x), grouped_eo, grouped_eo)
            else
                grouped_eo = [[i] for i in eo if i in sub_vs]
            end

            tree = eo2ct(grouped_eo, incidence_list, [1.0 for _ in 1:length(eo)])
            push!(trees, tree)
        end
    end
    tree = reduce((x,y) -> ContractionTree(x, y), trees)
    code = parse_eincode(incidence_list, tree, vertices = collect(1:length(ixs))) # this code is OMEinsumContractionOrders.NestedEinsum, not OMEinsum.NestedEinsum

    return decorate(code)
end

function update_code(g_new::SimpleGraph{Int}, code_old::NestedEinsum, vmap::Vector{Int})
    eo_old = eincode2order(code_old)
    eo_new = update_order(eo_old, vmap)
    return order2eincode(g_new, eo_new)
end

function ein2contraction_tree(code::NestedEinsum)
    @assert is_binary(code)
    return _ein2contraction_tree(code)
end

function _ein2contraction_tree(code)
    return isleaf(code) ? code.tensorindex : ContractionTree(_ein2contraction_tree(code.args[1]), _ein2contraction_tree(code.args[2]))
end