using OMEinsum: getixsv, getiyv, LeafString, flatten, _flatten, isleaf

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

function _log2_einsize(eincode::ET, size_dict::Dict{LT, Int}) where {ET, LT}
    return foldl((x, y) -> x + log2(size_dict[y]), eincode.iy, init = 0.0)
end

function get_subtree_pre(code::CT, size_dict::Dict{LT, Int}, threshold::T) where {CT, LT, T}
    for subtree in PreOrderDFS(code)
        (subtree isa LeafString) && continue
        if _log2_einsize(subtree.eins, size_dict) ≤ threshold
            return subtree
        end
    end
    # if no subtree larger than threshold, return the whole code
    return code
end

function get_subtree_post(code::CT, size_dict::Dict{LT, Int}, threshold::T) where {CT, LT, T}
    for subtree in PostOrderDFS(code)
        (subtree isa LeafString) && continue
        if _log2_einsize(subtree.eins, size_dict) ≥ threshold
            return subtree
        end
    end
    # if no subtree larger than threshold, return the whole code
    return code
end

function list_subtree(code::CT, size_dict::Dict{LT, Int}, threshold::T) where {CT, LT, T}
    subtrees = Vector{CT}()
    for subtree in PostOrderDFS(code)
        (subtree isa LeafString) && continue
        if _log2_einsize(subtree.eins, size_dict) ≥ threshold
            push!(subtrees, subtree)
        end
    end
    return subtrees
end

function most_label_subtree(code::CT, size_dict::Dict{LT, Int}, threshold::T) where {CT, LT, T}
    subtrees = list_subtree(code, size_dict, threshold)
    num_of_labels = [length(OMEinsum.uniquelabels(subtree)) for subtree in subtrees]
    return subtrees[findmax(num_of_labels)[2]]
end

# remove tensors from a network and reformulate the contraction tree

function remove_tensors!(code::DynamicNestedEinsum{LT}, tids::Vector{Int}) where LT
    _remove_tensors!(code, code.eins.iy, Set(tids))
    reform_tree!(code)
    return code
end

function remove_tensors(code::Union{DynamicNestedEinsum{LT}, SlicedEinsum{LT}}, tids::Vector{Int}) where LT
    ccode = true_eincode(deepcopy(code))
    remove_tensors!(ccode, tids)
    return ccode
end

function is_removed(code::CT, tids) where CT
    return hasfield(CT, :tensorindex) ? (code.tensorindex ∈ tids) : false
end

function _remove_tensors!(code::DynamicNestedEinsum{LT}, iy::Vector{LT}, tids::Set{Int}) where LT
    OMEinsum.isleaf(code) && return is_removed(code, tids) ? (true, LT[]) : (false, iy)
    
    dels = Int[]
    for (i, ix) in enumerate(code.eins.ixs)
        flag, new_ix = _remove_tensors!(code.args[i], ix, tids)
        if flag
            push!(dels, i)
        else
            code.eins.ixs[i] = new_ix
        end
    end
    deleteat!(code.eins.ixs, dels)
    deleteat!(code.args, dels)

    niy = isempty(code.eins.ixs) ? LT[] : union(code.eins.ixs...)
    dely = findall(x -> !(x ∈ niy), code.eins.iy)
    deleteat!(code.eins.iy, dely)
    return isempty(code.eins.ixs), code.eins.iy
end

function tensors_removed(code::Union{DynamicNestedEinsum{LT}, SlicedEinsum{LT}}, vs::Vector{LT}) where LT
    tids = Int[]
    ixd = Dict(OMEinsum._flatten(code))
    for tid in keys(ixd)
        if any(x -> x ∈ vs, ixd[tid])
            push!(tids, tid)
        end
    end
    return tids
end

function is_binary(code::DynamicNestedEinsum{LT}) where LT
    for node in PreOrderDFS(code)
        (node isa LeafString) && continue
        length(node.eins.ixs) != 2 && return false
    end
    return true
end

@inline vec_replace!(old::Vector{T}, new::Vector{T}) where T = append!(empty!(old), new)

# reformulate tree as binary tree
function reform_tree!(code::DynamicNestedEinsum{LT}) where LT
    idx = Dict(OMEinsum._flatten(code))
    isempty(idx) && return code

    # a special case, when the root of the tree is non binary
    reformed_code, _ = _reform_tree!(code, idx)

    if reformed_code.eins.ixs != code.eins.ixs
        vec_replace!(code.eins.ixs, reformed_code.eins.ixs)
        vec_replace!(code.args, reformed_code.args)
    end

    @assert is_binary(code)

    return code
end

function _reform_tree!(code::DynamicNestedEinsum{LT}, idx::Dict{Int, Vector{LT}}) where LT
    OMEinsum.isleaf(code) && return (code, idx[code.tensorindex])

    if length(code.args) == 1
        return _reform_tree!(code.args[1], idx)
    end
    
    for (i, arg) in enumerate(code.args)
        ncode, nix = _reform_tree!(arg, idx)
        code.args[i] = ncode
        code.eins.ixs[i] = nix
    end
    return code, code.eins.iy
end

function unsafe_flatten(code::DynamicNestedEinsum{LT}) where LT
    ixd = Dict(OMEinsum._flatten(code))
    DynamicEinCode([haskey(ixd, i) ? ixd[i] : LT[] for i=1:maximum(keys(ixd))], collect(OMEinsum.getiy(code.eins)))
end

function rethermalize(code::Union{DynamicNestedEinsum{LT}, SlicedEinsum{LT}}, size_dict::Dict{LT, Int}, βs::StepRangeLen, ntrials::Int, niters::Int, sc_target::Int) where LT
    return optimize_code(code, size_dict, TreeSA(initializer = :specified, βs=βs, ntrials=ntrials, niters=niters, sc_target=sc_target)).eins
end

# function reconfig_code(code::Union{DynamicNestedEinsum{LT}, SlicedEinsum{LT}}, vs::Vector{LT}, size_dict::Dict{LT, Int}, do_rethermalize::Bool, βs::StepRangeLen, ntrials::Int, niters::Int, sc_target::Int) where LT
#     tids = tensors_removed(code, vs)
#     new_code = remove_tensors(code, tensors_removed(code, vs))
#     return do_rethermalize ? rethermalize(new_code, size_dict, βs, ntrials, niters, sc_target) : new_code
# end

# this part is about reindex the tree with a vertex map

function inverse_vmap(vmap::Vector{Int})
    ivmap = zeros(Int, maximum(vmap))
    for (i, v) in enumerate(vmap)
        ivmap[v] = i
    end
    return ivmap
end

# reindex the tree with a vertex map
function reindex_tree!(code::DynamicNestedEinsum{LT}, vmap::Vector{Int}) where LT
    _reindex_tree!(code, inverse_vmap(vmap))
    return code
end
function _reindex_tree!(code::DynamicNestedEinsum{LT}, ivmap::Vector{Int}) where LT
    OMEinsum.isleaf(code) && return nothing
    
    # notice that eins.ixs[i] is actually eins.args[i].eins.iy, the same in memory
    for (i, ix) in enumerate(code.eins.ixs)
        for (j, ixi) in enumerate(ix)
            ix[j] = ivmap[ixi]
        end
        _reindex_tree!(code.args[i], ivmap)
    end
    
    return nothing
end

# Personally, I think the design of SlicedEinsum is terrible, same function, different output type
function true_eincode(code::Union{DynamicNestedEinsum{LT}, SlicedEinsum{LT}}) where LT
    return code isa SlicedEinsum ? code.eins : code
end


# TODO: map the contraction order of the original graph to the new graph by using the same elimination order, if we want to kernelize the graph after each branching step