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

function rethermalize(code::Union{DynamicNestedEinsum{LT}, SlicedEinsum{LT}}, size_dict::Dict{LT, Int}, βs::IT, ntrials::Int, niters::Int, sc_target::Int) where {LT, IT}
    return optimize_code(code, size_dict, TreeSA(initializer = :specified, βs=βs, ntrials=ntrials, niters=niters, sc_target=sc_target)).eins
end


# this part is about reindex the tree with a vertex map
function inverse_vmap(vmap::Vector{Int})
    ivmap = zeros(Int, maximum(vmap))
    for (i, v) in enumerate(vmap)
        ivmap[v] = i
    end
    return ivmap
end
function inverse_vmap_dict(vmap::Vector{Int})
    ivmap = Dict{Int, Int}()
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

function mis_complexity(code)
    return contraction_complexity(code, uniformsize(code, 2))
end

function auto_slicing(code, sc_target::Int)
    sc_init = Int(mis_complexity(code).sc)
    nslices = sc_init - sc_target
    
    local new_code
    while true
        new_code = optimize_code(code, uniformsize(code, 2), TreeSA(sc_target = sc_target, nslices = nslices))
        cc = mis_complexity(new_code)
        cc.sc <= sc_target && break
        nslices += max(1, Int(cc.sc - sc_target))
    end

    return new_code
end

function random_ksg(m::Int, n::Int, rho::Float64, seed::Int)
    Random.seed!(seed)
    ksg = SimpleGraph(GenericTensorNetworks.random_diagonal_coupled_graph(m, n, rho))
    return ksg
end

function contraction_peak_memory(code::NestedEinsum, size_dict)
    ixs = getixsv(code)
    tscs = Float64[]
    initial_sc = sum(prod(Float64(size_dict[i]) for i in ix) for ix in ixs)
    push!(tscs, initial_sc)
    _tsc!(tscs, code, size_dict, ixs)
    return maximum(log2.(tscs))
end

function _tsc!(tscs, code, size_dict, ixs)
    isleaf(code) && return zero(Float64)

    freed_size = zero(Float64)
    for subcode in code.args
        freed_size += _tsc!(tscs, subcode, size_dict, ixs)
    end

    future_freed_size = isempty(code.eins.ixs) ? 0.0 : sum(isempty(ix) ? 1.0 : prod(Float64(size_dict[i]) for i in ix) for ix in code.eins.ixs)
    allocated_size = isempty(code.eins.iy) ? 1.0 : prod(Float64(size_dict[i]) for i in code.eins.iy)
    new_size = tscs[end] + allocated_size - freed_size
    push!(tscs, new_size)
    return future_freed_size
end

function contraction_all_memory(code::NestedEinsum, size_dict)
    return log2(_ssc!(code, size_dict))
end

function _ssc!(code, size_dict)
    isleaf(code) && return zero(Float64)
    t = isempty(code.eins.iy) ? 1.0 : prod(Float64(size_dict[i]) for i in code.eins.iy)
    return t + sum(_ssc!(subcode, size_dict) for subcode in code.args)
end

function show_status(scs, sc_target, num_unfinished, num_finished)
    @info "current num of unfinished slices: $num_unfinished, finished slices: $num_finished"
    counts = zeros(Int, Int(maximum(scs) - minimum(scs) + 1))
    for sc in scs
        counts[Int(sc - minimum(scs) + 1)] += 1
    end
    println(barplot(Int(minimum(scs)):Int(maximum(scs)), counts, xlabel = "num of slices", ylabel = "sc, target = $(sc_target)"))
end