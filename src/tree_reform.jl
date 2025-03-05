export remove_tensors, remove_tensors!, tensors_removed, unsafe_flatten, rethermalize

function remove_tensors!(code::DynamicNestedEinsum{LT}, tids::Vector{Int}) where LT
    _remove_tensors!(code, code.eins.iy, Set(tids))
    reform_tree!(code)
    return code
end

function remove_tensors(code::Union{DynamicNestedEinsum{LT}, SlicedEinsum{LT}}, tids::Vector{Int}) where LT
    ccode = deepcopy(code)
    (ccode isa SlicedEinsum) && (ccode = ccode.eins)
    remove_tensors!(ccode, tids)
    reform_tree!(ccode)
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

# reformulate tree as binary tree
function reform_tree!(code::DynamicNestedEinsum{LT}) where LT
    idx = Dict(OMEinsum._flatten(code))
    _reform_tree!(code, idx)
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

function rethermalize(code::Union{DynamicNestedEinsum{LT}, SlicedEinsum{LT}}, size_dict::Dict{LT, Int}; βs = 14:0.05:15, ntrials = 1) where LT
    return optimize_code(code, size_dict, TreeSA(initializer = :specified, βs=βs, ntrials=ntrials)).eins
end
