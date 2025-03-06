export get_subtree_pre, get_subtree_post, list_subtree, most_label_subtree

function _log2_einsize(eincode::ET, size_dict::Dict{LT, Int}) where {ET, LT}
    return foldl((x, y) -> x + log2(size_dict[y]), eincode.iy, init = 0.0)
end

function get_subtree_pre(code::CT, size_dict::Dict{LT, Int}, threshold::T) where {CT, LT, T}
    for subtree in PreOrderDFS(code)
        (subtree isa OMEinsum.LeafString) && continue
        if _log2_einsize(subtree.eins, size_dict) ≤ threshold
            return subtree
        end
    end
    # if no subtree larger than threshold, return the whole code
    return code
end

function get_subtree_post(code::CT, size_dict::Dict{LT, Int}, threshold::T) where {CT, LT, T}
    for subtree in PostOrderDFS(code)
        (subtree isa OMEinsum.LeafString) && continue
        if _log2_einsize(subtree.eins, size_dict) ≥ threshold
            return subtree
        end
    end
    # if no subtree larger than threshold, return the whole code
    return code
end

function list_subtree(code::CT, size_dict::Dict{LT, Int}, threshold::T) where {CT, LT, T}
    subtrees = Vector{CT}()
    _list_subtree!(subtrees, code, size_dict, threshold)
    return subtrees
end

function _list_subtree!(subtrees::Vector{CT}, code::CT, size_dict::Dict{LT, Int}, threshold::T) where {CT, LT, T}
    OMEinsum.isleaf(code) && return nothing
    if _log2_einsize(code.eins, size_dict) ≤ threshold
        push!(subtrees, code)
        return nothing
    end
    for child in code.args
        _list_subtree!(subtrees, child, size_dict, threshold)
    end
    return nothing
end

function most_label_subtree(code::CT, size_dict::Dict{LT, Int}, threshold::T) where {CT, LT, T}
    subtrees = list_subtree(code, size_dict, threshold)
    num_of_labels = [length(OMEinsum.uniquelabels(subtree)) for subtree in subtrees]
    return subtrees[findmax(num_of_labels)[2]]
end