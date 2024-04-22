# clustering using the bit distance and hierarchical clustering
# not good, hamming distance is not transitive
function bithclust(bitstrings::Vector{BitStr{N, T}}) where {N, T}
    dm = [bdistance(b1, b2) for b1 in bitstrings, b2 in bitstrings]
    return hclust(dm, linkage = :average)
end

function clustering(bitstrings::Vector{BitStr{N, T}}, k::Int) where {N, T}
    tree = cutree(bithclust(bitstrings), k = k)
    clustered_bitstrings = Vector{Vector{BitStr{N, T}}}()
    for i in 1:k
        mask = tree .== i
        push!(clustered_bitstrings, bitstrings[mask])
    end
    return clustered_bitstrings
end

# each time gather the two clauses with the smallest bit distance
function gather2(c1::Clause{N, T}, c2::Clause{N, T}) where {N, T}
    b1 = c1.val & c1.mask
    b2 = c2.val & c2.mask
    mask = (b1 ⊻ flip_all(b2)) & c1.mask & c2.mask
    val = b1 & mask
    return Clause(mask, val)
end

function gather(clauses::AbstractVector{Clause{N, T}}) where {N, T}
    nc = length(clauses)
    @assert nc > 1
    dm = zeros(nc, nc)
    fill!(dm, Inf)
    for i in 1:nc - 1, j in i + 1:nc
        dm[i, j] = bdistance(clauses[i], clauses[j])
    end
    min_d = minimum(dm)
    ids = findall(dm .== min_d)

    # find the smallest complexity
    for id in ids
        gc = gather2(clauses[id[1]], clauses[id[2]])
        n_clauses = [clauses[i] for i in 1:nc if i ∉ Tuple(id)] ∪ [gc]
        γ_n = complexity(sbranches(n_clauses))
        γ = complexity(sbranches(clauses))
        if γ_n < γ
            return n_clauses
        end
    end

    #if failed, try a random id
    id = rand(ids)
    gc = gather2(clauses[id[1]], clauses[id[2]])
    return [clauses[i] for i in 1:nc if i ∉ Tuple(id)] ∪ [gc]
end