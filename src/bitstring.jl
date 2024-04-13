struct Clause{N, T}
    mask::BitStr{N, T}
    val::BitStr{N, T}
end

# returns true if bit string a is covered by bit string b
function covered_by(a, b, mask)
    return (a & mask) == (b & mask)
end
covered_by(a, clause::Clause) = covered_by(a, clause.val, clause.mask)

function covered_items(bitstrings, clause::Clause)
    return [k for (k, b) in enumerate(bitstrings) if covered_by(b, clause)]
end

function flip_all(b::BitStr{N, T}) where{N, T}
    return flip(b, bmask(BitStr{N, T}, 1:N))
end

function maxmask(bitstrings::AbstractVector{BitStr{N, T}}) where {N, T}
    mask = bmask(BitStr{N, T}, 1:N)
    for i in 1:length(bitstrings) - 1
        mask &= bitstrings[i] ⊻ flip_all(bitstrings[i+1])
    end
    val = bitstrings[1] & mask
    return Clause(mask, val)
end