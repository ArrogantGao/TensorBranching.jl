struct Clause{N, T}
    mask::BitStr{N, T}
    val::BitStr{N, T}
end

Clause(b::BitStr{N, T}) where {N, T} = Clause(bmask(BitStr{N, T}, 1:N), b)

Base.show(io::IO, c::Clause{N, T}) where{N, T} = print(io, "Clause{$N, $(count_ones(c.mask)), $T}: mask: $(c.mask), val: $(c.val)")
booleans(n::Int) = [Clause(bmask(BitStr{n, Int64}, i), bmask(BitStr{n, Int64}, i)) for i=1:n]
GenericTensorNetworks.:∧(x::Clause, xs::Clause...) = Clause(reduce(|, getfield.(xs, :mask); init=x.mask), reduce(|, getfield.(xs, :val); init=x.val))
GenericTensorNetworks.:¬(x::Clause) = Clause(x.mask, flip(x.val, x.mask))

function BitBasis.bdistance(c1::Clause{N, T}, c2::Clause{N, T}) where{N, T}
    b1 = c1.val & c1.mask & c2.mask
    b2 = c2.val & c1.mask & c2.mask
    return bdistance(b1, b2)
end

function BitBasis.bdistance(c::Clause{N, T}, b::BitStr{N, T}) where{N, T}
    b1 = b & c.mask
    c1 = c.val & c.mask
    return bdistance(b1, c1)
end

# returns true if bit string a is covered by bit string b
function covered_by(a, b, mask)
    return (a & mask) == (b & mask)
end
covered_by(a, clause::Clause) = covered_by(a, clause.val, clause.mask)

function covered_by(as::AbstractArray, clause::Clause)
    return [covered_by(a, clause) for a in as]
end

function covered_items(bitstrings, clause::Clause)
    return [k for (k, b) in enumerate(bitstrings) if covered_by(b, clause)]
end

function flip_all(b::BitStr{N, T}) where{N, T}
    return flip(b, bmask(BitStr{N, T}, 1:N))
end

function clause(bitstrings::AbstractVector{BitStr{N, T}}) where {N, T}
    mask = bmask(BitStr{N, T}, 1:N)
    for i in 1:length(bitstrings) - 1
        mask &= bitstrings[i] ⊻ flip_all(bitstrings[i+1])
    end
    val = bitstrings[1] & mask
    return Clause(mask, val)
end

function clauses(clustered_bs)
    return [clause(bitstrings) for bitstrings in clustered_bs]
end