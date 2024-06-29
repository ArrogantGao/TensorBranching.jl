struct Clause{N, T}
    mask::BitStr{N, T}
    val::BitStr{N, T}
    function Clause(mask::BitStr{N, T}, val::BitStr{N, T}) where{N, T}
        new{N, T}(mask, val & mask)
    end
end

Clause(b::BitStr{N, T}) where {N, T} = Clause(bmask(BitStr{N, T}, 1:N), b)

Base.show(io::IO, c::Clause{N, T}) where{N, T} = print(io, "Clause{$N, $(count_ones(c.mask)), $T}: mask: $(c.mask), val: $(c.val)")
booleans(n::Int) = [Clause(bmask(BitStr{n, Int64}, i), bmask(BitStr{n, Int64}, i)) for i=1:n]
GenericTensorNetworks.:∧(x::Clause, xs::Clause...) = Clause(reduce(|, getfield.(xs, :mask); init=x.mask), reduce(|, getfield.(xs, :val); init=x.val))
GenericTensorNetworks.:¬(x::Clause) = Clause(x.mask, flip(x.val, x.mask))

struct SubCover{N, T}
    ids::Set{Int}
    clause::Clause{N, T}
    n_rm::Int # number of vertices to remove
end

SubCover(ids::Vector{Int}, clause::Clause, n_rm::Int) = SubCover(Set(ids), clause, n_rm)

Base.show(io::IO, sc::SubCover{N, T}) where {N, T} = print(io, "SubCover{$N, $T}: ids: $(sc.ids), mask: $(sc.clause.mask), val: $(sc.clause.val), n_rm: $(sc.n_rm)")
Base.:(==)(sc1::SubCover{N, T}, sc2::SubCover{N, T}) where {N, T} = (sc1.ids == sc2.ids) && (sc1.clause == sc2.clause)
function Base.in(ids::Set{Int}, subcovers::AbstractVector{SubCover{N, T}}) where {N, T}
    for sc in subcovers
        if sc.ids == ids
            return true
        end
    end
    return false
end
Base.in(ids::Vector{Int}, subcovers::AbstractVector{SubCover{N, T}}) where {N, T} = Set(ids) ∈ subcovers
function Base.in(clause::Clause, subcovers::AbstractVector{SubCover{N, T}}) where {N, T}
    for sc in subcovers
        if sc.clause == clause
            return true
        end
    end
    return false
end

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
    return [k for (k, b) in enumerate(bitstrings) if any(covered_by(b, clause))]
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

function gather2(c1::Clause{N, T}, c2::Clause{N, T}) where {N, T}
    b1 = c1.val & c1.mask
    b2 = c2.val & c2.mask
    mask = (b1 ⊻ flip_all(b2)) & c1.mask & c2.mask
    val = b1 & mask
    return Clause(mask, val)
end