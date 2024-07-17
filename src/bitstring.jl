"""
    Clause{INT <: Integer}

A clause is a pair of bit strings, `mask` and `val`, where `mask` is a bit string that indicates the bits that are relevant to the clause, and `val` is a bit string that indicates the bits that must be satisfied. The clause is satisfied if and only if `val` is covered by the bit string and `mask`.
- `INT`: The number of integers as the storage.

### Examples
```jldoctest
julia> Clause(bit"1110", bit"1001")
Clause: mask: 1110 ₍₂₎, val: 1000 ₍₂₎
```

If some bit in mask is set as 0, then the corresponding bit in val must be 0.
"""
struct Clause{INT <: Integer}
    mask::INT
    val::INT
    function Clause(mask::INT, val::INT) where INT <: Integer
        new{INT}(mask, val & mask)
    end
end

Base.show(io::IO, c::Clause{INT}) where INT = print(io, "Clause{$INT}: mask: $(c.mask), val: $(c.val)")
function booleans(n::Int)
    C = (n + 63) ÷ 64
    INT = LongLongUInt{C}
    return [Clause(bmask(INT, i), bmask(INT, i)) for i=1:n]
end
GenericTensorNetworks.:∧(x::Clause, xs::Clause...) = Clause(reduce(|, getfield.(xs, :mask); init=x.mask), reduce(|, getfield.(xs, :val); init=x.val))
GenericTensorNetworks.:¬(x::Clause) = Clause(x.mask, flip(x.val, x.mask))

"""
    SubCover{INT <: Integer}

A subcover is a pair of a set of integers `ids`, a clause `clause` and a integer `n_rm`. The `ids` for the truth covered by the clause, and `n_rm` is the number of vertices to remove.
- `INT`: The number of integers as the storage.

### Examples
```jldoctest
julia> SubCover([1, 2], Clause(bit"1110", bit"1001"), 3)
SubCover{4, Int64}: ids: Set([2, 1]), mask: 1110 ₍₂₎, val: 1000 ₍₂₎, n_rm: 3
```

"""
struct SubCover{INT <: Integer}
    ids::Set{Int}
    clause::Clause{INT}
    n_rm::Int # number of vertices to remove
end

SubCover(ids::Vector{Int}, clause::Clause, n_rm::Int) = SubCover(Set(ids), clause, n_rm)

Base.show(io::IO, sc::SubCover{INT}) where INT = print(io, "SubCover{$INT}: ids: $(sc.ids), mask: $(sc.clause.mask), val: $(sc.clause.val), n_rm: $(sc.n_rm)")
Base.:(==)(sc1::SubCover{INT}, sc2::SubCover{INT}) where {INT} = (sc1.ids == sc2.ids) && (sc1.clause == sc2.clause)
function Base.in(ids::Set{Int}, subcovers::AbstractVector{SubCover{INT}}) where {INT}
    for sc in subcovers
        if sc.ids == ids
            return true
        end
    end
    return false
end
Base.in(ids::Vector{Int}, subcovers::AbstractVector{SubCover{INT}}) where {INT} = Set(ids) ∈ subcovers
function Base.in(clause::Clause, subcovers::AbstractVector{SubCover{INT}}) where INT <: Integer
    for sc in subcovers
        if sc.clause == clause
            return true
        end
    end
    return false
end

function BitBasis.bdistance(c1::Clause{INT}, c2::Clause{INT}) where INT <: Integer
    b1 = c1.val & c1.mask & c2.mask
    b2 = c2.val & c1.mask & c2.mask
    return bdistance(b1, b2)
end

function BitBasis.bdistance(c::Clause{INT}, b::INT) where INT <: Integer
    b1 = b & c.mask
    c1 = c.val & c.mask
    return bdistance(b1, c1)
end

"""
    covered_by(a::LongLongUInt, b::LongLongUInt, mask::LongLongUInt)

Check if `a` is covered by `b` with `mask`. The function returns `true` if and only if `a` and `b` are the same when masked by `mask`.
"""
function covered_by(a, b, mask)
    return (a & mask) == (b & mask)
end

"""
    covered_by(a::LongLongUInt, clause::Clause)

Check if `a` is covered by the clause. The function returns `true` if and only if `a` and `clause.val` are the same when masked by `clause.mask`.
"""
covered_by(a, clause::Clause) = covered_by(a, clause.val, clause.mask)

function covered_by(as::AbstractArray, clause::Clause)
    return [covered_by(a, clause) for a in as]
end

"""
    covered_items(bitstrings, clause::Clause)

Return the indices of the bit strings that are covered by the clause.
"""
function covered_items(bitstrings, clause::Clause)
    return [k for (k, b) in enumerate(bitstrings) if any(covered_by(b, clause))]
end

function flip_all(n::Int, b::INT) where INT <: Integer
    return flip(b, bmask(INT, 1:n))
end

"""
    clause(n::Int, bitstrings::AbstractVector{INT})

Return a clause that covers all the bit strings.
"""
function clause(n::Int, bitstrings::AbstractVector{INT}) where INT
    mask = bmask(INT, 1:n)
    for i in 1:length(bitstrings) - 1
        mask &= bitstrings[i] ⊻ flip_all(n, bitstrings[i+1])
    end
    val = bitstrings[1] & mask
    return Clause(mask, val)
end

function clauses(n::Int, clustered_bs)
    return [clause(n, bitstrings) for bitstrings in clustered_bs]
end

function gather2(n::Int, c1::Clause{INT}, c2::Clause{INT}) where INT
    b1 = c1.val & c1.mask
    b2 = c2.val & c2.mask
    mask = (b1 ⊻ flip_all(n, b2)) & c1.mask & c2.mask
    val = b1 & mask
    return Clause(mask, val)
end

function BitBasis.LongLongUInt(sv::StaticBitVector)
    @assert length(sv.data) == 1 "bit string too long!"
    return LongLongUInt(sv.data)
end