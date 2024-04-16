function setcover(bitstrings::AbstractVector{BitStr{N, T}}, clauses = [Clause{N, T}[]]) where {N, T}
    # generate longest clauses covering 1 bitstring
    for (k, b) in enumerate(bitstrings)
        clause = Clause{N, T}(bmask(BitStr{N, T}, 1:length(b)), b)
        clauses[1][clause] = [k]
    end
    # generate clauses covering 2 bitstrings
    for kcover in 2:length(bitstrings)
        for c in clauses[kcover-1]
            # try to combine c with a bitstring
            for b in bitstrings
                if !covered_by(c.val, b, c.mask)
                    push!(clauses[kcover], Clause{N, T}(c.mask | bmask(BitStr{N, T}, 1:length(b)), c.val | b))
                end
            end
        end
    end
end

function all_clauses_naive(bs::Vector{BitStr{N, T}}) where{N, T}
    allclauses = Vector{Clause{N, T}}()
    for ids in Iterators.product([0:1 for i in 1:length(bs)]...)
        masks = Bool.([ids...])
        cbs = bs[masks]
        if length(cbs) > 0
            ccbs = clause(cbs)
            if !(ccbs in allclauses)
                push!(allclauses, ccbs)
            end
        end
    end
    return allclauses
end