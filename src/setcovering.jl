function all_clauses_naive(bs::Vector{BitStr{N, T}}) where{N, T}
    allclauses = Vector{Clause{N, T}}()
    for ids in Iterators.product([0:1 for i in 1:length(bs)]...)
        masks = Bool.([ids...])
        cbs = bs[masks]
        if length(cbs) > 0
            ccbs = clause(cbs)
            if !(ccbs in allclauses) && (ccbs.mask != 0)
                push!(allclauses, ccbs)
            end
        end
    end
    return allclauses
end

function subcovers_naive(bs::Vector{BitStr{N, T}}) where{N, T}
    allclauses = all_clauses_naive(bs)
    allcovers = Vector{SubCover{N, T}}()
    for (i, c) in enumerate(allclauses)
        ids = Vector{Int}()
        for (j, b) in enumerate(bs)
            if covered_by(b, c)
                push!(ids, j)
            end
        end
        push!(allcovers, SubCover(ids, c))
    end
    return allcovers
end

function subcovers(bs::Vector{BitStr{N, T}}) where{N, T}
    allcovers = [Vector{SubCover{N, T}}() for i in 1:length(bs)]
    allcovers[1] = [SubCover([i], Clause(bmask(BitStr{N, T}, 1:length(bs[i])), bs[i])) for i in 1:length(bs)]
    for i in 2:length(bs) - 1
        temp_clauses = Vector{Clause{N, T}}()
        for c in allcovers[i-1]
            for b in bs
                c_new = gather2(c.clause, Clause(bmask(BitStr{N, T}, 1:length(b)), b))
                if !(c_new in temp_clauses)
                    push!(temp_clauses, c_new)
                end
            end
        end
        for c in temp_clauses
            if c.mask != 0
                ids = Vector{Int}()
                for (j, b) in enumerate(bs)
                    if covered_by(b, c)
                        push!(ids, j)
                    end
                end
                j = length(ids)
                if ids ∉ allcovers[j]
                    push!(allcovers[j], SubCover(ids, c))
                end
            end
        end
    end
    return vcat(allcovers...)
end