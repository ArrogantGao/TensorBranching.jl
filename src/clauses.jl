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

function all_clauses_naive(bss::AbstractVector{Vector{BitStr{N, T}}}) where{N, T}
    allclauses = Vector{Clause{N, T}}()
    for ids in Iterators.product([0:length(bss[i]) for i in 1:length(bss)]...)
        masks = [ids...]
        cbs = [bss[i][masks[i]] for i in 1:length(bss) if masks[i] != 0]
        if length(cbs) > 0
            ccbs = clause(cbs)
            if !(ccbs in allclauses) && (ccbs.mask != 0)
                push!(allclauses, ccbs)
            end
        end
    end
    return allclauses
end

function subcovers_naive(bs::Union{Vector{BitStr{N, T}}, AbstractVector{Vector{BitStr{N, T}}}}) where{N, T}
    allclauses = all_clauses_naive(bs)
    allcovers = Vector{SubCover{N, T}}()
    for (i, c) in enumerate(allclauses)
        ids = covered_items(bs, c)
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
                ids = covered_items(bs, c)
                j = length(ids)
                if ids ∉ allcovers[j]
                    push!(allcovers[j], SubCover(ids, c))
                end
            end
        end
    end
    return vcat(allcovers...)
end

function subcovers(bss::AbstractVector{Vector{BitStr{N, T}}}) where{N, T}
    allcovers = [Vector{SubCover{N, T}}() for i in 1:length(bss)]
    # 1-clauses
    for i in 1:length(bss)
        for b in bss[i]
            push!(allcovers[1], SubCover([i], Clause(bmask(BitStr{N, T}, 1:length(b)), b)))
        end
    end

    # 2-clauses and above
    for i in 2:length(bss) - 1
        temp_clauses = Vector{Clause{N, T}}()
        for c in allcovers[i-1]
            for j in 1:length(bss)
                if j ∉ c.ids
                    for b in bss[j]
                        c_new = gather2(c.clause, Clause(bmask(BitStr{N, T}, 1:length(b)), b))
                        if !(c_new in temp_clauses)
                            push!(temp_clauses, c_new)
                        end
                    end
                end
            end
        end
        for c in temp_clauses
            if c.mask != 0
                ids = covered_items(bss, c)
                j = length(ids)
                if c ∉ allcovers[j]
                    push!(allcovers[j], SubCover(ids, c))
                end
            end
        end
    end

    return vcat(allcovers...)
end

function subcovers(tbl::BranchingTable{N, C}) where{N, C}
    return subcovers(Tbl2BitStrs(tbl))
end