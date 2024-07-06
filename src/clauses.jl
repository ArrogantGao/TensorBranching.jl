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

function subcovers_naive(bs::Union{Vector{BitStr{N, T}}, AbstractVector{Vector{BitStr{N, T}}}}, vertices::Vector{Int}, g::SimpleGraph) where{N, T}
    allclauses = all_clauses_naive(bs)
    allcovers = Vector{SubCover{N, T}}()
    for (i, c) in enumerate(allclauses)
        ids = covered_items(bs, c)
        n_rm = length(rv(vertices, g, c))
        # n_rm = count_ones(c.mask)
        push!(allcovers, SubCover(ids, c, n_rm))
    end
    return allcovers
end

function subcovers(bs::Vector{BitStr{N, T}}, vertices::Vector{Int}, g::SimpleGraph) where{N, T}
    all_clauses = Set{Clause{N, T}}()
    temp_clauses = [Clause(bmask(BitStr{N, T}, 1:length(bs[i])), bs[i]) for i in 1:length(bs)]
    while !isempty(temp_clauses)
        c = pop!(temp_clauses)
        if !(c in all_clauses)
            push!(all_clauses, c)
            for b in bs
                c_new = gather2(c, Clause(bmask(BitStr{N, T}, 1:length(b)), b))
                if (c_new != c) && c_new.mask != 0
                    push!(temp_clauses, c_new)
                end
            end
        end
    end
    allcovers = [SubCover(covered_items(bs, c), c, length(rv(vertices, g, c))) for c in all_clauses]
    # allcovers = [SubCover(covered_items(bs, c), c, count_ones(c.mask)) for c in all_clauses]
    return allcovers
end

function subcovers(bss::AbstractVector{Vector{BitStr{N, T}}}, vertices::Vector{Int}, g::SimpleGraph) where{N, T}
    bs = vcat(bss...)
    all_clauses = Set{Clause{N, T}}()
    temp_clauses = [Clause(bmask(BitStr{N, T}, 1:length(bs[i])), bs[i]) for i in 1:length(bs)]
    while !isempty(temp_clauses)
        c = pop!(temp_clauses)
        if !(c in all_clauses)
            push!(all_clauses, c)
            idc = Set(covered_items(bss, c))
            for i in 1:length(bss)
                if i ∉ idc                
                    for b in bss[i]
                        c_new = gather2(c, Clause(bmask(BitStr{N, T}, 1:length(b)), b))
                        if (c_new != c) && c_new.mask != 0
                            push!(temp_clauses, c_new)
                        end
                    end
                end
            end
        end
    end
    allcovers = [SubCover(covered_items(bss, c), c, length(rv(vertices, g, c))) for c in all_clauses]
    # allcovers = [SubCover(covered_items(bss, c), c, count_ones(c.mask)) for c in all_clauses]
    return allcovers
end

function subcovers(tbl::BranchingTable{N, C}, vertices::Vector{Int}, g::SimpleGraph) where{N, C}
    return subcovers(Tbl2BitStrs(tbl), vertices, g)
end