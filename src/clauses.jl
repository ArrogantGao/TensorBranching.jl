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

"""
    subcovers_naive(bs, vertices, g, use_rv)

Compute the subcovers for a given set of bit strings by search all possible clauses.

# Arguments
- `bs::Union{Vector{BitStr{N, T}}, AbstractVector{Vector{BitStr{N, T}}}}`: The set of bit strings.
- `vertices::Vector{Int}`: The vertices.
- `g::SimpleGraph`: The graph.
- `use_rv::Bool`: A flag indicating whether to use the number of vertices removed.

# Returns
- `allcovers::Vector{SubCover{N, T}}`: The computed subcovers.

"""
function subcovers_naive(bs::Union{Vector{BitStr{N, T}}, AbstractVector{Vector{BitStr{N, T}}}}, vertices::Vector{Int}, g::SimpleGraph, use_rv::Bool) where{N, T}
    allclauses = all_clauses_naive(bs)
    allcovers = Vector{SubCover{N, T}}()
    for (i, c) in enumerate(allclauses)
        ids = covered_items(bs, c)

        if use_rv
            n_rm = length(rv(vertices, g, c))
        else
            n_rm = count_ones(c.mask)
        end

        push!(allcovers, SubCover(ids, c, n_rm))
    end
    return allcovers
end

function subcovers(bs::Vector{BitStr{N, T}}, vertices::Vector{Int}, g::SimpleGraph, use_rv::Bool) where{N, T}
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

    if use_rv
        allcovers = [SubCover(covered_items(bs, c), c, length(rv(vertices, g, c))) for c in all_clauses]
    else
        allcovers = [SubCover(covered_items(bs, c), c, count_ones(c.mask)) for c in all_clauses]
    end

    return allcovers
end

"""
    subcovers(bss::AbstractVector{Vector{BitStr{N, T}}}, vertices::Vector{Int}, g::SimpleGraph, use_rv::Bool) where {N, T}

Compute the subcovers of a set of bit strings by iteratively gathering clauses.

# Arguments
- `bss::AbstractVector{Vector{BitStr{N, T}}}`: A vector of vectors of bit strings.
- `vertices::Vector{Int}`: A vector of integers representing vertices.
- `g::SimpleGraph`: A simple graph.
- `use_rv::Bool`: A boolean indicating whether to use the number of vertices removed.

# Returns
- `allcovers::Vector{SubCover}`: A vector of `SubCover` objects representing the subcovers.

"""
function subcovers(bss::AbstractVector{Vector{BitStr{N, T}}}, vertices::Vector{Int}, g::SimpleGraph, use_rv::Bool) where{N, T}
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

    if use_rv
        allcovers = [SubCover(covered_items(bss, c), c, length(rv(vertices, g, c))) for c in all_clauses]
    else
        allcovers = [SubCover(covered_items(bss, c), c, count_ones(c.mask)) for c in all_clauses]
    end

    return allcovers
end

function subcovers(tbl::BranchingTable{N, C}, vertices::Vector{Int}, g::SimpleGraph, use_rv::Bool) where{N, C}
    return subcovers(Tbl2BitStrs(tbl), vertices, g, use_rv)
end