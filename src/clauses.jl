function all_clauses_naive(n::Int, bss::AbstractVector{Vector{INT}}) where INT
    allclauses = Vector{Clause{INT}}()
    for ids in Iterators.product([0:length(bss[i]) for i in 1:length(bss)]...)
        masks = [ids...]
        cbs = [bss[i][masks[i]] for i in 1:length(bss) if masks[i] != 0]
        if length(cbs) > 0
            ccbs = clause(n::Int, cbs)
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
- `bs::Union{Vector{INT}, AbstractVector{Vector{INT}}}`: The set of bit strings.
- `vertices::Vector{Int}`: The vertices.
- `g::SimpleGraph`: The graph.
- `use_rv::Bool`: A flag indicating whether to use the number of vertices removed.

# Returns
- `allcovers::Vector{SubCover{INT}}`: The computed subcovers.

"""
function subcovers_naive(n::Int, bs::Union{Vector{INT}, AbstractVector{Vector{INT}}}, vertices::Vector{Int}, g::SimpleGraph, measurement::AbstractMeasure) where {INT}
    allclauses = all_clauses_naive(n, bs)
    allcovers = Vector{SubCover{INT}}()
    for (i, c) in enumerate(allclauses)
        ids = covered_items(bs, c)

        n_rm = size_reduced(g, vertices, c, measurement)

        push!(allcovers, SubCover(ids, c, n_rm))
    end
    return allcovers
end

"""
    subcovers(bss::AbstractVector{Vector{INT}}, vertices::Vector{Int}, g::SimpleGraph, use_rv::Bool) where {INT}

Compute the subcovers of a set of bit strings by iteratively gathering clauses.

# Arguments
- `bss::AbstractVector{Vector{INT}}`: A vector of vectors of bit strings.
- `vertices::Vector{Int}`: A vector of integers representing vertices.
- `g::SimpleGraph`: A simple graph.
- `use_rv::Bool`: A boolean indicating whether to use the number of vertices removed.

# Returns
- `allcovers::Vector{SubCover}`: A vector of `SubCover` objects representing the subcovers.

"""
function subcovers(bss::AbstractVector{Vector{INT}}, vertices::Vector{Int}, g::SimpleGraph, measurement::AbstractMeasure) where {INT}
    n = nv(g)
    bs = vcat(bss...)
    all_clauses = Set{Clause{INT}}()
    temp_clauses = [Clause(bmask(INT, 1:n), bs[i]) for i in 1:length(bs)]
    while !isempty(temp_clauses)
        c = pop!(temp_clauses)
        if !(c in all_clauses)
            push!(all_clauses, c)
            idc = Set(covered_items(bss, c))
            for i in 1:length(bss)
                if i ∉ idc                
                    for b in bss[i]
                        c_new = gather2(n, c, Clause(bmask(INT, 1:n), b))
                        if (c_new != c) && c_new.mask != 0
                            push!(temp_clauses, c_new)
                        end
                    end
                end
            end
        end
    end

    allcovers = [SubCover(covered_items(bss, c), c, size_reduced(g, vertices, c, measurement)) for c in all_clauses]

    return allcovers
end

function subcovers(tbl::BranchingTable{INT}, vertices::Vector{Int}, g::SimpleGraph, measurement::AbstractMeasure) where{INT}
    return subcovers(tbl.table, vertices, g, measurement)
end