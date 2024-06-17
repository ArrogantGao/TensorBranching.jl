abstract type BranchingStrategy end
struct NaiveBranching <: BranchingStrategy end
struct SetCoverBranching <: BranchingStrategy 
    max_itr::Int
    SetCoverBranching() = new(1)
    SetCoverBranching(max_itr::Int) = new(max_itr)
end

function missolve(g::SimpleGraph; show_count::Bool = false, strategy::BranchingStrategy = NaiveBranching())
    if nv(g) == 0 || nv(g) == 1
        return show_count ? (nv(g), 1) : nv(g)
    end
    vertices, openvertices, dnf = optimal_branching_dnf(g, strategy = strategy)
    @assert !isempty(vertices)
    branch_count = 0
    max_mis = 0
    for clause in dnf.clauses
        removed_vertices = Int[]
        for (k, v) in enumerate(vertices)
            if readbit(clause.mask, k) == 1
                push!(removed_vertices, v)
                if readbit(clause.val, k) == 1
                    append!(removed_vertices, neighbors(g, v))
                end
            end
        end
        gi = copy(g)
        rem_vertices!(gi, removed_vertices)
        @assert !isempty(removed_vertices)
        if show_count
            mis, count = missolve(gi, show_count = true)
            branch_count += count
        else
            mis = missolve(gi)
        end
        max_mis = max(mis + count_ones(clause.val), max_mis)
    end
    return show_count ? (max_mis, branch_count) : max_mis
end

function optimal_branching_dnf(g::SimpleGraph; strategy::BranchingStrategy = NaiveBranching())
    # vertices up to second nearest neighbors
    # v = rand(vertices(g))
    # choose vertex with maximum degree
    v = argmax(degree(g))
    nn = neighbors(g, v)
    nnn = union(nn, [v], [neighbors(g, u) for u in nn]...)
    openvertices = [v for v in nnn if !all(u -> u ∈ nnn, neighbors(g, v))]
    subg, vmap = induced_subgraph(g, nnn)
    tbl = reduced_alpha_configs(subg, Int[findfirst(==(v), nnn) for v in openvertices])
    if strategy isa NaiveBranching
        return nnn, openvertices, naive_strategy(tbl)
    elseif strategy isa SetCoverBranching
        return nnn, openvertices, setcover_strategy(tbl, max_itr = strategy.max_itr)
    end
end

function naive_strategy(tbl::BranchingTable{N}) where N
    return DNF([Clause(bmask(BitStr{N, Int}, 1:N), BitStr(first(x))) for x in tbl.table])
end

function BitBasis.BitStr(sv::StaticBitVector)
    @assert length(sv.data) == 1 "bit string too long!"
    return BitBasis.BitStr{length(sv), Int}(sv.data[1])
end