abstract type BranchingStrategy end
struct NaiveBranching <: BranchingStrategy end
struct SetCoverBranching <: BranchingStrategy 
    max_itr::Int
    SetCoverBranching() = new(1)
    SetCoverBranching(max_itr::Int) = new(max_itr)
end

function missolve(g::SimpleGraph; show_count = false, strategy::BranchingStrategy = NaiveBranching(), kneighbor::Int = 2, use_rv::Bool = true)
    mis = mis_solver(g, strategy, kneighbor, use_rv)
    return show_count ? (mis.mis, mis.count) : mis.mis
end

function rv(vertices::Vector{Int}, g::SimpleGraph, clause::Clause{N}) where N
    removed_vertices = Int[]
    for (k, v) in enumerate(vertices)
        if readbit(clause.mask, k) == 1
            push!(removed_vertices, v)
            if readbit(clause.val, k) == 1
                append!(removed_vertices, neighbors(g, v))
            end
        end
    end
    return unique!(removed_vertices)
end

function mis_solver(g::SimpleGraph, strategy::BranchingStrategy, kneighbor::Int, use_rv::Bool)
    dg = degree(g)
    if nv(g) == 0 || nv(g) == 1
        return CountingMIS(nv(g))
    elseif (0 ∈ dg) || (1 ∈ dg)
        v = findfirst(x -> (x==0)||(x==1), dg)
        return 1 + mis_solver(induced_subgraph(g, setdiff(1:nv(g), v ∪ neighbors(g, v)))[1], strategy, kneighbor, use_rv)
    elseif maximum(dg) ≥ 6
        v = argmax(dg)
        return max(1 + mis_solver(induced_subgraph(g, setdiff(1:nv(g), v ∪ neighbors(g, v)))[1], strategy, kneighbor, use_rv), mis_solver(induced_subgraph(g, setdiff(1:nv(g), v))[1], strategy, kneighbor, use_rv))
    else
        vertices, openvertices, dnf = optimal_branching_dnf(g, strategy, kneighbor, use_rv)
        @assert !isempty(vertices)

        mis_count = Vector{CountingMIS}(undef, length(dnf.clauses))
        
        @threads for i in 1:length(dnf.clauses)
            clause = dnf.clauses[i]
            removed_vertices = rv(vertices, g, clause)
            gi = copy(g)
            rem_vertices!(gi, removed_vertices)
            @assert !isempty(removed_vertices)
            mis_count[i] = mis_solver(gi, strategy, kneighbor, use_rv) + count_ones(clause.val)
        end
        max_mis = max(mis_count...)

        return max_mis
    end
end

function neighbor_cover(g::SimpleGraph, v::Int, k::Int)
    @assert k >= 0
    vertices = [v]
    for _ = 1:k
        vertices = union(vertices, (neighbors(g, w) for w in vertices)...)
    end
    openvertices = [v for v in vertices if !all(x->x ∈ vertices, neighbors(g, v))]
    return vertices, openvertices
end

function optimal_branching_dnf(g::SimpleGraph, strategy::BranchingStrategy, kneighbor::Int, use_rv::Bool)
    # reference: Exaxt Exponential Algorithms by Fomin and Kratsch, chapter 2.3
    degree_g = degree(g)

    # use minimum boundary instead
    vs_min = Int[]
    ovs_min = [1:nv(g)...]
    for v in 1:nv(g)
        vs, ovs = neighbor_cover(g, v, kneighbor)
        if length(ovs) < length(ovs_min)
            vs_min = vs
            ovs_min = ovs
        end
    end
    vertices = vs_min
    openvertices = ovs_min

    subg, vmap = induced_subgraph(g, vertices)
    tbl = reduced_alpha_configs(subg, Int[findfirst(==(v), vertices) for v in openvertices])
    filtered_tbl = env_filter(tbl, g, vertices, openvertices)
    return vertices, openvertices, impl_strategy(filtered_tbl, strategy, vertices, g, use_rv)
end

function impl_strategy(tbl::BranchingTable{N}, strategy::NaiveBranching, vertices, g, use_rv) where N
    return DNF([Clause(bmask(BitStr{N, Int}, 1:N), BitStr(first(x))) for x in tbl.table])
end

function impl_strategy(tbl::BranchingTable{N}, strategy::SetCoverBranching, vertices::Vector{Int}, g::SimpleGraph, use_rv::Bool) where N
    return setcover_strategy(tbl, vertices, g, max_itr = strategy.max_itr, use_rv)
end

function BitBasis.BitStr(sv::StaticBitVector)
    @assert length(sv.data) == 1 "bit string too long!"
    return BitBasis.BitStr{length(sv), Int}(sv.data[1])
end