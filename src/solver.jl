function missolve(g::SimpleGraph)
    if nv(g) == 0
        return 0
    elseif nv(g) == 1
        return 1
    end
    vertices, openvertices, dnf = optimal_branching_dnf(g)
    @assert !isempty(vertices)
    return maximum(dnf.clauses) do clause
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
        missolve(gi) + count_ones(clause.val)
    end
end

function optimal_branching_dnf(g::SimpleGraph)
    # vertices up to second nearest neighbors
    v = rand(vertices(g))
    nn = neighbors(g, v)
    nnn = union(nn, [v], [neighbors(g, u) for u in nn]...)
    openvertices = [v for v in nnn if !all(u -> u ∈ nnn, neighbors(g, v))]
    subg, vmap = induced_subgraph(g, nnn)
    tbl = reduced_alpha_configs(subg, Int[findfirst(==(v), nnn) for v in openvertices])
    return nnn, openvertices, naive_strategy(tbl)
end

function naive_strategy(tbl::BranchingTable{N}) where N
    return DNF([Clause(bmask(BitStr{N, Int}, 1:N), BitStr(first(x))) for x in tbl.table])
end
function BitBasis.BitStr(sv::StaticBitVector)
    @assert length(sv.data) == 1 "bit string too long!"
    return BitBasis.BitStr{length(sv), Int}(sv.data[1])
end