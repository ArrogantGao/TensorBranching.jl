using OptimalBranching.OptimalBranchingCore: candidate_clauses, covered_items
using OptimalBranching.OptimalBranchingMIS: removed_vertices

function ob_region(g::SimpleGraph{Int}, code::DynamicNestedEinsum{Int}, slicer::ContractionTreeSlicer, selector::MaxIntersectRS, size_dict::Dict{Int, Int}, verbose::Int)

    large_tensors = list_subtree(code, size_dict, slicer.sc_target)
    large_tensors_iys = [t.eins.iy for t in large_tensors]
    unique_large_tensors_iys = unique!(vcat(large_tensors_iys...))

    best_region = Vector{Int}()
    best_loss = 0

    for i in unique_large_tensors_iys
        region_i = OptimalBranchingMIS.select_region(g, i, selector.n_max, selector.strategy)
        if selector.loss == :num_uniques
            loss = length(intersect(region_i, unique_large_tensors_iys))
        else
            error("Loss function $(selector.loss) not implemented")
        end

        (loss > best_loss) && (best_loss = loss; best_region = region_i)
    end

    (verbose ≥ 2) && (@info "best region: $best_region \n loss: $best_loss")

    return best_region, best_loss
end

# I notice that in some special cases, different candidates can have the same removed vertices
# in this case, I merge the covered items of these candidates and take the maximum fixed ones
function generate_subsets(g::SimpleGraph{Int}, tbl::BranchingTable, region::Vector{Int})

    candidates = candidate_clauses(tbl)

    dict_rvs = Dict{Vector{Int}, Vector{Int}}()
    for (i, clause) in enumerate(candidates)
        rv = sort!(removed_vertices(region, g, clause))
        if haskey(dict_rvs, rv)
            push!(dict_rvs[rv], i)
        else
            dict_rvs[rv] = [i]
        end
    end

    rvs = collect(keys(dict_rvs))
    subsets = Vector{Vector{Int}}()
    fixed_ones = zeros(Int, length(rvs))
    for (i, rv) in enumerate(rvs)
        clauses_ids = dict_rvs[rv]
        items = Vector{Int}()
        for j in clauses_ids
            append!(items, covered_items(tbl.table, candidates[j]))
            fixed_ones[i] = max(fixed_ones[i], count_ones(candidates[j].val & candidates[j].mask))
        end
        push!(subsets, unique!(items))
    end

    return subsets, rvs, fixed_ones
end

function optimal_branches(g::SimpleGraph{Int}, code::DynamicNestedEinsum{Int}, slicer::ContractionTreeSlicer, region::Vector{Int}, size_dict::Dict{Int, Int}, verbose::Int)

    p = MISProblem(g)
    tbl = branching_table(p, slicer.table_solver, region)
    subsets, rvs, fixed_ones = generate_subsets(g, tbl, region)

    (verbose ≥ 2) && (@info "generating candidates, table: $(length(tbl.table)), candidates: $(length(rvs))")

    # codes = [remove_tensors(code, tensors_removed(code, rv)) for rv in rvs]
    # losses = [slicer_loss(g, code, codes[i], slicer, size_dict) for i in 1:length(rvs)]

    losses = slicer_loss(g, code, rvs, slicer, size_dict)

    ## calculate the loss and select the best ones
    optimal_branches_ids = weighted_minimum_set_cover(slicer.setcover_solver, losses, subsets, length(tbl.table))

    (verbose ≥ 2) && (@info "length of optimal branches: $(length(optimal_branches_ids))")

    branches = Vector{SlicedBranch{Int}}()
    for i in optimal_branches_ids
        newcode = remove_tensors(code, tensors_removed(code, rvs[i]))
        newcode = rethermalize(newcode, size_dict, slicer.βs, slicer.ntrials, slicer.niters, slicer.sc_target)
        gi, vmapi = induced_subgraph(g, setdiff(1:nv(g), rvs[i]))
        !iszero(nv(gi)) && reindex_tree!(newcode, vmapi)
        push!(branches, SlicedBranch(gi, newcode, fixed_ones[i]))
    end

    return branches
end

function slicer_loss(g::SimpleGraph{Int}, code::DynamicNestedEinsum{Int}, rvs::Vector{Vector{Int}}, slicer::ContractionTreeSlicer, size_dict::Dict{Int, Int})
    if slicer.loss == :sc_score
        return sc_score(slicer.sc_target, code, rvs, size_dict)
    else
        error("Loss function $(slicer.loss) not implemented")
    end
end

function sc_score(sc_target::Int, code::DynamicNestedEinsum{Int}, rvs::Vector{Vector{Int}}, size_dict::Dict{Int, Int})
    scores = ones(Float64, length(rvs))

    large_tensors = list_subtree(code, size_dict, sc_target)
    large_tensors_iys = [Set(t.eins.iy) for t in large_tensors]

    for (i, rv) in enumerate(rvs)
        for lt_iy in large_tensors_iys
            t = max(0, length(lt_iy) - reduce((y, x) -> y += (x ∈ lt_iy), rv, init = 0) - sc_target)
            scores[i] += 2.0^t - 1.0
        end
    end

    return scores
end