using OptimalBranching.OptimalBranchingCore: candidate_clauses, covered_items
using OptimalBranching.OptimalBranchingMIS: removed_vertices, select_region, clause_size

function ob_region(g::SimpleGraph{Int}, code::DynamicNestedEinsum{Int}, slicer::ContractionTreeSlicer, selector::ScoreRS, size_dict::Dict{Int, Int}, verbose::Int)

    large_tensors = list_subtree(code, size_dict, slicer.sc_target)
    large_tensors_iys = [t.eins.iy for t in large_tensors]
    unique_large_tensors_iys = unique!(vcat(large_tensors_iys...))

    regions = [Vector{Int}() for _ in unique_large_tensors_iys]
    Threads.@threads for i in 1:length(unique_large_tensors_iys)
        iy = unique_large_tensors_iys[i]
        region_i = select_region(g, iy, selector.n_max, selector.strategy)
        regions[i] = region_i
    end

    if selector.loss == :sc_score
        losses = sc_score(slicer.sc_target, code, regions, size_dict)
        best_loss = minimum(losses)
        best_region = regions[argmin(losses)]
    elseif selector.loss == :bag_score
        losses = bag_score(slicer.sc_target, code, regions, size_dict)
        best_loss = minimum(losses)
        best_region = regions[argmin(losses)]
    elseif selector.loss == :num_uniques
        losses = zeros(Float64, length(unique_large_tensors_iys))
        Threads.@threads for i in 1:length(unique_large_tensors_iys)
            iy = unique_large_tensors_iys[i]
            losses[i] = length(intersect(regions[i], iy))
        end
        best_loss = maximum(losses)
        best_region = regions[argmax(losses)]
    else
        error("Loss function $(selector.loss) not implemented")
    end

    (verbose ≥ 2) && (@info "best region: $best_region \n loss: $best_loss")

    return best_region, best_loss
end


# I notice that in some special cases, different candidates can have the same removed vertices
# in this case, I merge the covered items of these candidates and take the maximum fixed ones
function generate_subsets(g::SimpleGraph{Int}, weights::VT, tbl::BranchingTable, region::Vector{Int}) where VT

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
    residuals = zeros(eltype(weights), length(rvs))
    for (i, rv) in enumerate(rvs)
        clauses_ids = dict_rvs[rv]
        items = Vector{Int}()
        for j in clauses_ids
            append!(items, covered_items(tbl.table, candidates[j]))
            residuals[i] = max(residuals[i], clause_size(weights, candidates[j].val & candidates[j].mask, region))
        end
        push!(subsets, unique!(items))
    end

    return subsets, rvs, residuals
end

function optimal_branches(p::MISProblem{INT, VT}, code::DynamicNestedEinsum{Int}, r::RT, slicer::ContractionTreeSlicer, reducer::AbstractReducer, region::Vector{Int}, size_dict::Dict{Int, Int}, verbose::Int) where {INT, VT, RT}

    cc = contraction_complexity(code, size_dict)
    (verbose ≥ 2) && (@info "solving g: $(nv(p.g)), $(ne(p.g)), code complexity: tc = $(cc.tc), sc = $(cc.sc)")

    tbl = branching_table(p, slicer.table_solver, region)
    (verbose ≥ 2) && (@info "table: $(length(tbl.table))")

    # special case: find reduction
    reduction = OptimalBranchingCore.intersect_clauses(tbl, :dfs)
    if !isempty(reduction)
        c = reduction[1]
        new_branch = generate_branch(p, code, sort!(removed_vertices(region, p.g, c)), clause_size(p.weights, c.val & c.mask, region), slicer, reducer, size_dict)
        return [add_r(new_branch, r)]
    end

    subsets, rvs, residuals = generate_subsets(p.g, p.weights, tbl, region)
    (verbose ≥ 2) && (@info "candidates: $(length(rvs))")

    losses = slicer_loss(p.g, code, rvs, slicer.brancher, slicer.sc_target, size_dict)

    ## calculate the loss and select the best ones
    optimal_branches_ids = set_cover(slicer.brancher, losses, subsets, length(tbl.table))

    (verbose ≥ 2) && (@info "length of optimal branches: $(length(optimal_branches_ids))")

    brs = Vector{SlicedBranch}() # brs for branches
    for i in optimal_branches_ids
        new_branch = generate_branch(p, code, rvs[i], residuals[i], slicer, reducer, size_dict)

        (verbose ≥ 2) && (@info "branching id = $i, g: $(nv(new_branch.p.g)), $(ne(new_branch.p.g)), rv = $(rvs[i])")
        (verbose ≥ 2) && (cc_ik = complexity(new_branch); @info "rethermalized code complexity: tc = $(cc_ik.tc), sc = $(cc_ik.sc)")

        push!(brs, add_r(new_branch, r))
    end

    return brs
end

function slicer_loss(g::SimpleGraph{Int}, code::DynamicNestedEinsum{Int}, rvs::Vector{Vector{Int}}, brancher::GreedyBrancher, sc_target::Int, size_dict::Dict{Int, Int})
    if brancher.loss == :sc_score
        return sc_score(sc_target, code, rvs, size_dict)
    elseif brancher.loss == :bag_score
        return bag_score(sc_target, code, rvs, size_dict)
    else
        error("Loss function $(brancher.loss) not implemented")
    end
end

function slicer_loss(g::SimpleGraph{Int}, code::DynamicNestedEinsum{Int}, rvs::Vector{Vector{Int}}, brancher::FixedPointBrancher, sc_target::Int, size_dict::Dict{Int, Int})
    if brancher.measure == :sc_measure
        return sc_measure(sc_target, code, rvs, size_dict)
    else
        error("Loss function $(brancher.loss) not implemented")
    end
end

# sc_score is the sum of the scores of the removed vertices, each of the tensors with size larger than sc_target contributes 2^(t - sc_target) - 1 to the score, where t is size of the tensor
function sc_score(sc_target::Int, code::DynamicNestedEinsum{Int}, rvs::Vector{Vector{Int}}, size_dict::Dict{Int, Int})
    scores = ones(Float64, length(rvs))

    large_tensors = list_subtree(code, size_dict, sc_target)
    large_tensors_iys = [Set(t.eins.iy) for t in large_tensors]

    Threads.@threads for i in 1:length(rvs)
        rv = rvs[i]
        for lt_iy in large_tensors_iys
            t = max(0, length(lt_iy) - reduce((y, x) -> y += (x ∈ lt_iy), rv, init = 0) - sc_target)
            scores[i] += 2.0^t - 1.0
        end
    end

    return scores
end

# tree bag score 
function bag_score(sc_target::Int, code::DynamicNestedEinsum{Int}, rvs::Vector{Vector{Int}}, size_dict::Dict{Int, Int})
    tree = decompose(code)
    scores = ones(Float64, length(rvs))
    
    intersects = Vector{Set{Int}}()
    for node in PostOrderDFS(tree)
        if !isempty(node.children)
            for child in node.children
                push!(intersects, intersect(node.bag, child.bag))
            end
        end
    end

    Threads.@threads for i in 1:length(rvs)
        rv = rvs[i]
        for lt_iy in intersects
            t = max(0, length(lt_iy) - reduce((y, x) -> y += (x ∈ lt_iy), rv, init = 0) - sc_target)
            scores[i] += 2.0^t - 1.0
        end
    end

    return scores
end

# sc_measure is similar to the D3 measure, where each tensor is counted by (t - sc_target)
function sc_measure(sc_target::Int, code::DynamicNestedEinsum{Int}, rvs::Vector{Vector{Int}}, size_dict::Dict{Int, Int})
    delta_rho = zeros(Float64, length(rvs))

    large_tensors = list_subtree(code, size_dict, sc_target)
    large_tensors_iys = [Set(t.eins.iy) for t in large_tensors]

    Threads.@threads for i in 1:length(rvs)
        rv = rvs[i]
        for lt_iy in large_tensors_iys
            delta_rho[i] += reduce((y, x) -> y += (x ∈ lt_iy), rv, init = 0)
        end
    end

    return delta_rho
end

function set_cover(solver::GreedyBrancher, losses::Vector{Float64}, subsets::Vector{Vector{Int}}, n_clauses::Int)
    return weighted_minimum_set_cover(solver.setcover_solver, losses, subsets, n_clauses)
end

function set_cover(solver::FixedPointBrancher, losses::Vector{Float64}, subsets::Vector{Vector{Int}}, n_clauses::Int)
    return fixed_point_set_cover(solver.setcover_solver, losses, subsets, n_clauses)
end

function fixed_point_set_cover(solver::AbstractSetCoverSolver, losses::Vector{Float64}, subsets::Vector{Vector{Int}}, n_clauses::Int)
    cx_old = cx = solver.γ0
    local picked_scs
    for i = 1:solver.max_itr
        weights = 1 ./ cx_old .^ losses
        picked_scs = weighted_minimum_set_cover(solver, weights, subsets, n_clauses)
        cx = OptimalBranching.OptimalBranchingCore.complexity_bv(losses[picked_scs])
        cx ≈ cx_old && break  # convergence
        cx_old = cx
    end
    return picked_scs
end

function generate_branch(p::MISProblem{INT, VT}, code::DynamicNestedEinsum{Int}, removed_vertices::Vector{Int}, r0::RT, slicer::ContractionTreeSlicer, reducer::AbstractReducer, size_dict::Dict{Int, Int}) where {INT, VT, RT}
    g = p.g
    weights = p.weights
    g_i, vmap_i = induced_subgraph(g, setdiff(1:nv(g), removed_vertices))
    weights_i = weights[vmap_i]

    res = kernelize(g_i, weights_i, reducer, vmap = vmap_i)
    g_ik = res.g
    weights_ik = res.weights
    r_ik = res.r
    vmap_ik = res.vmap

    nv(g_ik) == 0 && return SlicedBranch(MISProblem(g_ik, weights_ik), nothing, r_ik + r0)

    sc0 = contraction_complexity(code, size_dict).sc

    code_ik = update_code(g_ik, code, vmap_ik)        
    re_code_ik = refine(code_ik, size_dict, slicer.refiner, slicer.sc_target, sc0)

    return SlicedBranch(MISProblem(g_ik, weights_ik), re_code_ik, r_ik + r0)
end