function ob_region(g::SimpleGraph{Int}, code::DynamicNestedEinsum{Int}, slicer::ContractionTreeSlicer, selector::MaxIntersectRS, size_dict::Dict{Int, Int})

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

    return best_region, best_loss
end

function optimal_branches(g::SimpleGraph{Int}, code::DynamicNestedEinsum{Int}, slicer::ContractionTreeSlicer, region::Vector{Int}, size_dict::Dict{Int, Int})
    
    p = MISProblem(g)
    tbl = branching_table(p, slicer.table_solver, region)
    candidates = candidate_clauses(tbl)
    removed_vertices = [removed_vertices(region, g, clause) for clause in candidates]
    codes = [reconfig_code(code, removed_vertices[i], size_dict, slicer.βs, slicer.ntrials, slicer.niters, slicer.sc_target) for i in 1:length(candidates)]

    ## calculate the loss and select the best ones
    optimal_branches = minimize_loss() # should be a vector of Int

    branches = Vector{Tuple{SimpleGraph{Int}, DynamicNestedEinsum{Int}}}()
    for i in optimal_branches
        gi, vmapi = induced_subgraph(g, removed_vertices[i])
        reindex_tree!(codes[i], vmapi)
        push!(branches, (gi, codes[i]))
    end

    return branches
end

function minimize_loss()

end