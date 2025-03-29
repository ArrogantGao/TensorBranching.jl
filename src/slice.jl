# generate slices from the kernelized graph

function slice(g::SimpleGraph, code::DynamicNestedEinsum, r::Int, slicer::AbstractSlicer, reducer::AbstractReducer; verbose::Int = 0)
    branch = SlicedBranch(g, code, r)
    return slice(branch, slicer, reducer; verbose = verbose)
end

function slice(branch::SlicedBranch, slicer::AbstractSlicer, reducer::AbstractReducer; verbose::Int = 0)
    size_dict = uniformsize(branch.code, 2)

    (contraction_complexity(branch.code, size_dict).sc ≤ slicer.sc_target) && return [branch]

    if slicer.search_order == :dfs
        slices = Vector{SlicedBranch{Int}}()
        _slice_dfs!(slices, branch, slicer, reducer, size_dict, verbose)
    elseif slicer.search_order == :bfs
        unfinished_slices = Vector{Tuple{SlicedBranch{Int}, AbstractReducer}}()
        finished_slices = Vector{SlicedBranch{Int}}()
        push!(unfinished_slices, (branch, reducer))
        slices = _slice_bfs!(unfinished_slices, finished_slices, slicer, size_dict, verbose)
    else
        error("search_order must be :dfs or :bfs")
    end
    return slices
end

function _slice_dfs!(slices::Vector{SlicedBranch{Int}}, branch::SlicedBranch{Int}, slicer::AbstractSlicer, reducer::AbstractReducer, size_dict::Dict{Int, Int}, verbose::Int)

    if nv(branch.g) == 0
        push!(slices, branch)
        return nothing
    end

    if (contraction_complexity(branch.code, size_dict).sc ≤ slicer.sc_target)
        push!(slices, branch)
        (verbose ≥ 2) && (@info "current num of slices: $(length(slices))")
        return nothing
    end

    # res is a vector of (mask, code), each corresponding to a slice
    region, loss = ob_region(branch.g, branch.code, slicer, slicer.region_selector, size_dict, verbose)
    brs = optimal_branches(branch.g, branch.code, slicer, reducer, region, size_dict, verbose)

    for (new_branch, new_reducer) in brs
        _slice_dfs!(slices, SlicedBranch(new_branch.g, new_branch.code, branch.r + new_branch.r), slicer, new_reducer, size_dict, verbose)
    end

    return nothing
end

function _slice_bfs!(unfinished_slices::Vector{Tuple{SlicedBranch{Int}, AbstractReducer}}, finished_slices::Vector{SlicedBranch{Int}}, slicer::AbstractSlicer, size_dict::Dict{Int, Int}, verbose::Int)

    isempty(unfinished_slices) && return finished_slices

    if verbose ≥ 1
        scs = Int.([contraction_complexity(branch.code, size_dict).sc for (branch, reducer) in unfinished_slices])
        @info "current num of unfinished slices: $(length(unfinished_slices)), largest sc: $(maximum(scs)), smallest sc: $(minimum(scs))"
        @info "current num of finished slices: $(length(finished_slices))"
        counts = zeros(Int, maximum(scs) - minimum(scs) + 1)
        for sc in scs
            counts[sc - minimum(scs) + 1] += 1
        end
        println(barplot(minimum(scs):maximum(scs), counts, xlabel = "num of slices", ylabel = "sc"))
    end

    n = length(unfinished_slices)

    pb = (verbose ≥ 1) ? ProgressBar(1:n) : 1:n
    for _ in pb
        branch, reducer = popfirst!(unfinished_slices)
        region, loss = ob_region(branch.g, branch.code, slicer, slicer.region_selector, size_dict, verbose)
        brs = optimal_branches(branch.g, branch.code, slicer, reducer, region, size_dict, verbose)
        for (new_branch, new_reducer) in brs
            new_slice = SlicedBranch(new_branch.g, new_branch.code, branch.r + new_branch.r)
            if (complexity(new_slice).sc ≤ slicer.sc_target)
                push!(finished_slices, new_slice)
            else
                push!(unfinished_slices, (new_slice, new_reducer))
            end
        end
    end

    return _slice_bfs!(unfinished_slices, finished_slices, slicer, size_dict, verbose)
end