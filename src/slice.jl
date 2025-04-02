# generate slices from the kernelized graph

function slice(g::SimpleGraph, code::DynamicNestedEinsum, r::Int, slicer::AbstractSlicer, reducer::AbstractReducer; verbose::Int = 0)
    branch = SlicedBranch(g, code, r)
    return slice(branch, slicer, reducer; verbose = verbose)
end

function slice(branch::SlicedBranch, slicer::AbstractSlicer, reducer::AbstractReducer; verbose::Int = 0)
    size_dict = uniformsize(uncompress(branch.code), 2)

    (complexity(branch).sc ≤ slicer.sc_target) && return [branch]

    if slicer.search_order == :dfs
        slices = Vector{SlicedBranch{Int}}()
        _slice_dfs!(slices, branch, slicer, reducer, size_dict, verbose)
    elseif slicer.search_order == :bfs
        slices = _slice_bfs(branch, slicer, reducer, size_dict, verbose)
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

    if (complexity(branch).sc ≤ slicer.sc_target)
        push!(slices, branch)
        (verbose ≥ 2) && (@info "current num of slices: $(length(slices))")
        return nothing
    end

    # res is a vector of (mask, code), each corresponding to a slice
    uncompressed_code = uncompress(branch.code)
    region, loss = ob_region(branch.g, uncompressed_code, slicer, slicer.region_selector, size_dict, verbose)
    brs = optimal_branches(branch.g, uncompressed_code, branch.r, slicer, reducer, region, size_dict, verbose)

    for (new_branch, new_reducer) in brs
        _slice_dfs!(slices, new_branch, slicer, new_reducer, size_dict, verbose)
    end

    return nothing
end

function _slice_bfs(branch::SlicedBranch, slicer::AbstractSlicer, reducer::AbstractReducer, size_dict::Dict{Int, Int}, verbose::Int)

    unfinished_slices = [(branch, reducer)]
    finished_slices = SlicedBranch{Int}[]

    while true

        unfinished_slices, new_finished_slices, scs = __slice_bfs(unfinished_slices, slicer, size_dict, verbose)
        append!(finished_slices, new_finished_slices)
        isempty(unfinished_slices) && break

        if verbose ≥ 1
            @info "current num of unfinished slices: $(length(unfinished_slices)), finished slices: $(length(finished_slices))"
            counts = zeros(Int, maximum(scs) - minimum(scs) + 1)
            for sc in scs
                counts[sc - minimum(scs) + 1] += 1
            end
            println(barplot(minimum(scs):maximum(scs), counts, xlabel = "num of slices", ylabel = "sc, target = $(slicer.sc_target)"))
        end

        # force garbage collection
        Base.GC.gc()
    end

    return finished_slices
end

function __slice_bfs(unfinished_slices::Vector{ST}, slicer::AbstractSlicer, size_dict::Dict{Int, Int}, verbose::Int) where ST

    n = length(unfinished_slices)
    temp_slices = Vector{Vector{Tuple{SlicedBranch{Int}, AbstractReducer}}}(undef, n)
    nt = Threads.nthreads()
    chunks = collect(Iterators.partition(1:n, ceil(Int, n/nt)))

    time_start = time()
    Threads.@threads for chunk in chunks
        for i in chunk
            branch, reducer = unfinished_slices[i]
            uncompressed_code = uncompress(branch.code)
            region, loss = ob_region(branch.g, uncompressed_code, slicer, slicer.region_selector, size_dict, verbose)
            brs = optimal_branches(branch.g, uncompressed_code, branch.r, slicer, reducer, region, size_dict, verbose)
            temp_slices[i] = brs
        end
    end
    time_end = time()
    (verbose ≥ 1) && @info "time: $(time_end - time_start), average time: $((time_end - time_start) / n)"

    brss = vcat(temp_slices...)
    nbrss = length(brss)
    scs_0 = Vector{Int}(undef, nbrss)

    Threads.@threads for i in collect(Iterators.partition(1:nbrss, ceil(Int, nbrss/nt)))
        for j in i
            branch, reducer = brss[j]
            scs_0[j] = complexity(branch).sc
        end
    end

    unfinished_slices = Vector{ST}()
    new_finished_slices = Vector{SlicedBranch{Int}}()
    scs = Vector{Int}()
    
    for i in 1:nbrss
        if scs_0[i] ≤ slicer.sc_target
            push!(new_finished_slices, brss[i][1])
        else
            push!(unfinished_slices, brss[i])
            push!(scs, scs_0[i])
        end
    end

    return unfinished_slices, new_finished_slices, scs
end