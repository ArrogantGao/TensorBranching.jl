# generate slices from the kernelized graph

function slice(g::SimpleGraph, code::DynamicNestedEinsum, r::Int, slicer::AbstractSlicer, reducer::AbstractReducer; verbose::Int = 0, dirname = nothing)
    branch = SlicedBranch(g, code, r)
    return slice(branch, slicer, reducer; verbose = verbose, dirname = dirname)
end

function slice(branch::SlicedBranch, slicer::AbstractSlicer, reducer::AbstractReducer; verbose::Int = 0, dirname = nothing)
    size_dict = uniformsize(uncompress(branch.code), 2)

    (complexity(branch).sc ≤ slicer.sc_target) && return [branch]

    if slicer.search_order == :dfs
        slices = Vector{SlicedBranch{Int}}()
        _slice_dfs!(slices, branch, slicer, reducer, size_dict, verbose)
    elseif slicer.search_order == :bfs
        slices = _slice_bfs(branch, slicer, reducer, size_dict, verbose, dirname)
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

function _slice_bfs(branch::SlicedBranch, slicer::AbstractSlicer, reducer::AbstractReducer, size_dict::Dict{Int, Int}, verbose::Int, dirname)

    saveflag = !isnothing(dirname)
    if saveflag
        !isdir(dirname) && mkdir(dirname)
        df = CSV.write(joinpath(dirname, "slices.csv"), DataFrame(id = Int[], sc = Float64[], tc = Float64[], r = Int[]))
    end

    id = 1

    unfinished_slices = [(branch, reducer)]
    finished_slices = Vector{SlicedBranch{Int}}()
    cc = complexity(branch)
    scs = [cc.sc]
    
    while true

        new_slices, new_scs, new_tcs = __slice_bfs(unfinished_slices, slicer, size_dict, verbose)

        empty!(unfinished_slices)
        empty!(scs)

        saveflag && (ids = Int[])

        for (br, sc, tc) in zip(new_slices, new_scs, new_tcs)
            if sc ≤ slicer.sc_target
                if saveflag
                    CSV.write(df, DataFrame(id = id, sc = sc, tc = tc, r = br[1].r), append = true)
                    push!(ids, id)
                end
                push!(finished_slices, br[1])
                id += 1
            else
                push!(unfinished_slices, br)
                push!(scs, sc)
            end
        end

        if saveflag
            threaded_saveslices(dirname, finished_slices, ids)
            empty!(finished_slices)
            empty!(ids)
        end

        if isempty(unfinished_slices)
            verbose ≥ 1 && @info "all slices finished, $(id - 1) slices in total"
            break
        end

        if verbose ≥ 1
            @info "current num of unfinished slices: $(length(unfinished_slices)), finished slices: $(id - 1)"
            counts = zeros(Int, Int(maximum(scs) - minimum(scs) + 1))
            for sc in scs
                counts[Int(sc - minimum(scs) + 1)] += 1
            end
            println(barplot(Int(minimum(scs)):Int(maximum(scs)), counts, xlabel = "num of slices", ylabel = "sc, target = $(slicer.sc_target)"))
        end

        # force garbage collection
        Base.GC.gc()
    end

    return saveflag ? nothing : finished_slices
end

function __slice_bfs(unfinished_slices::Vector{ST}, slicer::AbstractSlicer, size_dict::Dict{Int, Int}, verbose::Int) where ST

    n = length(unfinished_slices)
    temp_slices = Vector{Tuple{SlicedBranch{Int}, AbstractReducer}}()
    nt = Threads.nthreads()
    chunks = collect(Iterators.partition(1:n, ceil(Int, n/nt)))

    temp_slices = Vector{Vector{Tuple{SlicedBranch{Int}, AbstractReducer}}}(undef, n)

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

    brss = Vector{Tuple{SlicedBranch{Int}, AbstractReducer}}()
    for i in 1:n
        brs = pop!(temp_slices)
        for (br, reducer) in brs
            push!(brss, (br, reducer))
        end
    end

    nbrss = length(brss)
    scs = Vector{Float64}(undef, nbrss)
    tcs = Vector{Float64}(undef, nbrss)

    Threads.@threads for i in collect(Iterators.partition(1:nbrss, ceil(Int, nbrss/nt)))
        for j in i
            branch, reducer = brss[j]
            cc = complexity(branch)
            scs[j] = cc.sc
            tcs[j] = cc.tc
        end
    end

    return brss, scs, tcs
end