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
    elseif slicer.search_order == :tree
        slices = _slice_tree(branch, slicer, reducer, size_dict, verbose)
    elseif slicer.search_order == :bfs_rw
        isnothing(dirname) && error("dirname must be provided for bfs_rw")
        !isdir(dirname) && mkdir(dirname)
        slices = _slice_bfs_rw(branch, slicer, reducer, size_dict, verbose, dirname)
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


function _slice_bfs_rw(branch::SlicedBranch, slicer::AbstractSlicer, reducer::AbstractReducer, size_dict::Dict{Int, Int}, verbose::Int, dirname)
    df = CSV.write(joinpath(dirname, "slices.csv"), DataFrame(id = Int[], sc = Float64[], tc = Float64[], r = Int[]))

    num_unfinished = 1
    num_finished = 0
    round = 1

    local flatten_unfinished
    local unfinished_dir

    time_start = time()
    while true
        old_unfinished_dir = joinpath(dirname, "unfinished_r$(round - 1)")
        unfinished_dir = joinpath(dirname, "unfinished_r$(round)")
        !isdir(unfinished_dir) && mkdir(unfinished_dir)

        info_finished = [Vector{Tuple{Int, Float64, Float64, Int}}() for _ in 1:num_unfinished]
        info_unfinished = [Vector{Tuple{Int, Float64, Float64, Int}}() for _ in 1:num_unfinished]

        nt = Threads.nthreads()
        chunks = collect(Iterators.partition(1:num_unfinished, ceil(Int, num_unfinished/nt)))

        time_start = time()
        Threads.@threads for chunk in chunks
            for i in chunk
                # load the branch
                (new_branch, new_reducer) = round == 1 ? (branch, reducer) : load_unfinished(old_unfinished_dir, flatten_unfinished[i][1])

                # generate new branches
                uncompressed_code = uncompress(new_branch.code)
                region, loss = ob_region(new_branch.g, uncompressed_code, slicer, slicer.region_selector, size_dict, verbose)
                brs = optimal_branches(new_branch.g, uncompressed_code, new_branch.r, slicer, new_reducer, region, size_dict, verbose)

                # update info, save finished or unfinished
                for (nb, nr) in brs
                    id = Int(time_ns())
                    cc = complexity(nb)
                    if cc.sc ≤ slicer.sc_target
                        push!(info_finished[i], (id, cc.sc, cc.tc, nb.r))
                        save_finished(dirname, nb, id)
                    else
                        push!(info_unfinished[i], (id, cc.sc, cc.tc, nb.r))
                        save_unfinished(unfinished_dir, nb, nr, id)
                    end
                end
            end
        end
        time_end = time()
        (verbose ≥ 1) && @info "round: $round, time: $(time_end - time_start), average time: $((time_end - time_start) / num_unfinished)"

        for i in 1:num_unfinished
            for e in info_finished[i]
                (id, sc, tc, r) = e
                CSV.write(df, DataFrame(id = id, sc = sc, tc = tc, r = r), append = true)
                num_finished += 1
            end
        end

        # clean the old unfinished directory
        round > 1 && rm(old_unfinished_dir, recursive = true)

        flatten_unfinished = vcat(info_unfinished...)
        isempty(flatten_unfinished) && break

        num_unfinished = length(flatten_unfinished)
        scs = [e[2] for e in flatten_unfinished]
        verbose ≥ 1 && show_status(scs, slicer.sc_target, num_unfinished, num_finished)

        round += 1
    end
    rm(unfinished_dir, recursive = true)

    return nothing
end

function slice_tree(g::SimpleGraph, code::DynamicNestedEinsum, r::Int, slicer::AbstractSlicer, reducer::AbstractReducer; verbose::Int = 0)
    branch = SlicedBranch(g, code, r)
    size_dict = uniformsize(code, 2)
    return _slice_tree(branch, slicer, reducer, size_dict, verbose)
end

function _slice_tree(branch::SlicedBranch, slicer::AbstractSlicer, reducer::AbstractReducer, size_dict::Dict{Int, Int}, verbose::Int)
    tree = SlicingTree(branch, SlicingTree[])
    if nv(branch.g) == 0 || (complexity(branch).sc ≤ slicer.sc_target) || isnothing(branch.code)
        return tree
    end

    uncompressed_code = uncompress(branch.code)
    region, loss = ob_region(branch.g, uncompressed_code, slicer, slicer.region_selector, size_dict, verbose)
    brs = optimal_branches(branch.g, uncompressed_code, branch.r, slicer, reducer, region, size_dict, verbose)

    for (new_branch, new_reducer) in brs
        push!(tree.children, _slice_tree(new_branch, slicer, new_reducer, size_dict, verbose))
    end

    return tree
end

