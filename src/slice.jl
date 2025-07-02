# generate slices from the kernelized graph
# two different ways are provided: bfs and bfs_rw

function slice(p::MISProblem{INT, VT}, code::DynamicNestedEinsum, r::RT, slicer::AbstractSlicer, reducer::AbstractReducer; verbose::Int = 0) where {INT, VT, RT}
    branch = SlicedBranch(p, code, r)
    if (complexity(branch).sc ≤ slicer.sc_target)
        return [branch] # no need to slice
    else
        return slice_bfs(branch, slicer, reducer, verbose) # slice the branch
    end
end

function slice_bfs(branch::SlicedBranch{INT, VT, RT}, slicer::AbstractSlicer, reducer::AbstractReducer, verbose::Int) where {INT, VT, RT}
    unfinished_slices = SlicedBranch[branch]
    finished_slices = SlicedBranch[]
    scs = [complexity(branch).sc]
    size_dict = uniformsize(uncompress(branch.code), 2)

    while !isempty(unfinished_slices)
        verbose ≥ 1 && show_status(scs, slicer.sc_target, length(unfinished_slices), length(finished_slices))
        new_slices, new_scs = _slice_bfs(unfinished_slices, slicer, reducer, size_dict, verbose)
        empty!(unfinished_slices)
        empty!(scs)
        for (slice, sc) in zip(new_slices, new_scs)
            if sc ≤ slicer.sc_target
                push!(finished_slices, slice)
            else
                push!(unfinished_slices, slice)
                push!(scs, sc)
            end
        end
    end

    return finished_slices
end

function _slice_bfs(unfinished_slices::Vector{ST}, slicer::AbstractSlicer, reducer::AbstractReducer, size_dict::Dict{Int, Int}, verbose::Int) where ST

    n = length(unfinished_slices)
    temp_slices = Vector{Vector{ST}}(undef, n)
    nt = Threads.nthreads()
    chunks = collect(Iterators.partition(1:n, ceil(Int, n/nt)))

    time_start = time()
    Threads.@threads for chunk in chunks
        for i in chunk
            branch = unfinished_slices[i]
            uncompressed_code = uncompress(branch.code)
            region, loss = ob_region(branch.p.g, uncompressed_code, slicer, slicer.region_selector, size_dict, verbose)
            branches = optimal_branches(branch.p, uncompressed_code, branch.r, slicer, reducer, region, size_dict, verbose)
            temp_slices[i] = branches
        end
    end
    time_end = time()
    (verbose ≥ 1) && @info "time: $(time_end - time_start), average time: $((time_end - time_start) / n)"

    new_slices = union(temp_slices...)
    new_scs = [complexity(slice).sc for slice in new_slices]

    return new_slices, new_scs
end

function slice_bfs_rw(branch::SlicedBranch{INT, VT, RT}, slicer::AbstractSlicer, reducer::AbstractReducer, dirname::String, verbose::Int) where {INT, VT, RT}
    size_dict = uniformsize(uncompress(branch.code), 2)

    df = CSV.write(joinpath(dirname, "slices.csv"), DataFrame(id = Int[], sc = Float64[], tc = Float64[], r = RT[]))

    if complexity(branch).sc ≤ slicer.sc_target
        cc = complexity(branch)
        CSV.write(df, DataFrame(id = 1, sc = cc.sc, tc = cc.tc, r = branch.r), append = true)
        save_finished(dirname, branch, 1)
        return nothing
    end

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

        info_finished = [Vector{Tuple{Int, Float64, Float64, RT}}() for _ in 1:num_unfinished]
        info_unfinished = [Vector{Tuple{Int, Float64, Float64, RT}}() for _ in 1:num_unfinished]

        nt = Threads.nthreads()
        chunks = collect(Iterators.partition(1:num_unfinished, ceil(Int, num_unfinished/nt)))

        time_start = time()
        Threads.@threads for chunk in chunks
            for i in chunk
                # load the branch
                new_branch = round == 1 ? branch : load_unfinished(old_unfinished_dir, flatten_unfinished[i][1])

                # generate new branches
                uncompressed_code = uncompress(new_branch.code)
                region, loss = ob_region(new_branch.p.g, uncompressed_code, slicer, slicer.region_selector, size_dict, verbose)
                branches = optimal_branches(new_branch.p, uncompressed_code, new_branch.r, slicer, reducer, region, size_dict, verbose)

                # update info, save finished or unfinished
                for new_branch in branches
                    id = Int(time_ns())
                    cc = complexity(new_branch)
                    if cc.sc ≤ slicer.sc_target
                        push!(info_finished[i], (id, cc.sc, cc.tc, new_branch.r))
                        save_finished(dirname, new_branch, id)
                    else
                        push!(info_unfinished[i], (id, cc.sc, cc.tc, new_branch.r))
                        save_unfinished(unfinished_dir, new_branch, id)
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