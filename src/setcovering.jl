function sbranches(clauses::Vector{Clause{N, T}}) where{N, T}
    return [count_ones(c.mask) for c in clauses]
end

function complexity(sbranches::Vector{Int})
    # solve x, where 1 = sum(x^(-i) for i in sbranches)
    f = x -> sum(x[1]^(-i) for i in sbranches) - 1.0
    sol = nlsolve(f, [1.0])
    return sol.zero[1]
end

function complexity(subcovers::AbstractVector{SubCover{N, T}}) where{N, T}
    return complexity([sc.n_rm for sc in subcovers])
end

function max_id(sub_covers::AbstractVector{SubCover{N, T}}) where{N, T}
    m0 = 1
    for sc in sub_covers
        m0 = max(m0, maximum(sc.ids))
    end
    return m0
end

function cover(sub_covers::AbstractVector{SubCover{N, T}}; max_itr::Int = 2, min_complexity::TF = 1.0) where{N, T, TF}
    n = max_id(sub_covers)
    γp = n^(1/N)
    scs_new = copy(sub_covers)
    for i =1:max_itr
        xs = LP_setcover(γp, scs_new, n)
        picked = random_pick(xs, sub_covers, n)
        cx = complexity(picked)
        @debug "Iteration $i: picked_num = $length(picked), complexity = $cx"
        if (cx < min_complexity) || (i == max_itr)
            return picked, cx
        end
    end
end

function LP_setcover(γp::TF, sub_covers::AbstractVector{SubCover{N, T}}, n::Int) where{N, T, TF}
    fixeds = [sc.n_rm for sc in sub_covers]
    γs = 1 ./ γp .^ fixeds
    nsc = length(sub_covers)

    sets_id = [Vector{Int}() for _=1:n]
    for i in 1:nsc
        for j in sub_covers[i].ids
            push!(sets_id[j], i)
        end
    end

    # LP by JuMP
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, 0 <= x[i = 1:nsc] <= 1)
    @objective(model, Min, sum(x[i] * γs[i] for i in 1:nsc))
    for i in 1:n
        @constraint(model, sum(x[j] for j in sets_id[i]) >= 1)
    end

    optimize!(model)
    @assert is_solved_and_feasible(model)

    return [value(x[i]) for i in 1:nsc]
end

function random_pick(xs::Vector{TF}, sub_covers::AbstractVector{SubCover{N, T}}, n::Int) where{N, T, TF}
    picked = Set{Int}()
    picked_ids = Set{Int}()
    nsc = length(sub_covers)
    flag = true
    while flag 
        for i in 1:nsc
            if (rand() < xs[i]) && !(i in picked)
                push!(picked, i)
                picked_ids = union(picked_ids, sub_covers[i].ids)
            end
            if length(picked_ids) == n
                flag = false
                break
            end
        end
    end

    return [sub_covers[i] for i in picked]
end

function setcover_strategy(tbl::BranchingTable{N, C}, vertices::Vector{Int}, g::SimpleGraph; max_itr::Int = 1) where{N, C}
    sub_covers = subcovers(tbl, vertices, g)
    cov, cx = cover(sub_covers, max_itr=max_itr)
    return DNF([c.clause for c in cov])
end