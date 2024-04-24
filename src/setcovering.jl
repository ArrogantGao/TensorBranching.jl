function max_id(sub_covers::AbstractVector{SubCover{N, T}}) where{N, T}
    m0 = 1
    for sc in sub_covers
        m0 = max(m0, maximum(sc.ids))
    end
    return m0
end

function cover(sub_covers::AbstractVector{SubCover{N, T}}; max_itr::Int = 10, min_complexity::TF = 1.0) where{N, T, TF}
    n = max_id(sub_covers)
    γp = n^(1/N)
    scs_new = copy(sub_covers)
    for i =1:max_itr
        xs = LP_setcover(γp, scs_new, n)
        picked = random_pick(xs, sub_covers, n)
        cx = complexity(picked)
        if (cx < min_complexity) || (i == max_itr)
            return picked, cx
        end
    end
end

function LP_setcover(γp::TF, sub_covers::AbstractVector{SubCover{N, T}}, n::Int) where{N, T, TF}
    fixeds = [count_ones(sc.clause.mask) for sc in sub_covers]
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