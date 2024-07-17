"""
    complexity(sbranches::Vector{Int})

Compute the complexity of a set of branches.

This function solves the equation 1 = sum(x^(-i) for i in sbranches) for x, where `sbranches` is a vector of integers representing the branches.

# Arguments
- `sbranches::Vector{Int}`: A vector of integers representing the branches.

# Returns
- `Float64`: The value of x that satisfies the equation.

"""
function complexity(sbranches::Vector{Int})
    # solve x, where 1 = sum(x^(-i) for i in sbranches)
    f = x -> sum(x[1]^(-i) for i in sbranches) - 1.0
    sol = nlsolve(f, [1.0])
    return sol.zero[1]
end

function complexity(subcovers::AbstractVector{SubCover{INT}}) where{INT}
    return complexity([sc.n_rm for sc in subcovers])
end

function max_id(sub_covers::AbstractVector{SubCover{INT}}) where{INT}
    m0 = 1
    for sc in sub_covers
        m0 = max(m0, maximum(sc.ids))
    end
    return m0
end

function γ0(sub_covers::AbstractVector{SubCover{INT}}) where{INT}
    n = max_id(sub_covers)
    max_rvs = Int[]
    for i in 1:n
        rvs = [sc.n_rm for sc in sub_covers if sc.ids == Set([i])]
        push!(max_rvs, maximum(rvs))
    end
    return complexity(max_rvs), n
end

"""
    cover(sub_covers::AbstractVector{SubCover{INT}}; max_itr::Int = 2, min_complexity::TF = 1.0) where{N, T, TF}

The `cover` function performs a set covering algorithm on a collection of sub-covers.

## Arguments
- `sub_covers::AbstractVector{SubCover{INT}}`: A collection of sub-covers.
- `max_itr::Int = 2`: The maximum number of iterations to perform.
- `min_complexity::TF = 1.0`: The minimum complexity threshold.

## Returns
- `picked`: The selected sub-covers.
- `cx`: The complexity of the selected sub-covers.

"""
function cover(sub_covers::AbstractVector{SubCover{INT}}, max_itr::Int; min_complexity::TF = 1.0) where{INT, TF}
    cx, n = γ0(sub_covers)
    scs_new = copy(sub_covers)
    for i =1:max_itr
        xs = LP_setcover(cx, scs_new, n)
        picked = random_pick(xs, sub_covers, n)
        cx = complexity(picked)
        @debug "Iteration $i: picked_num = $length(picked), complexity = $cx"
        if (cx < min_complexity) || (i == max_itr)
            return picked, cx
        end
    end
end

function LP_setcover(γp::TF, sub_covers::AbstractVector{SubCover{INT}}, n::Int) where{INT, TF}
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

function random_pick(xs::Vector{TF}, sub_covers::AbstractVector{SubCover{INT}}, n::Int) where{INT, TF}
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

function setcover_strategy(tbl::BranchingTable{INT}, vertices::Vector{Int}, g::SimpleGraph, max_itr::Int, measure::AbstractMeasure) where{INT}
    sub_covers = subcovers(tbl, vertices, g, measure)
    cov, cx = cover(sub_covers, max_itr)
    branches = Branches([Branch(sc.clause, vertices, g, measure) for sc in cov])
    return branches
end