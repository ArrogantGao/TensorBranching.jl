function ip_mwis(graph::SimpleGraph,weights::Vector{Float64}; optimizer = HiGHS.Optimizer)
    model = Model(optimizer)
    set_silent(model)
    nsc = nv(graph)
    @variable(model, 0 <= x[i = 1:nsc] <= 1)
    @objective(model, Max, sum(weights[i]*x[i] for i in 1:nsc))

    for e in edges(graph)
        i = src(e)
        j = dst(e)
        @constraint(model, x[i]+x[j] <= 1)
    end
    cliques = find_all_cliques(graph,3)
    for clique in cliques
        @constraint(model, sum(x[i] for i in clique) <= 1)
    end
    optimize!(model)

    return objective_bound(model)
end

function lp_mwis(graph::SimpleGraph, weights::Vector{Float64}; optimizer = HiGHS.Optimizer)
    model = Model(optimizer)
    set_silent(model)
    nsc = nv(graph)
    @variable(model, 0 <= x[i = 1:nsc] <= 1, Bin)
    @objective(model, Max, sum(weights[i]*x[i] for i in 1:nsc))

    for e in edges(graph)
        i = src(e)
        j = dst(e)
        @constraint(model, x[i]+x[j] <= 1)
    end
    cliques = find_all_cliques(graph,3)
    for clique in cliques
        @constraint(model, sum(x[i] for i in clique) <= 1)
    end
    optimize!(model)
    
    return objective_value(model)
end