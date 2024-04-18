function sbranches(clauses::Vector{Clause{N, T}}) where{N, T}
    return [count_ones(c.mask) for c in clauses]
end

function complexity(sbranches::Vector{Int})
    # solve x, where 1 = sum(x^(-i) for i in sbranches)
    f = x -> sum(x[1]^(-i) for i in sbranches) - 1.0
    sol = nlsolve(f, [1.0])
    return sol.zero
end