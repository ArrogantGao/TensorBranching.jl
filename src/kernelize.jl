# kernlize the graph
using OptimalBranchingCore: reduce_problem

function kernelize(g::SimpleGraph, reducer::AbstractReducer)
    p_old = MISProblem(g)
    while true
        rp, rv = reduce_problem(MaxSize, p_old, reducer)
        (p_old == rp) && return rp.g
        p_old = rp
    end
end