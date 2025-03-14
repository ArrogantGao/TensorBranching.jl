# kernlize the graph, via reduction
using OptimalBranching.OptimalBranchingCore: reduce_problem

function kernelize(g::SimpleGraph, reducer::AbstractReducer)
    p_old = MISProblem(g)
    r = MaxSize(0)
    while true
        rp, rv = reduce_problem(MaxSize, p_old, reducer)
        r *= rv
        (p_old == rp) && return (rp.g, TropicalF32(r.size))
        p_old = rp
    end
end