# dynamic optimal branching for the maximum independent set problem
# the algorithm has two phases, 1. kernelization 2. slicing via branching + tn contraction

function solve_mis(g::SimpleGraph; reducer::AbstractReducer = TensorNetworkReducer(), usecuda::Bool = false, loss::Symbol = :max_sc)
    # 1. kernelization
    g_kernelized = kernelize(g, reducer)

    # 2. slicing via branching
    g_slices = dynamic_ob_slicing(g_kernelized)

    # 3. contract the slices
    res = [solve(g_slice) for g_slice in g_slices]

    # 4. return the maximum results
    return sum(res)
end