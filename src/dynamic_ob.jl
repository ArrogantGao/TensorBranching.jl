# dynamic optimal branching for the maximum independent set problem
# the algorithm has two phases, 1. kernelization 2. slicing via branching + tn contraction

function solve_mis(g::SimpleGraph; reducer::AbstractReducer = TensorNetworkReducer(), slicer::AbstractSlicer, optimizer::CodeOptimizer = TreeSA(), usecuda::Bool = false)
    # 1. kernelization
    g_kernelized = kernelize(g, reducer)

    # 2. optimize the contraction tree
    net = GenericTensorNetwork(IndependentSet(g_kernelized), optimizer = optimizer)
    code = true_eincode(net.code)

    # 3. slicing via branching, each g_slice should be a tuple of (g_slice, contraction_tree)
    g_slices = slice(g_kernelized, code, slicer)

    # 3. contract the slices
    res = [solve(g_slice, usecuda) for g_slice in g_slices]

    # 4. return the maximum results
    return sum(res)
end