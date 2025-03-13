# dynamic optimal branching for the maximum independent set problem
# the algorithm has two phases, 1. kernelization 2. slicing via branching + tn contraction

function solve_mis(g::SimpleGraph; reducer::AbstractReducer = TensorNetworkReducer(), slicer::AbstractSlicer, optimizer::CodeOptimizer = TreeSA(), usecuda::Bool = false)
    # 1. kernelization
    g_kernelized, r = kernelize(g, reducer)

    # 2. optimize the contraction tree
    net = GenericTensorNetwork(IndependentSet(g_kernelized), optimizer = optimizer)
    code = true_eincode(net.code)
    # tensors = GenericTensorNetworks.generate_tensors(TropicalF32(1.0), IndependentSet(g_kernelized))

    # 3. slicing via branching, each g_slice should be a tuple of (sliced_graph, contraction_order)
    g_slices = slice(g_kernelized, code, slicer)

    # 3. contract the slices
    res = [solve(code, tensors, usecuda) for (g, code) in g_slices]

    # 4. return the maximum results
    return (sum(res) * r).n
end