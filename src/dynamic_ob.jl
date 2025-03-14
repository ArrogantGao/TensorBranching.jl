# dynamic optimal branching for the maximum independent set problem
# the algorithm has two phases, 1. kernelization 2. slicing via branching + tn contraction

function dynamic_ob_mis(g::SimpleGraph; reducer::AbstractReducer = TensorNetworkReducer(), slicer::AbstractSlicer = ContractionTreeSlicer(), optimizer::CodeOptimizer = TreeSA(), usecuda::Bool = false, kernelize::Bool = true)
    # 1. kernelization
    if kernelize
        g_kernelized, r = kernelize(g, reducer)
        r = r.n
    else
        g_kernelized = g
        r = 0
    end

    # 2. optimize the contraction tree
    net = GenericTensorNetwork(IndependentSet(g_kernelized), optimizer = optimizer)
    code = true_eincode(net.code)
    # tensors = GenericTensorNetworks.generate_tensors(TropicalF32(1.0), IndependentSet(g_kernelized))

    # 3. slicing via branching, each g_slice should be a tuple of (sliced_graph, contraction_order)
    branches = slice(g_kernelized, code, slicer)

    # 3. contract the slices
    res = [solve(branch.code, tensors, usecuda) + branch.r + r for branch in branches]

    # 4. return the maximum results
    return maximum(res)
end