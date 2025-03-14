# dynamic optimal branching for the maximum independent set problem
# the algorithm has two phases, 1. kernelization 2. slicing via branching + tn contraction

function dynamic_ob_mis(g::SimpleGraph; reducer::AbstractReducer = TensorNetworkReducer(), slicer::AbstractSlicer = ContractionTreeSlicer(), optimizer::CodeOptimizer = TreeSA(), usecuda::Bool = false, use_kernelize::Bool = true, tensor_type::Type = TropicalF32)
    # 1. kernelization
    if use_kernelize
        g_kernelized, r = kernelize(g, reducer)
        r = r.n
    else
        g_kernelized = g
        r = 0
    end

    nv(g_kernelized) == 0 && return r

    # 2. optimize the contraction tree
    net = GenericTensorNetwork(IndependentSet(g_kernelized), optimizer = optimizer)
    code = true_eincode(net.code)
    tensors = GenericTensorNetworks.generate_tensors(tensor_type(1.0), IndependentSet(g_kernelized))
    usecuda && (tensors = gpu_tensors(tensors))

    # 3. slicing via branching, each g_slice should be a tuple of (sliced_graph, contraction_order)
    branches = slice(g_kernelized, code, slicer)
    @show branches

    # 3. contract the slices
    res = [item(branch.code(tensors...)).n + branch.r + r for branch in branches]

    # 4. return the maximum results
    return Int(maximum(res))
end

item(t::Array{T, 0}) where {T} = t[]

function gpu_tensors(tensors) where {T}
    error("CUDA should be loaded as a external package.")
end