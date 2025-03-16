using OptimalBranching.OptimalBranchingMIS.GenericTensorNetworks: generate_tensors

# dynamic optimal branching for the maximum independent set problem
# the algorithm has two phases, 1. kernelization 2. slicing via branching + tn contraction

function dynamic_ob_mis(g::SimpleGraph; reducer::AbstractReducer = TensorNetworkReducer(), slicer::AbstractSlicer = ContractionTreeSlicer(), optimizer::CodeOptimizer = TreeSA(), usecuda::Bool = false, use_kernelize::Bool = true, tensor_type::Type = TropicalF32, verbose::Int = 0)
    # 1. kernelization
    if use_kernelize
        g_kernelized, r = kernelize(g, reducer)
        r = Int(r.n)
        (verbose ≥ 1) && (@info "kernelized graph: {$(nv(g_kernelized)), $(ne(g_kernelized))} simple graph")
        nv(g_kernelized) == 0 && return r
    else
        g_kernelized = g
        r = 0
    end

    # 2. optimize the contraction tree
    code = initialize_code(g_kernelized, optimizer)

    return dynamic_ob_mis(g_kernelized, code, r; slicer = slicer, usecuda = usecuda, tensor_type = tensor_type, verbose = verbose)
end

# if one already has the code of the kernelized graph, then one can directly use this function
function dynamic_ob_mis(g::SimpleGraph, code::DynamicNestedEinsum, r::Int; slicer::AbstractSlicer = ContractionTreeSlicer(), usecuda::Bool = false, tensor_type::Type = TropicalF32, verbose::Int = 0)

    (verbose ≥ 1) && (@info "initial code complexity: \n $(contraction_complexity(code, uniformsize(code, 2)))")

    tensors = initialize_tensors(g, usecuda, tensor_type)

    # 3. slicing via branching, each g_slice should be a tuple of (sliced_graph, contraction_order)
    branches = slice(SlicedBranch(g, code, r), slicer, verbose)

    (verbose ≥ 1) && (@info "branching result: \n num of branches: $(length(branches)). \n max sc: $(maximum(sc.(branches))). \n total tc: $(log2(sum(2 .^ tc.(branches)))).")

    # 4. contract the slices
    res = contract(branches, tensors)

    # 5. return the maximum results
    return maximum(res)
end

function initialize_code(g::SimpleGraph, optimizer::CodeOptimizer)
    net = GenericTensorNetwork(IndependentSet(g), optimizer = optimizer)
    code = true_eincode(net.code)
    return code
end

function initialize_tensors(g::SimpleGraph, usecuda::Bool, tensor_type::Type)
    tensors = generate_tensors(tensor_type(1.0), IndependentSet(g))
    usecuda && (tensors = gpu_tensors(tensors))
    return tensors
end

function contract(branches::Vector{SlicedBranch{T1}}, tensors::Vector{T2}) where {T1, T2}
    res = [Int(item(branch.code(tensors...)).n) + branch.r for branch in branches]
    return res
end

item(t::Array{T, 0}) where {T} = t[]

function gpu_tensors(x)
    error("CUDA should be loaded as a external package.")
end