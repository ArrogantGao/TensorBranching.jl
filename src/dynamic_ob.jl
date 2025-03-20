using OptimalBranching.OptimalBranchingMIS.GenericTensorNetworks: generate_tensors

# dynamic optimal branching for the maximum independent set problem
# the algorithm has two phases, 1. kernelization 2. slicing via branching + tn contraction

function dynamic_ob_mis(g::SimpleGraph; reducer::AbstractReducer = TensorNetworkReducer(), slicer::AbstractSlicer = ContractionTreeSlicer(), optimizer::CodeOptimizer = TreeSA(), usecuda::Bool = false, element_type::Type = Float32, verbose::Int = 0)
    # 1. kernelization
    # the first kernelize step is special, leads to tons of reductions
    g_kernelized, r, _ = kernelize(g, reducer, verbose = verbose)
    (verbose ≥ 1) && (@info "kernelized graph: {$(nv(g_kernelized)), $(ne(g_kernelized))} simple graph")
    nv(g_kernelized) == 0 && return r 

    # 2. optimize the contraction tree
    code = initialize_code(g_kernelized, optimizer)

    return dynamic_ob_mis(g_kernelized, code, r; reducer = reducer, slicer = slicer, usecuda = usecuda, element_type = element_type, verbose = verbose)
end

# if one already has the code of the kernelized graph, then one can directly use this function
function dynamic_ob_mis(g::SimpleGraph, code::DynamicNestedEinsum, r::Int; reducer::AbstractReducer = TensorNetworkReducer(), slicer::AbstractSlicer = ContractionTreeSlicer(), usecuda::Bool = false, element_type::Type = Float32, verbose::Int = 0)

    (verbose ≥ 1) && (@info "initial code complexity: \n $(contraction_complexity(code, uniformsize(code, 2)))")

    # 3. slicing via branching, each g_slice should be a tuple of (sliced_graph, contraction_order)
    branches = slice(SlicedBranch(g, code, r), slicer, reducer, verbose = verbose)

    (verbose ≥ 1) && (@info "branching result: \n num of branches: $(length(branches)). \n max sc: $(maximum(sc.(branches))). \n total tc: $(log2(sum(2 .^ tc.(branches)))).")

    # 4. contract the slices
    res = contract(branches, element_type, usecuda)

    # 5. return the maximum results
    return maximum(res)
end

function contract(branches::Vector{SlicedBranch{T1}}, element_type::Type, usecuda::Bool) where {T1}
    res = Int[]
    for branch in branches
        if nv(branch.g) == 0
            push!(res, branch.r)
        else
            net = GenericTensorNetwork(IndependentSet(branch.g), branch.code, Dict{Int, Int}())
            t = Array(solve(net, SizeMax(), T = element_type, usecuda = usecuda))[].n
            res_i = Int(t) + branch.r
            push!(res, res_i)
        end
    end
    return res
end

function initialize_code(g::SimpleGraph, optimizer::CodeOptimizer)
    net = GenericTensorNetwork(IndependentSet(g), optimizer = optimizer)
    code = true_eincode(net.code)
    return code
end