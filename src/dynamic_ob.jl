using OptimalBranching.OptimalBranchingMIS.GenericTensorNetworks: generate_tensors

# dynamic optimal branching for the maximum independent set problem
# the algorithm has two phases, 1. kernelization 2. slicing via branching + tn contraction

function dynamic_ob_mis(g::SimpleGraph, weights::VT; reducer::AbstractReducer = TensorNetworkReducer(), slicer::AbstractSlicer = ContractionTreeSlicer(), optimizer::CodeOptimizer = TreeSA(), usecuda::Bool = false, element_type::Type = Float32, verbose::Int = 0) where VT
    # 1. kernelization
    # the first kernelize step is special, leads to tons of reductions
    res = kernelize(g, weights, reducer, verbose = verbose)
    (verbose ≥ 1) && (@info "kernelized graph: {$(nv(res.g)), $(ne(res.g))} simple graph")
    nv(res.g) == 0 && return res.r 

    # 2. optimize the contraction tree
    code = initialize_code(res.g, optimizer)

    (verbose ≥ 1) && (@info "initial code complexity: \n $(contraction_complexity(code, uniformsize(code, 2)))")

    # 3. slicing via branching, each g_slice should be a tuple of (sliced_graph, contraction_order)
    branches = slice(MISProblem(res.g, res.weights), code, res.r, slicer, reducer, verbose = verbose)

    (verbose ≥ 1) && (@info "branching result: \n num of branches: $(length(branches)). \n max sc: $(maximum(sc.(branches))). \n total tc: $(log2(sum(2 .^ tc.(branches)))).")

    # 4. contract the slices
    res = contract_slices(branches, element_type, usecuda)

    # 5. return the maximum results
    return maximum(res)
end

function solve_slice(branch::SlicedBranch, element_type::Type, usecuda::Bool)
    net = GenericTensorNetwork(IndependentSet(branch.p.g, branch.p.weights), uncompress(branch.code), Dict{Int, Int}())
    res =  Array(solve(net, SizeMax(), T = element_type, usecuda = usecuda))[].n
    return res
end

function contract_slices(branches::Vector{SlicedBranch}, element_type::Type, usecuda::Bool)
    res = element_type[]
    for branch in branches
        if nv(branch.p.g) == 0
            push!(res, element_type(branch.r))
        else
            t = solve_slice(branch, element_type, usecuda)
            res_i = t + element_type(branch.r)
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