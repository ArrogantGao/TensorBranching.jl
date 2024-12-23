using TensorBranching
using OMEinsum, GenericTensorNetworks, TreeWidthSolver
using Graphs

using Random
Random.seed!(1234)

@testset "tree decomposition" begin
    g0 = random_regular_graph(100, 3)

    prob = GenericTensorNetwork(IndependentSet(g0); optimizer=TreeSA(; ntrials=5, niters=10))
    code = prob.code.eins

    cc = contraction_complexity(code, uniformsize(code, 2^10))
    tree = decompose(code)
    @test width(tree) + 1 ≈ floor(cc.tc / 10)

    @test is_treedecomposition(g0, tree)
end