using TensorBranching
using OMEinsum, OMEinsum.AbstractTrees, TreeWidthSolver
using Graphs
using GenericTensorNetworks, GenericTensorNetworks.TropicalNumbers

using Test
using Random
Random.seed!(1234)

# @testset "tree decomposition" begin
#     g0 = random_regular_graph(100, 3)

#     prob = GenericTensorNetwork(IndependentSet(g0); optimizer=TreeSA(; ntrials=5, niters=10))
#     code = prob.code.eins

#     cc = contraction_complexity(code, uniformsize(code, 2^10))
#     tree = decompose(code)
#     @test width(tree) + 1 ≈ floor(cc.tc / 10)

#     @test is_treedecomposition(g0, tree)

#     eincode = optein"ij, jk, kl, il-> "
#     tree = decompose(eincode)
#     cc = contraction_complexity(eincode, uniformsize(eincode, 2^10))
#     @test width(tree) + 1 ≈ floor(cc.tc / 10)
# end

@testset "tree reform" begin
    for i in 1:1000
        Random.seed!(i)
        n = rand(20:30)
        g = random_regular_graph(n, 4)
        net = GenericTensorNetwork(IndependentSet(g))
        order = net.code
        tensors = GenericTensorNetworks.generate_tensors(TropicalF32(1.0), net)

        remove = rand(1:n, rand(1:10))
        subg, vmap = induced_subgraph(g, setdiff(1:n, remove))
        mis = solve(GenericTensorNetwork(IndependentSet(subg)), SizeMax())[]

        tids = tensors_removed(order, remove)
        sub_order = remove_tensors(order, tids)
        @test sub_order(tensors...)[] == mis

        reindex_tree!(sub_order, vmap)
        @test sub_order(tensors...)[] == mis
    end
end

@testset "tree rethermalize" begin
    g = random_regular_graph(100, 3)
    net = GenericTensorNetwork(IndependentSet(g), optimizer=TreeSA())
    order = net.code.eins
    tensors = GenericTensorNetworks.generate_tensors(TropicalF32(1.0), net)
    
    n1 = neighbors(g, 1) ∪ [1]
    n2 = union([neighbors(g, x) for x in n1]...) ∪ n1
    subg, vmap = induced_subgraph(g, setdiff(1:100, n2))
    sub_order = remove_tensors(order, tensors_removed(order, n2))
    rt_order = rethermalize(sub_order, uniformsize(sub_order, 2))

    ri_order = reindex_tree!(deepcopy(sub_order), vmap)
    rt_ri_order = rethermalize(ri_order, uniformsize(ri_order, 2))
    
    @test rt_order(tensors...)[] == sub_order(tensors...)[] == ri_order(tensors...)[] == rt_ri_order(tensors...)[] == solve(GenericTensorNetwork(IndependentSet(subg)), SizeMax())[]
end