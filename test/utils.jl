using TensorBranching
using Graphs, TropicalNumbers, OMEinsum

using OptimalBranching
using OptimalBranching.OptimalBranchingCore, OptimalBranching.OptimalBranchingMIS
using OptimalBranching.OptimalBranchingMIS.EliminateGraphs

using GenericTensorNetworks
using GenericTensorNetworks: generate_tensors

using Test
using Random
Random.seed!(1234)

using TensorBranching: remove_tensors, remove_tensors!, tensors_removed, unsafe_flatten, rethermalize, reindex_tree!

@testset "tree reform" begin
    for i in 1:1000
        Random.seed!(i)
        n = 30
        g = random_regular_graph(n, 3)
        net = GenericTensorNetwork(IndependentSet(g))
        order = net.code
        tensors = generate_tensors(TropicalF32(1.0), net)

        # check if the result is still correct after removing for more than one rounds
        for _ in 1:4
            remove = rand(1:nv(g), rand(1:5))
            subg, vmap = induced_subgraph(g, setdiff(1:nv(g), remove))
            mis = mis2(EliminateGraph(subg))

            tids = tensors_removed(order, remove)
            sub_order = remove_tensors(order, tids)
            @test sub_order(tensors...)[].n ≈ mis

            ri_order = reindex_tree!(sub_order, vmap)
            @test ri_order(tensors...)[].n ≈ mis
            g = subg
            order = ri_order
        end
    end
end

@testset "tree rethermalize" begin
    for i in 1:10
        n = 60
        g = random_regular_graph(n, 3)
        net = GenericTensorNetwork(IndependentSet(g), optimizer=TreeSA())
        order = net.code
        tensors = GenericTensorNetworks.generate_tensors(TropicalF32(1.0), net)
        
        n1 = neighbors(g, 1) ∪ [1]
        n2 = union([neighbors(g, x) for x in n1]...) ∪ n1
        subg, vmap = induced_subgraph(g, setdiff(1:n, n2))
        sub_order = remove_tensors(order, tensors_removed(order, n2))
        rt_order = rethermalize(deepcopy(sub_order), uniformsize(sub_order, 2), 100.0:100.0, 1, 10, 25)

        ri_order = reindex_tree!(deepcopy(sub_order), vmap)
        rt_ri_order = rethermalize(deepcopy(ri_order), uniformsize(ri_order, 2), 100.0:100.0, 1, 10, 25)
        
        @test rt_order(tensors...)[].n ≈ sub_order(tensors...)[].n ≈ ri_order(tensors...)[].n ≈ rt_ri_order(tensors...)[].n ≈ mis2(EliminateGraph(subg))
    end
end