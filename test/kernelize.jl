using TensorBranching
using OptimalBranching.OptimalBranchingMIS
using OptimalBranching.OptimalBranchingMIS.EliminateGraphs
using Graphs
using GenericTensorNetworks

using Test
using Random
Random.seed!(1234)

@testset "kernelize" begin
    for n in 30:10:60
        g = random_regular_graph(n, 3)
        for reducer in [MISReducer(), XiaoReducer(), TensorNetworkReducer()]
            g_new, r, vmap, reducer_new = kernelize(g, reducer)
            mis_1 = mis2(EliminateGraph(g_new)) + r
            mis_2 = mis2(EliminateGraph(g))
            @test nv(g_new) ≤ nv(g)
            @test mis_1 == mis_2
        end
    end

    for g in [random_regular_graph(300, 3), SimpleGraph(GenericTensorNetworks.random_diagonal_coupled_graph(40, 40, 0.8))]
        reducer = TensorNetworkReducer()
        g_new, _, _, reducer_new = kernelize(g, reducer)
        g_new_new, _, _, reducer_new_new = kernelize(g_new, reducer_new)
        @test g_new == g_new_new
        @test length(reducer_new.region_list) == nv(g_new)
        @test length(reducer_new_new.region_list) == nv(g_new_new)
    end
end
