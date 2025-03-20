using TensorBranching
using OptimalBranching.OptimalBranchingMIS
using OptimalBranching.OptimalBranchingMIS.EliminateGraphs


using Test
using Random
Random.seed!(1234)

@testset "kernelize" begin
    for n in 30:10:60
        g = random_regular_graph(n, 3)
        for reducer in [MISReducer(), XiaoReducer(), TensorNetworkReducer()]
            g_new, r, vmap = kernelize(g, reducer)
            mis_1 = mis2(EliminateGraph(g_new)) + r
            mis_2 = mis2(EliminateGraph(g))
            @test nv(g_new) ≤ nv(g)
            @test mis_1 == mis_2
        end
    end
end
