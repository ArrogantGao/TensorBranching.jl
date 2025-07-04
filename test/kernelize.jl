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
        for ws in [UnitWeight(n), Float32.(1.0 .+ rand(n))]
            for reducer in [BasicReducer(), XiaoReducer(), TensorNetworkReducer()]
                @info "n = $n, ws = $ws, reducer = $reducer"
                res = kernelize(g, ws, reducer)
                mis_1 = solve(GenericTensorNetwork(IndependentSet(res.g, res.weights), optimizer = TreeSA()), SizeMax(), T = Float32)[].n + res.r
                mis_2 = solve(GenericTensorNetwork(IndependentSet(g, ws), optimizer = TreeSA()), SizeMax(), T = Float32)[].n
                @test nv(res.g) ≤ nv(g)
                @test mis_1 ≈ mis_2
            end
        end
    end
end
