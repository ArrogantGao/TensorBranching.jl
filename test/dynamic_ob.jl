using TensorBranching
using Graphs
using OptimalBranching.OptimalBranchingMIS
using OptimalBranching.OptimalBranchingMIS.EliminateGraphs

using Test
using Random
Random.seed!(1234)

@testset "dynamic_ob for mis, treesa refiner" begin
    for reducer in [MISReducer(), XiaoReducer(), TensorNetworkReducer()]
        for brancher in [GreedyBrancher(), FixedPointBrancher()]
            for refiner in [TreeSARefiner()]
                for search_order in [:bfs, :dfs]
                    slicer = ContractionTreeSlicer(sc_target = 10, brancher = brancher, refiner = refiner, search_order = search_order)
                    @info "reducer = $reducer, brancher = $brancher, refiner = $refiner, search_order = $search_order"
                    g = random_regular_graph(100, 3)
                    @test dynamic_ob_mis(g, slicer = slicer, reducer = reducer) == mis2(EliminateGraph(g))
                end
            end
        end
    end
end

@testset "dynamic_ob for mis, reoptimize refiner" begin
    reducer = XiaoReducer()
    brancher = GreedyBrancher()
    refiner = ReoptimizeRefiner()
    search_order = :bfs
    slicer = ContractionTreeSlicer(sc_target = 10, brancher = brancher, refiner = refiner, search_order = search_order)
    @info "reducer = $reducer, brancher = $brancher, refiner = $refiner, search_order = $search_order"
    g = random_regular_graph(100, 3)
    @test dynamic_ob_mis(g, slicer = slicer, reducer = reducer) == mis2(EliminateGraph(g))
end