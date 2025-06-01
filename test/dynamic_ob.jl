using TensorBranching
using Graphs
using OptimalBranching.OptimalBranchingMIS
using OptimalBranching.OptimalBranchingMIS.EliminateGraphs
using GenericTensorNetworks
using Test
using Random
using UnitDiskMapping
Random.seed!(1234)

@testset "dynamic_ob for mis, treesa refiner" begin
    for reducer in [MISReducer(), XiaoReducer()]
        for brancher in [GreedyBrancher(), FixedPointBrancher()]
            for refiner in [TreeSARefiner()]
                for selector in [ScoreRS(loss = :num_uniques), ScoreRS(loss = :sc_score), ScoreRS(loss = :bag_score)]
                    for search_order in [:bfs, :dfs]
                        slicer = ContractionTreeSlicer(sc_target = 12, brancher = brancher, refiner = refiner, search_order = search_order, region_selector = selector)
                        @info "reducer = $reducer, brancher = $brancher, refiner = $refiner, search_order = $search_order, region_selector = $selector"
                        g = random_regular_graph(100, 3)
                        @test dynamic_ob_mis(g, slicer = slicer, reducer = reducer) == mis2(EliminateGraph(g))
                    end
                end
            end
        end
    end
end

@testset "dynamic_ob for mis, reoptimize refiner" begin
    reducer = TensorNetworkReducer()
    brancher = GreedyBrancher()
    refiner = ReoptimizeRefiner()
    search_order = :bfs
    slicer = ContractionTreeSlicer(sc_target = 10, brancher = brancher, refiner = refiner, search_order = search_order)
    @info "reducer = $reducer, brancher = $brancher, refiner = $refiner, search_order = $search_order"
    g = random_regular_graph(100, 3)
    @test dynamic_ob_mis(g, slicer = slicer, reducer = reducer, verbose = 2) == mis2(EliminateGraph(g))
end

@testset "dynamic_ob for mwis, treesa refiner" begin
    for reducer in [MWISReducer(), TensorNetworkReducer()]
        for brancher in [GreedyBrancher(), FixedPointBrancher()]
            for refiner in [TreeSARefiner()]
                for selector in [ScoreRS(loss = :num_uniques), ScoreRS(loss = :sc_score), ScoreRS(loss = :bag_score)]
                    for search_order in [:bfs, :dfs]
                        slicer = ContractionTreeSlicer(sc_target = 12, brancher = brancher, refiner = refiner, search_order = search_order, region_selector = selector)
                        @info "reducer = $reducer, brancher = $brancher, refiner = $refiner, search_order = $search_order, region_selector = $selector"
                        g = random_regular_graph(100, 3)
                        weights = ones(Float64, nv(g))
                        problem = GenericTensorNetwork(IndependentSet(g, weights); optimizer = TreeSA())
                        mwis = solve(problem, SizeMax())[1].n
                        @test abs(dynamic_ob_mwis(g, weights, slicer = slicer, reducer = reducer) - mwis) < 1e-12
                    end
                end
            end
        end
    end
end

@testset "dynamic_ob for mwis, reoptimize refiner" begin
    reducer = TensorNetworkReducer()
    brancher = GreedyBrancher()
    refiner = ReoptimizeRefiner()
    search_order = :bfs
    slicer = ContractionTreeSlicer(sc_target = 10, brancher = brancher, refiner = refiner, search_order = search_order)
    @info "reducer = $reducer, brancher = $brancher, refiner = $refiner, search_order = $search_order"
    g = random_regular_graph(100, 3)
    weights = ones(Float64, nv(g))
    problem = GenericTensorNetwork(IndependentSet(g, weights); optimizer = TreeSA())
    mwis = solve(problem, SizeMax())[1].n
    @test abs(dynamic_ob_mwis(g, weights, slicer = slicer, reducer = reducer, verbose = 2) - mwis) < 1e-12
end