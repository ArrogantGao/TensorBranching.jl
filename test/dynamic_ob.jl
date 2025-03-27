using TensorBranching
using Graphs
using OptimalBranching.OptimalBranchingMIS
using OptimalBranching.OptimalBranchingMIS.EliminateGraphs
using GenericTensorNetworks
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
    reducer = GreedyBrancher()
    brancher = FixedPointBrancher()
    refiner = ReoptimizeRefiner()
    search_order = :bfs
    slicer = ContractionTreeSlicer(sc_target = 10, brancher = brancher, refiner = refiner, search_order = search_order)
    @info "reducer = $reducer, brancher = $brancher, refiner = $refiner, search_order = $search_order"
    g = random_regular_graph(100, 3)
    @test dynamic_ob_mis(g, slicer = slicer, reducer = reducer, verbose = 1) == mis2(EliminateGraph(g))
end

#It takes hours to run this testset.Orz
@testset "dynamic_ob for mis, comparison of original region selector and sc_score region selector" begin
    n_sqrt = 60
    filling = 0.8
    for seed in 1:5
        Random.seed!(seed)
        reducer = TensorNetworkReducer()
        brancher = FixedPointBrancher()
        refiner = TreeSARefiner()
        sc_target = 25
        search_order = :bfs
        slicer_maxintersect = ContractionTreeSlicer(sc_target = sc_target, brancher = brancher, refiner = refiner, search_order = search_order, region_selector = MaxIntersectRS())
        slicer_scscore = ContractionTreeSlicer(sc_target = sc_target, brancher = brancher, refiner = refiner, search_order = search_order, region_selector = ScScoreRS())

        graph = GenericTensorNetworks.random_diagonal_coupled_graph(n_sqrt, n_sqrt, filling)
        graph = SimpleGraph(graph)
        result_maxintersect = dynamic_ob_mis(graph, slicer = slicer_maxintersect, reducer = reducer, verbose = 1) 
        result_scscore = dynamic_ob_mis(graph, slicer = slicer_scscore, reducer = reducer, verbose = 1) 
        @test result_maxintersect == result_scscore
    end
end