using TensorBranching
using Graphs
using OptimalBranching.OptimalBranchingMIS
using OptimalBranching.OptimalBranchingMIS.EliminateGraphs
using GenericTensorNetworks, ProblemReductions
using Test
using Random
Random.seed!(123)

@testset "dynamic_ob for mis, treesa refiner" begin
    for reducer in [BasicReducer(), XiaoReducer()]
        for brancher in [GreedyBrancher(), FixedPointBrancher()]
            for refiner in [TreeSARefiner()]
                for selector in [ScoreRS(loss = :num_uniques), ScoreRS(loss = :sc_score), ScoreRS(loss = :bag_score)]
                    for ws in [UnitWeight(100), Float32.(1.0 .+ rand(100))]
                        slicer = ContractionTreeSlicer(sc_target = 10, brancher = brancher, refiner = refiner, region_selector = selector)
                        fw = ws isa UnitWeight ? "unweighted" : "weighted"
                        @info "reducer = $reducer, brancher = $brancher, refiner = $refiner, region_selector = $selector, weights = $fw"
                        g = random_regular_graph(100, 3)
                        @test dynamic_ob_mis(g, ws, slicer = slicer, reducer = reducer) ≈ solve(GenericTensorNetwork(IndependentSet(g, ws), optimizer = TreeSA()), SizeMax(), T = Float32)[].n
                    end
                end
            end
        end
    end
end

@testset "dynamic_ob for mis, tensor network reducer" begin
    # tensor network reducer
    reducer = TensorNetworkReducer()
    brancher = GreedyBrancher()
    refiner = TreeSARefiner()
    slicer = ContractionTreeSlicer(sc_target = 10, brancher = brancher, refiner = refiner)
    @info "reducer = $reducer, brancher = $brancher, refiner = $refiner"

    g = random_regular_graph(100, 3)
    for ws in [UnitWeight(100), Float32.(1.0 .+ rand(100))]
        @test dynamic_ob_mis(g, ws, slicer = slicer, reducer = reducer) ≈ solve(GenericTensorNetwork(IndependentSet(g, ws), optimizer = TreeSA()), SizeMax(), T = Float32)[].n
    end
end

@testset "dynamic_ob for mis, reoptimize refiner" begin
    reducer = XiaoReducer()
    brancher = FixedPointBrancher()
    refiner = ReoptimizeRefiner()
    slicer = ContractionTreeSlicer(sc_target = 10, brancher = brancher, refiner = refiner)
    @info "reducer = $reducer, brancher = $brancher, refiner = $refiner"

    g = random_regular_graph(100, 3)
    for ws in [UnitWeight(100), Float32.(1.0 .+ rand(100))]
        @test dynamic_ob_mis(g, ws, slicer = slicer, reducer = reducer, verbose = 1) ≈ solve(GenericTensorNetwork(IndependentSet(g, ws), optimizer = TreeSA()), SizeMax(), T = Float32)[].n
    end
end