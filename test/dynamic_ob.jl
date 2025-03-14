using TensorBranching
using Graphs
using OptimalBranching.OptimalBranchingMIS.EliminateGraphs

using Test
using Random
Random.seed!(1234)

@testset "dynamic_ob for mis" begin
    for i in 1:10
        g = random_regular_graph(50, 3)
        slicer = ContractionTreeSlicer(sc_target = 5)

        @test dynamic_ob_mis(g, slicer = slicer) == mis2(EliminateGraph(g))
    end
end


@testset "dynamic_ob large graph" begin
    g = random_regular_graph(100, 3)
    slicer = ContractionTreeSlicer(sc_target = 10)

    @test dynamic_ob_mis(g, slicer = slicer) == mis2(EliminateGraph(g))
end
