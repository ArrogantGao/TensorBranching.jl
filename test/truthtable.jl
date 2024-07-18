using Test, TensorBranching, TensorBranching.Graphs, GenericTensorNetworks
using TensorBranching: neighbor_cover

@testset "graphs from tuple" begin
    @testset "graph_from_tuples" begin
        g = graph_from_tuples(4, [(1, 2), (2, 3), (3, 4)])
        @test nv(g) == 4
        @test ne(g) == 3
    end
end

@testset "neighbor cover" begin
    g = smallgraph(:petersen)
    @test length(neighbor_cover(g, 1, 0)[1]) == 1
    @test length(neighbor_cover(g, 1, 0)[2]) == 1
    @test length(neighbor_cover(g, 1, 1)[1]) == 4
    @test length(neighbor_cover(g, 1, 1)[2]) == 3
    @test length(neighbor_cover(g, 1, 2)[1]) == 10
    @test length(neighbor_cover(g, 1, 2)[2]) == 0
end

@testset "reduced alpha" begin
    g = graph_from_tuples(3, [(1, 2), (2, 3), (3, 1)])
    @test reduced_alpha(g, [1, 2]) == Tropical.([1 -Inf; -Inf -Inf])

    cfgs = TensorBranching._reduced_alpha_configs(g, [1, 2])
    @test count(!iszero, cfgs) == 1
    @test collect_configs.(cfgs, Ref("abc")) == reshape([["c"], String[], String[], String[]], 2, 2)
    @test collect_configs.(cfgs) == reshape([[BitVector((0, 0, 1))], [], [], []], 2, 2)
    @test BranchingTable(cfgs) == BranchingTable(3, [[StaticElementVector(2, [0, 0, 1])]])
    println(BranchingTable(cfgs))
    @test reduced_alpha_configs(TensorNetworkSolver(), g, [1, 2]) == BranchingTable(cfgs)
end

@testset "satellite" begin
    graph_sat = graph_from_tuples(5, [(1, 2), (2,3), (2,4), (1,3), (3, 4), (4, 5), (2,5)])
    tbl = reduced_alpha_configs(TensorBranching.TensorNetworkSolver(), graph_sat, [1, 4, 5])
    @test tbl == BranchingTable(5, [
        [StaticElementVector(2, [0, 0, 1, 0, 0]), StaticElementVector(2, [0, 1, 0, 0, 0])],
        [StaticElementVector(2, [1, 0, 0, 1, 0])],
        [StaticElementVector(2, [0, 0, 1, 0, 1])]
    ])
    a, b, c, d, e = booleans(5)
    @test !covered_by(tbl, DNF(a ∧ ¬b))
    @test covered_by(tbl, DNF(a ∧ ¬b ∧ d ∧ ¬e, ∧(¬a, ¬b, c, ¬d)))
    @test covered_by(tbl, DNF(a ∧ ¬b ∧ d ∧ ¬e, ∧(¬a, ¬b, c, ¬d)))
    @test !covered_by(tbl, DNF(a ∧ ¬b ∧ d ∧ ¬e, ∧(¬a, ¬b, c, ¬d, e)))
    @test covered_by(tbl, DNF(a ∧ ¬b ∧ d ∧ ¬e, ∧(¬a, ¬b, c, ¬d, e), ∧(¬a, b, ¬c, ¬d, ¬e)))
end