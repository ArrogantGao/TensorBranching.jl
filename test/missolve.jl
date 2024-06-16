using Test, TensorBranching, TensorBranching.Graphs
using TensorBranching: BitStr, StaticElementVector, optimal_branching_dnf
using EliminateGraphs: mis1, EliminateGraph

@testset "branching strategy" begin
    @test BitStr(StaticElementVector(2, [0, 0, 1])) == BitStr{3, Int}(4)

    graph = graph_from_tuples(3, [(1, 2), (2, 3), (3, 1)])
    tbl = reduced_alpha_configs(graph, [1, 2])
    @test TensorBranching.naive_strategy(tbl) == DNF([Clause(bit"111", bit"100")])
    @test TensorBranching.setcover_strategy(tbl) == DNF([Clause(bit"111", bit"100")])
end

@testset "optimal_branching_dnf" begin
    petersen = smallgraph(:petersen)
    vertices, openvertices, dnf = optimal_branching_dnf(petersen)
    @test dnf isa DNF
    @test length(vertices) == 10
    @test length(openvertices) == 0
    @test count_ones(dnf.clauses[1].val) == 4

    vertices, openvertices, dnf = optimal_branching_dnf(petersen, strategy=SetCoverBranching())
    @test dnf isa DNF
    @test length(vertices) == 10
    @test length(openvertices) == 0
    @test count_ones(dnf.clauses[1].val) == 4
end

@testset "missolve" begin
    g = smallgraph(:tutte)
    @test missolve(g) == mis1(EliminateGraph(g))
    @test missolve(g, strategy=SetCoverBranching()) == mis1(EliminateGraph(g))
end