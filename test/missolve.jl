using Test, TensorBranching, TensorBranching.Graphs
using TensorBranching: BitStr, StaticElementVector, optimal_branching_dnf
using EliminateGraphs: mis2, EliminateGraph

@testset "branching strategy" begin
    @test BitStr(StaticElementVector(2, [0, 0, 1])) == BitStr{3, Int}(4)

    graph = graph_from_tuples(3, [(1, 2), (2, 3), (3, 1)])
    tbl = reduced_alpha_configs(graph, [1, 2])
    v = [1, 2, 3]
    @test TensorBranching.impl_strategy(graph, v, tbl, NaiveBranching(), NaiveMeasure())[1] == DNF([Clause(bit"111", bit"100")])
    @test TensorBranching.setcover_strategy(tbl, v, graph, 3, NaiveMeasure())[1] == DNF([Clause(bit"111", bit"100")])
end

@testset "optimal_branching_dnf" begin
    petersen = smallgraph(:petersen)
    vertices, openvertices, dnf = optimal_branching_dnf(petersen, SetCoverBranching(), D3Measure(), MinBoundSelector(2), EnvFilter())
    @test dnf isa DNF
    @test length(vertices) == 10
    @test length(openvertices) == 0
    @test count_ones(dnf.clauses[1].val) == 4

    vertices, openvertices, dnf = optimal_branching_dnf(petersen, SetCoverBranching(), D3Measure(), MinBoundSelector(2), EnvFilter())
    @test dnf isa DNF
    @test length(vertices) == 10
    @test length(openvertices) == 0
    @test count_ones(dnf.clauses[1].val) == 4
end

@testset "missolve" begin
    g = smallgraph(:tutte)
    mis = mis2(EliminateGraph(g))
    for strategy in [NaiveBranching(), SetCoverBranching()], vertex_select in [MinBoundSelector(2)], measurement in [NaiveMeasure(), D3Measure()], filter in [EnvFilter(), NoFilter()]
        @test missolve(g, strategy=strategy, vertex_select=vertex_select, measurement=measurement, filter=filter) == mis
        @test missolve(g, show_count=true, strategy=strategy, vertex_select=vertex_select, measurement=measurement, filter=filter)[1] == mis
    end
end