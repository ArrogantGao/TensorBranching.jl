using Test, TensorBranching, TensorBranching.Graphs
using TensorBranching: BitStr, StaticElementVector, optimal_branches
using EliminateGraphs: mis2, EliminateGraph

@testset "branching strategy" begin
    @test BitStr(StaticElementVector(2, [0, 0, 1])) == BitStr{3, Int}(4)

    graph = graph_from_tuples(3, [(1, 2), (2, 3), (3, 1)])
    tbl = reduced_alpha_configs(graph, [1, 2])
    v = [1, 2, 3]
    @test TensorBranching.impl_strategy(graph, v, tbl, NaiveBranching(), NumOfVertices()) == Branches{Int64}(Branch{Int64}[Branch{Int64}([1, 2, 3], 3, 1)])
    @test TensorBranching.setcover_strategy(tbl, v, graph, 3, NumOfVertices()) == Branches{Int64}(Branch{Int64}[Branch{Int64}([1, 2, 3], 3, 1)])
end

@testset "optimal_branches" begin
    petersen = smallgraph(:petersen)
    for strategy in [NaiveBranching(), SetCoverBranching()], vertex_select in [ManualSelector([1, 2, 3, 4]), MinBoundarySelector(2)], measurement in [NumOfVertices(), NumOfDegree()], filter in [EnvFilter(), NoFilter()]
        branches = optimal_branches(petersen, strategy, measurement, vertex_select, filter)
        @test branches isa Branches
    end
end

@testset "missolve" begin
    g = smallgraph(:tutte)
    mis = mis2(EliminateGraph(g))
    for strategy in [NaiveBranching(), SetCoverBranching()], vertex_select in [MinBoundarySelector(2)], measurement in [NumOfVertices(), NumOfDegree()], filter in [EnvFilter(), NoFilter()]
        @test missolve(g, strategy=strategy, vertex_select=vertex_select, measurement=measurement, filter=filter) == mis
        @test missolve(g, show_count=true, strategy=strategy, vertex_select=vertex_select, measurement=measurement, filter=filter)[1] == mis
    end
end