using Test, TensorBranching, TensorBranching.Graphs
using TensorBranching: LongLongUInt, StaticElementVector, optimal_branches
using EliminateGraphs: mis2, EliminateGraph

@testset "branching strategy" begin
    @test TensorBranching._vec2int(LongLongUInt, StaticElementVector(2, [0, 0, 1])) == LongLongUInt{1}(4)

    graph = graph_from_tuples(3, [(1, 2), (2, 3), (3, 1)])
    tbl = reduced_alpha_configs(TensorNetworkSolver(), graph, [1, 2])
    v = [1, 2, 3]
    @test TensorBranching.impl_strategy(graph, v, tbl, NaiveBranching(), NumOfVertices()) == Branches(Branch[Branch([1, 2, 3], 1)])
    @test TensorBranching.setcover_strategy(tbl, v, graph, 3, NumOfVertices(), IPSetCoverSolver(), false) == Branches(Branch[Branch([1, 2, 3], 1)])
    @test TensorBranching.setcover_strategy(tbl, v, graph, 3, NumOfVertices(), LPSetCoverSolver(), false) == Branches(Branch[Branch([1, 2, 3], 1)])
    @test TensorBranching.impl_strategy(graph, v, tbl, SetCoverBranching(), NumOfVertices()) == Branches(Branch[Branch([1, 2, 3], 1)])
end

@testset "optimal_branches" begin
    petersen = smallgraph(:petersen)
    for strategy in [NaiveBranching(), SetCoverBranching()], vertex_select in [MinBoundarySelector(2)], measure in [NumOfVertices(), D3Measure()], table_filter in [EnvFilter(), NoFilter()]
        vertices = TensorBranching.select_vertex(petersen, vertex_select)
        branches = optimal_branches(petersen, vertices, strategy; measure, table_filter)
        @test branches isa Branches
    end
end

@testset "optimal branches" begin
    g = smallgraph(:tutte)
    vertices = [1] ∪ neighbors(g, 1)
    cfg = SolverConfig(; branching_strategy=SetCoverBranching())
    branch = optimal_branches(g, vertices, cfg.branching_strategy; cfg.measure, cfg.table_filter)
    @test isapprox(TensorBranching.effective_γ(branch, g, cfg.measure), 1.2405351053760838)
end

@testset "missolve" begin
    g = smallgraph(:tutte)
    mis = mis2(EliminateGraph(g))
    for branching_strategy in [NaiveBranching(), SetCoverBranching()], vertex_selector in [MinBoundarySelector(2)], measure in [NumOfVertices(), D3Measure()], table_filter in [EnvFilter(), NoFilter()]
        @info "Testing missolve with: $branching_strategy, $vertex_selector, $measure, $table_filter"
        cfg = SolverConfig(; branching_strategy, vertex_selector, measure, table_filter)
        @test solve_mis(g, cfg) == mis
        @test count_mis(g, cfg)[1] == mis
    end
end
