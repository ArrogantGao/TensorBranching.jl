using Test, TensorBranching

@testset "contructing all subcovering from mis truth table" begin
    graph_sat = graph_from_tuples(5, [(1, 2), (2,3), (2,4), (1,3), (3, 4), (4, 5), (2,5)])
    vertices = [1, 2, 3, 4, 5]
    tbl = reduced_alpha_configs(TensorBranching.TensorNetworkSolver(), graph_sat, [1, 4, 5])
    bss = tbl2longlongint(tbl)
    for measurement in [NumOfVertices(), D3Measure()]
        scn = subcovers_naive(5, bss, vertices, graph_sat, measurement)
        sc = subcovers(bss, vertices, graph_sat, measurement)
        @test length(scn) == length(sc)
        for scni in scn
            @test scni ∈ sc
        end
    end
end
