using Test, TensorBranching

@testset "contructing all subcovers" begin
    bs = unique([rand(BitStr{10, Int64}) for i in 1:20])
    g = random_regular_graph(10, 3)
    v = [1:10...]
    scn = subcovers_naive(bs, v, g)
    sc = subcovers(bs, v, g)
    @test length(scn) == length(sc)
    for scni in scn
        @test scni ∈ sc
    end
end

@testset "contructing all subcovering from mis truth table" begin
    graph_sat = graph_from_tuples(5, [(1, 2), (2,3), (2,4), (1,3), (3, 4), (4, 5), (2,5)])
    vertices = [1, 2, 3, 4, 5]
    tbl = reduced_alpha_configs(graph_sat, [1, 4, 5])
    bss = Tbl2BitStrs(tbl)
    scn = subcovers_naive(bss, vertices, graph_sat)
    sc = subcovers(bss, vertices, graph_sat)
    @test length(scn) == length(sc)
    for scni in scn
        @test scni ∈ sc
    end
end