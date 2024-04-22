using Test, TensorBranching

@testset "contructing all subcovers" begin
    bs = unique([rand(BitStr{10, Int64}) for i in 1:20])
    scn = subcovers_naive(bs)
    sc = subcovers(bs)
    @test length(scn) == length(sc)
    for scni in scn
        @test scni ∈ sc
    end
end

@testset "contructing all subcovering from mis truth table" begin
    graph_sat = graph_from_tuples(5, [(1, 2), (2,3), (2,4), (1,3), (3, 4), (4, 5), (2,5)])
    tbl = reduced_alpha_configs(graph_sat, [1, 4, 5])
    bss = Tbl2BitStrs(tbl)
    scn = subcovers_naive(bss)
    sc = subcovers(bss)
    @test length(scn) == length(sc)
    for scni in scn
        @test scni ∈ sc
    end
end