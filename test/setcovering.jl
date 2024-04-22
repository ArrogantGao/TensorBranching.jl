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