using Test, TensorBranching

@testset "setcover by JuMP" begin
    bs = [bit"0000", bit"1000", bit"1010", bit"1110", bit"0001", bit"1001", bit"0101", bit"1101", bit"1011", bit"1111"]
    scs = subcovers(bs)
    cov, cx = cover(scs, length(bs))
    @test isapprox(cx, 1.618, atol = 1e-3)
end