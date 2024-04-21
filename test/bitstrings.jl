using TensorBranching, Test

@testset "clause and cover" begin
    c1 = Clause(bit"1110", bit"0000")
    c2 = Clause(bit"1110", bit"0001")
    c3 = Clause(bit"1110", bit"0010")
    c4 = Clause(bit"1100", bit"0001")
    @test c1 == c2
    @test c1 !== c3
    @test c1 !== c4
end