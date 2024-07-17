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

@testset "gather2" begin
    INT = LongLongUInt{1}
    mask = bmask(INT, 1:5)
    v1 = LongLongUInt{1}((0b00010,))
    v2 = LongLongUInt{1}((0b01001,))
    c1 = Clause(mask, v1)
    c2 = Clause(mask, v2)
    c3 = TensorBranching.gather2(5, c1, c2)
    @test c3 == Clause(LongLongUInt{1}((0b10100,)), LongLongUInt{1}((0b0,)))
end