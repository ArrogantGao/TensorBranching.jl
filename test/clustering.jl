using Test, TensorBranching

@testset "gathering algorithm" begin
    N = 30
    bs = [rand(BitStr{10, Int64}) for i in 1:30]
    ca = Clause.(bs)
    cca = copy(ca)
    for _=1:N - 1
        cca = gather(cca)
        @test foldl(.|, [covered_by(bs, ccai) for ccai in cca]; init = [false for i in 1:N]) == ones(N)
    end
end