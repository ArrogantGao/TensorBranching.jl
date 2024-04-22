using Test, TensorBranching

@testset "gathering algorithm" begin
    n = 30
    m = 10
    bs = [rand(BitStr{m, Int64}) for i in 1:n]
    ca = Clause.(bs)
    cca = copy(ca)
    for _=1:n - 1
        cca = gather(cca)
        @test foldl(.|, [covered_by(bs, ccai) for ccai in cca]; init = [false for i in 1:n]) == ones(n)
    end
end