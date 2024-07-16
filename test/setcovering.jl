using Test, TensorBranching

@testset "setcover by JuMP" begin
    bs = [bit"0000", bit"1000", bit"1010", bit"1110", bit"0001", bit"1001", bit"0101", bit"1101", bit"1011", bit"1111"]
    bbs = [[bsi] for bsi in bs]
    g = random_regular_graph(10, 3)
    v = [1:4...]
    scs = subcovers(bbs, v, g, NaiveMeasure())
    cov, cx = cover(scs, 4)
    picked = Set{Int}()
    for covi in cov
        picked = union(picked, covi.ids)
    end
    @test length(picked) == length(bs)
end