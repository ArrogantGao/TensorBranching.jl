using Test, TensorBranching
using Random
Random.seed!(1234)

@testset "setcover by JuMP" begin
    bs = [bit"0000", bit"1000", bit"1010", bit"1110", bit"0001", bit"1001", bit"0101", bit"1101", bit"1011", bit"1111"]
    bbs = [[bsi] for bsi in bs]
    g = smallgraph(:petersen)
    v = [1:4...]
    scs = subcovers(bbs, v, g, NumOfVertices())
    cov_lp, cx_lp = cover(copy(scs), 3, LPSetCoverSolver())
    cov_ip, cx_ip = cover(copy(scs), 3, IPSetCoverSolver())
    picked_lp = Set{Int}()
    picked_ip = Set{Int}()
    for covi in cov_lp
        picked_lp = union(picked_lp, covi.ids)
    end
    for covi in cov_ip
        picked_ip = union(picked_ip, covi.ids)
    end
    @test length(picked_lp) == length(picked_ip) == length(bs)
end