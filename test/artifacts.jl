@testset "graph_from_artifact" begin
    g = graph_from_artifact(1)
    @test nv(g) == 6160 
    @test ne(g) == 40207
end