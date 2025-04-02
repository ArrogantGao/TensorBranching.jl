using TensorBranching
using TensorBranching.GenericTensorNetworks
using Graphs, OMEinsum, TropicalNumbers

using TensorBranching: eincode2order, eincode2graph, order2eincode, rethermalize, update_order
using OptimalBranching.OptimalBranchingMIS

using Test
using Random

Random.seed!(1234)

@testset "converting graph and code" begin
    for n in 30:10:100
        g = random_regular_graph(n, 3)
        code = GenericTensorNetwork(IndependentSet(g)).code

        @test eincode2graph(code)[1] == g

        eo = eincode2order(code)
        code_new = order2eincode(g, eo)

        sc_target = Int(contraction_complexity(code, uniformsize(code, 2)).sc)
        re_code_new = rethermalize(code_new, uniformsize(code_new, 2), 100.0:100.0, 1, 10, sc_target)

        for T in (Float64, Tropical{Float64})
            tensors = GenericTensorNetworks.generate_tensors(T(1.0), IndependentSet(g))
            @test code(tensors...) ≈ code_new(tensors...) ≈ re_code_new(tensors...)
            # @show contraction_complexity(code, uniformsize(code, 2)).sc, contraction_complexity(code_new, uniformsize(code_new, 2)).sc, contraction_complexity(re_code_new, uniformsize(re_code_new, 2)).sc
        end
    end
end

@testset "reconstructing code" begin
    for n in 30:10:100
        for _ in 1:5
            g = random_regular_graph(n, 3)
            code = GenericTensorNetwork(IndependentSet(g)).code
            removed_vertices = unique!(rand(1:n, 10))

            g_new, vmap = induced_subgraph(g, setdiff(1:n, removed_vertices))

            for _ in 1:10
                i = rand(1:nv(g_new))
                j = rand(1:nv(g_new))

                add_edge!(g_new, i, j)
            end

            code = GenericTensorNetwork(IndependentSet(g)).code
            eo = eincode2order(code)
            eo_new = update_order(eo, vmap)
            code_new = order2eincode(g_new, eo_new)

            net_direct = GenericTensorNetwork(IndependentSet(g_new))
            net_new = GenericTensorNetwork(IndependentSet(g_new), code_new, Dict{Int, Int}())
            @test solve(net_direct, SizeMax()) ≈ solve(net_new, SizeMax())
        end
    end
end

@testset "corner case: disconnected graph have 1 vertex" begin
    g = random_regular_graph(30, 3)
    removed_vertices = neighbors(g, 1)
    g_new, vmap = induced_subgraph(g, setdiff(1:nv(g), removed_vertices))
    code = GenericTensorNetwork(IndependentSet(g)).code
    eo = eincode2order(code)
    eo_new = update_order(eo, vmap)
    code_new = order2eincode(g_new, eo_new)
    @test !is_connected(g_new)
    @test solve(GenericTensorNetwork(IndependentSet(g_new)), SizeMax()) ≈ solve(GenericTensorNetwork(IndependentSet(g_new), code_new, Dict{Int, Int}()), SizeMax())
end

@testset "corner case: disconnected graph have 2 vertices" begin
    g = random_regular_graph(30, 3)
    vs = [1, neighbors(g, 1)[1]]
    removed_vertices = OptimalBranchingMIS.open_neighbors(g, vs)
    g_new, vmap = induced_subgraph(g, setdiff(1:nv(g), removed_vertices))
    code = GenericTensorNetwork(IndependentSet(g)).code
    eo = eincode2order(code)
    eo_new = update_order(eo, vmap)
    code_new = order2eincode(g_new, eo_new)
    @test !is_connected(g_new)
    @test solve(GenericTensorNetwork(IndependentSet(g_new)), SizeMax()) ≈ solve(GenericTensorNetwork(IndependentSet(g_new), code_new, Dict{Int, Int}()), SizeMax())
end

@testset "compressing and uncompressing code" begin
    for n in 30:10:100
        g = random_regular_graph(n, 3)
        code = GenericTensorNetwork(IndependentSet(g)).code
        compressed_code = compress(code)
        uncompressed_code = uncompress(compressed_code)

        cc = mis_complexity(code)
        ccc = mis_complexity(uncompressed_code)

        @test cc.sc ≈ ccc.sc
        @test cc.tc ≈ ccc.tc
        @test cc.rwc ≈ ccc.rwc
    end
end