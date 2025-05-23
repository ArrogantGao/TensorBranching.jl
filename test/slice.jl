using TensorBranching
using OMEinsum, Graphs, AbstractTrees
using CSV, DataFrames

using Test
using Random
Random.seed!(1234)

@testset "slice tree" begin
    g = random_regular_graph(200, 3)
    code = initialize_code(g, TreeSA())
    slicer = ContractionTreeSlicer(sc_target = 20, brancher = FixedPointBrancher(), refiner = TreeSARefiner(), search_order = :bfs)
    reducer = XiaoReducer()
    tree = slice_tree(g, code, 10, slicer, reducer)
    
    for node in PostOrderDFS(tree)
        if TensorBranching.isslicingleaf(node)
            @test TensorBranching.sc(node.node) ≤ slicer.sc_target
        end
    end
end

@testset "edge case: smaller than sc_target" begin
    g = random_regular_graph(60, 3)
    code = initialize_code(g, TreeSA())
    reducer = XiaoReducer()
    for strategy in [:dfs, :bfs]
        slicer = ContractionTreeSlicer(sc_target = 100, brancher = FixedPointBrancher(), refiner = TreeSARefiner(), search_order = strategy)
        branches = slice(g, code, 10, slicer, reducer)
        @test length(branches) == 1
    end

    slicer_tree = ContractionTreeSlicer(sc_target = 100, brancher = FixedPointBrancher(), refiner = TreeSARefiner(), search_order = :tree)
    tree = slice(g, code, 10, slicer_tree, reducer)
    @test tree.node.r == 10
    @test isempty(tree.children)

    temp_dir = tempdir()
    slicer_bfs_rw = ContractionTreeSlicer(sc_target = 100, brancher = FixedPointBrancher(), refiner = TreeSARefiner(), search_order = :bfs_rw)
    slice(g, code, 10, slicer_bfs_rw, reducer, dirname = temp_dir)
    df = CSV.read(joinpath(temp_dir, "slices.csv"), DataFrame)
    @test length(df.id) == 1
end