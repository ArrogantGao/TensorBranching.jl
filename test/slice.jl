using TensorBranching
using OMEinsum, Graphs, AbstractTrees

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
