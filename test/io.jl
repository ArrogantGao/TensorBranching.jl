using TensorBranching
using Graphs, OMEinsum, OptimalBranching.OptimalBranchingMIS, TropicalNumbers
using Test

@testset "io" begin
    temp_dir = tempdir()
    filename = joinpath(temp_dir, "test.jld2")

    # create some slices and tensors
    g = random_regular_graph(30, 3)
    code = initialize_code(g, TreeSA())
    r = 3
    slices = slice(g, code, r, ContractionTreeSlicer(sc_target = 3), XiaoReducer())
    saveslices(filename, slices)

    # load the slices and tensors
    slices_loaded = loadslices(filename)
    @test length(slices_loaded) == length(slices)
    @test contract_slices(slices_loaded, TropicalF32, false) == contract_slices(slices, TropicalF32, false)
end
