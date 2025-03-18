using TensorBranching
using TropicalNumbers
using Graphs
using OMEinsum
using Test

@testset "io" begin
    temp_dir = tempdir()
    filename = joinpath(temp_dir, "test.jld2")

    # create some slices and tensors
    g = random_regular_graph(30, 3)
    tensors = initialize_tensors(g, false, TropicalF32)
    code = initialize_code(g, TreeSA())
    r = 3
    slices = slice(SlicedBranch(g, code, r), ContractionTreeSlicer(), 0)
    saveslices(slices, tensors, r, filename)

    # load the slices and tensors
    slices_loaded, tensors_loaded, r_loaded = loadslices(filename)
    @test length(slices_loaded) == length(slices)
    @test tensors_loaded == tensors
    @test r_loaded == r
end
