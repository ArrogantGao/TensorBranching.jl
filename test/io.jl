using TensorBranching
using Graphs, OMEinsum, OptimalBranching.OptimalBranchingMIS, TropicalNumbers
using CSV, DataFrames
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

    @test isnothing(slice(g, code, r, ContractionTreeSlicer(sc_target = 3), XiaoReducer(); dirname = temp_dir))
    df = CSV.read(joinpath(temp_dir, "slices.csv"), DataFrame)
    @test df isa DataFrame
end
