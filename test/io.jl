using TensorBranching
using Graphs, OMEinsum, OptimalBranching.OptimalBranchingMIS, TropicalNumbers
using GenericTensorNetworks, ProblemReductions
using CSV, DataFrames
using Test

using TensorBranching: save_weights, load_weights, save_finished, load_finished, save_unfinished, load_unfinished

@testset "io" begin
    temp_dir = tempdir()

    # create some slices and tensors
    g = random_regular_graph(100, 3)
    code = initialize_code(g, TreeSA())
    ws = UnitWeight(100)
    r = 10
    branch = SlicedBranch(MISProblem(g, ws), code, r)

    save_finished(temp_dir, branch, 1)
    g_loaded, code_loaded, ws_loaded = load_finished(temp_dir, 1)
    @test g_loaded == g
    @test code_loaded == code
    @test ws_loaded == ws

    save_unfinished(temp_dir, branch, 1)
    # load the slices and tensors
    branch_loaded = load_unfinished(temp_dir, 1)
    @test branch_loaded.p == branch.p
    @test uncompress(branch_loaded.code) == uncompress(branch.code)
    @test branch_loaded.r == branch.r

    for type in [Float64, Float32, Int64, Int32]
        ws = rand(type, 100)
        save_weights(joinpath(temp_dir, "weights.txt"), ws)
        ws_loaded = load_weights(joinpath(temp_dir, "weights.txt"))
        @test ws_loaded == ws
    end
end
