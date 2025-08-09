using TensorBranching
using OMEinsum, Graphs, AbstractTrees
using CSV, DataFrames
using GenericTensorNetworks, ProblemReductions
using TropicalGEMM

using Test
using Random
Random.seed!(1234)

using TensorBranching: load_all_finished

@testset "edge case: smaller than sc_target" begin
    g = random_regular_graph(60, 3)
    code = initialize_code(g, TreeSA())
    reducer = XiaoReducer()
    slicer = ContractionTreeSlicer(sc_target = 100, brancher = FixedPointBrancher(), refiner = TreeSARefiner())
    branches = slice(MISProblem(g, UnitWeight(nv(g))), code, 10, slicer, reducer)
    @test length(branches) == 1
end

@testset "bfs rw" begin
    temp_dir = tempdir()
    g = random_regular_graph(100, 3)
    code = initialize_code(g, TreeSA())
    reducer = XiaoReducer()
    slicer = ContractionTreeSlicer(sc_target = 10, brancher = FixedPointBrancher(), refiner = TreeSARefiner())

    for ws in [UnitWeight(nv(g)), Float32.(1.0 .+ rand(nv(g)))]
        slice_bfs_rw(SlicedBranch(g, ws, code, zero(eltype(ws))), slicer, reducer, temp_dir, 1)
        branches = load_all_finished(temp_dir)
        res = maximum(contract_slices(branches, Float32, false))
        @test res ≈ solve(GenericTensorNetwork(IndependentSet(g, ws), code, Dict{Int, Int}()), SizeMax(), T = Float32)[].n
    end
end

@testset "dfs lp" begin
    temp_dir = tempdir()
    g = random_regular_graph(100, 3)
    code = initialize_code(g, TreeSA())
    ws = ones(nv(g))
    reducer = XiaoReducer()
    slicer = ContractionTreeSlicer(sc_target = 10)
    slice_dfs_lp(SlicedBranch(g, ws, code, zero(eltype(ws))), slicer, reducer, temp_dir, Float32, false, 1)

    res = CSV.read(joinpath(temp_dir, "slices.csv"), DataFrame)
    @test maximum(res.solution) ≈ solve(GenericTensorNetwork(IndependentSet(g, ws), code, Dict{Int, Int}()), SizeMax(), T = Float32)[].n
end