using TensorBranching
using Test

@testset "TensorBranching.jl" begin
    include("bitstrings.jl")
    include("clustering.jl")
end

@testset "setcovering" begin
    include("setcovering.jl")
end

@testset "truth table" begin
    include("truthtable.jl")
end

@testset "missolve" begin
    include("missolve.jl")
end