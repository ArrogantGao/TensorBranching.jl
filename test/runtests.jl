using TensorBranching
using Test

@testset "TensorBranching.jl" begin
    include("bitstrings.jl")
end

@testset "setcovering" begin
    include("clauses.jl")
    include("setcovering.jl")
end

@testset "truth table" begin
    include("truthtable.jl")
end

@testset "missolve" begin
    include("missolve.jl")
end

@testset "artifacts" begin
    include("artifacts.jl")
end