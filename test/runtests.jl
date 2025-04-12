using TensorBranching
using Test

@testset "utils" begin
    include("utils.jl")
end

@testset "tree decomposition" begin
    include("decompose.jl")
end

@testset "kernelize" begin
    include("kernelize.jl")
end

@testset "slice" begin
    include("slice.jl")
end

@testset "dynamic ob" begin
    include("dynamic_ob.jl")
end

@testset "io" begin
    include("io.jl")
end
