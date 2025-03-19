using TensorBranching
using Test

@testset "tree decomposition" begin
    include("utils.jl")
    include("decompose.jl")
    include("dynamic_ob.jl")
    include("io.jl")
end
