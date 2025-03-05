module TensorBranching

using Graphs, AbstractTrees
using TreeWidthSolver
using OMEinsum, GenericTensorNetworks

using OMEinsum.AbstractTrees

export decompose

include("tree_decompose.jl")
include("tree_select.jl")
include("tree_reform.jl")

end
