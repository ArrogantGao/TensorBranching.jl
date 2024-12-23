module TensorBranching

using Graphs, AbstractTrees
using TreeWidthSolver
using OMEinsum, GenericTensorNetworks

export decompose

include("tree_decompose.jl")

end
