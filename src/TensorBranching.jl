module TensorBranching

using Graphs, AbstractTrees
using TreeWidthSolver
using OMEinsum, GenericTensorNetworks

using OMEinsum.AbstractTrees

using OptimalBranchingCore, OptimalBranchingMIS
using BitBasis
using OptimalBranchingMIS.KaHyPar, OptimalBranchingMIS.SparseArrays

export MaskedGraph

export decompose
export kernelize

include("types.jl")
include("tree_utils.jl")
include("graph_utils.jl")

# dynamic ob
include("dynamic_ob.jl")
include("kernelize.jl")
include("slice.jl")

end
