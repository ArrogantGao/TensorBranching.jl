module TensorBranching

using GenericTensorNetworks
using GenericTensorNetworks.Graphs
using GenericTensorNetworks.TropicalNumbers
using GenericTensorNetworks.OMEinsum, GenericTensorNetworks.OMEinsum.AbstractTrees
using GenericTensorNetworks.OMEinsum.OMEinsumContractionOrders.TreeWidthSolver

using OptimalBranching
using OptimalBranching.OptimalBranchingCore, OptimalBranching.OptimalBranchingMIS
using OptimalBranching.OptimalBranchingMIS.KaHyPar, OptimalBranching.OptimalBranchingMIS.SparseArrays

export decompose
export kernelize

include("types.jl")
include("tree_utils.jl")

# dynamic ob
include("dynamic_ob.jl")
include("kernelize.jl")
include("slice.jl")
include("branch.jl")

end
