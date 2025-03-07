module TensorBranching

using Graphs, AbstractTrees
using TreeWidthSolver
using OMEinsum, GenericTensorNetworks

using OMEinsum.AbstractTrees

using OptimalBranchingCore, OptimalBranchingMIS
using OptimalBranchingMIS.KaHyPar, OptimalBranchingMIS.SparseArrays

export decompose

export kernelize

include("types.jl")
include("utils.jl")

# function about tree decomposition and reformulation
include("tree_decompose.jl")
include("tree_select.jl")
include("tree_reform.jl")

# dynamic ob
include("kernelize.jl")
include("dynamic_ob.jl")

end
