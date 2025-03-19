module TensorBranching

using OptimalBranching
using OptimalBranching.OptimalBranchingCore, OptimalBranching.OptimalBranchingMIS

using OptimalBranching.OptimalBranchingMIS.GenericTensorNetworks
using Graphs, TropicalNumbers, OMEinsum, AbstractTrees, TreeWidthSolver

using JLD2

# types
export AbstractRegionSelector, MaxIntersectRS
export AbstractSlicer, ContractionTreeSlicer
export SlicedBranch

# tree decomposition
export decompose

# dynamic ob
export kernelize, initialize_code, initialize_tensors, slice, contract
export dynamic_ob_mis

#io
export saveslices, loadslices

include("types.jl")
include("utils.jl")
include("decompose.jl")

# dynamic ob
include("dynamic_ob.jl")
include("kernelize.jl")
include("slice.jl")
include("branch.jl")

# io functions
include("io.jl")

end
