module TensorBranching

using GenericTensorNetworks
using Graphs, TropicalNumbers, OMEinsum, AbstractTrees, TreeWidthSolver

using OptimalBranching
using OptimalBranching.OptimalBranchingCore, OptimalBranching.OptimalBranchingMIS

# types
export AbstractRegionSelector, MaxIntersectRS
export AbstractSlicer, ContractionTreeSlicer
export SlicedBranch

# tree decomposition
export decompose

# dynamic ob
export kernelize, slice
export dynamic_ob_mis

include("types.jl")
include("utils.jl")

# dynamic ob
include("dynamic_ob.jl")
include("kernelize.jl")
include("slice.jl")
include("branch.jl")

end
