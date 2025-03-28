module TensorBranching

using OptimalBranching
using OptimalBranching.OptimalBranchingCore, OptimalBranching.OptimalBranchingMIS

using GenericTensorNetworks
using Graphs, TropicalNumbers, OMEinsum, AbstractTrees, TreeWidthSolver
using Suppressor, JLD2, Random, ProgressBars

# types
export AbstractRegionSelector, MaxIntersectRS, ScScoreRS
export AbstractSlicer, ContractionTreeSlicer
export AbstractRefiner, TreeSARefiner, ReoptimizeRefiner
export AbstractBrancher, GreedyBrancher, FixedPointBrancher
export AbstractReducer, MISReducer, XiaoReducer, TensorNetworkReducer

export SlicedBranch

# tree decomposition
export decompose, max_bag
export order2eincode, eincode2order

# omeinsum interface
export mis_complexity, auto_slicing, random_ksg

# dynamic ob
export kernelize, initialize_code, slice, contract_slices
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
include("refine.jl")

# io functions
include("io.jl")

end
