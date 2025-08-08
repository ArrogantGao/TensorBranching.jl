module TensorBranching

using OptimalBranching
using OptimalBranching.OptimalBranchingCore, OptimalBranching.OptimalBranchingMIS

using GenericTensorNetworks, ProblemReductions
using Graphs, TropicalNumbers, OMEinsum, AbstractTrees, TreeWidthSolver, OMEinsumContractionOrders
using JuMP, SCIP
using JLD2, Random, UnicodePlots, CSV, DataFrames
using Base.Threads

# types
export AbstractRegionSelector, ScoreRS
export AbstractSlicer, ContractionTreeSlicer
export AbstractRefiner, TreeSARefiner, ReoptimizeRefiner
export AbstractBrancher, GreedyBrancher, FixedPointBrancher
export AbstractReducer, MISReducer, XiaoReducer, TensorNetworkReducer

export SlicedBranch, CompressedEinsum

# tree decomposition
export decompose, max_bag
export order2eincode, eincode2order
export compress, uncompress

# omeinsum interface
export mis_complexity, auto_slicing, random_ksg, contraction_peak_memory, contraction_all_memory

# dynamic ob
export kernelize, initialize_code, contract_slices
export slice, slice_bfs, slice_bfs_rw, slice_lpscore
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
