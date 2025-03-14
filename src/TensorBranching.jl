module TensorBranching

using GenericTensorNetworks
using Graphs, TropicalNumbers, OMEinsum, AbstractTrees, TreeWidthSolver

using OptimalBranching
using OptimalBranching.OptimalBranchingCore, OptimalBranching.OptimalBranchingMIS

export ContractionTreeSlicer, MaxIntersectRS
export SlicedBranch


export eincode2order, eincode2graph, decompose, max_bag
export kernelize, slice

export get_subtree_pre, get_subtree_post, list_subtree, most_label_subtree
export remove_tensors, remove_tensors!, tensors_removed, unsafe_flatten, rethermalize, reindex_tree!

include("types.jl")
include("utils.jl")

# dynamic ob
include("dynamic_ob.jl")
include("kernelize.jl")
include("slice.jl")
include("branch.jl")

end
