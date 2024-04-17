module TensorBranching

using Clustering, NLsolve
using Reexport
using GenericTensorNetworks, GenericTensorNetworks.Graphs

@reexport using BitBasis

export Clause, clause, clauses
export bithclust, clustering
export sbranches, complexity
export graph_from_tuples, reduced_alpha, reduced_alpha_configs, collect_configs, BranchingTable, DNF, booleans, covered_by

include("bitstring.jl")
include("clustering.jl")
include("setcovering.jl")
include("branching.jl")
include("truthtable.jl")

end
