module TensorBranching

using Clustering, NLsolve
using Reexport
using GenericTensorNetworks, GenericTensorNetworks.Graphs
import EliminateGraphs

@reexport using BitBasis

export Clause, clause, clauses
export bithclust, clustering, gather, gather2
export sbranches, complexity
export graph_from_tuples, reduced_alpha, reduced_alpha_configs, collect_configs, BranchingTable, DNF, booleans, covered_by
export missolve

include("bitstring.jl")
include("clustering.jl")
include("setcovering.jl")
include("branching.jl")
include("truthtable.jl")
include("solver.jl")

end
