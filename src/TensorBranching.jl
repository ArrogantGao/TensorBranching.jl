module TensorBranching

using Clustering, NLsolve, JuMP, HiGHS
using GenericTensorNetworks, GenericTensorNetworks.Graphs
import EliminateGraphs

using Reexport
@reexport using BitBasis

export Clause, clause, clauses, SubCover, Tbl2BitStrs
export subcovers, subcovers_naive, all_clauses_naive
export LP_setcover, random_pick, cover
export sbranches, complexity
export graph_from_tuples, reduced_alpha, reduced_alpha_configs, collect_configs, BranchingTable, DNF, booleans, covered_by
export BranchingStrategy, NaiveBranching, SetCoverBranching
export missolve


include("bitstring.jl")
include("truthtable.jl")
include("clauses.jl")
include("setcovering.jl")
include("branching.jl")
include("solver.jl")

end
