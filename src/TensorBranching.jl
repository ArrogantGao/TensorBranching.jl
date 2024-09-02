module TensorBranching

using NLsolve, JuMP, HiGHS, SCIP
using GenericTensorNetworks
import EliminateGraphs

using Base.Threads
using Reexport
@reexport using BitBasis, GenericTensorNetworks.Graphs

export Clause, clause, clauses, SubCover
export subcovers, subcovers_naive, all_clauses_naive
export LP_setcover, random_pick, cover
export sbranches, complexity
export graph_from_tuples, reduced_alpha, reduced_alpha_configs, collect_configs, BranchingTable, DNF, booleans, covered_by

export SolverConfig
export AbstractBranching, NaiveBranching, SetCoverBranching
export AbstractMeasure, NumOfVertices, D3Measure
export AbstractVertexSelector, MinBoundarySelector
export AbstractTruthFilter, NoFilter, EnvFilter
export AbstractMISSolver, TensorNetworkSolver
export AbstractSetCoverSolver, IPSetCoverSolver, LPSetCoverSolver

export Branch, Branches, effective_γ, optimal_branches

export solve_mis, count_mis, CountingMIS, branching_tree

export graph_from_artifact

include("artifacts.jl")
include("counting_mis.jl")
include("bitstring.jl")
include("truthtable.jl")
include("strategy.jl")
include("graphs.jl")
include("clauses.jl")
include("setcovering.jl")
include("solver.jl")
include("branching_tree.jl")

end
