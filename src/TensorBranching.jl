module TensorBranching

using Clustering, NLsolve
using Reexport

@reexport using BitBasis

export Clause, clause, clauses
export bithclust, clustering
export sbranches, complexity

include("bitstring.jl")
include("clustering.jl")
include("setcovering.jl")
include("branching.jl")

end
