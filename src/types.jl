abstract type AbstractRegionSelector end

# the maximum intersection region selector try to find the region with the maximum intersection with the largest tensors, where large means larger that the threshold
@kwdef struct ScoreRS <: AbstractRegionSelector
    n_max::Int = 20
    strategy::Symbol = :neighbor # :mincut or :neighbor
    loss::Symbol = :sc_score # num_uniques or sc_score
end

abstract type AbstractBrancher end

@kwdef struct GreedyBrancher <: AbstractBrancher
    loss::Symbol = :sc_score
    setcover_solver::AbstractSetCoverSolver = IPSolver()
end

@kwdef struct FixedPointBrancher <: AbstractBrancher
    measure::Symbol = :sc_measure
    setcover_solver::AbstractSetCoverSolver = IPSolver()
end

abstract type AbstractRefiner end

@kwdef struct TreeSARefiner{IT} <: AbstractRefiner
    βs::IT = 1.0:1.0:15.0 # range of βs for the rethermalization
    ntrials::Int = 3
    niters::Int = 30
    max_rounds::Int = 2
    reoptimize::Bool = true # setting this to true will reoptimize the code after each round of rethermalization if the resulting sc is larger than sc0
end

@kwdef struct ReoptimizeRefiner <: AbstractRefiner
    optimizer::CodeOptimizer = TreeSA()
end

abstract type AbstractSlicer end

@kwdef struct ContractionTreeSlicer <: AbstractSlicer
    sc_target::Int = 30
    search_order::Symbol = :bfs # :bfs or :dfs
    region_selector::AbstractRegionSelector = ScoreRS() # select the region to branch, what else methods?
    table_solver::AbstractTableSolver = TensorNetworkSolver()
    brancher::AbstractBrancher = GreedyBrancher()
    refiner::AbstractRefiner = TreeSARefiner()
end

struct SlicedBranch{T}
    g::SimpleGraph{T}
    code::Union{DynamicNestedEinsum{T}, Nothing}
    r::Int
end
function Base.show(io::IO, branch::SlicedBranch{T}) where T
    print(io, "SlicedBranch{$T}: ")
    print(io, "graph: {$(nv(branch.g)), $(ne(branch.g))} simple graph; ")
    cc = complexity(branch)
    print(io, "code complexity: sc: $(cc.sc), tc: $(cc.tc)")
    print(io, "; fixed ones: $(branch.r)")
end

function complexity(branch::SlicedBranch)
    code = branch.code
    isnothing(code) && return OMEinsum.OMEinsumContractionOrders.ContractionComplexity(0.0, 0.0, 0.0)
    return contraction_complexity(code, uniformsize(code, 2))
end
tc(branch::SlicedBranch) = complexity(branch).tc
sc(branch::SlicedBranch) = complexity(branch).sc
