abstract type AbstractRegionSelector end

# the maximum intersection region selector try to find the region with the maximum intersection with the largest tensors, where large means larger that the threshold
@kwdef struct MaxIntersectRS <: AbstractRegionSelector
    n_max::Int = 20
    strategy::Symbol = :mincut # :mincut or :neighbors
    loss::Symbol = :num_uniques # what else methods? may consider more complicated ones
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

@kwdef struct TreeSARefiner <: AbstractRefiner
    βs::StepRangeLen = 100.0:100.0 # range of βs for the rethermalization
    ntrials::Int = 1
    niters::Int = 100
end

abstract type AbstractSlicer end

@kwdef struct ContractionTreeSlicer <: AbstractSlicer
    sc_target::Int = 30
    region_selector::AbstractRegionSelector = MaxIntersectRS() # select the region to branch, what else methods?
    table_solver::AbstractTableSolver = TensorNetworkSolver()
    brancher::AbstractBrancher = GreedyBrancher()
    refiner::AbstractRefiner = TreeSARefiner()
end

struct SlicedBranch{T}
    g::SimpleGraph{T}
    code::DynamicNestedEinsum{T}
    r::Int
end
function Base.show(io::IO, branch::SlicedBranch{T}) where T
    print(io, "SlicedBranch{$T}: ")
    print(io, "graph: {$(nv(branch.g)), $(ne(branch.g))} simple graph; ")
    cc = contraction_complexity(branch.code, uniformsize(branch.code, 2))
    print(io, "code complexity: sc: $(cc.sc), tc: $(cc.tc)")
    print(io, "; fixed ones: $(branch.r)")
end

cc(branch::SlicedBranch) = contraction_complexity(branch.code, uniformsize(branch.code, 2))
tc(branch::SlicedBranch) = cc(branch).tc
sc(branch::SlicedBranch) = cc(branch).sc
