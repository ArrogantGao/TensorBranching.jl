abstract type AbstractRegionSelector end

# the maximum intersection region selector try to find the region with the maximum intersection with the largest tensors, where large means larger that the threshold
@kwdef struct MaxIntersectRS <: AbstractRegionSelector
    n_max::Int = 20
    strategy::Symbol = :mincut # :mincut or :neighbors
    loss::Symbol = :num_uniques # what else methods? may consider more complicated ones
end

# the sc_score region selector try to find the region with the minimum sc_score
@kwdef struct ScScoreRS <: AbstractRegionSelector
    n_max::Int = 20
    strategy::Symbol = :mincut # :mincut or :neighbors
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
    βs::StepRangeLen = 1000.0:1000.0 # range of βs for the rethermalization
    ntrials::Int = 10
    niters::Int = 200
    max_rounds::Int = 2
    reoptimize::Bool = true # setting this to true will reoptimize the code after each round of rethermalization if the resulting sc is larger than sc0
end

function refine(code::DynamicNestedEinsum{LT}, size_dict::Dict{LT, Int}, refiner::TreeSARefiner, sc_target::Int, sc0::Number) where LT
    refined_code = code
    for i in 1:refiner.max_rounds
        refined_code = @suppress rethermalize(refined_code, size_dict, refiner.βs, refiner.ntrials, refiner.niters, sc_target)
    end
    sc = contraction_complexity(refined_code, size_dict).sc
    if sc > sc0
        @warn "Refiner did not improve the code, original sc = $sc0, got $sc, reoptimizing = $(refiner.reoptimize)"
        if refiner.reoptimize
            refined_code = @suppress true_eincode(optimize_code(refined_code, size_dict, TreeSA(sc_target = sc0)))
            @info "Reoptimized the code, sc = $(contraction_complexity(refined_code, size_dict).sc)"
        end
    end
    
    return refined_code
end

@kwdef struct ReoptimizeRefiner <: AbstractRefiner
    optimizer::CodeOptimizer = TreeSA()
end

function refine(code::DynamicNestedEinsum{LT}, size_dict::Dict{LT, Int}, refiner::ReoptimizeRefiner, sc_target::Int, sc0::Number) where LT
    refined_code = @suppress true_eincode(optimize_code(code, size_dict, refiner.optimizer))
    (contraction_complexity(refined_code, size_dict).sc > sc0) && (@warn "Refiner did not improve the code, got $(contraction_complexity(refined_code, size_dict).sc) instead of $sc0")
    return refined_code
end

abstract type AbstractSlicer end

@kwdef struct ContractionTreeSlicer <: AbstractSlicer
    sc_target::Int = 30
    search_order::Symbol = :bfs # :bfs or :dfs
    region_selector::AbstractRegionSelector = MaxIntersectRS() # select the region to branch, what else methods?
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
    cc = contraction_complexity(branch.code, uniformsize(branch.code, 2))
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
