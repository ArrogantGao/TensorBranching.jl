using OMEinsum: getixsv, getiyv, LeafString, flatten, _flatten, isleaf, decorate
using OMEinsum.OMEinsumContractionOrders: IncidenceList, parse_eincode, eo2ct, ContractionTree
using TensorBranching.GenericTensorNetworks: rawcode
using OptimalBranching.OptimalBranchingCore: IPSolver, LPSolver

abstract type AbstractRegionSelector end

# the maximum intersection region selector try to find the region with the maximum intersection with the largest tensors, where large means larger that the threshold
@kwdef struct ScoreRS <: AbstractRegionSelector
    n_max::Int = 20
    strategy::Symbol = :neighbor # :neighbor
    loss::Symbol = :sc_score # num_uniques or sc_score ro bag_score
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
    region_selector::AbstractRegionSelector = ScoreRS() # select the region to branch, what else methods?
    table_solver::AbstractTableSolver = TensorNetworkSolver()
    brancher::AbstractBrancher = GreedyBrancher()
    refiner::AbstractRefiner = TreeSARefiner()
end

struct CompressedEinsum{LT}
    ixs::Vector{Vector{LT}}
    iy::Vector{LT}
    ct::ContractionTree
    function CompressedEinsum(ixs::Vector{Vector{LT}}, iy::Vector{LT}, ct::ContractionTree) where LT
        return new{LT}(ixs, iy, ct)
    end
end

AbstractTrees.nodevalue(ct::ContractionTree) = "-"
AbstractTrees.children(ct::ContractionTree) = [ct.left, ct.right]
Base.show(io::IO, ct::ContractionTree) = print_tree(io, ct)

function compress(code::NestedEinsum)
    ixs = getixsv(code)
    iy = getiyv(code)
    ct = ein2contraction_tree(code)
    return CompressedEinsum(ixs, iy, ct)
end

function compress(::Nothing)
    return nothing
end

function uncompress(ce::CompressedEinsum{LT}) where LT
    incidence_list = IncidenceList(Dict([i=>ix for (i, ix) in enumerate(ce.ixs)]))
    code = parse_eincode(incidence_list, ce.ct, vertices = collect(1:length(ce.ixs)))
    return decorate(code)
end

function uncompress(::Nothing)
    return nothing
end

struct SlicedBranch{INT, VT, RT}
    p::MISProblem{INT, VT}
    code::Union{CompressedEinsum, Nothing}
    r::RT
    function SlicedBranch(p::MISProblem{INT, VT}, ::Nothing, r::RT) where {INT, VT, RT}
        return new{INT, VT, RT}(p, nothing, r)
    end
    function SlicedBranch(g::SimpleGraph, weights::VT, code::DynamicNestedEinsum, r::RT) where {VT, RT}
        p = MISProblem(g, weights)
        INT = typeof(p).parameters[1]
        return new{INT, VT, RT}(p, compress(code), r)
    end
    function SlicedBranch(p::MISProblem{INT, VT}, code::CompressedEinsum, r::RT) where {INT, VT, RT}
        return new{INT, VT, RT}(p, code, r)
    end
    function SlicedBranch(p::MISProblem{INT, VT}, code::DynamicNestedEinsum, r::RT) where {INT, VT, RT}
        return new{INT, VT, RT}(p, compress(code), r)
    end
end
function Base.show(io::IO, branch::SlicedBranch{INT, VT, RT}) where {INT, VT, RT}
    print(io, "SlicedBranch{$INT, $VT, $RT}: ")
    gtype = VT isa UnitWeight ? "simple graph" : "weighted graph"
    print(io, "graph: {$(nv(branch.p.g)), $(ne(branch.p.g))} $gtype; ")
    cc = complexity(branch)
    print(io, "code complexity: sc: $(cc.sc), tc: $(cc.tc)")
    print(io, "; fixed weight: $(branch.r)")
end

add_r(branch::SlicedBranch{INT, VT, RT1}, r::RT2) where {INT, VT, RT1, RT2} = SlicedBranch(branch.p, branch.code, RT1(branch.r + r))

function complexity(branch::SlicedBranch)
    isnothing(branch.code) && return OMEinsum.OMEinsumContractionOrders.ContractionComplexity(0.0, 0.0, 0.0)
    code = uncompress(branch.code)
    return contraction_complexity(code, uniformsize(code, 2))
end
tc(branch::SlicedBranch) = complexity(branch).tc
sc(branch::SlicedBranch) = complexity(branch).sc
