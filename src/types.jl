abstract type AbstractRegionSelector end

# the maximum intersection region selector try to find the region with the maximum intersection with the largest tensors, where large means larger that the threshold
@kwdef struct MaxIntersectRS <: AbstractRegionSelector
    n_max::Int = 20
    strategy::Symbol = :mincut # :mincut or :neighbors
    loss::Symbol = :num_uniques # what else methods? may consider more complicated ones
end


abstract type AbstractSlicer end

@kwdef struct ContractionTreeSlicer <: AbstractSlicer
    sc_target::Int = 30
    loss_function::Symbol = :numofbags
    βs::StepRangeLen = 100.0:100.0 # range of βs for the rethermalization
    ntrials::Int = 1
    niters::Int = 200
    region_selector::AbstractRegionSelector = MaxIntersectRS() # select the region to branch, what else methods?
    table_solver::AbstractTableSolver = TensorNetworkSolver()
end