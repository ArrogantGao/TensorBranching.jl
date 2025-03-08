abstract type AbstractSlicer end

@kwdef struct ContractionTreeSlicer <: AbstractSlicer
    target_sc::Int = 30
    loss_function::Symbol = :numofbags
    rt_βs::StepRange = 100:0.05:101 # range of βs for the rethermalization
end
