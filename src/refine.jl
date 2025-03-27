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

function refine(code::DynamicNestedEinsum{LT}, size_dict::Dict{LT, Int}, refiner::ReoptimizeRefiner, sc_target::Int, sc0::Number) where LT
    refined_code = @suppress true_eincode(optimize_code(code, size_dict, refiner.optimizer))
    (contraction_complexity(refined_code, size_dict).sc > sc0) && (@warn "Refiner did not improve the code, got $(contraction_complexity(refined_code, size_dict).sc) instead of $sc0")
    return refined_code
end