function refine(code::DynamicNestedEinsum{LT}, size_dict::Dict{LT, Int}, refiner::TreeSARefiner, sc_target::Int, sc0::Number) where LT
    refined_code = code
    for i in 1:refiner.max_rounds
        refined_code = rethermalize(refined_code, size_dict, refiner.βs, refiner.ntrials, refiner.niters, sc_target)
    end
    sc = contraction_complexity(refined_code, size_dict).sc
    if sc > sc0
        # @info "Refiner did not improve the code, original sc = $sc0, got $sc, reoptimizing = $(refiner.reoptimize)"
        if refiner.reoptimize
            # refined_code = true_eincode(optimize_code(refined_code, size_dict, TreeSA(sc_target = sc_target)))
            reoptimized_code = rethermalize(refined_code, size_dict, 1.0:0.1:15.0, 5, 50, sc_target)
            resc = contraction_complexity(reoptimized_code, size_dict).sc
            refined_code = resc < sc ? reoptimized_code : refined_code
            # @info "Reoptimized the code, original sc = $sc0, refined sc = $sc, reoptimized sc = $resc"
        end
    end
    
    return refined_code
end

function refine(code::DynamicNestedEinsum{LT}, size_dict::Dict{LT, Int}, refiner::ReoptimizeRefiner, sc_target::Int, sc0::Number) where LT
    refined_code = true_eincode(optimize_code(code, size_dict, refiner.optimizer))
    (contraction_complexity(refined_code, size_dict).sc > sc0) && (@warn "Refiner did not improve the code, got $(contraction_complexity(refined_code, size_dict).sc) instead of $sc0")
    return refined_code
end