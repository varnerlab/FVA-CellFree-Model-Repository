function FVA_max_min(rxn_idx::Int64)
fva_max = readdlm("FVA_results/Flux_$rxn_idx"*"_max")[:,rxn_idx] # MAX, MIN were saved in the wrong places
fva_min = readdlm("FVA_results/Flux_$rxn_idx"*"_min")[:,rxn_idx]
return fva_max, fva_min
end

function FVA_max_min(rxn_idx::Array{Int64})
r1 = rxn_idx[1]
r2 = rxn_idx[2]
fva_max = readdlm("FVA_results/Flux_$r1"*"_max")[:,r1]-readdlm("FVA_results/Flux_$r1"*"_max")[:,r2] # MAX, MIN were saved in the wrong places
#fva_min = readdlm("FVA_results/Flux_1_max")[:,r1]-readdlm("FVA_results/Flux_1_max")[:,r2]
fva_min = readdlm("FVA_results/Flux_$r2"*"_max")[:,r1]-readdlm("FVA_results/Flux_$r2"*"_max")[:,r2]
return fva_max, fva_min
end
