function CalcFVA(rxn_idx::Int64)
fva_max = readdlm("FVA/base_case/$rxn_idx/FVA_min")[:,rxn_idx] # MAX, MIN were saved in the wrong places
fva_min = readdlm("FVA/base_case/$rxn_idx/FVA_max")[:,rxn_idx]
return fva_max, fva_min
end

function CalcFVA(rxn_idx::Array{Int64})
r1 = rxn_idx[1]
r2 = rxn_idx[2]
fva_max = readdlm("FVA/base_case/$r1/FVA_min")[:,r1]-readdlm("FVA/base_case/$r1/FVA_min")[:,r2] # MAX, MIN were saved in the wrong places
#fva_min = readdlm("FVA/base_case/$r1/FVA_max")[:,r1]-readdlm("FVA/base_case/$r1/FVA_max")[:,r2]
fva_min = readdlm("FVA/base_case/$r2/FVA_min")[:,r1]-readdlm("FVA/base_case/$r2/FVA_min")[:,r2]
return fva_max, fva_min
end
