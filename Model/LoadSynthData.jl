function LoadSynthData()

species_list = readdlm("./data/species_list")
#Mean_raw = readdlm("./data/Mean")
#Std_raw = readdlm("./data/Std")
Upper_raw = readdlm("./data/Upper")
Lower_raw = readdlm("./data/Lower")
#Mean = Dict()
#Std = Dict()
Upper = Dict()
Lower = Dict()
for i in 1:length(species_list)
	rxn_idx = species_list[i]
#	Mean[rxn_idx] = Mean_raw[:,i]
#	Std[rxn_idx] = Std_raw[:,i]
	Upper[rxn_idx] = Upper_raw[:,i]
	Lower[rxn_idx] = Lower_raw[:,i]
end

return Upper, Lower
end
