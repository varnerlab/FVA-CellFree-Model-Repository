include("DataDictionary.jl")
function CalcErrorbar()

data_dictionary = DataDictionary(0,0,0)
species_list = vcat(data_dictionary["list_of_metabolite_symbols"], "M_axp_c", "M_gxp_c", "M_cxp_c", "M_uxp_c", "M_xxp_c")

num_dir = Int64(readdlm("TXTL/num_dir")[1])
species_states = Dict()
for i in 1:length(species_list)
	species_states[species_list[i]] = zeros(301,num_dir)
end

for i in 1:num_dir
	time_state_array = readdlm("TXTL/$i/time_state_array")
	for j in 1:size(time_state_array,1)
		species_states[species_list[j]][:,i] = vec(time_state_array[j,:])
	end
end
species_states["M_axp_c"] =  species_states["M_atp_c"]+species_states["M_adp_c"]+species_states["M_amp_c"]
species_states["M_gxp_c"] = species_states["M_gtp_c"]+species_states["M_gdp_c"]+species_states["M_gmp_c"]
species_states["M_cxp_c"] = species_states["M_ctp_c"]+species_states["M_cdp_c"]+species_states["M_cmp_c"]
species_states["M_uxp_c"] = species_states["M_utp_c"]+species_states["M_udp_c"]+species_states["M_ump_c"]
species_states["M_xxp_c"] = species_states["M_axp_c"]+species_states["M_gxp_c"]+species_states["M_cxp_c"]+species_states["M_uxp_c"]

z_val_for_95 = 1.95996

species_upper = Dict()
species_lower = Dict()
#for i in 1:length(species_list)
#	species_upper[species_list[i]] = mean(species_states[species_list[i]],2)+z_val_for_95*std(species_states[species_list[i]],2)
#	species_lower[species_list[i]] = max(0,mean(species_states[species_list[i]],2)-z_val_for_95*std(species_states[species_list[i]],2))
#end

num_tps = size(species_states[species_list[1]],1)
num_samples = size(species_states[species_list[1]],2)
for i in 1:length(species_list)
	upper = zeros(num_tps,1)
	lower = zeros(num_tps,1)
	for j in 1:size(species_states[species_list[i]],1)
		upper_idx = Int64(round(num_samples*.975))
		lower_idx = Int64(round(num_samples*.025))
		tmp = species_states[species_list[i]][j,:]
		upper[j] = sort(tmp)[upper_idx]
		lower[j] = sort(tmp)[lower_idx]
	end
	species_upper[species_list[i]] = upper
	species_lower[species_list[i]] = lower
end

return species_upper, species_lower
end
