include("DataDictionary.jl")
include("CalcError.jl")
#include("FVA_max_min.jl")

data_dictionary = DataDictionary(0,0,0)
rxn_list = data_dictionary["list_of_reaction_strings"][1:191]
for i in 1:length(rxn_list)
	rxn_list[i] = split(rxn_list[i],"::")[1]
end

t0 = 0
tf = 3
tstep = 0.01
data_dictionary["tstep"] = 0.01
experimental_time = t0:tstep:tf

species_list = readdlm("./data/species_list")
Mean_raw = readdlm("./data/Mean")
Std_raw = readdlm("./data/Std")
Upper_raw = readdlm("./data/Upper")
Lower_raw = readdlm("./data/Lower")
Mean = Dict()
Std = Dict()
Upper = Dict()
Lower = Dict()
for i in 1:length(species_list)
	rxn_idx = species_list[i]
	Mean[rxn_idx] = Mean_raw[:,i]
	Std[rxn_idx] = Std_raw[:,i]
	Upper[rxn_idx] = Upper_raw[:,i]
	Lower[rxn_idx] = Lower_raw[:,i]
end


rxn_list_subset = [[1:51;];[62:75;]]


I = Int64[]
for i in 1:length(rxn_list)
	tmp = rxn_list[i]
	if length(tmp) > 7
		if tmp[end-7:end] == "_reverse"
			push!(I,i)
		end
	end
end

net_fluxes = sort(vcat(I-1,I))
idx = collect(1:length(rxn_list))
single_fluxes = Int64[]
for i in idx
	if !in(i,net_fluxes)
		push!(single_fluxes,i)
	end
end


single_flux_labels = rxn_list[single_fluxes]
net_flux_labels = rxn_list[net_fluxes[1:2:end]].*"_net"
#single fluxes
dataset_index = [15; 18; 22; 24; 27; 29; 31; 34; 38; 41; 42; 59; 60; 61; 66; 69; 70; 76; 79; 83; 84; 86; 88; 89; 106; 110; 114; 121; 123; 126; 128; 130; 132; 133; 134; 135; 140]

single_species_additions = [3 "M_13dpg_c"; 5 "M_2pg_c"; 6 "M_3pg_c"; 44 "M_dhap_c"; 48 "M_f6p_c"; 50 "M_fdp_c"; 55 "M_g3p_c"; 56 "M_g6p_c"; 105 "M_pep_c"; 4 "M_2ddg6p_c"; 11 "M_6pgc_c"; 12 "M_6pgl_c"; 46 "M_e4p_c"; 117 "M_r5p_c"; 118 "M_ru5p_D_c"; 119 "M_s7p_c"; 138 "M_xu5p_D_c"; 16 "M_accoa_c"; 21 "M_akg_c"; 36 "M_cit_c"; 54 "M_fum_c"; 65 "M_glx_c"; 78 "M_icit_c"; 102 "M_oaa_c"; 124 "M_succoa_c"]
combo_species_additions = [3 "M_13dpg_c"; 44 "M_dhap_c"; 105 "M_pep_c"; 11 "M_6pgc_c"; 138 "M_xu5p_D_c";  36 "M_cit_c"; 54 "M_fum_c"]

constraint_index_array = collect(1:num_constraint_sets)

Labels = ["base_case"]

for constraint_index in 1:1
	println(constraint_index)
	constraint_subset = Labels[constraint_index]
	FVA_MAX = zeros(300,194)
	FVA_MIN = zeros(300,194)
	for j in 1:length(single_fluxes)
		idx = single_fluxes[j]
		if in(idx,rxn_list_subset)
		fva_max = readdlm("FVA/$constraint_subset/$idx/FVA_min")
		fva_min = readdlm("FVA/$constraint_subset/$idx/FVA_max")
		FVA_MAX[:,j] = fva_max[:,idx]
		FVA_MIN[:,j] = fva_min[:,idx]
		end
	end
	for j in 2:2:length(net_fluxes)
		idx = single_fluxes[j-1:j]
		idx1 = idx[1]
		idx2 = idx[2]
		if in(idx1,rxn_list_subset)
			fva_max = readdlm("FVA/$constraint_subset/$idx1/FVA_min") #min and max are switched upon saving
			fva_min = readdlm("FVA/$constraint_subset/$idx2/FVA_min")
			FVA_MAX[:,j+length(single_fluxes)] = fva_max[:,idx[1]]-fva_max[:,idx[2]]
			FVA_MIN[:,j+length(single_fluxes)] = fva_min[:,idx[1]]-fva_min[:,idx[2]]
		end
	end
	FVA_norm = norm(FVA_MAX-FVA_MIN)

	# calculate accuracy
	time_state_array = readdlm("FVA/$constraint_subset/time_state_array.txt")
	error_vector = CalcError(Upper,Lower,Mean,experimental_time,time_state_array,data_dictionary)
	error = sum(error_vector)

	writedlm("FVA/$constraint_subset/Error_vector_centralcarbon",error_vector)
	writedlm("FVA/$constraint_subset/Error_centralcarbon",error)
	writedlm("FVA/$constraint_subset/FVA_Max_centralcarbon",FVA_MAX)
	writedlm("FVA/$constraint_subset/FVA_Min_centralcarbon",FVA_MIN)
	writedlm("FVA/$constraint_subset/FVA_norm_centralcarbon",FVA_norm)
	writedlm("FVA/$constraint_subset/Labels_centralcarbon",vcat(single_flux_labels,net_flux_labels))
end
