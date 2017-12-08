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

flux_subset = vcat(collect(1:51),collect(62:75))

Labels = ["Plus_M_13dpg_c_M_g3p_c"   
 "Plus_M_13dpg_c_M_g6p_c"   
 "Plus_M_13dpg_c_M_2ddg6p_c"
 "Plus_M_13dpg_c_M_6pgl_c"  
 "Plus_M_13dpg_c_M_e4p_c"   
 "Plus_M_13dpg_c_M_r5p_c"   
 "Plus_M_13dpg_c_M_ru5p_D_c"
 "Plus_M_13dpg_c_M_accoa_c" 
 "Plus_M_13dpg_c_M_akg_c"   
 "Plus_M_13dpg_c_M_glx_c"   
 "Plus_M_13dpg_c_M_icit_c"  
 "Plus_M_13dpg_c_M_oaa_c"   
 "Plus_M_13dpg_c_M_succoa_c"]

Labels = ["base_case"]

for constraint_index in 1:length(Labels)
	println(Labels[constraint_index])
	constraint_subset = Labels[constraint_index]
	FVA_MAX = zeros(300,194)
	FVA_MIN = zeros(300,194)
	for j in 1:length(single_fluxes)
		idx = single_fluxes[j]
		if in(idx,flux_subset)
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
		if in(idx1,flux_subset)
			fva_max = readdlm("FVA/$constraint_subset/$idx1/FVA_min") #min and max are switched upon saving
			fva_min = readdlm("FVA/$constraint_subset/$idx2/FVA_min")
			FVA_MAX[:,j+length(single_fluxes)] = fva_max[:,idx[1]]-fva_max[:,idx[2]]
			FVA_MIN[:,j+length(single_fluxes)] = fva_min[:,idx[1]]-fva_min[:,idx[2]]
		end
	end
	FVA_norm = norm(FVA_MAX-FVA_MIN)
	writedlm("FVA/$constraint_subset/FVA_norm_subset",FVA_norm)
#	time_state_array = readdlm("FVA/$constraint_subset/time_state_array.txt")
#	error_vector = CalcError(Upper,Lower,experimental_time,time_state_array,data_dictionary)
#	writedlm("FVA/$constraint_subset/error",sum(error_vector))
end

error = zeros(length(Labels),1)
FVA_norm_subset = zeros(length(Labels),1)
for constraint_index in 1:13
	println(Labels[constraint_index])
	constraint_subset = Labels[constraint_index]
	error[constraint_index] = readdlm("FVA/$constraint_subset/error")[1]
	FVA_norm_subset[constraint_index] = readdlm("FVA/$constraint_subset/FVA_norm_subset")[1]
end

tbl = [Labels error FVA_norm_subset]



