include("Include.jl")

data_dictionary, TXTL_dictionary, species_list, Mean, Upper, Lower = LoadDictionaries()

metabolite_list = data_dictionary["list_of_metabolite_symbols"]
species_indices = Dict()
for i in 1:length(metabolite_list)
	key = metabolite_list[i]
	species_indices[key] = i
end

species_to_vary = readdlm("species_to_vary")
indices_to_vary = Int64[]
for i in 1:length(species_to_vary)
	push!(indices_to_vary,species_indices[species_to_vary[i]])
end

L = zeros(146,0)
for i in 1:238
	L = [L readdlm("Ens_j/$i/species_boolean")]
end
writedlm("species_boolean_array",L)

m = round(mean(L,2)*100,1)
tbl = sortrows([m metabolite_list][indices_to_vary,:],rev=true)[1:63,[2;1]]
writedlm("most_common_metabolites",tbl)

error_bc = 3608115.9638814344
L_best = zeros(146,0)
for i in 1:238
	error = readdlm("Ens_j/$i/state_error")[1]
	if error/error_bc < 1e-3
		L_best = [L_best readdlm("Ens_j/$i/species_boolean")]
	end
end
writedlm("species_boolean_array_subset",L_best)

m_best = round(mean(L_best,2)*100,1)
tbl_best = sortrows([m_best metabolite_list][indices_to_vary,:],rev=true)[1:63,[2;1]]
writedlm("most_common_metabolites_subset",tbl_best)

