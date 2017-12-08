S = data_dictionary["stoichiometric_matrix"]
species_list = data_dictionary["list_of_metabolite_symbols"]
measurable_metabolites = readdlm("measurable_metabolites")[:,1]

U, S, V = svd(S,thin=false)

U_norm = abs(U)
for i in 1:size(U_norm,2)
	u = U_norm[:,i]
	u = (u-minimum(u))/(maximum(u)-minimum(u))
	U_norm[:,i] = u
end

s = (S-minimum(S))/(maximum(S)-minimum(S))

U_weighted = copy(U_norm)
for i in 1:size(U_weighted,2)
	u = U_weighted[:,i]
	u *= s[i]
	U_weighted[:,i] = u
end

importance = sum(U_weighted,2)

tmp = [importance species_list][measurable_metabolites,:]
importance_array = sortrows(tmp,rev=true)
writedlm("importance_array",importance_array)





