include("DataDictionary.jl")
include("FVA_max_min.jl")

data_dictionary = DataDictionary(0,0,0)
rxn_list = data_dictionary["list_of_reaction_strings"][1:191]
for i in 1:length(rxn_list)
	rxn_list[i] = split(rxn_list[i],"::")[1]
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

