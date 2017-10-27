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

Labels = String[]
FVA_MAX = zeros(300,165)
FVA_MIN = zeros(300,165)
for i in 1:length(single_fluxes)
	flux = single_fluxes[i]
	push!(Labels,rxn_list[flux])
	fva_max, fva_min = FVA_max_min(flux)
	FVA_MAX[:,i] = fva_max
	FVA_MIN[:,i] = fva_min
end
for i in 2:2:length(net_fluxes)
	fluxes = net_fluxes[i-1:i]
	key = rxn_list[fluxes[1]]*"_net"
	push!(Labels,key)
	fva_max, fva_min = FVA_max_min(net_fluxes[i-1:i])
	FVA_MAX[:,Int64(139+i/2)] = fva_max
	FVA_MIN[:,Int64(139+i/2)] = fva_min
end
FVA_SPREAD = FVA_MAX-FVA_MIN

fva_spread = zeros(size(FVA_SPREAD,2),1)
for i in 1:size(FVA_SPREAD,2)
	fva_spread[i] = norm(FVA_SPREAD[:,i])
end

sorted_fva_spread = sortrows([fva_spread Labels])

sorted_idx = sortrows([fva_spread Labels 1:165])[:,3]

if !isdir("raw_data")
	mkdir("raw_data")
end

#for i in sorted_idx
for i in 1:length(Labels)
	key = Labels[i]
	if !isdir("raw_data/$key")
		mkdir("raw_data/$key")
	end
	writedlm("raw_data/$key/FVA_max",FVA_MAX[:,i])
	writedlm("raw_data/$key/FVA_min",FVA_MIN[:,i])
end

