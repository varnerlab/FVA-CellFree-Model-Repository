include("DataDictionary.jl")
include("CalcFVA.jl")

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
FVA_MAX = zeros(3000,165)
FVA_MIN = zeros(3000,165)
for i in 1:length(single_fluxes)
	flux = single_fluxes[i]
	push!(Labels,rxn_list[flux])
	fva_max, fva_min = CalcFVA(flux)
	FVA_MAX[:,i] = fva_max
	FVA_MIN[:,i] = fva_min
end
for i in 2:2:length(net_fluxes)
	fluxes = net_fluxes[i-1:i]
	key = rxn_list[fluxes[1]]*"_net"
	push!(Labels,key)
	fva_max, fva_min = CalcFVA(net_fluxes[i-1:i])
	FVA_MAX[:,Int64(139+i/2)] = fva_max
	FVA_MIN[:,Int64(139+i/2)] = fva_min
end
FVA_SPREAD = FVA_MAX-FVA_MIN

#fva_spread = zeros(size(FVA_SPREAD,2),1)
#for i in 1:size(FVA_SPREAD,2)
#	fva_spread[i] = norm(FVA_SPREAD[:,i])
#end
#sorted_fva_spread = sortrows([fva_spread Labels])

fva_spread_phase_1 = zeros(size(FVA_SPREAD,2),1)
fva_spread_phase_2 = zeros(size(FVA_SPREAD,2),1)
for i in 1:size(FVA_SPREAD,2)
	fva_spread_phase_1[i] = norm(FVA_SPREAD[1:1000,i])
	fva_spread_phase_2[i] = norm(FVA_SPREAD[1001:3000,i])
end
sorted_fva_spread_phase_1 = sortrows([fva_spread_phase_1 Labels])
sorted_fva_spread_phase_2 = sortrows([fva_spread_phase_2 Labels])
writedlm("sorted_fva_phase_1",sorted_fva_spread_phase_1)
writedlm("sorted_fva_phase_2",sorted_fva_spread_phase_2)


