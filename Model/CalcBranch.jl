function CalcBranch(numerator_flux_index,denominator_flux_index)

if size(numerator_flux_index,2) < 1
	permutedims(numerator_flux_index,[2,1])
end
if size(denominator_flux_index,2) < 1
	permutedims(denominator_flux_index,[2,1])
end

tstep = 0.01
fraction = Float64[]
denominator_positive = Bool[]
for i in 1:194
	FVA = readdlm("FVA/base_case/$i/FVA_max")
	flux_array = FVA[:,numerator_flux_index]
	if size(flux_array,2) == 1
		net_flux_numerator = sum(flux_array)*tstep
	elseif size(flux_array,2) == 2
		net_flux_numerator = sum(flux_array[:,1]-flux_array[:,2])*tstep
	else
		throw("Wrong size")
	end
	flux_array = FVA[:,denominator_flux_index]
	if size(flux_array,2) == 1
		net_flux_denominator = sum(flux_array)*tstep
	elseif size(flux_array,2) == 2
		net_flux_denominator = sum(flux_array[:,1]-flux_array[:,2])*tstep
	else
		throw("Wrong size")
	end
	push!(fraction,net_flux_numerator/net_flux_denominator)
	push!(denominator_positive,net_flux_denominator>0)
	
	FVA = readdlm("FVA/base_case/$i/FVA_min")
	flux_array = FVA[:,numerator_flux_index]
	if size(flux_array,2) == 1
		net_flux_numerator = sum(flux_array)*tstep
	elseif size(flux_array,2) == 2
		net_flux_numerator = sum(flux_array[:,1]-flux_array[:,2])*tstep
	else
		throw("Wrong size")
	end
	flux_array = FVA[:,denominator_flux_index]
	if size(flux_array,2) == 1
		net_flux_denominator = sum(flux_array)*tstep
	elseif size(flux_array,2) == 2
		net_flux_denominator = sum(flux_array[:,1]-flux_array[:,2])*tstep
	else
		throw("Wrong size")
	end
	push!(fraction,net_flux_numerator/net_flux_denominator)
	push!(denominator_positive,net_flux_denominator>0)
end

return [minimum(fraction[denominator_positive]);maximum(fraction[denominator_positive])]
end
