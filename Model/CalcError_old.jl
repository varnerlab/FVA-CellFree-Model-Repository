function CalcError(Upper,Lower,experimental_time,time_state_array,data_dictionary)

species = data_dictionary["list_of_metabolite_symbols"]

tstart = 0.0
tstop = 3.0
tstep_sim = data_dictionary["tstep"]
tstep_exp = (tstop-tstart)/(length(Upper["GENE_CAT"])-1)
t_array_sim = tstart:tstep_sim:tstop
t_array_exp = tstart:tstep_exp:tstop

error = zeros(length(species),1)
for i in 1:length(species)
	key = species[i]
	u = Upper[key]
	l = Lower[key]
	sim = time_state_array[i,:]
	itp = interpolate((t_array_sim,),sim,Gridded(Linear()))
	sim_itp = itp[t_array_exp]
	
	error[i] = sum(max(0,sim_itp-u)+max(0,l-sim_itp))
end

return error
end
