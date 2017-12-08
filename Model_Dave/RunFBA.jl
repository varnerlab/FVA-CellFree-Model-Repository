function RunFBA(data_dictionary,TXTL_dictionary)

tic()
# set up initial conditions (use kinetic model data)
time_state_array = zeros(n_species,length_time)
time_state_array[:,1] = data_dictionary["initial_conditions"]
time_flux_array = zeros(n_rxn,length_time)
FBA = data_dictionary["default_flux_bounds_array"]
SBA = data_dictionary["species_bounds_array"]

Exit_flag = Int64[]
uptake_array = 0
flux_array = 0

species_constraint_index = copy(dataset_index)

# discretized dfba
time_array = experimental_time[1:end-1]
state_array = deepcopy(data_dictionary["initial_conditions"])
index = 1
Exit_flag = Int64[]

for time_index in time_array
	# constraint
	data_dictionary["objective_coefficient_array"] = zeros(number_of_fluxes)
	data_dictionary["objective_coefficient_array"][171] = -1
	data_dictionary = Bounds(data_dictionary,TXTL_dictionary,state_array)
	species_constraints = CalculateConstraints(data_dictionary, state_array, species_constraint_index, syn_data, index, tstep)
#	species_constraints = calculate_constraints(state_array, species_constraint_index, syn_data, time_index, tstep)
	data_dictionary["species_bounds_array"] = species_constraints
	data_dictionary["state_array"] = state_array

	# solve the lp problem -
	(objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(data_dictionary)
	push!(Exit_flag,exit_flag)

	# mass balance for concentration
	state_array = state_array + (uptake_array * tstep) # cell free mass balanace
	state_array[101] = data_dictionary["initial_conditions"][101]
	Z0 = objective_value
	
	index = index + 1
	time_state_array[:,index] = state_array
	time_flux_array[:,index] = flux_array

end # for time_index in time_array

writedlm("Exit_flag",Exit_flag)
error = CalcError(Upper,Lower,experimental_time,time_state_array,data_dictionary)
num_failed_tps = sum(5-Exit_flag)

time_elapsed = toc()

return time_state_array, error, num_failed_tps
end
