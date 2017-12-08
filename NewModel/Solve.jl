tic()
include("Include.jl")

data_dictionary, TXTL_dictionary, species_list, Mean, Upper, Lower = LoadDictionaries()

t0_exp = 0
tf_exp = 3
tstep_exp = 0.01
experimental_time_exp = [t0_exp:tstep_exp:tf_exp;]

metabolite_list = data_dictionary["list_of_metabolite_symbols"]
rxn_list = data_dictionary["list_of_reaction_strings"]

(n_species,n_rxn) = size(data_dictionary["stoichiometric_matrix"])
number_of_fluxes = find(rxn_list.=="tRNA_charging_M_val_L_c_CAT::16.0*M_val_L_c+16.0*M_atp_c+16.0*tRNA+16.0*M_h2o_c --> 16.0*M_val_L_c_tRNA+16.0*M_amp_c+16.0*M_ppi_c")[1]

t0 = 0
tf = 3
tstep = 0.01
experimental_time = t0:tstep:tf
length_time = length(experimental_time)
data_dictionary["tstep"] = tstep

# this determines which parameters are fitted to data
dataset_index = [15; 18; 22; 24; 27; 29; 31; 34; 38; 41; 42; 59; 60; 61; 66; 69; 70; 76; 79; 83; 84; 86; 88; 89; 106; 110; 114; 121; 123; 126; 128; 130; 132; 133; 134; 135; 140]
non_dataset_index = collect(1:n_species)
deleteat!(non_dataset_index,dataset_index)

species_list = data_dictionary["list_of_metabolite_symbols"]
species_list_dataset = species_list[dataset_index]
species_list_non_dataset = species_list[non_dataset_index]

species_constraint_index = copy(dataset_index)

plot_color = "orangered"

num_constraint_sets = 1#size(species_constraint_array,1)

# not currently used
parameter = ones(202,2)
parameter[:,1] = parameter[:,1]*110
parameter[:,2] = parameter[:,2]*0.03

initial_conditions = data_dictionary["initial_conditions"]

# set up initial conditions (use kinetic model data)
time_state_array = zeros(n_species,length_time)
time_state_array[:,1] = initial_conditions
time_flux_array = zeros(n_rxn,length_time)
FBA = data_dictionary["default_flux_bounds_array"]
SBA = data_dictionary["species_bounds_array"]

EXIT_FLAG = zeros(length_time-1,num_constraint_sets)

overall_state_array = zeros()
# loop for all constraint subsets
FVA_run = false

constraint_index_array = collect(1:num_constraint_sets)

Exit_flag = Int64[]
uptake_array = 0
flux_array = 0

species_constraint_array = zeros(37,36)
for i in 1:37
	tmp = deleteat!(collect(1:37),i)
	species_constraint_array[i,:] = dataset_index[tmp]
end

Labels = "Delta_".*species_list[dataset_index]

Labels = ["base_case"]

if !isdir("FVA")
	mkdir("FVA")
end

for constraint_index in 1:1#size(species_constraint_array,1) #;println(Labels[constraint_index])
	constraint_subset = Labels[constraint_index]
	if !isdir("FVA/$constraint_subset")
		mkdir("FVA/$constraint_subset")
	end
	
	# discretized dfba
	state_array = deepcopy(initial_conditions)
	index = 1
	
	species_constraint_index = copy(dataset_index)
	
	syn_data_upper = zeros(n_species,length_time)
	syn_data_lower = zeros(n_species,length_time)
	syn_mean = zeros(n_species,length_time)
	for syn_index in species_constraint_index
		itp = interpolate((experimental_time_exp,),Mean[metabolite_list[syn_index]],Gridded(Linear()))
		Mean_itp = itp[experimental_time]
		itp = interpolate((experimental_time_exp,),Lower[metabolite_list[syn_index]],Gridded(Linear()))
		Lower_itp = itp[experimental_time]
		itp = interpolate((experimental_time_exp,),Upper[metabolite_list[syn_index]],Gridded(Linear()))
		Upper_itp = itp[experimental_time]
		syn_data_upper[syn_index,:] = Upper_itp
		syn_data_lower[syn_index,:] = Lower_itp
		syn_mean[syn_index,:] = Mean_itp
	end
	syn_data = Dict{AbstractString,Any}()
	syn_data["mean"] = syn_mean
	syn_data["upper"] = syn_data_upper
	syn_data["lower"] = syn_data_lower
	
	Exit_flag = Int64[]
	uptake_array = 0
	flux_array = 0
	
	FVA_min = Dict()
	FVA_max = Dict()
	Percentage_failed_tps_min = Dict()
	Percentage_failed_tps_max = Dict()
	for rxn_idx in 1:191
		FVA_min[rxn_idx] = zeros(length(t0:tstep:tf-tstep),191) # Only record 191 fluxes
		FVA_max[rxn_idx] = zeros(length(t0:tstep:tf-tstep),191) # Only record 191 fluxes
	end
	
	#dfba
	fva_time = experimental_time[1:end-1]
	for t in 1:length(fva_time)
		time_index = fva_time[t] ;println(time_index)
		# constraint
		data_dictionary["objective_coefficient_array"] = zeros(number_of_fluxes)
		data_dictionary["objective_coefficient_array"][171] = -1
#		data_dictionary["objective_coefficient_array"][171:191] = -1
#		data_dictionary["objective_coefficient_array"] = readdlm("Obj/objective_coefficient_array")
		data_dictionary = Bounds(data_dictionary,TXTL_dictionary,state_array) #;println(data_dictionary["default_flux_bounds_array"][167,1])
		species_constraints = calculate_constraints(state_array, parameter, species_constraint_index, syn_data, index, tstep)
		data_dictionary["species_bounds_array"] = species_constraints
		#  data_dictionary["default_flux_bounds_array"] = flux_constraints
		data_dictionary["state_array"] = state_array

		# solve the lp problem -
		(objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(data_dictionary) #;println(flux_array[168])#-flux_array[169])
		push!(Exit_flag,exit_flag)

		# mass balance for concentration
		state_array = state_array + (uptake_array * tstep) # cell free mass balanace
		state_array[101] = initial_conditions[101]
		Z0 = objective_value
		if FVA_run && in(time_index,collect(0.0:0.1:3.0));println("TIME: ",time_index)
			for flux_index in 1:191
				#constraint objective flux to value determined
				data_dictionary["default_flux_bounds_array"][171] = 1*Z0
				data_dictionary["default_flux_bounds_array"][171] = Z0

				# minimize flux
				data_dictionary["objective_coefficient_array"] = zeros(number_of_fluxes)
				data_dictionary["objective_coefficient_array"][flux_index] = 1
				(objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(data_dictionary)
				
				rxn_idx = copy(flux_index)
				
				FVA_min[rxn_idx][t,:] = flux_array[1:191]
				Percentage_failed_tps_min[rxn_idx] = sum(5-exit_flag)/length(experimental_time[1:end-1])

				# maximize flux
				data_dictionary["objective_coefficient_array"] = zeros(number_of_fluxes)
				data_dictionary["objective_coefficient_array"][flux_index] = -1
				(objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(data_dictionary)
				
				FVA_max[rxn_idx][t,:] = flux_array[1:191]
				Percentage_failed_tps_max[rxn_idx] = sum(5-exit_flag)/length(experimental_time[1:end-1])
			end
		end # if FVA_run
		
		index = index + 1
		time_state_array[:,index] = state_array
		time_flux_array[:,index] = flux_array

	end # for time_index in [t0:tstep:tf-tstep;]
	#writedlm("FVA/$constraint_subset/time_state_array.txt",time_state_array)
	
	error = CalcError(Upper,Lower,Mean,experimental_time,time_state_array,data_dictionary)
	#writedlm("FVA/$constraint_subset/error",error)
	
	for rxn_idx in 1:191
		if !isdir("FVA/$constraint_subset/$rxn_idx")
			mkdir("FVA/$constraint_subset/$rxn_idx")
		end
		#writedlm("FVA/$constraint_subset/$rxn_idx/FVA_min",FVA_min[rxn_idx])
		#writedlm("FVA/$constraint_subset/$rxn_idx/FVA_max",FVA_max[rxn_idx])
		#writedlm("FVA/$constraint_subset/$rxn_idx/Percentage_failed_tps_min",Percentage_failed_tps_min[rxn_idx])
		#writedlm("FVA/$constraint_subset/$rxn_idx/Percentage_failed_tps_max",Percentage_failed_tps_max[rxn_idx])
	end

#	if !isdir("ef/$constraint_subset")
#		mkdir("ef/$constraint_subset")
#	end
	println(sum(5-Exit_flag))
	#writedlm("ef/$constraint_subset/Exit_flag",Exit_flag)
	
end # for constraint_index in constraint_index_array

include("PlotSeparate.jl")
plot_flag = true
if plot_flag
#	PlotErrorbar(data_dictionary, Upper, Lower, plot_color)
#	PlotMain(plot_color)
	PlotSeparate(data_dictionary, Upper, Lower, plot_color)
end

time_elapsed = toc()


