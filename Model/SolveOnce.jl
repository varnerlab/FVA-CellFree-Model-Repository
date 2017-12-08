include("Include.jl")

function dFBA(objective_coefficient_array)

# load the data dictionary -
data_dictionary = DataDictionary(0,0,0)
TXTL_dictionary = TXTLDictionary(0,0,0)
number_of_fluxes = length(data_dictionary["default_flux_bounds_array"][:,1])
number_of_species = length(data_dictionary["species_bounds_array"][:,1])

#Set objective reaction
data_dictionary["objective_coefficient_array"][171] = -1;

case = 1
if case == 1
  data_dictionary["AASyn"] = 100;
  data_dictionary["AAUptake"] = 30
  data_dictionary["AASecretion"] = 0;
end
if case == 2
  data_dictionary["AASyn"] = 0;
  data_dictionary["AAUptake"] = 30
  data_dictionary["AASecretion"] = 0;
end
if case == 3
  data_dictionary["AASyn"] = 100;
  data_dictionary["AAUptake"] = 0
  data_dictionary["AASecretion"] = 100;
end
#-------------------------------------------------------------------------------------------------
volume = TXTL_dictionary["volume"]
plasmid_concentration = 5;
gene_copy_number = (volume/1e9)*(6.02e23)*plasmid_concentration;
TXTL_dictionary["gene_copies"] = gene_copy_number
data_dictionary = DataDictionary(0,0,0)
#Set objective reaction
data_dictionary["objective_coefficient_array"][194] = -1;
data_dictionary["AASyn"] = 100;
data_dictionary["AAUptake"] = 30
data_dictionary["AASecretion"] = 0;
data_dictionary["GlcUptake"] = 30#*rand(1)[1];
data_dictionary["Oxygen"] = 30#+2*rand(1)[1];
TXTL_dictionary["RNAP_concentration_nM"] = 75 #+ (80-60)*rand(1)[1];
TXTL_dictionary["RNAP_elongation_rate"] = 25 #+ (30-20)*rand(1)[1];
TXTL_dictionary["RIBOSOME_concentration"] = 0.0016# + (0.0018-0.0012)*rand(1)[1];
TXTL_dictionary["RIBOSOME_elongation_rate"] = 2# + (3-1.5)*rand(1)[1];
data_dictionary["Oxygen"] = 100;
data_dictionary["AAUptake"] = 1 #+ (1.5-0.5)*rand(1)[1];
data_dictionary["GlcUptake"] = 30# + (38-28)*rand(1)[1];
data_dictionary["AC"] = 8# + (12-6)*rand(1)[1];
data_dictionary["SUCC"] = 2# + (3-1)*rand(1)[1];
data_dictionary["PYR"] = 7# + (10-5)*rand(1)[1];
data_dictionary["MAL"] = 3 #+ (5-2)*rand(1)[1];
data_dictionary["LAC"] = 8 #+ (12-8)*rand(1)[1];
data_dictionary["ENERGY"] = 1#1*rand(1)[1];
data_dictionary["ALA"] = 2 #+ (4-3)*rand(1)[1];
data_dictionary["ASP"] = 1#+rand(1)[1];
data_dictionary["GLN"] = 0.25#*rand(1)[1];
data_dictionary["LysUptake"] = 1.5#*rand(1)[1];
data_dictionary["LysSecretion"] = 0.5#1*rand(1)[1];

t0_exp = 0
tf_exp = 3
tstep_exp = 0.01
experimental_time_exp = [t0_exp:tstep_exp:tf_exp;]
Mean = Unpack("./data/Mean")
Std = Unpack("./data/Std")
Upper = Unpack("./data/Upper")
Lower = Unpack("./data/Lower")

#amplitude = 2, decay_rate = 1
data_dictionary["species_drift_amplitude"] = 2 # percent allowed to drift for unconstraint species for each iteration (0.01 hr)
data_dictionary["species_drift_decay_rate"] = 1
metabolite_list = data_dictionary["list_of_metabolite_symbols"]
rxn_list = data_dictionary["list_of_reaction_strings"]

(n_species,n_rxn) = size(data_dictionary["stoichiometric_matrix"])
t0 = 0
tf = 3
tstep = 0.01
experimental_time = [t0:tstep:tf;]
length_time = length(experimental_time)
data_dictionary["tstep"] = tstep

initial_conditions = zeros(n_species,1)
for species_index in 1:n_species
  initial_conditions[species_index] = Mean[metabolite_list[species_index]][1]
#  @show species_index,metabolite_list[species_index],Mean[metabolite_list[species_index]][1]
end

# this determines which parameters are fitted to data
dataset_index = [15; 18; 22; 24; 27; 29; 31; 34; 38; 41; 42; 59; 60; 61; 66; 69; 70; 76; 79; 83; 84; 86; 88; 89; 106; 110; 114; 121; 123; 126; 128; 130; 132; 133; 134; 135; 140]
non_dataset_index = collect(1:number_of_species)
deleteat!(non_dataset_index,dataset_index)

species_list = data_dictionary["list_of_metabolite_symbols"]
species_list_dataset = species_list[dataset_index]
species_list_non_dataset = species_list[non_dataset_index]

#tRNA_index = [23; 26; 28; 30; 43; 62; 64; 67; 78; 81; 86; 88; 91; 109; 113; 124; 130; 132; 134; 139]
tRNA_index = Int64[]
for i in 1:length(metabolite_list)
	key = metabolite_list[i]
	if length(key) > 5
		if key[end-4:end] == "_tRNA"
			push!(tRNA_index,i)
		end
	end
end
#species_constraint_index = vcat(dataset_index,tRNA_index)
species_constraint_index = copy(dataset_index)

CAT_index = find(metabolite_list.=="PROTEIN_CAT")
species_constraint_index = species_constraint_index[find(species_constraint_index.!=CAT_index)]

ATP_index = find(metabolite_list.=="M_atp_c")
species_constraint_index = species_constraint_index[find(species_constraint_index.!=ATP_index)]

num_constraint_sets = 1#size(species_constraint_array,1)

# not currently used
parameter = ones(202,2)
parameter[:,1] = parameter[:,1]*110
parameter[:,2] = parameter[:,2]*0.03

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
fva_flux_array_1_max = zeros(number_of_fluxes,number_of_fluxes,num_constraint_sets)
fva_flux_array_2_max = zeros(number_of_fluxes,number_of_fluxes,num_constraint_sets)
fva_flux_array_3_max = zeros(number_of_fluxes,number_of_fluxes,num_constraint_sets)
fva_flux_array_4_max = zeros(number_of_fluxes,number_of_fluxes,num_constraint_sets)
fva_flux_array_1_min = zeros(number_of_fluxes,number_of_fluxes,num_constraint_sets)
fva_flux_array_2_min = zeros(number_of_fluxes,number_of_fluxes,num_constraint_sets)
fva_flux_array_3_min = zeros(number_of_fluxes,number_of_fluxes,num_constraint_sets)
fva_flux_array_4_min = zeros(number_of_fluxes,number_of_fluxes,num_constraint_sets)

constraint_index_array = collect(1:num_constraint_sets)

	state_array = deepcopy(initial_conditions)
	index = 1
	
	syn_data_upper = zeros(n_species,length_time)
	syn_data_lower = zeros(n_species,length_time)
	syn_mean = zeros(n_species,length_time)
	for syn_index in species_constraint_index
		itp = interpolate((experimental_time_exp,),Mean[metabolite_list[syn_index]],Gridded(Linear()))
		Mean_itp = itp[experimental_time]
		itp = interpolate((experimental_time_exp,),Std[metabolite_list[syn_index]],Gridded(Linear()))
		Std_itp = itp[experimental_time]
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

	#dfba
	for time_index in [t0:tstep:tf-tstep;] #;println(time_index)
		# constraint
		data_dictionary["objective_coefficient_array"] = zeros(number_of_fluxes)
		data_dictionary["objective_coefficient_array"][171] = -1
		data_dictionary = Bounds(data_dictionary,TXTL_dictionary,state_array)
		species_constraints = calculate_constraints(state_array,parameter,species_constraint_index,syn_data,index,tstep)
		data_dictionary["species_bounds_array"] = species_constraints
		#  data_dictionary["default_flux_bounds_array"] = flux_constraints
		data_dictionary["state_array"] = state_array

		# solve the lp problem -
		(objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(data_dictionary)
		push!(Exit_flag,exit_flag)

		# mass balance for concentration
		state_array = state_array + (uptake_array * tstep) # cell free mass balanace
		state_array[101] = initial_conditions[101]
		Z0 = objective_value
		if FVA_run
			# run FVA at time = 0hr, 1hr, 2hr, 3hr
			if time_index in [0, 1, 2, 3-tstep]
				if time_index == 0
					fva_flag = "1"
				elseif time_index == 1
					fva_flag = "2"
				elseif time_index == 2
					fva_flag = "3"
				elseif time_index == 3-tstep
					fva_flag = "4"
				end
				for flux_index in [1:number_of_fluxes;]
					#constraint objective flux to value determined
					data_dictionary["default_flux_bounds_array"][171] = 1*Z0
					data_dictionary["default_flux_bounds_array"][171] = Z0

					# minimize flux
					data_dictionary["objective_coefficient_array"] = zeros(number_of_fluxes)
					data_dictionary["objective_coefficient_array"][flux_index] = -1
					(objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(data_dictionary)

#		            @show flux_index, constraint_index
#		            @show flux_array
#					eval(parse("fva_flux_array_"fva_flag*"_max"))[:,flux_index,constraint_index] = flux_array

					# maximize flux
					data_dictionary["objective_coefficient_array"] = zeros(number_of_fluxes)
					data_dictionary["objective_coefficient_array"][flux_index] = 1
					(objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(data_dictionary)
#					eval(parse("fva_flux_array_"fva_flag*"_min"))[:,flux_index,constraint_index] = flux_array
#					eval(parse("fva_flux_array_"fva_flag*"_diff = fva_flux_array_"fva_flag*"_max-fva_flux_array_"fva_flag*"_min"))
				end # for flux_index in [1:number_of_fluxes;]
			end # if time_index in [0, 1, 2, 3-tstep]
		end # if FVA_run

		index = index + 1
		time_state_array[:,index] = state_array
		time_flux_array[:,index] = flux_array

	end # for time_index in [t0:tstep:tf-tstep;]

	# evalution of performance

	# "Accuracy" - based on dFBA result
	error = CalcError(Upper,Lower,experimental_time,time_state_array,data_dictionary)

