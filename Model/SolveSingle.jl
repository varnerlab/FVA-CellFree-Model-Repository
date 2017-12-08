tic()
include("Include.jl")

# load the data dictionary -
data_dictionary = DataDictionary(0,0,0)
TXTL_dictionary = TXTLDictionary(0,0,0)
number_of_fluxes = length(data_dictionary["default_flux_bounds_array"][:,1])
number_of_species = length(data_dictionary["species_bounds_array"][:,1])
#Set objective reaction
data_dictionary["objective_coefficient_array"][171] = -1;
#-------------------------------------------------------------------------------------------------#Define case number
# 1 = Amino Acid Uptake & Synthesis
# 2 = Amino Acid Uptake w/o Synthesis
# 3 = Amino Acid Synthesis w/o Uptake
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
#TXTL_dictionary["RNAP_concentration_nM"] = 60 + (80-60)*rand(1)[1];
#TXTL_dictionary["RNAP_elongation_rate"] = 20 + (30-20)*rand(1)[1];
#TXTL_dictionary["RIBOSOME_concentration"] = 0.0012 + (0.0018-0.0012)*rand(1)[1];
#TXTL_dictionary["RIBOSOME_elongation_rate"] = 1.5 + (3-1.5)*rand(1)[1];
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
t0_exp = 0
tf_exp = 3
tstep_exp = 0.01
experimental_time_exp = [t0_exp:tstep_exp:tf_exp;]
Mean = Unpack("./data/Mean")
Std = Unpack("./data/Std")
Upper = Unpack("./data/Upper")
Lower = Unpack("./data/Lower")
#data_dictionary = DataDictionary(0,0,0)#148x301
data_dictionary["species_drift_constant"] = [2 2] # percent allowed to drift for unconstraint species for each iteration (0.01 hr)
metabolite_list = data_dictionary["list_of_metabolite_symbols"]
rxn_list = data_dictionary["list_of_reaction_strings"]

#data_dictionary["stoichiometric_matrix"][143:144,:] = 0;
#data_dictionary["stoichiometric_matrix"][1,:] = 0 # gene_cat
#mrna
#data_dictionary["stoichiometric_matrix"][145,180] = 1 # transcription
#data_dictionary["stoichiometric_matrix"][1,180] = 0 # transcription
#data_dictionary["stoichiometric_matrix"][145,182] = 0 # translation
(n_species,n_rxn) = size(data_dictionary["stoichiometric_matrix"])
t0 = 0
tf = 3
tstep = 0.001
experimental_time = [t0:tstep:tf;]
length_time = length(experimental_time)
data_dictionary["tstep"] = tstep

#notes: two possible implementations ---
#  1. select subset of species to be constraint to synthetic data
#    have everything else that are not constrained to data be open or a function of concentration
#  2. select subset of species to be constraint to synthetic data
#    have everything else that are not constrained to data constrained at steady state (sv=0)

initial_conditions = zeros(n_species,1)
for species_index in 1:n_species
  initial_conditions[species_index] = Mean[metabolite_list[species_index]][1]
#  @show species_index,metabolite_list[species_index],Mean[metabolite_list[species_index]][1]
end
initial_conditions[1]
initial_conditions[142]
initial_conditions[143]
initial_conditions[144]
initial_conditions[145]
initial_conditions[146]
#initial_conditions = ones(n_species,1)*0.1

# this determines which parameters are fitted to data
dataset_index = [15; 18; 22; 24; 27; 29; 31; 34; 38; 41; 42; 59; 60; 61; 66; 69; 70; 76; 79; 83; 84; 86; 88; 89; 106; 110; 114; 121; 123; 126; 128; 130; 132; 133; 134; 135; 140]
non_dataset_index = collect(1:number_of_species)
deleteat!(non_dataset_index,dataset_index)

species_list = data_dictionary["list_of_metabolite_symbols"]
species_list_dataset = species_list[dataset_index]
species_list_non_dataset = species_list[non_dataset_index]

#-------------------------------------------------------------------------------------------------
## Glycolysis
#glycolysis_group = ["M_glc_D_c"; "M_glc_D_c"; 
## Amino acids
#glu_group = ["M_ala_L_c"; "M_leu_L_c"; "M_phe_L_c"; "M_pro_L_c"; "M_tyr_L_c"; "M_val_L_c"]
#gln_group = ["M_asn_L_c"; "M_gln_L_c"; "M_his_L_c"; "M_trp_L_c"]
#asp_group = ["M_asp_L_c"; "M_cys_L_c"; "M_lys_L_c"; "M_met_L_c"]
#ser_group = ["M_gly_L_c"; "M_ile_L_c"; "M_ser_L_c"; "M_thr_L_c"]

#-------------------------------------------------------------------------------------------------

species_array = readdlm("Results/Other/species_array")
writedlm("Results/Other/species_array",species_array[2:end,:])

species_constraint_label = species_array[1,1]
println(species_constraint_label)
mkdir("Results/Other/$species_constraint_label")

species_constraint_index = species_array[1,2:end]
writedlm("Results/Other/$species_constraint_label/species_constraint_index", species_constraint_index)


#species_constraint_array = Array[]
#Labels = String[]
##for i in 1:length(dataset_index)
##	push!(Labels,"Delta_"*species_list_dataset[i])
##	tmp = copy(dataset_index)
##	deleteat!(tmp,i)
##	push!(species_constraint_array,tmp)
##end
#for i in 2:length(dataset_index)
#	for j in i+1:length(dataset_index)
#		push!(Labels,"Delta_"*species_list_dataset[i]*"_Delta_"*species_list_dataset[j])
#		tmp = copy(dataset_index)
#		deleteat!(tmp,[i;j])
#		push!(species_constraint_array,tmp)
#	end
#end
##for i in 1:length(non_dataset_index)
##	push!(Labels,"Plus_"*species_list_non_dataset[i])
##	tmp = copy(dataset_index)
##	push!(tmp,non_dataset_index[i])
##	push!(species_constraint_array,sort(tmp))
##end
##for i in 1:length(non_dataset_index)
##	for j in i+1:length(non_dataset_index)
##		push!(Labels,"Plus_"*species_list_non_dataset[i]*"_Plus_"*species_list_non_dataset[j])
##		tmp = copy(dataset_index)
##		push!(tmp,non_dataset_index[i])
##		push!(tmp,non_dataset_index[j])
##		push!(species_constraint_array,sort(tmp))
##	end
##end
#writedlm("species_constraint_labels",Labels)
#writedlm("species_constraint_array",species_constraint_array)

num_constraint_sets = 1 #size(species_constraint_array,1)

# not currently used
parameter = ones(202,2)
parameter[:,1] = parameter[:,1]*110
parameter[:,2] = parameter[:,2]*0.03

# import synthetic data from kinetic model
#initial_conditions = readdlm("IC_data.dat")

#data_dictionary["species_bounds_array"] = ones(202,2)
#data_dictionary["species_bounds_array"][:,1] = -1;
#data_dictionary["species_bounds_array"][:,2] = 1
# select subset of data to use as constraints
  # possible groups for analysis: glycolysis, TCA, PPP, AA synthesis, fatty acid synthesis, branching points

# set up initial conditions (use kinetic model data)
time_state_array = zeros(n_species,length_time)
time_state_array[:,1] = initial_conditions
time_flux_array = zeros(n_rxn,length_time)
FBA = data_dictionary["default_flux_bounds_array"]
SBA = data_dictionary["species_bounds_array"]

#EXIT_FLAG = zeros(length_time-1,num_constraint_sets)

overall_state_array = zeros()
# loop for all constraint subsets
FVA_run = true
fva_flux_array_1_max = zeros(number_of_fluxes,number_of_fluxes,num_constraint_sets)
fva_flux_array_2_max = zeros(number_of_fluxes,number_of_fluxes,num_constraint_sets)
fva_flux_array_3_max = zeros(number_of_fluxes,number_of_fluxes,num_constraint_sets)
fva_flux_array_4_max = zeros(number_of_fluxes,number_of_fluxes,num_constraint_sets)
fva_flux_array_1_min = zeros(number_of_fluxes,number_of_fluxes,num_constraint_sets)
fva_flux_array_2_min = zeros(number_of_fluxes,number_of_fluxes,num_constraint_sets)
fva_flux_array_3_min = zeros(number_of_fluxes,number_of_fluxes,num_constraint_sets)
fva_flux_array_4_min = zeros(number_of_fluxes,number_of_fluxes,num_constraint_sets)

#constraint_index_array = collect(1:num_constraint_sets)
constraint_index = 1

#Accuracy = zeros(num_constraint_sets,1)
#Precision = zeros(num_constraint_sets,1)
#for constraint_index in constraint_index_array;println(Labels[constraint_index])
	# discretized dfba
	state_array = deepcopy(initial_conditions)
	index = 1
	
#	species_constraint_index = species_constraint_array[constraint_index]
	
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

	#z_factor = 10
	#syn_data["upper"] = syn_mean+(syn_data_upper-syn_mean)*z_factor
	#syn_data["lower"] = max(0,syn_mean-(syn_mean-syn_data_lower)*z_factor)

	#for syn_index in [1:length_time]
	#  syn_data[60,syn_index] = 1*exp(-0.01*syn_index)
	#end
	
	Exit_flag = Int64[]
	Objective_value = Float64[]
	Flux_array = zeros(number_of_fluxes,0)
	Dual_array = zeros(number_of_fluxes,0)
	Uptake_array = zeros(number_of_species,0)
	
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
		push!(Objective_value,objective_value)
		Flux_array = [Flux_array flux_array]
		Dual_array = [Dual_array dual_array]
		Uptake_array = [Uptake_array uptake_array]

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
					eval(parse("fva_flux_array_"fva_flag*"_max"))[:,flux_index,constraint_index] = flux_array

					# maximize flux
					data_dictionary["objective_coefficient_array"] = zeros(number_of_fluxes)
					data_dictionary["objective_coefficient_array"][flux_index] = 1
					(objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(data_dictionary)
					eval(parse("fva_flux_array_"fva_flag*"_min"))[:,flux_index,constraint_index] = flux_array
					eval(parse("fva_flux_array_"fva_flag*"_diff = fva_flux_array_"fva_flag*"_max-fva_flux_array_"fva_flag*"_min"))
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
	accuracy = 1/sum(error)
#	Accuracy[constraint_index] = accuracy
	writedlm("Results/Other/$species_constraint_label/accuracy",accuracy)
	
	# "Precision" - based on FVA result
	fva_result = zeros(number_of_fluxes,0)
	for u in 1:4
	fva_flag = "$u"
	array = eval(parse("fva_flux_array_"fva_flag*"_diff"))
	diag_vector = diag(array[:,:,constraint_index])
	#		eval(parse("diag_vector = diag(fva_flux_array_"fva_flag*"_diff[:,:,constraint_index])"))
	fva_result = [fva_result diag_vector]
	#		eval(parse("fva_result_"fva_flag*" = diag_vector"))
	end
	precision = 1/norm(fva_result)
#	Precision[constraint_index] = precision
	writedlm("Results/Other/$species_constraint_label/precision",precision)

	# save all states
	
#	EXIT_FLAG[:,constraint_index] = Exit_flag
	writedlm("Results/Other/$species_constraint_label/Exit_flag",Exit_flag)
	writedlm("Results/Other/$species_constraint_label/Flux_array",Flux_array)
	writedlm("Results/Other/$species_constraint_label/Dual_array",Dual_array)
	writedlm("Results/Other/$species_constraint_label/Uptake_array",Uptake_array)
	
#end # for constraint_index in constraint_index_array

##writedlm("EXIT_FLAG",[collect(1:size(EXIT_FLAG,1)) EXIT_FLAG])
#writedlm("EXIT_FLAG",EXIT_FLAG)
#writedlm("Accuracy",Accuracy)
#writedlm("Precision",Precision)

#max_flux_array_1 = maximum(fva_flux_array_1_max,2)
#max_flux_array_2 = maximum(fva_flux_array_2_max,2)
#max_flux_array_3 = maximum(fva_flux_array_3_max,2)
#max_flux_array_4 = maximum(fva_flux_array_4_max,2)
#min_flux_array_1 = minimum(fva_flux_array_1_min,2)
#min_flux_array_2 = minimum(fva_flux_array_2_min,2)
#min_flux_array_3 = minimum(fva_flux_array_3_min,2)
#min_flux_array_4 = minimum(fva_flux_array_4_min,2)

#plot_flag = false
#if plot_flag
#	Plot()
#end

time_elapsed = toc()


