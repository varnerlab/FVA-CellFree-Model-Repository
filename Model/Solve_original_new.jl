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

species_list = readdlm("./data/species_list")
Mean_raw = readdlm("./data/Mean")
Std_raw = readdlm("./data/Std")
Upper_raw = readdlm("./data/Upper")
Lower_raw = readdlm("./data/Lower")
Mean = Dict()
Std = Dict()
Upper = Dict()
Lower = Dict()
for i in 1:length(species_list)
	rxn_idx = species_list[i]
	Mean[rxn_idx] = Mean_raw[:,i]
	Std[rxn_idx] = Std_raw[:,i]
	Upper[rxn_idx] = Upper_raw[:,i]
	Lower[rxn_idx] = Lower_raw[:,i]
end

data_dictionary["species_drift_amplitude"] = 2 # percent allowed to drift for unconstraint species for each iteration (0.01 hr)
data_dictionary["species_drift_decay_rate"] = 1
metabolite_list = data_dictionary["list_of_metabolite_symbols"]
rxn_list = data_dictionary["list_of_reaction_strings"]

(n_species,n_rxn) = size(data_dictionary["stoichiometric_matrix"])
t0 = 0
tf = 3
tstep = 0.01
experimental_time = t0:tstep:tf
length_time = length(experimental_time)
data_dictionary["tstep"] = tstep

tRNA_index = Int64[]
for i in 1:length(metabolite_list)
	rxn_idx = metabolite_list[i]
	if length(rxn_idx) > 5
		if rxn_idx[end-4:end] == "_tRNA"
			push!(tRNA_index,i)
		end
	end
end

initial_conditions = zeros(n_species,1)
for species_index in 1:n_species
	if in(species_index,tRNA_index)
		initial_conditions[species_index] = Mean[metabolite_list[species_index]][2]
	else
		initial_conditions[species_index] = Mean[metabolite_list[species_index]][1]
	end
end
initial_conditions[142] = .001
initial_conditions[144] = .0001

# this determines which parameters are fitted to data
dataset_index = [15; 18; 22; 24; 27; 29; 31; 34; 38; 41; 42; 59; 60; 61; 66; 69; 70; 76; 79; 83; 84; 86; 88; 89; 106; 110; 114; 121; 123; 126; 128; 130; 132; 133; 134; 135; 140]

species_list = data_dictionary["list_of_metabolite_symbols"]

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

Exit_flag = Int64[]
uptake_array = 0
flux_array = 0

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
time_array = experimental_time[1:end-1]
for time_index in time_array
	# constraint
	data_dictionary["objective_coefficient_array"] = zeros(number_of_fluxes)
	data_dictionary["objective_coefficient_array"][171] = -1
	data_dictionary = Bounds(data_dictionary,TXTL_dictionary,state_array)
	species_constraints = calculate_constraints(state_array,parameter,species_constraint_index,syn_data,index,tstep)
	data_dictionary["species_bounds_array"] = species_constraints
	data_dictionary["state_array"] = state_array

	# solve the lp problem -
	(objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(data_dictionary)
	push!(Exit_flag,exit_flag)

	# mass balance for concentration
	state_array = state_array + (uptake_array * tstep) # cell free mass balanace
	state_array[101] = initial_conditions[101]
	Z0 = objective_value
	
	index = index + 1
	time_state_array[:,index] = state_array
	time_flux_array[:,index] = flux_array

end # for time_index in time_array

error = CalcError(Upper,Lower,experimental_time,time_state_array,data_dictionary)
num_failed_tps = sum(5-Exit_flag)

time_elapsed = toc()


