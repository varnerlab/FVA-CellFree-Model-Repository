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
#Mean = Unpack("./data/Mean")
#Std = Unpack("./data/Std")
#Upper = Unpack("./data/Upper")
#Lower = Unpack("./data/Lower")

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

#Lower["PROTEIN_CAT"] = .01*collect(0.0:.01:3.0)/3.0

#data_dictionary = DataDictionary(0,0,0)#148x301

#amplitude = 2, decay_rate = 1
data_dictionary["species_drift_amplitude"] = 2 # percent allowed to drift for unconstraint species for each iteration (0.01 hr)
data_dictionary["species_drift_decay_rate"] = 1
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
tstep = 0.01
experimental_time = t0:tstep:tf
length_time = length(experimental_time)
data_dictionary["tstep"] = tstep

#notes: two possible implementations ---
#  1. select subset of species to be constraint to synthetic data
#    have everything else that are not constrained to data be open or a function of concentration
#  2. select subset of species to be constraint to synthetic data
#    have everything else that are not constrained to data constrained at steady state (sv=0)

#tRNA_index = [23; 26; 28; 30; 43; 62; 64; 67; 78; 81; 86; 88; 91; 109; 113; 124; 130; 132; 134; 139]
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
#for species_index in 1:n_species
#  initial_conditions[species_index] = Mean[metabolite_list[species_index]][1]
##  @show species_index,metabolite_list[species_index],Mean[metabolite_list[species_index]][1]
#end
#initial_conditions[tRNA_index] =  .01 #initial_conditions[145]*4/20
for species_index in 1:n_species
	if in(species_index,tRNA_index)
		initial_conditions[species_index] = Mean[metabolite_list[species_index]][2]
	else
		initial_conditions[species_index] = Mean[metabolite_list[species_index]][1]
	end
end
#initial_conditions[142] = .001
initial_conditions[143] = .001
initial_conditions[144] = .0001

data_dictionary["initial_conditions"] = initial_conditions

# Original values:
GENE_IC = 5.19e-6
RNAP_IC = 6.75e-5
RIBO_IC = 0.002
mRNA_IC = 0.0001
tRNA_IC = 1.6
TX_Vmax = 122.7
plasmid_saturation_coefficient = 3.5e-6
TL_Vmax = 163.6
mRNA_saturation_coefficient = 0.045
mRNA_degradation_rate = 20

params = vcat(GENE_IC, RNAP_IC, RIBO_IC, mRNA_IC, tRNA_IC, TX_Vmax, plasmid_saturation_coefficient, TL_Vmax, mRNA_saturation_coefficient, mRNA_degradation_rate)

params[8] = 40

#params = [2.37878e-5 
# 0.0434249  
# 0.679137   
# 1.59025e-5 
# 0.337003   
# 0.175074   
# 1.73004e-7 
# 4.5247     
# 0.000711719
# 8.2838]

RNAP_IC = 1e-3
RIBO_IC = 0.002
RNAP_elongation_rate = 25
plasmid_saturation_coefficient = 116e-6
polysome_amplification = 5
RIBOSOME_elongation_rate = 2
mRNA_saturation_coefficient = 0.045
mRNA_degradation_rate = 10

params = vcat(RNAP_IC, RIBO_IC, RNAP_elongation_rate, plasmid_saturation_coefficient, polysome_amplification, RIBOSOME_elongation_rate, mRNA_saturation_coefficient, mRNA_degradation_rate)

data_dictionary["initial_conditions"][143] = params[1]
data_dictionary["initial_conditions"][141] = params[2]
TXTL_dictionary["RNAP_elongation_rate"] = params[3]
TXTL_dictionary["plasmid_saturation_coefficient"] = params[4]
TXTL_dictionary["polysome_amplification"] = params[5]
TXTL_dictionary["RIBOSOME_elongation_rate"] = params[6]
TXTL_dictionary["mRNA_saturation_coefficient"] = params[7]
TXTL_dictionary["mRNA_degradation_rate"] = params[8]

# this determines which parameters are fitted to data
dataset_index = [15; 18; 22; 24; 27; 29; 31; 34; 38; 41; 42; 59; 60; 61; 66; 69; 70; 76; 79; 83; 84; 86; 88; 89; 106; 110; 114; 121; 123; 126; 128; 130; 132; 133; 134; 135; 140]
#dataset_index = [15; 18; 22; 24; 27; 29; 34; 38; 42; 59; 60; 61; 66; 69; 76; 79; 83; 84; 86; 88; 89; 106; 110; 114; 121; 123; 126; 128; 130; 132; 133; 135; 140]
non_dataset_index = collect(1:number_of_species)
deleteat!(non_dataset_index,dataset_index)

species_list = data_dictionary["list_of_metabolite_symbols"]
species_list_dataset = species_list[dataset_index]
species_list_non_dataset = species_list[non_dataset_index]

#species_constraint_index = vcat(dataset_index,tRNA_index)
species_constraint_index = copy(dataset_index)
#species_constraint_index = dataset_index[1:10]

plot_color = "pink"

num_constraint_sets = 1#size(species_constraint_array,1)

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

EXIT_FLAG = zeros(length_time-1,num_constraint_sets)

overall_state_array = zeros()
# loop for all constraint subsets
FVA_run = false

constraint_index_array = collect(1:num_constraint_sets)

Exit_flag = Int64[]
uptake_array = 0
flux_array = 0











species_list = data_dictionary["list_of_metabolite_symbols"]
rxn_list = data_dictionary["list_of_reaction_strings"]

species_indices = Dict()
for i in 1:length(species_list)
	species_indices[species_list[i]] = i
end

rxn_indices = Dict()
for i in 1:length(rxn_list)
	key = String(split(rxn_list[i],':')[1])
	rxn_indices[key] = i
end

tstart = 0.0
tstop = 3.0
tstep_errorbar = data_dictionary["tstep"]
tstep_exp = (tstop-tstart)/(length(Upper["GENE_CAT"])-1)
t_array_exp = tstart:tstep_exp:tstop
t_array_synthetic = collect(tstart:tstep_errorbar:tstop)

Upper_itp = Dict()
Lower_itp = Dict()
Mean_itp = Dict()
for key in collect(keys(Upper))
	itp = interpolate((t_array_exp,),Upper[key],Gridded(Linear()))
	Upper_itp[key] = itp[t_array_synthetic]
	itp = interpolate((t_array_exp,),Lower[key],Gridded(Linear()))
	Lower_itp[key] = itp[t_array_synthetic]
	itp = interpolate((t_array_exp,),Mean[key],Gridded(Linear()))
	Mean_itp[key] = itp[t_array_synthetic]
end

tRNA_index = Int64[]
for i in 1:length(metabolite_list)
	key = metabolite_list[i]
	if length(key) >= 4
		if key[end-3:end] == "tRNA"
			push!(tRNA_index,i)
		end
	end
end

time_state_array_bc = readdlm("SubstrateUtil/time_state_array_bc")
time_state_array_glycolysis = readdlm("SubstrateUtil/time_state_array_glycolysis")
time_state_array_PP = readdlm("SubstrateUtil/time_state_array_PP")
time_state_array_ED = readdlm("SubstrateUtil/time_state_array_ED")
time_state_array_none = readdlm("SubstrateUtil/time_state_array_none")

figure("main",figsize=(4,3))
species = ["PROTEIN_CAT"]
for i in 1:length(species)
	key = species[i]
	species_number = species_indices[key]
	fill_between(vec(t_array_synthetic),vec(1000*Lower_itp[key]),vec(1000*Upper_itp[key]),color="lightblue")
	plot(experimental_time,1000*time_state_array_none[species_number,:],"g",lw=1.5)
#	plot(experimental_time,1000*time_state_array_ED[species_number,:],"g",lw=1.5)
	plot(experimental_time,1000*time_state_array_PP[species_number,:],"r",lw=1.5)
	plot(experimental_time,1000*time_state_array_glycolysis[species_number,:],"b",lw=1.5)
	plot(experimental_time,1000*time_state_array_bc[species_number,:],"k",lw=1.5)
	a = collect(axis())
	a[3] = 0
	axis(a)
	xlabel("Time (h)")
	ylabel(L"CAT ($\mu$M)")
	xticks(collect(0.0:1.0:3.0))
end
tight_layout()
savefig("SubstrateUtil/fig8_SubstrateUtilization.pdf")





