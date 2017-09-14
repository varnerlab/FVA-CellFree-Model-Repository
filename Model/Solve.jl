include("Include.jl")
include("calculate_constraints.jl")
include("./data/Unpack.jl")
using PyPlot
using Interpolations


t0_exp = 0
tf_exp = 3
tstep_exp = 0.01
experimental_time_exp = [t0_exp:tstep_exp:tf_exp;]
Mean = Unpack("./data/Mean")
Std = Unpack("./data/Std")
Upper = Unpack("./data/Upper")
Lower = Unpack("./data/Lower")
data_dictionary = DataDictionary(0,0,0)#148x301
data_dictionary["species_drift_constant"] = 0.01 # percent allowed to drift for unconstraint species for each iteration (0.01 hr)
metabolite_list = data_dictionary["list_of_metabolite_symbols"]

# 182 translation_CAT
# 60 M_glc_D_c
# 142 PROTEIN_CAT
# load the data dictionary -

# error in generation of STM associated with RNAP and RIBOSOME since species appear on both side of reaction
# correction needed --- species taken out
data_dictionary["stoichiometric_matrix"][143:144,:] = 0;
data_dictionary["stoichiometric_matrix"][1,:] = 0 # gene_cat
#mrna
data_dictionary["stoichiometric_matrix"][145,180] = 1 # transcription
data_dictionary["stoichiometric_matrix"][1,180] = 0 # transcription
data_dictionary["stoichiometric_matrix"][145,182] = 0 # translation
(n_species,n_rxn) = size(data_dictionary["stoichiometric_matrix"])
t0 = 0
tf = 3
tstep = 0.001
experimental_time = [t0:tstep:tf;]
length_time = length(experimental_time)


#= notes: two possible implementations ---
  1. select subset of species to be constraint to synthetic data
    have everything else that are not constrained to data be open or a function of concentration
  2. select subset of species to be constraint to synthetic data
    have everything else that are not constrained to data constrained at steady state (sv=0)

=#
initial_conditions = ones(n_species,1)*0.1
for species_index in 1:n_species
  initial_conditions[species_index] = Mean[metabolite_list[species_index]][1]
  @show species_index,metabolite_list[species_index],Mean[metabolite_list[species_index]][1]
end
#initial_conditions = ones(n_species,1)*0.1
num_constraint_sets = 10
# this determines which parameters are fitted to data
species_constraint_index = [60,142,15,21,84,116,74,31,38,41,69,70,135,137,145]
tRNA_index = [23, 26, 28, 30, 43, 62, 64, 67, 78, 81, 86, 88, 91, 109, 113, 124, 130, 132, 134, 139]
AA_index = tRNA_index - 1
species_constraint_index = vcat(species_constraint_index,tRNA_index,AA_index)
#species_constraint_index = [1:146;]
#species_constraint_index = 1
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
#for syn_index in [1:length_time]
#  syn_data[60,syn_index] = 1*exp(-0.01*syn_index)
#end

# not currently used
parameter = ones(202,2)
parameter[:,1] = parameter[:,1]*110
parameter[:,2] = parameter[:,2]*0.03

# import synthetic data from kinetic model
#initial_conditions = readdlm("IC_data.dat")
#

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


overall_state_array = zeros()
# loop for all constraint subsets
for constraint_index = [1];
  # discretized dfba
  state_array = deepcopy(initial_conditions)
  index = 1

  #dfba
  for time_index in [t0:tstep:tf-tstep;]
    # constraint
      (species_constraints, flux_constraints) = calculate_constraints(state_array,parameter,species_constraint_index,syn_data,index,tstep)
      data_dictionary["species_bounds_array"] = species_constraints
      data_dictionary["default_flux_bounds_array"] = flux_constraints


      # solve the lp problem -
      (objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(data_dictionary)


      # mass balance for concentration
      state_array = state_array + (uptake_array * tstep) # cell free mass balanace
      state_array[103] = initial_conditions[103]


      index = index + 1
      time_state_array[:,index] = state_array
      time_flux_array[:,index] = flux_array



    end

  # evalution of performance

  # save all states

  # waters



end


# plotting simulations
figure("timecourse")
subplot(2,3,1)
species_number = 60 # gluc
a = plot(experimental_time,time_state_array[species_number,:],"k")
b = plot(experimental_time_exp,Mean[metabolite_list[species_number]],"blue")
fill_between(vec(experimental_time_exp),vec(Lower[metabolite_list[species_number]]),vec(Upper[metabolite_list[species_number]]),color="lightblue",lw=2)
title(metabolite_list[species_number])
subplot(2,3,2)
species_number = 15 # acetate
a = plot(experimental_time,time_state_array[species_number,:],"k")
b = plot(experimental_time_exp,Mean[metabolite_list[species_number]],"blue")
fill_between(vec(experimental_time_exp),vec(Lower[metabolite_list[species_number]]),vec(Upper[metabolite_list[species_number]]),color="lightblue",lw=2)
title(metabolite_list[species_number])
subplot(2,3,3)
species_number = 21 # akg
a = plot(experimental_time,time_state_array[species_number,:],"k")
b = plot(experimental_time_exp,Mean[metabolite_list[species_number]],"blue")
fill_between(vec(experimental_time_exp),vec(Lower[metabolite_list[species_number]]),vec(Upper[metabolite_list[species_number]]),color="lightblue",lw=2)
title(metabolite_list[species_number])
subplot(2,3,4)
species_number = 142 # cat
a = plot(experimental_time,time_state_array[species_number,:],"k")
b = plot(experimental_time_exp,Mean[metabolite_list[species_number]],"blue")
fill_between(vec(experimental_time_exp),vec(Lower[metabolite_list[species_number]]),vec(Upper[metabolite_list[species_number]]),color="lightblue",lw=2)
title(metabolite_list[species_number])
subplot(2,3,5)
species_number = 84 # lactate
a = plot(experimental_time,time_state_array[species_number,:],"k")
b = plot(experimental_time_exp,Mean[metabolite_list[species_number]],"blue")
fill_between(vec(experimental_time_exp),vec(Lower[metabolite_list[species_number]]),vec(Upper[metabolite_list[species_number]]),color="lightblue",lw=2)
title(metabolite_list[species_number])
subplot(2,3,6)
species_number = 116 # pyruvate
a = plot(experimental_time,time_state_array[species_number,:],"k")
b = plot(experimental_time_exp,Mean[metabolite_list[species_number]],"blue")
fill_between(vec(experimental_time_exp),vec(Lower[metabolite_list[species_number]]),vec(Upper[metabolite_list[species_number]]),color="lightblue",lw=2)
title(metabolite_list[species_number])

figure("misc")
subplot(2,3,1)
species_number = 143 # rnap
a = plot(experimental_time,time_state_array[species_number,:],"k")
b = plot(experimental_time_exp,Mean[metabolite_list[species_number]],"blue")
fill_between(vec(experimental_time_exp),vec(Lower[metabolite_list[species_number]]),vec(Upper[metabolite_list[species_number]]),color="lightblue",lw=2)
title(metabolite_list[species_number])

subplot(2,3,2)
species_number = 144 # rib
a = plot(experimental_time,time_state_array[species_number,:],"k")
b = plot(experimental_time_exp,Mean[metabolite_list[species_number]],"blue")
fill_between(vec(experimental_time_exp),vec(Lower[metabolite_list[species_number]]),vec(Upper[metabolite_list[species_number]]),color="lightblue",lw=2)
title(metabolite_list[species_number])
subplot(2,3,3)
species_number = 39 # co2
a = plot(experimental_time,time_state_array[species_number,:],"k")
b = plot(experimental_time_exp,Mean[metabolite_list[species_number]],"blue")
fill_between(vec(experimental_time_exp),vec(Lower[metabolite_list[species_number]]),vec(Upper[metabolite_list[species_number]]),color="lightblue",lw=2)
title(metabolite_list[species_number])
subplot(2,3,4)
species_number = 74 #
a = plot(experimental_time,time_state_array[species_number,:],"k")
b = plot(experimental_time_exp,Mean[metabolite_list[species_number]],"blue")
fill_between(vec(experimental_time_exp),vec(Lower[metabolite_list[species_number]]),vec(Upper[metabolite_list[species_number]]),color="lightblue",lw=2)
title(metabolite_list[species_number])
subplot(2,3,5)
species_number = 84 #
a = plot(experimental_time,time_state_array[species_number,:],"k")
b = plot(experimental_time_exp,Mean[metabolite_list[species_number]],"blue")
fill_between(vec(experimental_time_exp),vec(Lower[metabolite_list[species_number]]),vec(Upper[metabolite_list[species_number]]),color="lightblue",lw=2)
title(metabolite_list[species_number])
subplot(2,3,6)
species_number = 107 #
a = plot(experimental_time,time_state_array[species_number,:],"k")
b = plot(experimental_time_exp,Mean[metabolite_list[species_number]],"blue")
fill_between(vec(experimental_time_exp),vec(Lower[metabolite_list[species_number]]),vec(Upper[metabolite_list[species_number]]),color="lightblue",lw=2)
title(metabolite_list[species_number])

# 23 26 28 30 43 62 64 67 78 81 86 88 91109 113 124130 132 134 139
# plot tRNA charging
tRNA_index = [23 26 28 30 43 62 64 67 78 81 86 88 91 109 113 124 130 132 134 139]
count = 1
figure("tRNA")
for index in tRNA_index
  subplot(4,5,count)
  species_number = index # rib
  a = plot(experimental_time,time_state_array[species_number,:],"k")
  b = plot(experimental_time_exp,Mean[metabolite_list[species_number]],"blue")
  fill_between(vec(experimental_time_exp),vec(Lower[metabolite_list[species_number]]),vec(Upper[metabolite_list[species_number]]),color="lightblue",lw=2)
  title(metabolite_list[species_number])
  count = count+1
end

AA_index = [23 26 28 30 43 62 64 67 78 81 86 88 91 109 113 124 130 132 134 139] - 1
count = 1
figure("AA")
for index in AA_index
  subplot(4,5,count)
  species_number = index # rib
  a = plot(experimental_time,time_state_array[species_number,:],"k")
  b = plot(experimental_time_exp,Mean[metabolite_list[species_number]],"blue")
  fill_between(vec(experimental_time_exp),vec(Lower[metabolite_list[species_number]]),vec(Upper[metabolite_list[species_number]]),color="lightblue",lw=2)
  title(metabolite_list[species_number])
  count = count+1
end

#energetics
figure("energetics")
energetics_index = [24 31 38 41 69 70 135 137]
count = 1
for index in energetics_index
  subplot(2,4,count)
  species_number = index # rib
  a = plot(experimental_time,time_state_array[species_number,:],"k")
  b = plot(experimental_time_exp,Mean[metabolite_list[species_number]],"blue")
  fill_between(vec(experimental_time_exp),vec(Lower[metabolite_list[species_number]]),vec(Upper[metabolite_list[species_number]]),color="lightblue",lw=2)
  title(metabolite_list[species_number])
  count = count+1
end



figure("TXTL")
TXTL_index = [1 145 146 142]
# 1 GENE_CAT # 142 PROTEIN_CAT # 143 RIBOSOME # 144 RNAP# 145 mRNA_CAT # 146 tRNA
count = 1
for index in TXTL_index
  subplot(2,2,count)
  species_number = index # rib
  a = plot(experimental_time,time_state_array[species_number,:],"k")
  b = plot(experimental_time_exp,Mean[metabolite_list[species_number]],"blue")
  fill_between(vec(experimental_time_exp),vec(Lower[metabolite_list[species_number]]),vec(Upper[metabolite_list[species_number]]),color="lightblue",lw=2)
  title(metabolite_list[species_number])
  count = count+1
end
