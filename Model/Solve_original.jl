include("Include.jl")
include("calculate_constraints.jl")
include("./data/Unpack.jl")
include("Bounds.jl")
include("TXTLDictionary.jl")
using PyPlot
using Interpolations

include("Include.jl")
include("Bounds.jl")
include("TXTLDictionary.jl")

# load the data dictionary -
data_dictionary = DataDictionary(0,0,0)
TXTL_dictionary = TXTLDictionary(0,0,0)
number_of_fluxes = length(data_dictionary["default_flux_bounds_array"][:,1])
number_of_species = length(data_dictionary["species_bounds_array"][:,1])
#Set objective reaction
data_dictionary["objective_coefficient_array"][171] = -1;

#=============================Cases=========================================#
#Define case number
# 1 = Amino Acid Uptake & Synthesis
# 2 = Amino Acid Uptake w/o Synthesis
# 3 = Amino Acid Synthesis w/o Uptake
case = 3
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
#===========================================================================#
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




















t0_exp = 0
tf_exp = 3
tstep_exp = 0.01
experimental_time_exp = [t0_exp:tstep_exp:tf_exp;]
Mean = Unpack("./data/Mean")
Std = Unpack("./data/Std")
Upper = Unpack("./data/Upper")
Lower = Unpack("./data/Lower")
#data_dictionary = DataDictionary(0,0,0)#148x301#amplitude = 2, decay_rate = 1
data_dictionary["species_drift_amplitude"] = 10 # percent allowed to drift for unconstraint species for each iteration (0.01 hr)
data_dictionary["species_drift_decay_rate"] = 1
metabolite_list = data_dictionary["list_of_metabolite_symbols"]
rxn_list = data_dictionary["list_of_reaction_strings"]

# 182 translation_CAT
# 60 M_glc_D_c
# 142 PROTEIN_CAT
# load the data dictionary -

# error in generation of STM associated with RNAP and RIBOSOME since species appear on both side of reaction
# correction needed --- species taken out
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
experimental_time = [t0:tstep:tf;]
length_time = length(experimental_time)
data_dictionary["tstep"] = tstep

#= notes: two possible implementations ---
  1. select subset of species to be constraint to synthetic data
    have everything else that are not constrained to data be open or a function of concentration
  2. select subset of species to be constraint to synthetic data
    have everything else that are not constrained to data constrained at steady state (sv=0)

=#
initial_conditions = zeros(n_species,1)
for species_index in 1:n_species
  initial_conditions[species_index] = Mean[metabolite_list[species_index]][1]
  @show species_index,metabolite_list[species_index],Mean[metabolite_list[species_index]][1]
end
initial_conditions[1]
initial_conditions[142]
initial_conditions[143]
initial_conditions[144]
initial_conditions[145]
initial_conditions[146]
#initial_conditions = ones(n_species,1)*0.1
num_constraint_sets = 1
# this determines which parameters are fitted to data
species_constraint_index = [60,15,21,83,114,31,140,22,25,27]
tRNA_index = [23, 26, 28, 30, 43, 62, 64, 67, 77, 80, 85, 87, 90, 107, 111, 122, 127, 129, 131, 136]
AA_index = tRNA_index - 1

species_constraint_index = [15; 18; 22; 24; 27; 29; 31; 34; 38; 41; 42; 59; 60; 61; 66; 69; 70; 76; 79; 83; 84; 86; 88; 89; 106; 110; 114; 121; 123; 126; 128; 130; 132; 133; 134; 135; 140]
#species_constraint_index = vcat(species_constraint_index,AA_index)
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
FVA_run = false
fva_flux_array_1_max = zeros(number_of_fluxes,number_of_fluxes,num_constraint_sets)
fva_flux_array_2_max = zeros(number_of_fluxes,number_of_fluxes,num_constraint_sets)
fva_flux_array_3_max = zeros(number_of_fluxes,number_of_fluxes,num_constraint_sets)
fva_flux_array_4_max = zeros(number_of_fluxes,number_of_fluxes,num_constraint_sets)
fva_flux_array_1_min = zeros(number_of_fluxes,number_of_fluxes,num_constraint_sets)
fva_flux_array_2_min = zeros(number_of_fluxes,number_of_fluxes,num_constraint_sets)
fva_flux_array_3_min = zeros(number_of_fluxes,number_of_fluxes,num_constraint_sets)
fva_flux_array_4_min = zeros(number_of_fluxes,number_of_fluxes,num_constraint_sets)
for constraint_index = [1];
  # discretized dfba
  state_array = deepcopy(initial_conditions)
  index = 1

  #dfba
  for time_index in [t0:tstep:tf-tstep;]
      # constraint
      data_dictionary["objective_coefficient_array"] = zeros(number_of_fluxes)
      data_dictionary["objective_coefficient_array"][171] = -1
      data_dictionary = Bounds(data_dictionary,TXTL_dictionary,state_array);
      (species_constraints) = calculate_constraints(state_array,parameter,species_constraint_index,syn_data,index,tstep)
      data_dictionary["species_bounds_array"] = species_constraints
      #  data_dictionary["default_flux_bounds_array"] = flux_constraints
      data_dictionary["state_array"] = state_array

      # solve the lp problem -
      (objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(data_dictionary)


      # mass balance for concentration
      state_array = state_array + (uptake_array * tstep) # cell free mass balanace
      state_array[101] = initial_conditions[101]
      Z0 = objective_value
      if FVA_run == true
        # run FVA at time = 0hr, 1hr, 2hr, 3hr
        if time_index in [0, 1, 2, 3-tstep]
          if time_index == 0
            fva_flag = "1"
          end
          if time_index == 1
            fva_flag = "2"
          end
          if time_index == 2
            fva_flag = "3"
          end
          if time_index == 3-tstep
            fva_flag = "4"
          end
          for flux_index in [1:number_of_fluxes;]
            #constraint objective flux to value determined
            data_dictionary["default_flux_bounds_array"][171] = 1*Z0
            data_dictionary["default_flux_bounds_array"][171] = Z0


            # minimuze flux
            data_dictionary["objective_coefficient_array"] = zeros(number_of_fluxes)
            data_dictionary["objective_coefficient_array"][flux_index] = -1
            (objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(data_dictionary)

            #@show flux_index, constraint_index
            #@show flux_array
            eval(parse("fva_flux_array_"fva_flag*"_max"))[:,flux_index,constraint_index] = flux_array

            # maximize flux
            data_dictionary["objective_coefficient_array"] = zeros(number_of_fluxes)
            data_dictionary["objective_coefficient_array"][flux_index] = 1
            (objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(data_dictionary)
            eval(parse("fva_flux_array_"fva_flag*"_min"))[:,flux_index,constraint_index] = flux_array
          end
        end
      end



      index = index + 1
      time_state_array[:,index] = state_array
      time_flux_array[:,index] = flux_array



    end

  # evalution of performance

  # save all states

  # waters



end
max_flux_array_1 = maximum(fva_flux_array_1_max,2)
max_flux_array_2 = maximum(fva_flux_array_2_max,2)
max_flux_array_3 = maximum(fva_flux_array_3_max,2)
max_flux_array_4 = maximum(fva_flux_array_4_max,2)
min_flux_array_1 = minimum(fva_flux_array_1_min,2)
min_flux_array_2 = minimum(fva_flux_array_2_min,2)
min_flux_array_3 = minimum(fva_flux_array_3_min,2)
min_flux_array_4 = minimum(fva_flux_array_4_min,2)
plot_flag = true
# plotting simulations
if plot_flag == true
  figure("timecourse")
  subplot(2,3,1)
  species_number = 60 # gluc
  xcolor = "lightblue"
if species_number in species_constraint_index
  xcolor = "salmon"
end
a = plot(experimental_time,time_state_array[species_number,:],"k")
  b = plot(experimental_time_exp,Mean[metabolite_list[species_number]],"blue")
  fill_between(vec(experimental_time_exp),vec(Lower[metabolite_list[species_number]]),vec(Upper[metabolite_list[species_number]]),color=xcolor,lw=2)
  title(metabolite_list[species_number])
  subplot(2,3,2)
  species_number = 15 # acetate
  xcolor = "lightblue"
if species_number in species_constraint_index
  xcolor = "salmon"
end
a = plot(experimental_time,time_state_array[species_number,:],"k")
  b = plot(experimental_time_exp,Mean[metabolite_list[species_number]],"blue")
  fill_between(vec(experimental_time_exp),vec(Lower[metabolite_list[species_number]]),vec(Upper[metabolite_list[species_number]]),color=xcolor,lw=2)
  title(metabolite_list[species_number])
  subplot(2,3,3)
  species_number = 21 # akg
  xcolor = "lightblue"
if species_number in species_constraint_index
  xcolor = "salmon"
end
a = plot(experimental_time,time_state_array[species_number,:],"k")
  b = plot(experimental_time_exp,Mean[metabolite_list[species_number]],"blue")
  fill_between(vec(experimental_time_exp),vec(Lower[metabolite_list[species_number]]),vec(Upper[metabolite_list[species_number]]),color=xcolor,lw=2)
  title(metabolite_list[species_number])
  subplot(2,3,4)
  species_number = 140 # cat
  xcolor = "lightblue"
if species_number in species_constraint_index
  xcolor = "salmon"
end
a = plot(experimental_time,time_state_array[species_number,:],"k")
  b = plot(experimental_time_exp,Mean[metabolite_list[species_number]],"blue")
  fill_between(vec(experimental_time_exp),vec(Lower[metabolite_list[species_number]]),vec(Upper[metabolite_list[species_number]]),color=xcolor,lw=2)
  title(metabolite_list[species_number])
  subplot(2,3,5)
  species_number = 83 # lactate
  xcolor = "lightblue"
if species_number in species_constraint_index
  xcolor = "salmon"
end
a = plot(experimental_time,time_state_array[species_number,:],"k")
  b = plot(experimental_time_exp,Mean[metabolite_list[species_number]],"blue")
  fill_between(vec(experimental_time_exp),vec(Lower[metabolite_list[species_number]]),vec(Upper[metabolite_list[species_number]]),color=xcolor,lw=2)
  title(metabolite_list[species_number])
  subplot(2,3,6)
  species_number = 114 # pyruvate
  xcolor = "lightblue"
if species_number in species_constraint_index
  xcolor = "salmon"
end
a = plot(experimental_time,time_state_array[species_number,:],"k")
  b = plot(experimental_time_exp,Mean[metabolite_list[species_number]],"blue")
  fill_between(vec(experimental_time_exp),vec(Lower[metabolite_list[species_number]]),vec(Upper[metabolite_list[species_number]]),color=xcolor,lw=2)
  title(metabolite_list[species_number])

  figure("misc")
  subplot(2,3,1)
  species_number = 143 # rnap
  xcolor = "lightblue"
if species_number in species_constraint_index
  xcolor = "salmon"
end
a = plot(experimental_time,time_state_array[species_number,:],"k")
  b = plot(experimental_time_exp,Mean[metabolite_list[species_number]],"blue")
  fill_between(vec(experimental_time_exp),vec(Lower[metabolite_list[species_number]]),vec(Upper[metabolite_list[species_number]]),color=xcolor,lw=2)
  title(metabolite_list[species_number])

  subplot(2,3,2)
  species_number = 141 # rib
  xcolor = "lightblue"
if species_number in species_constraint_index
  xcolor = "salmon"
end
a = plot(experimental_time,time_state_array[species_number,:],"k")
  b = plot(experimental_time_exp,Mean[metabolite_list[species_number]],"blue")
  fill_between(vec(experimental_time_exp),vec(Lower[metabolite_list[species_number]]),vec(Upper[metabolite_list[species_number]]),color=xcolor,lw=2)
  title(metabolite_list[species_number])
  subplot(2,3,3)
  species_number = 144 # mRNA
  xcolor = "lightblue"
if species_number in species_constraint_index
  xcolor = "salmon"
end
a = plot(experimental_time,time_state_array[species_number,:],"k")
  b = plot(experimental_time_exp,Mean[metabolite_list[species_number]],"blue")
  fill_between(vec(experimental_time_exp),vec(Lower[metabolite_list[species_number]]),vec(Upper[metabolite_list[species_number]]),color=xcolor,lw=2)
  title(metabolite_list[species_number])
  subplot(2,3,4)
  species_number = 139 #
  xcolor = "lightblue"
if species_number in species_constraint_index
  xcolor = "salmon"
end
a = plot(experimental_time,time_state_array[species_number,:],"k")
  b = plot(experimental_time_exp,Mean[metabolite_list[species_number]],"blue")
  fill_between(vec(experimental_time_exp),vec(Lower[metabolite_list[species_number]]),vec(Upper[metabolite_list[species_number]]),color=xcolor,lw=2)
  title(metabolite_list[species_number])
  subplot(2,3,5)
  species_number = 142 #
  xcolor = "lightblue"
if species_number in species_constraint_index
  xcolor = "salmon"
end
a = plot(experimental_time,time_state_array[species_number,:],"k")
  b = plot(experimental_time_exp,Mean[metabolite_list[species_number]],"blue")
  fill_between(vec(experimental_time_exp),vec(Lower[metabolite_list[species_number]]),vec(Upper[metabolite_list[species_number]]),color=xcolor,lw=2)
  title(metabolite_list[species_number])
  subplot(2,3,6)
  species_number = 1 #
  xcolor = "lightblue"
if species_number in species_constraint_index
  xcolor = "salmon"
end
a = plot(experimental_time,time_state_array[species_number,:],"k")
  b = plot(experimental_time_exp,Mean[metabolite_list[species_number]],"blue")
  fill_between(vec(experimental_time_exp),vec(Lower[metabolite_list[species_number]]),vec(Upper[metabolite_list[species_number]]),color=xcolor,lw=2)
  title(metabolite_list[species_number])

  # 23 26 28 30 43 62 64 67 78 81 86 88 91109 113 124130 132 134 139
  # plot tRNA charging
  tRNA_index = [23, 26, 28, 30, 43, 62, 64, 67, 77, 80, 85, 87, 90, 107, 111, 122, 127, 129, 131, 136]
  count = 1
  figure("tRNA")
  for index in tRNA_index
    subplot(4,5,count)
    species_number = index # rib
    xcolor = "lightblue"
if species_number in species_constraint_index
  xcolor = "salmon"
end
a = plot(experimental_time,time_state_array[species_number,:],"k")
    b = plot(experimental_time_exp,Mean[metabolite_list[species_number]],"blue")
    fill_between(vec(experimental_time_exp),vec(Lower[metabolite_list[species_number]]),vec(Upper[metabolite_list[species_number]]),color=xcolor,lw=2)
    title(metabolite_list[species_number])
    count = count+1
  end

  AA_index = tRNA_index - 1
  count = 1
  figure("AA")
  for index in AA_index
    subplot(4,5,count)
    species_number = index # rib
    xcolor = "lightblue"
if species_number in species_constraint_index
  xcolor = "salmon"
end
a = plot(experimental_time,time_state_array[species_number,:],"k")
    b = plot(experimental_time_exp,Mean[metabolite_list[species_number]],"blue")
    fill_between(vec(experimental_time_exp),vec(Lower[metabolite_list[species_number]]),vec(Upper[metabolite_list[species_number]]),color=xcolor,lw=2)
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
    xcolor = "lightblue"
if species_number in species_constraint_index
  xcolor = "salmon"
end
a = plot(experimental_time,time_state_array[species_number,:],"k")
    b = plot(experimental_time_exp,Mean[metabolite_list[species_number]],"blue")
    fill_between(vec(experimental_time_exp),vec(Lower[metabolite_list[species_number]]),vec(Upper[metabolite_list[species_number]]),color=xcolor,lw=2)
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
    xcolor = "lightblue"
if species_number in species_constraint_index
  xcolor = "salmon"
end
a = plot(experimental_time,time_state_array[species_number,:],"k")
    b = plot(experimental_time_exp,Mean[metabolite_list[species_number]],"blue")
    fill_between(vec(experimental_time_exp),vec(Lower[metabolite_list[species_number]]),vec(Upper[metabolite_list[species_number]]),color=xcolor,lw=2)
    title(metabolite_list[species_number])
    count = count+1
  end
end
