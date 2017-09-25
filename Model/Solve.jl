include("Include.jl")
include("Bounds.jl")
include("TXTLDictionary.jl")

# load the data dictionary -
data_dictionary = DataDictionary(0,0,0)
TXTL_dictionary = TXTLDictionary(0,0,0)
number_of_fluxes = length(data_dictionary["default_flux_bounds_array"][:,1])
number_of_species = length(data_dictionary["species_bounds_array"][:,1])
#Set objective reaction
data_dictionary["objective_coefficient_array"][194] = -1;

#=============================Cases=========================================#
#Define case number
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
#===========================================================================#
volume = TXTL_dictionary["volume"]

plasmid_concentration = 5;
gene_copy_number = (volume/1e9)*(6.02e23)*plasmid_concentration;
TXTL_dictionary["gene_copies"] = gene_copy_number
runs = number_of_fluxes
Protein_array = Float64[]
flux_ensemble = zeros(number_of_fluxes,runs)
for i = 1:runs
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


  # solve the lp problem -
  data_dictionary = Bounds(data_dictionary,TXTL_dictionary);
  flux_bounds = data_dictionary["default_flux_bounds_array"]
  flux_bounds[i,1:2] = 0
  data_dictionary["default_flux_bounds_array"] = flux_bounds
  (objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(data_dictionary)

  push!(Protein_array,abs(objective_value)*1e3)
  flux_ensemble[:,i] = flux_array
end


efficiency = (flux_ensemble[194,:]*2196)./(flux_ensemble[1,:]*21)
include("yield.jl")
res = zeros(1,3)
res[1,1] = Y*100
res[1,2] = mean(Protein_array)
res[1,3] = mean(efficiency)*100


#writedlm("Flux/Yield_Prod_$case.txt",res)
