function calculate_constraints(state_array, parameters, constraint_index, syn_data,time_index,tstep)
  alpha = 1.1 # looseness/tightness of constraints
  beta = 10
  n_species = length(state_array[:,1])
  syn_upper = syn_data["upper"]
  syn_lower = syn_data["lower"]
  syn_mean = syn_data["mean"]
  state_array = max(0,state_array)
#  # state array can only vary by a factor of 2 per iteration
#  species_drift_array = data_dictionary["species_drift_constant"]
#  if time_index*tstep < .15
#  	species_drift_constant = species_drift_array[1]
#  else
#  	species_drift_constant = species_drift_array[2]
#  end
  
  species_drift_amplitude = data_dictionary["species_drift_amplitude"]
  species_drift_decay_rate = data_dictionary["species_drift_decay_rate"]
  species_drift = species_drift_amplitude*exp(-time_index*tstep*species_drift_decay_rate)
  
  lower_bound = 0
  upper_bound = 100000000
  species_lower_bound = max(min(state_array*(-species_drift)/tstep,-0.00000000001/tstep),(lower_bound-state_array)/tstep);
  species_upper_bound = min(max(state_array*(species_drift)/tstep,0.00000000001/tstep),(upper_bound-state_array)/tstep);
  
#  species_lower_bound = (lower_bound-state_array)/tstep;
#  species_upper_bound = (upper_bound-state_array)/tstep;
  
#  println(mean(state_array)*(-species_drift)/tstep)
#  println((lower_bound-mean(state_array))/tstep)
#  println(mean(species_lower_bound))

  #species_lower_bound = species_lower_bound*0
  #species_upper_bound = species_lower_bound*0

  #@show syn_data,state_array
  #@show size(syn_data), size(state_array),constraint_index
  #@show syn_data[1,time_index+1]
  #@show state_array[1],constraint_index
  # constrain the species to synthetic data
  for c_index in constraint_index
    #if syn_data[c_index,time_index+1]-state_array[c_index] >= 0
      species_lower_bound[c_index] = (syn_lower[c_index,time_index+1]-state_array[c_index]) / tstep
      species_upper_bound[c_index] = (syn_upper[c_index,time_index+1]-state_array[c_index]) / tstep
    #else
      #species_lower_bound[c_index] = (syn_lower[c_index,time_index+1]-state_array[c_index])
    #  species_upper_bound[c_index] = (syn_upper[c_index,time_index+1]-state_array[c_index])
    #end

  end
  # oxygen
  species_lower_bound[101] = -10000
  species_upper_bound[101] = 10000
  #h2o
  species_lower_bound[72] = -10000
  species_upper_bound[72] = 10000

  # TXTL machinery
  #species_lower_bound[1] = 0 # GENE_CAT
  #species_upper_bound[1] = 0 # GENE_CAT

#  species_lower_bound[24] = -100
#  species_lower_bound[31] = -100
#  species_lower_bound[38] = -100
#  species_lower_bound[41] = -100
#  species_lower_bound[69] = -100
#  species_lower_bound[70] = -100
#  species_lower_bound[135] = -100
#  species_lower_bound[137] = -100
#  species_upper_bound[24] = 100
#  species_upper_bound[31] = 100
#  species_upper_bound[38] = 100
#  species_upper_bound[41] = 100
#  species_upper_bound[69] = 100
#  species_upper_bound[70] = 100
#  species_upper_bound[135] = 100
#  species_upper_bound[137] = 100

  # flux constraints

  # just start with Vmax and 0 for bounds
  #flux_lower_bound = 0*ones(202) #parameters[:,1].*state_array./(parameter[:,2]+state_array) / beta
  #flux_upper_bound = 1000*ones(202) #parameters[:,1].*state_array./(parameter[:,2]+state_array) * beta

  #flux_lower_bound[173:179] = 00
  #flux_upper_bound[173:179] = 10000


  # constraints on transcription_CAT
  n = 1
  k = 0.001
  #flux_lower_bound[180] = state_array[1]^n/(state_array[1]^n + k^n)/10 * 0
  #flux_upper_bound[180] = state_array[1]^n/(state_array[1]^n + k^n)


  # constraints for translation -- 145
  n = 1
  k = 0.00001

  #flux_lower_bound[182] = state_array[145]^n/(state_array[145]^n + k^n)/10 * 0
  #flux_upper_bound[182] = state_array[145]^n/(state_array[145]^n + k^n)
  # 173 M_o2_c --> []
  # 174 M_co2_c --> []
  # 175 M_pi_c --> []
  # 176 M_nh4_c --> []
  # 177 M_so4_c --> []
  # 178 M_h2o_c --> []
  # 179 M_h_e --> M_h_c

  #0.0	0.0	;	# 74 M_h_c
  #0.0	0.0	;	# 75 M_h_e

  species_constraints = [species_lower_bound species_upper_bound]
  #flux_constraints = [flux_lower_bound flux_upper_bound]

  return species_constraints#, flux_constraints
end
