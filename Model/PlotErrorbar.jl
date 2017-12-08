include("CalcErrorbar.jl")
using LaTeXStrings

function PlotErrorbar(data_dictionary, Upper, Lower, plot_color)

#data_dictionary = DataDictionary(0,0,0)

species = [["M_glc_D_c" 100 10 40 "Glucose (mM)"],
["PROTEIN_CAT" 98 30 120 L"CAT ($\mu$M)"],
["M_pyr_c" 102 5 20 "Pyruvate (mM)"],
["M_lac_D_c" 104 5 20 "Lactate (mM)"],
["M_ac_c" 103 20 60 "Acetate (mM)"],
["M_succ_c" 118 2 6 "Succinate (mM)"],
["M_mal_L_c" 105 5 10 "Malate (mM)"],
["M_xxp_c" (18,24,31,59,69,70,34,38,41,132,133,134) 5 10 "Energy Total (mM)"],

["M_atp_c" 106 .5 2 "ATP (mM)"],
["M_adp_c" 107 .5 1 "ADP (mM)"],
["M_amp_c" 108 .5 1 "AMP (mM)"],
["M_axp_c" (18,24,31) .5 2.5 L"$\Sigma$ AXP (mM)"],
["M_gtp_c" 109 .5 2 "GTP (mM)"],
["M_gdp_c" 110 .5 1 "GDP (mM)"],
["M_gmp_c" 111 .5 1 "GMP (mM)"],
["M_gxp_c" (59,69,70) .5 2.5 L"$\Sigma$ GXP (mM)"],
["M_ctp_c" 115 .5 2 "CTP (mM)"],
["M_cdp_c" 116 .5 1 "CDP (mM)"],
["M_cmp_c" 117 .5 1 "CMP (mM)"],
["M_cxp_c" (34,38,41) .5 2.5 L"$\Sigma$ CXP (mM)"],
["M_utp_c" 112 .5 2 "UTP (mM)"],
["M_udp_c" 113 .5 1 "UDP (mM)"],
["M_ump_c" 114 .5 1 "UMP (mM)"],
["M_uxp_c" (132,133,134) .5 2.5 L"$\Sigma$ UXP (mM)"],

["M_ala_L_c" 126 2 8 "ALA (mM)"],
["M_arg_L_c" 36 1 4 "ARG (mM)"],
["M_asn_L_c" 122 1 4 "ASN (mM)"],
["M_asp_L_c" 119 1 4 "ASP (mM)"],
["M_cys_L_c" 123 1 4 "CYS (mM)"],
["M_gln_L_c" 137 1 5 "GLN (mM)"],
["M_glu_L_c" 136 50 200 "GLU (mM)"],
["M_gly_L_c" 120 1 4 "GLY (mM)"],
["M_his_L_c" 125 1 4 "HIS (mM)"],
["M_ile_L_c" 121 1 4 "ILE (mM)"],
["M_leu_L_c" 135 1 4 "LEU (mM)"],
["M_lys_L_c" 124 1 4 "LYS (mM)"],
["M_met_L_c" 134 1 4 "MET (mM)"],
["M_phe_L_c" 127 1 4 "PHE (mM)"],
["M_pro_L_c" 128 1 4 "PRO (mM)"],
["M_ser_L_c" 129 1 4 "SER (mM)"],
["M_thr_L_c" 130 1 4 "THR (mM)"],
["M_trp_L_c" 131 1 4 "TRP (mM)"],
["M_tyr_L_c" 132 1 4 "TYR (mM)"],
["M_val_L_c" 133 1 4 "VAL (mM)"]]

carbon_species = permutedims(reshape(species[1:8,1],(4,2)),(2,1))
energy_species = reshape(species[9:24,1],(4,4))
amino_species = permutedims(reshape(species[25:44,1],(4,5)),(2,1))

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
#Mean_itp = Dict()
for key in collect(keys(Upper))
	itp = interpolate((t_array_exp,),Upper[key],Gridded(Linear()))
	Upper_itp[key] = itp[t_array_synthetic]
	itp = interpolate((t_array_exp,),Lower[key],Gridded(Linear()))
	Lower_itp[key] = itp[t_array_synthetic]
#	itp = interpolate((t_array_exp,),Mean[key],Gridded(Linear()))
#	Mean_itp[key] = itp[t_array_synthetic]
end
Upper_itp["M_axp_c"] = Upper_itp["M_atp_c"]+Upper_itp["M_adp_c"]+Upper_itp["M_amp_c"]
Upper_itp["M_gxp_c"] = Upper_itp["M_gtp_c"]+Upper_itp["M_gdp_c"]+Upper_itp["M_gmp_c"]
Upper_itp["M_cxp_c"] = Upper_itp["M_ctp_c"]+Upper_itp["M_cdp_c"]+Upper_itp["M_cmp_c"]
Upper_itp["M_uxp_c"] = Upper_itp["M_utp_c"]+Upper_itp["M_udp_c"]+Upper_itp["M_ump_c"]
Upper_itp["M_xxp_c"] = Upper_itp["M_axp_c"]+Upper_itp["M_gxp_c"]+Upper_itp["M_cxp_c"]+Upper_itp["M_uxp_c"]
Lower_itp["M_axp_c"] = Lower_itp["M_atp_c"]+Lower_itp["M_adp_c"]+Lower_itp["M_amp_c"]
Lower_itp["M_gxp_c"] = Lower_itp["M_gtp_c"]+Lower_itp["M_gdp_c"]+Lower_itp["M_gmp_c"]
Lower_itp["M_cxp_c"] = Lower_itp["M_ctp_c"]+Lower_itp["M_cdp_c"]+Lower_itp["M_cmp_c"]
Lower_itp["M_uxp_c"] = Lower_itp["M_utp_c"]+Lower_itp["M_udp_c"]+Lower_itp["M_ump_c"]
Lower_itp["M_xxp_c"] = Lower_itp["M_axp_c"]+Lower_itp["M_gxp_c"]+Lower_itp["M_cxp_c"]+Lower_itp["M_uxp_c"]
#Mean_itp["M_axp_c"] = Mean_itp["M_atp_c"]+Mean_itp["M_adp_c"]+Mean_itp["M_amp_c"]
#Mean_itp["M_gxp_c"] = Mean_itp["M_gtp_c"]+Mean_itp["M_gdp_c"]+Mean_itp["M_gmp_c"]
#Mean_itp["M_cxp_c"] = Mean_itp["M_ctp_c"]+Mean_itp["M_cdp_c"]+Mean_itp["M_cmp_c"]
#Mean_itp["M_uxp_c"] = Mean_itp["M_utp_c"]+Mean_itp["M_udp_c"]+Mean_itp["M_ump_c"]
#Mean_itp["M_xxp_c"] = Mean_itp["M_axp_c"]+Mean_itp["M_gxp_c"]+Mean_itp["M_cxp_c"]+Mean_itp["M_uxp_c"]

species_upper, species_lower = CalcErrorbar()
experimental_time = collect(0:3.0/(length(species_upper["PROTEIN_CAT"])-1):3.0)

# Carbon
f_Carbon,axarr_Carbon = plt[:subplots](size(carbon_species,1),size(carbon_species,2))
f_Carbon[:set_figheight](3.5)
f_Carbon[:set_figwidth](10)
for i in 1:size(carbon_species,1)
	for j in 1:size(carbon_species,2)
		key = carbon_species[i,j][1]
		if in(key,keys(species_indices))
			species_index = species_indices[key]
		else
			species_index = collect(carbon_species[i,j][2])
		end
		y_step = carbon_species[i,j][3]
		y_max = carbon_species[i,j][4]
		y_label = carbon_species[i,j][5]
		
		# Check if CAT
		if key == "PROTEIN_CAT" # Convert from mM to μM
			conversion_factor = 1000
		else
			conversion_factor = 1
		end
		
		if i < size(carbon_species,1)
			axarr_Carbon[i,j][:xaxis][:set_major_formatter](plt[:NullFormatter]())
		else
			axarr_Carbon[i,j][:set_xlabel]("Time (h)")
		end
		
		axarr_Carbon[i,j][:set_ylabel](y_label)
		axarr_Carbon[i,j][:axis]([0,tstop,0,y_max])
		axarr_Carbon[i,j][:set_xticks](collect(0:1:tstop))
		axarr_Carbon[i,j][:set_yticks](collect(0:y_step:y_max))
		
		# Plot synthetic data
		shade_color = "lightblue"
		axarr_Carbon[i,j][:fill_between](vec(t_array_synthetic),vec(conversion_factor*Lower_itp[key]),vec(conversion_factor*Upper_itp[key]),color=shade_color)
		
#		# Plot simulations
#		if typeof(species_index) == Int64
#		    dfba = time_state_array[species_index,:]
#		elseif typeof(species_index) == Array{Int64,1}
#		    dfba = sum(time_state_array[species_index,:],1)
#		else
#		    throw("Wrong type")
#		end
#		axarr_Carbon[i,j][:plot](experimental_time,vec(conversion_factor*dfba),plot_color,lw=1.5)
		axarr_Carbon[i,j][:fill_between](experimental_time,vec(conversion_factor*species_upper[key]),vec(species_lower[key]),color=plot_color,alpha=0.4)
	end
end
tight_layout()
savefig("Carbon.pdf")
close()

# Amino
f_Amino,axarr_Amino = plt[:subplots](size(amino_species,1),size(amino_species,2))
f_Amino[:set_figheight](7.9)
f_Amino[:set_figwidth](10)
for i in 1:size(amino_species,1)
	for j in 1:size(amino_species,2)
		key = amino_species[i,j][1]
		if in(key,keys(species_indices))
			species_index = species_indices[key]
		else
			species_index = collect(amino_species[i,j][2])
		end
		y_step = amino_species[i,j][3]
		y_max = amino_species[i,j][4]
		y_label = amino_species[i,j][5]
		
		# Check if CAT
		if key == "PROTEIN_CAT" # Convert from mM to μM
			conversion_factor = 1000
		else
			conversion_factor = 1
		end
		
		if i < size(amino_species,1)
			axarr_Amino[i,j][:xaxis][:set_major_formatter](plt[:NullFormatter]())
		else
			axarr_Amino[i,j][:set_xlabel]("Time (h)")
		end
		
		axarr_Amino[i,j][:set_ylabel](y_label)
		axarr_Amino[i,j][:axis]([0,tstop,0,y_max])
		axarr_Amino[i,j][:set_xticks](collect(0:1:tstop))
		axarr_Amino[i,j][:set_yticks](collect(0:y_step:y_max))
		
		# Plot synthetic data
		shade_color = "lightblue"
		axarr_Amino[i,j][:fill_between](vec(t_array_synthetic),vec(conversion_factor*Lower_itp[key]),vec(conversion_factor*Upper_itp[key]),color=shade_color)
		
#		# Plot simulations
#		if typeof(species_index) == Int64
#		    dfba = time_state_array[species_index,:]
#		elseif typeof(species_index) == Array{Int64,1}
#		    dfba = sum(time_state_array[species_index,:],1)
#		else
#		    throw("Wrong type")
#		end
#		axarr_Amino[i,j][:plot](experimental_time,vec(conversion_factor*dfba),plot_color,lw=1.5)
		axarr_Amino[i,j][:fill_between](experimental_time,vec(conversion_factor*species_upper[key]),vec(species_lower[key]),color=plot_color,alpha=0.4)
	end
end
tight_layout()
savefig("Amino.pdf")
close()

# Energy
f_Energy,axarr_Energy = plt[:subplots](size(energy_species,1),size(energy_species,2))
f_Energy[:set_figheight](6.5)
f_Energy[:set_figwidth](10)
for i in 1:size(energy_species,1)
	for j in 1:size(energy_species,2)
		key = energy_species[i,j][1]
		if in(key,keys(species_indices))
			species_index = species_indices[key]
		else
			species_index = collect(energy_species[i,j][2])
		end
		y_step = energy_species[i,j][3]
		y_max = energy_species[i,j][4]
		y_label = energy_species[i,j][5]
		
		# Check if CAT
		if key == "PROTEIN_CAT" # Convert from mM to μM
			conversion_factor = 1000
		else
			conversion_factor = 1
		end
		
		if i < size(energy_species,1)
			axarr_Energy[i,j][:xaxis][:set_major_formatter](plt[:NullFormatter]())
		else
			axarr_Energy[i,j][:set_xlabel]("Time (h)")
		end
		
		axarr_Energy[i,j][:set_ylabel](y_label)
		axarr_Energy[i,j][:axis]([0,tstop,0,y_max])
		axarr_Energy[i,j][:set_xticks](collect(0:1:tstop))
		axarr_Energy[i,j][:set_yticks](collect(0:y_step:y_max))
		
		# Plot synthetic data
		shade_color = "lightblue"
		axarr_Energy[i,j][:fill_between](vec(t_array_synthetic),vec(conversion_factor*Lower_itp[key]),vec(conversion_factor*Upper_itp[key]),color=shade_color)
		
#		# Plot simulations
#		if typeof(species_index) == Int64
#		    dfba = time_state_array[species_index,:]
#		elseif typeof(species_index) == Array{Int64,1}
#		    dfba = sum(time_state_array[species_index,:],1)
#		else
#		    throw("Wrong type")
#		end
#		axarr_Energy[i,j][:plot](experimental_time,vec(conversion_factor*dfba),plot_color,lw=1.5)
		axarr_Energy[i,j][:fill_between](experimental_time,vec(conversion_factor*species_upper[key]),vec(species_lower[key]),color=plot_color,alpha=0.4)
	end
end
tight_layout()
savefig("Energy.pdf")
close()

return nothing
end
