function Plot(plot_color)

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

figure("main",figsize=(23.5,12.5))
species = ["M_glc_D_c"; "PROTEIN_CAT"; "M_pyr_c"; "M_lac_D_c"; "M_ac_c"; "M_succ_c"; "M_mal_L_c"; "M_atp_c"; "M_adp_c"; "M_amp_c"; "M_gtp_c"; "M_gdp_c"; "M_gmp_c"; "M_ctp_c"; "M_cdp_c"; "M_cmp_c"; "M_utp_c"; "M_udp_c"; "M_ump_c"; "M_ala_L_c"; "M_arg_L_c" ; "M_asn_L_c"; "M_asp_L_c"; "M_cys_L_c"; "M_gln_L_c"; "M_glu_L_c"; "M_gly_L_c"; "M_his_L_c"; "M_ile_L_c"; "M_leu_L_c"; "M_lys_L_c"; "M_met_L_c"; "M_phe_L_c"; "M_pro_L_c"; "M_ser_L_c"; "M_thr_L_c"; "M_trp_L_c"; "M_tyr_L_c"; "M_val_L_c"]
subplot_index = [1; 2; 3; 4; 9; 10; 11; 17; 25; 33; 18; 26; 34; 19; 27; 35; 20; 28; 36; 5; 6; 7; 8; 13; 14; 15; 16; 21; 22; 23; 24; 29; 30; 31; 32; 37; 38; 39; 40]
for i in 1:length(species)
	key = species[i]
	subplot(5,8,subplot_index[i])
	species_number = species_indices[key]
	shade_color = (species_number in species_constraint_index)?"lightblue":"lightgray"
	fill_between(vec(t_array_synthetic),vec(Lower_itp[key]),vec(Upper_itp[key]),color=shade_color)
	plot(experimental_time,time_state_array[species_number,:],plot_color,lw=1.5)
	title(key)
	a = collect(axis())
	a[3] = 0
	axis(a)
end

#subplot(5,8,12)
#plot(experimental_time,vec(sum(time_state_array[tRNA_index,:],1)),plot_color,lw=1.5)
#a = collect(axis())
#a[3] = 0
#axis(a)

tight_layout()

figure("TXTL",figsize=(23.5,12.5))
species = ["M_nad_c"; "M_nadh_c"; "M_nadp_c"; "M_nadph_c"; "M_q8_c"; "M_q8h2_c"; "M_mql8_c"; "M_mqn8_c"; "M_h_c"; "M_h_e"; "M_pi_c"; "M_ppi_c"; "GENE_CAT"; "RNAP"; "OPEN_GENE_CAT"; "mRNA_CAT"; "tRNA"; "RIBOSOME"; "RIBOSOME_START_CAT"; "PROTEIN_CAT"; "M_ala_L_c_tRNA"; "M_arg_L_c_tRNA"; "M_asn_L_c_tRNA"; "M_asp_L_c_tRNA"; "M_cys_L_c_tRNA"; "M_gln_L_c_tRNA"; "M_glu_L_c_tRNA"; "M_gly_L_c_tRNA"; "M_his_L_c_tRNA"; "M_ile_L_c_tRNA"; "M_leu_L_c_tRNA"; "M_lys_L_c_tRNA"; "M_met_L_c_tRNA"; "M_phe_L_c_tRNA"; "M_pro_L_c_tRNA"; "M_ser_L_c_tRNA"; "M_thr_L_c_tRNA"; "M_trp_L_c_tRNA"; "M_tyr_L_c_tRNA"; "M_val_L_c_tRNA"]
subplot_index = [1; 2; 3; 4; 9; 10; 11; 12; 17; 18; 19; 20; 25; 26; 27; 28; 33; 34; 35; 36; 5; 6; 7; 8; 13; 14; 15; 16; 21; 22; 23; 24; 29; 30; 31; 32; 37; 38; 39; 40]
for i in 1:length(species)
	key = species[i]
	subplot(5,8,subplot_index[i])
	species_number = species_indices[key]
	shade_color = (species_number in species_constraint_index)?"lightblue":"lightgray"
	fill_between(vec(t_array_synthetic),vec(Lower_itp[key]),vec(Upper_itp[key]),color=shade_color)
	plot(experimental_time,time_state_array[species_number,:],plot_color,lw=1.5)
	title(key)
	a = collect(axis())
	a[3] = 0
	axis(a)
end
tight_layout()

figure("Glycolysis_PPP_TCA",figsize=(23.5,12.5))
species = ["M_glc_D_c"; "M_g6p_c"; "M_f6p_c"; "M_fdp_c"; "M_dhap_c"; "M_g3p_c"; "M_13dpg_c"; "M_3pg_c"; "M_2pg_c"; "M_pep_c"; "M_pyr_c"; "M_6pgl_c"; "M_6pgc_c"; "M_ru5p_D_c"; "M_r5p_c"; "M_xu5p_D_c"; "M_s7p_c"; "M_e4p_c"; "M_2ddg6p_c"; "M_accoa_c"; "M_cit_c"; "M_icit_c"; "M_akg_c"; "M_succoa_c"; "M_succ_c"; "M_fum_c"; "M_mal_L_c"; "M_oaa_c"]
subplot_index = collect(1:28)
for i in 1:length(species)
	key = species[i]
	subplot(4,7,subplot_index[i])
	species_number = species_indices[key]
	shade_color = (species_number in species_constraint_index)?"lightblue":"lightgray"
	fill_between(vec(t_array_synthetic),vec(Lower_itp[key]),vec(Upper_itp[key]),color=shade_color)
	plot(experimental_time,time_state_array[species_number,:],plot_color,lw=1.5)
	title(key)
	a = collect(axis())
	a[3] = 0
	axis(a)
end
tight_layout()

return nothing
end
