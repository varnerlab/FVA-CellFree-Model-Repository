function Bounds(DF,TXTL,state_array)

FB = DF["default_flux_bounds_array"]

  FB[195,2] = DF["Oxygen"]; #[]-> O2
  FB[207,1] = 0.95*DF["GlcUptake"]
  FB[207,2] = DF["GlcUptake"]; #[]-> GLC
  FB[211,2] = 0; #[]-> PYR

  FB[214:217,2]= 0; #succ, Mal, fum, etoh, mglx

  FB[264,2] = 0; #ATP ->[]
  FB[265,2] = 0; #[]-> ADP

  #METABOLITES
  FB[210,1] = DF["PYR"]
  FB[210,2] = 1.2*DF["PYR"]
  FB[212,1] = DF["AC"]
  FB[212,2] = 1.2*DF["AC"]
  FB[213,1] = DF["LAC"]
  FB[213,2] = 1.2*DF["LAC"]
  FB[214,1] = DF["SUCC"]
  FB[214,2] = 1.2*DF["SUCC"]
  FB[215,1] = DF["MAL"]
  FB[215,2] = 1.2*DF["MAL"]
  FB[216:218,2]= 0; # fum, etoh, mglx
  #Amino Acids
  FB[224:2:262,2] = DF["AAUptake"]; #AA UPTAKE
  FB[225:2:263,2] = DF["AASecretion"]; #AA ->[]
  FB[224,2] = 0; #ALA
  FB[225,1] = 0.1*DF["ALA"]
  FB[225,2] = DF["ALA"]
  FB[230,2] = DF["ASP"]
  FB[231,2] = 0;
  FB[236,2] = 0;
  FB[237,2] = DF["GLN"]
  FB[246,2] = DF["LysUptake"]
  FB[247,2] = DF["LysSecretion"]

  #ENERGY
  FB[264:287,2] = DF["ENERGY"];
  FB[283,2] = 0;

  AASyn = TXTL["AA_synthesis_rxn"]
  AADeg = TXTL["AA_degradation_rxn"]
  FB[AASyn,2] = DF["AASyn"];
  FB[AADeg,2] = 0;

  FB[1:191,1] = 0
  FB[1:191,2] = 100000
  FB[192:end,1] = 0
  FB[192:end,2] = 0
  FB[197,2] = 100 # hydrogen transport

#==============================================TXTL=====================================================#
  RNAP_concentration_nM = state_array[143] * 10^6 # conversion to nM #TXTL["RNAP_concentration_nM"];
  RNAP_elongation_rate = TXTL["RNAP_elongation_rate"];
  RIBOSOME_concentration = state_array[141]#TXTL["RIBOSOME_concentration"];
  RIBOSOME_elongation_rate = TXTL["RIBOSOME_elongation_rate"];
  kd = TXTL["mRNA_degradation_rate"]*2;
  mRNA_length = TXTL["mRNA_length"];
  protein_length = TXTL["protein_length"];
  gene_copies = TXTL["gene_copies"];
  volume = TXTL["volume"];
  polysome_amplification = TXTL["polysome_gain"];
  plasmid_saturation_coefficient = TXTL["plasmid_saturation_coefficient"];
  mRNA_saturation_coefficient = TXTL["mRNA_saturation_coefficient"]
  Promoter = TXTL["Promoter"]
  inducer = TXTL["inducer"]
  ala_tRNA_saturation_coefficient = .001
  arg_tRNA_saturation_coefficient = .001
  asn_tRNA_saturation_coefficient = .0001
  asp_tRNA_saturation_coefficient = .001
  cys_tRNA_saturation_coefficient = .001
  gln_tRNA_saturation_coefficient = .001
  glu_tRNA_saturation_coefficient = .001
  gly_tRNA_saturation_coefficient = .001
  his_tRNA_saturation_coefficient = .001
  ile_tRNA_saturation_coefficient = .001
  leu_tRNA_saturation_coefficient = .001
  lys_tRNA_saturation_coefficient = .001
  met_tRNA_saturation_coefficient = .001
  phe_tRNA_saturation_coefficient = .001
  pro_tRNA_saturation_coefficient = .0001
  ser_tRNA_saturation_coefficient = .001
  thr_tRNA_saturation_coefficient = .001
  trp_tRNA_saturation_coefficient = .0001
  tyr_tRNA_saturation_coefficient = .0001
  val_tRNA_saturation_coefficient = .001
  
  tRNA_coeff = [15; 5; 10; 12; 5; 12; 13; 10; 12; 9; 13; 12; 9; 20; 7; 10; 13; 5; 11; 16]
#====================================Transcription===================================================#
#Compute the promoter strength P -
  K1 = Promoter[1]
  P = 0.8*(K1)/(1+K1);
  gene_concentration = state_array[1] * 10^6#gene_copies*(1e9/6.02e23)*(1/volume);
  saturation_term = (gene_concentration)/(plasmid_saturation_coefficient+gene_concentration);
  RNAP_concentration = RNAP_concentration_nM/1e6; #nM to mM
  TX = (RNAP_elongation_rate*(1/mRNA_length)*(RNAP_concentration)*(saturation_term)*3600)*P;
#====================================Translation===================================================#
  mRNA_state = state_array[144]
  
  ala_tRNA_state = state_array[23]
  arg_tRNA_state = state_array[26]
  asn_tRNA_state = state_array[28]
  asp_tRNA_state = state_array[30]
  cys_tRNA_state = state_array[43]
  gln_tRNA_state = state_array[62]
  glu_tRNA_state = state_array[64]
  gly_tRNA_state = state_array[67]
  his_tRNA_state = state_array[77]
  ile_tRNA_state = state_array[80]
  leu_tRNA_state = state_array[85]
  lys_tRNA_state = state_array[87]
  met_tRNA_state = state_array[90]
  phe_tRNA_state = state_array[107]
  pro_tRNA_state = state_array[111]
  ser_tRNA_state = state_array[122]
  thr_tRNA_state = state_array[127]
  trp_tRNA_state = state_array[129]
  tyr_tRNA_state = state_array[131]
  val_tRNA_state = state_array[136]
  
  translation_rate_constant = polysome_amplification*(3*RIBOSOME_elongation_rate)*(1/mRNA_length)*3600;
  TL = translation_rate_constant*RIBOSOME_concentration*mRNA_state/(mRNA_saturation_coefficient+mRNA_state)#*(ala_tRNA_state/(ala_tRNA_saturation_coefficient+ala_tRNA_state)*arg_tRNA_state/(arg_tRNA_saturation_coefficient+arg_tRNA_state)*asn_tRNA_state/(asn_tRNA_saturation_coefficient+asn_tRNA_state)*asp_tRNA_state/(asp_tRNA_saturation_coefficient+asp_tRNA_state)*cys_tRNA_state/(cys_tRNA_saturation_coefficient+cys_tRNA_state)*gln_tRNA_state/(gln_tRNA_saturation_coefficient+gln_tRNA_state)*glu_tRNA_state/(glu_tRNA_saturation_coefficient+glu_tRNA_state)*gly_tRNA_state/(gly_tRNA_saturation_coefficient+gly_tRNA_state)*his_tRNA_state/(his_tRNA_saturation_coefficient+his_tRNA_state)*ile_tRNA_state/(ile_tRNA_saturation_coefficient+ile_tRNA_state)*leu_tRNA_state/(leu_tRNA_saturation_coefficient+leu_tRNA_state)*lys_tRNA_state/(lys_tRNA_saturation_coefficient+lys_tRNA_state)*met_tRNA_state/(met_tRNA_saturation_coefficient+met_tRNA_state)*phe_tRNA_state/(phe_tRNA_saturation_coefficient+phe_tRNA_state)*pro_tRNA_state/(pro_tRNA_saturation_coefficient+pro_tRNA_state)*ser_tRNA_state/(ser_tRNA_saturation_coefficient+ser_tRNA_state)*thr_tRNA_state/(thr_tRNA_saturation_coefficient+thr_tRNA_state)*trp_tRNA_state/(trp_tRNA_saturation_coefficient+trp_tRNA_state)*tyr_tRNA_state/(tyr_tRNA_saturation_coefficient+tyr_tRNA_state)*val_tRNA_state/(val_tRNA_saturation_coefficient+val_tRNA_state))
  DEG = mRNA_state * kd
#===================================================================================================#
  FB[167,1] = TX #transcriptional initiation
  FB[168,1] = TX #0 #transcriptional elongation
  FB[167,2] = TX #transcriptional initiation
  FB[168,2] = TX #transcriptional elongation
  FB[169,1] = 0 #mRNA_degradation
  FB[169,2] = DEG #mRNA_degradation
  FB[170,1] = 0 #translation initiation
  FB[171,1] = 0 #translation elongation
  FB[170,2] = TL #translation initiation
  FB[171,2] = TL #translation elongation
  
#  charging_lb = (state_array[140]<.005)?.5:.5
  charging_lb = state_array[140]*20
  
  FB[172:191,1] = 0#TL/5 #tRNA charging
  FB[172:191,2] = 50*TL #tRNA charging
  FB = max(0,FB)
  DF["default_flux_bounds_array"] = FB
  return DF
end
