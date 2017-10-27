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
  mRNA_saturation_coefficient = TXTL["mRNA_saturation_coefficient"];
  Promoter = TXTL["Promoter"]
  inducer = TXTL["inducer"]
#====================================Transcription===================================================#
#Compute the promoter strength P -
  K1 = Promoter[1]
  P = 0.8*(K1)/(1+K1);
  gene_concentration = state_array[1] * 10^6#gene_copies*(1e9/6.02e23)*(1/volume);
  saturation_term = (gene_concentration)/(plasmid_saturation_coefficient+gene_concentration);
  RNAP_concentration = RNAP_concentration_nM/1e6; #nM to mM
  TX = (RNAP_elongation_rate*(1/mRNA_length)*(RNAP_concentration)*(saturation_term)*3600)*P;

#====================================Translation===================================================#
  mRNA_steady_state = (TX/kd);
  mRNA_state = state_array[144]
  translation_rate_constant = polysome_amplification*(3*RIBOSOME_elongation_rate)*(1/mRNA_length)*3600;
  TL = translation_rate_constant*RIBOSOME_concentration*mRNA_state/(mRNA_saturation_coefficient+mRNA_state);
  DEG = mRNA_state * kd
#===================================================================================================#
  FB[167,2] = TX #transcriptional initiation
  FB[168,2] = TX #transcriptional elongation
  FB[167,1] = TX #transcriptional initiation
  FB[168,1] = TX#0 #transcriptional elongation
  FB[169,1:2] = DEG #mRNA_degradation
  FB[170,2] = TL #translations initiation
  FB[171,2] = TL #translations elongation
  FB[170,1] = 0 #translations initiation
  FB[171,1] = TL #translations elongation
  FB = max(0,FB)
  DF["default_flux_bounds_array"] = FB
  return DF
end
