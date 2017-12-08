using PyPlot
percent_appearance = ["GLC" 89.9 "substrate"; "G6P" 89.9 "substrate"; "CIT" 85.7 "TCA"; "ICIT" 84.9 "TCA"; "FUM" 84.5 "TCA"; "F6P" 83.6 "glycolysis"; "6PGL" 81.5 "PP"; "S7P" 79.4 "PP"; "ALA" 77.7 "amino"; "GTP" 77.3 "energy"; "MAL" 75.6 "TCA"; "RU5P" 74.8 "PP"; "E4P" 73.1 "PP"; "ADP" 72.3 "energy"; "AKG" 71.8 "TCA"; "UDP" 70.2 "energy"; "SUCC" 69.7 "TCA"; "CMP" 69.3 "energy"; "GDP" 67.2 "energy"; "6PGC" 67.2 "PP"; "ARG" 66.8 "amino"; "R5P" 66.4 "PP"; "MET" 65.5 "amino"; "GLX" 65.5 "TCA"; "GLN" 63.9 "amino"; "PHE" 63.4 "amino"; "VAL" 62.6 "amino"; "G3P" 62.2 "glycolysis"; "AMP" 62.2 "energy"; "PRO" 60.1 "amino"; "FDP" 60.1 "glycolysis"; "DHAP" 60.1 "glycolysis"; "HIS" 59.7 "amino"; "GLY" 59.7 "amino"; "OAA" 58.4 "TCA"; "2DDG6P" 58.4 "PP"; "CAT" 56.3 "CAT"; "CYS" 55.5 "amino"; "AC" 54.2 "TCA"; "SUCCOA" 53.8 "TCA"; "UMP" 52.9 "energy"; "TRP" 52.5 "amino"; "LAC" 52.5 "glycolysis"; "UTP" 52.1 "energy"; "ASP" 51.7 "amino"; "GMP" 50.8 "energy"; "ASN" 50.4 "amino"; "CDP" 50.0 "energy"; "PEP" 48.3 "glycolysis"; "3PG" 47.9 "glycolysis"; "2PG" 47.1 "glycolysis"; "LYS" 43.3 "amino"; "THR" 42.0 "amino"; "GLU" 38.7 "amino"; "TYR" 37.0 "amino"; "ATP" 37.0 "energy"; "XU5P" 35.7 "PP"; "ACCOA" 34.9 "TCA"; "CTP" 33.2 "energy"; "SER" 31.9 "amino"; "ILE" 29.8 "amino"; "LEU" 26.5 "amino"; "PYR" 21.8 "glycolysis"]

# Plot the bar graph with fixed width of the bar
width = 1
a = [1:4:4*size(percent_appearance,1);]
threshold = 0.05*ones(a)
PyPlot.matplotlib[:rc]("text",usetex=true)
fig,ax = plt[:subplots](figsize=(15,7))

color_dict = Dict()
color_dict["substrate"] = "blue"
color_dict["glycolysis"] = "green"
color_dict["PP"] = "red"
color_dict["TCA"] = "yellow"
color_dict["energy"] = "orange"
color_dict["amino"] = "purple"
color_dict["CAT"] = "pink"

color_array = percent_appearance[:,3]
for i in 1:length(color_array)
	color_array[i] = color_dict[color_array[i]]
end

# Set up the format for presenting each dataset
for i in 1:size(percent_appearance,1)
	b = ax[:bar](a[i]+(i-1)*width, percent_appearance[i,2], width, align="center", color=color_array[i]) #, ecolor="black")
end

#b_control = ax[:bar](a,Sensitivity[:,1],yerr=Error[:,1],width,align="center",color="black",ecolor="black")
#b_woAA = ax[:bar](a+width,Sensitivity[:,2],yerr=Error[:,2],width,align="center",color="#858585 ",ecolor="black")
#b_woSyn = ax[:bar](a+2*width,Sensitivity[:,3],yerr=Error[:,3],width,align="center",color="lightgrey",ecolor="black")

## Set up the characteristic of the figure
#ax[:set_xticks]([])
#plt[:margins](0.002)
#plt[:ylim]([0,.8])
#plt[:yticks]([0:.2:.8],fontsize=20)
#plt[:ylabel]("Sensitivity Index",fontsize=20)
#plt[:tick_params](axis="x",which="both",bottom="off",top="off")

##plt[:xlabel]("GLC OX RNAP rate RIBO rate ALA ARG ASN ASP CYS GLU GLN GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL",fontsize=14.5)
#x = vec(["GLC" "OX" "RNAP" "rate" "RIBO" "rate" "ALA" "ARG" "ASN" "ASP" "CYS" "GLU" "GLN" "GLY" "HIS" "ILE" "LEU" "LYS" "MET" "PHE" "PRO" "SER" "THR" "TRP" "TYR" "VAL"]')

#a = [1:4:4*length(x);]
#ax[:set_xticks](a)
#plt[:tick_params](axis="x",which="both",bottom="off",top="off")

#ax[:set_xticklabels](x,rotation=70)

