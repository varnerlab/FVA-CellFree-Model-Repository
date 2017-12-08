using PyPlot
using LaTeXStrings

#data = readdlm("updated_FVA/data",'\t')
#data_13dpg = readdlm("updated_FVA/combined_data_combo_13dpg")
#data = vcat(data,data_13dpg)

data = readdlm("updated_FVA/data_corrected",'\t')

#s7p = find(data[:,1].=="+S7P")[1]
#data = data[[1:s7p-1;s7p+1:size(data,1)],:]

labels = data[:,1]
error = data[:,2]
FVA_norm = data[:,3]

species_indices = Dict()
for i in 1:length(labels)
	species_indices[labels[i]] = i
end

error /= error[1]
FVA_norm /= FVA_norm[1]

m = 5
a = .3

num_singles = 24
num_doubles = 15
num_deltas = 37

figure(figsize=(5,3.5))
plot(error[1],FVA_norm[1],"ko",markersize=m,markeredgecolor="none")
plot(error[2:num_singles+1],FVA_norm[2:num_singles+1],"go",markersize=m,markeredgecolor="none",alpha=a)
plot(error[num_singles+2:num_singles+num_doubles+1],FVA_norm[26:num_singles+num_doubles+1],"ro",markersize=m,markeredgecolor="none",alpha=a)
plot(error[num_singles+num_doubles+2:end],FVA_norm[num_singles+num_doubles+2:end],"bo",markersize=m,markeredgecolor="none",alpha=a)
xlabel("error")
ylabel("FVA_norm")

ax = collect(axis())
#plot(vec(ax[1:2]),vec([0 0]),"k")
#plot(vec([0 0]),vec(ax[3:4]),"k")

key = "–GLC"
annotate(key,xy=(error[species_indices[key]],FVA_norm[species_indices[key]]))

key = "Base case"
annotate(key,xy=(error[species_indices[key]]+50,FVA_norm[species_indices[key]]+10000))

#for i in 1:length(labels)
#	key = labels[i]
#	annotate(key,xy=(error[species_indices[key]],FVA_norm[species_indices[key]]))
#end

#savefig("updated_FVA/Main.pdf")

for key in labels
	if key[1:4] == "+13D"
		annotate(key,xy=(error[species_indices[key]],FVA_norm[species_indices[key]]-.0005))
	else
		annotate(key,xy=(error[species_indices[key]],FVA_norm[species_indices[key]]))
	end
end

ax = collect(axis())
ax[3] = 0
axis(ax)

#savefig("updated_FVA/Cluster1.pdf")

#ax[1] = -14200
#ax[2] = -13600
#ax[3] = -1690000
#ax[4] = -1650000
#axis(ax)

##key = "+13DPG"
##annotate(key,xy=(error[species_indices[key]],FVA_norm[species_indices[key]]-700))

##key = "+DHAP+FUM"
##annotate(key,xy=(error[species_indices[key]],FVA_norm[species_indices[key]]+200))

##for key in ["+DHAP+PEP"; "+DHAP+6PGC"; "+DHAP+XU5P"; "+DHAP+CIT"]
##	annotate(key,xy=(error[species_indices[key]],FVA_norm[species_indices[key]]+200))
##end

#savefig("updated_FVA/Cluster2.pdf")









