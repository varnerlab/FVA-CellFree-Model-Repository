using PyPlot
using LaTeXStrings

#data = readdlm("updated_FVA/data",'\t')
#data_13dpg = readdlm("updated_FVA/combined_data_combo_13dpg")
#data = vcat(data,data_13dpg)
#data = readdlm("updated_FVA/data_corrected",'\t')
#combined_data_combo_new = readdlm("updated_FVA/combined_data_combo_new")
#combined_data_combo_13dpg = readdlm("updated_FVA/13dpg")
#data = vcat(data,combined_data_combo_new,combined_data_combo_13dpg)

data = readdlm("updated_FVA/data_final_flux_subset",'\t')

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
#num_doubles = 48
num_13dpg = 23
num_other_doubles = 25
num_deltas = 37

labels_singles = labels[2:num_singles+1]
#labels_doubles = labels[num_singles+1:num_singles+num_doubles]
labels_13dpg = labels[num_singles+2:num_singles+num_13dpg+1]
labels_other_doubles = labels[num_singles+num_13dpg+2:num_singles+num_13dpg+num_other_doubles+1]
labels_deltas = labels[num_singles+num_13dpg+num_other_doubles+2:end]

figure(figsize=(22.5,11.5))
#figure(figsize=(8,6))
plot(error[num_singles+2:num_singles+num_13dpg+1],FVA_norm[26:num_singles+num_13dpg+1],"go",markersize=m,markeredgecolor="none",alpha=a)
#plot(error[num_singles+num_13dpg+2:num_singles+num_13dpg+num_other_doubles+1],FVA_norm[num_singles+num_13dpg+2:num_singles+num_13dpg+num_other_doubles+1],"yo",markersize=m,markeredgecolor="none",alpha=a)
plot(error[num_singles+num_13dpg+num_other_doubles+2:end-1],FVA_norm[num_singles+num_13dpg+num_other_doubles+2:end-1],"ro",markersize=m,markeredgecolor="none",alpha=a)
plot(error[end],FVA_norm[end],"yo",markersize=m,markeredgecolor="none",alpha=a)
plot(error[2:num_singles+1],FVA_norm[2:num_singles+1],"bo",markersize=m,markeredgecolor="none",alpha=a)
plot(error[1],FVA_norm[1],"ko",markersize=m,markeredgecolor="none")
xlabel("error / base case error")
ylabel("flux uncertainty / base case flux uncertainty")

ax = collect(axis())
#plot(vec(ax[1:2]),vec([0 0]),"k")
#plot(vec([0 0]),vec(ax[3:4]),"k")

#key = "–GLC"
#annotate(key,xy=(error[species_indices[key]],FVA_norm[species_indices[key]]))

#key = "Base case"
#annotate(key,xy=(error[species_indices[key]]+50,FVA_norm[species_indices[key]]+10000))

#for i in 1:length(labels_singles)
#	key = labels_singles[i]
#	annotate(key,xy=(error[species_indices[key]],FVA_norm[species_indices[key]]))
#end
#for i in 1:length(labels_13dpg)
#	key = labels_13dpg[i]
#	annotate(key,xy=(error[species_indices[key]],FVA_norm[species_indices[key]]))
#end
#for i in 1:length(labels_other_doubles)
#	key = labels_other_doubles[i]
#	annotate(key,xy=(error[species_indices[key]],FVA_norm[species_indices[key]]))
#end
for i in 1:length(labels_deltas)
	key = labels_deltas[i]
	if key == "–LYS"
		annotate(key,xy=(error[species_indices[key]]-.003,FVA_norm[species_indices[key]]))
	else
		annotate(key,xy=(error[species_indices[key]],FVA_norm[species_indices[key]]))
	end
end

#ax = collect(axis())
#ax[1] = 0.9
#ax[2] = 1.3
#ax[3] = 0.6
#ax[4] = 1.8
#axis(ax)
#xticks(collect(ax[1]:0.1:ax[2]))

#savefig("updated_FVA/Main.pdf")

#for key in labels
#	if key[1:4] == "+13D"
#		annotate(key,xy=(error[species_indices[key]],FVA_norm[species_indices[key]]-.0005))
#	else
#		annotate(key,xy=(error[species_indices[key]],FVA_norm[species_indices[key]]))
#	end
#end

#ax = collect(axis())
#ax[1] = 0.175
#ax[2] = 0.24
#ax[3] = 0.295
#ax[4] = 0.39
#axis(ax)

#savefig("updated_FVA/Cluster1.pdf")

#ax = collect(axis())
#ax[1] = 0.95
#ax[2] = 1.04
#ax[3] = 0.8
#ax[4] = 1.3
#axis(ax)

#savefig("updated_FVA/Cluster2.pdf")


