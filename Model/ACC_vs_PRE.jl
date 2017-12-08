using PyPlot
using LaTeXStrings

#base_case = ["Base case" readdlm("FVA/base_case/error")[1] readdlm("FVA/base_case/fva_norm")[1]]
#single_FVA_data = readdlm("updated_FVA/combined_data_singles")
#combo_FVA_data = readdlm("updated_FVA/combined_data_combo")
#KO_FVA_data = readdlm("updated_FVA/combined_data_KO")
#data = vcat(base_case, single_FVA_data, combo_FVA_data, KO_FVA_data)
data = readdlm("updated_FVA/data",'\t')

s7p = find(data[:,1].=="+S7P")[1]
data = data[[1:s7p-1;s7p+1:size(data,1)],:]

labels = data[:,1]
error = data[:,2]
FVA_norm = data[:,3]

species_indices = Dict()
for i in 1:length(labels)
	species_indices[labels[i]] = i
end

error -= error[1]
FVA_norm -= FVA_norm[1]

#error *= -1
#FVA_norm *= -1

#single_FVA_data = readdlm("single_FVA_data")
##error = single_FVA_data[:,2]
##FVA_norm = single_FVA_data[:,3]
#labels = single_FVA_data[[1:15;17:25],1] # Ignore 13dpg and s7p
#error = single_FVA_data[[1:15;17:25],2]
#FVA_norm = single_FVA_data[[1:15;17:25],3]

#combo_FVA_data = readdlm("combo_FVA_data")
#labels = [labels;combo_FVA_data[:,1]]
#error = [error;combo_FVA_data[:,2]]
#FVA_norm = [FVA_norm;combo_FVA_data[:,3]]

#error_itp = error-minimum(error)
#error_itp /= maximum(error_itp)
#error_itp *= .8
#error_itp += .1
#FVA_norm_itp = FVA_norm-minimum(FVA_norm)
#FVA_norm_itp /= maximum(FVA_norm_itp)
#FVA_norm_itp *= .8
#FVA_norm_itp += .1

m = 8
a = .3

figure(figsize=(23.5,12.5))
#plot(error_itp,FVA_norm_itp,"wo")
#plot(error_itp[1:24],FVA_norm_itp[1:24],"go")
#plot(error_itp[25:end],FVA_norm_itp[25:end],"ro")
plot(error[1],FVA_norm[1],"ko",markersize=m,markeredgecolor="none")
plot(error[2:25],FVA_norm[2:25],"go",markersize=m,markeredgecolor="none",alpha=a)
plot(error[26:40],FVA_norm[26:40],"ro",markersize=m,markeredgecolor="none",alpha=a)
plot(error[41:end],FVA_norm[41:end],"bo",markersize=m,markeredgecolor="none",alpha=a)
#ax = collect(axis())
#ax[1] = 0
#ax[3] = 0
#axis(ax)
xlabel("error")
ylabel("FVA_norm")

ax = collect(axis())
plot(vec(ax[1:2]),vec([0 0]),"k")
plot(vec([0 0]),vec(ax[3:4]),"k")

key = "–GLC"
annotate(key,xy=(error[species_indices[key]],FVA_norm[species_indices[key]]))

key = "Base case"
annotate(key,xy=(error[species_indices[key]]+50,FVA_norm[species_indices[key]]+10000))

#for i in 1:length(labels)
#	key = labels[i]
#	annotate(key,xy=(error[species_indices[key]],FVA_norm[species_indices[key]]))
#end




