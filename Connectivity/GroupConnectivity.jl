include("SpeciesIndices.jl")

function GroupConnectivity(idx::Array{Int64})

sort!(idx)
Degree = readdlm("Degree.dat")
sub_array = Degree[idx,idx]
group_connectivity_index = maximum(sub_array)
if group_connectivity_index != Inf
	group_connectivity_index = Int64(group_connectivity_index)
end

return group_connectivity_index
end

function GroupConnectivity(species_list::Array{String})

species_indices = SpeciesIndices()

idx = Int64[]
for key in species_list
	push!(idx,species_indices[key])
end

sort!(idx)
Degree = readdlm("Degree.dat")
sub_array = Degree[idx,idx]
group_connectivity_index = maximum(sub_array)
if group_connectivity_index != Inf
	group_connectivity_index = Int64(group_connectivity_index)
end

return group_connectivity_index
end
