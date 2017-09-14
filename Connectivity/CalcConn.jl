include("SpeciesIndices.jl")

# Define function to calculate connected species
function Connectivity(species_idx)

sub_array = S[species_idx,:]
if ndims(sub_array) == 1
	sub_array = sub_array'
end
rxn_inx = find(sum(sub_array,1))
connectivity_array = S[:,rxn_inx]
connected_species = find(sum(connectivity_array,2))

return connected_species
end

# Read the stoichiometric network
S = readdlm("Network.dat")

# Make sure that TXTL stoichiometry is correct
S[1,:] = 0 # GENE_CAT
S[143:144,:] = 0 # RNAP, Ribosome
S[145,180] = 1 # mRNA (transcription)
S[145,180] = -1 # mRNA (degradation)
S[145,182] = 0 # mRNA (translation)

# Replace all non-zeros with 1
S = Array{Int64}(round(abs(S).^.00000001,1))

# Choose species/rxns to ignore
species_to_ignore = [1;39;47;71:76;101:103;110:111;125;143:144] # GENE_CAT, co2, etoh, h2o2, h2o, h2s, h_c, h_e, hco3, nh3, nh4, o2, pi, ppi, so4, RIBOSOME, RNAP
S[species_to_ignore,:] = 0

# Define number of species
num_species = size(S,1)

# Determine connected species
connected_species = Int64[]
for i in 1:num_species
	stoich_vector = S[i,:]
	rxn_inx = find(stoich_vector)
	if !isempty(rxn_inx)
		push!(connected_species,i)
	end
end

# Update array to contain connected species only
#S = S[connected_species,:]

# Calculate first-order connectivity
FirstOrderConn = Array{Int64}(zeros(num_species,1))
for i in connected_species
	tmp = Connectivity(i)
	FirstOrderConn[i] = length(tmp)-1
end

# Calculate second-order connectivity
SecondOrderConn = Array{Int64}(zeros(num_species,1))
for i in connected_species
	tmp = Connectivity(Connectivity(i))
	SecondOrderConn[i] = length(tmp)-1
end

# Calculate third-order connectivity
ThirdOrderConn = Array{Int64}(zeros(num_species,1))
for i in connected_species
	tmp = Connectivity(Connectivity(Connectivity(i)))
	ThirdOrderConn[i] = length(tmp)-1
end

# Check if third-order connectivity covers all connections
minimum(ThirdOrderConn[find(ThirdOrderConn)]) == length(FirstOrderConn[find(FirstOrderConn)])-1

# Calculate fourth-order connectivity
FourthOrderConn = Array{Int64}(zeros(num_species,1))
for i in connected_species
	tmp = Connectivity(Connectivity(Connectivity(Connectivity(i))))
	FourthOrderConn[i] = length(tmp)-1
end

# Check if fourth-order connectivity covers all connections
minimum(FourthOrderConn[find(FourthOrderConn)]) == length(FirstOrderConn[find(FirstOrderConn)])-1

# Calculate degree-of-connectivity matrix; set infinity as the default value
Degree = Array{Any}(Inf*ones(num_species,num_species))

# All species have zero-degree connectivity with themselves
for i in 1:num_species
	Degree[i,i] = 0
end

# Calculate all other connectivity indices
for i in 1:num_species
	for j in i+1:num_species
		if in(j,Connectivity(i))
			Degree[i,j] = 1
			Degree[j,i] = 1
		elseif in(j,Connectivity(Connectivity(i)))
			Degree[i,j] = 2
			Degree[j,i] = 2
		elseif in(j,Connectivity(Connectivity(Connectivity(i))))
			Degree[i,j] = 3
			Degree[j,i] = 3
		elseif in(j,Connectivity(Connectivity(Connectivity(Connectivity(i)))))
			Degree[i,j] = 4
			Degree[j,i] = 4
		end
	end
end

writedlm("Degree.dat",Degree)

type Species
	model_index::Int64
	first_order_connnectivity::Int64
	second_order_connnectivity::Int64
	third_order_connnectivity::Int64
	fourth_order_connnectivity::Int64
	connectivity_array::Array
	function Species()
		this = new()
	end
end

species_indices = SpeciesIndices()
SpeciesConnectivity = Dict()
for key in keys(species_indices)
	SpeciesConnectivity[key] = Species()
	SpeciesConnectivity[key].model_index = species_indices[key]
	SpeciesConnectivity[key].first_order_connnectivity = FirstOrderConn[SpeciesConnectivity[key].model_index]
	SpeciesConnectivity[key].second_order_connnectivity = SecondOrderConn[SpeciesConnectivity[key].model_index]
	SpeciesConnectivity[key].third_order_connnectivity = ThirdOrderConn[SpeciesConnectivity[key].model_index]
	SpeciesConnectivity[key].fourth_order_connnectivity = FourthOrderConn[SpeciesConnectivity[key].model_index]
	SpeciesConnectivity[key].connectivity_array = Degree[SpeciesConnectivity[key].model_index,:]
end













