tic()

include("Include.jl")

if !isdir("Obj")
	mkdir("Obj")
end

data_dictionary = DataDictionary(0,0,0)
TXTL_dictionary = TXTLDictionary(0,0,0)
number_of_fluxes = length(data_dictionary["default_flux_bounds_array"][:,1])
number_of_species = length(data_dictionary["species_bounds_array"][:,1])
data_dictionary["objective_coefficient_array"][194] = -1

data_dictionary["species_drift_amplitude"] = 2 # percent allowed to drift for unconstraint species for each iteration (0.01 hr)
data_dictionary["species_drift_decay_rate"] = 1

if isfile("Obj/objective_coefficient_array_best")
	objective_coefficient_array_best = readdlm("Obj/objective_coefficient_array_best")
else
	objective_coefficient_array_best = copy(data_dictionary["objective_coefficient_array"])
	objective_coefficient_array_best = max(-1,objective_coefficient_array_best-.001)
	writedlm("Obj/objective_coefficient_array_best",objective_coefficient_array_best)
end

if isfile("Obj/objective_coefficient_array")
	objective_coefficient_array = readdlm("Obj/objective_coefficient_array")
else
	objective_coefficient_array = copy(objective_coefficient_array_best)
	writedlm("Obj/objective_coefficient_array",objective_coefficient_array)
end

if isfile("Obj/cost_best")
	cost_best = readdlm("Obj/cost_best")[1]
else
	error_vector, Exit_flag = dFBA(objective_coefficient_array)
	failed_tps = sum(5 - Exit_flag)
	println("Number of failed tps: ",failed_tps)
	cost_best = sum(error_vector)
	writedlm("Obj/cost_best",cost_best)
end

if isfile("Obj/cost")
	cost = readdlm("Obj/cost")[1]
else
	cost = copy(cost_best)
	writedlm("Obj/cost",cost)
end

alpha = 25
variance_min = .5
variance_max = 50

for i in 1:10
variance = exp(log(variance_min)+rand()*(log(variance_max)-log(variance_min)))
num_to_vary = randperm(length(objective_coefficient_array))[1]
indices_to_vary = sort(randperm(length(objective_coefficient_array))[1:num_to_vary])

# Create perturbation vector
perturb = ones(size(objective_coefficient_array))
tmp = exp(variance*randn(length(objective_coefficient_array)))
perturb[indices_to_vary] = tmp[indices_to_vary]

# Perturb objective coefficient array
objective_coefficient_array_new = copy(objective_coefficient_array)
objective_coefficient_array_new = objective_coefficient_array.*perturb
objective_coefficient_array_new = max(-1,objective_coefficient_array_new)
objective_coefficient_array_new = min(1,objective_coefficient_array_new)

# Run dFBA
error_vector_new, Exit_flag = dFBA(objective_coefficient_array_new)

failed_tps = sum(5 - Exit_flag)
println("Number of failed tps: ",failed_tps)
cost_new = sum(error_vector_new)

acc_prob = exp(alpha*(cost-cost_new)/cost)

cost_round = round(cost,2)
cost_new_round = round(cost_new,2)
cost_best_round = round(cost_best,2)
acc_prob_round = round(acc_prob,2)
var_round = round(variance,2)

# If a new best is achieved, overwrite parameters and cost and save the relevant information to Best directory
if cost_new < cost_best
    # Save to Best directory
    objective_coefficient_array_best = copy(objective_coefficient_array_new)
    cost_best = copy(cost_new)
    writedlm("Obj/objective_coefficient_array",objective_coefficient_array_best)
    writedlm("Obj/cost_best",cost_new)
    println("cost_new = $cost_new_round, cost = $cost_round, best = $cost_best_round, acc_prob = $acc_prob_round, var = $var_round, num_to_vary = $num_to_vary, NEW BEST")
else
    println("cost_new = $cost_new_round, cost = $cost_round, best = $cost_best_round, acc_prob = $acc_prob_round, var = $var_round, num_to_vary = $num_to_vary")
end

# If new cost is lower than previous cost, choose new params as reference point for perturbation
if rand() < acc_prob # Accept all better sets and some worse sets
    # Create directory for next sample and save the relevant information
#            mkdir("Ensemble/$i")
#            writedlm("Ensemble/$i/rate_constant.dat",rate)
#            writedlm("Ensemble/$i/saturation_constant.dat",sat)
#            writedlm("Ensemble/$i/control_constant.dat",cont)
#            writedlm("Ensemble/$i/Tsim",Tsim)
#            writedlm("Ensemble/$i/X",X)
#            writedlm("Ensemble/$i/cost",cost_new)
    # Overwrite variables
    objective_coefficient_array = copy(objective_coefficient_array_new)
    cost = copy(cost_new)
    writedlm("Obj/cost",cost_new)
end
end
time_elapsed = toc()

