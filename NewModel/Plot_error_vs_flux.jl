num_dir = Int64(readdlm("Ens/num_dir")[1])

L = zeros(146,0)
for i in 1:num_dir
	if isdir("Ens/$i")
		L = [L readdlm("Ens/$i/species_boolean")]
	end
end

Error = Float64[]
for i in 1:num_dir
	if isdir("Ens/$i")
		push!(Error,readdlm("Ens/$i/state_error")[1])
	end
end

Flux = Float64[]
for i in 1:num_dir
	if isdir("Ens/$i")
		push!(Flux,readdlm("Ens/$i/flux_uncertainty")[1])
	end
end

m = 8
a = .3

using PyPlot

r = collect(0:1/316:1)
g = collect(0:1/316:1)
b = collect(0:1/316:1)

figure(figsize=(12,8))
for i in 1:length(Error)
	semilogx(Error[i],Flux[i],color=vec([r[i] g[i] b[i]]),"o",markersize=m,markeredgecolor="black",alpha=a)
end
xlabel("error")
ylabel("flux uncertainty")

j = 0
for i in 1:num_dir;println(i)
	if isdir("Ens/$i")
		j += 1
		cp("Ens/$i","Ens_j/$j")
	end
end




