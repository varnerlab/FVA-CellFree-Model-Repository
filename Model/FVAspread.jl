include("CalcFVA.jl")

FVA_MAX = zeros(300,194)
FVA_MIN = zeros(300,194)
for i in 1:194;println(i)
	FVA_MAX[:,i] = readdlm("FVA/base_case/$i/FVA_max")[:,i]
	FVA_MIN[:,i] = readdlm("FVA/base_case/$i/FVA_min")[:,i]
end
FVA_spread = FVA_MIN-FVA_MAX # MAX, MIN were saved in the wrong places
FVA_norm = norm(FVA_spread)
writedlm("FVA/base_case/FVA_MAX",FVA_MAX)
writedlm("FVA/base_case/FVA_MIN",FVA_MIN)
writedlm("FVA/base_case/FVA_norm",FVA_norm)

#FVA_spread = FVA_spread[:,find(sum(FVA_spread,1).!=0)]

using PyPlot
figure()
for i in 1:194
	tmp = FVA_spread[:,i]
	plot(collect(0.01:0.01:3.0),tmp)
end

fva_max, fva_min = CalcFVA(1)

fva_max, fva_min = CalcFVA([2 3])




