using PyPlot

param_bounds = [8e-4 12e-4
				0.0015 0.003
				20 30
				116e-7 116e-5
				5 15
				1.5 3
				0.0045 0.45
				8 12]

num_dir = Int64(readdlm("TXTL/num_dir")[1])
Params = zeros(8,num_dir)
for i in 1:num_dir
	Params[:,i] = readdlm("TXTL/$i/params")
end
mean_std = [mean(Params,2) std(Params,2)]
max_min = [minimum(Params,2) maximum(Params,2)]

#[minimum(Params,2) mean(Params,2) maximum(Params,2)]

max_min_interp = zeros(size(max_min))
for i in 1:size(param_bounds,1)
	tmp = max_min[i,:]
	bounds = param_bounds[i,:]
	tmp_interp = tmp-bounds[1]
	tmp_interp = tmp_interp/(bounds[2]-bounds[1])
	max_min_interp[i,:] = tmp_interp
end

figure()
for i in 1:size(max_min_interp,1)
	fill_between(vec([i-.5 i]),vec([max_min_interp[i,1] max_min_interp[i,1]]),vec([max_min_interp[i,2] max_min_interp[i,2]]))
end

