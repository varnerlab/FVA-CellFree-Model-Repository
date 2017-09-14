function Unpack(path::String)

raw_data = readdlm(path)

unpacked_dictionary = Dict()
for i in 1:size(raw_data,1)
	key = raw_data[i,1]
	data = eval(parse(raw_data[i,2]))
	unpacked_dictionary[key] = data
end

return unpacked_dictionary
end

function Unpack(raw_data::Array)

if size(raw_data,2)!=2
	throw("Array is of wrong size")
end

unpacked_dictionary = Dict()
for i in 1:size(raw_data,1)
	key = raw_data[i,1]
	data = eval(parse(raw_data[i,2]))
	unpacked_dictionary[key] = data
end

return unpacked_dictionary
end
