using PyPlot

Accuracy = readdlm("Accuracy")
Precision = readdlm("Precision")
a = 1./Accuracy
p = 1./Precision

EXIT_FLAG = 5-readdlm("EXIT_FLAG")
ef = sum(EXIT_FLAG,1)/300

Labels = readdlm("species_constraint_labels")

tmp = sortrows([ef' a p Labels])
ef = tmp[:,1]
a = tmp[:,2]
p = tmp[:,3]
Labels = tmp[:,4]

figure()
plot(collect(1:length(ef)),ef)
ylabel("Exit flag")

figure()
plot(collect(1:length(Accuracy)),a)
ylabel("1/Accuracy")

figure()
plot(collect(1:length(Precision)),p)
ylabel("1/Precision")









