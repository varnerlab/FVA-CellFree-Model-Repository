include("Include.jl")

data_dictionary, TXTL_dictionary, species_list, Mean, Upper, Lower = LoadDictionaries()

time_state_array = readdlm("Ens/68/time_state_array")

plot_color = "green"

#PlotErrorbar(data_dictionary, Upper, Lower, plot_color)
PlotMain(plot_color)

