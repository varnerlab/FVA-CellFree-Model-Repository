include("Include.jl")
data_dictionary, TXTL_dictionary = LoadDictionaries()

Upper, Lower = LoadSynthData()

plot_color = "orangered"
PlotErrorbar(data_dictionary, Upper, Lower, plot_color)



