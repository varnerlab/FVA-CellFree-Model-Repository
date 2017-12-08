include("DataDictionary.jl")
include("LoadDictionaries.jl")
include("LoadSynthData.jl")
include("FluxDriver.jl")
include("calculate_constraints.jl")
include("Bounds.jl")
include("TXTLDictionary.jl")
include("Utility.jl")
include("CalcError.jl")
include("PlotMain.jl")
include("PlotErrorbar.jl")
include("CalcErrorbar.jl")
using PyPlot
using Interpolations
using LaTeXStrings
using GLPK
term_out(GLPK.OFF)
