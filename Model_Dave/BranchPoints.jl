include("CalcBranch.jl")

glc_uptake = 1
glycolysis = [2;3]
PP = [23;24]
PP_ru5p = 26
ED = [37;38]

glycolysis_fraction = CalcBranch(glycolysis,glc_uptake)
PP_fraction = CalcBranch(PP,glc_uptake)

PP_to_ru5p_fraction = CalcBranch(PP_ru5p,PP)
PP_to_ED_fraction = CalcBranch(ED,PP)












