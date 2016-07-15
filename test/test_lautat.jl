# workspace()
# include("UNSflow.jl")
# using UNSflow

alphadef = EldUpDef(45,0.2,0.8)
hdef = ConstDef(0.)
udef = ConstDef(1.)
full_kinem = KinemDef(alphadef, hdef, udef)

lesp_crit = [50;]

surf = TwoDSurf(1., 1., "sd7003_fine.dat", 0.35, 70, 35, "Prescribed", full_kinem, lesp_crit)
#surf = TwoDSurf(1., 1., "FlatPlate", 0.35, 70, 35, "Prescribed", full_kinem)

curfield = TwoDFlowField()

lautat(surf, curfield)

