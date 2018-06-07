include("../../src/UNSflow.jl")
using UNSflow

include("types.jl")
include("calcs.jl")
        
alphadef = EldRampReturnDef(25,0.11,11)
hdef = ConstDef(0.)
udef = ConstDef(1.)
full_kinem = KinemDef(alphadef, hdef, udef)

pvt = 0.0 #leading edge

surf = TwoDSurfPanel("sd7003_fine.dat", pvt, full_kinem)

curfield = TwoDFlowField()
