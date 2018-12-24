using UNSflow

using Test

alphadef = UNSflow.ConstDef(5. *pi/180)
hdef = UNSflow.ConstDef(0.)
udef = UNSflow.ConstDef(1.)
full_kinem = UNSflow.KinemDef(alphadef, hdef, udef)
pvt = 0.25
geometry = "FlatPlate"
lespc = [10.15;]

#The value of leading-edge radius (rho) must be specified according to geometry
surf = UNSflow.TwoDSurf(geometry, pvt, full_kinem, lespc, rho=0.016)

@test typeof(surf)==TwoDSurf
