push!(LOAD_PATH,"../../src/")
import UNSflow


alphadef = UNSflow.ConstDef(10. *pi/180)

hdef = UNSflow.ConstDef(0.)

udef = UNSflow.ConstDef(1.)

full_kinem = UNSflow.KinemDef(alphadef, hdef, udef)

pvt = 0.0

geometry = "FlatPlate"

lespc = [10.15;]

#The value of leading-edge radius (rho) must be specified according to geometry
surf = UNSflow.TwoDSurf(geometry, pvt, full_kinem, lespc, rho=0.02085)

curfield = UNSflow.TwoDFlowField()

dtstar = 0.015 #UNSflow.find_tstep(alphadef)

t_tot = 5.

nsteps =Int(round(t_tot/dtstar))+1

startflag = 0

writeflag = 1

writeInterval = t_tot/10.

#delvort = delSpalart(500, 12, 1e-5)
delvort = UNSflow.delNone()

mat, surf, curfield = UNSflow.ldvmLin(surf, curfield, nsteps, dtstar,startflag, writeflag, writeInterval, delvort)

#UNSflow.makeVortPlots2D()

#UNSflow.makeForcePlots2D()

#UNSflow.cleanWrite()
