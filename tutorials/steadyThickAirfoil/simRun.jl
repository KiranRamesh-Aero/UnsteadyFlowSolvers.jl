push!(LOAD_PATH,"../../src/")
import UNSflow

alphadef = UNSflow.ConstDef(10. *pi/180)

hdef = UNSflow.ConstDef(0.)

udef = UNSflow.ConstDef(1.)

full_kinem = UNSflow.KinemDef(alphadef, hdef, udef)

pvt = 0.25

geometry = "NACA0012"

surf = UNSflow.TwoDSurfThick(geometry, pvt, full_kinem, ndiv=140, naterm=138)

curfield = UNSflow.TwoDFlowField()

dtstar = 0.015#UNSflow.find_tstep(alphadef)

t_tot = 5.

nsteps =Int(round(t_tot/dtstar))+1

startflag = 0

writeflag = 1

writeInterval = t_tot/10.

#delvort = delSpalart(500, 12, 1e-5)
delvort = UNSflow.delNone()

mat, surf, curfield = UNSflow.lautat(surf, curfield, nsteps, dtstar,startflag, writeflag, writeInterval, delvort)

#UNSflow.makeVortPlots2D()

#UNSflow.makeForcePlots2D()

#UNSflow.makeInfoPlots2D()

UNSflow.cleanWrite()
