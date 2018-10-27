push!(LOAD_PATH,"../../src/")
import UNSflow

AR = 6. 

alphadef = UNSflow.ConstDef(10. *pi/180)

hdef = UNSflow.ConstDef(0.)

udef = UNSflow.ConstDef(1.)

full_kinem = UNSflow.KinemDef(alphadef, hdef, udef)

pvt = 0.25

geometry = "FlatPlate"

surf = UNSflow.ThreeDSurfSimple(AR, full_kinem, geometry, pvt)

curfield = UNSflow.ThreeDFieldSimple()

dtstar = 0.015

t_tot = 8.

nsteps = Int(round(t_tot/dtstar))+1

startflag = 0

writeflag = 1

writeInterval = t_tot/8.

#delvort = delSpalart(500, 12, 1e-5)
delvort = UNSflow.delNone()

mat, surf, curfield = UNSflow.QSLLT_lautat(surf, curfield, nsteps, dtstar,startflag, writeflag, writeInterval, delvort)

UNSflow.makeVortPlots3D()

UNSflow.makeForcePlots3D()

UNSflow.makeInfoPlots3D()

UNSflow.cleanWrite()
