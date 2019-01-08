push!(LOAD_PATH,"../../src/")
import UNSflow

AR = 6. 

alphadef = UNSflow.ConstDef(5. *pi/180)

hdef = UNSflow.SinDef(0., 0.5, 1.0, 0.)

udef = UNSflow.ConstDef(1.)

full_kinem = UNSflow.KinemDef(alphadef, hdef, udef)

pvt = 0.25

geometry = "FlatPlate"

lespcrit = [10.1;]

surf = UNSflow.ThreeDSurfSimple(AR, full_kinem, geometry, pvt, lespcrit)

curfield = UNSflow.ThreeDFieldSimple()

dtstar = 0.015

t_tot = 5. *pi/hdef.k

nsteps = Int(round(t_tot/dtstar))+1

println("nsteps ", nsteps)

startflag = 0

writeflag = 1

writeInterval = t_tot/16.

#delvort = UNSflow.delSpalart(500, 12, 1e-5)
delvort = UNSflow.delNone()

mat, surf, curfield = UNSflow.QSLLT_lautat(surf, curfield, nsteps, dtstar,startflag, writeflag, writeInterval, delvort)

#UNSflow.makeVortPlots3D()

#UNSflow.makeForcePlots3D()

UNSflow.makeInfoPlots3D()

UNSflow.cleanWrite()
