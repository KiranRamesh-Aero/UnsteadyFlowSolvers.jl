push!(LOAD_PATH,"../../src/")
import UNSflow

alphadef = UNSflow.CosDef(0., 1. *pi/180, 0.05, 0.)

hdef = UNSflow.ConstDef(0.)

udef = UNSflow.ConstDef(1.)

full_kinem = UNSflow.KinemDef(alphadef, hdef, udef)

pvt = 0.25

#geometry = "FlatPlate0821"
geometry = "NACA0018"

lespcrit = [0.1;]

surf = UNSflow.TwoDSurfThick(geometry, pvt, full_kinem, lespcrit, ndiv=140, naterm=137)

curfield = UNSflow.TwoDFlowField()

dtstar = UNSflow.find_tstep(alphadef)/5. 

t_tot = 5*pi/0.1

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
