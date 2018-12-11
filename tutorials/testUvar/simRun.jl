push!(LOAD_PATH,"../../src/")
import UNSflow


alphadef = UNSflow.ConstDef(6. *pi/180)

hdef = UNSflow.ConstDef(0.)

udef = UNSflow.SinDef(0., 0.5, 0.25, 0.)

full_kinem = UNSflow.KinemDef(alphadef, hdef, udef)

pvt = 0.25

geometry = "FlatPlate"

lespcrit = [0.2;]

surf = UNSflow.TwoDSurf(geometry, pvt, full_kinem, lespcrit)

curfield = UNSflow.TwoDFlowField(UNSflow.ConstDef(1.), UNSflow.ConstDef(0.))

dtstar = 0.05

t_tot = 60.

nsteps = Int(round(t_tot/dtstar))+1

startflag = 0

writeflag = 1

writeInterval = t_tot/5.

#delvort = delSpalart(500, 12, 1e-5)
delvort = UNSflow.delNone()

mat, surf, curfield = UNSflow.ldvmLin(surf, curfield, nsteps, dtstar,startflag, writeflag, writeInterval, delvort)

UNSflow.makeVortPlots2D()

UNSflow.makeForcePlots2D()

UNSflow.cleanWrite()
