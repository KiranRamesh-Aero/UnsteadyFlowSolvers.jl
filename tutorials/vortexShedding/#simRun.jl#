push!(LOAD_PATH,"../../src/")
import UNSflow


alphadef = UNSflow.EldUpDef(45. *pi/180, 0.4, 0.8)

hdef = UNSflow.ConstDef(0.)

udef = UNSflow.ConstDef(1.)

full_kinem = UNSflow.KinemDef(alphadef, hdef, udef)

pvt = 0.0

geometry = "FlatPlate"

lespcrit = [0.11;]

surf = UNSflow.TwoDSurf(geometry, pvt, full_kinem, lespcrit)

curfield = UNSflow.TwoDFlowField()

dtstar = UNSflow.find_tstep(alphadef)

t_tot = 9.

nsteps =Int(round(t_tot/dtstar))+1

startflag = 0

writeflag = 1

writeInterval = t_tot/18.

#delvort = delSpalart(500, 12, 1e-5)
delvort = UNSflow.delNone()

mat, surf, curfield = UNSflow.ldvm(surf, curfield, nsteps, dtstar,startflag, writeflag, writeInterval, delvort)

UNSflow.makeVortPlots2D()

UNSflow.makeForcePlots()

UNSflow.cleanWrite()
