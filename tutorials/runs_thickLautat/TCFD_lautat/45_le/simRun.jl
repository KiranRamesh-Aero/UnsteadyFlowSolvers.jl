push!(LOAD_PATH,"../../../UNSflow/src/")
using UNSflow

k = 0.11
amp = 25.*pi/180
a = 11.

alphadef = EldRampReturnDef(amp, k, a)

hdef = ConstDef(0.)

udef = ConstDef(1.)

kinem = KinemDef(alphadef, hdef, udef)

pvt = 0.0

geometry = "FlatPlate0417"

surf = TwoDSurfThick(geometry, pvt, kinem)

curfield = TwoDFlowField()

dtstar = 0.015#find_tstep(alphadef)

t_tot = 7.

nsteps =Int(round(t_tot/dtstar))+1

startflag = 0

writeflag = 1

writeInterval = t_tot/10.

#delvort = delSpalart(500, 12, 1e-5)
delvort = delNone()

mat, surf, curfield = lautat(surf, curfield, nsteps, dtstar,startflag, writeflag, writeInterval, delvort)

makeVortPlots2D()

makeForcePlots()

cleanWrite()
