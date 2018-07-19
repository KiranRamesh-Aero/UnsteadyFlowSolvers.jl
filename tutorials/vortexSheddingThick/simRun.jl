push!(LOAD_PATH,"../../src/")
using UNSflow

alphadef = ConstDef(45.*pi/180)

hdef = ConstDef(0.)

udef = ConstDef(1.)

kinem = KinemDef(alphadef, hdef, udef)

pvt = 0.25

geometry = "FlatPlate0417"

lespcrit = [0.1;]

surf = TwoDSurfThick(geometry, pvt, kinem, lespcrit)

curfield = TwoDFlowField()

dtstar = 0.015

t_tot = 10.

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
