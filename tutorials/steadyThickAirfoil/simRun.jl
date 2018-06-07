push!(LOAD_PATH,"../../src/")
using UNSflow

alphadef = ConstDef(0.*pi/180)

hdef = ConstDef(0.)

udef = ConstDef(1.)

kinem = KinemDef(alphadef, hdef, udef)

pvt = 0.25

geometry = "NACA0012"

surf = TwoDSurfThick(geometry, pvt, kinem)

curfield = TwoDFlowField()

dtstar = find_tstep(alphadef)

t_tot = 10.

nsteps =Int(round(t_tot/dtstar))+1

startflag = 0

writeflag = 0

writeInterval = t_tot/10.

#delvort = delSpalart(500, 12, 1e-5)
delvort = delNone()

mat, surf, curfield = lautat(surf, curfield, nsteps, dtstar,startflag, writeflag, writeInterval, delvort)

#makeVortPlots2D()

makeForcePlots()

cleanWrite()
