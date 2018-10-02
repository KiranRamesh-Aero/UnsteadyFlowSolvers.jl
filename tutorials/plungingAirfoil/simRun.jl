push!(LOAD_PATH,"../../src/")
using UNSflow

alphadef = ConstDef(4. *pi/180)

hdef = SinDef(0., 0.05, 3.93, 0.)

udef = ConstDef(1.)

full_kinem = KinemDef(alphadef, hdef, udef)

pvt = 0.25

geometry = "sd7003.dat"

surf = TwoDSurf(geometry, pvt, full_kinem)

curfield = TwoDFlowField()

dtstar = find_tstep(hdef)

t_tot = 5. *pi/hdef.k

nsteps =Int(round(t_tot/dtstar))+1

startflag = 0

writeflag = 1

writeInterval = t_tot/20.

#delvort = delSpalart(500, 12, 1e-5)
delvort = delNone()

mat, surf, curfield = lautatRoll(surf, curfield, nsteps, dtstar,startflag, writeflag, writeInterval, delvort)

makeVortPlots2D()

makeForcePlots()

cleanWrite()
