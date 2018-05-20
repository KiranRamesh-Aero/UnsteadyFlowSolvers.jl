push!(LOAD_PATH,"../../src/")
using UNSflow

alphadef = ConstDef(5.*pi/180)

hdef = ConstDef(0.)

udef = ConstDef(1.)

kinem = KinemDef3D(alphadef, hdef, udef)

AR = 10.

pvt = 0.25

geometry = "FlatPlate"

surf = ThreeDSurfSimple(AR, kinem, geometry, pvt)

field = ThreeDFieldSimple()

dtstar = find_tstep(alphadef)

t_tot = 10.0

nsteps =Int(round(t_tot/dtstar))+1

startflag = 0

writeflag = 1

writeInterval = t_tot/12.

#delvort = delSpalart(500, 12, 1e-5)
delvort = delNone()

mat, surf, field = QSLLTlautatRoll(surf, field, nsteps, dtstar, startflag,
writeflag, writeInterval, delvort)

#maxwrite = 100; nround=6

makeForcePlots3Dstrip()

makeVortPlots3Dstrip()

makeTevstrPlots3Dstrip()

cleanWrite()
