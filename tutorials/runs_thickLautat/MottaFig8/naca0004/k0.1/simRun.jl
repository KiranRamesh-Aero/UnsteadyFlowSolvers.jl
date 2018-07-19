push!(LOAD_PATH,"../../../UNSflow/src/")
using UNSflow

k = 0.1
amp = 1.*pi/180

alphadef = SinDef(0., amp, k, 0.)

hdef = ConstDef(0.)

udef = ConstDef(1.)

kinem = KinemDef(alphadef, hdef, udef)

pvt = 0.25

geometry = "NACA0004"

surf = TwoDSurfThick(geometry, pvt, kinem)

curfield = TwoDFlowField()

dtstar = 0.15#find_tstep(alphadef)

T = pi/k
t_tot = 5.*T

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
