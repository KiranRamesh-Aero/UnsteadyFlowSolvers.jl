push!(LOAD_PATH,"../../src/")
using UnsteadyFlowSolvers

#alphadef = ConstDef(0. * pi/180)
alphadef = CosDef(0. *pi/180, 30. *pi/180, 0.1, 0.)

hdef = ConstDef(0.)

udef = ConstDef(1.)

betadef = ConstDef(60. *pi/180)
#betadef = SinDef(0. *pi/180, 10. *pi/180, 0.1, 0.)

full_kinem = KinemDefFlap(alphadef, hdef, udef, betadef)

pvt = 0.0

hinge = 0.9

geometry = "FlatPlate"

lespcrit = [0.11;]

surf = TwoDSurfFlap(geometry, pvt, hinge, full_kinem, lespcrit)

curfield = TwoDFlowField()

dtstar = find_tstep(betadef)

t_tot = 40.

nsteps = Int(round(t_tot/dtstar))+1

startflag = 0

writeflag = 1

writeInterval = t_tot/80.

#delvort = delSpalart(500, 12, 1e-5)
delvort = delNone()

mat, surf, curfield = ldvmFlap(surf, curfield, nsteps, dtstar,startflag, writeflag, writeInterval, delvort)

makeVortPlots2D()

makeForcePlots2DFlap()

cleanWrite()
