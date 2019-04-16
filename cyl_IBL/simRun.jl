push!(LOAD_PATH,"../src/")
using UnsteadyFlowSolvers


alphadef = ConstDef(0. *pi/180)

hdef = ConstDef(0.)

udef = ConstDef(1.)

full_kinem = KinemDef(alphadef, hdef, udef)

pvt = 0.25

geometry = "Cylinder"

lespcrit = [10.25;]

surf = TwoDSurfThick(geometry, pvt, full_kinem, ndiv=140, naterm=136)

curfield = TwoDFlowField()

dtstar = 0.005

t_tot = 1.

nsteps = Int(round(t_tot/dtstar))+1

println("nsteps ", nsteps)

startflag = 0

writeflag = 1

writeInterval = t_tot/10.

#delvort = delSpalart(500, 12, 1e-5)
delvort = delNone()

mat, surf, curfield, del, E, quf, qu,thick_orig, thick_orig_slope = transpCoupled(surf, curfield, 200,  nsteps, dtstar, startflag, writeflag, writeInterval, delvort)


