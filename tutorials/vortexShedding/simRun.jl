push!(LOAD_PATH,"../../src/")
using UnsteadyFlowSolvers
using PyPlot

alphadef = ConstDef(5. * pi/180)
#alphadef = SinDef(0. *pi/180, 10. *pi/180, 0.2, 0.)

hdef = ConstDef(0.)
#hdef = SinDef(0. ,0.5 , 0.1, 0.)

udef = ConstDef(1.)

#betadef = ConstDef(22. *pi/180)
betadef = CosDef(0. *pi/180, 20. *pi/180, 0.2, 0.)

full_kinem = KinemDefFlap(alphadef, hdef, udef, betadef)

pvt = 0.0

hinge = 0.8

geometry = "FlatPlate"

lespcrit = [0.11;]

surf = TwoDSurfFlap(geometry, pvt, hinge, full_kinem, lespcrit)

curfield = TwoDFlowField()

#dtstar = find_tstep(betadef)
dtstar = 0.015

t_tot = 2. *pi/betadef.k
#t_tot = 10.

nsteps = Int(round(t_tot/dtstar))+1
println("nsteps ", nsteps)

startflag = 0

writeflag = 1

writeInterval = t_tot/10.

#delvort = delSpalart(500, 12, 1e-5)
delvort = delNone()

mat, surf, curfield = ldvmFlap(surf, curfield, nsteps, dtstar,startflag, writeflag, writeInterval, delvort)

makeVortPlots2D()

makeForcePlots2DFlap()

#plot(surf.x_star, surf.cam_slope)

cleanWrite()
