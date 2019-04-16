push!(LOAD_PATH,"../src/")
using UnsteadyFlowSolvers
using Debugger
using Revise

alphadef = ConstDef(0. *pi/180)

hdef = ConstDef(0.)

udef = ConstDef(1.)

full_kinem = KinemDef(alphadef, hdef, udef)

pvt = 0.25

geometry = "Cylinder"

lespcrit = [10.25;]

surfvis = TwoDSurfThick(geometry, pvt, full_kinem, ndiv=140, naterm=136)

surfdup = TwoDSurfThick(geometry, pvt, full_kinem, ndiv=140, naterm=136)


curfieldvis = TwoDFlowField()

dtstar = 0.015

t_tot = 7.

nsteps = Int(round(t_tot/dtstar))+1

println("nsteps ", nsteps)

startflag = 0

writeflag = 1

writeInterval = t_tot/10.

#delvort = delSpalart(500, 12, 1e-5)
delvort = delNone()

#mat, surf, curfield = lautat(surf, curfield, nsteps, dtstar,startflag, writeflag, writeInterval, delvort)

#qu, ql = calc_edgeVel(surf, [curfield.u[1], curfield.w[1]])
istep = 0
#mat, surf, curfield, del, E, quf, qu,thick_orig, thick_orig_slope = IBLThickCoupled(surfvis, curfieldvis, 200, 1000, dtstar, startflag, writeflag, writeInterval, delvort)
includet("../src/thickCoupled/IBLThickCoupled.jl")
includet("../src/thickCoupled/IBLFV.jl")
breakpoint("../src/thickCoupled/IBLThickCoupled.jl",56)
breakpoint("../src/thickCoupled/IBLThickCoupled.jl",134)

breakpoint("../src/thickCoupled/IBLFV.jl",10)
@run IBLThickCoupled(surfvis, curfieldvis, 200, false,300, dtstar, startflag, writeflag, writeInterval, delvort)
#UNSflow.makeVortPlots2D()

#UNSflow.makeForcePlots2D()

#getEndCycle(mat, alphadef.k)

#UNSflow.cleanWrite()
