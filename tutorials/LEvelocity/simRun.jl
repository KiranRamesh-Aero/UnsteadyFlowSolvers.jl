push!(LOAD_PATH,"../../src/")
import UNSflow


alphadef = UNSflow.ConstDef(5. *pi/180)

hdef = UNSflow.ConstDef(0.)

udef = UNSflow.ConstDef(1.)

full_kinem = UNSflow.KinemDef(alphadef, hdef, udef)

pvt = 0.25

geometry = "FlatPlate"

lespc = [0.15;]

#The value of leading-edge radius (rho) must be specified according to geometry
surf = UNSflow.TwoDSurf(geometry, pvt, full_kinem, lespc, rho=0.016)

curfield = UNSflow.TwoDFlowField()

dtstar = UNSflow.find_tstep(alphadef)

t_tot = 3.

nsteps =Int(round(t_tot/dtstar))+1

startflag = 0

writeflag = 1

writeInterval = t_tot/10.

#delvort = delSpalart(500, 12, 1e-5)
delvort = UNSflow.delNone()

mat, surf, curfield = UNSflow.ldvmLin(surf, curfield, nsteps, dtstar,startflag, writeflag, writeInterval, delvort)

UNSflow.makeVortPlots2D()

UNSflow.makeForcePlots2D()

UNSflow.cleanWrite()

#Edge velocity at end of simulation
q_u, q_l = UNSflow.calc_edgeVel(surf, [curfield.u[1], curfield.w[1]])

using PyPlot
plot(surf.x, q_u)
plot(surf.x, q_l)
