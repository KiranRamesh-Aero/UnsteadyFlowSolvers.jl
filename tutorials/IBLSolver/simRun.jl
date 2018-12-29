push!(LOAD_PATH,"../../src")
import UNSflow

alphadef = UNSflow.ConstDef(5. *pi/180)

hdef = UNSflow.ConstDef(0.)

udef = UNSflow.ConstDef(1.)

full_kinem = UNSflow.KinemDef(alphadef, hdef, udef)

pvt = 0.25

geometry = "FlatPlate"

lespc = [10.15;]

numDiv =500
#The value of leading-edge radius (rho) must be specified c according to geometry
surf = UNSflow.TwoDSurf(geometry, pvt, full_kinem, lespc, ndiv=numDiv,rho=0.02085)

curfield = UNSflow.TwoDFlowField()

dtstar = UNSflow.find_tstep(alphadef)

t_tot = 3.

nsteps =Int(round(t_tot/dtstar))+1

startflag = 0

writeflag = 1

writeInterval = t_tot/10.

#delvort = delSpalart(500, 12, 1e-5)
delvort = UNSflow.delNone()

# define number of timesteps for the boundary layer problem.
# the current implementation assumes numDiv=blDiv,

blDiv=numDiv

# define the viscous time-step value.
# setting the value to viscousTimeStep =0.0 implies time-step calculation based on CFL value.
# automatic time-steppping is recommedned for CFL<0.2

#viscousTimeStep =0.0004
viscousTimeStep =0.0

# setting up the operations conditions for the boundary layer problem opCond=[cfl,Reynolnds number]
opCond =[0.2, 10000.0]

mat, surf, curfield, soln, fluxSplit ,invis = UNSflow.invisicViscousCoupledSolver(surf, curfield, blDiv, opCond[1], opCond[2], userdefinedDt=viscousTimeStep)

#mat, surf, curfield = UNSflow.ldvmLin(surf, curfield, nsteps, dtstar,startflag, writeflag, writeInterval, delvort)
#q_u, q_l = UNSflow.calc_edgeVel(surf, [curfield.u[1], curfield.w[1]])


#UNSflow.makeVortPlots2D()

#UNSflow.makeForcePlots2D()

#UNSflow.cleanWrite()

#Edge velocity at end of simulation

#using PyPlot
#display(plot(surf.x, q_u))
#display(plot(surf.x, q_l))
