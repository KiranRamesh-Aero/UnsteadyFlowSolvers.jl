include("src/UnsteadyFlowSolvers.jl")

## Define Geometry and Kinematics
# Kinematics
alphadef = UnsteadyFlowSolvers.SinDef(0,5,.4,0)
hdef = UnsteadyFlowSolvers.ConstDef(0)
udef = UnsteadyFlowSolvers.ConstDef(1)
full_kinem = UnsteadyFlowSolvers.KinemDef(alphadef, hdef, udef)
# Geometry
pvt = 0.25 ;
geometry = "bin/airfoil.dat"
lespcrit = [0.11;]
surf = UnsteadyFlowSolvers.TwoDSurf(geometry,pvt,full_kinem,lespcrit)
curfield = UnsteadyFlowSolvers.TwoDFlowField()
# Iteration Parameters
dtstar = UnsteadyFlowSolvers.find_tstep(alphadef)
t_tot = 9
nsteps = Int(round(t_tot/dtstar))+1

newsurf, test = UnsteadyFlowSolvers.LVE(surf,curfield,nsteps,dtstar)
