include("src/UnsteadyFlowSolvers.jl")

## Define Geometry and Kinematics
# Kinematics
alphadef = UnsteadyFlowSolvers.SinDef(0,20,.4,0)
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
t_tot = 2
nsteps = Int(round(t_tot/dtstar))+1

newsurf, frames, newfield, glo_field, test = UnsteadyFlowSolvers.LVE(surf,curfield,nsteps,dtstar)
#aniStep = UnsteadyFlowSolvers.LVE(surf,curfield,nsteps,dtstar)

using Plots
pyplot()
mesh = frames[end]
x = mesh.x[1,:]
y = mesh.z[:,1]
contour(x,y, mesh.velMag,fill=true,levels = 200,c = :lightrainbow_r)
#scatter!(map(q -> q.x, glo_field.tev), map(q -> q.z, glo_field.tev), color = :black)
plot!(mesh.camX,mesh.camZ,color = :black, linewidth = 3, aspect_ratio = :equal,legend = :none,axis = false)
