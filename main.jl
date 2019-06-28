include("src/UnsteadyFlowSolvers.jl")

## Define Geometry and Kinematics
# Kinematics
case = 3 # kinematics case to run
if case == 1
    amp = 35 #degrees
    k = .005
    tstart = 60 # non-dimensional time
elseif case == 2
    amp = 45
    k = .05
    tstart = 20
elseif case == 3
    amp = 90
    k = .4
    tstart = 10
elseif case == 4
    amp = 25
    k = .11
    a = 11
    tstart = 1
end
#amp = amp * pi / 180
if case != 4
a = (pi^2 * k*180 )/(2*amp*pi *(1-0.1))
end

alphadef = UnsteadyFlowSolvers.EldUptstartDef(amp,k,a,0)#tstart)
hdef = UnsteadyFlowSolvers.ConstDef(0)
udef = UnsteadyFlowSolvers.ConstDef(1)
full_kinem = UnsteadyFlowSolvers.KinemDef(alphadef, hdef, udef)
# Geometry
pvt = 0.25 ;
geometry = "FlatPlate"#bin/airfoil.dat"
lespcrit = [0.11;]
surf = UnsteadyFlowSolvers.TwoDSurf(geometry,pvt,full_kinem,lespcrit)
curfield = UnsteadyFlowSolvers.TwoDFlowField()
# Iteration Parameters
dtstar = UnsteadyFlowSolvers.find_tstep(alphadef)
#=Solve for time needed to acheive maximum AoA
global t = 0
global alpha = alphadef(t)
P = .5

while alpha <= amp || alpha >= amp + 1 # Timestep is within 1 degree
    global alpha = alphadef(t)
    error = amp - alpha
    global t += P*error
    println(alpha )#* 180 / pi)
    if t > 100
        break
    end
end
=#
t_tot = 10
nsteps = Int(round(t_tot/dtstar))+1

newsurf, frames, newfield, glo_field, mat, test = UnsteadyFlowSolvers.LVE(surf,curfield,nsteps,dtstar,40)
#u,w,n = UnsteadyFlowSolvers.LVE(surf,curfield,nsteps,dtstar)

# Velocity plot
using Plots
pyplot()
mesh = frames[end]
x = mesh.x[1,:] ;
y = mesh.z[:,1] ;
contour(x,y, mesh.velMag,fill=true,levels = 200,c = :lightrainbow_r)
#scatter!(map(q -> q.x, glo_field.tev), map(q -> q.z, glo_field.tev), color = :black)
plot!(mesh.camX,mesh.camZ,color = :black, linewidth = 3, aspect_ratio = :equal,legend = :none,axis = false)

# Force plots
#UnsteadyFlowSolvers.makeForcePlots2D()
