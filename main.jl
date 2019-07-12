include("src/UnsteadyFlowSolvers.jl")
using Revise
## Define Geometry and Kinematics
# Kinematics
case = 2 # kinematics case to run
if case == 1
    amp = 35 #degrees
    k = .005
    tstart = 60 # non-dimensional time
elseif case == 2
    amp = 45
    k = .4#05
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
a = .1#(pi^2 * k*180 )/(2*amp*pi *(1-0.1))
end

alphadef = UnsteadyFlowSolvers.EldUptstartDef(amp,k,a,1)#tstart)
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

#Solve for time needed to acheive maximum AoA with a proportional loop
global t = tstart
global alpha = alphadef(t)
P = 10
while alpha <= amp - .05 || alpha >= amp + .05 # Timestep is outside .05 degrees
    global alpha = alphadef(t)
    error = amp - alpha
    global t += P*error
    #println(alpha )#* 180 / pi)
    if t > 10000
        break
    end
end
t_tot = ceil(t)
resol = 50 # Step resolution from calculation (x/100%)
nsteps = round(t_tot/dtstar)+1
nsteps = Int(ceil(nsteps * resol / 100))
dtstar = dtstar * 100 / resol
#Calculate time to run
#nsteps = 20000 #DEBUG

if nsteps >= 1900 && nsteps <= 23000
    T = 2*10^-6 * nsteps^2 - .0164 * nsteps + 31.279 # in minutes
    hour = Int(floor(T/60))
    min = Int(floor(T%60))
    if hour == 0
        println("Approximate time to run is $min minutes.")
    else
        println("Approximate time to run is $hour hours, $min minutes.")
    end
end

println("Running LVE")
frames, glo_field, mat, test = UnsteadyFlowSolvers.LVE(surf,curfield,nsteps,dtstar,40)
#frames = UnsteadyFlowSolvers.LVE(surf,curfield,nsteps,dtstar)

#=
println("Running LDVM")
startflag = 0
writeflag = 1
writeInterval = t_tot/18.
delvort = UnsteadyFlowSolvers.delNone()

mat2, surf2, curfield2 = UnsteadyFlowSolvers.ldvm(surf, curfield, nsteps, dtstar,startflag, writeflag, writeInterval, delvort)
=#
# mat = [t , alpha , h , u , LESP , cl , cd , cm] for each timestep

# Velocity plot
using Plots
pyplot()
mesh = frames[end]

contour(mesh.x[1,:],mesh.z[:,1], mesh.velMag,fill=true,levels = 200,c = :lightrainbow_r)
scatter!(mesh.vorX, mesh.vorZ,markersize = 1.5, markerstrokestyle = :dot, markerstrokecolor = :white)
plot!(mesh.camX,mesh.camZ,color = :black, linewidth = 3, aspect_ratio = :equal,legend = :none,axis = false)

#= Plot circ vs alpha
circ = map(q -> q.circ[1],frames)
alpha = map(q -> q.alpha,frames)
plot(alpha,circ)
circ1 = frames[1].circ
x = (1:length(circ1))/length(circ1)
plot(x,circ1)


# Force plots
#UnsteadyFlowSolvers.makeForcePlots2D()

plot(mat[:,2],mat[:,6]) # cl vs AoA
=#
