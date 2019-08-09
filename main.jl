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

alpha = 10
alphadef = UnsteadyFlowSolvers.ConstDef(alpha *pi/180)#EldUptstartDef(amp,k,a,1)#tstart)
hdef = UnsteadyFlowSolvers.ConstDef(0)
udef = UnsteadyFlowSolvers.ConstDef(1)
full_kinem = UnsteadyFlowSolvers.KinemDef(alphadef, hdef, udef)
# Geometry
pvt = 0.25 ;
geometry = "FlatPlate"#"bin/sd7003.dat"
lespcrit = [10000;]
surf = UnsteadyFlowSolvers.TwoDSurf(geometry,pvt,full_kinem,lespcrit; ndiv = 70, camberType = "linear")
curfield = UnsteadyFlowSolvers.TwoDFlowField()
# Iteration Parameters
dtstar = UnsteadyFlowSolvers.find_tstep(alphadef)

#=Solve for time needed to acheive maximum AoA with a proportional loop
global t = tstart
global alpha = alphadef(t) *180/pi
P = 10
while alpha <= amp - .05 || alpha >= amp + .05 # Timestep is outside .05 degrees
    global alpha = alphadef(t) *180/pi
    error = amp - alpha
    global t += P*error
    #println(alpha )#* 180 / pi)
    if t > 10000
        break
    end
end
t_tot = ceil(t)
=#
t_tot = 9
nsteps = Int(round(t_tot/dtstar))
#nsteps = 3 #DEBUG

println("Running LVE")
frames, ifr_field, mat, test = UnsteadyFlowSolvers.LVE(surf,curfield,nsteps,dtstar,40,"longwake")
#circ,a,RHS = UnsteadyFlowSolvers.LVE(surf,curfield,nsteps,dtstar)

#
println("Running LDVM")
startflag = 0
writeflag = 0
writeInterval = t_tot/18.
surf2 = UnsteadyFlowSolvers.TwoDSurf(geometry,pvt,full_kinem,lespcrit; ndiv = 4)
curfield2 = UnsteadyFlowSolvers.TwoDFlowField()
delvort = UnsteadyFlowSolvers.delNone()

mat2, surf2, curfield2 = UnsteadyFlowSolvers.ldvm(surf2, curfield2, nsteps, dtstar,startflag, writeflag, writeInterval, delvort)
#mat2 = UnsteadyFlowSolvers.ldvm(surf2, curfield2, nsteps, dtstar,startflag, writeflag, writeInterval, delvort)
#
# mat = [t , alpha , h , u , LESP , cl , cd , cm] for each timestep

# Velocity plot
using Plots
pyplot()
single = "none" ;
if single == "true"
    mesh = frames[end]

    contour(mesh.x[1,:],mesh.z[:,1], mesh.velMag,fill=true,levels = 200,c = :lightrainbow_r)
    scatter!(mesh.vorX, mesh.vorZ,markersize = 1.5, markerstrokestyle = :dot, markerstrokecolor = :white)
    plot!(mesh.camX,mesh.camZ,color = :black, linewidth = 3, aspect_ratio = :equal,legend = :none,axis = false)
elseif single == "false" # Make full list of images in folder
    dir = "case_$case steps_$nsteps"
    mkdir(dir)
    for i = 1:length(frames)
        mesh = frames[i]

        contour(mesh.x[1,:],mesh.z[:,1], mesh.velMag,fill=true,levels = 200,c = :lightrainbow_r)
        scatter!(mesh.vorX, mesh.vorZ,markersize = 1.5, markerstrokestyle = :dot, markerstrokecolor = :white)
        plot!(mesh.camX,mesh.camZ,color = :black, linewidth = 3, aspect_ratio = :equal,legend = :none,axis = false)

        savefig("$dir/$i")
    end
end
#
#= Plot circ vs alpha
circ = map(q -> q.circ[1],frames)
alpha = map(q -> q.alpha,frames)
plot(alpha,circ)
circ1 = frames[1].circ
x = (1:length(circ1))/length(circ1)
plot(x,circ1)
=#

#= Force plots
#UnsteadyFlowSolvers.makeForcePlots2D()

plot(mat[:,1],mat[:,6],color = :black) # cl vs AoA
plot!(mat2[:,1],mat2[:,6],color = :red)
=#
alpha = Int(alpha)
#plot(mat2[:,1],mat2[:,2], xlabel = "t*", ylabel = "AoA")
#savefig("RameshData_Motion_$alpha")

UnsteadyFlowSolvers.subPlotForces(mat,1,"t*","MyData_$alpha")

plot(mat[:,1],test[1:end-1], xlabel = "t*", ylabel = "TEV Strength")
