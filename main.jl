include("src/UnsteadyFlowSolvers.jl")
## Define Geometry and Kinematics
# Kinematics
case = 6 # kinematics case to run
if case == 1 # ~ 8 hours to run, Ramshes @ 4 hrs
    amp = 35 #degrees
    k = .005
    tstart = 60 # non-dimensional time
    t_tot = 277 # 138.525 for up only, 277 for full
elseif case == 2
    amp = 45
    k = .05
    tstart = 20
    t_tot = 60.15 # 29.5 for up only, 60.15 for full
elseif case == 3
    amp = 90
    k = .4
    tstart = 10
    t_tot = 25 # 12.5 for up only, 25 for full
elseif case == 4
    amp = 25
    k = .11
    a = 11
    tstart = 1
    t_tot = 7.09 # 3.55 for up only, 7.09 for full
elseif case == 5
    amp = 45
    k = 0.4
    a = 11
    tstart = 1
    t_tot = 4.52 # 2.26 for up only, 4.52 for full
elseif case == 6
    mean = 0
    amp = 30
    k = 0.1
    phi = 0
    t_tot = 31.4 #roughly
elseif case == 7
    mean = 0
    amp = 30
    k = 0.4
    phi = 0
    t_tot = 7.85
end
#amp = amp * pi / 180
if case < 4
a = (pi^2 * k*180 )/(2*amp*pi *(1-0.1))
end

#alpha = 10
if case < 6
    alphadef = UnsteadyFlowSolvers.EldRampReturntstartDef(amp*pi/180,k,a,tstart) # ConstDef(alpha*pi/180)#
else
    alphadef = UnsteadyFlowSolvers.SinDef(mean,amp*pi/180,k,phi)
end
hdef = UnsteadyFlowSolvers.ConstDef(0)
udef = UnsteadyFlowSolvers.ConstDef(1)
full_kinem = UnsteadyFlowSolvers.KinemDef(alphadef, hdef, udef)
# Geometry
pvt = 0.25 ;
geometry = "FlatPlate"#"bin/sd7003.dat"
lespcrit = [.2;]

surf = UnsteadyFlowSolvers.TwoDSurf(geometry,pvt,full_kinem,lespcrit; ndiv = 101, camberType = "linear", rho = 1.225)
curfield = UnsteadyFlowSolvers.TwoDFlowField()
# Iteration Parameters
dtstar = .015 #UnsteadyFlowSolvers.find_tstep(alphadef)

nsteps = Int(round(t_tot/dtstar))
#nsteps = 1 #DEBUG

println("Running LVE")
frames, IFR_field, mat, test = UnsteadyFlowSolvers.LVE(surf,curfield,nsteps,dtstar,40,"longwake")

# Velocity plot
using Plots
pyplot()
single = "false" ;
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
