include("src/UnsteadyFlowSolvers.jl")
using DelimitedFiles

dir0 = "Full Run"
mkdir(dir0)
for case = 2:7 # For each case
    dir1 = "$dir0/Case $case"
    mkdir(dir1)
    println("Case $case")
    for levToggle = 1:2 # for LEV on and off
        if levToggle == 1 # LEV on
            dir2 = "$dir1/LEV on"
            lespcrit = [.2;]
        else # LEV off
            dir2 = "$dir1/LEV off"
            lespcrit = [10000;]
        end
        mkdir(dir2)
        for pvt = [0,.25,1] # for each pivot point
            if pvt == 0
                pvtType = "LE"
            elseif pvt == .25
                pvtType = "QC"
            else
                pvtType = "TE"
            end
            dir3 = "$dir2/$pvtType"
            #mkdir(dir3)

            ## Define Geometry and Kinematics
            # Kinematics
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
                tstart = 0
            end
            hdef = UnsteadyFlowSolvers.ConstDef(0)
            udef = UnsteadyFlowSolvers.ConstDef(1)
            full_kinem = UnsteadyFlowSolvers.KinemDef(alphadef, hdef, udef)
            # Geometry
            geometry = "FlatPlate"#"bin/sd7003.dat"
            global surf = UnsteadyFlowSolvers.TwoDSurf(geometry,pvt,full_kinem,lespcrit; ndiv = 101, camberType = "linear", rho = 1.225)
            global curfield = UnsteadyFlowSolvers.TwoDFlowField()
            # Iteration Parameters
            dtstar = .015 #UnsteadyFlowSolvers.find_tstep(alphadef)

            nsteps = Int(round(t_tot/dtstar))
            #nsteps = 150 #DEBUG

            println("Running LVE")
            global frames, IFR_field, mat, test1, test2 = UnsteadyFlowSolvers.LVE(surf,curfield,nsteps,dtstar,40,"longwake")
            writedlm("$dir3/myMat.txt",mat)

            #=
            println("Running LDVM")
            startflag = 0
            writeflag = 0
            writeInterval = t_tot/18.
            surf2 = UnsteadyFlowSolvers.TwoDSurf(geometry,pvt,full_kinem,lespcrit; ndiv = 101,rho = 1.225)
            curfield2 = UnsteadyFlowSolvers.TwoDFlowField()
            delvort = UnsteadyFlowSolvers.delNone()

            global mat2, surf2, curfield2 = UnsteadyFlowSolvers.ldvm(surf2, curfield2, nsteps, dtstar,startflag, writeflag, writeInterval, delvort)
            # mat = [t , alpha , h , u , LESP , cl , cd , cm] for each timestep
            #writedlm("$dir3/ramMat.txt",mat2)
            =#
        end
    end
end
