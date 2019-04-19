function transpCoupled(surf::TwoDSurfThickBL, curfield::TwoDFlowField, ncell::Int64, nsteps::Int64 = 300, dtstar::Float64 = 0.015, startflag = 0, writeflag = 0, writeInterval = 1000., delvort = delNone(); maxwrite = 50, nround=6)

    # If a restart directory is provided, read in the simulation data
    if startflag == 0
        mat = zeros(0, 12)
        t = 0.

    elseif startflag == 1
        dirvec = readdir()
        dirresults = map(x->(v = tryparse(Float64,x); typeof(v) == Nothing ? 0.0 : v),dirvec)
        latestTime = maximum(dirresults)
        mat = DelimitedFiles.readdlm("resultsSummary")
        t = mat[end,1]
    else
        throw("invalid start flag, should be 0 or 1")
    end
    mat = mat'

    dt = dtstar*surf.c/surf.uref

    # if writeflag is on, determine the timesteps to write at
    if writeflag == 1
        writeArray = Int64[]
        tTot = nsteps*dt
        for i = 1:maxwrite
            tcur = writeInterval*real(i)
            if t > tTot
                break
            else
                push!(writeArray, Int(round(tcur/dt)))
            end
        end
    end

    vcore = 0.02*surf.c

    int_wax = zeros(surf.ndiv)
    int_c = zeros(surf.ndiv)
    int_t = zeros(surf.ndiv)

    for istep = 1:nsteps

        t = t + dt

        #Update kinematic parameters
        update_kinem(surf, t)

        #Update flow field parameters if any
        update_externalvel(curfield, t)

        #Update bound vortex positions
        update_boundpos(surf, dt)

        #Update incduced velocities on airfoil
        update_indbound(surf, curfield)

        #Set up the matrix problem
        surf, xloc_tev, zloc_tev = update_thickLHS(surf, curfield, dt, vcore)

        #Construct RHS vector
        update_thickRHS(surf, curfield)

        #Solve inviscid problem
        invsoln = surf.LHS[1:surf.ndiv*2-2, 1:surf.naterm*2+1] \ surf.RHS[1:surf.ndiv*2-2]



        #Assign the solution
        for i = 1:surf.naterm
            surf.aterm[i] = invsoln[i]
            surf.bterm[i] = invsoln[i+surf.naterm]
        end
        tevstr = invsoln[2*surf.naterm+1]*surf.uref*surf.c
        push!(curfield.tev, TwoDVort(xloc_tev, zloc_tev, tevstr, vcore, 0., 0.))

        avisc = zeros(surf.naterm)
        bvisc = zeros(surf.naterm)

        #Update induced velocities to include effect of last shed vortex
        update_indbound(surf, curfield)

        res = 1.

        qux = zeros(surf.ndiv)
        qlx = zeros(surf.ndiv)
        qut = zeros(surf.ndiv)
        qlt = zeros(surf.ndiv)

        iter_delu = zeros(surf.ndiv)
        iter_dell = zeros(surf.ndiv)
        iter_Eu = zeros(surf.ndiv)
        iter_El = zeros(surf.ndiv)
        wtu = zeros(surf.ndiv)
        wtl = zeros(surf.ndiv)

        #Iterate for viscous solution and interaction

        while res > 1e-6

            surf.aterm[:] = invsoln[1:surf.naterm] .+ avisc[:]
            surf.bterm[:] = invsoln[surf.naterm+1:2*surf.naterm] .+ bvisc[:]

            surf.qu[:], surf.ql[:], qux[:], qlx[:] = calc_edgeVel(surf, [curfield.u[1], curfield.w[1]])

            if istep == 1
                surf.quprev[:] = surf.qu[:]
                surf.qlprev[:] = surf.ql[:]
            end

            qut[:] = (surf.qu[:] .- surf.quprev[:])./dt
            qlt[:] = (surf.ql[:] .- surf.qlprev[:])./dt

            w0 = [surf.delu surf.delu.*(surf.Eu .+ 1)]
            w, j1 ,j2 = FVMIBL(w0, surf.qu, qut, qux, surf.x, dt)
            iter_delu = w[:,1]
            iter_Eu = (w[:,2]./w[:,1]) .- 1.0
            iter_dell[:] = surf.delu[:]
            iter_El[:] = surf.Eu[:]

            wtu[2:end] = (1/100)*diff(surf.qu.*iter_delu)./diff(surf.x)
            wtu[1] = 2*wtu[2] - wtu[3]
            wtl[:] = wtu[:]

            RHStransp = zeros(surf.ndiv*2-2)

            #Add transpiration velocity to RHS
            for i = 2:surf.ndiv-1
                RHStransp[i-1] = 0.5*(sqrt(1 + (surf.cam_slope[i] + surf.thick_slope[i])^2)*wtu[i]
                                       + sqrt(1 + (surf.cam_slope[i] - surf.thick_slope[i])^2)*wtl[i])

                RHStransp[surf.ndiv+i-3] = 0.5*(sqrt(1 + (surf.cam_slope[i] + surf.thick_slope[i])^2)*wtu[i]
                                                      - sqrt(1 + (surf.cam_slope[i] - surf.thick_slope[i])^2)*wtl[i])
            end

            #Solve viscous problem
            viscsoln = surf.LHS[1:surf.ndiv*2-2, 1:surf.naterm*2] \ RHStransp[1:surf.ndiv*2-2]

            res = sqrt(sum((viscsoln[:] .- [avisc;bvisc]).^2))

            println(res)

            avisc[:] = viscsoln[1:surf.naterm]
            bvisc[:] = viscsoln[surf.naterm+1:end]


            error("here")

        end


        #Calculate adot
        update_atermdot(surf, dt)

        #Set previous values of aterm to be used for derivatives in next time step
        surf.a0prev[1] = surf.a0[1]
        for ia = 1:3
            surf.aprev[ia] = surf.aterm[ia]
        end

        surf.quprev[:] = surf.qu[:]
        surf.qlprev[:] = surf.ql[:]

        #Calculate bound vortex strengths
        update_bv_src(surf)

        #Add effect of transpiration to sources and vortices

        #Wake rollup
        wakeroll(surf, curfield, dt)

        #Force calculation
        cnc, cnnc, cn, cs, cl, cd, int_wax, int_c, int_t = calc_forces(surf, int_wax, int_c, int_t, dt)



        vle = surf.qu[1]

        if vle > 0.
            qspl = Spline1D(surf.x, surf.ql)
            stag = try
                roots(qspl, maxn=1)[1]
OA            catch
                0.
            end
        else
            qspl = Spline1D(surf.x, surf.qu)
            stag = try
                roots(qspl, maxn=1)[1]
            catch
                0.
            end
        end

        mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, vle,
        cl, cd, cnc, cnnc, cn, cs, stag])

        println("here")
    end

    mat = mat'



    f = open("resultsSummary", "w")
    Serialization.serialize(f, ["#time \t", "alpha (deg) \t", "h/c \t", "u/uref \t", "A0 \t", "Cl \t", "Cd \t", "Cm \n"])
    DelimitedFiles.writedlm(f, mat)
    close(f)

    return mat, surf, curfield

end
