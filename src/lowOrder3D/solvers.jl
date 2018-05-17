function QSLLTlautat(surf :: ThreeDSurfSimple, field :: ThreeDFieldSimple, nsteps :: Int64, dtstar :: Float64, startflag = 0,
    writeflag = 0, writeInterval = 1000., delvort = delNone(); maxwrite = 100, nround=6)

    # If a restart directory is provided, read in the simulation data
    if startflag == 0
        mat = Array{Float64}(0, 5)
        t = 0.
    elseif startflag == 1
        dirvec = readdir()
        dirresults = map(x->(v = tryparse(Float64,x); isnull(v) ? 0.0 : get(v)),dirvec)
        latestTime = maximum(dirresults)
        mat = readdlm("resultsSummary")
        t = mat[end,1]
    else
        throw("invalid start flag, should be 0 or 1")
    end
    mat = mat'

    # if writeflag is on, determine the timesteps to write at
    if writeflag == 1
        writeArray = Int64[]
        tTot = nsteps*dtstar
        for i = 1:maxwrite
            tcur = writeInterval*real(i)
            if tcur > tTot
                break
            else
                push!(writeArray, Int(round(tcur/dtstar)))
            end
        end
    end

    dt = dtstar*surf.cref/surf.uref

    cl = zeros(surf.nspan)
    cd = zeros(surf.nspan)
    cm = zeros(surf.nspan)

    for istep = 1:nsteps
        #Udpate current time
        t = t + dt

        for i = 1:surf.nspan
            #Define the flow field
            push!(field.f2d, TwoDFlowField())
            #Update kinematic parameters
            update_kinem(surf.s2d[i], t)

            #Update flow field parameters if any
            update_externalvel(field.f2d[i], t)

            #Update bound vortex positions
            update_boundpos(surf.s2d[i], dt)

            #Add a TEV with dummy strength
            place_tev(surf.s2d[i], field.f2d[i], dt)
        end

        kelv = KelvinConditionLLT(surf, field)

        #Solve for TEV strength to satisfy Kelvin condition

        soln = nlsolve(not_in_place(kelv), [-0.01*ones(surf.nspan); zeros(surf.nspan)])

        for i = 1:surf.nspan
            field.f2d[i].tev[length(field.f2d[i].tev)].s = soln.zero[i]

            #Update incduced velocities on airfoil
            update_indbound(surf.s2d[i], field.f2d[i])
            #Calculate downwash
            update_downwash(surf.s2d[i], [field.f2d[i].u[1],field.f2d[i].w[1]])

            #Calculate first two fourier coefficients
            update_a0anda1(surf.s2d[i])
            surf.bc[i] = surf.s2d[i].a0[1] + 0.5*surf.s2d[i].aterm[1]
        end

        for i = 1:surf.nspan
            surf.a03d[i] = 0
            for n = 1:surf.nspan
                nn = 2*n - 1
                surf.a03d[i] = surf.a03d[i] - real(nn)*soln.zero[n+surf.nspan]*sin(nn*surf.psi[i])/sin(surf.psi[i])
            end

        end


        for i = 1:surf.nspan
            #Update 3D effect on A0
            surf.s2d[i].a0[1] = surf.s2d[i].a0[1] + surf.a03d[i]

            #Update rest of Fourier terms
            update_a2toan(surf.s2d[i])

            #Update derivatives of Fourier coefficients
            update_adot(surf.s2d[i],dt)

            #Set previous values of aterm to be used for derivatives in next time step
            surf.s2d[i].a0prev[1] = surf.s2d[i].a0[1]
            for ia = 1:3
                surf.s2d[i].aprev[ia] = surf.s2d[i].aterm[ia]
            end

            #Calculate bound vortex strengths
            update_bv(surf.s2d[i])

            wakeroll(surf.s2d[i], field.f2d[i], dt)

            cl[i], cd[i], cm[i] = calc_forces(surf.s2d[i])
            end

        end

        cl3d = 0
        cd3d = 0
        cm3d = 0

        for i = 1:surf.nspan-1
            cl3d = cl3d + 0.5*(cl[i] + cl[i+1])*sin(0.5*(surf.psi[i] + surf.psi[i+1]))*(surf.psi[i+1] - surf.psi[i])/2
            cd3d = cd3d + 0.5*(cd[i] + cd[i+1])*sin(0.5*(surf.psi[i] + surf.psi[i+1]))*(surf.psi[i+1] - surf.psi[i])/2
            cm3d = cm3d + 0.5*(cm[i] + cm[i+1])*sin(0.5*(surf.psi[i] + surf.psi[i+1]))*(surf.psi[i+1] - surf.psi[i])/2
        end

        mat = hcat(mat, [t, maximum(map(q->q.a0[1],surf.s2d)), cl3d, cd3d, cm3d])
    end
    mat = mat'
    mat, surf, field

end
