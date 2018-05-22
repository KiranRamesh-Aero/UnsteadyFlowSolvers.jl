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

        soln = nlsolve(not_in_place(kelv), -0.01*ones(surf.nspan))

        for i = 1:surf.nspan
            field.f2d[i].tev[end].s = soln.zero[i]

            #Update incduced velocities on airfoil
            update_indbound(surf.s2d[i], field.f2d[i])
            #Calculate downwash
            update_downwash(surf.s2d[i], [field.f2d[i].u[1],field.f2d[i].w[1]])

            #Calculate first two fourier coefficients
            update_a0anda1(surf.s2d[i])
            surf.bc[i] = surf.s2d[i].a0[1] + 0.5*surf.s2d[i].aterm[1]
        end

        calc_a0a13d(surf)
        calc_a2toan3d(surf)

        for i = 1:surf.nspan
            #Update 3D effect on A0 and A1
            surf.s2d[i].a0[1] += surf.a03d[i]
            surf.s2d[i].aterm[1] += surf.aterm3d[1,i]

            #Update rest of Fourier terms
            update_a2toan(surf.s2d[i])
            #Update 3D effect on An
            for ia = 2:surf.naterm
                surf.s2d[i].aterm[ia] += surf.aterm3d[ia,i]
            end

            #Update derivatives of Fourier coefficients
            update_adot(surf.s2d[i],dt)

            #Set previous values of aterm to be used for derivatives in next time step
            surf.s2d[i].a0prev[1] = surf.s2d[i].a0[1]
            for ia = 1:3
                surf.s2d[i].aprev[ia] = surf.s2d[i].aterm[ia]
            end

            #Calculate bound vortex strengths
            update_bv(surf.s2d[i])

            #wakeroll(surf.s2d[i], field.f2d[i], dt)

        end

        cl3d, cd3d, cm3d = calc_forces(surf)



        # write flow details if required
        if writeflag == 1
            if istep in writeArray
                dirname = "$(round(t,nround))"
                writeStamp(dirname, t, surf, field)
            end
        end

        mat = hcat(mat, [t, maximum(map(q->q.a0[1],surf.s2d)), cl3d, cd3d, cm3d])
    end

    mat = mat'

    f = open("resultsSummary", "w")
    write(f, ["#time \t", "CL \t", "CD \t", "CM \n"])
    writedlm(f, mat)
    close(f)

    mat, surf, field

end


function QSLLTlautatRoll(surf :: ThreeDSurfSimple, field :: ThreeDFieldSimple, nsteps :: Int64, dtstar :: Float64, startflag = 0,
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

        soln = nlsolve(not_in_place(kelv), -0.01*ones(surf.nspan))

        for i = 1:surf.nspan
            field.f2d[i].tev[end].s = soln.zero[i]

            #Update incduced velocities on airfoil
            update_indbound(surf.s2d[i], field.f2d[i])
            #Calculate downwash
            update_downwash(surf.s2d[i], [field.f2d[i].u[1],field.f2d[i].w[1]])

            #Calculate first two fourier coefficients
            update_a0anda1(surf.s2d[i])
            surf.bc[i] = surf.s2d[i].a0[1] + 0.5*surf.s2d[i].aterm[1]
        end

        calc_a0a13d(surf)
        calc_a2toan3d(surf)

        for i = 1:surf.nspan
            #Update 3D effect on A0 and A1
            surf.s2d[i].a0[1] += surf.a03d[i]
            surf.s2d[i].aterm[1] += surf.aterm3d[1,i]

            #Update rest of Fourier terms
            update_a2toan(surf.s2d[i])
            #Update 3D effect on An
            for ia = 2:surf.naterm
                surf.s2d[i].aterm[ia] += surf.aterm3d[ia,i]
            end

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

        end

        cl3d, cd3d, cm3d = calc_forces(surf)



        # write flow details if required
        if writeflag == 1
            if istep in writeArray
                dirname = "$(round(t,nround))"
                writeStamp(dirname, t, surf, field)
            end
        end

        mat = hcat(mat, [t, maximum(map(q->q.a0[1],surf.s2d)), cl3d, cd3d, cm3d])
    end

    mat = mat'

    f = open("resultsSummary", "w")
    write(f, ["#time \t", "CL \t", "CD \t", "CM \n"])
    writedlm(f, mat)
    close(f)

    mat, surf, field

end

function QSLLTldvm(surf :: ThreeDSurfSimple, field :: ThreeDFieldSimple, nsteps :: Int64, dtstar :: Float64, startflag = 0,
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

        soln = nlsolve(not_in_place(kelv), -0.01*ones(surf.nspan))

        for i = 1:surf.nspan
            field.f2d[i].tev[end].s = soln.zero[i]

            #Update incduced velocities on airfoil
            update_indbound(surf.s2d[i], field.f2d[i])
            #Calculate downwash
            update_downwash(surf.s2d[i], [field.f2d[i].u[1],field.f2d[i].w[1]])

            #Calculate first two fourier coefficients
            update_a0anda1(surf.s2d[i])
            surf.bc[i] = surf.s2d[i].a0[1] + 0.5*surf.s2d[i].aterm[1]
            #Update rest of Fourier terms
            update_a2toan(surf.s2d[i])
        end

        calc_a0a13d(surf)
        calc_a2toan3d(surf)

        #Update 3D effect on A0-An and calculate derivtives
        for i = 1:surf.nspan
            surf.s2d[i].a0[1] += surf.a03d[i]
            for ia = 1:surf.naterm
                surf.s2d[i].aterm[ia] += surf.aterm3d[ia,i]
            end

            #Update derivatives of Fourier coefficients
            update_adot(surf.s2d[i],dt)
        end

        #Check for LEV formation with the LESP criterion
        nshed = Int(0)
        for i = 1:surf.nspan
            #Condition if LESP_crit is exceeded --> shedding tev + lev
            if abs(surf.s2d[i].a0[1]) > surf.s2d[i].lespcrit[1]   #condition lesp exceeded
                #Remove the previous tev
                pop!(field.f2d[i].tev)

                #Add a TEV with dummy strength
                place_tev(surf.s2d[i], field.f2d[i], dt)

                #Add a LEV with dummy strength
                place_lev(surf.s2d[i], field.f2d[i], dt)

                surf.s2d[i].levflag[1] = 1   #shedding mode "on"
                nshed += 1   #counter increased
            else
                surf.s2d[i].levflag[1] = 0   #shedding mode "off"
            end
        end

        if nshed > 0   #if there is lev shed        ding
            #Kelvin + Kutta Conditions with TEV + LEV
            kelv = KelvinKuttaLLT(surf, field, nshed)

            #Solver for TEV and LEV to satisfy Kelvin + Kutta conditions
            soln = nlsolve(not_in_place(kelv), [-0.01*ones(surf.nspan); 0.01*ones(nshed);])

            cntr = surf.nspan + 1

            for i = 1:surf.nspan
                field.f2d[i].tev[end].s = soln.zero[i]
            end
            for i = 1:surf.nspan
                if surf.s2d[i].levflag[1] == 1
                    field.f2d[i].lev[end].s = soln.zero[cntr]
                    cntr += 1
                end
            end
            for i = 1:surf.nspan
                #Update induced velocities on airfoil
                update_indbound(surf.s2d[i], field.f2d[i])

                #Calculate downwash
                update_downwash(surf.s2d[i], [field.f2d[i].u[1],field.f2d[i].w[1]])

                #Calculate first two fourier coefficients
                update_a0anda1(surf.s2d[i])
                surf.bc[i] = surf.s2d[i].a0[1] + 0.5*surf.s2d[i].aterm[1]
            end

            calc_a0a13d(surf)
            calc_a2toan3d(surf)

            for i = 1:surf.nspan
                #Update 3D effect on A0 and A1
                surf.s2d[i].a0[1] += surf.a03d[i]
                surf.s2d[i].aterm[1] += surf.aterm3d[1,i]

                #Update rest of Fourier terms
                update_a2toan(surf.s2d[i])
                #Update 3D effect on An
                for ia = 2:surf.naterm
                    surf.s2d[i].aterm[ia] += surf.aterm3d[ia,i]
                end
            end
        end

        for i = 1:surf.nspan
            #Set previous values of aterm to be used for derivatives in next time step
            surf.s2d[i].a0prev[1] = surf.s2d[i].a0[1]
            for ia = 1:3
                surf.s2d[i].aprev[ia] = surf.s2d[i].aterm[ia]
            end


            #Calculate bound vortex strengths
            update_bv(surf.s2d[i])

            wakeroll(surf.s2d[i], field.f2d[i], dt)
        end

        cl3d, cd3d, cm3d = calc_forces(surf)

        # write flow details if required
        if writeflag == 1
            if istep in writeArray
                dirname = "$(round(t,nround))"
                writeStamp(dirname, t, surf, field)
            end
        end

        mat = hcat(mat, [t, maximum(map(q->q.a0[1],surf.s2d)), cl3d, cd3d, cm3d])
    end

    mat = mat'

    f = open("resultsSummary", "w")
    write(f, ["#time \t", "CL \t", "CD \t", "CM \n"])
    writedlm(f, mat)
    close(f)

    mat, surf, field
end

function QSLLTlautatLin(surf :: ThreeDSurfSimple, field :: ThreeDFieldSimple, nsteps :: Int64, dtstar :: Float64, startflag = 0,
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
    vcore = 0.02

    T1 = zeros(surf.ndiv)
    T2 = zeros(surf.ndiv)
    T3 = zeros(surf.ndiv)

    for istep = 1:nsteps
        #Udpate current time
        t = t + dt

        LHS = zeros(2*surf.nspan, 2*surf.nspan)
        RHS = zeros(2*surf.nspan)

        for i = 1:surf.nspan
            #Define the flow field
            push!(field.f2d, TwoDFlowField())
            #Update kinematic parameters
            update_kinem(surf.s2d[i], t)

            #Update flow field parameters if any
            update_externalvel(field.f2d[i], t)

            #Update bound vortex positions
            update_boundpos(surf.s2d[i], dt)

            #Update incduced velocities on airfoil
            update_indbound(surf.s2d[i], field.f2d[i])

            #Calculate downwash
            update_downwash(surf.s2d[i], [field.f2d[i].u[1],field.f2d[i].w[1]])

            T1[:] = surf.s2d[i].downwash[:]
            I1 = surf.s2d[i].c*simpleTrapz(T1.*(cos.(surf.s2d[i].theta) - 1), surf.s2d[i].theta)
            J1 = -simpleTrapz(T1,surf.s2d[i].theta)/(surf.s2d[i].uref*pi)

            # T2 depends on recenetly shed TEV
            ntev = length(field.f2d[i].tev)

            if ntev == 0
                xloc_tev = surf.s2d[i].bnd_x[end] + 0.5*surf.s2d[i].kinem.u*dt
                zloc_tev = surf.s2d[i].bnd_z[end]
            else
                xloc_tev = surf.s2d[i].bnd_x[end]+(1./3.)*(field.f2d[i].tev[ntev].x - surf.s2d[i].bnd_x[end])
                zloc_tev = surf.s2d[i].bnd_z[end]+(1./3.)*(field.f2d[i].tev[ntev].z - surf.s2d[i].bnd_z[end])
            end

            for ib = 1:surf.s2d[i].ndiv
                xdist = surf.s2d[i].bnd_x[ib] - xloc_tev
                zdist = surf.s2d[i].bnd_z[ib] - zloc_tev
                distsq = xdist*xdist + zdist*zdist
                T2[ib] = (surf.s2d[i].cam_slope[ib]*zdist + xdist)/(2*pi*sqrt(distsq^2 + vcore^4))
            end

            sig_prev = sum(map(q->q.s, field.f2d[i].tev)) + sum(map(q->q.s, field.f2d[i].lev))

            I2 = simpleTrapz(T2.*(cos.(surf.s2d[i].theta) - 1.), surf.s2d[i].theta)

            LHS[i,i] = 1. + I2
            for n = 1:surf.nspan
                nn = 2*n - 1
                LHS[i,surf.nspan+n] = -real(nn)*sin(nn*surf.psi[i])/sin(surf.psi[i])
            end

            integ = simpleTrapz(surf.s2d[i].cam_slope, surf.s2d[i].theta)
            -simpleTrapz(surf.s2d[i].cam_slope.*cos.(surf.s2d[i].theta), surf.s2d[i].theta)
            for n = 1:surf.nspan
                nn = 2*n - 1
                LHS[i+surf.nspan,n+surf.nspan] = sin(nn*surf.psi[i])*(sin(surf.psi[i]) + (nn*pi/(2*surf.AR))
                *(cos(surf.s2d[i].kinem.alpha) + sin(surf.s2d[i].kinem.alpha)*integ/pi))
            end

            LHS[i+surf.nspan, i] = -sin(surf.psi[i])*I1/(2*surf.AR)

            RHS[i] = -I1 - sig_prev
            RHS[i+surf.nspan] = sin(surf.psi[i])*I1/(2*surf.AR)
        end

        soln = LHS \ RHS

        #Assign and update

        for i = 1:surf.nspan
            tevstr = soln[i]

            ntev = length(field.f2d[i].tev)
            if ntev == 0
                xloc_tev = surf.s2d[i].bnd_x[end] + 0.5*surf.s2d[i].kinem.u*dt
                zloc_tev = surf.s2d[i].bnd_z[end]
            else
                xloc_tev = surf.s2d[i].bnd_x[end]+(1./3.)*(field.f2d[i].tev[ntev].x - surf.s2d[i].bnd_x[end])
                zloc_tev = surf.s2d[i].bnd_z[end]+(1./3.)*(field.f2d[i].tev[ntev].z - surf.s2d[i].bnd_z[end])
            end
            push!(field.f2d[i].tev, TwoDVort(xloc_tev, zloc_tev, tevstr, vcore, 0., 0.))

            surf.bcoeff[i] = soln[i+surf.nspan]
        end

        for i = 1:surf.nspan
            #Update incduced velocities on airfoil
            update_indbound(surf.s2d[i], field.f2d[i])
            #Calculate downwash
            update_downwash(surf.s2d[i], [field.f2d[i].u[1],field.f2d[i].w[1]])

            #Calculate first two fourier coefficients
            update_a0anda1(surf.s2d[i])
            surf.bc[i] = surf.s2d[i].a0[1] + 0.5*surf.s2d[i].aterm[1]
        end

        calc_a0a13d(surf)
        calc_a2toan3d(surf)

        for i = 1:surf.nspan
            #Update 3D effect on A0 and A1
            surf.s2d[i].a0[1] += surf.a03d[i]
            surf.s2d[i].aterm[1] += surf.aterm3d[1,i]

            #Update rest of Fourier terms
            update_a2toan(surf.s2d[i])
            #Update 3D effect on An
            for ia = 2:surf.naterm
                surf.s2d[i].aterm[ia] += surf.aterm3d[ia,i]
            end

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
        end

        cl3d, cd3d, cm3d = calc_forces(surf)

        # write flow details if required
        if writeflag == 1
            if istep in writeArray
                dirname = "$(round(t,nround))"
                writeStamp(dirname, t, surf, field)
            end
        end

        mat = hcat(mat, [t, maximum(map(q->q.a0[1],surf.s2d)), cl3d, cd3d, cm3d])
    end

    mat = mat'

    f = open("resultsSummary", "w")
    write(f, ["#time \t", "CL \t", "CD \t", "CM \n"])
    writedlm(f, mat)
    close(f)

    mat, surf, field
end
