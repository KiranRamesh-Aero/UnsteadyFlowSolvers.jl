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
