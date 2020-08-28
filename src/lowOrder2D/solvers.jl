"""
    lautat(surf, curfield, nsteps=500, dtstar=0.015, startflag=0,
writeflag=0, writeInterval=1000., delvort=delNone(), maxwrite=100,
nround=6)

Simulates potential flow for an airfoil undergoing unsteady motion
    using Large-Angle Unsteady Thin-Airfoil Theory.

    Inputs required:
    - `TwoDSurf` : based on airfoil geometry and kinematics
    - `TwoDFlowField` : initial flowfield

    `dtstar` should be in accordance with the kinematics; see `?find_tstep`.

    `nsteps` should be calculated according to total run time required and
    dtstar.

    `startflag=1` can be used to continue a simulation using the
    *resultsSummary* file in current directory.

    `writeflag=1` along with a suitable writeInterval is used to write flow
    data at intervals. These can be used to create plots and movies of the
    flowfield.

    `delvort` is used to control the vortex count. The simulation cost
    increases as O(n``^2``) with number of vortices in the flow field and so
    long run times can get very computationally expensive. Look at
    `?delVortDef` for algorthims to control vortex count.

    `maxwrite` is the maximum number of timestamps/directories that may be
    written.

    `nround` is the number of digits to which the output timestamps are
    rounded off to.

    The Large-Angle Unsteady Thin-Airfoil theory is described in: Ramesh,
    K. et al., "An unsteady airfoil theory applied to pitching motions
    validated against experiment and computation", Theor. Comput. Fluid
    Dyn. (2013) 27: 843. [Weblink](https://doi.org/10.1007/s00162-012-0292-8)

"""
function lautat(surf::TwoDSurf, curfield::TwoDFlowField, nsteps::Int64 = 500, dtstar::Float64 = 0.015, startflag = 0, writeflag = 0, writeInterval = 1000., delvort = delNone(); maxwrite = 50, nround=6, wakerollup=1)

    # If a restart directory is provided, read in the simulation data
    if startflag == 0
        mat = zeros(0, 8)
        t = 0.
    elseif startflag == 1
        dirvec = readdir()
        dirresults = map(x->(v = tryparse(Float64,x); typeof(v) == Nothing ? 0.0 : v),dirvec)
        latestTime = maximum(dirresults)
        mat = DelimitedFiles.readdlm("resultsSummary",header = true)
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

    #Intialise flowfield
    for istep = 1:nsteps
        #Udpate current time
        t = t + dt

        #Update kinematic parameters
        update_kinem(surf, t)

        #Update bound vortex positions
        update_boundpos(surf, dt)

        #Add a TEV with dummy strength
        place_tev(surf,curfield,dt)

        #Solve for TEV strength to satisfy Kelvin condition
        kelv = KelvinCondition(surf,curfield)
        soln = nlsolve(not_in_place(kelv), [-0.01])
        curfield.tev[length(curfield.tev)].s = soln.zero[1]

        #Update adot
        update_a2a3adot(surf,dt)

        #Set previous values of aterm to be used for derivatives in next time step
        surf.a0prev[1] = surf.a0[1]
        for ia = 1:3
            surf.aprev[ia] = surf.aterm[ia]
        end

        #Update rest of Fourier terms
        update_a2toan(surf)

        #Calculate bound vortex strengths
        update_bv(surf)

        # Delete or merge vortices if required
        controlVortCount(delvort, surf.bnd_x[Int(round(surf.ndiv/2))], surf.bnd_z[Int(round(surf.ndiv/2))], curfield)

        #Wake rollup if flag on
        if wakerollup == 1
            wakeroll(surf, curfield, dt)
        end

        # Calculate force and moment coefficients
        cl, cd, cm = calc_forces(surf, [curfield.u[1], curfield.w[1]])

        # write flow details if required
        if writeflag == 1
            if istep in writeArray
                dirname = "$(round(t, sigdigits=nround))"
                writeStamp(dirname, t, surf, curfield)
            end
        end

        # for writing in resultsSummary
        mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, surf.a0[1], cl, cd, cm])

    end

    mat = mat'

    f = open("resultsSummary", "w")
    println(f,"#time \t","alpha(deg) \t","h/c \t", "u/uref \t", "A0 \t", "Cl \t", "Cd \t", "Cm ")
    DelimitedFiles.writedlm(f, mat)
    close(f)

    mat, surf, curfield

end

function ldvm(surf::TwoDSurf, curfield::TwoDFlowField, nsteps::Int64 = 500, dtstar::Float64 = 0.015, startflag = 0, writeflag = 0, writeInterval = 1000., delvort = delNone(); maxwrite = 50, nround=6)

    # If a restart directory is provided, read in the simulation data
    if startflag == 0
        mat = zeros(0, 8)
        t = 0.
    elseif startflag == 1
        dirvec = readdir()
        dirresults = map(x->(v = tryparse(Float64,x); typeof(v) == Nothing ? 0.0 : v),dirvec)
        latestTime = maximum(dirresults)
        mat = DelimitedFiles.readdlm("resultsSummary",header = true)
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

    # time stepping
    for istep = 1:nsteps
        #Udpate current time
        t = t + dt

        #Update external flowfield
        update_externalvel(curfield, t)

        #Update kinematic parameters
        update_kinem(surf, t)

        #Update bound vortex positions
        update_boundpos(surf, dt)

        #Add a TEV with dummy strength
        place_tev(surf,curfield,dt)

        #Update incduced velocities on airfoil
        update_indbound(surf, curfield)

        #Solve for TEV strength to satisfy Kelvin condition
        kelv = KelvinCondition(surf,curfield)
        soln = nlsolve(not_in_place(kelv), [-0.01])
        curfield.tev[length(curfield.tev)].s = soln.zero[1]

        #Update adot
        update_a2a3adot(surf,dt)

        lesp = surf.a0[1]

        #Check for LESP condition
        if abs(lesp)>surf.lespcrit[1]
            #2D iteration if LESP_crit is exceeded

            #Add a LEV with dummy strength
            place_lev(surf,curfield,dt)

            #Solve for TEV and LEV strengths to satisfy Kelvin condition and Kutta condition at leading edge
            kelvkutta = KelvinKutta(surf,curfield)
            soln = nlsolve(not_in_place(kelvkutta), [-0.01; 0.01])
            (curfield.tev[length(curfield.tev)].s, curfield.lev[length(curfield.lev)].s) = soln.zero[1], soln.zero[2]

            #set flag for levshedding=on
            surf.levflag[1] = 1
        else
            surf.levflag[1] = 0
        end

        #Update rest of Fourier terms
        update_a2toan(surf)

        #Set previous values of aterm to be used for derivatives in next time step
        surf.a0prev[1] = surf.a0[1]
        for ia = 1:3
            surf.aprev[ia] = surf.aterm[ia]
        end

        #Calculate bound vortex strengths
        update_bv(surf)

        # Delete or merge vortices if required
        controlVortCount(delvort, surf.bnd_x[Int(round(surf.ndiv/2))], surf.bnd_z[Int(round(surf.ndiv/2))], curfield)

        # free wake rollup
        wakeroll(surf, curfield, dt)

        # Calculate force and moment coefficients
        cl, cd, cm, cn = calc_forces(surf, [curfield.u[1], curfield.w[1]])

        # write flow details if required
        if writeflag == 1
            if istep in writeArray
                dirname = "$(round(t, sigdigits=nround))"
                writeStamp(dirname, t, surf, curfield)
            end
        end

        # for writing in resultsSummary
        mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, surf.a0[1], cl, cd, cm])

    end

    mat = mat'

    f = open("resultsSummary", "w")
    println(f,"#time \t","alpha(deg) \t","h/c \t", "u/uref \t", "A0 \t", "Cl \t", "Cd \t", "Cm ")
    DelimitedFiles.writedlm(f, mat)
    close(f)

    mat, surf, curfield

end

function ldvm(surf::Vector{TwoDSurf}, curfield::TwoDFlowField, nsteps::Int64 = 500, dtstar::Float64 = 0.015, startflag = 0, writeflag = 0, writeInterval = 1000., delvort = delNone(); maxwrite = 50, nround=6)

    nsurf = length(surf)

    # If a restart directory is provided, read in the simulation data
    if startflag == 0
        mat = zeros(0, 7*nsurf+1)
        t = 0.
    elseif startflag == 1
        dirvec = readdir()
        dirresults = map(x->(v = tryparse(Float64,x); typeof(v) == Nothing ? 0.0 : v),dirvec)
        latestTime = maximum(dirresults)
        mat = DelimitedFiles.readdlm("resultsSummary",header = true)
        t = mat[end,1]
    else
        throw("invalid start flag, should be 0 or 1")
    end
    mat = mat'

    dt = dtstar*minimum(map(q->q.c/q.uref, surf))

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

    lesp = zeros(nsurf)
    cl = zeros(nsurf)
    cd = zeros(nsurf)
    cm = zeros(nsurf)

    # time stepping
    for istep = 1:nsteps
        #Udpate current time
        t = t + dt

        #Update external flowfield
        update_externalvel(curfield, t)

        for is = 1:nsurf
            #Update kinematic parameters
            update_kinem(surf[is], t)

            #Update bound vortex positions
            update_boundpos(surf[is], dt)
        end

        #Add TEVs with dummy strength
        place_tev(surf,curfield,dt)

        for is = 1:nsurf
            #Update induced velocities on airfoil
            update_indbound(surf[is], curfield)
            for js = 1:nsurf
                if is != js
                    #add influence of other surfaces to downwash
                    add_indbound_b(surf[is], surf[js])
                end
            end
        end


        #Solve for TEV strength to satisfy Kelvin condition
        kelv = KelvinConditionMult(surf,curfield)
        soln = nlsolve(not_in_place(kelv), ones(nsurf)*-0.01)

        for i = 1:nsurf
            curfield.tev[end-nsurf+i].s = soln.zero[i]

            #Update adot
            update_a2a3adot(surf[i],dt)

            lesp[i] = surf[i].a0[1]
        end

        shed_ind = Int[]

        #Check for LESP condition
        for i = 1:nsurf
            if (abs(lesp[i])>surf[i].lespcrit[1])

                push!(shed_ind, i)
            end
        end
        nshed = length(shed_ind)

        if nshed > 0
            #2D iteration if LESP_crit is exceeded
            #Add LEVs with dummy strength
            place_lev(surf,curfield,dt,shed_ind)

            #Solve for TEV and LEV strengths to satisfy Kelvin condition and Kutta condition at leading edge
            kelvkutta = KelvinKuttaMult(surf,curfield,shed_ind)
            soln = nlsolve(not_in_place(kelvkutta), -0.01*ones(nshed+nsurf))

            for i = 1:nsurf
                if i in shed_ind
                    curfield.tev[end-nsurf+i].s = soln.zero[i]
                    curfield.lev[end-nsurf+shed_ind[i]].s = soln.zero[nsurf+i]
                    surf[i].levflag[1] = 1
                else
                    curfield.tev[end-nsurf+i].s = soln.zero[i]
                    surf[i].levflag[1] = 0
                end
            end
        else
            for i = 1:nsurf
                surf[i].levflag[1] = 0
            end
        end

        for i = 1:nsurf
            #Set previous values of aterm to be used for derivatives in next time step
            surf[i].a0prev[1] = surf[i].a0[1]
            for ia = 1:3
                surf[i].aprev[ia] = surf[i].aterm[ia]
            end
        end

        # Delete or merge vortices if required
        mean_bndx = sum(map(q->q.bnd_x[Int(round(q.ndiv/2))], surf))/nsurf
        mean_bndz = sum(map(q->q.bnd_z[Int(round(q.ndiv/2))], surf))/nsurf
        controlVortCount(delvort, mean_bndx, mean_bndz, curfield)

        # free wake rollup
        wakeroll(surf, curfield, dt)

        # Calculate force and moment coefficients
        for i = 1:nsurf
            cl[i], cd[i], cm[i] = calc_forces(surf[i], [curfield.u[1], curfield.w[1]])
        end

        # write flow details if required
        if writeflag == 1
            if istep in writeArray
                dirname = "$(round(t, sigdigits=nround))"
                writeStamp(dirname, t, surf, curfield)
            end
        end

        # for writing in resultsSummary
        matvect = [t;]
        for is = 1:nsurf
            matvect = [matvect;[surf[is].kinem.alpha, surf[is].kinem.h,
                                surf[is].kinem.u, surf[is].a0[1], cl[is], cd[is], cm[is]]]
        end
        mat = hcat(mat,matvect)

    end

mat = mat'

f = open("resultsSummary", "w")
println(f,"#time \t","alpha(deg) \t","h/c \t", "u/uref \t", "A0 \t", "Cl \t", "Cd \t", "Cm ")
DelimitedFiles.writedlm(f, mat)
close(f)

mat, surf, curfield

end


function ldvmLin(surf::TwoDSurf, curfield::TwoDFlowField, nsteps::Int64 = 500, dtstar::Float64 = 0.015, startflag = 0, writeflag = 0, writeInterval = 1000., delvort = delNone(); maxwrite = 50, nround=6)

    # If a restart directory is provided, read in the simulation data
    if startflag == 0
        mat = zeros(0, 8)
        t = 0.
    elseif startflag == 1
        dirvec = readdir()
        dirresults = map(x->(v = tryparse(Float64,x); typeof(v) == Nothing ? 0.0 : v),dirvec)
        latestTime = maximum(dirresults)
        mat = DelimitedFiles.readdlm("resultsSummary",header=true)
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

    T1 = zeros(surf.ndiv)
    T2 = zeros(surf.ndiv)
    T3 = zeros(surf.ndiv)

    # time stepping
    for istep = 1:nsteps
        #Udpate current time
        t = t + dt

        #Update external flowfield
        update_externalvel(curfield, t)

        #Update kinematic parameters
        update_kinem(surf, t)

        #Update bound vortex positions
        update_boundpos(surf, dt)

        #Update induced velocities on airfoil
        update_indbound(surf, curfield)

        #Calculate downwash
        update_downwash(surf, [curfield.u[1],curfield.w[1]])

        #The integrals I1 and J1 are based on this downwash which doesnt
        #include vortives shed at current step (equal to term T1)
        #I1 has units if circulation, J1 is dimensionless
        T1[:] = surf.downwash[:]
        I1 = surf.c*simpleTrapz(T1.*(cos.(surf.theta) .- 1. ), surf.theta)
        J1 = -simpleTrapz(T1,surf.theta)/(surf.uref*pi)

        # T2 depends on recenetly shed TEV
        ntev = length(curfield.tev)

        if ntev == 0
            xloc_tev = surf.bnd_x[surf.ndiv] + 0.5*surf.kinem.u*dt
            zloc_tev = surf.bnd_z[surf.ndiv]
        else
            xloc_tev = surf.bnd_x[surf.ndiv]+(1. /3.)*(curfield.tev[ntev].x - surf.bnd_x[surf.ndiv])
            zloc_tev = surf.bnd_z[surf.ndiv]+(1. /3.)*(curfield.tev[ntev].z - surf.bnd_z[surf.ndiv])
        end

        for ib = 1:surf.ndiv
            xdist = surf.bnd_x[ib] - xloc_tev
            zdist = surf.bnd_z[ib] - zloc_tev
            distsq = xdist*xdist + zdist*zdist
            T2[ib] = (surf.cam_slope[ib]*zdist + xdist)/(2*pi*sqrt(distsq^2 + vcore^4))
        end

        #sig_prev = sum(map(q->q.s, curfield.tev)) + sum(map(q->q.s, curfield.lev))
        sig_prev = -surf.uref*surf.c*pi*(surf.a0prev[1] + surf.aprev[1]/2. )

        I2 = simpleTrapz(T2.*(cos.(surf.theta) .- 1. ), surf.theta)
        J2 = -simpleTrapz(T2, surf.theta)/(pi*surf.uref)

        tevstr = -(I1 + sig_prev)/(1 + I2)

        #Calc first 3 fourier coefficients and derivatives
        surf.a0[1] = J1 + J2*tevstr
        for ia = 1:surf.naterm
            surf.aterm[ia] = 2. *(simpleTrapz(T1.*cos.(ia*surf.theta), surf.theta) + tevstr*simpleTrapz(T2.*cos.(ia*surf.theta), surf.theta))/(pi*surf.uref)
        end

        #Calculate adot
        surf.a0dot[1] = (surf.a0[1] - surf.a0prev[1])/dt
        for ia = 1:surf.naterm
            surf.adot[ia] = (surf.aterm[ia]-surf.aprev[ia])/dt
        end

        #Check if LEV shedding is true
        if abs(surf.a0[1]) > surf.lespcrit[1]
            if (surf.a0[1] >= 0.)
                lesp_cond = surf.lespcrit[1]
            else
                lesp_cond = -surf.lespcrit[1]
            end

            # T3 depends on recenetly shed LEV
            nlev = length(curfield.lev)
            if surf.levflag[1] == 0
                le_vel_x = surf.kinem.u - surf.kinem.alphadot*sin(surf.kinem.alpha)*surf.pvt*surf.c + surf.uind[1]
                le_vel_z = -surf.kinem.alphadot*cos(surf.kinem.alpha)*surf.pvt*surf.c- surf.kinem.hdot + surf.wind[1]
                xloc_lev = surf.bnd_x[1] + 0.5*le_vel_x*dt
                zloc_lev = surf.bnd_z[1] + 0.5*le_vel_z*dt
            else
                xloc_lev = surf.bnd_x[1] + (1. /3.)*(curfield.lev[nlev].x - surf.bnd_x[1])
                zloc_lev = surf.bnd_z[1]+(1. /3.)*(curfield.lev[nlev].z - surf.bnd_z[1])
            end

            for ib = 1:surf.ndiv
                xdist = surf.bnd_x[ib] - xloc_lev
                zdist = surf.bnd_z[ib] - zloc_lev
                distsq = xdist*xdist + zdist*zdist
                T3[ib] = (surf.cam_slope[ib]*zdist + xdist)/(2*pi*sqrt(distsq^2 + vcore^4))
            end
            I3 = simpleTrapz(T3.*(cos.(surf.theta) .- 1. ), surf.theta)
            J3 = -simpleTrapz(T3, surf.theta)/(pi*surf.uref)

            det = J3*(I2 + 1. ) - J2*(I3 + 1. )

            tevstr = (-J3*(I1 + sig_prev) + (I3 + 1)*(J1 - lesp_cond))/det
            levstr = (J2*(I1 + sig_prev) - (I2 + 1)*(J1 - lesp_cond))/det

            #Recalculate required fourier terms
            surf.a0[1] = J1 + J2*tevstr + J3*levstr
            for ia = 1:3
                surf.aterm[ia] = 2. *(simpleTrapz(T1.*cos.(ia*surf.theta), surf.theta) +
                                      tevstr*simpleTrapz(T2.*cos.(ia*surf.theta), surf.theta) +
                                      levstr*simpleTrapz(T3.*cos.(ia*surf.theta), surf.theta))/(pi*surf.uref)
            end

            push!(curfield.tev, TwoDVort(xloc_tev, zloc_tev, tevstr, vcore, 0., 0.))
            push!(curfield.lev, TwoDVort(xloc_lev, zloc_lev, levstr, vcore, 0., 0.))

            for ia = 4:surf.naterm
                surf.aterm[ia] = 2. *(simpleTrapz(T1.*cos.(ia*surf.theta), surf.theta) +
                                      tevstr*simpleTrapz(T2.*cos.(ia*surf.theta), surf.theta) +
                                      levstr*simpleTrapz(T3.*cos.(ia*surf.theta), surf.theta))/(pi*surf.uref)
            end

            surf.levflag[1] = 1
        else
            push!(curfield.tev, TwoDVort(xloc_tev, zloc_tev, tevstr, vcore, 0., 0.))

            for ia = 4:surf.naterm
                surf.aterm[ia] = 2. *(simpleTrapz(T1.*cos.(ia*surf.theta), surf.theta) +
                                      tevstr*simpleTrapz(T2.*cos.(ia*surf.theta), surf.theta))/(pi*surf.uref)
            end

            surf.levflag[1] = 0
        end

#Set previous values of aterm to be used for derivatives in next time step
surf.a0prev[1] = surf.a0[1]
for ia = 1:3
    surf.aprev[ia] = surf.aterm[ia]
end

#Calculate bound vortex strengths
update_bv(surf)

# Delete or merge vortices if required
controlVortCount(delvort, surf.bnd_x[Int(round(surf.ndiv/2))], surf.bnd_z[Int(round(surf.ndiv/2))], curfield)

# free wake rollup
wakeroll(surf, curfield, dt)

# Calculate force and moment coefficients
cl, cd, cm = calc_forces(surf, [curfield.u[1], curfield.w[1]])

# write flow details if required
if writeflag == 1
    if istep in writeArray
        dirname = "$(round(t,sigdigits=nround))"
        writeStamp(dirname, t, surf, curfield)
    end
end

# for writing in resultsSummary
mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, surf.a0[1], cl, cd, cm])

end

mat = mat'

f = open("resultsSummary", "w")
println(f,"#time \t","alpha(deg) \t","h/c \t", "u/uref \t", "A0 \t", "Cl \t", "Cd \t", "Cm ")
DelimitedFiles.writedlm(f, mat)
close(f)

mat, surf, curfield

end

function ldvm2DOF(surf::TwoDSurf, curfield::TwoDFlowField, strpar::TwoDOFPar, kinem::KinemPar2DOF, nsteps::Int64 = 500, dtstar::Float64 = 0.015, startflag = 0, writeflag = 0, writeInterval = 1000., delvort = delNone(); maxwrite = 50, nround=6)

    # If a restart directory is provided, read in the simulation data
    if startflag == 0
        mat = zeros(0, 8)
        t = 0.
    elseif startflag == 1
        dirvec = readdir()
        dirresults = map(x->(v = tryparse(Float64,x); typeof(v) == Nothing ? 0.0 : v),dirvec)
        latestTime = maximum(dirresults)
        mat = DelimitedFiles.readdlm("resultsSummary",header=true)
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



    cl = 0.
    cm = 0.

    #Intialise flowfield

    for istep = 1:nsteps
        #Udpate current time
        t = t + dt

        #Update external flowfield
        update_externalvel(curfield, t)

        #Update kinematic parameters (based on 2DOF response)
        update_kinem2DOF(surf, strpar, kinem, dt, cl, cm)

        #Update bound vortex positions
        update_boundpos(surf, dt)

        #Add a TEV with dummy strength
        place_tev(surf,curfield,dt)

        #Update incduced velocities on airfoil
        update_indbound(surf, curfield)

        #Solve for TEV strength to satisfy Kelvin condition
        kelv = KelvinCondition(surf,curfield)
        soln = nlsolve(not_in_place(kelv), [-0.01])
        curfield.tev[length(curfield.tev)].s = soln.zero[1]

        #Update adot
        #update_a2a3adot(surf,dt)
        update_atermdot(surf, dt)

        lesp = surf.a0[1]

        #Check for LESP condition
        if (abs(lesp)>surf.lespcrit[1])
            #2D iteration if LESP_crit is exceeded
            #Remove the previous tev

            #Add a LEV with dummy strength
            place_lev(surf,curfield,dt)

            #Solve for TEV and LEV strengths to satisfy Kelvin condition and Kutta condition at leading edge
            kelvkutta = KelvinKutta(surf,curfield)
            soln = nlsolve(not_in_place(kelvkutta), [-0.01; 0.01])
            (curfield.tev[length(curfield.tev)].s, curfield.lev[length(curfield.lev)].s) = soln.zero[1], soln.zero[2]

            #set flag for levshedding=on
            surf.levflag[1] = 1
        else
            surf.levflag[1] = 0
        end

        #Update rest of Fourier terms
        update_a2toan(surf)

        #Set previous values of aterm to be used for derivatives in next time step
        surf.a0prev[1] = surf.a0[1]
        for ia = 1:surf.naterm
            surf.aprev[ia] = surf.aterm[ia]
        end

        #Calculate bound vortex strengths
        update_bv(surf)

        # Delete or merge vortices if required
        controlVortCount(delvort, surf.bnd_x[Int(round(surf.ndiv/2))], surf.bnd_z[Int(round(surf.ndiv/2))], curfield)

        # free wake rollup
        wakeroll(surf, curfield, dt)

        # Calculate force and moment coefficients
        cl, cd, cm = calc_forces(surf, [curfield.u[1], curfield.w[1]])

        # write flow details if required
        if writeflag == 1
            if istep in writeArray
                dirname = "$(round(t,sigdigits=nround))"
                writeStamp(dirname, t, surf, curfield)
            end
        end

        # for writing in resultsSummary
        mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, surf.a0[1], cl, cd, cm])

    end

    mat = mat'

    f = open("resultsSummary", "w")
    println(f,"#time \t","alpha(deg) \t","h/c \t", "u/uref \t", "A0 \t", "Cl \t", "Cd \t", "Cm ")
    DelimitedFiles.writedlm(f, mat)
    close(f)

    mat, surf, curfield
end
