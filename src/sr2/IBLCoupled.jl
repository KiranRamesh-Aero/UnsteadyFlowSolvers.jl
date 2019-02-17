function IBLCoupled(surf::TwoDSurf, curfield::TwoDFlowField, ncell::Int64, nsteps::Int64 = 500, dtstar::Float64 = 0.015, startflag = 0, writeflag = 0, writeInterval = 1000., delvort = delNone(); maxwrite = 50, nround=6, wakerollup=1)
    # If a restart directory is provided, read in the simulation data
    if startflag == 0
        mat = zeros(0, 8)
        t = 0.
        tv = 0.
        del = zeros(ncell-1)
        E   = zeros(ncell-1)

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

    # initial momentum and energy shape factors

    del, E, x = initViscous(ncell-1)

    # time loop
    for istep = 1:nsteps
        t = t + dt
        mat, surf, curfield = lautat(surf, curfield, t, dt, istep, startflag, writeflag, writeInterval, delvort)
        qu, ql = calc_edgeVel(surf, [curfield.u[1], curfield.w[1]])
        w0u, Uu,  Utu, Uxu = inviscidInterface(qu, dt)
        #w0l, Ul, Utl, Uxl = inviscidInterface(ql, dt)

        while tv < t

            w, dtv, j1 ,j2 = FVMIBL(w0u, Uu, Utu, Uxu);
            ww = w
            del = ww[:,1]
            E = (ww[:,2]./ww[:,1]) .- 1.0

        # the plots of del and E

        #display(plot(x/pi,[del, E], xticks = 0:0.1:1, layout=(2,1), legend = false))

        # convergence of the Eigen-value
            p1 = plot(x[1:end-1]/pi,j1)
            p2 = plot(x[1:end-1]/pi,j2)
            p3 = plot(x/pi,del)
            p4 = plot(x/pi,E)
            display(plot(p1, p2, p3, p4, xticks = 0:0.1:1, layout=(2,2), legend = false))

        #display(plot(sep, xticks = 0:10:200, legend = false))
            sleep(0.05)
            w0u = w;

            @printf("viscous Time :%1.10f , viscous Time step size %1.10f \n", tv, dtv);
        #map(x -> @sprintf("Seperation occure at :%1.10f  \n",x), xtrunc[seperation]./pi);

            tv = tv + dtv

            sols = w
        end
    end

    mat = mat'

    f = open("resultsSummary", "w")
    Serialization.serialize(f, ["#time \t", "alpha (deg) \t", "h/c \t", "u/uref \t", "A0 \t", "Cl \t", "Cd \t", "Cm \n"])
    DelimitedFiles.writedlm(f, mat)
    close(f)

mat, surf, curfield

end

function lautat(surf::TwoDSurf, curfield::TwoDFlowField, t::Float64, dt::Float64, istep::Int64, startflag = 0, writeflag = 0, writeInterval = 1000., delvort = delNone(); maxwrite = 50, nround=6, wakerollup=1)


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

        return  mat, surf, curfield
    end
end

function inviscidInterface(del::Array{Float64,1}, E::Array{Float64,1}, q::Array{Float64,1}, dt::Float64)

    n =length(q)
    #del, E, F ,B = init(n-1)


    U0 = q[1:n]
    #x =x[1:n]
    U0 = U0[1:n]

    Ut = temporalDerivates(q, qu0, dt)
    Ux = spatialDerivates(q)

    w1 = del
    w2 = del.*(E.+1.0)
    w0 = hcat(w1,w2)


    return U0, Ut, Ux
end

function spatialDerivates(qu::Array{Float64,1})

    n = length(qu)
    dx = 1./n

    dudx = qu[2:end]-qu[1:end-1]./dx

    return dudx
end

function temporalDerivates(qu::Array{Float64,1}, qu0::Array{Float64}, dt:Float64)

    return qu[2:end]-qu0[2:end]./dt

end

function initViscous(ncell::Int64)


    del, E = init(ncell-1)
    x = collect(0:ncell-1)./(ncell-1)


    return del, E, x



end
