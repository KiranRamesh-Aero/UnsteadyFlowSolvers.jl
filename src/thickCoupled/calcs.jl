function initDelE(n)

    E = 0.4142 * ones(n)
    B = 131.9*E.^3 - 167.32*E.^2 + 76.642.*E .- 11.068
    del = sqrt.(B.*0.005);

    #Same for upper and lower surface at initialisation
    return del, del, E, E

end


function correlate(w)

    # the model correlations for quantities F, B, S, dedf based on del and E

    del = w[:,1]
    E = w[:,2]./del .- 1

    F = 4.8274*E.^4 - 5.9816*E.^3 + 4.0274*E.^2 + 0.23247.*E .+ 0.15174

    B = 0.5*(-225.86.*E.^3 - 3016.6.*E.^2 - 208.68*E .- 17.915
             + 131.9*E.^3 - 167.32*E.^2 + 76.642.*E .- 11.068)
    B[E.<-0.0616] = -225.86.*E[E.<-0.0616].^3 - 3016.6.*E[E.<-0.0616].^2 - 208.68*E[E.<-0.0616] .- 17.915
    B[E.>-0.0395] =  131.9*E[E.>-0.0395].^3 - 167.32*E[E.>-0.0395].^2 + 76.642.*E[E.>-0.0395] .- 11.068

    S = 0.5*(451.55*E.^3 + 2010*E.^2 + 138.96*E .+ 11.296
             - 96.739*E.^3 + 117.74*E.^2 - 46.432*E .+ 6.8074)
    S[E.<-0.0582] = 451.55*E[E.<-0.0582].^3 + 2010*E[E.< -0.0582].^2 + 138.96*E[E.< -0.0582] .+ 11.296
    S[E.>-0.042]  = -96.739*E[E.>-0.042].^3 + 117.74*E[E.> -0.042].^2 - 46.432*E[E.> -0.042] .+ 6.8074

    dfde = 4*4.8274*E.^3 - 3*5.9816*E.^2 + 2*4.0274*E .+ 0.23247

    return del, E, F ,B, S, dfde
end

function smoothEdgeVelocity(qu::Array{Float64,1}, theta::Array{Float64,1}, ncell::Int64, xCoarse::Int64)


    thetacoarse = collect(range(0, stop = pi, length = xCoarse))
    thetafine = collect(range(0, stop = pi, length = ncell))
    quCoarseInter = Spline1D(theta, qu)
    quCoarse = evaluate(quCoarseInter, thetacoarse)
    quFineInt = Spline1D(thetacoarse, quCoarse)
    quFine = evaluate(quFineInt, thetafine)

    return quFine, thetafine

end

function FVMIBL(w::Array{Float64,2}, U::Array{Float64,1}, Ut::Array{Float64,1}, Ux::Array{Float64,1}, x::Array{Float64,1}, dt::Float64)

    n = length(x)

    dx = zeros(n)
    dx[2:end] = diff(x)
    dx[1] = dx[2]

    # correlate the unknown values from the del and E values
    del, E, FF ,B, S, dfde = correlate(w)

    fL, fR, UipL ,UipR, FFipL, FFipR, dfdeipL, dfdeipR, wipL, wipR = fluxReconstruction(w , U, FF, dfde, del, E)

    #lamb1L ,lamb2L = eigenlamb(UipL, dfdeipL, FFipL, wipL)
    #lamb1R ,lamb2R = eigenlamb(UipR, dfdeipR, FFipR, wipR)

    #lamb1 ,lamb2 = eigenlamb(U, dfde, FF, w)
    lamb1 ,lamb2 = calc_eigenjac(E, FF, dfde, U)
    dt = calc_Dt(lamb1 ,lamb2, 0.5, dx)

    #if t_cur + dt > t_tot
    #    dt = t_tot - t_cur
    #end

    # two step forward Euler methods for adding the source term to the right hand-side
    # of the transport equations.
    #z = RHSSource(U,B, del,Ut, Ux, FF, E, S)

    # step 1 : assuming this as a homogeneous equation and advanced half a step
    w1 = w.+ ((fL - fR).*((dt/2)./(dx)))

    del , E, FF ,B, S, dfde = correlate(w1)

    fL, fR, UipL ,UipR, FFipL, FFipR, dfdeipL, dfdeipR, wipL, wipR = fluxReconstruction(w1 , U, FF, dfde, del, E)

    z = RHSSourcejac(U, B, del, Ut, Ux, FF, E, S)

    #dtL = calc_Dt(UipL, dfdeipL, FFipL, wipL, 0.8, dx)
    #dtR = calc_Dt(UipR, dfdeipR, FFipR, wipR, 0.8, dx)

    #dt = calc_Dt(lamb1 ,lamb2, 0.6, dx)

    j1, j2 = separationJ(lamb1, lamb2, dt, dx)

    # step 2 : by considering source terms advanced a full step using 2nd order midpoint rule
    w2 = (w) + (fL - fR).* ((dt)./(dx)) .+ (dt).*z


    return w2, dt, lamb1, lamb2
end

function FVMIBLorig(w::Array{Float64,2}, U::Array{Float64,1}, Ut::Array{Float64,1}, Ux::Array{Float64,1}, x::Array{Float64,1}, t::Float64, t_tot::Float64)

    n = length(x)

    dx = Float64((x[end]-x[1])/n)

    csep = zeros(n)

    while t < t_tot
        # correlate the unknown values from the del and E values
        del, E, FF ,B, S, dfde = correlate(w)

        fL, fR, UipL ,UipR, FFipL, FFipR, dfdeipL, dfdeipR, wipL, wipR = fluxReconstruction(w , U, FF, dfde, del, E)

        #lamb1L ,lamb2L = eigenlamb(UipL, dfdeipL, FFipL, wipL)
        #lamb1R ,lamb2R = eigenlamb(UipR, dfdeipR, FFipR, wipR)

        #lamb1 ,lamb2 = eigenlamb(U, dfde, FF, w)
        lamb1 ,lamb2 = calc_eigen(E, FF, dfde, U)
        dt = calc_Dt(lamb1 ,lamb2, 0.5, ones(length(x)).*dx)

        if t + dt > t_tot
            dt = t_tot - t
        end

        #if t_cur + dt > t_tot
        #    dt = t_tot - t_cur
        #end

        # two step forward Euler methods for adding the source term to the right hand-side
        # of the transport equations.
        #z = RHSSource(U,B, del,Ut, Ux, FF, E, S)

        # step 1 : assuming this as a homogeneous equation and advanced half a step
        w1 = w.+ ((fL - fR).*((dt/2)./(dx)))

        del , E, FF ,B, S, dfde = correlate(w1)

        fL, fR, UipL ,UipR, FFipL, FFipR, dfdeipL, dfdeipR, wipL, wipR = fluxReconstruction(w1 , U, FF, dfde, del, E)

        z = RHSSource(U, B, del, Ut, Ux, FF, E, S)

        #dtL = calc_Dt(UipL, dfdeipL, FFipL, wipL, 0.8, dx)
        #dtR = calc_Dt(UipR, dfdeipR, FFipR, wipR, 0.8, dx)

        #dt = calc_Dt(lamb1 ,lamb2, 0.6, dx)

        j1, j2 = separationJ(lamb1, lamb2, dt, ones(length(x)).*dx)

        # step 2 : by considering source terms advanced a full step using 2nd order midpoint rule
        w2 = (w) + (fL - fR).* ((dt)./(dx)) .+ (dt).*z

        w[:,:] = w2[:,:]
        t += dt

        for i = 20:n-20
            csep[i] = (w[i+1,1] - w[i,1])/(x[i+1,1] - w[i,1])/pi
        end
        #println(t, "   ", maximum(csep))

    end

    return x, w, t
end


function update_indbound(surf::TwoDSurfThickBL, curfield::TwoDFlowField)
    surf.uind_u[1:surf.ndiv], surf.wind_u[1:surf.ndiv] = ind_vel([curfield.tev;
                                                                  curfield.lev; curfield.extv], surf.bnd_x_u, surf.bnd_z_u)

    surf.uind_l[1:surf.ndiv], surf.wind_l[1:surf.ndiv] = ind_vel([curfield.tev;
                                                                  curfield.lev; curfield.extv], surf.bnd_x_l, surf.bnd_z_l)

    return surf
end

function calc_Dt(lamb1::Array{Float64,1}, lamb2::Array{Float64,1}, cfl::Float64, dx::Array{Float64})

    # calculate time step values based on eigenvalues

    dti = cfl.*(dx./(abs.(lamb1+lamb2)))
    dt = minimum(abs.(dti))

    #@printf(" Max l1: %1.5f, Min l1: %1.5f, Max l2: %1.5f, Min l2: %1.5f \n", maximum(lamb1), minimum(lamb1),maximum(lamb2), minimum(lamb2));

    if dt< 0.00001
        println("Very low dt, singularity     , i_s=", argmin(dti))
        #dt =0.0005
    end

    return dt
end

function calc_Dtjac(lamb1, lamb2, cfl, dx)

    # calculate time step values based on eigenvalues

    dti = cfl.*(dx./(abs.(lamb1+lamb2)))
    dt = minimum(abs.(dti))

    return dt
end


function update_kinem(surf::TwoDSurfThickBL, t)

    # Pitch kinematics
    if (typeof(surf.kindef.alpha) == EldUpDef)
        surf.kinem.alpha = surf.kindef.alpha(t)
        surf.kinem.alphadot = ForwardDiff.derivative(surf.kindef.alpha,t)*surf.uref/surf.c
    elseif (typeof(surf.kindef.alpha) == EldUptstartDef)
        surf.kinem.alpha = surf.kindef.alpha(t)
        surf.kinem.alphadot = ForwardDiff.derivative(surf.kindef.alpha,t)*surf.uref/surf.c
    elseif (typeof(surf.kindef.alpha) == EldRampReturnDef)
        surf.kinem.alpha = surf.kindef.alpha(t)
        surf.kinem.alphadot = ForwardDiff.derivative(surf.kindef.alpha,t)*surf.uref/surf.c
    elseif (typeof(surf.kindef.alpha) == ConstDef)
        surf.kinem.alpha = surf.kindef.alpha(t)
        surf.kinem.alphadot = 0.
    elseif (typeof(surf.kindef.alpha) == SinDef)
        surf.kinem.alpha = surf.kindef.alpha(t)
        surf.kinem.alphadot = ForwardDiff.derivative(surf.kindef.alpha,t)*surf.uref/surf.c
    elseif (typeof(surf.kindef.alpha) == CosDef)
        surf.kinem.alpha = surf.kindef.alpha(t)
        surf.kinem.alphadot = ForwardDiff.derivative(surf.kindef.alpha,t)*surf.uref/surf.c
    elseif (typeof(surf.kindef.alpha) == VAWTalphaDef)
        surf.kinem.alpha = surf.kindef.alpha(t)
        surf.kinem.alphadot = ForwardDiff.derivative(surf.kindef.alpha,t)*surf.uref/surf.c
    end
    # ---------------------------------------------------------------------------------------------

    # Plunge kinematics
    if (typeof(surf.kindef.h) == EldUpDef)
        surf.kinem.h = surf.kindef.h(t)*surf.c
        surf.kinem.hdot = ForwardDiff.derivative(surf.kindef.h,t)*surf.uref
    elseif (typeof(surf.kindef.h) == EldUptstartDef)
        surf.kinem.h = surf.kindef.h(t)*surf.c
        surf.kinem.hdot = ForwardDiff.derivative(surf.kindef.h,t)*surf.uref
    elseif (typeof(surf.kindef.h) == EldUpIntDef)
        surf.kinem.h = surf.kindef.h(t)*surf.c
        surf.kinem.hdot = ForwardDiff.derivative(surf.kindef.h,t)*surf.uref
    elseif (typeof(surf.kindef.h) == EldUpInttstartDef)
        surf.kinem.h = surf.kindef.h(t)*surf.c
        surf.kinem.hdot = ForwardDiff.derivative(surf.kindef.h,t)*surf.uref
    elseif (typeof(surf.kindef.h) == EldRampReturnDef)
        surf.kinem.h = surf.kindef.h(t)*surf.c
        surf.kinem.hdot = ForwardDiff.derivative(surf.kindef.h,t)*surf.uref
    elseif (typeof(surf.kindef.h) == ConstDef)
        surf.kinem.h = surf.kindef.h(t)*surf.c
        surf.kinem.hdot = 0.
    elseif (typeof(surf.kindef.h) == SinDef)
        surf.kinem.h = surf.kindef.h(t)*surf.c
        surf.kinem.hdot = ForwardDiff.derivative(surf.kindef.h,t)*surf.uref
    elseif (typeof(surf.kindef.h) == CosDef)
        surf.kinem.h = surf.kindef.h(t)*surf.c
        surf.kinem.hdot =  ForwardDiff.derivative(surf.kindef.h,t)*surf.uref
    elseif (typeof(surf.kindef.h) == VAWThDef)
        surf.kinem.h = surf.kindef.h(t)*surf.c
        surf.kinem.hdot = ForwardDiff.derivative(surf.kindef.h,t)*surf.uref
    end
    # ---------------------------------------------------------------------------------------------

    # Forward velocity
    if (typeof(surf.kindef.u) == EldUpDef)
        surf.kinem.u = surf.kindef.u(t)*surf.uref
        surf.kinem.udot = ForwardDiff.derivative(surf.kindef.u,t)*surf.uref*surf.uref/surf.c
    elseif (typeof(surf.kindef.u) == EldRampReturnDef)
        surf.kinem.u, surf.kinem.udot = surf.kindef.u(t)
        surf.kinem.u = surf.kinem.u*surf.uref
        surf.kinem.udot = surf.kinem.udot*surf.uref*surf.uref/surf.c
    elseif (typeof(surf.kindef.u) == ConstDef)
        surf.kinem.u = surf.kindef.u(t)*surf.uref
        surf.kinem.udot = 0.
    elseif (typeof(surf.kindef.u) == SinDef)
        surf.kinem.u = surf.kindef.u(t)*surf.uref
        surf.kinem.udot = ForwardDiff.derivative(surf.kindef.u,t)*surf.uref*surf.uref/surf.c
    elseif (typeof(surf.kindef.u) == CosDef)
        surf.kinem.u = surf.kindef.u(t)*surf.uref
        surf.kinem.udot = ForwardDiff.derivative(surf.kindef.u,t)*surf.uref*surf.uref/surf.c
    elseif (typeof(surf.kindef.u) == LinearDef)
        surf.kinem.u = surf.kindef.u(t)*surf.uref
        surf.kinem.udot = ForwardDiff.derivative(surf.kindef.u,t)*surf.uref*surf.uref/surf.c
    elseif (typeof(surf.kindef.u) == VAWTuDef)
        surf.kinem.u = surf.kindef.u(t)*surf.uref
        surf.kinem.udot = ForwardDiff.derivative(surf.kindef.u,t)*surf.uref*surf.uref/surf.c
    end
    # ---------------------------------------------------------------------------------------------
    return surf
end

function wakeroll(surf::TwoDSurfThickBL, curfield::TwoDFlowField, dt)

    nlev = length(curfield.lev)
    ntev = length(curfield.tev)
    nextv = length(curfield.extv)

    #Clean induced velocities
    for i = 1:ntev
        curfield.tev[i].vx = 0
        curfield.tev[i].vz = 0
    end

    for i = 1:nlev
        curfield.lev[i].vx = 0
        curfield.lev[i].vz = 0
    end

    for i = 1:nextv
        curfield.extv[i].vx = 0
        curfield.extv[i].vz = 0
    end

    #Velocities induced by free vortices on each other
    mutual_ind([curfield.tev; curfield.lev; curfield.extv])

    #Add the influence of velocities induced by bound vortices
    utemp = zeros(ntev + nlev + nextv)
    wtemp = zeros(ntev + nlev + nextv)
    utemp, wtemp = ind_vel(surf.bv, [map(q -> q.x, curfield.tev); map(q -> q.x, curfield.lev); map(q -> q.x, curfield.extv)], [map(q -> q.z, curfield.tev); map(q -> q.z, curfield.lev); map(q -> q.z, curfield.extv) ])

    for i = 1:ntev
        curfield.tev[i].vx += utemp[i]
        curfield.tev[i].vz += wtemp[i]
    end
    for i = ntev+1:ntev+nlev
        curfield.lev[i-ntev].vx += utemp[i]
        curfield.lev[i-ntev].vz += wtemp[i]
    end
    for i = ntev+nlev+1:ntev+nlev+nextv
        curfield.extv[i-ntev-nlev].vx += utemp[i]
        curfield.extv[i-ntev-nlev].vz += wtemp[i]
    end

    #Add effect of sources
    utemp, wtemp = ind_vel_src(surf.src, [map(q -> q.x, curfield.tev);
                                          map(q -> q.x, curfield.lev); map(q -> q.x, curfield.extv)],
                               [map(q -> q.z, curfield.tev); map(q -> q.z, curfield.lev);
                                map(q -> q.z, curfield.extv)])

    for i = 1:ntev
        curfield.tev[i].vx += utemp[i]
        curfield.tev[i].vz += wtemp[i]
    end
    for i = ntev+1:ntev+nlev
        curfield.lev[i-ntev].vx += utemp[i]
        curfield.lev[i-ntev].vz += wtemp[i]
    end
    for i = ntev+nlev+1:ntev+nlev+nextv
        curfield.extv[i-ntev-nlev].vx += utemp[i]
        curfield.extv[i-ntev-nlev].vz += wtemp[i]
    end

    #Convect free vortices with their induced velocities
    for i = 1:ntev
        curfield.tev[i].x += dt*curfield.tev[i].vx
        curfield.tev[i].z += dt*curfield.tev[i].vz
    end
    for i = 1:nlev
        curfield.lev[i].x += dt*curfield.lev[i].vx
        curfield.lev[i].z += dt*curfield.lev[i].vz
    end
    for i = 1:nextv
        curfield.extv[i].x += dt*curfield.extv[i].vx
        curfield.extv[i].z += dt*curfield.extv[i].vz
    end

    return curfield
end

function update_boundpos(surf::TwoDSurfThickBL, dt::Float64)
    for i = 1:surf.ndiv
        surf.bnd_x_u[i] = surf.bnd_x_u[i] + dt*((surf.pvt*surf.c - surf.x[i])*
                                                sin(surf.kinem.alpha)*surf.kinem.alphadot - surf.kinem.u + (surf.cam[i] +
                                                                                                            surf.thick[i])*cos(surf.kinem.alpha)*surf.kinem.alphadot)

        surf.bnd_z_u[i] = surf.bnd_z_u[i] + dt*(surf.kinem.hdot + (surf.pvt*surf.c -
                                                                   surf.x[i])*cos(surf.kinem.alpha)*surf.kinem.alphadot - (surf.cam[i] +
                                                                                                                           surf.thick[i])*sin(surf.kinem.alpha)*surf.kinem.alphadot)

        surf.bnd_x_l[i] = surf.bnd_x_l[i] + dt*((surf.pvt*surf.c - surf.x[i])*
                                                sin(surf.kinem.alpha)*surf.kinem.alphadot - surf.kinem.u + (surf.cam[i] -
                                                                                                            surf.thick[i])*cos(surf.kinem.alpha)*surf.kinem.alphadot)

        surf.bnd_z_l[i] = surf.bnd_z_l[i] + dt*(surf.kinem.hdot + (surf.pvt*
                                                                   surf.c - surf.x[i])*cos(surf.kinem.alpha)*surf.kinem.alphadot - (surf.cam[i] -
                                                                                                                                    surf.thick[i])*sin(surf.kinem.alpha)*surf.kinem.alphadot)

        surf.bnd_x_chord[i] = surf.bnd_x_chord[i] + dt*((surf.pvt*surf.c - surf.x[i])*
                                                        sin(surf.kinem.alpha)*surf.kinem.alphadot - surf.kinem.u)
        surf.bnd_z_chord[i] = surf.bnd_z_chord[i] + dt*(surf.kinem.hdot + (surf.pvt*surf.c -
                                                                           surf.x[i])*cos(surf.kinem.alpha)*surf.kinem.alphadot)
    end

    return surf
end


function update_thickLHS(surf::TwoDSurfThickBL, curfield::TwoDFlowField, dt::Float64, vcore::Float64)
    #Calculate the missing column in LHS that depends on last shed vortex location

    ntev = length(curfield.tev)

    if ntev == 0
        xloc_tev = surf.bnd_x_chord[surf.ndiv] + 0.5*surf.kinem.u*dt
        zloc_tev = surf.bnd_z_chord[surf.ndiv]
    else
        xloc_tev = surf.bnd_x_chord[surf.ndiv] + (1. /3.)*(curfield.tev[ntev].x - surf.bnd_x_chord[surf.ndiv])
        zloc_tev = surf.bnd_z_chord[surf.ndiv] + (1. /3.)*(curfield.tev[ntev].z - surf.bnd_z_chord[surf.ndiv])
    end

    dummyvort = TwoDVort(xloc_tev, zloc_tev, 1., vcore, 0., 0.)

    uu, wu = ind_vel([dummyvort], surf.bnd_x_u, surf.bnd_z_u)
    ul, wl = ind_vel([dummyvort], surf.bnd_x_l, surf.bnd_z_l)

    for i = 2:surf.ndiv-1


        #Sweep all rows (corresponding to ndiv) for lifting equation

        #Sweep columns for aterms
        for n = 1:surf.naterm
           surf.LHS[i-1,n] = cos(n*surf.theta[i]) - surf.thick_slope[i]*sin(n*surf.theta[i])
        end

        #Sweep columns for bterm
        for n = 1:surf.naterm
           surf.LHS[i-1,n+surf.naterm] = surf.cam_slope[i]*cos(n*surf.theta[i])
        end

        #TEV term must be updated in the loop after its location is known
        #Sweep all rows (corresponding to ndiv) for nonlifting equation

        for n = 1:surf.naterm
           surf.LHS[surf.ndiv+i-3,n]  = -surf.cam_slope[i]*sin(n*surf.theta[i])
        end
        for n = 1:surf.naterm
           surf.LHS[surf.ndiv+i-3,surf.naterm+n] = sin(n*surf.theta[i]) + surf.thick_slope[i]*cos(n*surf.theta[i])
        end


        wlz = 0.5*(wu[i]*cos(surf.kinem.alpha) + uu[i]*sin(surf.kinem.alpha) +
                   wl[i]*cos(surf.kinem.alpha) + ul[i]*sin(surf.kinem.alpha))

        wtz = 0.5*(wu[i]*cos(surf.kinem.alpha) + uu[i]*sin(surf.kinem.alpha) -
                   wl[i]*cos(surf.kinem.alpha) - ul[i]*sin(surf.kinem.alpha))

        wtx = 0.5*(uu[i]*cos(surf.kinem.alpha) - wu[i]*sin(surf.kinem.alpha) +
                   ul[i]*cos(surf.kinem.alpha) - wl[i]*sin(surf.kinem.alpha))

        wlx = 0.5*(uu[i]*cos(surf.kinem.alpha) - wu[i]*sin(surf.kinem.alpha) -
                   ul[i]*cos(surf.kinem.alpha) + wl[i]*sin(surf.kinem.alpha))

        surf.LHS[i-1,1+2*surf.naterm] = -surf.cam_slope[i]*wtx -
            surf.thick_slope[i]*wlx + wlz

        surf.LHS[surf.ndiv+i-3,1+2*surf.naterm] = -surf.cam_slope[i]*wlx -
            surf.thick_slope[i]*wtx + wtz

    end

    #surf.LHS[2*surf.ndiv-1,2+2*surf.naterm] = -wu[1]*cos(surf.kinem.alpha) -
    #uu[1]*sin(surf.kinem.alpha)


    return surf, xloc_tev, zloc_tev
end

function update_thickLHS2V(surf::TwoDSurfThickBL, curfield::TwoDFlowField, dt::Float64, vcore::Float64, i_xsl)
    #Calculate the missing column in LHS that depends on last shed vortex location

    ntev = length(curfield.tev)

    vx = surf.qu[i_xsl]/sqrt(1. + (-surf.thick_slope[i_xsl] + surf.cam_slope[i_xsl])^2)
    vz = surf.qu[i_xsl]*(-surf.thick_slope[i_xsl] + surf.cam_slope[i_xsl])/sqrt(1. + (-surf.thick_slope[i_xsl] + surf.cam_slope[i_xsl])^2)
    
    alpha = surf.kinem.alpha
    R = [cos(alpha) -sin(alpha); sin(alpha) cos(alpha)]
    
    vx1, vz1 = R*[vx; 0]
    vx2, vz2 = R*[0; vz]
    vx = vx1 + vx2
    vz = vz1 + vz2

    
    if ntev == 0
        xloc_tev = surf.bnd_x_l[i_xsl] + 0.5*vx*dt
        zloc_tev = surf.bnd_z_l[i_xsl] + 0.5*vz*dt
    else
        xloc_tev = surf.bnd_x_l[i_xsl] + (1. /3.)*(curfield.tev[ntev].x - surf.bnd_x_l[i_xsl])
        zloc_tev = surf.bnd_z_l[i_xsl] + (1. /3.)*(curfield.tev[ntev].z - surf.bnd_z_l[i_xsl])
    end

    dummylsv = TwoDVort(xloc_tev, zloc_tev, 1., vcore, 0., 0.)
    #dummylsv = TwoDVort(xloc_tev, zloc_tev, 1., vcore, 0., 0.)

    uu, wu = ind_vel([dummylsv], surf.bnd_x_u, surf.bnd_z_u)
    ul, wl = ind_vel([dummylsv], surf.bnd_x_l, surf.bnd_z_l)

    for i = 2:surf.ndiv-1

        
        #Sweep all rows (corresponding to ndiv) for lifting equation
        
        # #Sweep columns for aterms
        # for n = 1:surf.naterm
        #     surf.LHS[i-1,n] = cos(n*surf.theta[i]) - surf.thick_slope[i]*sin(n*surf.theta[i])
        # end

        # #Sweep columns for bterm
        # for n = 1:surf.naterm
        #    surf.LHS[i-1,n+surf.naterm] = surf.cam_slope[i]*cos(n*surf.theta[i])
        # end

        # #TEV term must be updated in the loop after its location is known
        # #Sweep all rows (corresponding to ndiv) for nonlifting equation

        # for n = 1:surf.naterm
        #    surf.LHS[surf.ndiv+i-3,n]  = -surf.cam_slope[i]*sin(n*surf.theta[i])
        # end
        # for n = 1:surf.naterm
        #    surf.LHS[surf.ndiv+i-3,surf.naterm+n] = sin(n*surf.theta[i]) + surf.thick_slope[i]*cos(n*surf.theta[i])
        # end


        wlz = 0.5*(wu[i]*cos(surf.kinem.alpha) + uu[i]*sin(surf.kinem.alpha) +
                   wl[i]*cos(surf.kinem.alpha) + ul[i]*sin(surf.kinem.alpha))

        wtz = 0.5*(wu[i]*cos(surf.kinem.alpha) + uu[i]*sin(surf.kinem.alpha) -
                   wl[i]*cos(surf.kinem.alpha) - ul[i]*sin(surf.kinem.alpha))

        wtx = 0.5*(uu[i]*cos(surf.kinem.alpha) - wu[i]*sin(surf.kinem.alpha) +
                   ul[i]*cos(surf.kinem.alpha) - wl[i]*sin(surf.kinem.alpha))

        wlx = 0.5*(uu[i]*cos(surf.kinem.alpha) - wu[i]*sin(surf.kinem.alpha) -
                   ul[i]*cos(surf.kinem.alpha) + wl[i]*sin(surf.kinem.alpha))

        surf.LHS[i-1,1+2*surf.naterm] = -surf.cam_slope[i]*wtx -
            surf.thick_slope[i]*wlx + wlz

        surf.LHS[surf.ndiv+i-3,1+2*surf.naterm] = -surf.cam_slope[i]*wlx -
            surf.thick_slope[i]*wtx + wtz

    end

    return surf, xloc_tev, zloc_tev
end

function update_thickRHS(surf::TwoDSurfThickBL, curfield::TwoDFlowField)
    for i = 2:surf.ndiv-1
        #RHS for lifting equation

        wlz = 0.5*(surf.wind_u[i]*cos(surf.kinem.alpha) + surf.uind_u[i]*sin(surf.kinem.alpha) +
                   surf.wind_l[i]*cos(surf.kinem.alpha) + surf.uind_l[i]*sin(surf.kinem.alpha))

        wtz = 0.5*(surf.wind_u[i]*cos(surf.kinem.alpha) + surf.uind_u[i]*sin(surf.kinem.alpha) -
                   surf.wind_l[i]*cos(surf.kinem.alpha) - surf.uind_l[i]*sin(surf.kinem.alpha))

        wtx = 0.5*(surf.uind_u[i]*cos(surf.kinem.alpha) - surf.wind_u[i]*sin(surf.kinem.alpha) +
                   surf.uind_l[i]*cos(surf.kinem.alpha) - surf.wind_l[i]*sin(surf.kinem.alpha))

        wlx = 0.5*(surf.uind_u[i]*cos(surf.kinem.alpha) - surf.wind_u[i]*sin(surf.kinem.alpha) -
                   surf.uind_l[i]*cos(surf.kinem.alpha) + surf.wind_l[i]*sin(surf.kinem.alpha))

        surf.RHS[i-1] = -(surf.kinem.u + curfield.u[1])*sin(surf.kinem.alpha)/surf.uref -
            surf.kinem.alphadot*(surf.x[i] - surf.pvt*surf.c)/surf.uref +
            (surf.kinem.hdot - curfield.w[1])*cos(surf.kinem.alpha)/surf.uref - wlz/surf.uref +
            surf.cam_slope[i]*((surf.kinem.u + curfield.u[1])*cos(surf.kinem.alpha) +
                               (surf.kinem.hdot - curfield.w[1])*sin(surf.kinem.alpha) + wtx -
                               surf.kinem.alphadot*surf.cam[i]) +
                               surf.thick_slope[i]*(wlx - surf.kinem.alphadot*surf.thick[i])

        surf.RHS[surf.ndiv+i-3] = surf.cam_slope[i]*(wlx - surf.kinem.alphadot*surf.thick[i]) +
            surf.thick_slope[i]*((surf.kinem.u + curfield.u[1])*cos(surf.kinem.alpha) +
                                 (surf.kinem.hdot - curfield.w[1])*sin(surf.kinem.alpha) + wtx -
                                 surf.kinem.alphadot*surf.cam[i]) - wtz
    end

    #RHS for Kelvin condition (negative strength of all previously shed vortices)
    surf.RHS[2*surf.ndiv-3] = -100*(sum(map(q->q.s, curfield.tev)) + sum(map(q->q.s, curfield.lev)))/(surf.uref*surf.c)

    #Kutta condition
    surf.RHS[2*surf.ndiv-2] = 0.
    #surf.RHS[2*surf.ndiv-3] = pi*(surf.a0prev[1] + 0.5*surf.aprev[1])

        #RHS for LE velocity condition
    #surf.RHS[2*surf.ndiv-1] = sin(surf.kinem.alpha) - surf.kinem.hdot*cos(surf.kinem.alpha)/surf.uref -
    #surf.kinem.alphadot*surf.pvt*surf.c/surf.uref +
    #   (surf.wind_u[1]*cos(surf.kinem.alpha) + surf.uind_u[1]*sin(surf.kinem.alpha))/surf.uref


   # surf.RHS[2*surf.ndiv-2] = -(surf.kinem.u + curfield.u[1])*cos(surf.kinem.alpha) - (surf.kinem.hdot - curfield.w[1])*sin(surf.kinem.alpha)

    #surf.RHS[2*surf.ndiv-1] = -(surf.kinem.u + curfield.u[1])*cos(surf.kinem.alpha) - (surf.kinem.hdot - curfield.w[1])*sin(surf.kinem.alpha)

end


function fluxReconstruction(w::Array{Float64,2}, U::Array{Float64,1}, FF::Array{Float64,1}, dfde::Array{Float64,1}, del::Array{Float64,1} , E::Array{Float64,1})

    delF = U.*(w[:,2] - w[:,1])
    EF = U.* FF.* w[:,1]
    F = hcat(delF, EF)
    # first-order approximation of left and right side of the i+1/2
    # cell interface

    # 1) approximation of edge velocity
    UR = U
    UL = U

    #UR = limiter(UR, U, [U[2:end]; U[1]])
    #UL = limiter(UL, U, [U[end]; U[1:end-1]])

    #Why is this kind of cyclic condition used?

    UipL = UR
    UipR = [UL[2:end]; UL[1]]



    # 2) approximation of FF
    FFR = FF
    FFL = FF
    FFipL = FFR
    FFipR = [FFL[2:end]; FFL[1]]

    # 3) approximation of dfde
    dfdeR = dfde
    dfdeL = dfde
    dfdeipL = dfdeR
    dfdeipR = [dfdeL[2:end]; dfdeL[1]]

    # 4) approximation w (a two-dimensional array)
    wR = w
    wL = w
    wipL = wR
    wipR = [wL[2:end,:]; wL[1:1,:]]

    # calculating left and right flux based on calculated left and right quantities
    delFluxL = UipL.*(wipL[:,2] - wipL[:,1])
    delFluxR = UipR.*(wipR[:,2] - wipR[:,1])
    EFluxL = UipL.* FFipL.* wipL[:,1]
    EFluxR = UipR.* FFipR.* wipR[:,1]

    fpL = hcat(delFluxL, EFluxL)
    fpR = hcat(delFluxR, EFluxR)

    # calculating wave speed at the left and right interaces
    wsL = maxWaveSpeed(UipL, wipL, dfdeipL, FFipL)
    wsR = maxWaveSpeed(UipR, wipR, dfdeipR, FFipR)
    # selecting the maximum wave speed
    ws = max.(wsL,wsR)
    #ws = max.(wsL+wsR)
    ww = hcat(ws,ws)
    
    # flux reconstruction of left and right side of the i+1/2 interface
    fR = 0.5*((fpL + fpR) + ww.* (wipL - wipR))
    fL = [fR[end:end,:];fR[1:end-1,:]]

    # specifying the boundary conditions using internal extrapolation from the calculated flux of the
    # neighbouring cell centers
    #fL[1,:] = [F[1,1]; F[1,2]]
    fL[1,:] = 0.5*((F[1,:]) - wsR[1,:].* (wipR[1,:]))
    #fL[1,:] = [0; 0]
    
    #fR[end,:] = ((F[end,:]) - wsR[end,:].* (wipR[end,:]))
    fR[end,:] = [F[end,1];F[end,2]]
    #fL[end,:] .= 0

    #fR[end,:] = [0.0;0.0]
    

    return fL, fR, UipL ,UipR, FFipL, FFipR, dfdeipL, dfdeipR, wipL, wipR

end


function limiter(wExtrapolated::Array{Float64,1}, wCell::Array{Float64,1}, wNeighbor::Array{Float64,1})

    w = wExtrapolated;
    wMax = max.(wCell, wNeighbor);
    wMin = min.(wCell, wNeighbor);
    w[w.>wMax] = wMax[w.>wMax];
    w[w.<wMin] = wMin[w.<wMin];

    return w

end


function maxWaveSpeed(Uip::Array{Float64,1}, wip::Array{Float64,2}, dfdeip::Array{Float64,1}, FFip::Array{Float64,1})

    # calculate wave speed at the interfaces of the cell
    wsP = abs.(0.5*Uip).*((dfdeip .- 1.0)
        + sqrt.(1.0 .+ 4.0*FFip .- 2.0.* dfdeip .- (4.0*((wip[:,2]./wip[:,1]) .-1.0).*dfdeip) .+ dfdeip.^2))

    wsN = abs.(0.5*Uip).*((dfdeip .- 1.0)
            - sqrt.(1.0 .+ 4.0*FFip .- 2.0.* dfdeip .- (4.0*((wip[:,2]./wip[:,1]) .-1.0).*dfdeip) .+ dfdeip.^2))

    return max(wsP, wsN)
end


function calc_eigen(E::Array{Float64}, F::Array{Float64},
                    dfde::Array{Float64}, ue::Array{Float64})

    ncell = length(E)
    lamb1 = zeros(ncell); lamb2 = zeros(ncell)

    for i = 1:ncell
        a_q = 1.
        b_q = -ue[i] * (dfde[i] - 1)
        c_q = ue[i] * ue[i] * (E[i]*dfde[i] - F[i])
        lamb1[i] = (-b_q + sqrt(b_q*b_q - 4*a_q*c_q))/(2*a_q)
        lamb2[i] = (-b_q  -sqrt(b_q*b_q - 4*a_q*c_q))/(2*a_q)

        #Always have lamb1 > lamb2
        if lamb2[i] > lamb1[i]
            temp  = lamb2[i]
            lamb2[i] = lamb1[i]
            lamb1[i] = temp
        end
    end

    return lamb1, lamb2
end

function RHSSource(U::Array{Float64,1} ,B::Array{Float64,1}, del::Array{Float64,1},Ut::Array{Float64,1}, Ux::Array{Float64,1}, FF::Array{Float64,1}, E::Array{Float64,1}, S::Array{Float64,1} )

    z1 = B./(2.0*del) .- del.* (Ut./U) .- (E.+ 1.0).*del.*Ux
    z2 = S./del .- 2.0*E.*del.* (Ut./U) .- 2.0*FF.*del.*Ux
    z = hcat(z1,z2)

    return z

end


function separationJ(lamb1::Array{Float64,1}, lamb2::Array{Float64,1}, dt::Float64, dx::Array{Float64,1})

    N1 = length(lamb1)
    lambj1 = sum(lamb1)/N1
    J1Sep = (dt)./(dx[2:end]) .* (lamb1[1:end-1] - lamb1[2:end])./ (lambj1*N1)

    for i = 2:N1-1
        if lamb1[i] .< 0.0
            J1Sep[i] =  (dt)/(dx[i]) .* (lamb1[i+1] - lamb1[i-1])./ (lambj1*N1)
        end
    end

    N2 = length(lamb2)
    lambj2 = sum(lamb2)/N2
    J2Sep = (dt)/(dx[2:end]) .* (lamb2[1:end-1] - lamb2[2:end])./ (lambj2*N2)
    
    for i = 2:N2-1
        if lamb2[i] .< 0.0
            J2Sep[i] =  (dt)/(dx[i]) * (lamb2[i+1] - lamb2[i-1]) / (lambj2*N2)
        end
    end

    return J1Sep, J2Sep

end

function update_atermdot(surf::TwoDSurfThickBL,dt)
    surf.a0dot[1] = (surf.a0[1] - surf.a0prev[1])/dt
    for ia = 1:surf.naterm
        surf.adot[ia] = (surf.aterm[ia]-surf.aprev[ia])/dt
    end
    return surf
end

function update_bv_src(surf::TwoDSurfThickBL)
    gamma = zeros(surf.ndiv)
    src = zeros(surf.ndiv)
    for ib = 1:surf.ndiv
        gamma[ib] = surf.a0[1]*(1 + cos(surf.theta[ib]))
        for ia = 1:surf.naterm
            gamma[ib] += surf.aterm[ia]*sin(ia*surf.theta[ib])*sin(surf.theta[ib])
            src[ib] += surf.bterm[ia]*sin(ia*surf.theta[ib])*sin(surf.theta[ib])
        end
        gamma[ib] = gamma[ib]*surf.uref*surf.c
        src[ib] = src[ib]*surf.uref*surf.c
    end

    for ib = 2:surf.ndiv
        surf.bv[ib-1].s = (gamma[ib]+gamma[ib-1])*(surf.theta[2]-surf.theta[1])/2.
        surf.bv[ib-1].x = (surf.bnd_x_chord[ib] + surf.bnd_x_chord[ib-1])/2.
        surf.bv[ib-1].z = (surf.bnd_z_chord[ib] + surf.bnd_z_chord[ib-1])/2.
        surf.src[ib-1].s = (src[ib] + src[ib-1])*(surf.theta[2]-surf.theta[1])/2.
        surf.src[ib-1].x = surf.bv[ib-1].x
        surf.src[ib-1].z = surf.bv[ib-1].z
    end

    for ib = surf.ndiv:2*surf.ndiv-2
        cn = ib-surf.ndiv+1
        surf.src[ib].s = (surf.wtu[cn] + surf.wtu[cn+1])/2*(surf.su[cn+1] - surf.su[cn])
        surf.src[ib].x = 0.5*(surf.bnd_x_u[cn] + surf.bnd_x_u[cn+1])
        surf.src[ib].z = 0.5*(surf.bnd_z_u[cn] + surf.bnd_z_u[cn+1])
    end
    for ib = 2*surf.ndiv-1:3*surf.ndiv-3
        cn = ib-2*surf.ndiv+2
        surf.src[ib].s = (surf.wtl[cn] + surf.wtl[cn+1])/2*(surf.sl[cn+1] - surf.sl[cn])
        surf.src[ib].x = 0.5*(surf.bnd_x_l[cn] + surf.bnd_x_l[cn+1])
        surf.src[ib].z = 0.5*(surf.bnd_z_l[cn] + surf.bnd_z_l[cn+1])
    end
    
end


function smoothEdges!(q::Array{Float64}, nsm::Int)
    ndiv = length(q)
    
    for i = 1:nsm
        q[nsm-i+1] = 2*q[nsm+2-i] - q[nsm+3-i]
        q[ndiv-nsm+i] = 2*q[ndiv-nsm+i-1] - q[ndiv-nsm+i-2]
    end

    return q
end



function transResidual(x, naterm, uref, theta, xtev, ztev, vctev, bnd_x_u, bnd_z_u, bnd_x_l, bnd_z_l, uind_u, wind_u, uind_l, wind_l, vels, alpha, alphadot, u, hdot, cam, thick, cam_slope, thick_slope, pvt, c, su, delstart, Estart, LHS, RHS, qucprev, dt, t, Re)

    ndiv = length(theta)

    aterm = x[1:naterm]
    bterm = x[naterm+1:2*naterm]
    stev = x[2*naterm+1:2*naterm+1]
    
    quc = zeros(ndiv-1)
    quxc = zeros(ndiv-1)
    qutc = zeros(ndiv-1)
    dsu = zeros(ndiv-1)
    suc = zeros(ndiv-1)

    
    wtu = zeros(ndiv)
    wtl = zeros(ndiv)

    F = zeros(2*naterm+1)
    
    #println(typeof(aterm))
    #println(typeof(stev))

    # surf.aterm[:] = x[1:surf.naterm]
    # surf.bterm[:] = x[surf.naterm+1:2*surf.naterm]
    # #curfield.tev[end].s = x[2*surf.naterm+1]
    # #tev hasnt been added yet
    # del = x[2*surf.naterm+1:2*surf.naterm+surf.nfvm]
    # E = x[2*surf.naterm+surf.nfvm+1:2*surf.naterm+2*surf.nfvm]

    qu, ql = calc_edgeVelIBL(uref, theta, xtev, ztev, stev, vctev, bnd_x_u, bnd_z_u, bnd_x_l, bnd_z_l, aterm, bterm, uind_u, wind_u, uind_l, wind_l, vels, alpha, alphadot, u, hdot, cam, thick, cam_slope, thick_slope, pvt, c)

    #Solve the FV problem at cell centres
    for i = 1:ndiv-1
        quc[i] = (qu[i] + qu[i+1])/2
        suc[i] = (su[i] + su[i+1])/2
    end
    
    quxc[2:end] = diff(quc)./diff(suc)
    quxc[1] = 2*quxc[2] - quxc[3]
    
    #smoothEdges!(qux, 10)

    qutc[:] .= (quc[:] .- qucprev[:])./dt

    dsu = diff(su)
    
    w0 = [delstart delstart.*(Estart .+ 1)]
    
    w0 = FVMIBLgridfixed(w0, quc, qutc, quxc, dsu, dt)
    
    delu = w0[:,1]
    Eu = w0[:,2]./w0[:,1] .- 1.

    delu = [delu; delu[end]]
    Eu = [Eu; Eu[end]]
    #smoothEdges!(delu, 5)
    
    xf = c/2 .* (1. .- cos.(theta))

    wtu[2:end] = (1/sqrt(Re))*diff(qu.*delu)./diff(xf)
    wtu[1] = wtu[2]

    wtl[:] = -wtu[:]
    
    #smoothEdges!(wtu, 5)

    RHStransp = zeros(2*ndiv-2)

    for i = 2:ndiv-1
        RHStransp[i-1] = 0.5*(wtu[i] + wtl[i])
        RHStransp[ndiv+i-3] = 0.5*(wtu[i] - wtl[i])
    end

    figure(1)
    plot(xf, wtu)

    x_solved = LHS[1:ndiv*2-2, 1:naterm*2+1] \ (RHS[1:ndiv*2-2] + RHStransp[:])
    
    #Residual for inviscid problem
    #F[1:2*ndiv-2] = LHS[1:ndiv*2-2, 1:naterm*2+1]*x[1:2*naterm+1] - (RHS[1:ndiv*2-2] + RHStransp[:])
    F[1:2*naterm+1] = abs.(x_solved .- x[1:2*naterm+1])
    #Residual for viscous problem
    println(sum(wtu))
    #F[2*naterm+2:end] = abs.([delu[1:end-1]; Eu[1:end-1]] .- x[2*naterm+2:end])

    return F
end

function FVMIBLgridfixed(w, U, Ut, Ux, dx, dt)


    #while t < t_tot
        # correlate the unknown values from the del and E values
    del, E, FF ,B, S, dfde = correlate(w)

    fL, fR, UipL ,UipR, FFipL, FFipR, dfdeipL, dfdeipR, wipL, wipR = fluxReconstruction(w,  U, FF, dfde, del, E)

        #lamb1 ,lamb2 = calc_eigenjac(E, FF, dfde, U)
        #dt = calc_Dtjac(lamb1 ,lamb2, 0.5, ones(length(x)).*dx)

        #if t + dt > t_tot
     #       dt = t_tot - t
      #  end

    # step 1 : assuming this as a homogeneous equation and advanced half a step
    w1 = w .+ ((fL - fR).*((dt/2)./(dx)))

    del , E, FF ,B, S, dfde = correlate(w1)

    fL, fR, UipL ,UipR, FFipL, FFipR, dfdeipL, dfdeipR, wipL, wipR = fluxReconstruction(w1, U, FF, dfde, del, E)

    z = RHSSource(U, B, del, Ut, Ux, FF, E, S)
    
    # step 2 : by considering source terms advanced a full step using 2nd order midpoint rule
    w2 = w + (fL - fR).* ((dt)./(dx)) .+ (dt).*z
    
    # t += dt
    
    # for i = 20:n-20
    #     csep[i] = (wfin[i+1,1] - w[i,1])/(x[i+1,1] - w[i,1])/pi
    # end

  #  end

    return w2
end

function FVMIBLgridvar(w, U, Ut, Ux, dx, t, t_tot)

    n = length(dx)

    i_s = n

    j1 = zeros(n-1)
    j2 = zeros(n-1)
    
    while t < t_tot

        # correlate the unknown values from the del and E values
        del, E, FF ,B, S, dfde = correlate(w)
        
        fL, fR, UipL ,UipR, FFipL, FFipR, dfdeipL, dfdeipR, wipL, wipR = fluxReconstruction(w,  U, FF, dfde, del, E)
        
        lamb1 ,lamb2 = calc_eigen(E, FF, dfde, U)
        dt = calc_Dt(lamb1 ,lamb2, 0.5, dx)
        
        if t + dt > t_tot
            dt = t_tot - t
        end
        
        # step 1 : assuming this as a homogeneous equation and advanced half a step
        w1 = w .+ ((fL - fR).*((dt/2)./(dx)))
        
        del , E, FF ,B, S, dfde = correlate(w1)
        
        fL, fR, UipL ,UipR, FFipL, FFipR, dfdeipL, dfdeipR, wipL, wipR = fluxReconstruction(w1, U, FF, dfde, del, E)
        
        z = RHSSource(U, B, del, Ut, Ux, FF, E, S)
        
        # step 2 : by considering source terms advanced a full step using 2nd order midpoint rule
        w = w + (fL - fR).* ((dt)./(dx)) .+ (dt).*z
        
        t += dt

        for i = 10:n-10
            csep[i] = (w[i+1,1] - w[i,1])/(w[i,1] - w[i-1,1])
        end
        
        #csepmax = maximum(csep)
        #println(t, "   ", maximum(csepmax), "    ", argmax(csep))
        # if csepmax > 2.
        #     println("Separation identified", "    csep=$csepmax", "   i_s=$(argmax(csep))")
        #     i_s = argmax(csep) 
        # end
        
        j1, j2 = separationJ(lamb1, lamb2, dt, dx)
        
    end
    
      
    return w, i_s, j1*1000, j2*1000
end


function newton(f, x1, x2, maxiter, tol)

    f1 = f(x1)
    n = 2
    l = length(x1)
    fl = length(f1)
    
    while n < maxiter+1
        f2 = f(x2)
        df = zeros(fl,l)
        for i = 2:fl
            for j = 2:l
                df[i,j] = (f2[i] - f1[i])/(x2[j] - x1[j])
            end
        end
        df[1,:] = df[2,:]
        df[:,1] = df[:,2]
        
        x1 = x2
        f1 = f2
        x2 = x1 - df*f1
        res = sqrt(mean((x2 .-x1).^2))
        println(res)
        if res < tol
            break
        end
        n += 1
    end
    return x2
end
