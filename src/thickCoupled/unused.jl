
function transpSimul(surf::TwoDSurfThickBL, curfield::TwoDFlowField, ncell::Int64, nsteps::Int64 = 300, dtstar::Float64 = 0.015, startflag = 0, writeflag = 0, writeInterval = 1000., delvort = delNone(); maxwrite = 50, nround=6)

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

    Re = 1000


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

    qucprev = zeros(surf.ndiv-1)


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

        #Solve for all unknowns together through a newton iteration
        #start with an inviscid solution and BL from prev time step
        xinv = surf.LHS[1:surf.ndiv*2-2, 1:surf.naterm*2+1] \ surf.RHS[1:surf.ndiv*2-2]

        #Starting solution for del and E
        #srange = collect(range(0, stop=surf.su[end], length=surf.nfvm))
        #delInter = Spline1D(surf.su, surf.delu)
        #EInter = Spline1D(surf.su, surf.Eu)
        #delstart = evaluate(delInter, srange)
        #Estart = evaluate(EInter, srange)

        xinit = zeros(2*surf.naterm+1)

        xinit[1:2*surf.naterm+1] = xinv[:]
        #xinit[2*surf.naterm+2:2*surf.naterm+surf.ndiv] = surf.delu[:]
        #xinit[2*surf.naterm+surf.ndiv+1:end] = surf.Eu[:]


        #cache1=DiffCache(surf.theta)

        # IBLsimul!(Fvec, xvec) = transResidual!(Fvec, xvec, surf.naterm, surf.uref, surf.theta, xloc_tev, zloc_tev, vcore, surf.bnd_x_u, surf.bnd_z_u, surf.bnd_x_l, surf.bnd_z_l, surf.uind_u, surf.wind_u, surf.uind_l, surf.wind_l, [curfield.u[1]; curfield.w[1]], surf.kinem.alpha, surf.kinem.alphadot, surf.kinem.u, surf.kinem.hdot, surf.cam, surf.thick, surf.cam_slope, surf.thick_slope, surf.pvt, surf.c, surf.su, surf.delu, surf.Eu, surf.LHS, surf.RHS, quprev, dt, t, Re)

        resfn(xvec) = transResidual(xvec, surf.naterm, surf.uref, surf.theta, xloc_tev, zloc_tev, vcore, surf.bnd_x_u, surf.bnd_z_u, surf.bnd_x_l, surf.bnd_z_l, surf.uind_u, surf.wind_u, surf.uind_l, surf.wind_l, [curfield.u[1]; curfield.w[1]], surf.kinem.alpha, surf.kinem.alphadot, surf.kinem.u, surf.kinem.hdot, surf.cam, surf.thick, surf.cam_slope, surf.thick_slope, surf.pvt, surf.c, surf.su, surf.delu, surf.Eu, surf.LHS, surf.RHS, qucprev, dt, t, Re)

        tol = 1e-6
        xsoln = newton(resfn, xinit, xinit .+ 0.01, 5, tol)
        
        #soln = nlsolve(resfn, xinit, iterations=5, show_trace=true)
        
        println("soln")
        error("here")

        #Iterate for viscous solution and interaction

        # if istep ==70 && iter == 1
        #     figure(1)
        #     plot(surf.x, iter_delu)
        #     figure(2)
        #     plot(srange, quxf)
        #     figure(3)
        #     plot(surf.x, wtu)
        #     error("here")
        # end
        
        #Assign inviscid soln
        #surf.aterm[]
        
        
        #Assign bl
        surf.delu[:] = iter_delu[:]
        #surf.delu[end] = surf.delu[end-1]
        #surf.Eu[end] = surf.Eu[end-1]

        surf.Eu[:] = iter_Eu[:]

        #Calculate adot
        update_atermdot(surf, dt)

        #Set previous values of aterm to be used for derivatives in next time step
        surf.a0prev[1] = surf.a0[1]
        for ia = 1:3
            surf.aprev[ia] = surf.aterm[ia]
        end

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



function dtfunjac(w, U, Ut, Ux, x)

    n = Int(length(w)/2)

    del = w[1:n]
    E = (w[n+1:end]./w[1:n]) .- 1.

    dx = Float64((x[end]-x[1])/n)

    # correlate the unknown values from the del and E values
    FF = 4.8274*E.^4 - 5.9816*E.^3 + 4.0274*E.^2 + 0.23247.*E .+ 0.15174

    B = 0.5*(-225.86.*E.^3 - 3016.6.*E.^2 - 208.68*E .- 17.915
             + 131.9*E.^3 - 167.32*E.^2 + 76.642.*E .- 11.068)
    B[E.<-0.0616] = -225.86.*E[E.<-0.0616].^3 - 3016.6.*E[E.<-0.0616].^2 - 208.68*E[E.<-0.0616] .- 17.915
    B[E.>-0.0395] =  131.9*E[E.>-0.0395].^3 - 167.32*E[E.>-0.0395].^2 + 76.642.*E[E.>-0.0395] .- 11.068

    S = 0.5*(451.55*E.^3 + 2010*E.^2 + 138.96*E .+ 11.296
             - 96.739*E.^3 + 117.74*E.^2 - 46.432*E .+ 6.8074)
    S[E.<-0.0582] = 451.55*E[E.<-0.0582].^3 + 2010*E[E.< -0.0582].^2 + 138.96*E[E.< -0.0582] .+ 11.296
    S[E.>-0.042]  = -96.739*E[E.>-0.042].^3 + 117.74*E[E.> -0.042].^2 - 46.432*E[E.> -0.042] .+ 6.8074

    dfde = 4*4.8274*E.^3 - 3*5.9816*E.^2 + 2*4.0274*E .+ 0.23247

    fL, fR, UipL ,UipR, FFipL, FFipR, dfdeipL, dfdeipR, wipL, wipR = dtfunfluxReconjac(del, E , U, FF, dfde)

    z = RHSSourcejac(U, B, del, Ut, Ux, FF, E, S)

    dfdt = ([fL[:,1]; fL[:,2]] - [fR[:,1]; fR[:,2]])./ones(2*n)*dx .+ [z[:,1]; z[:,2]]

    return dfdt
end


function dtfunFVM(w::Array{Float64,1}, U::Array{Float64,1}, Ut::Array{Float64,1}, Ux::Array{Float64,1}, x::Array{Float64,1})

    n = Int(length(w)/2)

    del = w[1:n]
    E = (w[n+1:end]./w[1:n]) .- 1.

    dx = Float64((x[end]-x[1])/n)

    # correlate the unknown values from the del and E values
    FF = 4.8274*E.^4 - 5.9816*E.^3 + 4.0274*E.^2 + 0.23247.*E .+ 0.15174

    B = 0.5*(-225.86.*E.^3 - 3016.6.*E.^2 - 208.68*E .- 17.915
             + 131.9*E.^3 - 167.32*E.^2 + 76.642.*E .- 11.068)
    B[E.<-0.0616] = -225.86.*E[E.<-0.0616].^3 - 3016.6.*E[E.<-0.0616].^2 - 208.68*E[E.<-0.0616] .- 17.915
    B[E.>-0.0395] =  131.9*E[E.>-0.0395].^3 - 167.32*E[E.>-0.0395].^2 + 76.642.*E[E.>-0.0395] .- 11.068

    S = 0.5*(451.55*E.^3 + 2010*E.^2 + 138.96*E .+ 11.296
             - 96.739*E.^3 + 117.74*E.^2 - 46.432*E .+ 6.8074)
    S[E.<-0.0582] = 451.55*E[E.<-0.0582].^3 + 2010*E[E.< -0.0582].^2 + 138.96*E[E.< -0.0582] .+ 11.296
    S[E.>-0.042]  = -96.739*E[E.>-0.042].^3 + 117.74*E[E.> -0.042].^2 - 46.432*E[E.> -0.042] .+ 6.8074

    dfde = 4*4.8274*E.^3 - 3*5.9816*E.^2 + 2*4.0274*E .+ 0.23247

    fL, fR, UipL ,UipR, FFipL, FFipR, dfdeipL, dfdeipR, wipL, wipR = dtfunfluxRecon(del, E , U, FF, dfde)

    z = RHSSource(U, B, del, Ut, Ux, FF, E, S)

    dfdt = ([fL[:,1]; fL[:,2]] - [fR[:,1]; fR[:,2]])./ones(2*n)*dx .+ [z[:,1]; z[:,2]]

    return dfdt
end


function RHSSourcejac(U ,B, del,Ut, Ux, FF, E, S)

    z1 = B./(2.0*del) .- del.* (Ut./U) .- (E.+ 1.0).*del.*Ux
    z2 = S./del .- 2.0*E.*del.* (Ut./U) .- 2.0*FF.*del.*Ux
    z = hcat(z1,z2)

    return z

end


function calc_eigenjac(E, F,
                    dfde, ue)

    ncell = length(E)
    lamb1 = zeros(Number, ncell); lamb2 = zeros(Number, ncell)

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

function maxWaveSpeedjac(Uip, wip, dfdeip, FFip)

    # calculate wave speed at the interfaces of the cell
    wsP = abs.(0.5*Uip).*((dfdeip .- 1.0)
        + sqrt.(1.0 .+ 4.0*FFip .- 2.0.* dfdeip .- (4.0*((wip[:,2]./wip[:,1]) .-1.0).*dfdeip) .+ dfdeip.^2))

    wsN = abs.(0.5*Uip).*((dfdeip .- 1.0)
            - sqrt.(1.0 .+ 4.0*FFip .- 2.0.* dfdeip .- (4.0*((wip[:,2]./wip[:,1]) .-1.0).*dfdeip) .+ dfdeip.^2))

    return max(wsP, wsN)
end

function limiterjac(wExtrapolated, wCell, wNeighbor)

    w = wExtrapolated;
    wMax = max.(wCell, wNeighbor);
    wMin = min.(wCell, wNeighbor);
    w[w.>wMax] = wMax[w.>wMax];
    w[w.<wMin] = wMin[w.<wMin];

    return w

end

function dtfunfluxReconjac(delF, Ef, U, FF, dfde)

    F = hcat(delF, Ef)
    # first-order approximation of left and right side of the i+1/2
    # cell interface

    # 1) approximation of edge velocity
    UR = U
    UL = U

    UR = limiterjac(UR, U, [U[2:end]; U[1]])
    UL = limiterjac(UL, U, [U[end]; U[1:end-1]])

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
    wR = [delF delF.*(Ef .+ 1.)]
    wL = [delF delF.*(Ef .+ 1.)]
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
    wsL = maxWaveSpeedjac(UipL, wipL, dfdeipL, FFipL)
    wsR = maxWaveSpeedjac(UipR, wipR, dfdeipR, FFipR)
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


function dtfunfluxRecon(delF::Array{Float64,1}, Ef::Array{Float64,1}, U::Array{Float64,1}, FF::Array{Float64,1}, dfde::Array{Float64,1})

    F = hcat(delF, Ef)
    # first-order approximation of left and right side of the i+1/2
    # cell interface

    # 1) approximation of edge velocity
    UR = U
    UL = U

    UR = limiter(UR, U, [U[2:end]; U[1]])
    UL = limiter(UL, U, [U[end]; U[1:end-1]])

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
    wR = [delF delF.*(Ef .+ 1.)]
    wL = [delF delF.*(Ef .+ 1.)]
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


function FVMIBLgrid(w, U, Ut, Ux, x, dt)

    n = length(x)

    dx = zeros(n)
    dx[2:end] = diff(x)
    dx[1] = dx[2]

    #csep = zeros(n)

    wfin = zeros(Number, size(w))

    wfin[:,:] = w[:,:]

    #while t < t_tot
        # correlate the unknown values from the del and E values
        del, E, FF ,B, S, dfde = correlate(w)

        fL, fR, UipL ,UipR, FFipL, FFipR, dfdeipL, dfdeipR, wipL, wipR = dtfunfluxReconjac(del, E, U, FF, dfde)

        #lamb1 ,lamb2 = calc_eigenjac(E, FF, dfde, U)
        #dt = calc_Dtjac(lamb1 ,lamb2, 0.5, ones(length(x)).*dx)

        #if t + dt > t_tot
     #       dt = t_tot - t
      #  end

        # step 1 : assuming this as a homogeneous equation and advanced half a step
        w1 = wfin.+ ((fL - fR).*((dt/2)./(dx)))

        del , E, FF ,B, S, dfde = correlate(w1)

        fL, fR, UipL ,UipR, FFipL, FFipR, dfdeipL, dfdeipR, wipL, wipR = dtfunfluxReconjac(del, E, U, FF, dfde)

        z = RHSSourcejac(U, B, del, Ut, Ux, FF, E, S)

        # step 2 : by considering source terms advanced a full step using 2nd order midpoint rule
        w2 = (wfin) + (fL - fR).* ((dt)./(dx)) .+ (dt).*z

        wfin[:,:] = w2[:,:]
       # t += dt

        # for i = 20:n-20
        #     csep[i] = (wfin[i+1,1] - w[i,1])/(x[i+1,1] - w[i,1])/pi
        # end

  #  end

    return x, wfin
end

