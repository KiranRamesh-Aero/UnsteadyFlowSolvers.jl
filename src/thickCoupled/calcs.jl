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
        + sqrt.(abs.(1.0 .+ 4.0*FFip .- 2.0.* dfdeip .- (4.0*((wip[:,2]./wip[:,1]) .-1.0).*dfdeip)) .+ dfdeip.^2))

    wsN = abs.(0.5*Uip).*((dfdeip .- 1.0)
            - sqrt.(abs.(1.0 .+ 4.0*FFip .- 2.0.* dfdeip .- (4.0*((wip[:,2]./wip[:,1]) .-1.0).*dfdeip) .+ dfdeip.^2)))

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

        #csepmax = maximum(csep)
        #println(t, "   ", maximum(csepmax), "    ", argmax(csep))
        # if csepmax > 2.
        #     println("Separation identified", "    csep=$csepmax", "   i_s=$(argmax(csep))")
        #     i_s = argmax(csep) 
        # end
        
        j1, j2 = separationJ(lamb1, lamb2, dt, dx)
        jmax = maximum(j1)
        if jmax > 1e-4
            println("Separation identified", "    jsep=$jmax", "   i_s=$(argmax(j1))")
            i_s = argmax(j1) 
        else
            i_s = 0
        end
        
    end
    
    return w, i_s
end


function find_nacaCoef(surf::TwoDSurfThick, thick::Array{Float64}, bstart)

    th = parse(Int, surf.coord_file[7:8])/100.
    
    b1 = 0.2969
    
    @. nacath(x, b) = 5*th*(b1*sqrt(x) + b[1]*x + b[2]*x^2 + b[3]*x^3 + b[4]*x^4)
    
    fit = curve_fit(nacath, surf.x[2:end], thick[2:end], bstart)
    
    b = coef(fit)

    return b
end


    
