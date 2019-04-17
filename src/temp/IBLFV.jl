function FVMIBL(w::Array{Float64,2}, U::Array{Float64,1}, Ut::Array{Float64,1}, Ux::Array{Float64,1}, x::Array{Float64,1})

    n = Int(length(w)/2)
    #dx = Float64((x[end]-x[1])/n)
    dx = Float64((x[end]-x[1])/n)
    #dx_half = (x[end]-x[end-1])/2
    #dx = [x[2:end]-x[1:end-1]; dx_half]

    # correlate the unknown values from the del and E values
    del , E, FF ,B, S, dfde = correlate(w)
    #del , E, FF ,B, S, dfde  = calc_shapes(n,w)


    fL, fR, UipL ,UipR, FFipL, FFipR, dfdeipL, dfdeipR, wipL, wipR = fluxReconstruction(w , U, FF, dfde, del, E)

    #lamb1L ,lamb2L = eigenlamb(UipL, dfdeipL, FFipL, wipL)
    #lamb1R ,lamb2R = eigenlamb(UipR, dfdeipR, FFipR, wipR)

    #lamb1 ,lamb2 = eigenlamb(U, dfde, FF, w)
    lamb1 ,lamb2 = calc_eigen(E, FF, dfde, U)
    dt = calc_Dt(lamb1 ,lamb2, 0.5, dx)

    # two step forward Euler methods for adding the source term to the right hand-side
    # of the transport equations.
    #z = RHSSource(U,B, del,Ut, Ux, FF, E, S)



    # step 1 : assuming this as a homogeneous equation and advanced half a step
    w1 = w.+ ((fL - fR).*((dt/2)./(dx)))

    del , E, FF ,B, S, dfde = correlate(w1)

    fL, fR, UipL ,UipR, FFipL, FFipR, dfdeipL, dfdeipR, wipL, wipR = fluxReconstruction(w1 , U, FF, dfde, del, E)

    z = RHSSource(U,B, del, Ut, Ux, FF, E, S)

    #dtL = calc_Dt(UipL, dfdeipL, FFipL, wipL, 0.8, dx)
    #dtR = calc_Dt(UipR, dfdeipR, FFipR, wipR, 0.8, dx)

    #dt = calc_Dt(lamb1 ,lamb2, 0.6, dx)

    j1, j2 = sepeartionJ(lamb1, lamb2, dt, dx)

    # step 2 : by considering source terms advanced a full step using 2nd order midpoint rule
    w2 = (w) + (fL - fR).* ((dt)./(dx)) .+ (dt).*z


    return w2, dt, lamb1, lamb2
end

function init(n)

    # initial conditions of IBVP
    E_init = 0.4142 * ones(n)
    E = E_init
    B = 131.9*E_init.^3 - 167.32*E_init.^2 + 76.642.*E_init .- 11.068
    del = sqrt.(B.*0.005);
    F = 4.8274.*E_init.^4 - 5.9816*E_init.^3 + 4.0274*E_init.^2 + 0.23247.*E_init .+ 0.15174;

    return del, E, F ,B
end


function eigenlamb(U::Array{Float64,1}, dfde::Array{Float64,1}, FF::Array{Float64,1}, w::Array{Float64,2})
    lamb1 =  0.5*U.*((dfde .- 1.0)
            + sqrt.(1.0 .+ 4.0*FF .- 2.0.* dfde .- (4.0*((w[:,2]./w[:,1]) .-1.0).*dfde) .+ dfde.^2))
    lamb2 = 0.5*U.*((dfde .- 1.0)
            - sqrt.(1.0 .+ 4.0*FF .- 2.0.* dfde .- (4.0*((w[:,2]./w[:,1]) .-1.0).*dfde) .+ dfde.^2))

    return lamb1, lamb2
end

function calc_Dt(lamb1::Array{Float64,1}, lamb2::Array{Float64,1}, cfl::Float64, dx::Array{Float64})

    # calculate time step values based on eigenvalues

    dti = cfl.*(dx./(abs.(lamb1+lamb2)))
    dt = minimum(abs.(dti))

    #@printf(" Max l1: %1.5f, Min l1: %1.5f, Max l2: %1.5f, Min l2: %1.5f \n", maximum(lamb1), minimum(lamb1),maximum(lamb2), minimum(lamb2));

    if dt< 0.00001
        println("dt modified for max dt :$dt")
        dt =0.0005
    end

    return dt
end

function calc_Dt(lamb1::Array{Float64,1}, lamb2::Array{Float64,1}, cfl::Float64, dx::Float64)

    # calculate time step values based on eigenvalues

    dti = cfl.*(dx./(abs.(lamb1+lamb2)))
    dt = minimum(abs.(dti))

    #@printf(" Max l1: %1.5f, Min l1: %1.5f, Max l2: %1.5f, Min l2: %1.5f \n", maximum(lamb1), minimum(lamb1),maximum(lamb2), minimum(lamb2));

    if dt< 0.00001
        println("dt modified for max dt :$dt")
        dt =0.0005
    end

    return dt
end




function fluxReconstruction(w::Array{Float64,2}, U::Array{Float64,1}, FF::Array{Float64,1}, dfde::Array{Float64,1}, del::Array{Float64,1} , E::Array{Float64,1})

    delF = U.*(w[:,2]- w[:,1])
    EF = U.* FF.* w[:,1]
    F = hcat(delF, EF)
    # first-order approximation of left and right side of the i+1/2
    # cell interface

    # 1) approximation of edge velocity
    UR = U
    UL = U

    UR = limiter(UR, U, [U[2:end]; U[1]])
    UL = limiter(UL, U, [U[end]; U[1:end-1]])

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

    #fR[end,:] = [0.0;0.0]


    return fL, fR, UipL ,UipR, FFipL, FFipR, dfdeipL, dfdeipR, wipL, wipR

end

function maxWaveSpeed(Uip::Array{Float64,1}, wip::Array{Float64,2}, dfdeip::Array{Float64,1}, FFip::Array{Float64,1})

    # calculate wave speed at the interfaces of the cell
    wsP = abs.(0.5*Uip).*((dfdeip .- 1.0)
        + sqrt.(1.0 .+ 4.0*FFip .- 2.0.* dfdeip .- (4.0*((wip[:,2]./wip[:,1]) .-1.0).*dfdeip) .+ dfdeip.^2))

    wsN = abs.(0.5*Uip).*((dfdeip .- 1.0)
            - sqrt.(1.0 .+ 4.0*FFip .- 2.0.* dfdeip .- (4.0*((wip[:,2]./wip[:,1]) .-1.0).*dfdeip) .+ dfdeip.^2))

    return max(wsP, wsN)
end

function RHSSource(U::Array{Float64,1} ,B::Array{Float64,1}, del::Array{Float64,1},Ut::Array{Float64,1}, Ux::Array{Float64,1}, FF::Array{Float64,1}, E::Array{Float64,1}, S::Array{Float64,1} )

    z1 = B./(2.0*del) .- del.* (Ut./U) .- (E.+ 1.0).*del.*Ux
    z2 = S./del .- 2.0*E.*del.* (Ut./U) .- 2.0*FF.*del.*Ux
    z = hcat(z1,z2)

    return z

end


function calc_shapes(ncell::Int64, sol::Array{Float64})

    E = zeros(ncell); B = zeros(ncell); F = zeros(ncell); S = zeros(ncell)
    dfde = zeros(ncell); del = zeros(ncell)

    for i = 1:ncell
        del[i] = sol[i,1]
        E[i] = sol[i,2]/del[i] - 1.
        F[i] = 4.8274*E[i]^4 - 5.9816*E[i]^3 + 4.0274*E[i]^2 + 0.23247*E[i] + 0.15174
        if E[i] < -0.0616
            B[i] = -225.86*E[i]^3 - 3016.6*E[i]^2 - 208.68*E[i] - 17.915
        elseif E[i] > -0.0395
            B[i] = 131.9*E[i]^3 - 167.32*E[i]^2 + 76.642*E[i] - 11.068
        else
            B[i] = 0.5*(-225.86*E[i]^3 - 3016.6*E[i]^2 - 208.68*E[i] - 17.915
                        + 131.9*E[i]^3 - 167.32*E[i]^2 + 76.642*E[i] - 11.068)
        end
        if E[i] < -0.0582
            S[i] = 451.55*E[i]^3 + 2010*E[i]^2 + 138.96*E[i] + 11.296
        elseif E[i] > -0.042
            S[i] = -96.739*E[i]^3 + 117.74*E[i]^2 - 46.432*E[i] + 6.8074
        else
            S[i] = 0.5*(451.55*E[i]^3 + 2010*E[i]^2 + 138.96*E[i] + 11.296
                        - 96.739*E[i]^3 + 117.74*E[i]^2 - 46.432*E[i] + 6.8074)
        end
        dfde[i] = 4*4.8274*E[i]^3 - 3*5.9816*E[i]^2 + 2*4.0274*E[i] + 0.23247
    end

    return del, E, F ,B, S, dfde
end


function sepeartionJ(lamb1::Array{Float64,1}, lamb2::Array{Float64,1}, dt::Float64, dx::Float64)

    N1 = length(lamb1)
    lambj1 = sum(lamb1)/N1
    J1Sep = (dt)./(dx) .* (lamb1[1:end-1] - lamb1[2:end])./ (lambj1*N1)

    for i=2:N1-1
    if lamb1[i].<0.0
    J1Sep[i] =  (dt)/(dx) .* (lamb1[i+1] - lamb1[i-1])./ (lambj1*N1)
    end
    end

    N2 = length(lamb2)
    lambj2 = sum(lamb2)/N2
    J2Sep = (dt)/(dx) .* (lamb2[1:end-1] - lamb2[2:end])./ (lambj2*N2)

    for i=2:N2-1
    if lamb2[i].<0.0
    J2Sep[i] =  (dt)/(dx) * (lamb2[i+1] - lamb2[i-1]) / (lambj2*N2)
    end
    end

    return J1Sep, J2Sep

end


function initDt(w::Array{Float64,2}, U::Array{Float64,1})

    del , E, FF ,B, S, dfde = correlate(w)
    lamb1 ,lamb2 = eigenlamb(U, dfde, FF, w)
    n = Int(length(w)/2)
    dx = Float64(pi/n)
    dt = calc_Dt(lamb1 ,lamb2, 0.3, dx)

    return dt

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
        if lamb2[i]>lamb1[i]
            temp  = lamb2[i]
            lamb2[i] = lamb1[i]
            lamb1[i] = temp
        end
    end

    return lamb1, lamb2
end

function limiter(wExtrapolated::Array{Float64,1}, wCell::Array{Float64,1}, wNeighbor::Array{Float64,1})

    w = wExtrapolated;
    wMax = max.(wCell, wNeighbor);
    wMin = min.(wCell, wNeighbor);
    w[w.>wMax] = wMax[w.>wMax];
    w[w.<wMin] = wMin[w.<wMin];

    return w

end
