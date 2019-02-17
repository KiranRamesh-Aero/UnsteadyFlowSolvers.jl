function FVMIBL(w::Array{Float64,2}, U::Array{Float64,1}, Ut::Array{Float64,1}, Ux::Array{Float64,1})

    n = Int(length(w)/2)
    dx = Float64(1.0/n)

    # correlate the unknown values from the del and E values
    del , E, FF ,B, S, dfde = correlate(w)


    fL, fR, UipL ,UipR, FFipL, FFipR, dfdeipL, dfdeipR, wipL, wipR = fluxReconstruction(w , U, FF, dfde, del, E)

    dtL = calc_Dt(UipL, dfdeipL, FFipL, wipL, 0.5, dx)
    dtR = calc_Dt(UipR, dfdeipR, FFipR, wipR, 0.5, dx)

    dt = max(dtL,dtR)

    # two step forward Euler methods for adding the source term to the right hand side
    # of the transport equations.

    z = RHSSource(U,B, del,Ut, Ux, FF, E, S)

    # step 1 : assuming this as a homogeneous equation
    w1 = w .+ ((fL - fR) .* ((dt/2)/(dx)))

    del , E, FF ,B, S, dfde = correlate(w1)

    #fL, fR, UipL ,UipR, FFipL, FFipR, dfdeipL, dfdeipR, wipL, wipR = fluxReconstruction(w1 , U, FF, dfde, del, E)

    z = RHSSource(U,B, del,Ut, Ux, FF, E, S)

    #dtL = calc_Dt(UipL, dfdeipL, FFipL, wipL, 0.8, dx)
    #dtR = calc_Dt(UipR, dfdeipR, FFipR, wipR, 0.8, dx)

    #dt = max(dtL,dtR)

    w2 = w1 .+ (dt).*z
    #w2 = w + z.*dt
    return w2, dt
end

function init(n)

    E_init = 0.4142 * ones(n,1)
    E = E_init
    B = 131.9*E_init.^3 - 167.32*E_init.^2 + 76.642.*E_init .- 11.068
    del = sqrt.(B.*0.005);
    F= 4.8274.*E_init.^4 - 5.9816*E_init.^3 + 4.0274*E_init.^2 + 0.23247.*E_init .+ 0.15174;

    return del, E, F ,B
end


function calc_Dt(U::Array{Float64,1}, dfde::Array{Float64,1}, FF::Array{Float64,1}, w::Array{Float64,2}, cfl::Float64, dx::Float64 )

    lamb1 =  U.*((dfde .- 1.0)
            + sqrt.(1.0 .+ 4.0*FF .- 2.0.* dfde .- (4.0*((w[:,2]./w[:,1]) .-1.0).*dfde) .+ dfde.^2))
    lamb2 = U.*((dfde .- 1.0)
            - sqrt.(1.0 .+ 4.0*FF .- 2.0.* dfde .- (4.0*((w[:,2]./w[:,1]) .-1.0).*dfde) .+ dfde.^2))

    dti = cfl.*(dx./(max.(lamb1,lamb2)))
    dt = minimum(dti)
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
    ww = hcat(ws,ws)

    # flux reconstruction of left and right side of the i+1/2 interface
    fR = 0.5*((fpL + fpR) + ww.* (wipL - wipR))
    fL = [fR[end:end,:];fR[1:end-1,:]]

    # specifying the boundary conditions using internal extrapolation from the calculated flux of the
    # neighbouring cell centers
    fL[1,:] = [F[1,1]; F[1,2]]
    #fL[1,:] = [0; 0]
    fR[end,:] = [F[end,1]; F[end,2]]


    return fL, fR, UipL ,UipR, FFipL, FFipR, dfdeipL, dfdeipR, wipL, wipR

end

function maxWaveSpeed(Uip::Array{Float64,1}, wip::Array{Float64,2}, dfdeip::Array{Float64,1}, FFip::Array{Float64,1})

    ws = abs.(Uip).*((dfdeip .- 1.0)
        + sqrt.(1.0 .+ 4.0*FFip .- 2.0.* dfdeip .- (4.0*((wip[:,2]./wip[:,1]) .-1.0).*dfdeip) .+ dfdeip.^2))

    return ws
end

function RHSSource(U::Array{Float64,1} ,B::Array{Float64,1}, del::Array{Float64,1},Ut::Array{Float64,1}, Ux::Array{Float64,1}, FF::Array{Float64,1}, E::Array{Float64,1}, S::Array{Float64,1} )

    z1 = B./(2.0*del) .- del.* (Ut./U) .- (E.+ 1.0).*del.*Ux
    z2 = S./del .- 2.0*E.*del.* (Ut./U) .- 2.0*FF.*del.*Ux
    z = hcat(z1,z2)

    return z

end
