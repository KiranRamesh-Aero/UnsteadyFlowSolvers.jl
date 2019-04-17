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


function mappingAerofoilToFVGrid(qu::Array{Float64,1}, theta::Array{Float64,1}, ncell::Int16)

    xfvm = range(0,stop=pi,length=ncell)

    n = length(xfvm)
    dx = (xfvm[end]-xfvm[1])/(n)
    quFine = smoothEdgeVelocity(qu, theta, n, 20)

    return quFine, xfvm
end

function inviscidInterface(del::Array{Float64,1}, E::Array{Float64,1}, q::Array{Float64,1}, qu0::Array{Float64,1}, dt::Float64, surf::TwoDSurfThick, xfvm::Array{Float64,1})

    U0 = q[1:end]
    #x =x[1:n]
    #U0 = U0[1:n]

    U00 = qu0[1:end]

    #What is this?
    U0[U0.< 0.0] .= 1e-8
    U00[U00.< 0.0] .= 1e-8
    #println("finding the length",length(U00))
    #println("finding the m ", m)

    UxInter = Spline1D(xfvm, U0)
    Ux = derivative(UxInter, xfvm)  #diff1(xfvm, U0)
    Ut = temporalDeriv(U0, U00, dt)


    w1 = zeros(length(del))
    w2 = zeros(length(del))
    w1[1:end] = del[1:end]
    w2[1:end] = del[1:end].*(E[1:end] .+ 1.0)
    w0 = hcat(w1,w2)

    return w0, U0, Ut, Ux
end
