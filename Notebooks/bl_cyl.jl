

function FfromE(E::Float64)
      F = 4.8274*E^4 - 5.9816*E^3 + 4.0274*E^2 + 0.23247*E + 0.15174
end

function BfromE(E::Float64)
    if E < -0.0616
        B = -225.86*E^3 - 3016.6*E^2 - 208.68*E - 17.915
    elseif E > -0.0395
        B = 131.9*E^3 - 167.32*E^2 + 76.642*E - 11.068
    else
        B = 0.5*(-225.86*E^3 - 3016.6*E^2 - 208.68*E - 17.915 + 131.9*E^3 - 167.32*E^2 + 76.642*E - 11.068)
    end
end

function SfromE(E::Float64)
    if E < -0.0582
        S = 451.55*E^3 + 2010.*E^2 + 138.96*E + 11.296
    elseif E > -0.042
        S = -96.739*E^3 + 117.74*E^2 - 46.432*E + 6.8074
    else
        S = 0.5*(451.55*E^3 + 2010.*E^2 + 138.96*E + 11.296 - 96.739*E^3 + 117.74*E^2 - 46.432*E + 6.8074)
    end
end

function dfdefromE(E::Float64)
    dfde = 4*4.8274*E^3 - 3*5.9816*E^2 + 2*4.0274*E + 0.23247
end

function calcdt(dx::Float64, cfl::Float64, lamb::Array{Float64,2})
    dt = 10000
    for ic = 2:size(lamb,2)-1
        dti = cfl*(dx/(abs(lamb[1,ic]) + abs(lamb[2,ic])))
        if dti < dt
            dt=dti
        end
    end
    return dt
end

function calcEigen(ue::Vector{Float64}, E::Vector{Float64}, F::Vector{Float64}, dfde::Vector{Float64})
    ncell = length(E) - 2
    lamb = zeros(2,ncell+2)
    for ic = 1:ncell+2
        aq = 1.
        bq = -ue[ic]*(dfde[ic] - 1.)
        cq = ue[ic]*ue[ic]*(E[ic]*dfde[ic] - F[ic])
        lamb[1,ic] = (-bq + sqrt(bq*bq - 4*aq*cq))/(2*aq)
        lamb[2,ic] = (-bq - sqrt(bq*bq - 4*aq*cq))/(2*aq)

        #Always have lamb1 > lamb2
        if lamb[2,ic] > lamb[1,ic]
            temp  = lamb[2,ic]
            lamb[2,ic] = lamb[1,ic]
            lamb[1,ic] = temp
        end
    end
    return lamb
end


function calc_flux(lamb::Array{Float64,2}, ue::Vector{Float64}, E::Vector{Float64}, del::Vector{Float64}, F::Vector{Float64})
    ncell = length(ue) - 2
    flux = zeros(2,2,ncell+2)
    Apos = zeros(2,2)
    Aneg = zeros(2,2)

    for ic = 1:ncell+2
        if lamb[1,ic] >= 0. && lamb[2,ic] >= 0.
            flux[1,1,ic] = ue[ic]*E[ic]*del[ic]
            flux[1,2,ic] = ue[ic]*F[ic]*del[ic]
        elseif lamb[1,ic] < 0. && lamb[2,ic] < 0.
            flux[1,:,ic] = 0.
            else
            Apos[1,1] = (ue[ic]*lamb[1,ic]/(lamb[1,ic] - lamb[2,ic]))*(-1. - lamb[2,ic]/ue[ic])
            Apos[1,2] = ue[ic]*lamb[1,ic]/(lamb[1,ic] - lamb[2,ic])
            Apos[2,1] = -(ue[ic]*lamb[1,ic]/(lamb[1,ic] - lamb[2,ic]))*(1 + lamb[1,ic]/ue[ic])*(1 + lamb[2,ic]/ue[ic])
            Apos[2,2] = (ue[ic]*lamb[1,ic]/(lamb[1,ic] - lamb[2,ic]))*(1 + lamb[1,ic]/ue[ic])
            flux[1,1,ic] =  Apos[1,1]*del[ic] + Apos[1,2]*(E[ic] + 1.)*del[ic]
            flux[1,2,ic] = Apos[2,1]*del[ic] + Apos[2,2]*(E[ic] + 1.)*del[ic]
        end
    end
    for ic = 1:ncell+2
        if lamb[1,ic] >= 0. && lamb[2,ic] >= 0.
            flux[2,:,ic] = 0.
         elseif lamb[1,ic] < 0. && lamb[2,ic] < 0.
            flux[2,1,ic] = ue[ic]*E[ic]*del[ic]
            flux[2,2,ic] = ue[ic]*F[ic]*del[ic]
         else
            Aneg[1,1] = (ue[ic]*lamb[2,ic]/(lamb[1,ic] - lamb[2,ic]))*(1. + lamb[1,ic]/ue[ic])
            Aneg[1,2] = -ue[ic]*lamb[2,ic]/(lamb[1,ic] - lamb[2,ic])
            Aneg[2,1] = (ue[ic]*lamb[2,ic]/(lamb[1,ic] - lamb[2,ic]))*(1. + lamb[1,ic]/ue[ic])*(1. + lamb[2,ic]/ue[ic])
            Aneg[2,2] = (ue[ic]*lamb[2,ic]/(lamb[1,ic] - lamb[2,ic]))*(-1. - lamb[2,ic]/ue[ic])
            flux[2,1,ic] = Aneg[1,1]*del[ic] + Aneg[1,2]*(E[ic] + 1.)*del[ic]
            flux[2,2,ic] = Aneg[2,1]*del[ic] + Aneg[2,2]*(E[ic] + 1.)*del[ic]
         end
    end
    return flux
end


function IBLm_cyl(ntime = 1000, cfl = 0.5, ncell = 64, ttime = 1.4)

    sepflag = Int(0)
    dx = pi/(real(ncell) + 1.)
    t = 0.005

    x = zeros(ncell+2)
    ue = zeros(ncell+2)
    uex = zeros(ncell+2)
    uet = zeros(ncell+2)

    #Set sources
    for i = 1:ncell+2
        x[i] = real(i-1)*dx
        ue[i] = 2*sin(x[i])
        uex[i] = 2*cos(x[i])
    end

    uet[1:ncell+2] = 0

    E = zeros(ncell+2)
    B = zeros(ncell+2)
    del = zeros(ncell+2)
    F = zeros(ncell+2)
    dfde = zeros(ncell+2)
    S = zeros(ncell+2)
    unk = zeros(2,ncell+2)

    lamb = zeros(2,ncell+2)
    unkh = zeros(2,ncell+2)
    rhs = zeros(2,ncell+2)
    flux = zeros(2,2,ncell+2)
    crit = zeros(ncell+2)
    Apos = zeros(2,2)
    Aneg = zeros(2,2)

    #Set initial conditions
    E[:] = 0.4142
    for ic = 1:ncell+2
        B[ic] = BfromE(E[ic])
        F[ic] = FfromE(E[ic])
    end
    del = sqrt(B*t)

    unk[1,:] = del
    unk[2,:] = del.*(E + 1.)

    #Main loop over time steps
    for i = 1:ntime
        #i = 1
        unk[:,1] = 2*unk[:,2] - unk[:,3]
        unk[:,ncell+2] = 2*unk[:,ncell+1] - unk[:,ncell]
        #Calculate derived quantities
        for ic = 1:ncell+2
            del[ic] = unk[1,ic]
            E[ic] = (unk[2,ic]./del[ic]) - 1.
            F[ic] = FfromE(E[ic])
            B[ic] = BfromE(E[ic])
            S[ic] = SfromE(E[ic])
            dfde[ic] = dfdefromE(E[ic])
        end

        #Compute eigenvalues
        lamb = calcEigen(ue, E, F, dfde)

        #Compute timestep
        dt = calcdt(dx, cfl, lamb)

        #Compute fluxes
        flux = calc_flux(lamb, ue, E, del, F)

        #compute rhs
        for ic = 2:ncell+1
        rhs[1,ic] = B[ic]/(2*del[ic]) - del[ic]*uet[ic]/ue[ic] - (E[ic] + 1.)*del[ic]*uex[ic]
            rhs[2,ic] = S[ic]/del[ic] - 2*E[ic]*del[ic]*uet[ic]/ue[ic] - 2*F[ic]*del[ic]*uex[ic]
        end

        for ic = 2:ncell+1
            unkh[1,ic] = unk[1,ic] - (dt/dx)*(flux[1,1,ic] - flux[1,1,ic-1] + flux[2,1,ic+1]
                                              - flux[2,1,ic]) + dt*rhs[1,ic]
            unkh[2,ic] = unk[2,ic] - (dt/dx)*(flux[1,2,ic] - flux[1,2,ic-1] + flux[2,2,ic+1]
                                              - flux[2,2,ic]) + dt*rhs[2,ic]
        end
        ic = ncell+1
        #unkh[1,ic] = unk[1,ic] - (dt/dx)*(flux[1,1,ic] - flux[1,1,ic-1]) + dt*rhs[1,ic]
        #unkh[2,ic] = unk[2,ic] - (dt/dx)*(flux[1,2,ic] - flux[1,2,ic-1]) + dt*rhs[2,ic]

        #Update ghost cells
        unkh[:,1] = 2*unkh[:,2] - unkh[:,3]
        unkh[:,ncell+2] = 2*unkh[:,ncell+1] - unkh[:,ncell]

        #Calculate derived quantities
        for ic = 1:ncell+2
            del[ic] = unkh[1,ic]
            E[ic] = (unkh[2,ic]./del[ic]) - 1.
            F[ic] = FfromE(E[ic])
            B[ic] = BfromE(E[ic])
            S[ic] = SfromE(E[ic])
            dfde[ic] = dfdefromE(E[ic])
        end

        #Compute eigenvalues
        lamb = calcEigen(ue, E, F, dfde)

        #Compute fluxes
        fluxhalf = calc_flux(lamb, ue, E, del, F)

        #compute rhs
        for ic = 2:ncell+1
            rhs[1,ic] = B[ic]/(2*del[ic]) - del[ic]*uet[ic]/ue[ic] - (E[ic] + 1.)*del[ic]*uex[ic]
            rhs[2,ic] = S[ic]/del[ic] - 2*E[ic]*del[ic]*uet[ic]/ue[ic] - 2*F[ic]*del[ic]*uex[ic]
        end

        ic = 2

        unk[1,ic] = 0.5*(unk[1,ic] + unkh[1,ic]) -
        (0.5*dt/dx)*(flux[1,1,ic] - flux[1,1,ic-1] - flux[2,1,ic] +
        2*flux[2,1,ic+1] - flux[2,1,ic+2] + fluxhalf[1,1,ic] -
        fluxhalf[1,1,ic-1] + fluxhalf[2,1,ic+1] - fluxhalf[2,1,ic]) +
        0.5*dt*rhs[1,ic]

        unk[2,ic] = 0.5*(unk[2,ic] + unkh[2,ic]) -
        (0.5*dt/dx)*(flux[1,2,ic] - flux[1,2,ic-1] - flux[2,2,ic] +
        2*flux[2,2,ic+1] - flux[2,2,ic+2] + fluxhalf[1,2,ic] -
        fluxhalf[1,2,ic-1] + fluxhalf[2,2,ic+1] - fluxhalf[2,2,ic]) +
        0.5*dt*rhs[2,ic]

        for ic = 3:ncell

            unk[1,ic] = 0.5*(unk[1,ic] + unkh[1,ic]) -
            (0.5*dt/dx)*(flux[1,1,ic] - 2*flux[1,1,ic-1] +
            flux[1,1,ic-2] - flux[2,1,ic] + 2*flux[2,1,ic+1] -
            flux[2,1,ic+2] + fluxhalf[1,1,ic] -fluxhalf[1,1,ic-1] +
            fluxhalf[2,1,ic+1] - fluxhalf[2,1,ic]) + 0.5*dt*rhs[1,ic]

            unk[2,ic] = 0.5*(unk[2,ic] + unkh[2,ic]) -
            (0.5*dt/dx)*(flux[1,2,ic] - 2*flux[1,2,ic-1] +
            flux[1,2,ic-2] - flux[2,2,ic] + 2*flux[2,2,ic+1] -
            flux[2,2,ic+2] + fluxhalf[1,2,ic] - fluxhalf[1,2,ic-1] +
            fluxhalf[2,2,ic+1] - fluxhalf[2,2,ic]) + 0.5*dt*rhs[2,ic]

        end
        ic = ncell+1

        unk[1,ic] = 0.5*(unk[1,ic] + unkh[1,ic]) -
        (0.5*dt/dx)*(flux[1,1,ic] - 2*flux[1,1,ic-1] + flux[1,1,ic-2]
        +fluxhalf[1,1,ic] - fluxhalf[1,1,ic-1]) + 0.5*dt*rhs[1,ic]

        unk[2,ic] = 0.5*(unk[2,ic] + unkh[2,ic]) -
        (0.5*dt/dx)*(flux[1,2,ic] - 2*flux[1,2,ic-1] + flux[1,2,ic-2]
        +fluxhalf[1,2,ic] - fluxhalf[1,2,ic-1]) + 0.5*dt*rhs[2,ic]

        t = t + dt
        if t > ttime
            break
        end

        for ic = 2:ncell+1
            crit[ic] = abs((del[ic+1] - del[ic])/(del[ic] - del[ic-1]))
            if abs(crit[ic]) > 10.
                sepflag = 1
                println(t, " ",x[ic]," ", ue[ic]," ", crit[ic])
                break
            end
        end
        if t > ttime || sepflag == 1
            break
        end
    end
t, x, del, E
end





