function calc_shapes(ncell::Int64, sol::Array{Float64})

    E = zeros(ncell+2); B = zeros(ncell+2); F = zeros(ncell+2); S = zeros(ncell+2)
    dfde = zeros(ncell+2); del = zeros(ncell+2)

    for i = 1:ncell+2
        del[i] = sol[1,i]
        E[i] = sol[2,i]/del[i] - 1.
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

    return E, F, B, S, dfde, del
end

function calc_eigen(ncell::Int64, E::Array{Float64}, F::Array{Float64},
                    dfde::Array{Float64}, ue::Array{Float64})

    lamb1 = zeros(ncell+2); lamb2 = zeros(ncell+2)

    for i = 1:ncell+2
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

function calc_dt(ncell::Int64, cfl::Float64, dx::Float64, lamb1::Array{Float64},
                 lamb2::Array{Float64})

    dt = 10000
    for i = 2:ncell+1
        dti = cfl*(dx/(abs(lamb1[i]) + abs(lamb2[i])))
        if dti < dt
            dt = dti
        end
    end

    return dt
end

function calc_fluxes(ncell::Int64, E::Array{Float64}, F::Array{Float64},
                     del::Array{Float64}, ue::Array{Float64}, lamb1::Array{Float64}, lamb2::Array{Float64})

    Apos = zeros(2,2); Aneg = zeros(2,2)
    flux = zeros(2,2,ncell+2)

    for i = 1:ncell+2
        if lamb1[i] >= 0. && lamb2[i] >= 0.
            flux[1,1,i] = ue[i]*E[i]*del[i]
            flux[1,2,i] = ue[i]*F[i]*del[i]
        elseif lamb1[i] <= 0. && lamb2[i] <= 0.
            flux[1,:,i] = 0.
        else
            Apos[1,1] = ue[i]*lamb1[i]/(lamb1[i] - lamb2[i])*(-1. - lamb2[i]/ue[i])
            Apos[1,2] = ue[i]*lamb1[i]/(lamb1[i] - lamb2[i])
            Apos[2,1] = -ue[i]*lamb1[i]/(lamb1[i] - lamb2[i])*(1. + lamb1[i]/ue[i])*(1. + lamb2[i]/ue[i])
            Apos[2,2] = ue[i]*lamb1[i]/(lamb1[i] - lamb2[i])*(1. + lamb1[i]/ue[i])
            flux[1,1,i] = Apos[1,1]*del[i] + Apos[1,2]*(E[i] + 1.)*del[i]
            flux[1,2,i] = Apos[2,1]*del[i] + Apos[2,2]*(E[i] + 1.)*del[i]
        end
    end

    for i = 1:ncell+2
        if lamb1[i] >= 0. && lamb2[i] >= 0.
            flux[2,1,i] = 0
            flux[2,2,i] = 0
        elseif lamb1[i] < 0. && lamb2[i] < 0.
            flux[2,1,i]= ue[i]*E[i]*del[i]
            flux[2,2,i]= ue[i]*F[i]*del[i]
        else
            Aneg[1,1] = ue[i]*lamb2[i]/(lamb1[i]-lamb2[i])*(1. + lamb1[i]/ue[i])
            Aneg[1,2] = -ue[i]*lamb2[i]/(lamb1[i]-lamb2[i])
            Aneg[2,1] = ue[i]*lamb2[i]/(lamb1[i]-lamb2[i])*(1. + lamb1[i]/ue[i])*(1. + lamb2[i]/ue[i])
            Aneg[2,2] = ue[i]*lamb2[i]/(lamb1[i]-lamb2[i])*(-1. - lamb2[i]/ue[i])
            flux[2,1,i] = Aneg[1,1]*del[i] + Aneg[1,2]*(E[i]+1)*del[i]
            flux[2,2,i] = Aneg[2,1]*del[i] + Aneg[2,2]*(E[i]+1)*del[i]
        end
    end

    return flux
end

function matuCyl(tTot::Float64, ncell::Int64 = 200, cfl::Float64 = 0.8, nstage::Int64 = 2)

    dx = pi/ncell

    #Timestep determined dynamically
    dt = 0

    uinf = 1.; c = 1.; re = 1000.

    x = zeros(ncell+2); ue = zeros(ncell+2); uex = zeros(ncell+2); uet = zeros(ncell+2)
    E = zeros(ncell+2); B = zeros(ncell+2); F = zeros(ncell+2); S = zeros(ncell+2); dfde = zeros(ncell+2);
    del = zeros(ncell+2);

    sol = zeros(2,ncell+2); solt = zeros(2,ncell+2)
    lamb1 = zeros(ncell+2); lamb2 = zeros(ncell+2)
    Jsep = zeros(2,ncell+2)
    Csep = zeros(ncell+2)

    flux = zeros(2,2,ncell+2)
    cf = zeros(ncell+2)

    #Set sources (cell centred)
    for i = 2:ncell+1
        x[i] = (i-1.5)*dx
        ue[i] = 2*sin(x[i])
        uex[i] = 2*cos(x[i])
    end
    x[1] = 2*x[2] - x[3]
    x[ncell+2] = 2*x[ncell+1] - x[ncell]
    uet[:] = zeros(ncell+2)
    ue[1] = 0.
    ue[ncell+2] = 0.
    uex[1] = cos(0.)
    uex[ncell+2] = cos(pi)

    #Set initial conditions
    E_init = 0.4142
    for i = 1:ncell+2
        E[i] = E_init
        B[i] = 131.9*E_init^3 - 167.32*E_init^2 + 76.642*E_init - 11.068
        del[i] = sqrt.(B[i]*0.005)
        F[i] = 4.8274*E_init^4 - 5.9816*E_init^3 + 4.0274*E_init^2 + 0.23247*E_init + 0.15174
    end

    t = 0.005

    sol[1,:] = del[:]
    sol[2,:] = del[:].*(E[:] .+ 1)

    rhs = zeros(2, ncell+2)

    istep = 0
    while t < tTot
        istep += 1

        #Update ghost cells
        sol[:,1] = 2*sol[:,2] - sol[:,3]
        sol[:,ncell+2] = 2*sol[:,ncell+1] - sol[:,ncell]

        #Calculate derived quantities
        E, F, B, S, dfde, del = calc_shapes(ncell, sol)

        #Compute eigenvalues
        lamb1, lamb2 = calc_eigen(ncell, E, F, dfde, ue)

        #Compute timestep
        dt = calc_dt(ncell, cfl, dx, lamb1, lamb2)

        flux = calc_fluxes(ncell, E, F, del, ue, lamb1, lamb2)

        rhs[:,:] = zeros(2, ncell+2)

        #The ghost rhs are never used - dont matter
        for i = 2:ncell+1
            rhs[1,i] = B[i]/(2*del[i]) - del[i]*uet[i]/ue[i] - (E[i]+1)*del[i]*uex[i]
            rhs[2,i] = S[i]/del[i] - 2*E[i]*del[i]*uet[i]/ue[i] - 2*F[i]*del[i]*uex[i]
        end

        for i = 2:ncell
            solt[1,i] = sol[1,i] - (dt/dx)*(flux[1,1,i] - flux[1,1,i-1] + flux[2,1,i+1]
                                            -flux[2,1,i]) + dt*rhs[1,i]
            solt[2,i] = sol[2,i] - (dt/dx)*(flux[1,2,i] - flux[1,2,i-1] + flux[2,2,i+1]
                                            -flux[2,2,i]) + dt*rhs[2,i]
        end

        i = ncell+1
        solt[1,i] = sol[1,i] - (dt/dx)*(flux[1,1,i] - flux[1,1,i-1]) + dt*rhs[1,i]
        solt[2,i] = sol[2,i] - (dt/dx)*(flux[1,2,i] - flux[1,2,i-1]) + dt*rhs[2,i]

        #Update ghost cells
        solt[:,1] = 2*solt[:,2] - solt[:,3]
        solt[:,ncell+2] = 2*solt[:,ncell+1] - solt[:,ncell]

        #Calculate derived quantities
        E, F, B, S, dfde, del = calc_shapes(ncell, solt)

        #Compute eigenvalues
        lamb1, lamb2 = calc_eigen(ncell, E, F, dfde, ue)

        fluxt = calc_fluxes(ncell, E, F, del, ue, lamb1, lamb2)

        rhs[:,:] = zeros(2, ncell+2)
        #The ghost rhs are never used - dont matter
        for i = 2:ncell+1
            rhs[1,i] = B[i]/(2*del[i]) - del[i]*uet[i]/ue[i] - (E[i]+1)*del[i]*uex[i]
            rhs[2,i] = S[i]/del[i] - 2*E[i]*del[i]*uet[i]/ue[i] - 2*F[i]*del[i]*uex[i]
        end

        i = 2
        sol[1,i] = 0.5*(sol[1,i] + solt[1,i]) - (0.5*dt/dx)*(flux[1,1,i] - flux[1,1,i-1]
                                                             - flux[2,1,i] + 2*flux[2,1,i+1] - flux[2,1,i+2] + fluxt[1,1,i] - fluxt[1,1,i-1]
                                                             + fluxt[2,1,i+1] - fluxt[2,1,i]) + 0.5*dt*rhs[1,i]
        sol[2,i] = 0.5*(sol[2,i] + solt[2,i]) - (0.5*dt/dx)*(flux[1,2,i] - flux[1,2,i-1]
                                                             - flux[2,2,i] + 2*flux[2,2,i+1]-flux[2,2,i+2]+fluxt[1,2,i] - fluxt[1,2,i-1]
                                                             + fluxt[2,2,i+1] - fluxt[2,2,i])+0.5*dt*rhs[2,i]
        for i = 3:ncell
            sol[1,i] = 0.5*(sol[1,i] + solt[1,i]) - (0.5*dt/dx)*(flux[1,1,i] - 2*flux[1,1,i-1] + flux[1,1,i-2]
                                                                 - flux[2,1,i] + 2*flux[2,1,i+1] - flux[2,1,i+2] + fluxt[1,1,i] - fluxt[1,1,i-1]
                                                                 + fluxt[2,1,i+1] - fluxt[2,1,i]) + 0.5*dt*rhs[1,i]
            sol[2,i] = 0.5*(sol[2,i] + solt[2,i]) - (0.5*dt/dx)*(flux[1,2,i] - 2*flux[1,2,i-1] + flux[1,2,i-2]
                                                                 - flux[2,2,i] + 2*flux[2,2,i+1] - flux[2,2,i+2] + fluxt[1,2,i] - fluxt[1,2,i-1]
                                                                 + fluxt[2,2,i+1] - fluxt[2,2,i]) + 0.5*dt*rhs[2,i]
        end
        i = ncell+1
        sol[1,i] = 0.5*(sol[1,i] + solt[1,i]) - (0.5*dt/dx)*(flux[1,1,i] - 2*flux[1,1,i-1] + flux[1,1,i-2]
                                                             +fluxt[1,1,i] - fluxt[1,1,i-1]) + 0.5*dt*rhs[1,i]
        sol[2,i] = 0.5*(sol[2,i] + solt[2,i]) - (0.5*dt/dx)*(flux[1,2,i] - 2*flux[1,2,i-1] + flux[1,2,i-2]
                                                             +fluxt[1,2,i] - fluxt[1,2,i-1]) + 0.5*dt*rhs[2,i]

        t = t + dt
        if t > tTot
            break
        end

        for i = 2:ncell+1
            cf[i] = B[i]*ue[i]*c/(del[i]*uinf*re)
            Csep[i] = (del[i+1] - del[i])/(del[i] - del[i-1])
            if Csep[i] > 10.
                println("singularity (separation) detected at t=$t, x=$(x[i])")
            end
        end
        #println("$istep \t", "$dt \t", "$t \t", "$(maximum(Csep)) \t", "$(minimum(Csep)) \n")
        println("$t \t")
        display(plot(x/pi,[del, E] , layout=(2,1), legend = false))
        sleep(0.1)

    end

mat = zeros(0, 8)
mat = mat'
for i = 2:ncell+1
    mat = hcat(mat,[x[i], ue[i], uex[i], uet[i], del[i], E[i], cf[i], Csep[i]])
end
mat = mat'

return mat
end

# f = open("t$istep", "w")
# write(f, ["#x \t", "#ue \t", "#uex \t", "#uet \t", "#del \t", "#E \t", "#cf \t", "#Csep \n"])
# writedlm(f, mat)
# close(f)
