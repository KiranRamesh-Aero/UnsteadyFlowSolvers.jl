function calcEigenValues(soln::Solutions, fluxsplit::FluxSplittingParameters, ue::Array{Float64},ncell::Int64)

    for i = 1:ncell+2
        a_q = 1.
        b_q = -ue[i]*(fluxsplit.dfde[i] - 1)
        c_q = ue[i] * ue[i] * (fluxsplit.E[i]*fluxsplit.dfde[i] - fluxsplit.F[i])
        soln.lamb1[i] = (-b_q + sqrt(b_q*b_q - 4*a_q*c_q))/(2*a_q)
        soln.lamb2[i] = (-b_q  -sqrt(b_q*b_q - 4*a_q*c_q))/(2*a_q)

        #Always have lamb1 > lamb2
        if soln.lamb2[i]>soln.lamb1[i]
            temp  = sol.lamb2[i]
            soln.lamb2[i] = soln.lamb1[i]
            soln.lamb1[i] = temp
        end
    end
end


function calcDt(opCond::OperationalConditions, sol::Solutions, ncell::Int64)

    dt = 10000.0
    for i = 1:ncell+2
        dti = opCond.cfl*(opCond.dx/(abs(sol.lamb1[i]) + abs(sol.lamb2[i])))
        if dti < dt
            dt = dti
        end
    end

    if(dt==10000.0)
        error("Problem diverges")
    end
    return dt
end

function calcDtCirc(cfl::Float64, dx::Float64, sol::Solutions)

    dt = 10000.0
    for i = 1:length(sol.lamb1)
        dti = cfl*(dx/(abs(sol.lamb1[i]) + abs(sol.lamb2[i])))
        if (dti < dt)
            dt = dti
        end
    end
    return dt
end


function calcFluxes(fluxsplit::FluxSplittingParameters, flux::Array{Float64,3} ,ue::Array{Float64,1}, sol::Solutions, ncell::Int64)

    Apos = zeros(2,2); Aneg = zeros(2,2)

    for i = 1:ncell+2
        if sol.lamb1[i] >= 0. && sol.lamb2[i] >= 0.
            flux[1,1,i] = ue[i]*fluxsplit.E[i]*fluxsplit.del[i]
            flux[1,2,i] = ue[i]*fluxsplit.F[i]*fluxsplit.del[i]
        elseif sol.lamb1[i] <= 0. && sol.lamb2[i] <= 0.
            flux[1,:,i] = zeros(length(flux[1,:,i]))

        else
            Apos[1,1] = ue[i]*sol.lamb1[i]/(sol.lamb1[i] - sol.lamb2[i])*(-1. - sol.lamb2[i]/ue[i])
            Apos[1,2] = ue[i]*sol.lamb1[i]/(sol.lamb1[i] - sol.lamb2[i])
            Apos[2,1] = -ue[i]*sol.lamb1[i]/(sol.lamb1[i] - sol.lamb2[i])*(1. + sol.lamb1[i]/ue[i])*(1. + sol.lamb2[i]/ue[i])
            Apos[2,2] = ue[i]*sol.lamb1[i]/(sol.lamb1[i] - sol.lamb2[i])*(1. + sol.lamb1[i]/ue[i])
            flux[1,1,i] = Apos[1,1]*fluxsplit.del[i] + Apos[1,2]*(fluxsplit.E[i] + 1.)*fluxsplit.del[i]
            flux[1,2,i] = Apos[2,1]*fluxsplit.del[i] + Apos[2,2]*(fluxsplit.E[i] + 1.)*fluxsplit.del[i]
        end
    end

    for i = 1:ncell+2
        if sol.lamb1[i] >= 0. && sol.lamb2[i] >= 0.
              flux[2,:,i] = zeros(length(flux[2,:,i]))
        elseif sol.lamb1[i] < 0. && sol.lamb2[i] < 0.
            flux[2,1,i]= ue[i]*fluxsplit.E[i]*fluxsplit.del[i]
            flux[2,2,i]= ue[i]*fluxsplit.F[i]*fluxsplit.del[i]
        else
            Aneg[1,1] = ue[i]*sol.lamb2[i]/(sol.lamb1[i]-sol.lamb2[i])*(1. + sol.lamb1[i]/ue[i])
            Aneg[1,2] = -ue[i]*sol.lamb2[i]/(sol.lamb1[i]-sol.lamb2[i])
            Aneg[2,1] = ue[i]*sol.lamb2[i]/(sol.lamb1[i]-sol.lamb2[i])*(1. + sol.lamb1[i]/ue[i])*(1. + sol.lamb2[i]/ue[i])
            Aneg[2,2] = ue[i]*sol.lamb2[i]/(sol.lamb1[i]-sol.lamb2[i])*(-1. - sol.lamb2[i]/ue[i])
            flux[2,1,i] = Aneg[1,1]*fluxsplit.del[i] + Aneg[1,2]*(fluxsplit.E[i]+1)*fluxsplit.del[i]
            flux[2,2,i] = Aneg[2,1]*fluxsplit.del[i] + Aneg[2,2]*(fluxsplit.E[i]+1)*fluxsplit.del[i]
        end
    end

    return flux
end

# Setting tha parameters for the right-hand-side of the equations

function rhs(sol::Solutions,fluxSp::FluxSplittingParameters,uet::Array{Float64,1},uex::Array{Float64,1},ue::Array{Float64,1},ncell::Int64)

    rhs =zeros(2,ncell+2)
    for i = 2:ncell+1
        rhs[1,i] = fluxSp.B[i]/(2*fluxSp.del[i]) - (fluxSp.del[i]*uet[i])/ue[i] - (fluxSp.E[i]+1.0)*fluxSp.del[i]*uex[i]
        rhs[2,i] = fluxSp.S[i]/fluxSp.del[i] - (2.0*fluxSp.E[i]*fluxSp.del[i]*uet[i])/ue[i] - (2.0*fluxSp.F[i]*fluxSp.del[i]*uex[i])
    end
# The boundary values do not need
    boundaryCorrection(rhs[1,:])
    boundaryCorrection(rhs[2,:])

    return rhs
end

# First-order (intermediate) flux spliiting function

function fluxSplittingSchema(flux::Array{Float64,3},sol::Solutions,rhs::Array{Float64,2},dt::Float64,dx::Float64,ncell::Int64)

    for i = 2:ncell
        sol.solt[1,i] = sol.sol[1,i] - (dt/dx)*(flux[1,1,i] - flux[1,1,i-1] + flux[2,1,i+1]
        -flux[2,1,i]) + dt*rhs[1,i]
        sol.solt[2,i] = sol.sol[2,i] - (dt/dx)*(flux[1,2,i] - flux[1,2,i-1] + flux[2,2,i+1]
        -flux[2,2,i]) + dt*rhs[2,i]
     end
    i = ncell+1
        sol.solt[1,i] = sol.sol[1,i] - (dt/dx)*(flux[1,1,i] - flux[1,1,i-1]) + dt*rhs[1,i]
        sol.solt[2,i] = sol.sol[2,i] - (dt/dx)*(flux[1,2,i] - flux[1,2,i-1]) + dt*rhs[2,i]

end

# Second-order (final) flux spliiting function

function fluxSplittingSchema(flux::Array{Float64,3},fluxt::Array{Float64,3},sol::Solutions,rhs::Array{Float64,2},dt::Float64,dx::Float64, ncell::Int64)

    i = 2
    sol.sol[1,i] = 0.5*(sol.sol[1,i] + sol.solt[1,i]) - (0.5*dt/dx)*(flux[1,1,i] - flux[1,1,i-1]
    - flux[2,1,i] + 2*flux[2,1,i+1] - flux[2,1,i+2] + fluxt[1,1,i] - fluxt[1,1,i-1]
    + fluxt[2,1,i+1] - fluxt[2,1,i]) + 0.5*dt*rhs[1,i]
    sol.sol[2,i] = 0.5*(sol.sol[2,i] + sol.solt[2,i]) - (0.5*dt/dx)*(flux[1,2,i] - flux[1,2,i-1]
    - flux[2,2,i] + 2*flux[2,2,i+1]-flux[2,2,i+2]+fluxt[1,2,i] - fluxt[1,2,i-1]
    + fluxt[2,2,i+1] - fluxt[2,2,i])+0.5*dt*rhs[2,i]

    for i = 3:ncell
        sol.sol[1,i] = 0.5*(sol.sol[1,i] + sol.solt[1,i]) - (0.5*dt/dx)*(flux[1,1,i] - 2*flux[1,1,i-1] + flux[1,1,i-2]
        - flux[2,1,i] + 2*flux[2,1,i+1] - flux[2,1,i+2] + fluxt[1,1,i] - fluxt[1,1,i-1]
        + fluxt[2,1,i+1] - fluxt[2,1,i]) + 0.5*dt*rhs[1,i]
        sol.sol[2,i] = 0.5*(sol.sol[2,i] + sol.solt[2,i]) - (0.5*dt/dx)*(flux[1,2,i] - 2*flux[1,2,i-1] + flux[1,2,i-2]
        - flux[2,2,i] + 2*flux[2,2,i+1] - flux[2,2,i+2] + fluxt[1,2,i] - fluxt[1,2,i-1]
        + fluxt[2,2,i+1] - fluxt[2,2,i]) + 0.5*dt*rhs[2,i]
    end
    i = ncell+1
    sol.sol[1,i] = 0.5*(sol.sol[1,i] + sol.solt[1,i]) - (0.5*dt/dx)*(flux[1,1,i] - 2*flux[1,1,i-1] + flux[1,1,i-2]
    +fluxt[1,1,i] - fluxt[1,1,i-1]) + 0.5*dt*rhs[1,i]
    sol.sol[2,i] = 0.5*(sol.sol[2,i] + sol.solt[2,i]) - (0.5*dt/dx)*(flux[1,2,i] - 2*flux[1,2,i-1] + flux[1,2,i-2]
    +fluxt[1,2,i] - fluxt[1,2,i-1]) + 0.5*dt*rhs[2,i]

end
