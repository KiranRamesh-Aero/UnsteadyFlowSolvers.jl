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
