function identifyFlowSeperation(surf::TwoDSurf,fluxSplit::FluxSplittingParameters,soln::Solutions,ue::Array{Float64,1},x::Array{Float64,1},re::Float64, ncell::Int64)


    for i = 2:ncell+1
        #soln.cf[i] = fluxSplit.B[i]*ue[i]*surf.c/(fluxSplit.del[i]*surf.uref*re)
        soln.Csep[i] = abs(fluxSplit.del[i+1] - fluxSplit.del[i])/abs(fluxSplit.del[i] - fluxSplit.del[i-1])
        if abs(soln.Csep[i]) > 10.
            println("singularity (separation) detected at x=$(x[i])")
        end
    end

end

function boundaryCorrection(u::Array{Float64,1})

    n_bl = length(u)
    u[1]= 2*u[2]-u[3]
    u[n_bl] = 2*u[n_bl-1]-u[n_bl-2]

end

function boundaryCorrection(u::Array{Float64,2})

    for i =1:ndims(u)
        n_bl = length(u[i,:])
        u[i,1]= 2*u[i,2]-u[i,3]
        u[i,n_bl] = 2*u[i,n_bl-1]-u[i,n_bl-2]
    end
end

function derivativesViscous(dA::Array{Float64,1}, dB::Array{Float64,1},ncell::Int64)

    der =zeros(ncell+2)
    if ((ndims(dA)==1) && (ndims(dB)==1))
        if (length(dA)==length(dB))
            for i=2:ncell+1
                der[i] = (dA[i]-dA[i-1])/(dB[i]-dB[i-1])
            end
        else
            error("Input arrays should have identical dimensions")
        end
    else
        error("Function only acccepts one-dimensional arrays")
    end

   der[1]=2*der[2]-der[3]
   der[ncell+2] = 2*der[ncell+1]-der[ncell]
   return der
end

function derivativesViscous(dA::Array{Float64,1},dA0::Array{Float64,1}, dt::Float64, ncell::Int64)

    #n_ele = length(dA)
    der =zeros(ncell+2)
        if ((ndims(dA)==1) && (ndims(dA0)==1))
            if (length(dA)==length(dA0))
                for i=2:ncell+1
                    der[i] = (dA[i]-dA0[i])/(dt)
                end
            else
                error("Input arrays should have identical dimensions")
            end
        else
            error("Function only acccepts one-dimensional arrays")
        end

    der[1]=2*der[2]-der[3]
    der[ncell+2] = 2*der[ncell+1]-der[ncell]
    return der
end

function correlateFunction(fluxSplitPara::FluxSplittingParameters, sol::Array{Float64,2},ncell::Int64)

    for i=1:ncell+2
        fluxSplitPara.del[i] = sol[1,i]
        fluxSplitPara.E[i] = sol[2,i]/(fluxSplitPara.del[i]) - 1.0
        fluxSplitPara.F[i] = 4.8274*fluxSplitPara.E[i]^4 - 5.9816*fluxSplitPara.E[i]^3 + 4.0274*fluxSplitPara.E[i]^2 + 0.23247*fluxSplitPara.E[i] + 0.15174

        if fluxSplitPara.E[i] < -0.0616
            fluxSplitPara.B[i] = -225.86*fluxSplitPara.E[i]^3 - 3016.6*fluxSplitPara.E[i]^2 - 208.68*fluxSplitPara.E[i] - 17.915
        elseif fluxSplitPara.E[i] > -0.0395
            fluxSplitPara.B[i] = 131.9*fluxSplitPara.E[i]^3 - 167.32*fluxSplitPara.E[i]^2 + 76.642*fluxSplitPara.E[i] - 11.068
        else
            fluxSplitPara.B[i] = 0.5*(-225.86*fluxSplitPara.E[i]^3 - 3016.6*fluxSplitPara.E[i]^2 - 208.68*fluxSplitPara.E[i] - 17.915
            + 131.9*fluxSplitPara.E[i]^3 - 167.32*fluxSplitPara.E[i]^2 + 76.642*fluxSplitPara.E[i] - 11.068)
        end
        if fluxSplitPara.E[i] < -0.0582
            fluxSplitPara.S[i] = 451.55*fluxSplitPara.E[i]^3 + 2010*fluxSplitPara.E[i]^2 + 138.96*fluxSplitPara.E[i] + 11.296
        elseif fluxSplitPara.E[i] > -0.042
            fluxSplitPara.S[i] = -96.739*fluxSplitPara.E[i]^3 + 117.74*fluxSplitPara.E[i]^2 - 46.432*fluxSplitPara.E[i] + 6.8074
        else
            fluxSplitPara.S[i] = 0.5*(451.55*fluxSplitPara.E[i]^3 + 2010*fluxSplitPara.E[i]^2 + 138.96*fluxSplitPara.E[i] + 11.296
            - 96.739*fluxSplitPara.E[i]^3 + 117.74*fluxSplitPara.E[i]^2 - 46.432*fluxSplitPara.E[i] + 6.8074)
        end

        fluxSplitPara.dfde[i] = 4*4.8274*fluxSplitPara.E[i]^3 - 3*5.9816*fluxSplitPara.E[i]^2 + 2*4.0274*fluxSplitPara.E[i] + 0.23247
    end

end
