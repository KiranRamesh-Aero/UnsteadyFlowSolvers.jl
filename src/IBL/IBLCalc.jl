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
