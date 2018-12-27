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
