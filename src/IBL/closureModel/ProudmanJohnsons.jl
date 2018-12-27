function correlate(ncell::Int64, sol::Array{Float64})

    E = zeros(ncell+2); B = zeros(ncell+2); F = zeros(ncell+2); S = zeros(ncell+2)
    dfde = zeros(ncell+2); del = zeros(ncell+2)

    for i = 1:ncell+2
        sol.del[i] = sol[1,i]
        E[i] = sol[2,i]/del[i] - 1.
        F[i] = 4.8274*E[i]^4 - 5.9816*E[i]^3 + 4.0274*E[i]^2 + 0.23247*E[i] + 0.15174
        if E[i] < -0.0616
            B[i] = -225.86*E[i]^3 - 3016.6*E[i]^2 - 208.68*E[i] - 17.915
        elseif E[i] >-0.0395
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

end
