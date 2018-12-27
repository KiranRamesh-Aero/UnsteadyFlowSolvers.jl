mutable struct InvisicidTransport

    ue_us    :: Array{Float64,1}
    ue_ls    :: Array{Float64,1}
    ue_us_t0 :: Array{Float64,1}
    ue_ls_t0 :: Array{Float64,1}

    function InvisicidTransport(ncell::Int64)

        ue_us = zeros(ncell+2)
        ue_ls = zeros(ncell+2)
        ue_us_t0 = zeros(ncell+2)
        ue_ls_t0 = zeros(ncell+2)
        new(ue_us,ue_ls,ue_us_t0,ue_ls_t0)
    end

end


mutable struct FluxSplittingParameters

    E :: Array{Float64,1}
    B :: Array{Float64,1}
    F :: Array{Float64,1}
    S :: Array{Float64,1}
    dfde :: Array{Float64,1}
    del :: Array{Float64,1}

    function FluxSplittingParameters(ncell::Int64)

        E=ones(ncell+2)*0.4142
        B=ones(ncell+2)*((131.9*0.4142^3)-(167.32*0.4142^2)+(76.642*0.4142)-11.068)
        F=ones(ncell+2)*(4.8274*0.4142^4 - 5.9816*0.4142^3 + 4.0274*0.4142^2 + 0.23247*0.4142 + 0.15174)
        S=zeros(ncell+2)
        dfde=zeros(ncell+2)
        del=ones(ncell+2).*sqrt.(B[:]*0.005)
        new(E,B,F,S,dfde,del)
    end

end

mutable struct AerofoilThickness

xCoordinates :: Array{Float64,1}
yCoordinates :: Array{Float64,1}

 end

mutable struct Solutions

    sol ::   Array{Float64,2}
    solt ::  Array{Float64,2}
    lamb1 :: Array{Float64,1}
    lamb2 :: Array{Float64,1}
    Jsep ::  Array{Float64,2}
    Csep  :: Array{Float64,1}
    cf :: Array{Float64,1}

    function Solutions(ncell::Int64, fluxSplit::FluxSplittingParameters)

        sol = zeros(2,ncell+2)
        sol[1,:] = fluxSplit.del[:]
        sol[2,:] = fluxSplit.del[:].*(fluxSplit.E[:]+ ones(ncell+2))
        solt = zeros(2,ncell+2)
        lamb1 = zeros(ncell+2)
        lamb2 = zeros(ncell+2)
        Jsep = zeros(2,ncell+2)
        Csep = zeros(ncell+2)
        cf = zeros(ncell+2)
        new(sol,solt,lamb1,lamb2,Jsep,Csep,cf)
    end

end


mutable struct OperationalConditions

    cfl::Float64
    Re::Float64
    dx::Float64
    x:: Array{Float64,1}

    function OperationalConditions(cfl::Float64, Re::Float64,ncell::Int64)
        cfl=cfl
        Re = Re
        dx = 0.0
        x= zeros(ncell+2)
        new(cfl, Re, dx, x)
    end
end
