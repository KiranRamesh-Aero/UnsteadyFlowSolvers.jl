type KinemDef3D
    alpha :: MotionDef
    h :: MotionDef
    u :: MotionDef
    vartype :: String
    vary ::Int8
    add :: Array{Float64,2}

    function KinemDef3D(alpha :: MotionDef, h::MotionDef, u::MotionDef, vartype = "Constant", vary = 0, add = zeros(1,1))
        new(alpha, h, u, vartype, vary, add)
    end
end

immutable ThreeDFieldSimple
    f2d :: Vector{TwoDFlowField}
    function ThreeDFieldSimple()
        f2d = TwoDFlowField[]
    new(f2d)
    end
end

immutable ThreeDSurfSimple
    cref :: Float64
    AR :: Float64
    uref :: Float64 
    pvt :: Float64
    lespcrit :: Vector{Float64}
    coord_file :: String
    ndiv :: Int8
    nspan :: Int8
    naterm :: Int8
    kindef :: KinemDef3D
    psi :: Vector{Float64}
    yle :: Vector{Float64}
    s2d :: Vector{TwoDSurf}
    a03d :: Vector{Float64}
    bc :: Vector{Float64}
    nshed :: Vector{Float64}
    bcoeff :: Vector{Float64}
    a03dprev :: Vector{Float64}
    a03ddot :: Vector{Float64}
    levstr :: Vector{Float64}
    
    function ThreeDSurfSimple(AR, kindef, coord_file, pvt, lespcrit = [10.;]; nspan = 10, cref = 1., uref=1., ndiv=70, naterm=35)

        bref = AR*cref

        psi = zeros(nspan)
        yle = zeros(nspan)

        s2d = TwoDSurf[]

	for i = 1:nspan
            psi[i] = real(i)*(pi/2)/nspan
            yle[i] = -bref*cos(psi[i])/2.
        end

        #This code should be made more general to allow more motion types and combinations
        if typeof(kindef.h) == BendingDef
            for i = 1:nspan
                h_amp = evaluate(kindef.h.spl, yle[i])*kindef.h.scale 
                h2d = CosDef(0., h_amp, kindef.h.k, kindef.h.phi)
                kinem2d = KinemDef(kindef.alpha, h2d, kindef.u)
                push!(s2d, TwoDSurf(coord_file, pvt,  kinem2d, lespcrit, c=cref, uref=uref, ndiv=ndiv, naterm=naterm))
            end
        else
            for i = 1:nspan
                kinem2d = KinemDef(kindef.alpha, kindef.h, kindef.u)
                lespc = lespcrit[1]
                push!(s2d, TwoDSurf(coord_file, pvt,  kinem2d, [lespc;], c=cref, uref=uref, ndiv=ndiv, naterm=naterm))
            end
        end
        
        a03d = zeros(nspan)
        a03ddot = zeros(nspan)
        a03dprev = zeros(nspan)
        bc = zeros(nspan)
        nshed = [0.;]
        bcoeff = zeros(nspan)
        levstr = zeros(nspan)

        new(cref, AR, uref, pvt, lespcrit, coord_file,  ndiv, nspan, naterm, kindef, psi, yle, s2d, a03d, bc, nshed, bcoeff, a03dprev, a03ddot, levstr)
end
end

