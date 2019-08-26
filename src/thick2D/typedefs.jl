#Source distribution to model nonlifting/thickness solution
mutable struct TwoDSource
    x :: Float64
    z :: Float64
    s :: Float64
end

mutable struct TwoDSurfThick
    c :: Float64
    uref :: Float64
    coord_file :: String
    pvt :: Float64
    ndiv :: Int16
    naterm :: Int16
    kindef :: KinemDef
    cam :: Vector{Float64}
    cam_slope :: Vector{Float64}
    thick :: Vector{Float64}
    thick_slope :: Vector{Float64}
    theta :: Vector{Float64}
    x :: Vector{Float64}
    kinem :: KinemPar
    bnd_x_u :: Vector{Float64}
    bnd_z_u :: Vector{Float64}
    bnd_x_l :: Vector{Float64}
    bnd_z_l :: Vector{Float64}
    bnd_x_chord :: Vector{Float64}
    bnd_z_chord :: Vector{Float64}
    uind_u :: Vector{Float64}
    uind_l :: Vector{Float64}
    wind_u :: Vector{Float64}
    wind_l :: Vector{Float64}
    downwash :: Vector{Float64}
    a0 :: Vector{Float64}
    aterm :: Vector{Float64}
    a0dot :: Vector{Float64}
    adot :: Vector{Float64}
    a0prev :: Vector{Float64}
    aprev :: Vector{Float64}
    bterm :: Vector{Float64}
    bv :: Vector{TwoDVort}
    src :: Vector{TwoDSource}
    lespcrit :: Vector{Float64}
    levflag :: Vector{Int8}
    initpos :: Vector{Float64}
    rho :: Float64
    LHS :: Array{Float64}
    RHS :: Vector{Float64}

    function TwoDSurfThick(coord_file, pvt, kindef,lespcrit=zeros(1); c=1., uref=1., ndiv=70, naterm=35, initpos = [0.; 0.])
        theta = zeros(ndiv); x = zeros(ndiv); cam = zeros(ndiv); cam_slope = zeros(ndiv)
        thick = zeros(ndiv); thick_slope = zeros(ndiv); bnd_x_u = zeros(ndiv); bnd_z_u = zeros(ndiv)
        bnd_x_l = zeros(ndiv); bnd_z_l = zeros(ndiv); bnd_x_chord = zeros(ndiv); bnd_z_chord = zeros(ndiv)

        kinem = KinemPar(0, 0, 0, 0, 0, 0)

        dtheta = pi/(ndiv-1)
        for ib = 1:ndiv
            theta[ib] = real(ib-1.)*dtheta
            x[ib] = c/2. *(1-cos(theta[ib]))
        end

        thick, thick_slope, rho, cam, cam_slope = camber_thick_calc(x, coord_file)

        kinem = KinemPar(0., 0., 0., 0., 0., 0.)

        if (typeof(kindef.alpha) == EldUpDef)
            kinem.alpha = kindef.alpha(0.)
            kinem.alphadot = ForwardDiff.derivative(kindef.alpha,0.)*uref/c
        elseif (typeof(kindef.alpha) == EldRampReturnDef)
            kinem.alpha = kindef.alpha(0.)
            kinem.alphadot = ForwardDiff.derivative(kindef.alpha,0.)*uref/c
        elseif (typeof(kindef.alpha) == ConstDef)
            kinem.alpha = kindef.alpha(0.)
            kinem.alphadot = 0.
        elseif (typeof(kindef.alpha) == SinDef)
            kinem.alpha = kindef.alpha(0.)
            kinem.alphadot = ForwardDiff.derivative(kindef.alpha,0.)*uref/c
        elseif (typeof(kindef.alpha) == CosDef)
            kinem.alpha = kindef.alpha(0.)
            kinem.alphadot = ForwardDiff.derivative(kindef.alpha,0.)*uref/c
        end
        if (typeof(kindef.h) == EldUpDef)
            kinem.h = kindef.h(0.)*c
            kinem.hdot = ForwardDiff.derivative(kindef.h,0.)*uref
        elseif (typeof(kindef.h) == EldUpIntDef)
            kinem.h = kindef.h(0.)*c
            kinem.hdot = ForwardDiff.derivative(kindef.h,0.)*uref
        elseif (typeof(kindef.h) == EldRampReturnDef)
            kinem.h= kindef.h(0.)*c
            kinem.hdot = ForwardDiff.derivative(kindef.h,0.)*uref
        elseif (typeof(kindef.h) == ConstDef)
            kinem.h = kindef.h(0.)*c
            kinem.hdot = 0.
        elseif (typeof(kindef.h) == SinDef)
            kinem.h = kindef.h(0.)*c
            kinem.hdot = ForwardDiff.derivative(kindef.h,0.)*uref
        elseif (typeof(kindef.h) == CosDef)
            kinem.h = kindef.h(0.)*c
            kinem.hdot = ForwardDiff.derivative(kindef.h,0.)*uref
        end

        if (typeof(kindef.u) == EldUpDef)
            kinem.u = kindef.u(0.)*uref
            kinem.udot = ForwardDiff.derivative(kindef.u,0.)*uref*uref/c
        elseif (typeof(kindef.u) == EldRampReturnDef)
            kinem.u, kinem.udot = kindef.u(0.)
            kinem.u = kinem.u*uref
            kinem.udot = kinem.udot*uref*uref/c
        elseif (typeof(kindef.u) == ConstDef)
            kinem.u = kindef.u(0.)*uref
            kinem.udot = 0.
        end

        for i = 1:ndiv
            #Upper and lower are used to find the induced velocities
            #of wake. Chord is used for placing singular distributions

            bnd_x_u[i] = -((c - pvt*c)+((pvt*c - x[i])*cos(kinem.alpha))) + ((cam[i]+thick[i])*sin(kinem.alpha)) + initpos[1]
            bnd_z_u[i] = kinem.h + ((pvt*c - x[i])*sin(kinem.alpha))+((cam[i]+thick[i])*cos(kinem.alpha)) + initpos[2]
            bnd_x_l[i] = -((c - pvt*c)+((pvt*c - x[i])*cos(kinem.alpha))) + ((cam[i]-thick[i])*sin(kinem.alpha)) + initpos[1]
            bnd_z_l[i] = kinem.h + ((pvt*c - x[i])*sin(kinem.alpha))+((cam[i]-thick[i])*cos(kinem.alpha)) + initpos[2]
            bnd_x_chord[i] = -((c - pvt*c)+((pvt*c - x[i])*cos(kinem.alpha))) + initpos[1]
            bnd_z_chord[i] = kinem.h + ((pvt*c - x[i])*sin(kinem.alpha)) + initpos[2]
        end

        uind_u = zeros(ndiv); uind_l = zeros(ndiv); wind_u = zeros(ndiv); wind_l = zeros(ndiv)
        downwash = zeros(ndiv); a0 = zeros(1); a0dot = zeros(1); aterm = zeros(naterm)
        adot = zeros(naterm); a0prev = zeros(1); aprev = zeros(naterm); bterm = zeros(naterm);

        bv = TwoDVort[]
        src = TwoDSource[]

        for i = 1:ndiv-1
            push!(bv,TwoDVort(0,0,0,0.02*c,0,0))
        end

        for i = 1:ndiv-1
            xsrc = 0.5*(x[i] + x[i+1])
            dx = x[i+1] - x[i]
            thder = 0.5*(thick_slope[i] + thick_slope[i+1])
            push!(src, TwoDSource(xsrc, 0, 2*uref*thder*dx))
        end

        LHS = zeros(2*ndiv-1,naterm*2+2)
        RHS = zeros(2*ndiv-1)

        #Construct constant columns in LHS (all except the last one involving shed vortex)
        for i = 2:ndiv-1
            
             #Sweep all rows (corresponding to ndiv) for lifting equation
            
            #Sweep columns for aterms
            for n = 1:naterm
                LHS[i-1,n] = cos(n*theta[i]) - thick_slope[i]*sin(n*theta[i]) 
            end

            #Sweep columns for bterm
            for n = 1:naterm
                LHS[i-1,n+naterm] = cam_slope[i]*cos(n*theta[i]) 
            end

            #TEV term must be updated in the loop after its location is known
            #Sweep all rows (corresponding to ndiv) for nonlifting equation
         
            for n = 1:naterm
                LHS[ndiv+i-3,n]  = -cam_slope[i]*sin(n*theta[i])
            end
            for n = 1:naterm
                LHS[ndiv+i-3,naterm+n] = sin(n*theta[i]) + thick_slope[i]*cos(n*theta[i]) 
            end
        end

        #Terms for Kelvin condition
        LHS[2*ndiv-3,1] = 100*pi/2

        LHS[2*ndiv-3,2*naterm+1] = 100.
        #LHS[2*ndiv-3,2*naterm+3] = 1.   #FOR LEV

        # #Kutta
        for n = 1:naterm
            LHS[2*ndiv-2,n] = ((-1)^n)*100.
        end

        #Boundary conditions
        for n = 1:naterm
            LHS[2*ndiv-1,n+naterm] = ((-1)^n)*100.
        end
        
        
        # LHS[2*ndiv-1,1] = sqrt(2. /rho) + 1.
        # for n = 1:naterm
        #     LHS[2*ndiv-1,n+1] = -1.
        #     LHS[2*ndiv-1,n+naterm+1] = 1.
        # end
        levflag = [0;]

        new(c, uref, coord_file, pvt, ndiv, naterm, kindef, cam, cam_slope, thick, thick_slope, theta, x, kinem, bnd_x_u, bnd_z_u, bnd_x_l, bnd_z_l, bnd_x_chord, bnd_z_chord, uind_u, uind_l, wind_u, wind_l, downwash, a0, aterm, a0dot, adot, a0prev, aprev, bterm, bv, src, lespcrit, levflag, initpos, rho, LHS, RHS)
    end
end
