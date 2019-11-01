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
    ate :: Vector{Float64}
    aterm :: Vector{Float64}
    atedot :: Vector{Float64}
    adot :: Vector{Float64}
    ateprev :: Vector{Float64}
    aprev :: Vector{Float64}
    bterm :: Vector{Float64}
    bv :: Vector{TwoDVort}
    src :: Vector{TwoDSource}
    lespcrit :: Vector{Float64}
    levflag :: Vector{Int8}
    initpos :: Vector{Float64}
    LHS :: Array{Float64}
    RHS :: Vector{Float64}

    function TwoDSurfThick(coord_file, pvt, kindef,lespcrit=zeros(1); c=1., uref=1., ndiv=140, naterm=138, initpos = [0.; 0.])
        theta = zeros(ndiv); x = zeros(ndiv); cam = zeros(ndiv); cam_slope = zeros(ndiv)
        thick = zeros(ndiv); thick_slope = zeros(ndiv); bnd_x_u = zeros(ndiv); bnd_z_u = zeros(ndiv)
        bnd_x_l = zeros(ndiv); bnd_z_l = zeros(ndiv); bnd_x_chord = zeros(ndiv); bnd_z_chord = zeros(ndiv)
        
        kinem = KinemPar(0, 0, 0, 0, 0, 0)

        dtheta = pi/(ndiv-1)
        for ib = 1:ndiv
            theta[ib] = real(ib-1.)*dtheta
            x[ib] = c/2. *(1-cos(theta[ib]))
        end

        thick, thick_slope, _, cam, cam_slope = camber_thick_calc(x, coord_file)

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

        if ndiv != naterm+2
            println("Warning : 'ndiv=naterm+2' is required, modifying naterm from $(naterm)to $(ndiv-2)")
            naterm = ndiv-2
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
        downwash = zeros(ndiv); ate = zeros(1); atedot = zeros(1); aterm = zeros(naterm)
        adot = zeros(naterm); ateprev = zeros(1); aprev = zeros(naterm); bterm = zeros(naterm);

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

        LHS = zeros(2*ndiv-2,naterm*2+2)
        RHS = zeros(2*ndiv-2)

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
            LHS[i-1,2*naterm+2] = 1. - thick_slope[i]*tan(theta[i]/2)
            
            #TEV term must be updated in the loop after its location is known
            #Sweep all rows (corresponding to ndiv) for nonlifting equation
         
            for n = 1:naterm
                LHS[ndiv+i-3,n]  = -cam_slope[i]*sin(n*theta[i])
            end
            for n = 1:naterm
                LHS[ndiv+i-3,naterm+n] = sin(n*theta[i]) + thick_slope[i]*cos(n*theta[i])            
            end
            LHS[ndiv+i-3,2*naterm+2] = -cam_slope[i]*tan(theta[i]/2)
        end
        
        #Kelvin condition
        for n = 1:naterm
            s1 = 0; s2 = 0; s3 = 0; s4 = 0; s5 = 0; s6 = 0; s7 = 0; s8 = 0
            for i = 2:ndiv-1
                ds = sqrt((x[i] - x[i-1])^2 + (cam[i] + thick[i] - cam[i-1] - thick[i-1])^2)
                den_u = sqrt(1. + (cam_slope[i] + thick_slope[i])^2)
                den_u_p = sqrt(1. + (cam_slope[i-1] + thick_slope[i-1])^2)
                den_l = sqrt(1. + (cam_slope[i] - thick_slope[i])^2)
                den_l_p = sqrt(1. + (cam_slope[i-1] - thick_slope[i-1])^2)
                
                s1 += 0.5*ds*(sin(n*theta[i])/den_u + sin(n*theta[i-1])/den_u_p)
                s2 += 0.5*ds*((cam_slope[i] + thick_slope[i])*cos(n*theta[i])/den_u + (cam_slope[i-1] + thick_slope[i-1])*cos(n*theta[i-1])/den_u_p)
                s5 += 0.5*ds*(-cos(n*theta[i])/den_u - cos(n*theta[i-1])/den_u_p)
                s6 += 0.5*ds*((cam_slope[i] + thick_slope[i])*sin(n*theta[i])/den_u + (cam_slope[i-1] + thick_slope[i-1])*sin(n*theta[i-1])/den_u_p)
            
                ds = sqrt((x[i] - x[i-1])^2 + (cam[i] - thick[i] - cam[i-1] + thick[i-1])^2)
                s3 += 0.5*ds*(-sin(n*theta[i])/den_l - sin(n*theta[i-1])/den_l_p)
                s4 += 0.5*ds*((cam_slope[i] - thick_slope[i])*cos(n*theta[i])/den_l + (cam_slope[i-1] - thick_slope[i-1])*cos(n*theta[i-1])/den_l_p)
                s7 += 0.5*ds*(-cos(n*theta[i])/den_l - cos(n*theta[i-1])/den_l_p)
                s8 += 0.5*ds*(-(cam_slope[i] - thick_slope[i])*sin(n*theta[i])/den_l - (cam_slope[i-1] - thick_slope[i-1])*sin(n*theta[i-1])/den_l_p)
            end
            LHS[2*ndiv-3,n] = uref*(s1 + s2 - s3 - s4)
            LHS[2*ndiv-3,n+naterm] = uref*(s5 + s6 - s7 - s8)
        end
        s9 = 0; s10 = 0; s11 = 0; s12 = 0
        for i = 2:ndiv-1
            ds = sqrt((x[i] - x[i-1])^2 + (cam[i] + thick[i] - cam[i-1] - thick[i-1])^2)
            den_u = sqrt(1. + (cam_slope[i] + thick_slope[i])^2)
            den_u_p = sqrt(1. + (cam_slope[i-1] + thick_slope[i-1])^2)
            den_l = sqrt(1. + (cam_slope[i] - thick_slope[i])^2)
            den_l_p = sqrt(1. + (cam_slope[i-1] - thick_slope[i-1])^2)
            s9 += 0.5*ds*(tan(theta[i]/2)/den_u + tan(theta[i-1]/2)/den_u_p)
            s10 += 0.5*ds*((cam_slope[i] + thick_slope[i])/den_u + (cam_slope[i-1] + thick_slope[i-1])/den_u_p)
            ds = sqrt((x[i] - x[i-1])^2 + (cam[i] - thick[i] - cam[i-1] + thick[i-1])^2)     
            s11 += 0.5*ds*(-tan(theta[i]/2)/den_l - tan(theta[i-1]/2)/den_l_p)
            s12 += 0.5*ds*((cam_slope[i] - thick_slope[i])/den_l + (cam_slope[i-1] - thick_slope[i-1])/den_l_p)
        end
        LHS[2*ndiv-3,2*naterm+2] = uref*(s9 + s10 - s11 - s12)
        
        levflag = [0;]
        
        new(c, uref, coord_file, pvt, ndiv, naterm, kindef, cam, cam_slope, thick, thick_slope, theta, x, kinem, bnd_x_u, bnd_z_u, bnd_x_l, bnd_z_l, bnd_x_chord, bnd_z_chord, uind_u, uind_l, wind_u, wind_l, downwash, ate, aterm, atedot, adot, ateprev, aprev, bterm, bv, src, lespcrit, levflag, initpos, LHS, RHS)
    end
end


struct KelvinThick
    surf :: TwoDSurfThick
    field :: TwoDFlowField
    bc_prev :: Float64
    dt :: Float64
    aorig :: Array{Float64}
    borig :: Array{Float64}
end

function (kelv::KelvinThick)(tev_iter::Array{Float64})

    nlev = length(kelv.field.lev)
    ntev = length(kelv.field.tev)

    #Update the TEV strength and induced velocities
    minus_indbound_lasttev(kelv.surf, kelv.field)
    kelv.field.tev[ntev].s = tev_iter[1]
    add_indbound_lasttev(kelv.surf, kelv.field)

    update_LHSRHS_kutta(kelv.surf, kelv.field, kelv.dt, tev_iter[1], kelv.aorig, kelv.borig)
    
    #Update fourier coefficients
    soln = kelv.surf.LHS \ kelv.surf.RHS
    kelv.surf.aterm[:] = soln[1:kelv.surf.naterm]
    kelv.surf.bterm[:] = soln[kelv.surf.naterm+1:2*kelv.surf.naterm]
    kelv.surf.ate[1] = soln[2*kelv.surf.naterm+1]
    
    kelv.surf.RHS[2*kelv.surf.ndiv-3] += tev_iter[1]
    
    phi_u, phi_l = calc_phi(kelv.surf)
    bc = phi_u - phi_l
    
    val = bc - kelv.bc_prev + tev_iter[1]

    return val
end
