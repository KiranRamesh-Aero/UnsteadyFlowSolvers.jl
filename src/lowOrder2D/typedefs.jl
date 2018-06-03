"""
TwoDVort(x,z,s,vc,vx,vz)

Defines a 2D vortex at `(x,z)` with vorticity `s` and vortex-core
radius `vc`.

`vx` and `vz` are induced velocity at centroid used in vortex
interaction calculations

"""
type TwoDVort
    x :: Float64
    z :: Float64
    s :: Float64
    vc :: Float64
    vx :: Float64
    vz :: Float64
end

#Source distribution to model nonlifting/thickness solution
type TwoDSource
    x :: Float64
    z :: Float64
    s :: Float64
end

immutable TwoDFlowField
    velX :: MotionDef
    velZ :: MotionDef
    u :: Vector{Float64}
    w :: Vector{Float64}
    tev :: Vector{TwoDVort}
    lev :: Vector{TwoDVort}
    extv :: Vector{TwoDVort}
    function TwoDFlowField(velX = ConstDef(0.), velZ = ConstDef(0.))
        u = [0;]
        w = [0;]
        tev = TwoDVort[]
        lev = TwoDVort[]
        extv = TwoDVort[]
        new(velX, velZ, u, w, tev, lev, extv)
    end
end

immutable TwoDSurf
    c :: Float64
    uref :: Float64
    coord_file :: String
    pvt :: Float64
    ndiv :: Int8
    naterm :: Int8
    kindef :: KinemDef
    cam :: Vector{Float64}
    cam_slope :: Vector{Float64}
    theta :: Vector{Float64}
    x :: Vector{Float64}
    kinem :: KinemPar
    bnd_x :: Vector{Float64}
    bnd_z :: Vector{Float64}
    uind :: Vector{Float64}
    wind :: Vector{Float64}
    downwash :: Vector{Float64}
    a0 :: Vector{Float64}
    aterm :: Vector{Float64}
    a0dot :: Vector{Float64}
    adot :: Vector{Float64}
    a0prev :: Vector{Float64}
    aprev :: Vector{Float64}
    bv :: Vector{TwoDVort}
    lespcrit :: Vector{Float64}
    levflag :: Vector{Int8}
    initpos :: Vector{Float64}

    function TwoDSurf(coord_file, pvt, kindef,lespcrit=zeros(1); c=1., uref=1., ndiv=70, naterm=35, initpos = [0.; 0.])
        theta = zeros(ndiv)
        x = zeros(ndiv)
        cam = zeros(ndiv)
        cam_slope = zeros(ndiv)
        bnd_x = zeros(ndiv)
        bnd_z = zeros(ndiv)
        kinem = KinemPar(0, 0, 0, 0, 0, 0)

        dtheta = pi/(ndiv-1)
        for ib = 1:ndiv
            theta[ib] = real(ib-1.)*dtheta
            x[ib] = c/2.*(1-cos(theta[ib]))
        end
        if (coord_file != "FlatPlate")
            cam, cam_slope = camber_calc(x, coord_file)
        end

        #camspl = Spline1D(theta, cam)
        #for i = 1:ndiv
        #    cam_slope[i] = Dierckx.derivative(camspl, theta[i])
        #end

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
            bnd_x[i] = -((c - pvt*c)+((pvt*c - x[i])*cos(kinem.alpha))) + (cam[i]*sin(kinem.alpha)) + initpos[1]
            bnd_z[i] = kinem.h + ((pvt*c - x[i])*sin(kinem.alpha))+(cam[i]*cos(kinem.alpha)) + initpos[2]
        end
        uind = zeros(ndiv)
        wind = zeros(ndiv)
        downwash = zeros(ndiv)
        a0 = zeros(1)
        a0dot = zeros(1)
        aterm = zeros(naterm)
        adot = zeros(3)
        a0prev = zeros(1)
        aprev = zeros(3)
        bv = TwoDVort[]
        for i = 1:ndiv-1
            push!(bv,TwoDVort(0,0,0,0.02*c,0,0))
        end
        levflag = [0]
        new(c, uref, coord_file, pvt, ndiv, naterm, kindef, cam, cam_slope, theta, x, kinem, bnd_x, bnd_z, uind, wind, downwash, a0, aterm, a0dot, adot, a0prev, aprev, bv, lespcrit, levflag, initpos)
    end
end



immutable TwoDSurfThick
    c :: Float64
    uref :: Float64
    coord_file :: String
    pvt :: Float64
    ndiv :: Int8
    naterm :: Int8
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
            x[ib] = c/2.*(1-cos(theta[ib]))
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
        adot = zeros(3); a0prev = zeros(1); aprev = zeros(3); bterm = zeros(naterm);
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

        LHS = zeros(ndiv*2-3,naterm*2+2)
        RHS = zeros(ndiv*2-3)

        #Construct constant columns in LHS (all except the last one involving shed vortex)
        for i = 2:ndiv-1

            #Sweep all rows (corresponding to ndiv) for lifting equation
            #A0 term
            LHS[i-1,1] = uref*(-1 + thick[i]/(c*sin(theta[i])*sin(theta[i]/2)^2) - thick_slope[i]*cot(theta[i]/2))

            #Sweep columns for aterms
            for n = 1:naterm
                LHS[i-1,n+1] = uref*(cos(n*theta[i]) - thick_slope[i]*sin(n*theta[i]) + 2*thick[i]*n*cos(n*theta[i])/(c*sin(theta[i])))
            end

            #Sweep columns for bterm
            for n = 1:naterm
                LHS[i-1,n+naterm+1] = uref*(cam_slope[i]*cos(n*theta[i]) - 2*cam[i]*n*sin(n*theta[i])/(c*sin(theta[i])))
            end
            #TEV term must be updated in the loop after its location is known
            #Sweep all rows (corresponding to ndiv) for nonlifting equation
            LHS[ndiv+i-3,1]  = uref*(cam[i]/(c*sin(theta[i])*sin(theta[i]/2)^2) - cam_slope[i]*cot(theta[i]/2))
            for n = 1:naterm
                LHS[ndiv+i-3,1+n]  = -uref*(cam_slope[i]*sin(n*theta[i]) + 2*cam[i]*n*cos(n*theta[i])/(c*sin(theta[i])))
            end
            for n = 1:naterm
                LHS[ndiv+i-3,1+naterm+n] = uref*(sin(n*theta[i]) + thick_slope[i]*cos(n*theta[i]) - 2*thick[i]*n*sin(n*theta[i])/(c*sin(theta[i])))
            end
        end

        #Terms for Kelvin condition
        LHS[2*ndiv-3,1] = uref*c*pi
        LHS[2*ndiv-3,2] = uref*c*pi/2
        LHS[2*ndiv-3,2*naterm+2] = 1.

        levflag = [0]
        new(c, uref, coord_file, pvt, ndiv, naterm, kindef, cam, cam_slope, thick, thick_slope, theta, x, kinem, bnd_x_u, bnd_z_u, bnd_x_l, bnd_z_l, bnd_x_chord, bnd_z_chord, uind_u, uind_l, wind_u, wind_l, downwash, a0, aterm, a0dot, adot, a0prev, aprev, bterm, bv, src, lespcrit, levflag, initpos, rho, LHS, RHS)
    end
end

type TwoDOFPar
    x_alpha :: Float64
    r_alpha :: Float64
    kappa :: Float64
    w_alpha :: Float64
    w_h :: Float64
    w_alphadot :: Float64
    w_hdot :: Float64
    cubic_h_1 :: Float64
    cubic_h_3 :: Float64
    cubic_alpha_1 :: Float64
    cubic_alpha_3 :: Float64
end

#For airfoil with 2DOF in pitch and plunge
type KinemPar2DOF
    alpha :: Float64
    h :: Float64
    alphadot :: Float64
    hdot :: Float64
    u :: Float64
    udot :: Float64
    alphaddot :: Float64
    hddot :: Float64
    alpha_pr :: Float64
    h_pr :: Float64
    alphadot_pr :: Float64
    alphadot_pr2 :: Float64
    alphadot_pr3 :: Float64
    hdot_pr :: Float64
    hdot_pr2 :: Float64
    hdot_pr3 :: Float64
    alphaddot_pr :: Float64
    alphaddot_pr2 :: Float64
    alphaddot_pr3 :: Float64
    hddot_pr :: Float64
    hddot_pr2 :: Float64
    hddot_pr3 :: Float64

    function KinemPar2DOF(alpha, h, alphadot, hdot, u)
        new(alpha, h, alphadot, hdot, u, 0., 0., 0., alpha, h, alphadot, alphadot, alphadot, hdot, hdot, hdot, 0., 0., 0., 0., 0., 0.)
    end
end

immutable TwoDSurf2DOF
    c :: Float64
    uref :: Float64
    coord_file :: String
    pvt :: Float64
    ndiv :: Int8
    naterm :: Int8
    cam :: Vector{Float64}
    cam_slope :: Vector{Float64}
    theta :: Vector{Float64}
    x :: Vector{Float64}
    strpar :: TwoDOFPar
    kinem :: KinemPar2DOF
    bnd_x :: Vector{Float64}
    bnd_z :: Vector{Float64}
    uind :: Vector{Float64}
    wind :: Vector{Float64}
    downwash :: Vector{Float64}
    a0 :: Vector{Float64}
    aterm :: Vector{Float64}
    a0dot :: Vector{Float64}
    adot :: Vector{Float64}
    a0prev :: Vector{Float64}
    aprev :: Vector{Float64}
    bv :: Vector{TwoDVort}
    lespcrit :: Vector{Float64}
    levflag :: Vector{Int8}

    function TwoDSurf2DOF(c, uref, coord_file, pvt, strpar, kinem ,lespcrit=zeros(1), ndiv=70, naterm=45)
        theta = zeros(ndiv)
        x = zeros(ndiv)
        cam = zeros(ndiv)
        cam_slope = zeros(ndiv)
        bnd_x = zeros(ndiv)
        bnd_z = zeros(ndiv)

        dtheta = pi/(ndiv-1)
        for ib = 1:ndiv
            theta[ib] = real(ib-1.)*dtheta
            x[ib] = c/2.*(1-cos(theta[ib]))
        end
        if (coord_file != "FlatPlate")
            cam, cam_slope = camber_calc(x, coord_file)
        end

        for i = 1:ndiv
            bnd_x[i] = -((c - pvt*c)+((pvt*c - x[i])*cos(kinem.alpha))) + (cam[i]*sin(kinem.alpha))
            bnd_z[i] = kinem.h + ((pvt*c - x[i])*sin(kinem.alpha))+(cam[i]*cos(kinem.alpha))
        end
        uind = zeros(ndiv)
        wind = zeros(ndiv)
        downwash = zeros(ndiv)
        a0 = zeros(1)
        a0dot = zeros(1)
        aterm = zeros(naterm)
        adot = zeros(3)
        a0prev = zeros(1)
        aprev = zeros(3)
        bv = TwoDVort[]
        for i = 1:ndiv-1
            push!(bv,TwoDVort(0,0,0,0.02*c,0,0))
        end
        levflag = [0]
        new(c, uref, coord_file, pvt, ndiv, naterm, cam, cam_slope, theta, x, strpar, kinem, bnd_x, bnd_z, uind, wind, downwash, a0, aterm, a0dot, adot, a0prev, aprev, bv,lespcrit,levflag)
    end
end

immutable KelvinCondition
    surf :: TwoDSurf
    field :: TwoDFlowField
end

function (kelv::KelvinCondition)(tev_iter::Array{Float64})
    #Update the TEV strength
    nlev = length(kelv.field.lev)
    ntev = length(kelv.field.tev)
    kelv.field.tev[ntev].s = tev_iter[1]

    #Update incduced velocities on airfoil
    update_indbound(kelv.surf, kelv.field)

    #Calculate downwash
    update_downwash(kelv.surf, [kelv.field.u[1],kelv.field.w[1]])

    #Calculate first two fourier coefficients
    update_a0anda1(kelv.surf)

    val = kelv.surf.uref*kelv.surf.c*pi*(kelv.surf.a0[1] + kelv.surf.aterm[1]/2.)

    for iv = 1:ntev
        val = val + kelv.field.tev[iv].s
    end
    for iv = 1:nlev
        val = val + kelv.field.lev[iv].s
    end

    #Add kelv_enforced if necessary - merging will be better
    # val is the value of Gam_b + sum Gam_tev + Gam_lev which will equal zero
    # if the condition is satified

    return val
end

immutable KelvinCondition2DOF
    surf :: TwoDSurf2DOF
    field :: TwoDFlowField
end

function (kelv::KelvinCondition2DOF)(tev_iter::Array{Float64})
    #Update the TEV strength
    nlev = length(kelv.field.lev)
    ntev = length(kelv.field.tev)
    kelv.field.tev[ntev].s = tev_iter[1]

    #Update incduced velocities on airfoil
    update_indbound(kelv.surf, kelv.field)

    #Calculate downwash
    update_downwash(kelv.surf, [kelv.field.u[1],kelv.field.w[1]])

    #Calculate first two fourier coefficients
    update_a0anda1(kelv.surf)

    val = kelv.surf.uref*kelv.surf.c*pi*(kelv.surf.a0[1] + kelv.surf.aterm[1]/2.)

    for iv = 1:ntev
        val = val + kelv.field.tev[iv].s
    end
    for iv = 1:nlev
        val = val + kelv.field.lev[iv].s
    end

    return val
end

immutable KelvinKutta
    surf :: TwoDSurf
    field :: TwoDFlowField
end

function (kelv::KelvinKutta)(v_iter::Array{Float64})
    val = zeros(2)

    #Update the TEV and LEV strengths
    nlev = length(kelv.field.lev)
    ntev = length(kelv.field.tev)
    kelv.field.tev[ntev].s = v_iter[1]
    kelv.field.lev[nlev].s = v_iter[2]

    #Update incduced velocities on airfoil
    update_indbound(kelv.surf, kelv.field)

    #Calculate downwash
    update_downwash(kelv.surf ,[kelv.field.u[1],kelv.field.w[1]])

    #Calculate first two fourier coefficients
    update_a0anda1(kelv.surf)

    val[1] = kelv.surf.uref*kelv.surf.c*pi*(kelv.surf.a0[1] + kelv.surf.aterm[1]/2.)

    for iv = 1:ntev
        val[1] = val[1] + kelv.field.tev[iv].s
    end
    for iv = 1:nlev
        val[1] = val[1] + kelv.field.lev[iv].s
    end

    if (kelv.surf.a0[1] > 0)
        lesp_cond = kelv.surf.lespcrit[1]
    else
        lesp_cond = -kelv.surf.lespcrit[1]
    end
    val[2] = kelv.surf.a0[1]-lesp_cond

    #Add kelv_enforced if necessary - merging will be better
    return val
end

immutable KelvinKutta2DOF
    surf :: TwoDSurf2DOF
    field :: TwoDFlowField
end

function (kelv::KelvinKutta2DOF)(v_iter::Array{Float64})
    val = zeros(2)

    #Update the TEV and LEV strengths
    nlev = length(kelv.field.lev)
    ntev = length(kelv.field.tev)
    kelv.field.tev[ntev].s = v_iter[1]
    kelv.field.lev[nlev].s = v_iter[2]

    #Update incduced velocities on airfoil
    update_indbound(kelv.surf, kelv.field)

    #Calculate downwash
    update_downwash(kelv.surf ,[kelv.field.u[1],kelv.field.w[1]])

    #Calculate first two fourier coefficients
    update_a0anda1(kelv.surf)

    val[1] = kelv.surf.uref*kelv.surf.c*pi*(kelv.surf.a0[1] + kelv.surf.aterm[1]/2.)

    for iv = 1:ntev
        val[1] = val[1] + kelv.field.tev[iv].s
    end
    for iv = 1:nlev
        val[1] = val[1] + kelv.field.lev[iv].s
    end

    if (kelv.surf.a0[1] > 0)
        lesp_cond = kelv.surf.lespcrit[1]
    else
        lesp_cond = -kelv.surf.lespcrit[1]
    end
    val[2] = kelv.surf.a0[1]-lesp_cond

    #Add kelv_enforced if necessary - merging will be better
    return val
end
