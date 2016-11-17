#Constants

const accl_g = 9.80665

# ---------------------------------------------------------------------------------------------
# Theodorsen solver input
type TheoDef
    alpha_amp :: Float64
    h_amp :: Float64
    alpha_mean :: Float64
    alpha_zl :: Float64
    k :: Float64
    phi :: Float64
    pvt :: Float64
end

type DelVortDef
    flag :: Int8
    limit :: Int16
    dist :: Float64
end

type TheoDefwFlap
    alpha_amp :: Float64
    h_amp :: Float64
    alpha_mean :: Float64
    alpha_zl :: Float64
    k :: Float64
    phi :: Float64
    pvt :: Float64
    beta_amp :: Float64
    xf :: Float64
    psi :: Float64
end
# ---------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------
# Data structures used to represent the airfoil kinematics
type KinemPar
    alpha :: Float64
    h :: Float64
    alphadot :: Float64
    hdot :: Float64
    u :: Float64
    udot :: Float64
end

# Added by Laura Merchant 2016
type KinemParwFlap
    alpha :: Float64
    h :: Float64
    alphadot :: Float64
    hdot :: Float64
    u :: Float64
    udot :: Float64
    n :: Float64
    ndot :: Float64
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
    alpha_pr2 :: Float64
    alpha_pr3 :: Float64
    h_pr :: Float64
    h_pr2 :: Float64
    h_pr3 :: Float64
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
end

type KinemPar2DFree
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
    u_pr :: Float64
    udot_pr :: Float64
    udot_pr2 ::Float64
    udot_pr3 :: Float64
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

type TwoDFreePar
    r_g :: Float64
    x_g :: Float64
    kappa :: Float64
end


# ---------------------------------------------------------------------------------------------
# Various motion kinematics definitions
abstract MotionDef

immutable EldUpDef <: MotionDef
    amp :: Float64
    K :: Float64
    a :: Float64
end

immutable EldUptstartDef <: MotionDef
    amp :: Float64
    K :: Float64
    a :: Float64
    tstart :: Float64
end

immutable EldRampReturnDef <: MotionDef
    amp :: Float64
    K :: Float64
    a :: Float64
end

immutable ConstDef <: MotionDef
    amp :: Float64
end


function (eld::EldUpDef)(t)
    sm = pi*pi*eld.K/(2*(eld.amp)*(1 - eld.a))
    t1 = 1.
    t2 = t1 + ((eld.amp)/(2*eld.K))
    ((eld.K/sm)*log(cosh(sm*(t - t1))/cosh(sm*(t - t2))))+(eld.amp/2)
end

function (eld::EldUptstartDef)(t)
    sm = pi*pi*eld.K/(2*(eld.amp)*(1 - eld.a))
    t1 = eld.tstart
    t2 = t1 + ((eld.amp)/(2*eld.K))
    ((eld.K/sm)*log(cosh(sm*(t - t1))/cosh(sm*(t - t2))))+(eld.amp/2)
end

function (cons::ConstDef)(t)
    cons.amp
end

function (eld::EldRampReturnDef)(tt)
    fr = eld.K/(pi*abs(eld.amp));
    t1 = 1.
    t2 = t1 + (1./(2*pi*fr));
    t3 = t2 + ((1/(4*fr)) - (1/(2*pi*fr)));
    t4 = t3 + (1./(2*pi*fr));
    t5 = t4+1.;

    nstep = round(Int,t5/0.015) + 1
    g = zeros(nstep)
    t = zeros(nstep)

    for i = 1:nstep
        t[i] = (i-1.)*0.015
        g[i] = log((cosh(eld.a*(t[i] - t1))*cosh(eld.a*(t[i] - t4)))/(cosh(eld.a*(t[i] - t2))*cosh(eld.a*(t[i] - t3))))
    end
    maxg = maximum(g);

    gg = log((cosh(eld.a*(tt - t1))*cosh(eld.a*(tt - t4)))/(cosh(eld.a*(tt - t2))*cosh(eld.a*(tt - t3))))

    return eld.amp*gg/(maxg);

end

immutable SinDef <: MotionDef
  mean :: Float64
  amp :: Float64
  k :: Float64
  phi :: Float64
end

immutable CosDef <: MotionDef
  mean :: Float64
  amp :: Float64
  k :: Float64
  phi :: Float64
end

function (kin::SinDef)(t)
  (kin.mean) + (kin.amp)*sin(2*kin.k*t + kin.phi)
end

function (kin::CosDef)(t)
  (kin.mean) + (kin.amp)*cos(2*kin.k*t + kin.phi)
end
# ---------------------------------------------------------------------------------------------

# Added by Laura Merchant 2016
type KinemDefwFlap
    alpha :: MotionDef
    h :: MotionDef
    u :: MotionDef
    n :: MotionDef
end

type KinemDef
    alpha :: MotionDef
    h :: MotionDef
    u :: MotionDef
end

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

# ---------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------
#Definition of a 2D vortex blob
type TwoDVort
    x :: Float64
    z :: Float64
    s :: Float64
    vc :: Float64
    vx :: Float64
    vz :: Float64
end
# ---------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------
#Definition of a 3D vortex blob
type ThreeDVort
    x :: Vector{Float64}
    s :: Vector{Float64}
    vc :: Float64
    vx :: Float64
    vy :: Float64
    vz :: Float64
end
# ---------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------
#Definition of flowfield in 2D
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

#Definition of flowfield in 3D
immutable ThreeDFlowField
    velX :: MotionDef
    velY :: MotionDef
    velZ :: MotionDef
    u :: Vector{Float64}
    v :: Vector{Float64}
    w :: Vector{Float64}
    tev :: Vector{ThreeDVort}
    lev :: Vector{ThreeDVort}
    extv :: Vector{ThreeDVort}
    function ThreeDFlowField(velX = ConstDef(0.), velY = ConstDef(0.), velZ = ConstDef(0.))
        u = [0;]
        v = [0;]
        w = [0;]
        tev = ThreeDVort[]
        lev = ThreeDVort[]
        extv = ThreeDVort[]
        new(velX, velY, velZ, u, v, w, tev, lev, extv)
    end
end


immutable TwoDFlowData
    tev :: Vector{TwoDVort}
    lev :: Vector{TwoDVort}
    extv :: Vector{TwoDVort}
    bv :: Vector{TwoDVort}
end

# ---------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------
# BEGIN TwoDSurfwFlap
# Added by Laura Merchant 2016
immutable TwoDSurfwFlap
    c :: Float64
    uref :: Float64
    coord_file :: String
    pvt :: Float64
    ndiv :: Int8
    naterm :: Int8
    dynamics_type :: String
    kindef :: KinemDefwFlap
    cam_af :: Vector{Float64}
    cam :: Vector{Float64}
    cam_slope :: Vector{Float64}
    cam_tder :: Vector{Float64}
    theta :: Vector{Float64}
    x :: Vector{Float64}
    x_b :: Vector{Float64}
    kinem :: KinemParwFlap
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
    bv_prev :: Vector{TwoDVort}
    lespcrit :: Vector{Float64}
    levflag :: Vector{Int8}

    function TwoDSurfwFlap(c, uref, coord_file, pvt, ndiv, naterm, dynamics_type, kindef, x_b=1., lespcrit=zeros(1))
        x_b[1] = x_b[1]*c #Dimensional value
        theta = zeros(ndiv)
        x = zeros(ndiv)
        cam = zeros(ndiv)
        cam_af = zeros(ndiv)
        cam_slope = zeros(ndiv)
        cam_tder = zeros(ndiv)
        bnd_x = zeros(ndiv)
        bnd_z = zeros(ndiv)
        kinem = KinemParwFlap(0, 0, 0, 0, 0, 0, 0, 0)

        dtheta = pi/(ndiv-1)
        for ib = 1:ndiv
            theta[ib] = real(ib-1.)*dtheta
            x[ib] = c/2.*(1-cos(theta[ib]))
        end
        if (coord_file != "FlatPlate")
            cam, cam_slope = camber_calc(x, coord_file)
        end

        # Calculates the initial condition for airfoil pitch
        # Various options for prescribing the airfoil pitch kinematics
        if (typeof(kindef.alpha) == EldUpDef)
            kinem.alpha = kindef.alpha(0.)
            kinem.alphadot = ForwardDiff.derivative(kindef.alpha,0.)*uref/c
        elseif (typeof(kindef.alpha) == EldUptstartDef)
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
        # ---------------------------------------------------------------------------------------------

        # Calculates the initial condition for airfoil plunge
        # Various options for prescribing the airfoil plunge kinematics
        if (typeof(kindef.h) == EldUpDef)
            kinem.h = kindef.h(0.)*c
            kinem.hdot = ForwardDiff.derivative(kindef.h,0.)*uref
        elseif (typeof(kindef.h) == EldUptstartDef)
            kinem.h = kindef.h(0.)*c
            kinem.hdot = ForwardDiff.derivative(kindef.h,0.)*uref
        elseif (typeof(kindef.h) == EldUpIntDef)
            kinem.h = kindef.h(0.)*c
            kinem.hdot = ForwardDiff.derivative(kindef.h,0.)*uref
        elseif (typeof(kindef.h) == EldUpInttstartDef)
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
        # ---------------------------------------------------------------------------------------------

        # Calculates the initial condition for forward velocity of airfoil
        # Various options for prescribing the airfoil forward motion
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
        # ---------------------------------------------------------------------------------------------

        # ---------------------------------------------------------------------------------------------
        # Calculates the initial condition for flap deflection beta (n)
        # Various options for prescribing the flap kinematics
        if (typeof(kindef.n) == EldUpDef)
            kinem.n = kindef.n(0.)
            kinem.ndot = ForwardDiff.derivative(kindef.n,0.)*uref/c
        elseif (typeof(kindef.n) == EldRampReturnDef)
            kinem.n = kindef.n(0.)
            kinem.n = ForwardDiff.derivative(kindef.n,0.)*uref/c
        elseif (typeof(kindef.n) == ConstDef)
            kinem.n = kindef.n(0.)
            kinem.n = 0.
        elseif (typeof(kindef.n) == SinDef)
            kinem.n = kindef.n(0.)
            kinem.ndot = ForwardDiff.derivative(kindef.n,0.)*uref/c
        elseif (typeof(kindef.n) == CosDef)
            kinem.n = kindef.n(0.)
            kinem.ndot = ForwardDiff.derivative(kindef.n,0.)*uref/c
        end
# ---------------------------------------------------------------------------------------------


# Defines the aerofoil as a single flat plate
# Calls the calculation of the camber and it's derivatives (spatial and temporal)
#Deflection of flap is modelled as camber variation
if (coord_file != "FlatPlate")
    (cam_af, cam_slope) = camber_calc(x, coord_file)
end

for i = 1:ndiv
    if x[i] < x_b[1]
        cam[i] = cam_af[i];
    else
        cam[i] = cam_af[i] + (x_b[1] - x[i])*tan(kinem.n);
        cam_tder[i] = (x_b[1] - x[i])*kinem.ndot*sec(kinem.n)*sec(kinem.n);
    end
end
# ---------------------------------------------------------------------------------------------

# Populates the arrays for bnd_x, bnd_z based on initial conditions
for i = 1:ndiv
    bnd_x[i] = -((c - pvt*c)+((pvt*c - x[i])*cos(kinem.alpha))) + (cam[i]*sin(kinem.alpha))
    bnd_z[i] = kinem.h + ((pvt*c - x[i])*sin(kinem.alpha))+(cam[i]*cos(kinem.alpha))
end
# ---------------------------------------------------------------------------------------------

# Defines the arrays for various parameters such that they can be later populated
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
bv_prev = TwoDVort[]
for i = 1:ndiv-1
    push!(bv,TwoDVort(0,0,0,0.02*c,0,0))
    push!(bv_prev,TwoDVort(0,0,0,0.02*c,0,0))
end
# ---------------------------------------------------------------------------------------------
levflag = [0]

new(c, uref, coord_file, pvt, ndiv, naterm, dynamics_type, kindef, cam_af, cam, cam_slope, cam_tder, theta, x, x_b, kinem, bnd_x, bnd_z, uind, wind, downwash, a0, aterm, a0dot, adot, a0prev, aprev, bv, bv_prev, lespcrit,levflag)
end
end
# ---------------------------------------------------------------------------------------------
# END TwoDSurfwFlap

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

    function TwoDSurf(coord_file, pvt, kindef,lespcrit=zeros(1), c=1., uref=1., ndiv=70, naterm=35)
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
        new(c, uref, coord_file, pvt, ndiv, naterm, kindef, cam, cam_slope, theta, x, kinem, bnd_x, bnd_z, uind, wind, downwash, a0, aterm, a0dot, adot, a0prev, aprev, bv,lespcrit,levflag)
    end
end


immutable TwoDSurf_2DOF
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

    function TwoDSurf_2DOF(c, uref, coord_file, pvt, ndiv, naterm, strpar, kinem ,lespcrit=zeros(1))
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

immutable TwoDFreeSurf
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
    strpar :: TwoDFreePar
    kinem :: KinemPar2DFree
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

    function TwoDFreeSurf(c, uref, coord_file, pvt, ndiv, naterm, strpar, kinem ,lespcrit=zeros(1))
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


#----------------------------------------------------------------------------------------------
# Definition of a 3D surface
immutable patch
    x :: Float64
    y :: Float64
    z :: Float64
    pvt :: Float64
    coord_file :: String
    c :: Float64
    twist :: Float64
    lc :: Float64
    nspan :: Int8
end
    
immutable ThreeDSurf
    cref :: Float64
    bref :: Float64
    sref :: Float64
    uref :: Float64
    patchdata :: Vector{patch}
    ndiv :: Int8
    nspan :: Int8
    naterm :: Int8
    nbterm :: Int8
    kindef :: KinemDef3D
    c :: Vector{Float64}
    pvt :: Vector{Float64}
    cam :: Array{Float64,2}
    cam_slope :: Array{Float64,2}
    theta :: Vector{Float64}
    psi :: Vector{Float64}
    x :: Array{Float64,2}
    xle :: Vector{Float64}
    yle :: Vector{Float64}
    zle :: Vector{Float64}
    kinem :: Vector{KinemPar}
    bnd_x :: Array{Float64,2}
    bnd_z :: Array{Float64,2}
    uind :: Array{Float64,2}
    vind :: Array{Float64,2}
    wind :: Array{Float64,2}
    downwash :: Array{Float64,2}
    a0 :: Vector{Float64}
    aterm :: Array{Float64,2}
    a0dot :: Vector{Float64}
    adot :: Array{Float64,2}
    a0prev :: Vector{Float64}
    aprev :: Array{Float64,2}
    bv :: Array{ThreeDVort,2}
    lespcrit :: Vector{Float64}
    levflag :: Vector{Int8}
    
    function ThreeDSurf(cref, bref, sref, patchdata, kindef, uref=1., ndiv=70, naterm=35, nbterm = 70)
        nspan = 0
        for i = 1:length(patchdata)
            nspan += patchdata[i].nspan
        end
        nspan = nspan + 1
        lespcrit = zeros(nspan)

        psi = zeros(nspan)
        cam = zeros(nspan,ndiv)
        cam_slope = zeros(nspan,ndiv)
        cam1 = zeros(ndiv)
        cam_slope1 = zeros(ndiv)
        cam2 = zeros(ndiv)
        cam_slope2 = zeros(ndiv)
        bnd_x = zeros(nspan,ndiv)
        bnd_z = zeros(nspan,ndiv)
        c = zeros(nspan)
        pvt = zeros(nspan)
	yle = zeros(nspan)
	xle = zeros(nspan)
	zle = zeros(nspan)	    
       	twist = zeros(nspan)

        kinem = KinemPar[]
        
	theta = zeros(ndiv)
        x = zeros(nspan,ndiv)

        dtheta = pi/(ndiv-1)
        for ib = 1:ndiv
            theta[ib] = real(ib-1.)*dtheta
	    for is = 1:nspan	      
  	        x[is,ib] = c[is]/2.*(1-cos(theta[ib]))
	    end
	end
        
	nspan = 0
        for i = 1:length(patchdata) - 1
            startpsi = (patchdata[i].y - patchdata[1].y)*pi/(patchdata[length(patchdata)].y - patchdata[1].y)
            lenpsi = (patchdata[i+1].y - patchdata[i].y)*pi/(patchdata[length(patchdata)].y - patchdata[1].y)
            dpsi =  lenpsi/patchdata[i].nspan
#            divspan = (patchdata[i+1].y - patchdata[i].y)/patchdata[i].nspan
            for j = 1:patchdata[i].nspan
                nspan += 1
                psi[nspan] = startpsi + (j-1)*dpsi
                #                yle[nspan] = patchdata[i].y + (i-1)*divspan
                yle[nspan] = -bref*cos(psi[nspan])/2.
                xle[nspan] = interp(patchdata[i].y, patchdata[i+1].y, patchdata[i].x, patchdata[i+1].x, yle[nspan])
                zle[nspan] = interp(patchdata[i].y, patchdata[i+1].y, patchdata[i].z, patchdata[i+1].z, yle[nspan])
                c[nspan] = interp(patchdata[i].y, patchdata[i+1].y, patchdata[i].c, patchdata[i+1].c, yle[nspan])
                twist[nspan] = interp(patchdata[i].y, patchdata[i+1].y, patchdata[i].twist, patchdata[i+1].twist, yle[nspan])
                if (patchdata[i].coord_file != "FlatPlate")
                    cam1, cam_slope1 = camber_calc(0.5*(1-cos(theta)), patchdata[i].coord_file)
                end
                if (patchdata[i+1].coord_file != "FlatPlate")
                    cam2, cam_slope2 = camber_calc(0.5*(1-cos(theta)), patchdata[i+1].coord_file)
                end
                for k = 1:ndiv
                    cam[nspan,k] = interp(patchdata[i].y, patchdata[i+1].y, cam1[k], cam2[k], yle[nspan])*c[nspan]
                    cam_slope[nspan,k] = interp(patchdata[i].y, patchdata[i+1].y, cam_slope1[k], cam_slope2[k], yle[nspan])*c[nspan]
                end
                lespcrit[nspan] = interp(patchdata[i].y, patchdata[i+1].y, patchdata[i].lc, patchdata[i+1].lc, yle[nspan])
                pvt[nspan] = interp(patchdata[i].y, patchdata[i+1].y, patchdata[i].pvt, patchdata[i+1].pvt, yle[nspan])
            end
        end
        i = length(patchdata)
        nspan += 1
        psi[nspan] = pi
        yle[nspan] = patchdata[i].y
        xle[nspan] = patchdata[i].x
        zle[nspan] = patchdata[i].z
        c[nspan] = patchdata[i].c
        twist[nspan] = patchdata[i].twist
        pvt[nspan] = patchdata[i].pvt
        lespcrit[nspan] = patchdata[i].lc
        if (patchdata[i].coord_file != "FlatPlate")
            cam[nspan,:], cam_slope[nspan,:] = camber_calc(x[nspan,:], patchdata[i].coord_file)
        end
        
        inkim = KinemPar(0, 0, 0, 0, 0, 0)
        for i = 1:nspan
            push!(kinem, inkim)
        end


        if (kindef.vartype == "Constant")
            i = 1
            if (typeof(kindef.alpha) == EldUpDef)
                kinem[i].alpha = kindef.alpha(0.)
                kinem[i].alphadot = ForwardDiff.derivative(kindef.alpha,0.)*uref/cref
            elseif (typeof(kindef.alpha) == EldRampReturnDef)
                kinem[i].alpha = kindef.alpha(0.)
                kinem[i].alphadot = ForwardDiff.derivative(kindef.alpha,0.)*uref/cref
            elseif (typeof(kindef.alpha) == ConstDef)
                kinem[i].alpha = kindef.alpha(0.)
                kinem[i].alphadot = 0.
            elseif (typeof(kindef.alpha) == SinDef)
                kinem[i].alpha = kindef.alpha(0.)
                kinem[i].alphadot = ForwardDiff.derivative(kindef.alpha,0.)*uref/cref
            elseif (typeof(kindef.alpha) == CosDef)
                kinem[i].alpha = kindef.alpha(0.)
                kinem[i].alphadot = ForwardDiff.derivative(kindef.alpha,0.)*uref/cref
            end
                
            if (typeof(kindef.h) == EldUpDef)
                kinem[i].h = kindef.h(0.)*cref
                kinem[i].hdot = ForwardDiff.derivative(kindef.h,0.)*uref
            elseif (typeof(kindef.h) == EldUpIntDef)
                kinem[i].h = kindef.h(0.)*cref
                kinem[i].hdot = ForwardDiff.derivative(kindef.h,0.)*uref
            elseif (typeof(kindef.h) == EldRampReturnDef)
                kinem[i].h = kindef.h(0.)*cref
                kinem[i].hdot = ForwardDiff.derivative(kindef.h,0.)*uref
            elseif (typeof(kindef.h) == ConstDef)
                kinem[i].h = kindef.h(0.)*cref
                kinem[i].hdot = 0.
            elseif (typeof(kindef.h) == SinDef)
                kinem[i].h = kindef.h(0.)*cref
                kinem[i].hdot = ForwardDiff.derivative(kindef.h,0.)*uref
            elseif (typeof(kindef.h) == CosDef)
                kinem[i].h = kindef.h(0.)*cref
                kinem[i].hdot = ForwardDiff.derivative(kindef.h,0.)*uref
            end
            
            if (typeof(kindef.u) == EldUpDef)
                kinem[i].u = kindef.u(0.)*uref
                kinem[i].udot = ForwardDiff.derivative(kindef.u,0.)*uref*uref/cref
            elseif (typeof(kindef.u) == EldRampReturnDef)
                kinem[i].u, kinem[i].udot = kindef.u(0.)
                kinem[i].u = kinem[i].u*uref
                kinem[i].udot = kinem[i].udot*uref*uref/cref
            elseif (typeof(kindef.u) == ConstDef)
                kinem[i].u = kindef.u(0.)*uref
                kinem[i].udot = 0.
            end
            
            for i = 2:nspan
                kinem[i] = kinem[1]
            end
        elseif (typeof(kindef.vartype) == "Read")
            # Provide the max amplitude as the input, and the read array goes from 0 to 1 as a function of yle
            i = nspan
            if (typeof(kindef.alpha) == EldUpDef)
                kinem[i].alpha = kindef.alpha(0.)
                kinem[i].alphadot = ForwardDiff.derivative(kindef.alpha,0.)*uref/cref
            elseif (typeof(kindef.alpha) == EldRampReturnDef)
                kinem[i].alpha = kindef.alpha(0.)
                kinem[i].alphadot = ForwardDiff.derivative(kindef.alpha,0.)*uref/cref
            elseif (typeof(kindef.alpha) == ConstDef)
                kinem[i].alpha = kindef.alpha(0.)
                kinem[i].alphadot = 0.
            elseif (typeof(kindef.alpha) == SinDef)
                kinem[i].alpha = kindef.alpha(0.)
                kinem[i].alphadot = ForwardDiff.derivative(kindef.alpha,0.)*uref/cref
            elseif (typeof(kindef.alpha) == CosDef)
                kinem[i].alpha = kindef.alpha(0.)
                kinem[i].alphadot = ForwardDiff.derivative(kindef.alpha,0.)*uref/cref
            end
            
            if (typeof(kindef.h) == EldUpDef)
                kinem[i].h = kindef.h(0.)*cref
                kinem[i].hdot = ForwardDiff.derivative(kindef.h,0.)*uref
            elseif (typeof(kindef.h) == EldUpIntDef)
                kinem[i].h = kindef.h(0.)*cref
                kinem[i].hdot = ForwardDiff.derivative(kindef.h,0.)*uref
            elseif (typeof(kindef.h) == EldRampReturnDef)
                kinem[i].h = kindef.h(0.)*cref
                kinem[i].hdot = ForwardDiff.derivative(kindef.h,0.)*uref
            elseif (typeof(kindef.h) == ConstDef)
                kinem[i].h = kindef.h(0.)*cref
                kinem[i].hdot = 0.
            elseif (typeof(kindef.h) == SinDef)
                kinem[i].h = kindef.h(0.)*cref
                kinem[i].hdot = ForwardDiff.derivative(kindef.h,0.)*uref
            elseif (typeof(kindef.h) == CosDef)
                kinem[i].h = kindef.h(0.)*cref
                kinem[i].hdot = ForwardDiff.derivative(kindef.h,0.)*uref
            end
            
            if (typeof(kindef.u) == EldUpDef)
                kinem[i].u = kindef.u(0.)*uref
                kinem[i].udot = ForwardDiff.derivative(kindef.u,0.)*uref*uref/cref
            elseif (typeof(kindef.u) == EldRampReturnDef)
                kinem[i].u, kinem[i].udot = kindef.u(0.)
                kinem[i].u = kinem[i].u*uref
                kinem[i].udot = kinem[i].udot*uref*uref/cref
            elseif (typeof(kindef.u) == ConstDef)
                kinem[i].u = kindef.u(0.)*uref
                kinem[i].udot = 0.
            end

addspl = Spline1D(kindef.add[:,1], kindef.add[:,2])

if (kindef.vary == 1)
    
    for i = 1:nspan
        kinem.h[i] = kinem.h[nspan]*evaluate(addspl,yle[i])
    end
elseif (kindef.vary == 2)
    addspl = Spline1D(kinem.add[:,1], kinem.add[:,2])
    for i = 1:nspan
        kinem.alpha[i] = kinem.alpha[nspan]*evaluate(addspl,yle[i])
    end 
end

end



for j = 1:nspan
    for i = 1:ndiv
        bnd_x[j,i] = -((c[j] - pvt[j]*c[j])+((pvt[j]*c[j] - x[j,i])*cos(kinem[j].alpha))) + (cam[j,i]*sin(kinem[j].alpha))
        bnd_z[j,i] = kinem[j].h + ((pvt[j]*c[j] - x[j,i])*sin(kinem[j].alpha))+(cam[j,i]*cos(kinem[j].alpha))
    end
end

uind = zeros(nspan,ndiv)
vind = zeros(nspan,ndiv)
wind = zeros(nspan,ndiv)
downwash = zeros(nspan,ndiv)
a0 = zeros(nspan)
a0dot = zeros(nspan)
aterm = zeros(nspan,naterm)
adot = zeros(nspan,3)
a0prev = zeros(nspan)
aprev = zeros(nspan,3)
bv = Array(ThreeDVort,nspan,ndiv-1)

for j = 1:nspan
    for i = 1:ndiv-1
        bv[j,i] = ThreeDVort([0; 0; 0;], [0; 0; 0;], 0.02*cref, 0, 0, 0)
    end
end
levflag = zeros(nspan)

new(cref, bref, sref, uref, patchdata, ndiv, nspan, naterm, nbterm, kindef, c, pvt, cam, cam_slope, theta, psi, x, xle, yle, zle, kinem, bnd_x, bnd_z, uind, vind, wind, downwash, a0, aterm, a0dot, adot, a0prev, aprev, bv,lespcrit,levflag)
end
end
 


# ---------------------------------------------------------------------------------------------
# BEGIN EldUpIntDef
immutable EldUpIntDef <: MotionDef
    amp :: Float64
    K :: Float64
    a :: Float64
end

immutable EldUpInttstartDef <: MotionDef
    amp :: Float64
    K :: Float64
    a :: Float64
    tstart :: Float64
end


function (eld::EldUpIntDef)(t)

    dt = 0.015
    nsteps = t/dt + 1
    nsteps = round(Int,nsteps)
    dt = t/(nsteps-1)
    sm = pi*pi*eld.K/(2*eld.amp*(1 - eld.a))
    t1 = 1.
    t2 = t1 + ((eld.amp)/(2*eld.K))


    prev_h = 0
    amp = 0
    for i = 1:nsteps
      tmpt = (i-1)*dt
      if (eld.amp == 0.)
      	 hdot = 0.
      else
         hdot = ((eld.K/sm)*log(cosh(sm*(tmpt - t1))/cosh(sm*(tmpt - t2))))+(eld.amp/2.)
      end
      amp = prev_h + hdot*dt
      prev_h = amp
    end
    if (nsteps == 1)
      amp = 0.
    end
    amp
end

function (eld::EldUpInttstartDef)(t)

    dt = 0.015
    nsteps = t/dt + 1
    nsteps = round(Int,nsteps)
    dt = t/(nsteps-1)
    sm = pi*pi*eld.K/(2*eld.amp*(1 - eld.a))
    t1 = eld.tstart
    t2 = t1 + ((eld.amp)/(2*eld.K))


    prev_h = 0
    amp = 0
    for i = 1:nsteps
      tmpt = (i-1)*dt
      if (eld.amp == 0.)
      	 hdot = 0.
      else
         hdot = ((eld.K/sm)*log(cosh(sm*(tmpt - t1))/cosh(sm*(tmpt - t2))))+(eld.amp/2.)
      end
      amp = prev_h + hdot*dt
      prev_h = amp
    end
    if (nsteps == 1)
      amp = 0.
    end
    amp
end

# END EldUpIntDef
# ---------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------
# BEGIN KelvinCondition
immutable KelvinCondition
    surf :: TwoDSurf
    field :: TwoDFlowField
end

immutable KelvinConditionwFlap
    surf :: TwoDSurfwFlap
    field :: TwoDFlowField
end

immutable KelvinCondition2DOF
    surf :: TwoDSurf_2DOF
    field :: TwoDFlowField
    kelv_enf :: Float64
end

immutable KelvinCondition2DFree
    surf :: TwoDFreeSurf
    field :: TwoDFlowField
    kelv_enf :: Float64
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



function (kelv::KelvinCondition2DOF)(tev_iter::Array{Float64})
    #Update the TEV strength
    nlev = length(kelv.field.lev)
    ntev = length(kelv.field.tev)
    kelv.field.tev[ntev].s = tev_iter[1]

    #Update incduced velocities on airfoil
    update_indbound(kelv.surf, kelv.field)

    #Calculate downwash
    update_downwash(kelv.surf)

    #Calculate first two fourier coefficients
    update_a0anda1(kelv.surf)

    val = kelv.surf.uref*kelv.surf.c*pi*(kelv.surf.a0[1] + kelv.surf.aterm[1]/2.)

    for iv = 1:ntev
        val = val + kelv.field.tev[iv].s
    end
    for iv = 1:nlev
        val = val + kelv.field.lev[iv].s
    end
    val = val + kelv.kelv_enf

    #Add kelv_enforced if necessary - merging will be better
    # val is the value of Gam_b + sum Gam_tev + Gam_lev which will equal zero
    # if the condition is satified

    return val
end

function (kelv::KelvinCondition2DFree)(tev_iter::Array{Float64})
    #Update the TEV strength
    nlev = length(kelv.field.lev)
    ntev = length(kelv.field.tev)
    kelv.field.tev[ntev].s = tev_iter[1]

    #Update incduced velocities on airfoil
    update_indbound(kelv.surf, kelv.field)

    #Calculate downwash
    update_downwash(kelv.surf)

    #Calculate first two fourier coefficients
    update_a0anda1(kelv.surf)

    val = kelv.surf.uref*kelv.surf.c*pi*(kelv.surf.a0[1] + kelv.surf.aterm[1]/2.)

    for iv = 1:ntev
        val = val + kelv.field.tev[iv].s
    end
    for iv = 1:nlev
        val = val + kelv.field.lev[iv].s
    end
    val = val + kelv.kelv_enf

    #Add kelv_enforced if necessary - merging will be better
    # val is the value of Gam_b + sum Gam_tev + Gam_lev which will equal zero
    # if the condition is satified

    return val
end

function (kelv::KelvinConditionwFlap)(tev_iter::Array{Float64})
    #Update the TEV strength
    nlev = length(kelv.field.lev)
    ntev = length(kelv.field.tev)
    kelv.field.tev[ntev].s = tev_iter[1]

    #Update incduced velocities on airfoil
    update_indbound(kelv.surf, kelv.field)

    #Calculate downwash
    update_downwash(kelv.surf)

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

# END KelvinCondition
# ---------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------
# BEGIN KelvinKutta
immutable KelvinKutta
    surf :: TwoDSurf
    field :: TwoDFlowField
end

immutable KelvinKutta2DOF
    surf :: TwoDSurf_2DOF
    field :: TwoDFlowField
    kelv_enf :: Float64
end

immutable KelvinKutta2DFree
    surf :: TwoDFreeSurf
    field :: TwoDFlowField
    kelv_enf :: Float64
end

immutable KelvinKuttawFlap
    surf :: TwoDSurfwFlap
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
    update_downwash(kelv.surf)

    #Calculate first two fourier coefficients
    update_a0anda1(kelv.surf)

    val[1] = kelv.surf.uref*kelv.surf.c*pi*(kelv.surf.a0[1] + kelv.surf.aterm[1]/2.)

    for iv = 1:ntev
        val[1] = val[1] + kelv.field.tev[iv].s
    end
    for iv = 1:nlev
        val[1] = val[1] + kelv.field.lev[iv].s
    end
    val[1] = val[1] + kelv.kelv_enf

    if (kelv.surf.a0[1] > 0)
        lesp_cond = kelv.surf.lespcrit[1]
    else
        lesp_cond = -kelv.surf.lespcrit[1]
    end
    val[2] = kelv.surf.a0[1]-lesp_cond

    #Add kelv_enforced if necessary - merging will be better
    return val
end

function (kelv::KelvinKutta2DFree)(v_iter::Array{Float64})
    val = zeros(2)

    #Update the TEV and LEV strengths
    nlev = length(kelv.field.lev)
    ntev = length(kelv.field.tev)
    kelv.field.tev[ntev].s = v_iter[1]
    kelv.field.lev[nlev].s = v_iter[2]

    #Update incduced velocities on airfoil
    update_indbound(kelv.surf, kelv.field)

    #Calculate downwash
    update_downwash(kelv.surf)

    #Calculate first two fourier coefficients
    update_a0anda1(kelv.surf)

    val[1] = kelv.surf.uref*kelv.surf.c*pi*(kelv.surf.a0[1] + kelv.surf.aterm[1]/2.)

    for iv = 1:ntev
        val[1] = val[1] + kelv.field.tev[iv].s
    end
    for iv = 1:nlev
        val[1] = val[1] + kelv.field.lev[iv].s
    end
    val[1] = val[1] + kelv.kelv_enf

    if (kelv.surf.a0[1] > 0)
        lesp_cond = kelv.surf.lespcrit[1]
    else
        lesp_cond = -kelv.surf.lespcrit[1]
    end
    val[2] = kelv.surf.a0[1]-lesp_cond

    #Add kelv_enforced if necessary - merging will be better
    return val
end

function (kelv::KelvinKuttawFlap)(v_iter::Array{Float64})
    val = zeros(2)

    #Update the TEV and LEV strengths
    nlev = length(kelv.field.lev)
    ntev = length(kelv.field.tev)
    kelv.field.tev[ntev].s = v_iter[1]
    kelv.field.lev[nlev].s = v_iter[2]

    #Update incduced velocities on airfoil
    update_indbound(kelv.surf, kelv.field)

    #Calculate downwash
    update_downwash(kelv.surf)

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

# END KelvinKutta
# ---------------------------------------------------------------------------------------------
