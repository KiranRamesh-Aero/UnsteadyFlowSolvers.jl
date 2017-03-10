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
    v :: Vector{Float64}
    ds :: Vector{Float64}
end
# ---------------------------------------------------------------------------------------------
#Panel used in the lumped vorte method
type TwoDLVPanel
    xv :: Float64
    zv :: Float64
    xv_I :: Float64
    zv_I :: Float64
    s :: Float64
    l :: Float64
    xc :: Float64
    zc :: Float64
    xc_I :: Float64
    zc_I :: Float64
    nx :: Float64
    nz :: Float64
    tx :: Float64
    tz :: Float64
    uind :: Float64
    wind :: Float64
    uind_I :: Float64
    wind_I :: Float64
end


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

immutable TwoDFlowFieldMultSurf
    nsurf :: Int
    velX :: MotionDef
    velZ :: MotionDef
    u :: Vector{Float64}
    w :: Vector{Float64}
    tev :: Vector{Vector{TwoDVort}}
    lev :: Vector{Vector{TwoDVort}}
    extv :: Vector{TwoDVort}
    function TwoDFlowFieldMultSurf(nsurf, velX = ConstDef(0.), velZ = ConstDef(0.))
        u = [0;]
        w = [0;]
        tev = TwoDVort[]
        lev = TwoDVort[]
        extv = TwoDVort[]
        new(nsurf, velX, velZ, u, w, tev, lev, extv)
    end
end


#Definition of flowfield in 3D
immutable ThreeDFieldSimple
    f2d :: Vector{TwoDFlowField}
    function ThreeDFieldSimple()
        f2d = TwoDFlowField[]
    new(f2d)
    end
end


# immutable ThreeDFlowField
#     velX :: MotionDef
#     velY :: MotionDef
#     velZ :: MotionDef
#     u :: Vector{Float64}
#     v :: Vector{Float64}
#     w :: Vector{Float64}
#     tev :: Vector{ThreeDVort}
#     lev :: Vector{ThreeDVort}
#     extv :: Vector{ThreeDVort}
#     function ThreeDFlowField(velX = ConstDef(0.), velY = ConstDef(0.), velZ = ConstDef(0.))
#         u = [0;]
#         v = [0;]
#         w = [0;]
#         tev = ThreeDVort[]
#         lev = ThreeDVort[]
#         extv = ThreeDVort[]
#         new(velX, velY, velZ, u, v, w, tev, lev, extv)
#     end
# end


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

#TwoD surface represented by Lumped Vortex Panels
immutable TwoDSurfLV
    c :: Float64
    uref :: Float64
    coord_file :: String
    pvt :: Float64
    npanel :: Int8
    distype :: String
    kindef :: KinemDef
    cam :: Vector{Float64}
    cam_slope :: Vector{Float64}
    x :: Vector{Float64}
    kinem :: KinemPar
    X0 :: Vector{Float64}
    tevloc :: Float64
    a0 :: Vector{Float64}
    a0dot :: Vector{Float64}
    a0prev :: Vector{Float64}
    lv :: Vector{TwoDLVPanel}
    lespcrit :: Vector{Float64}
    levflag :: Vector{Int8}
    IC :: Array{Float64}
    rhs :: Vector{Float64}
    gamma :: Vector{Float64}
    gamma_prev :: Vector{Float64}
    tgl :: Array{Float64}
    tlg :: Array{Float64}
    
    function TwoDSurfLV(coord_file, pvt, kindef;lespcrit=zeros(1), c=1., uref=1., npanel=70, distype="Cosine", tevloc=0.5)

        x = zeros(npanel+1)
        cam = zeros(npanel+1)
        cam_slope = zeros(npanel+1)
        # bnd_x = zeros(ndiv)
        # bnd_z = zeros(ndiv)
        kinem = KinemPar(0, 0, 0, 0, 0, 0)

        if (distype == "linear")
            dc = c/npanel
            for i = 1:npanel+1
                x[i] = real(i-1.)*dc
            end

            if (coord_file != "FlatPlate")
                cam, cam_slope = camber_calc(x, coord_file)
            end
            
            camspl = Spline1D(x, cam)
            cam_slopespl = Spline1D(x, cam_slope)
            
            lv = TwoDLVPanel[]
            for i = 1:npanel
                push!(lv, TwoDLVPanel(0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.) )
                #Length of panel
                lv[i].l = sqrt((x[i+1] - x[i])^2 + (cam[i+1] - cam[i])^2)
                # Vortex Points (1/4c) (body frame)
                lv[i].xv = x[i] + dc/4; #linspace(0,c-dc,N) + dc/4;
                lv[i].zv = evaluate(camspl, lv[i].xv)
                # collocation Points (3/4c) (body frame)
                lv[i].xc = x[i] + 3*dc/4; #linspace(0,c-dc,N) + 3*dc/4;
                lv[i].zc = evaluate(camspl, lv[i].xc)
                deta = evaluate(cam_slopespl, lv[i].xc)
                den = sqrt(deta^2 + 1)
                # Normal and Tangent Vectors at Collocation Points (ref to body frame)
                lv[i].nx = -deta/ den
                lv[i].nz = 1/den
                lv[i].tx = lv[i].nz
                lv[i].tz = -lv[i].nx
            end

        elseif (distype == "Cosine")
            theta = zeros(npanel+1)
            dtheta = pi/(npanel)
            for i = 1:npanel+1
                theta[i] = real(i-1.)*dtheta
                x[i] = c/2.*(1-cos(theta[i]))
            end
            
            if (coord_file != "FlatPlate")
                cam, cam_slope = camber_calc(x, coord_file)
            end
            
            camspl = Spline1D(x, cam)
            cam_slopespl = Spline1D(x, cam_slope)

            lv = TwoDLVPanel[]
            for i = 1:npanel
                push!(lv, TwoDLVPanel(0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.) )
                #Length of panel
                lv[i].l = sqrt((x[i+1] - x[i])^2 + (cam[i+1] - cam[i])^2)
                # Vortex Points (1/4c) (body frame)
                dc = x[i+1] - x[i]
                lv[i].xv = x[i] + dc/4; #linspace(0,c-dc,N) + dc/4;
                lv[i].zv = evaluate(camspl, lv[i].xv) 
                # collocation Points (3/4c) (body frame)
                lv[i].xc = x[i] + 3*dc/4; #linspace(0,c-dc,N) + 3*dc/4;
                lv[i].zc = evaluate(camspl, lv[i].xc)
                deta = evaluate(cam_slopespl, lv[i].xc)
                den = sqrt(deta^2 + 1.)
                # Normal and Tangent Vectors at Collocation Points (ref to body frame)
                lv[i].nx = -deta/ den
                lv[i].nz = 1./den
                lv[i].tx = lv[i].nz
                lv[i].tz = -lv[i].nx
            end
        else
            error("Invalid panel distribution type")
        end
        
        # Self-Induced Velocity is constant, as long as the airfoil doesn't change geometry
        #Influence coefficient is a constant
        IC = zeros(npanel+1,npanel+1)
        uw = zeros(2)
        for i = 1:npanel
            for j = 1:npanel+1
                if j > npanel
                    #For varying time steps this column be updated at every time step
                    dtref = 0.015*c/uref
                    xloc = x[npanel+1] + tevloc*uref*dtref
                    zloc = cam[npanel+1]
                    dummytev = TwoDVort(xloc, zloc, 1.0, 0.0, 0.0, 0.0)
                    ui, wi = ind_vel([dummytev;], lv[i].xc, lv[i].zc)
                    uw[1] = ui[1]
                    uw[2] = wi[1]
                    IC[i,j] = dot(uw,[lv[i].nx; lv[i].nz])
                else
                    dummyj = TwoDVort(lv[j].xv, lv[j].zv, 1.0, 0.0, 0.0, 0.0) 
                    ui, wi =ind_vel([dummyj;], lv[i].xc, lv[i].zc)
                    uw[1] = ui[1]
                    uw[2] = wi[1]
                    IC[i,j] = dot(uw, [lv[i].nx; lv[i].nz])
                end
            end
        end
        i = npanel+1
        for j = 1:npanel+1
            IC[i,j]= 1.
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

tlg = [cos(kinem.alpha) -sin(kinem.alpha); sin(kinem.alpha) cos(kinem.alpha)]
tgl = transpose(tlg)

X0 = zeros(2) 
X0[2] = kinem.h

rhs = zeros(npanel+1)
a0 = [0.;]
a0dot = [0.;]
a0prev = [0.;]
levflag = [0]
gamma = [0.]
gamma_prev = [0.]

new(c, uref, coord_file, pvt, npanel, distype, kindef, cam, cam_slope, x, kinem, X0, tevloc, a0, a0dot, a0prev, lv, lespcrit,levflag, IC, rhs, gamma, gamma_prev, tgl, tlg)
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

# immutable ThreeDSurf
#     cref :: Float64
#     bref :: Float64
#     sref :: Float64
#     uref :: Float64
#     patchdata :: Vector{patch}
#     ndiv :: Int8
#     nspan :: Int8
#     naterm :: Int8
#     kindef :: KinemDef3D
#     c :: Vector{Float64}
#     pvt :: Vector{Float64}
#     cam :: Array{Float64,2}
#     cam_slope :: Array{Float64,2}
#     theta :: Vector{Float64}
#     psi :: Vector{Float64}
#     x :: Array{Float64,2}
#     xle :: Vector{Float64}
#     yle :: Vector{Float64}
#     zle :: Vector{Float64}
#     kinem :: Vector{KinemPar}
#     bnd_x :: Array{Float64,2}
#     bnd_z :: Array{Float64,2}
#     uind :: Array{Float64,2}
#     vind :: Array{Float64,2}
#     wind :: Array{Float64,2}
#     downwash :: Array{Float64,2}
#     a0 :: Vector{Float64}
#     aterm :: Array{Float64,2}
#     a0dot :: Vector{Float64}
#     adot :: Array{Float64,2}
#     a0prev :: Vector{Float64}
#     aprev :: Array{Float64,2}
#     bv :: Array{ThreeDVort,2}
#     lespcrit :: Vector{Float64}
#     levflag :: Vector{Int8}

#     function ThreeDSurf(cref, bref, sref, patchdata, kindef, uref=1., ndiv=70, naterm=35)
#         nspan = 0
#         for i = 1:length(patchdata)-1
#             nspan += patchdata[i].nspan
#         end

#         lespcrit = zeros(nspan)

#         psi = zeros(nspan)
#         cam = zeros(nspan,ndiv)
#         cam_slope = zeros(nspan,ndiv)
#         cam1 = zeros(ndiv)
#         cam_slope1 = zeros(ndiv)
#         cam2 = zeros(ndiv)
#         cam_slope2 = zeros(ndiv)
#         bnd_x = zeros(nspan,ndiv)
#         bnd_z = zeros(nspan,ndiv)
#         c = zeros(nspan)
#         pvt = zeros(nspan)
# 	yle = zeros(nspan)
# 	xle = zeros(nspan)
# 	zle = zeros(nspan)
#        	twist = zeros(nspan)

#         kinem = KinemPar[]

# 	theta = zeros(ndiv)
#         x = zeros(nspan,ndiv)

#         dtheta = pi/(ndiv-1)
#         for ib = 1:ndiv
#             theta[ib] = real(ib-1.)*dtheta
# 	    for is = 1:nspan
#   	        x[is,ib] = c[is]/2.*(1-cos(theta[ib]))
# 	    end
# 	end

# 	nspan = 0
#         for i = 1:length(patchdata) - 1
#             if i != 1
#                 nspan += 1
#                 psi[nspan] = (patchdata[i].y - patchdata[1].y)*pi/(patchdata[length(patchdata)].y - patchdata[1].y)
#                 yle[nspan] = patchdata[i].y
#                 xle[nspan] = patchdata[i].x
#                 zle[nspan] = patchdata[i].z
#                 c[nspan] = patchdata[i].c
#                 twist[nspan] = patchdata[i].twist
#                 if (patchdata[i].coord_file != "FlatPlate")
#                     cam[nspan,:], cam_slope[nspan,:] = camber_calc(0.5*(1-cos(theta)), patchdata[i].coord_file)
#                 end
#                 lespcrit[nspan] = interp(patchdata[i].y, patchdata[i+1].y, patchdata[i].lc, patchdata[i+1].lc, yle[nspan])
#                 pvt[nspan] = interp(patchdata[i].y, patchdata[i+1].y, patchdata[i].pvt, patchdata[i+1].pvt, yle[nspan])
#             end
#             startpsi = (patchdata[i].y - patchdata[1].y)*pi/(patchdata[length(patchdata)].y - patchdata[1].y)
#             lenpsi = (patchdata[i+1].y - patchdata[i].y)*pi/(patchdata[length(patchdata)].y - patchdata[1].y)
#             dpsi =  lenpsi/(patchdata[i].nspan + 1.)
# #            divspan = (patchdata[i+1].y - patchdata[i].y)/patchdata[i].nspan
#             for j = 1:patchdata[i].nspan
#                 nspan += 1
#                 psi[nspan] = startpsi + real(j)*dpsi
#                 #                yle[nspan] = patchdata[i].y + (i-1)*divspan
#                 yle[nspan] = -bref*cos(psi[nspan])/2.
#                 xle[nspan] = interp(patchdata[i].y, patchdata[i+1].y, patchdata[i].x, patchdata[i+1].x, yle[nspan])
#                 zle[nspan] = interp(patchdata[i].y, patchdata[i+1].y, patchdata[i].z, patchdata[i+1].z, yle[nspan])
#                 c[nspan] = interp(patchdata[i].y, patchdata[i+1].y, patchdata[i].c, patchdata[i+1].c, yle[nspan])
#                 twist[nspan] = interp(patchdata[i].y, patchdata[i+1].y, patchdata[i].twist, patchdata[i+1].twist, yle[nspan])
#                 if (patchdata[i].coord_file != "FlatPlate")
#                     cam1, cam_slope1 = camber_calc(0.5*(1-cos(theta)), patchdata[i].coord_file)
#                 end
#                 if (patchdata[i+1].coord_file != "FlatPlate")
#                     cam2, cam_slope2 = camber_calc(0.5*(1-cos(theta)), patchdata[i+1].coord_file)
#                 end
#                 for k = 1:ndiv
#                     cam[nspan,k] = interp(patchdata[i].y, patchdata[i+1].y, cam1[k], cam2[k], yle[nspan])*c[nspan]
#                     cam_slope[nspan,k] = interp(patchdata[i].y, patchdata[i+1].y, cam_slope1[k], cam_slope2[k], yle[nspan])*c[nspan]
#                 end
#                 lespcrit[nspan] = interp(patchdata[i].y, patchdata[i+1].y, patchdata[i].lc, patchdata[i+1].lc, yle[nspan])
#                 pvt[nspan] = interp(patchdata[i].y, patchdata[i+1].y, patchdata[i].pvt, patchdata[i+1].pvt, yle[nspan])
#             end
#         end

#         inkim = KinemPar(0, 0, 0, 0, 0, 0)
#         for i = 1:nspan
#             push!(kinem, inkim)
#         end


#         if (kindef.vartype == "Constant")
                  #             i = 1
#             if (typeof(kindef.alpha) == EldUpDef)
#                 kinem[i].alpha = kindef.alpha(0.)
#                 kinem[i].alphadot = ForwardDiff.derivative(kindef.alpha,0.)*uref/cref
#             elseif (typeof(kindef.alpha) == EldRampReturnDef)
#                 kinem[i].alpha = kindef.alpha(0.)
#                 kinem[i].alphadot = ForwardDiff.derivative(kindef.alpha,0.)*uref/cref
#             elseif (typeof(kindef.alpha) == ConstDef)
#                 kinem[i].alpha = kindef.alpha(0.)
#                 kinem[i].alphadot = 0.
#             elseif (typeof(kindef.alpha) == SinDef)
#                 kinem[i].alpha = kindef.alpha(0.)
#                 kinem[i].alphadot = ForwardDiff.derivative(kindef.alpha,0.)*uref/cref
#             elseif (typeof(kindef.alpha) == CosDef)
#                 kinem[i].alpha = kindef.alpha(0.)
#                 kinem[i].alphadot = ForwardDiff.derivative(kindef.alpha,0.)*uref/cref
#             end

#             if (typeof(kindef.h) == EldUpDef)
#                 kinem[i].h = kindef.h(0.)*cref
#                 kinem[i].hdot = ForwardDiff.derivative(kindef.h,0.)*uref
#             elseif (typeof(kindef.h) == EldUpIntDef)
#                 kinem[i].h = kindef.h(0.)*cref
#                 kinem[i].hdot = ForwardDiff.derivative(kindef.h,0.)*uref
#             elseif (typeof(kindef.h) == EldRampReturnDef)
#                 kinem[i].h = kindef.h(0.)*cref
#                 kinem[i].hdot = ForwardDiff.derivative(kindef.h,0.)*uref
#             elseif (typeof(kindef.h) == ConstDef)
#                 kinem[i].h = kindef.h(0.)*cref
#                 kinem[i].hdot = 0.
#             elseif (typeof(kindef.h) == SinDef)
#                 kinem[i].h = kindef.h(0.)*cref
#                 kinem[i].hdot = ForwardDiff.derivative(kindef.h,0.)*uref
#             elseif (typeof(kindef.h) == CosDef)
#                 kinem[i].h = kindef.h(0.)*cref
#                 kinem[i].hdot = ForwardDiff.derivative(kindef.h,0.)*uref
#             end

#             if (typeof(kindef.u) == EldUpDef)
#                 kinem[i].u = kindef.u(0.)*uref
#                 kinem[i].udot = ForwardDiff.derivative(kindef.u,0.)*uref*uref/cref
#             elseif (typeof(kindef.u) == EldRampReturnDef)
#                 kinem[i].u, kinem[i].udot = kindef.u(0.)
#                 kinem[i].u = kinem[i].u*uref
#                 kinem[i].udot = kinem[i].udot*uref*uref/cref
#             elseif (typeof(kindef.u) == ConstDef)
#                 kinem[i].u = kindef.u(0.)*uref
#                 kinem[i].udot = 0.
#             end

#             for i = 2:nspan
#                 kinem[i] = kinem[1]
#             end
#         elseif (typeof(kindef.vartype) == "Read")
#             # Provide the max amplitude as the input, and the read array goes from 0 to 1 as a function of yle
#             i = nspan
#             if (typeof(kindef.alpha) == EldUpDef)
#                 kinem[i].alpha = kindef.alpha(0.)
#                 kinem[i].alphadot = ForwardDiff.derivative(kindef.alpha,0.)*uref/cref
#             elseif (typeof(kindef.alpha) == EldRampReturnDef)
#                 kinem[i].alpha = kindef.alpha(0.)
#                 kinem[i].alphadot = ForwardDiff.derivative(kindef.alpha,0.)*uref/cref
#             elseif (typeof(kindef.alpha) == ConstDef)
#                 kinem[i].alpha = kindef.alpha(0.)
#                 kinem[i].alphadot = 0.
#             elseif (typeof(kindef.alpha) == SinDef)
#                 kinem[i].alpha = kindef.alpha(0.)
#                 kinem[i].alphadot = ForwardDiff.derivative(kindef.alpha,0.)*uref/cref
#             elseif (typeof(kindef.alpha) == CosDef)
#                 kinem[i].alpha = kindef.alpha(0.)
#                 kinem[i].alphadot = ForwardDiff.derivative(kindef.alpha,0.)*uref/cref
#             end

#             if (typeof(kindef.h) == EldUpDef)
#                 kinem[i].h = kindef.h(0.)*cref
#                 kinem[i].hdot = ForwardDiff.derivative(kindef.h,0.)*uref
#             elseif (typeof(kindef.h) == EldUpIntDef)
#                 kinem[i].h = kindef.h(0.)*cref
#                 kinem[i].hdot = ForwardDiff.derivative(kindef.h,0.)*uref
#             elseif (typeof(kindef.h) == EldRampReturnDef)
#                 kinem[i].h = kindef.h(0.)*cref
#                 kinem[i].hdot = ForwardDiff.derivative(kindef.h,0.)*uref
#             elseif (typeof(kindef.h) == ConstDef)
#                 kinem[i].h = kindef.h(0.)*cref
#                 kinem[i].hdot = 0.
#             elseif (typeof(kindef.h) == SinDef)
#                 kinem[i].h = kindef.h(0.)*cref
#                 kinem[i].hdot = ForwardDiff.derivative(kindef.h,0.)*uref
#             elseif (typeof(kindef.h) == CosDef)
#                 kinem[i].h = kindef.h(0.)*cref
#                 kinem[i].hdot = ForwardDiff.derivative(kindef.h,0.)*uref
#             end

#             if (typeof(kindef.u) == EldUpDef)
#                 kinem[i].u = kindef.u(0.)*uref
#                 kinem[i].udot = ForwardDiff.derivative(kindef.u,0.)*uref*uref/cref
#             elseif (typeof(kindef.u) == EldRampReturnDef)
#                 kinem[i].u, kinem[i].udot = kindef.u(0.)
#                 kinem[i].u = kinem[i].u*uref
#                 kinem[i].udot = kinem[i].udot*uref*uref/cref
#             elseif (typeof(kindef.u) == ConstDef)
#                 kinem[i].u = kindef.u(0.)*uref
#                 kinem[i].udot = 0.
#             end

# addspl = Spline1D(kindef.add[:,1], kindef.add[:,2])

# if (kindef.vary == 1)

#     for i = 1:nspan
#         kinem.h[i] = kinem.h[nspan]*evaluate(addspl,yle[i])
#     end
# elseif (kindef.vary == 2)
#     addspl = Spline1D(kinem.add[:,1], kinem.add[:,2])
#     for i = 1:nspan
#         kinem.alpha[i] = kinem.alpha[nspan]*evaluate(addspl,yle[i])
#     end
# end

# end



# for j = 1:nspan
#     for i = 1:ndiv
#         bnd_x[j,i] = -((c[j] - pvt[j]*c[j])+((pvt[j]*c[j] - x[j,i])*cos(kinem[j].alpha))) + (cam[j,i]*sin(kinem[j].alpha))
#         bnd_z[j,i] = kinem[j].h + ((pvt[j]*c[j] - x[j,i])*sin(kinem[j].alpha))+(cam[j,i]*cos(kinem[j].alpha))
#     end
# end

# uind = zeros(nspan,ndiv)
# vind = zeros(nspan,ndiv)
# wind = zeros(nspan,ndiv)
# downwash = zeros(nspan,ndiv)
# a0 = zeros(nspan)
# a0dot = zeros(nspan)
# aterm = zeros(nspan,naterm)
# adot = zeros(nspan,3)
# a0prev = zeros(nspan)
# aprev = zeros(nspan,3)
# bv = Array(ThreeDVort,nspan,ndiv-1)

# for j = 1:nspan
#     for i = 1:ndiv-1
#         bv[j,i] = ThreeDVort([0; 0; 0;], [0; 0; 0;], 0.02*cref, [0; 0; 0;], [0; 0; 0;])
#     end
# end
# levflag = zeros(nspan)

# new(cref, bref, sref, uref, patchdata, ndiv, nspan, naterm, nbterm, kindef, c, pvt, cam, cam_slope, theta, psi, x, xle, yle, zle, kinem, bnd_x, bnd_z, uind, vind, wind, downwash, a0, aterm, a0dot, adot, a0prev, aprev, bv,lespcrit,levflag)
# end
# end


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
    kindef :: KinemDef
    psi :: Vector{Float64}
    yle :: Vector{Float64}
    s2d :: Vector{TwoDSurf}
    a03d :: Vector{Float64}
    bc :: Vector{Float64}
    nshed :: Vector{Float64}

    function ThreeDSurfSimple(AR, kindef, coord_file, pvt, lespcrit = [10.;], nspan = 20, cref = 1., uref=1., ndiv=70, naterm=35)

        bref = AR*cref

        psi = zeros(nspan)
        yle = zeros(nspan)

        s2d = TwoDSurf[]

	for i = 1:nspan
            psi[i] = real(i)*(pi/2)/nspan
            yle[i] = -bref*cos(psi[i])/2.
        end

        for i = 1:nspan
            push!(s2d, TwoDSurf(coord_file, pvt,  kindef, lespcrit, cref, uref, ndiv, naterm))
        end

        a03d = zeros(nspan)
        bc = zeros(nspan)
        nshed = [0.;]

new(cref, AR, uref, pvt, lespcrit, coord_file,  ndiv, nspan, naterm, kindef, psi, yle, s2d, a03d, bc, nshed)
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

immutable KelvinConditionMultSurf
    surf :: Vector{TwoDSurf}
    field :: TwoDFlowFieldMultSurf
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

immutable KelvinConditionLLTldvm
    surf :: ThreeDSurfSimple
    field :: ThreeDFieldSimple
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

function (kelv::KelvinConditionMultSurf)(tev_iter::Array{Float64})

    val = zeros(kelv.field.nsurf)
    ntev = length(kelv.field.tev)
    nlev = length(kelv.field.lev)
    
    #Update the TEV strength
    for i = 1:kelv.field.nsurf
        kelv.field.tev[ntev][i].s = tev_iter[i]
    end
    
    
    #Update induced velocities on airfoil
    update_indbound(kelv.surf, kelv.field)
    
    #Calculate downwash
    update_downwash(kelv.surf, [kelv.field.u[1],kelv.field.w[1]])
        
    #Calculate first two fourier coefficients
    update_a0anda1(kelv.surf)
    
    for i = 1:kelv.field.nsurf
        val[i] = kelv.surf[i].uref*kelv.surf[i].c*pi*(kelv.surf[i].a0[1] + kelv.surf[i].aterm[1]/2.)

        for iv = 1:ntev
            val[i] = val[i] + kelv.field.tev[iv][i].s
        end
        for iv = 1:nlev
            val[i] = val[i] + kelv.field.lev[iv][i].s
        end
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

function (kelv::KelvinConditionLLTldvm)(tev_iter::Array{Float64})
    val = zeros(2*kelv.surf.nspan)
    lhs = zeros(kelv.surf.nspan,kelv.surf.nspan)
    rhs = zeros(kelv.surf.nspan)

    #Assume symmetry condition for now
    for i = 1:kelv.surf.nspan
        nlev = length(kelv.field.f2d[i].lev)
        ntev = length(kelv.field.f2d[i].tev)
        kelv.field.f2d[i].tev[ntev].s = tev_iter[i]

        #Update incduced velocities on airfoil
        update_indbound(kelv.surf.s2d[i], kelv.field.f2d[i])

        #Calculate downwash
        update_downwash(kelv.surf.s2d[i], [kelv.field.f2d[i].u[1], kelv.field.f2d[i].w[1]])

        #Calculate first two fourier coefficients
        update_a0anda1(kelv.surf.s2d[i])

        kelv.surf.bc[i] = kelv.surf.s2d[i].a0[1] + 0.5*kelv.surf.s2d[i].aterm[1]

        kelv.surf.a03d[i] = 0

        for n = 1:kelv.surf.nspan
            nn = 2*n - 1
            kelv.surf.a03d[i] = kelv.surf.a03d[i] - real(nn)*tev_iter[n+kelv.surf.nspan]*sin(nn*kelv.surf.psi[i])/sin(kelv.surf.psi[i])
        end

        val[i] = kelv.surf.s2d[i].uref*kelv.surf.s2d[i].c*pi*(kelv.surf.bc[i] + kelv.surf.a03d[i])

        for iv = 1:ntev
            val[i] = val[i] + kelv.field.f2d[i].tev[iv].s
        end
        for iv = 1:nlev
            val[i] = val[i] + kelv.field.f2d[i].lev[iv].s
        end
    end

    for i = 1:kelv.surf.nspan
        for n = 1:kelv.surf.nspan
            nn = 2*n - 1
            lhs[i,n] = sin(nn*kelv.surf.psi[i])*(sin(kelv.surf.psi[i]) + (nn*pi/(2*kelv.surf.AR)))
        end
        rhs[i] = pi*sin(kelv.surf.psi[i])*kelv.surf.bc[i]/(2*kelv.surf.AR)
    end

    val[kelv.surf.nspan+1:2*kelv.surf.nspan] = lhs*tev_iter[kelv.surf.nspan+1:2*kelv.surf.nspan] - rhs

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

immutable KelvinKuttaLLTldvm
    surf :: ThreeDSurfSimple
    field :: Vector{TwoDFlowField}
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

function (kelv::KelvinKuttaLLTldvm)(tev_iter::Array{Float64})
    val = zeros(2*kelv.surf3d.nspan + kelv.nshed)
    bc = zeros(kelv.surf3d.nspan)
    AR = kelv.surf3d.bref/kelv.surf3d.cref
    lhs = zeros(kelv.surf3d.nspan,kelv.surf3d.nspan)
    rhs = zeros(kelv.surf3d.nspan)
    bcoeff = zeros(kelv.surf3d.nspan)
    a03d = zeros(kelv.surf3d.nspan)

    cntr = kelv.surf3d.nspan + 1

    #Assume symmetry condition for now
    for i = 1:kelv.surf3d.nspan
        nlev = length(kelv.field[i].lev)
        ntev = length(kelv.field[i].tev)
        kelv.field[i].tev[ntev].s = tev_iter[i]
        if kelv.surf[i].levflag[1] == 1
            kelv.field[i].lev[nlev].s = tev_iter[cntr]
            cntr += 1
        end

        #Update incduced velocities on airfoil
        update_indbound(kelv.surf[i], kelv.field[i])

        #Calculate downwash
        update_downwash(kelv.surf[i], [kelv.field[i].u[1],kelv.field[i].w[1]])

        #Calculate first two fourier coefficients
        update_a0anda1(kelv.surf[i])

        bc[i] = kelv.surf[i].a0[1] + 0.5*kelv.surf[i].aterm[1]

        for i = 1:kelv.surf3d.nspan
        a03d[i] = 0
            for n = 1:kelv.surf3d.nspan
                nn = 2*n - 1
                a03d[i] = a03d[i] - real(nn)*tev_iter[n+kelv.surf3d.nspan+kelv.nshed]*sin(nn*kelv.surf3d.psi[i])/sin(kelv.surf3d.psi[i])
            end
        end

        val[i] = kelv.surf[i].uref*kelv.surf[i].c*pi*(bc[i] + a03d[i])

        for iv = 1:ntev
            val[i] = val[i] + kelv.field[i].tev[iv].s
        end
        for iv = 1:nlev
            val[i] = val[i] + kelv.field[i].lev[iv].s
        end

        if kelv.surf[i].levflag[1] == 1
            if kelv.surf[i].a0[1] > 0
                lesp_cond = kelv.surf[i].lespcrit[1]
            else
                lesp_cond = -kelv.surf[i].lespcrit[1]
            end
            val[cntr-1] = kelv.surf[i].a0[1] + a03d[i] - lesp_cond
        end
    end

    for i = 1:kelv.surf3d.nspan
        for n = 1:kelv.surf3d.nspan
            nn = 2*n - 1
            lhs[i,n] = sin(nn*kelv.surf3d.psi[i])*(sin(kelv.surf3d.psi[i]) + (nn*pi/(2*AR)))
        end
        rhs[i] = pi*sin(kelv.surf3d.psi[i])*bc[i]/(2*AR)
    end

    val[kelv.surf3d.nspan+kelv.nshed+1:2*kelv.surf3d.nspan+kelv.nshed] = lhs*tev_iter[kelv.surf3d.nspan+kelv.nshed+1:2*kelv.surf3d.nspan+kelv.nshed] - rhs

    return val
end

# END KelvinKutta
# ---------------------------------------------------------------------------------------------
#Types for 3D vortex ring method
type ThreeDSurfVRingGrid
    xv :: Float64
    yv :: Float64
    zv :: Float64
    xv_I :: Float64
    yv_I :: Float64
    zv_I :: Float64
end

type  ThreeDSurfVRingPanel
    s :: Float64
    xc :: Float64
    yc :: Float64
    zc :: Float64
    xc_I :: Float64
    yc_I :: Float64
    zc_I :: Float64
    nx :: Float64
    ny :: Float64
    nz :: Float64
end



# ---------------------------------------------------------------------------------------------
immutable ThreeDSurfVR
    cref :: Float64
    AR :: Float64
    uref :: Float64 
    pvt :: Float64
    lespcrit :: Vector{Float64}
    coord_file :: String
    nspan :: Int
    nchord :: Int
    kinem :: KinemPar
    kindef :: KinemDef
    a0 :: Vector{Float64}
    vr_g :: Vector{ThreeDSurfVRingGrid}
    vr_p :: Vector{ThreeDSurfVRingPanel}
    distype :: String
    IC :: Array{Float64}
    IC_W :: Array{Float64}
    rhs:: Vector{Float64}
    tgl :: Array{Float64}
    tlg :: Array{Float64}
    X0 :: Vector{Float64}
    tevloc :: Float64
    x :: Array{Float64}
    y :: Array{Float64}
    z :: Array{Float64}
    
    function ThreeDSurfVR(AR, kindef, coord_file, pvt, lespcrit = [10.;]; nspan = 20, nchord = 10, cref = 1., uref=1., distype = "Linear", tevloc = 0.2)
        
        bref = AR*cref
        sref = cref*bref

        kinem = KinemPar(0, 0, 0, 0, 0, 0)

        #Rectangular wing is assumed - to be updated
        if distype == "Linear" || distype == "linear"
            x = zeros(nchord+1, nspan+1)
            y = zeros(nchord+1, nspan+1)
            z = zeros(nchord+1, nspan+1)

            # Panel edges
            for i = 1:nchord+1
                for j = 1:nspan+1
                    x[i,j] = real(i - 1.)*cref/nchord
                    y[i,j] = 0.5*real(j - 1.)*bref/nspan
                    z[i,j] = 0.
                end
            end
            
            #Votex rings definitions
            vr_g = ThreeDSurfVRingGrid[]
            vr_p = ThreeDSurfVRingPanel[]            
            
            sp = sref/(nspan*nchord) # Panel area
            bp = bref/nspan # Panels span
            cp = cref/nchord # Panels chord

            #Set panel properties of the rings
            for i = 1:nchord
                for j = 1:nspan
                    
                    # Collocation Points
                    xc = 3.*cp/4. + x[i,j]
                    yc = bp/2. + y[i,j]
                    zc = z[i,j]
                    
                    # Normal and Tangent Vectors at Collocation Points
                    a_k = [x[i+1,j+1] - x[i,j]; y[i+1,j+1] - y[i,j]; z[i+1,j+1] - z[i,j]]
                    b_k = [x[i,j+1] - x[i+1,j]; y[i,j+1] - y[i+1,j]; z[i,j+1] - z[i+1,j]]
                    (nx, ny, nz) = cross(a_k,b_k)/norm(cross(a_k,b_k))
                    
                    push!(vr_p, ThreeDSurfVRingPanel(0., xc, yc, zc, 0., 0., 0., nx, ny, nz)) 
                end
            end

            #Set grid properties of the rings
            for i = 1:nchord
                for j = 1:nspan + 1
                    # Vortex Ring Vertice Position
                    xv = cp/4. + x[i,j]
                    yv = y[i,j]
                    zv = z[i,j]
                    
                    push!(vr_g, ThreeDSurfVRingGrid(xv, yv, zv, 0., 0., 0.)) 
                end
            end
            
            #Set the wake properties
            i = nchord + 1 #chord+1 is the wake vortex used in IC calculation
            dtref = 0.015*cref/uref
            for j = 1:nspan+1
                xv = x[i,j] + tevloc*uref*dtref
                yv = y[i,j]
                zv = z[i,j]
                push!(vr_g, ThreeDSurfVRingGrid(xv, yv, zv, 0., 0., 0.))
            end                      
            
        elseif distype == "Cosine" || distype == "cosine"
            # Panel edges
            thetax = zeros(nchord+1)
            dthetax = pi/nchord
            thetay = zeros(nspan+1)
            dthetay = pi/nspan
            x = zeros(nchord+1, nspan+1)
            y = zeros(nchord+1, nspan+1)
            z = zeros(nchord+1, nspan+1)
            
            for i = 1:nchord+1
                for j = 1:nspan+1
                    thetax[i] = real(i - 1.)*dthetax
                    thetay[j] = real(j - 1.)*dthetay
                    x[i,j] = (cref/2.)*(1 - cos(thetax[i]))
                    y[i,j] = (bref/4.)*(1 - cos(thetay[j]))
                    z[i,j] = 0.
                end
            end
            
            #Votex rings definitions
            vr_g = ThreeDSurfVRingGrid[]
            vr_p = ThreeDSurfVRingPanel[]            
            
            for i = 1:nchord
                for j = 1:nspan
                  
                    bp = y[i,j+1] - y[i,j]
                    cp = x[i+1,j] - x[i,j]
                    
                    # Collocation Points
                    xc = 3.*cp/4. + x[i,j]
                    yc = bp/2. + y[i,j]
                    zc = z[i,j]
                    
                    # Normal and Tangent Vectors at Collocation Points
                    a_k = [x[i+1,j+1] - x[i,j]; y[i+1,j+1] - y[i,j]; z[i+1,j+1] - z[i,j]]
                    b_k = [x[i,j+1] - x[i+1,j]; y[i,j+1] - y[i+1,j]; z[i,j+1] - z[i+1,j]]
                    (nx, ny, nz) = cross(a_k,b_k)/norm(cross(a_k,b_k))
                    
                    push!(vr_p, ThreeDSurfVRingPanel( 0., xc, yc, zc, 0., 0., 0., nx, ny, nz)) 
                end
            end

            # Vortex Ring Vertice Position
            for i = 1:nchord
                for j = 1:nspan + 1
                    xv = cp/4. + x[i,j]
                    yv = y[i,j]
                    zv = z[i,j]
                    push!(vr_g, ThreeDSurfVRingGrid(xv, yv, zv, 0., 0., 0.))
                end
            end

            #Set the wake properties
            i = nchord + 1 #chord+1 is the wake vortex used in IC calculation
            dtref = 0.015*cref/uref
            for j = 1:nspan+1
                xv = x[i,j] + tevloc*uref*dtref
                yv = y[i,j]
                zv = z[i,j]
                push!(vr_g, ThreeDSurfVRingGrid(xv, yv, zv, 0., 0., 0.))
            end
                
        else
            error("Invalid distribution type")
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
            kinem.h = kindef.h(0.)*cref
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

phi = 0.
theta = kinem.alpha
psi = 0.

p = 0.
q = kinem.alphadot
R = 0.

tlg = [1. 0. 0.; 0 cos(phi) sin(phi); 0 -sin(phi) cos(phi)]*[cos(theta) 0 -sin(theta); 0 1 0; sin(theta) 0 cos(theta)]*[cos(psi) sin(psi) 0; -sin(psi) cos(psi) 0; 0 0 1]
tgl = transpose(tlg)
                    
rhs = zeros(nchord*nspan)
levflag = [9]
a0 = zeros(nspan)

X0 = zeros(3)
X0[3] = kinem.h

IC = zeros(nchord*nspan, nspan*nchord)
IC_W = zeros(nchord*nspan, nspan*nchord)

new(cref, AR, uref, pvt, lespcrit, coord_file,  nspan, nchord, kinem, kindef, a0, vr_g, vr_p, distype, IC, IC_W, rhs, tgl, tlg, X0, tevloc, x, y, z)
end
end
# ---------------------------------------------------------------------------------------------

type ThreeDWakeVRingGrid
    xv_I :: Float64
    yv_I :: Float64
    zv_I :: Float64
end    

immutable ThreeDFlowFieldVR
    nspan :: Int
    velX :: MotionDef
    velY :: MotionDef
    velZ :: MotionDef
    u :: Vector{Float64}
    v :: Vector{Float64}
    w :: Vector{Float64}
    tev :: Vector{ThreeDWakeVRingGrid}
    tev_s :: Vector{Float64}
    lev :: Vector{ThreeDWakeVRingGrid}
    lev_s :: Vector{Float64}
    extv :: Vector{ThreeDWakeVRingGrid}
    extv_s :: Vector{Float64}
    
    # Wake rings are stored as ordered vectors instead of arrays for ease of dynamic memory allocation
    
    function ThreeDFlowFieldVR(nspan, velX = ConstDef(0.), velY = ConstDef(0.), velZ = ConstDef(0.))
        u = [0;]
        v = [0;]
        w = [0;]
        tev = ThreeDWakeVRingGrid[]
        lev = ThreeDWakeVRingGrid[]
        extv = ThreeDWakeVRingGrid[]
        tev_s = Float64[]
        lev_s = Float64[]
        extv_s = Float64[]
        new(nspan, velX, velY, velZ, u, v, w, tev, tev_s, lev, lev_s, extv, extv_s)
    end
end
# ---------------------------------------------------------------------------------------------
