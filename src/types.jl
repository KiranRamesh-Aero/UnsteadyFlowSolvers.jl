type KinemPar
    alpha :: Float64
    h :: Float64
    alphadot :: Float64
    hdot :: Float64
    u :: Float64
    udot :: Float64
end


abstract MotionDef

immutable EldUpDef <: MotionDef
    amp :: Float64
    K :: Float64
    a :: Float64
end

immutable EldRampReturnDef <:MotionDef
    amp :: Float64
    K :: Float64
    a :: Float64
end

immutable ConstDef <: MotionDef
    amp :: Float64
end


function call(eld::EldUpDef, t)
    sm = pi*pi*eld.K/(2*(eld.amp)*(1 - eld.a))
    t1 = 1.
    t2 = t1 + ((eld.amp)/(2*eld.K))
    ((eld.K/sm)*log(cosh(sm*(t - t1))/cosh(sm*(t - t2))))+(eld.amp/2)
end

function call(cons::ConstDef, t)
    cons.amp
end

function call(eld::EldRampReturnDef, tt)
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
  w :: Float64
  phi :: Float64
end

immutable CosDef <: MotionDef
  mean :: Float64
  amp :: Float64
  w :: Float64
  phi :: Float64
end

function call(kin::SinDef, t)
  (kin.mean) + (kin.amp)*sin(kin.w*t + kin.phi)
end

function call(kin::CosDef, t)
  (kin.mean) + (kin.amp)*cos(kin.w*t + kin.phi)
end


type TwoDVort
    x :: Float64
    z :: Float64
    s :: Float64
    vc :: Float64
    vx :: Float64
    vz :: Float64
end

immutable TwoDFlowField
    velX :: Float64
    velZ :: Float64
    tev :: Vector{TwoDVort}
    lev :: Vector{TwoDVort}
    function TwoDFlowField()
        velX = 0
        velZ = 0
        tev = TwoDVort[]
        lev = TwoDVort[]
        new(velX, velZ, tev, lev)
    end
end


type KinemDef
    alpha :: MotionDef
    h :: MotionDef
    u :: MotionDef
end


immutable TwoDSurf
    c :: Float64
    uref :: Float64
    coord_file :: ASCIIString
    pvt :: Float64
    ndiv :: Int8
    naterm :: Int8
    dynamics_type :: ASCIIString
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

    function TwoDSurf(c, uref, coord_file, pvt, ndiv, naterm, dynamics_type, kindef,lespcrit=zeros(1))
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
            kinem.alpha, kinem.alphadot = kindef.alpha(0.)
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
            kinem.h, kinem.hdot = kindef.h(0.)
            kinem.h = kinem.h*c
            kinem.hdot = kinem.hdot*uref
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
        new(c, uref, coord_file, pvt, ndiv, naterm, dynamics_type, kindef, cam, cam_slope, theta, x, kinem, bnd_x, bnd_z, uind, wind, downwash, a0, aterm, a0dot, adot, a0prev, aprev, bv,lespcrit,levflag)
    end
end






immutable EldUpIntDef <: MotionDef
    amp :: Float64
    K :: Float64
    a :: Float64
end

function call(eld::EldUpIntDef, t)

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



immutable KelvinCondition
    surf :: TwoDSurf
    field :: TwoDFlowField
end

function call(kelv::KelvinCondition, tev_iter::Array{Float64})
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
    return val
end

immutable KelvinKutta
    surf :: TwoDSurf
    field :: TwoDFlowField
end

function call(kelv::KelvinKutta, v_iter::Array{Float64})
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
