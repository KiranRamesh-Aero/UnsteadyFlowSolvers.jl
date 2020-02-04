"""
    TwoDVort(x,z,s,vc,vx,vz)

Defines a 2D vortex at `(x,z)` with vorticity `s` and vortex-core
radius `vc`.

`vx` and `vz` are induced velocity at centroid used in vortex
interaction calculations

"""
mutable struct  TwoDVort
    x :: Float64
    z :: Float64
    s :: Float64
    vc :: Float64
    vx :: Float64
    vz :: Float64
end

struct TwoDFlowField
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

struct TwoDSurf
    c :: Float64
    uref :: Float64
    coord_file :: String
    pvt :: Float64
    ndiv :: Int16
    naterm :: Int16
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
    rho :: Float64

    function TwoDSurf(coord_file, pvt, kindef,lespcrit=zeros(1); c=1., uref=1., ndiv=70, naterm=35, initpos = [0.; 0.], rho = 0.04, camberType = "radial")
        theta = zeros(ndiv)
        x = zeros(ndiv)
        cam = zeros(ndiv)
        cam_slope = zeros(ndiv)
        bnd_x = zeros(ndiv)
        bnd_z = zeros(ndiv)
        kinem = KinemPar(0, 0, 0, 0, 0, 0)

        if camberType == "radial"
            dtheta = pi/(ndiv-1)
            for ib = 1:ndiv
                theta[ib] = (ib-1.)*dtheta
                x[ib] = c/2. *(1-cos(theta[ib]))
            end
        elseif camberType == "linear"
            dx = c / (ndiv-1)
            for ib = 2:ndiv
                x[ib] = x[ib-1] + dx
            end
        end

        if (coord_file != "FlatPlate")
            cam, cam_slope = camber_calc(x, coord_file)
        end
        #camspl = Spline1D(theta, cam)
        #for i = 1:ndiv
        #    cam_slope[i] = Dierckx.derivative(camspl, theta[i])
        #end

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
        elseif (typeof(kindef.alpha) == FileDef)
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
        elseif (typeof(kindef.h) == FileDef)
            kinem.h = kindef.h(0.)
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
        elseif (typeof(kindef.alpha) == FileDef)
            kinem.u = kindef.u(0.)*uref
            kinem.udot = ForwardDiff.derivative(kindef.u,0.)*uref*uref/c
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
        adot = zeros(naterm)
        a0prev = zeros(1)
        aprev = zeros(naterm)
        bv = TwoDVort[]
        for i = 1:ndiv-1
            push!(bv,TwoDVort(0,0,0,0.02*c,0,0))
        end
        levflag = [0]
        new(c, uref, coord_file, pvt, ndiv, naterm, kindef, cam, cam_slope, theta, x, kinem, bnd_x, bnd_z, uind, wind, downwash, a0, aterm, a0dot, adot, a0prev, aprev, bv, lespcrit, levflag, initpos, rho)
    end
end

mutable struct TwoDOFPar
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
mutable struct KinemPar2DOF
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


struct KelvinCondition
    surf :: TwoDSurf
    field :: TwoDFlowField
end

function (kelv::KelvinCondition)(tev_iter::Array{Float64})

    nlev = length(kelv.field.lev)
    ntev = length(kelv.field.tev)

    uprev, wprev = ind_vel([kelv.field.tev[ntev]], kelv.surf.bnd_x, kelv.surf.bnd_z)

    #Update the TEV strength
    kelv.field.tev[ntev].s = tev_iter[1]

    unow, wnow = ind_vel([kelv.field.tev[ntev]], kelv.surf.bnd_x, kelv.surf.bnd_z)

    kelv.surf.uind[:] = kelv.surf.uind[:] .- uprev .+ unow
    kelv.surf.wind[:] = kelv.surf.wind[:] .- wprev .+ wnow

    #Calculate downwash
    update_downwash(kelv.surf, [kelv.field.u[1],kelv.field.w[1]])

    #Calculate first two fourier coefficients
    update_a0anda1(kelv.surf)

    val = kelv.surf.uref*kelv.surf.c*pi*(kelv.surf.a0[1] + kelv.surf.aterm[1]/2.) -
        kelv.surf.uref*kelv.surf.c*pi*(kelv.surf.a0prev[1] + kelv.surf.aprev[1]/2.) +
        kelv.field.tev[ntev].s

    return val
end

struct KelvinConditionMult
    surf :: Vector{TwoDSurf}
    field :: TwoDFlowField
end

function (kelv::KelvinConditionMult)(tev_iter::Array{Float64})

    nsurf = length(kelv.surf)

    val = zeros(nsurf)

    nlev = length(kelv.field.lev)
    ntev = length(kelv.field.tev)

    tev_list = kelv.field.tev[ntev-nsurf+1:ntev]
    for i = 1:nsurf

        bv_list = TwoDVort[]
        for j = 1:nsurf
            if i != j
                bv_list = [bv_list; kelv.surf[j].bv]
            end
        end
        uprev, wprev = ind_vel([tev_list; bv_list], kelv.surf[i].bnd_x, kelv.surf[i].bnd_z)

        #Update the TEV strength
        kelv.field.tev[ntev-nsurf+i].s = tev_iter[i]

        unow, wnow = ind_vel([tev_list; bv_list], kelv.surf[i].bnd_x, kelv.surf[i].bnd_z)

        kelv.surf[i].uind[:] = kelv.surf[i].uind[:] .- uprev .+ unow
        kelv.surf[i].wind[:] = kelv.surf[i].wind[:] .- wprev .+ wnow

        #Calculate downwash
        update_downwash(kelv.surf[i], [kelv.field.u[1],kelv.field.w[1]])
    end

    for i = 1:nsurf
        #Update Fourier coefficients and bv strength
        update_a0anda1(kelv.surf[i])
        update_a2toan(kelv.surf[i])
        update_bv(kelv.surf[i])


        val[i] = kelv.surf[i].uref*kelv.surf[i].c*pi*(kelv.surf[i].a0[1] + kelv.surf[i].aterm[1]/2.) -
            kelv.surf[i].uref*kelv.surf[i].c*pi*(kelv.surf[i].a0prev[1] + kelv.surf[i].aprev[1]/2.) +
            kelv.field.tev[ntev-nsurf+i].s
    end
    return val
end


struct KelvinKutta
    surf :: TwoDSurf
    field :: TwoDFlowField
end

function (kelv::KelvinKutta)(v_iter::Array{Float64})
    val = zeros(2)

    nlev = length(kelv.field.lev)
    ntev = length(kelv.field.tev)

    uprev, wprev = ind_vel([kelv.field.tev[ntev]; kelv.field.lev[nlev]], kelv.surf.bnd_x, kelv.surf.bnd_z)

    #Update the TEV and LEV strengths
    kelv.field.tev[ntev].s = v_iter[1]
    kelv.field.lev[nlev].s = v_iter[2]

    unow, wnow = ind_vel([kelv.field.tev[ntev]; kelv.field.lev[nlev]], kelv.surf.bnd_x, kelv.surf.bnd_z)

    kelv.surf.uind[:] = kelv.surf.uind[:] .- uprev .+ unow
    kelv.surf.wind[:] = kelv.surf.wind[:] .- wprev .+ wnow

    #Calculate downwash
    update_downwash(kelv.surf ,[kelv.field.u[1],kelv.field.w[1]])

    #Calculate first two fourier coefficients
    update_a0anda1(kelv.surf)

    val[1] = kelv.surf.uref*kelv.surf.c*pi*(kelv.surf.a0[1] + kelv.surf.aterm[1]/2.) -
        kelv.surf.uref*kelv.surf.c*pi*(kelv.surf.a0prev[1] + kelv.surf.aprev[1]/2.) +
        kelv.field.tev[ntev].s + kelv.field.lev[nlev].s

    if (kelv.surf.a0[1] > 0)
        lesp_cond = kelv.surf.lespcrit[1]
    else
        lesp_cond = -kelv.surf.lespcrit[1]
    end
    val[2] = kelv.surf.a0[1]-lesp_cond

    return val
end

struct KelvinKuttaMult
    surf :: Vector{TwoDSurf}
    field :: TwoDFlowField
    shed_ind :: Vector{Int}
end

function (kelv::KelvinKuttaMult)(v_iter::Array{Float64})
    nsurf = length(kelv.surf)
    nshed = length(kelv.shed_ind)
    val = zeros(nsurf+nshed)

    nlev = length(kelv.field.lev)
    ntev = length(kelv.field.tev)

    tev_list = kelv.field.tev[ntev-nsurf+1:ntev]
    lev_list = kelv.field.lev[nlev-nsurf.+kelv.shed_ind]

    levcount = 0
    for i = 1:nsurf
        bv_list = TwoDVort[]
        for j = 1:nsurf
            if i != j
                bv_list = [bv_list; kelv.surf[j].bv]
            end
        end
        uprev, wprev = ind_vel([tev_list; lev_list; bv_list], kelv.surf[i].bnd_x, kelv.surf[i].bnd_z)

        #Update the shed vortex strengths
        kelv.field.tev[ntev-nsurf+i].s = v_iter[i]

        if i in kelv.shed_ind
            levcount += 1
            kelv.field.lev[nlev-nsurf+i].s = v_iter[nsurf+levcount]
        end

        unow, wnow = ind_vel([tev_list; lev_list; bv_list], kelv.surf[i].bnd_x, kelv.surf[i].bnd_z)

        kelv.surf[i].uind[:] = kelv.surf[i].uind[:] .- uprev .+ unow
        kelv.surf[i].wind[:] = kelv.surf[i].wind[:] .- wprev .+ wnow

        #Calculate downwash
        update_downwash(kelv.surf[i], [kelv.field.u[1],kelv.field.w[1]])
    end

    levcount = 0
    for i = 1:nsurf
        #Update Fourier coefficients and bv strength
        update_a0anda1(kelv.surf[i])
        update_a2toan(kelv.surf[i])
        update_bv(kelv.surf[i])

        if i in kelv.shed_ind
            val[i] = kelv.surf[i].uref*kelv.surf[i].c*pi*(kelv.surf[i].a0[1] + kelv.surf[i].aterm[1]/2.) -
                kelv.surf[i].uref*kelv.surf[i].c*pi*(kelv.surf[i].a0prev[1] + kelv.surf[i].aprev[1]/2.) +
                kelv.field.tev[ntev-nsurf+i].s + kelv.field.lev[nlev-nsurf+i].s
            levcount += 1
            if (kelv.surf[i].a0[1] > 0)
                lesp_cond = kelv.surf[i].lespcrit[1]
            else
                lesp_cond = -kelv.surf[i].lespcrit[1]
            end

            val[levcount+nsurf] = kelv.surf[i].a0[1] - lesp_cond
        else
            val[i] = kelv.surf[i].uref*kelv.surf[i].c*pi*(kelv.surf[i].a0[1] + kelv.surf[i].aterm[1]/2.) -
                kelv.surf[i].uref*kelv.surf[i].c*pi*(kelv.surf[i].a0prev[1] + kelv.surf[i].aprev[1]/2.) +
                kelv.field.tev[ntev-nsurf+i].s
        end
    end
    return val
end

mutable struct meshgrid
    t :: Float64
    alpha :: Float64
    tev :: Array{Float64}
    lev :: Array{Float64}
    camX :: Vector{Float64}
    camZ :: Vector{Float64}
    x :: Array{Float64}
    z :: Array{Float64}
    uMat :: Array{Float64}
    wMat :: Array{Float64}
    velMag :: Array{Float64}

    function meshgrid(surf::TwoDSurf,tevs::Vector{TwoDVort},levs::Vector{TwoDVort},offset,t,width::Int64 = 100,view::String = "square")
        farBnd = surf.x[end] + surf.c*offset
        nearBnd = surf.x[1] - surf.c*offset
        zBnd = ( farBnd - nearBnd ) / 2
        if view == "wake"
            farBnd = 2*farBnd
        elseif view == "largewake"
            farBnd = 3*farBnd
            zBnd = 2*zBnd
        elseif view == "longwake"
            farBnd = 5*farBnd
            zBnd = 2*zBnd
        elseif view == "UI Window"
            farBnd = 2*farBnd
            zBnd = 3/4*farBnd
        end

        # Global frame translation
        X0 = -surf.kindef.u(t)*t
        Z0 = surf.kindef.h(t)

        lowX = nearBnd + X0 - surf.pvt*surf.c
        uppX = farBnd + X0 - surf.pvt*surf.c
        lowZ = -zBnd + Z0
        uppZ = zBnd + Z0

        # Finding global camber line postion
        camX, camZ = IFR(surf,surf.x,surf.cam,t)

        # Creating meshgrid
        step = (farBnd-nearBnd)/(width-1)
        range = lowX:step:uppX
        height = zBnd*2/step + 1
        x = [ j for i = 1:height, j = range]

        range = uppZ:-step:lowZ
        z = [ i for i = range , j = 1:size(x,2)]

        uMat = 0 .* x .+ surf.kinem.u
        wMat = 0 .* x .+ surf.kinem.hdot
        velMag = 0 .* x

        t = 0
        alpha = surf.kinem.alpha
        circ = zeros(surf.ndiv-1)

        tevX = []
        tevZ = []
        for i = 1:length(tevs)
            vx = tevs[i].x
            vz = tevs[i].z
            if vx >= lowX && vx <= uppX
                if vz >= lowZ && vz <= uppZ
                    tevX = [tevX;vx]
                    tevZ = [tevZ;vz]
                end
            end
        end
        tev = hcat(tevX,tevZ)
        levX = []
        levZ = []
        if length(levs) > 0
            for i = 1:length(levs)
                vx = levs[i].x
                vz = levs[i].z
                if vx >= lowX && vx <= uppX
                    if vz >= lowZ && vz <= uppZ
                        levX = [levX;vx]
                        levZ = [levZ;vz]
                    end
                end
            end
        end
        lev = hcat(levX,levZ)

        new(t,alpha,tev,lev,camX, camZ, x, z, uMat, wMat, velMag)
    end
end
