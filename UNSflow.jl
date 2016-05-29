module UNSflow

using Dierckx
export Spline1D, derivative, evaluate

using ForwardDiff
export derivative

#using Roots
#export secant_method
#using Roots

using PyPlot
export plot, scatter, figure

using Debug

using NLsolve
export nlsolve, not_in_place

export camber_calc, update_boundpos, update_kinem, update_indbound, update_downwash, update_a0anda1, place_tev, update_a2toan, mutual_ind, trapz, update_a2a3adot, update_bv, ind_vel, view_vorts, wakeroll, lautat, lautat_wakeroll, ldvm

export KinemPar, MotionDef, KinemDef, EldUpDef, ConstDef, TwoDSurf, TwoDVort, TwoDFlowField, KelvinCondition, KelvinKutta


function camber_calc(x::Vector,airfoil::ASCIIString)
    #Determine camber and camber slope on airfoil from airfoil input file

    ndiv = length(x);
    c = x[ndiv];

    cam = Array(Float64,ndiv)
    cam_slope = Array(Float64,ndiv)

    in_air = readdlm(airfoil);
    xcoord = in_air[:,1];
    ycoord = in_air[:,2];
    ncoord = length(xcoord);
    xcoord_sum = Array(Float64,ncoord);
    xcoord_sum[1] = 0;
    for i = 1:ncoord-1
        xcoord_sum[i+1] = xcoord_sum[i] + abs(xcoord[i+1]-xcoord[i]);
    end
    y_spl = Spline1D(xcoord_sum,ycoord);
    y_ans = Array(Float64,2*ndiv);

    for i=1:ndiv
        y_ans[i] = evaluate(y_spl,x[i]/c);
    end
    for i=ndiv+1:2*ndiv
        y_ans[i] = evaluate(y_spl,(x[ndiv]/c) + (x[i-ndiv]/c));
    end
    cam[1:ndiv] = [(y_ans[i] + y_ans[(2*ndiv) + 1 - i])*c/2 for i = ndiv:-1:1];
    cam[1] = 0;
    cam_spl = Spline1D(x,cam);
    cam_slope[1:ndiv] = Dierckx.derivative(cam_spl,x);
    return cam, cam_slope
end


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

immutable ConstDef <: MotionDef
    amp :: Float64
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

    function TwoDSurf(c, uref, coord_file, pvt, ndiv, naterm, dynamics_type, kindef)
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

        kinem.alpha = kindef.alpha(0.)
        kinem.alphadot = ForwardDiff.derivative(kindef.alpha,0.)*uref/c
        kinem.h = kindef.h(0.)*c
        kinem.hdot = ForwardDiff.derivative(kindef.h,0.)*uref
        kinem.u = kindef.u(0.)*uref
        kinem.udot = ForwardDiff.derivative(kindef.u,0.)*uref*uref/c

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
        lespcrit = zeros(1)
        levflag = [0]
        new(c, uref, coord_file, pvt, ndiv, naterm, dynamics_type, kindef, cam, cam_slope, theta, x, kinem, bnd_x, bnd_z, uind, wind, downwash, a0, aterm, a0dot, adot, a0prev, aprev, bv,lespcrit,levflag)
    end
end


function update_kinem(surf::TwoDSurf, t)
    #Dimensional values of kinemtic parameters
    surf.kinem.alpha = surf.kindef.alpha(t)
    surf.kinem.alphadot = ForwardDiff.derivative(surf.kindef.alpha,t)*surf.uref/surf.c
    surf.kinem.h = surf.kindef.h(t)*surf.c
    surf.kinem.hdot = ForwardDiff.derivative(surf.kindef.h,t)*surf.uref
    surf.kinem.u = surf.kindef.u(t)*surf.uref
    surf.kinem.udot = ForwardDiff.derivative(surf.kindef.u,t)*surf.uref*surf.uref/surf.c
    return surf
end

function update_bv(surf::TwoDSurf)
    gamma = zeros(surf.ndiv)
    for ib = 1:surf.ndiv
        gamma[ib] = (surf.a0[1]*(1 + cos(surf.theta[ib])))
        for ia = 1:surf.naterm
            gamma[ib] = gamma[ib] + surf.aterm[ia]*sin(ia*surf.theta[ib])*sin(surf.theta[ib])
        end
        gamma[ib] = gamma[ib]*surf.uref*surf.c
    end

    for ib = 2:surf.ndiv
        surf.bv[ib-1].s = (gamma[ib]+gamma[ib-1])*(surf.theta[2]-surf.theta[1])/2.
        surf.bv[ib-1].x = (surf.bnd_x[ib] + surf.bnd_x[ib-1])/2.
        surf.bv[ib-1].z = (surf.bnd_z[ib] + surf.bnd_z[ib-1])/2.
    end
end

function update_indbound(surf::TwoDSurf, curfield::TwoDFlowField)
    surf.uind[1:surf.ndiv], surf.wind[1:surf.ndiv] = ind_vel([curfield.tev; curfield.lev], surf.bnd_x, surf.bnd_z)
    return surf
end

function update_downwash(surf::TwoDSurf)
    for ib = 1:surf.ndiv
        surf.downwash[ib] = -surf.kinem.u*sin(surf.kinem.alpha) - surf.uind[ib]*sin(surf.kinem.alpha) + surf.kinem.hdot*cos(surf.kinem.alpha) - surf.wind[ib]*cos(surf.kinem.alpha) - surf.kinem.alphadot*(surf.x[ib] - surf.pvt*surf.c) + surf.cam_slope[ib]*(surf.uind[ib]*cos(surf.kinem.alpha) + surf.kinem.u*cos(surf.kinem.alpha) + surf.kinem.hdot*sin(surf.kinem.alpha) - surf.wind[ib]*sin(surf.kinem.alpha))
    end
    return surf
end

function update_a0anda1(surf::TwoDSurf)
    surf.a0[1] = trapz(surf.downwash,surf.theta)
    surf.aterm[1] = trapz(surf.downwash.*cos(surf.theta),surf.theta)
    surf.a0[1] = -surf.a0[1]/(surf.uref*pi)
    surf.aterm[1] = 2.*surf.aterm[1]/(surf.uref*pi)
    return surf
end


function update_a2toan(surf::TwoDSurf)

    for ia = 2:surf.naterm
        surf.aterm[ia] = trapz(surf.downwash.*cos(ia*surf.theta),surf.theta)
        surf.aterm[ia] = 2.*surf.aterm[ia]/(surf.uref*pi)
    end
    return surf
end

function update_a2a3adot(surf::TwoDSurf,dt)
    for ia = 2:3
        surf.aterm[ia] = trapz(surf.downwash.*cos(ia*surf.theta),surf.theta)
        surf.aterm[ia] = 2.*surf.aterm[ia]/(surf.uref*pi)
    end
    surf.a0dot[1] = (surf.a0[1] - surf.a0prev[1])/dt
    for ia = 1:3
        surf.adot[ia] = (surf.aterm[ia]-surf.aprev[ia])/dt
    end
    return surf
end



function trapz{T<:Real}(y::Vector{T}, x::Vector{T})
    local len = length(y)
    if (len != length(x))
        error("Vectors must be of same length")
    end
    r = 0.0
    for i in 2:len
        r += (x[i] - x[i-1]) * (y[i] + y[i-1])
    end
    r/2.0
end


immutable KelvinCondition
    surf :: TwoDSurf
    field :: TwoDFlowField
end

function call(kelv::KelvinCondition, tev_iter::Float64)
    #Update the TEV strength
    nlev = length(kelv.field.lev)
    ntev = length(kelv.field.tev)
    kelv.field.tev[ntev].s = tev_iter

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


function view_vorts(surf::TwoDSurf, field::TwoDFlowField)
    scatter(map(q->q.x, field.tev),map(q->q.z,field.tev))
    scatter(map(q->q.x, field.lev),map(q->q.z,field.lev))
    plot(map(q->q.x, surf.bv),map(q->q.z,surf.bv))
end


function lautat(surf::TwoDSurf)
    outfile = open("results.dat", "w")

    dtstar = 0.015
    dt = dtstar*surf.c/surf.uref
    nsteps = 500
    t = 0.

    curfield = TwoDFlowField()

    #Intialise flowfield
    for istep = 1:nsteps
        #Udpate current time
        t = t + dt

        #Update kinematic parameters
        update_kinem(surf, t)

        #Update bound vortex positions
        update_boundpos(surf, dt)

        #Add a TEV with dummy strength
        place_tev(surf,curfield,dt)

        kelv = KelvinCondition(surf,curfield)
        #Solve for TEV strength to satisfy Kelvin condition
        #curfield.tev[length(curfield.tev)].s = secant_method(kelv, 0., -0.01)
        soln = nlsolve(not_in_place(kelv), -0.01)
        curfield.tev[length(curfield.tev)].s = soln.zeros

        #Update adot
        update_a2a3adot(surf,dt)

        #Check for LEV and shed if yes

        #Update rest of Fourier terms
        #update_a2toan(surf)

        #Calculate bound vortex strengths
        #update_bv(surf)

        #wakeroll(surf, curfield)

        write(outfile, join((t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, surf.a0[1])," "), "\n")

    end
    close(outfile)

    #Plot flowfield viz and A0 history
    figure(0)
    view_vorts(surf, curfield)

    figure(1)
    data = readdlm("results.dat")
    plot(data[:,1],data[:,5])
end

function lautat_wakeroll(surf::TwoDSurf)
    outfile = open("results.dat", "w")

    dtstar = 0.015
    dt = dtstar*surf.c/surf.uref
    nsteps = 500
    t = 0.

    curfield = TwoDFlowField()

    #Intialise flowfield
    for istep = 1:nsteps
        #Udpate current time
        t = t + dt

        #Update kinematic parameters
        update_kinem(surf, t)

        #Update bound vortex positions
        update_boundpos(surf, dt)

        #Add a TEV with dummy strength
        place_tev(surf,curfield,dt)

        kelv = KelvinCondition(surf,curfield)
        #Solve for TEV strength to satisfy Kelvin condition
        curfield.tev[length(curfield.tev)].s = secant_method(kelv, 0., -0.01)

        #Update adot
        update_a2a3adot(surf,dt)

        #Check for LEV and shed if yes

        #Update rest of Fourier terms
        update_a2toan(surf)

        #Calculate bound vortex strengths
        update_bv(surf)

        wakeroll(surf, curfield, dt)

        write(outfile, join((t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, surf.a0[1])," "), "\n")

    end
    close(outfile)

    #Plot flowfield viz and A0 history
    figure(0)
    view_vorts(surf, curfield)

    figure(1)
    data = readdlm("results.dat")
    plot(data[:,1],data[:,5])
end

function ldvm(surf::TwoDSurf, curfield::TwoDFlowField)
    outfile = open("results.dat", "w")

    dtstar = 0.015
    dt = dtstar*surf.c/surf.uref
    nsteps = 500
    t = 0.

    #Intialise flowfield
    for istep = 1:nsteps
        #Udpate current time
        t = t + dt

        #Update kinematic parameters
        update_kinem(surf, t)

        #Update bound vortex positions
        update_boundpos(surf, dt)

        #Add a TEV with dummy strength
        place_tev(surf,curfield,dt)

        kelv = KelvinCondition(surf,curfield)
        #Solve for TEV strength to satisfy Kelvin condition
        curfield.tev[length(curfield.tev)].s = secant_method(kelv, 0., -0.01)

        #Check for LESP condition
        #Update values with converged value of shed tev
        #Update incduced velocities on airfoil
        update_indbound(kelv.surf, kelv.field)

        #Calculate downwash
        update_downwash(kelv.surf)

        #Calculate first two fourier coefficients
        update_a0anda1(kelv.surf)

        lesp = surf.a0[1]

        #2D iteration if LESP_crit is exceeded
        if (abs(lesp)>surf.lespcrit[1])
            #Add a TEV with dummy strength
            place_tev(surf,curfield,dt)

            #Add a LEV with dummy strength
            place_lev(surf,curfield,dt)

            kelvkutta = KelvinKutta(surf,curfield)
            #Solve for TEV and LEV strengths to satisfy Kelvin condition and Kutta condition at leading edge

            soln = nlsolve(not_in_place(kelvkutta), [-0.01; 0.01])
            (curfield.tev[length(curfield.tev)].s, curfield.lev[length(curfield.lev)].s) = soln.zero

            surf.levflag[1] = 1
        else
            surf.levflag[1] = 0
        end

        #Update adot
        update_a2a3adot(surf,dt)

        #Check for LEV and shed if yes

        #Update rest of Fourier terms
        update_a2toan(surf)

        #Calculate bound vortex strengths
        update_bv(surf)

        wakeroll(surf, curfield, dt)
        write(outfile, join((t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, surf.a0[1])," "), "\n")
    end

    close(outfile)
end


function wakeroll(surf::TwoDSurf, curfield::TwoDFlowField, dt)
    #Clean induced velocities
    for i = 1:length(curfield.tev)
        curfield.tev[i].vx = 0
        curfield.tev[i].vz = 0
    end

    for i = 1:length(curfield.lev)
        curfield.lev[i].vx = 0
        curfield.lev[i].vz = 0
    end

    #Velocities induced by free vortices on each other
    mutual_ind([curfield.tev; curfield.lev])

    #Add the influence of velocities induced by bound vortices
    utemp = zeros(length(curfield.tev)+length(curfield.lev))
    wtemp = zeros(length(curfield.tev)+length(curfield.lev))
    utemp, wtemp = ind_vel(surf.bv, [map(q -> q.x, curfield.tev); map(q -> q.x, curfield.lev)], [map(q -> q.z, curfield.tev); map(q -> q.z, curfield.lev)])

    for i = 1:length(curfield.tev)
        curfield.tev[i].vx += utemp[i]
        curfield.tev[i].vz += wtemp[i]
    end
    for i = length(curfield.tev)+1:length(utemp)
        curfield.lev[i-length(curfield.tev)].vx += utemp[i]
        curfield.lev[i-length(curfield.tev)].vz += wtemp[i]
    end

    #Convect free vortices with their induced velocities
    for i = 1:length(curfield.tev)
        curfield.tev[i].x += dt*curfield.tev[i].vx
        curfield.tev[i].z += dt*curfield.tev[i].vz
    end
    for i = 1:length(curfield.lev)
            curfield.lev[i].x += dt*curfield.lev[i].vx
        curfield.lev[i].z += dt*curfield.lev[i].vz
    end
    return curfield
end

function place_tev(surf::TwoDSurf,field::TwoDFlowField,dt)
    ntev = length(field.tev)
    if ntev == 0
        xloc = surf.bnd_x[surf.ndiv] + 0.5*surf.kinem.u*dt
        zloc = surf.bnd_z[surf.ndiv]
        else
        xloc = surf.bnd_x[surf.ndiv]+(1./3.)*(field.tev[ntev].x - surf.bnd_x[surf.ndiv])

        zloc = surf.bnd_z[surf.ndiv]+(1./3.)*(field.tev[ntev].z - surf.bnd_z[surf.ndiv])
    end
    push!(field.tev,TwoDVort(xloc,zloc,0.,0.02*surf.c,0.,0.))
    return field
end

function place_lev(surf::TwoDSurf,field::TwoDFlowField,dt)
    nlev = length(field.lev)

    le_vel_x = surf.kinem.u - surf.kinem.alphadot*sin(surf.kinem.alpha)*surf.pvt*surf.c + surf.uind[1]
    le_vel_z = -surf.kinem.alphadot*cos(surf.kinem.alpha)*surf.pvt*surf.c- surf.kinem.hdot + surf.wind[1]

    if (surf.levflag[1] == 0) then
        xloc = surf.bnd_x[1] + 0.5*le_vel_x*dt
        zloc = surf.bnd_z[1] + 0.5*le_vel_z*dt
    else
        xloc = surf.bnd_x[1]+(1./3.)*(field.lev[nlev].x - surf.bnd_x[1])
        zloc = surf.bnd_z[1]+(1./3.)*(field.lev[nlev].z - surf.bnd_z[1])
    end

    push!(field.lev,TwoDVort(xloc,zloc,0.,0.02*surf.c,0.,0.))

    return field
end


function call(eld::EldUpDef, t)
    sm = pi*pi*eld.K/(2*(eld.amp*pi/180)*(1 - eld.a))
    t1=1.
    t2=t1 + ((eld.amp*pi/180)/(2*eld.K))
    ((eld.K/sm)*log(cosh(sm*(t - t1))/cosh(sm*(t - t2))))+(eld.amp*pi/360)
end

function call(cons::ConstDef, t)
    cons.amp
end


function eld_fn(t::Vector;K=0.2,amp=45,t1=1.,tf=1.,a=11.)
    fr = K/(pi*abs(amp)*pi/180);
    t2=t1+(1./(2*pi*fr));
    t3 = t2+((1/(4*fr))-(1/(2*pi*fr)));
    t4 = t3+(1./(2*pi*fr));
    t5 = t4+tf;

    nstep = length(t)
    g = Array(Float64,nstep)
    res = Array(Float64,nstep)

    for i = 1:nstep
        g[i] = log((cosh(a*(t[i] - t1))*cosh(a*(t[i] - t4)))/(cosh(a*(t[i] - t2))*cosh(a*(t[i] - t3))))
    end
    maxg = maximum(g);
    return res = amp*g/maxg;
end

function start_bound(alpha,h,pvt,ndiv,c,x,cam)
    bnd_x = Array(Float64,ndiv)
    bnd_z = Array(Float64,ndiv)

    for i = 1:ndiv
        bnd_x[i] = -((c - pvt*c)+((pvt*c - x[i])*cos(alpha))) + (cam[i]*sin(alpha))
        bnd_z[i] = h + ((pvt*c - x[i])*sin(alpha))+(cam[i]*cos(alpha))
    end
    return bnd_x, bnd_z
end

function mutual_ind(vorts::Vector{TwoDVort})
    for i = 1:length(vorts)
        for j = i+1:length(vorts)
            dx = vorts[i].x - vorts[j].x
            dz = vorts[i].z - vorts[j].z
            #source- tar
            dsq = dx*dx + dz*dz

            magitr = 1./(2*pi*sqrt(vorts[j].vc*vorts[j].vc*vorts[j].vc*vorts[j].vc + dsq*dsq))
            magjtr = 1./(2*pi*sqrt(vorts[i].vc*vorts[i].vc*vorts[i].vc*vorts[i].vc + dsq*dsq))

            vorts[j].vx -= dz * vorts[i].s * magjtr
            vorts[j].vz += dx * vorts[i].s * magjtr

            vorts[i].vx += dz * vorts[j].s * magitr
            vorts[i].vz -= dx * vorts[j].s * magitr
        end
    end
    return vorts
end



function update_boundpos(surf::TwoDSurf, dt::Float64)
    for i = 1:surf.ndiv
        surf.bnd_x[i] = surf.bnd_x[i] + dt*((surf.pvt*surf.c - surf.x[i])*sin(surf.kinem.alpha)*surf.kinem.alphadot - surf.kinem.u + surf.cam[i]*cos(surf.kinem.alpha)*surf.kinem.alphadot)
        surf.bnd_z[i] = surf.bnd_z[i] + dt*(surf.kinem.hdot + (surf.pvt*surf.c - surf.x[i])*cos(surf.kinem.alpha)*surf.kinem.alphadot - surf.cam[i]*sin(surf.kinem.alpha)*surf.kinem.alphadot)
    end
    return surf
end

function ind_vel(src::Vector{TwoDVort},t_x,t_z)
    #'s' stands for source and 't' for target
    uind = zeros(length(t_x))
    wind = zeros(length(t_x))

    for itr = 1:length(t_x)
        for isr = 1:length(src)
            xdist = src[isr].x - t_x[itr]
            zdist = src[isr].z - t_z[itr]
            distsq = xdist*xdist + zdist*zdist
            uind[itr] = uind[itr] - src[isr].s*zdist/(2*pi*sqrt(src[isr].vc*src[isr].vc*src[isr].vc*src[isr].vc + distsq*distsq))
            wind[itr] = wind[itr] + src[isr].s*xdist/(2*pi*sqrt(src[isr].vc*src[isr].vc*src[isr].vc*src[isr].vc + distsq*distsq))
        end
    end
    return uind, wind
end

function convect_vort(v_x,v_z,uind,wind,dt)
    for iv = 1:length(v_x)
        v_x[iv] = v_x[iv] + dt*uind[iv]
        v_z[iv] = v_z[iv] + dt*wind[iv]
    end
    return v_x, v_z
end


function force_calc(u,alpha,hdot,u_ref,c,a0,aterm,a0dot, adot,uind,wind,ndiv,bv_s,x,cm_pvt)

    cnc = 2*pi*(u*cos(alpha)/u_ref + hdot*sin(alpha)/u_ref)*(a0 + aterm[1]/2)
    cnnc = 2*pi*(3*c*a0dot/(4*u_ref) + c*adot[1]/(4*u_ref) + c*adot[2]/(8*u_ref))
    cs = 2*pi*a0*a0
    #The components of normal force and moment from induced velocities are calulcated in dimensional units and nondimensionalized later
    nonl = 0
    nonl_m = 0
    for ib = 2:ndiv
        nonl = nonl + (uind[ib]*cos(alpha) - wind[ib]*sin(alpha))*bv_s[ib]
        nonl_m = nonl_m + (uind[ib]*cos(alpha) - wind[ib]*sin(alpha))*x[ib]*bv_s[ib]
    end
    nonl = nonl*2/(u_ref*u_ref*c)
    nonl_m = nonl_m*2/(u_ref*u_ref*c*c)

    cn = cnc+cnnc+nonl
    cl = cn*cos(alpha) + cs*sin(alpha)
    cd = cn*sin(alpha) - cs*cos(alpha)
    #Pitching moment is clockwise or nose up positive
    cm = cn*cm_pvt - 2*pi*((u*cos(alpha)/u_ref + hdot*sin(alpha)/u_ref)*(a0/4. + aterm[1]/4. - aterm[2]/8.) + c/u_ref*(7*a0dot/16 + 3*adot[1]/16 + adot[2]/16 - adot[3]/64)) - nonl_m
    return cl, cd, cm
end

function calc_bndvorts(a0,aterm,naterm,theta,dtheta,ndiv,bnd_x,bnd_z,u_ref,c)
    gamma = Array(Float64,ndiv)
    bv_s = Array(Float64,ndiv)
    bv_x = Array(Float64,ndiv)
    bv_z = Array(Float64,ndiv)
    bv_vc = Array(Float64,ndiv)

    for ib = 1:ndiv
        gamma[ib] = (a0*(1 + cos(theta[ib])))
        for ia = 1:naterm
            gamma[ib] = gamma[ib] + aterm[ia]*sin(ia*theta[ib])*sin(theta[ib])
        end
        gamma[ib] = gamma[ib]*u_ref*c
    end

    for ib = 2:ndiv
        bv_s[ib] = (gamma[ib]+gamma[ib-1])*dtheta/2.
        bv_x[ib] = (bnd_x[ib] + bnd_x[ib-1])/2.
        bv_z[ib] = (bnd_z[ib] + bnd_z[ib-1])/2.
        bv_vc[ib] = 0.02*c
    end
    return bv_s, bv_x, bv_z, bv_vc
end


# function calc_lesp_constU()
#     lesp=aterm()0
# end

# function twod_iter()
#     if (abs(lesp)>lesp_crit)
#         if (lesp>0)
#             lesp_cond=lesp_crit
#         else
#             lesp_cond=-lesp_crit
#         end
#         add_lev()
#         add_tev()
#         nlsolve(kelva0, [0.;0.], autodiff = true)
#     else
#         levflag = 0
#     end
# end

# function kelva0()
#     push!(tev_s, tev_iter)
#     push!(lev_s, lev_iter)
#     kelv = kelv_enf
#     for i = 1:nlev
#         kelv = kelv + lev_s(i)
#     end

#     for i = 1:ntev
#         kelv = kelv + tev_s(i)
#     end
#     gam = bcirc()
#     kelv = kelv + gam
#     kelva0_vec[1] = kelv
#     kelva0_vec[1] = aterm(0)-lesp_cond
# end

# function LDVM()
# end



end
