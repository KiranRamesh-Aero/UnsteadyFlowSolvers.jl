immutable TwoDSurfPanel
    c :: Float64
    uref :: Float64
    coord_file :: String
    pvt :: Float64
    ndiv :: Int64
    kindef :: KinemDef
    kinem :: KinemPar
    bnd_x :: Vector{Float64}
    bnd_z :: Vector{Float64}
    uind :: Vector{Float64}
    wind :: Vector{Float64}
    downwash :: Vector{Float64}
    a0 :: Vector{Float64}
    lespcrit :: Vector{Float64}
    levflag :: Vector{Int8}
    initpos :: Vector{Float64}
    x1p :: Vector{Float64}
    z1p :: Vector{Float64}
    x2p :: Vector{Float64}
    z2p :: Vector{Float64}
    alpha :: Vector{Float64}
    normx :: Vector{Float64}
    normz :: Vector{Float64}
    gam1 :: Vector{Float64}
    gam2 :: Vector{Float64}
    IC :: Array{Float64}

    function TwoDSurfPanel(coord_file, pvt, kindef, lespcrit=zeros(1); c=1., uref=1., ndiv=200, initpos = [0.; 0.])

        bnd_x = zeros(ndiv)
        bnd_z = zeros(ndiv)
        kinem = KinemPar(0, 0, 0, 0, 0, 0)

        in_air = readdlm(coord_file)
        xcoord = in_air[:,1]
        zcoord = in_air[:,2]
        ncoord = length(xcoord)

        s = zeros(ncoord)
        for i = 1:ncoord-1
            s[i+1] = s[i] + sqrt((xcoord[i+1] - xcoord[i])^2 + (zcoord[i+1] - zcoord[i])^2)
        end

        x_s = Spline1D(s, xcoord)
        z_s = Spline1D(s, zcoord)

        #Form panels

        x1p = zeros(ndiv-1)
        z1p = zeros(ndiv-1)
        x2p = zeros(ndiv-1)
        z2p = zeros(ndiv-1)
        alpha = zeros(ndiv-1)
        normx = zeros(ndiv-1)
        normz = zeros(ndiv-1)

        for i = 1:ndiv-1
            s1 = (i-1)*s[ncoord]/(ndiv-1)
            s2 = (i)*s[ncoord]/(ndiv-1)

            x1p[i] = Dierckx.evaluate(x_s, s1)
            x2p[i] = Dierckx.evaluate(x_s, s2)
            z1p[i] = Dierckx.evaluate(z_s, s1)
            z2p[i] = Dierckx.evaluate(z_s, s2)
            alpha[i] = atan((z2p[i] - z1p[i])/(x2p[i] - x1p[i]))
            normx[i] = -(z2p[i] - z1p[i])/sqrt((z2p[i] - z1p[i])^2 + (x2p[i] - x1p[i])^2)
            normz[i] = (x2p[i] - x1p[i])/sqrt((z2p[i] - z1p[i])^2 + (x2p[i] - x1p[i])^2)
        end

        #Form influence coefficients
        IC = zeros(ndiv, ndiv)

        for i = 1:ndiv-1
            for j = 1:ndiv-1
                if j == 1
                    _, _, ua, wa, _, _  = vor2dl(1., 1., 0.5*(x1p[i] + x2p[i]),
                     0.5*(z1p[i] + z2p[i]), x1p[j], z1p[j], x2p[j], z2p[j])

                    IC[i,j] = dot([ua; wa], [normx[i]; normz[i]])

                elseif j == ndiv-1
                    _, _, _, _, ub, wb = vor2dl(1., 1., 0.5*(x1p[i] + x2p[i]),
                    0.5*(z1p[i] + z2p[i]), x1p[j], z1p[j], x2p[j], z2p[j])

                    IC[i,j] = dot([ub; wb], [normx[i]; normz[i]])

                    #ua comes from j and ub comes from j-1
                    _, _, ua, wa, _, _ = vor2dl(1., 1., 0.5*(x1p[i] + x2p[i]), 0.5*(z1p[i] + z2p[i]), x1p[j], z1p[j], x2p[j], z2p[j])
                    _, _, _, _, ub, wb = vor2dl(1., 1., 0.5*(x1p[i] + x2p[i]), 0.5*(z1p[i] + z2p[i]), x1p[j-1], z1p[j-1], x2p[j-1], z2p[j-1])

                    IC[i,j] = dot([ua + ub; wa + wb], [normx[i], normz[i]])
                else
                    _, _, ua, wa, _, _ = vor2dl(1., 1., 0.5*(x1p[i] + x2p[i]), 0.5*(z1p[i] + z2p[i]), x1p[j], z1p[j], x2p[j], z2p[j])
                    _, _, _, _, ub, wb = vor2dl(1., 1., 0.5*(x1p[i] + x2p[i]), 0.5*(z1p[i] + z2p[i]), x1p[j-1], z1p[j-1], x2p[j-1], z2p[j-1])
                    IC[i,j] = dot([ua + ub; wa + wb], [normx[i], normz[i]])
                end
            end
        end
        IC[ndiv,1] = 1.
        IC[ndiv,ndiv] = 1.

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

        #for i = 1:ndiv
        #    bnd_x[i] = -((c - pvt*c)+((pvt*c - x[i])*cos(kinem.alpha))) + (cam[i]*sin(kinem.alpha)) + initpos[1]
        #    bnd_z[i] = kinem.h + ((pvt*c - x[i])*sin(kinem.alpha))+(cam[i]*cos(kinem.alpha)) + initpos[2]
        #end
        uind = zeros(ndiv)
        wind = zeros(ndiv)
        downwash = zeros(ndiv)
        a0 = zeros(1)

levflag = [0]

gam1 = zeros(ndiv)
gam2 = zeros(ndiv)

new(c, uref, coord_file, pvt, ndiv, kindef, kinem, bnd_x, bnd_z, uind, wind, downwash, a0, lespcrit, levflag, initpos, x1p, z1p, x2p, z2p, alpha, normx, normz, gam1, gam2, IC)
    end
end
