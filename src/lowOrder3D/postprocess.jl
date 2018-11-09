function writeStamp(dirname::String, t::Float64, surf::ThreeDSurfSimple, curfield::ThreeDFieldSimple)
    dirvec = readdir()
    if dirname in dirvec
        rm(dirname, recursive=true)
    end
    mkdir(dirname)
    cd(dirname)
    
    f = open("time", "w")
    Serialization.Serialization.serialize(f, ["#time \t"])
    DelimitedFiles.writedlm(f, [t]')
    close(f)

    for i = 1:surf.nspan
        f = open("kinem-$i", "w")
        Serialization.Serialization.serialize(f, ["alpha (deg) \t", "h/c \t", "u/uref \t"])
        DelimitedFiles.writedlm(f, [surf.s2d[i].kinem.alpha, surf.s2d[i].kinem.h, surf.s2d[i].kinem.u]')
        close(f)

        f = open("FourierCoeffs-$i", "w")
        Serialization.serialize(f, ["#Fourier coeffs (0-n) \t", "d/dt (Fourier coeffs) \n"])
        matfour = zeros(surf.s2d[i].naterm+1, 2)
        matfour[:,1] = [surf.s2d[i].a0[1]; surf.s2d[i].aterm[:]]
        matfour[1:4,2] = [surf.s2d[i].a0dot[1];surf.s2d[i].adot[:]]
        DelimitedFiles.writedlm(f, matfour)
        close(f)

        f = open("tev-$i", "w")
        Serialization.serialize(f, ["#strength \t", "x-position \t", "z-position \n"])
        tevmat = zeros(length(curfield.f2d[i].tev), 3)
        for iv = 1:length(curfield.f2d[i].tev)
            tevmat[iv,:] = [curfield.f2d[i].tev[iv].s curfield.f2d[i].tev[iv].x curfield.f2d[i].tev[iv].z]
        end
        DelimitedFiles.writedlm(f, tevmat)
        close(f)

        f = open("lev-$i", "w")
        Serialization.serialize(f, ["#strength \t", "x-position \t", "z-position \n"])
        levmat = zeros(length(curfield.f2d[i].lev), 3)
        for iv = 1:length(curfield.f2d[i].lev)
            levmat[iv,:] = [curfield.f2d[i].lev[iv].s curfield.f2d[i].lev[iv].x curfield.f2d[i].lev[iv].z]
        end
        DelimitedFiles.writedlm(f, levmat)
        close(f)

        f = open("boundv-$i", "w")
        Serialization.serialize(f, ["#strength \t", "x-position \t", "z-position \n"])
        bvmat = zeros(length(surf.s2d[i].bv), 3)
        for iv = 1:length(surf.s2d[i].bv)
            bvmat[iv,:] = [surf.s2d[i].bv[iv].s surf.s2d[i].bv[iv].x surf.s2d[i].bv[iv].z]
        end
        DelimitedFiles.writedlm(f, bvmat)
        close(f)
    end

    #Write spanwise properties of interest in a file
    f = open("spanwise-var", "w")
    Serialization.serialize(f, ["#y_le \t", "shed TEV str \t", "shed LEV str \t", "Gamma \t", "A0 2D \t", "A03D \t", "A0 tot \n"])
    mat = zeros(0,7)
    mat = [mat; [-surf.AR/2, 0. , 0. , 0. , 0. , 0. , 0. ]']
    for i = 1:surf.nspan
        yle = surf.yle[i]
        tev_s = curfield.f2d[i].tev[end].s
        if surf.s2d[i]. levflag == 1
            lev_s = curfield.f2d[i].lev[end].s
        else
            lev_s = 0.
        end
        gam = surf.s2d[i].uref*surf.s2d[i].c*pi*(surf.s2d[i].a0[1] + 0.5*surf.s2d[i].aterm[1])
        a02d = surf.s2d[i].a0[1] - surf.a03d[i]

        mat = [mat; [yle, tev_s, lev_s, gam, a02d, surf.a03d[i], surf.s2d[i].a0[1]]']
    end
    DelimitedFiles.writedlm(f, mat)
    close(f)

    cd("..")
end

function calc_forces(surf::ThreeDSurfSimple, field :: ThreeDFieldSimple, dt :: Float64)

    cl = zeros(surf.nspan)
    cd = zeros(surf.nspan)
    cm = zeros(surf.nspan)

    lev_s = zeros(surf.nspan)
    for i = 1:surf.nspan
        if surf.s2d[i].levflag[1] == 1
            lev_s[i] = field.f2d[i].lev[end].s
        end
    end
    
    for i = 1:surf.nspan
        # First term in eqn (2.30) Ramesh et al. in coefficient form
        wi = 0
        for n = 1:surf.nspan
            nn = 2*n -2
            wi += real(nn)*surf.bcoeff[n]*sin(nn*surf.psi[i])/sin(surf.psi[i])
        end

        cnc = 2*pi*(surf.s2d[i].kinem.u*cos(surf.s2d[i].kinem.alpha)/surf.s2d[i].uref + surf.s2d[i].kinem.hdot*sin(surf.s2d[i].kinem.alpha)/surf.s2d[i].uref + wi*sin(surf.s2d[i].kinem.alpha))*(surf.s2d[i].a0[1] + surf.s2d[i].aterm[1]/2.)

        # Second term in eqn (2.30) Ramesh et al. in coefficient form
        cnnc = 2*pi*(3*surf.s2d[i].c*surf.s2d[i].a0dot[1]/(4*surf.s2d[i].uref) + surf.s2d[i].c*surf.s2d[i].adot[1]/(4*surf.s2d[i].uref) + surf.s2d[i].c*surf.s2d[i].adot[2]/(8*surf.s2d[i].uref)) + (2*lev_s[i]/(dt*surf.s2d[i].uref^2))
        
        # Suction force given in eqn (2.31) Ramesh et al.
        cs = 2*pi*surf.s2d[i].a0[1]*surf.s2d[i].a0[1]

        #The components of normal force and moment from induced velocities are calulcated in dimensional units and nondimensionalized later
        nonl=0
        nonl_m=0
        for ib = 1:surf.s2d[i].ndiv-1
            nonl = nonl + (surf.s2d[i].uind[ib]*cos(surf.s2d[i].kinem.alpha) - surf.s2d[i].wind[ib]*sin(surf.s2d[i].kinem.alpha))*surf.s2d[i].bv[ib].s
            nonl_m = nonl_m + (surf.s2d[i].uind[ib]*cos(surf.s2d[i].kinem.alpha) - surf.s2d[i].wind[ib]*sin(surf.s2d[i].kinem.alpha))*surf.s2d[i].x[ib]*surf.s2d[i].bv[ib].s
        end
        nonl = nonl*2. /(surf.s2d[i].uref*surf.s2d[i].uref*surf.s2d[i].c)
        nonl_m = nonl_m*2. /(surf.s2d[i].uref*surf.s2d[i].uref*surf.s2d[i].c*surf.s2d[i].c)
        
        # Normal force coefficient
        cn = cnc + cnnc + nonl
        
        # Lift and drag coefficients
        cl[i] = cn*cos(surf.s2d[i].kinem.alpha) + cs*sin(surf.s2d[i].kinem.alpha)
        cd[i] = cn*sin(surf.s2d[i].kinem.alpha) - cs*cos(surf.s2d[i].kinem.alpha)

        #Pitching moment is clockwise or nose up positive
        cm[i] = cn*surf.s2d[i].pvt - 2*pi*((surf.s2d[i].kinem.u*cos(surf.s2d[i].kinem.alpha)/surf.s2d[i].uref + surf.s2d[i].kinem.hdot*sin(surf.s2d[i].kinem.alpha)/surf.s2d[i].uref)*(surf.s2d[i].a0[1]/4. + surf.s2d[i].aterm[1]/4. - surf.s2d[i].aterm[2]/8.) + (surf.s2d[i].c/surf.s2d[i].uref)*(7. *surf.s2d[i].a0dot[1]/16. + 3. *surf.s2d[i].adot[1]/16. + surf.s2d[i].adot[2]/16. - surf.s2d[i].adot[3]/64. )) - nonl_m
    end

    cl3d = 0
    cd3d = 0
    cm3d = 0
    
    for is = 1:surf.nspan
        if is == 1
            yl = -0.5*surf.AR*surf.cref
            yu = 0.5*(surf.yle[is] + surf.yle[is+1])
        elseif is == surf.nspan
            yl = 0.5*(surf.yle[is] + surf.yle[is-1])
            yu = 0. 
        else
            yl = 0.5*(surf.yle[is] + surf.yle[is-1])
            yu = 0.5*(surf.yle[is] + surf.yle[is+1])
        end
        cl3d += cl[is]*(yu-yl) 
        cd3d += cd[is]*(yu-yl)
        cm3d += cm[is]*(yu-yl)
    end
    cl3d = cl3d/(0.5*surf.AR*surf.cref)
    cd3d = cd3d/(0.5*surf.AR*surf.cref)
    cm3d = cm3d/(0.5*surf.AR*surf.cref)
    
    return cl3d, cd3d, cm3d, cl, cd, cm
end
