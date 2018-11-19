function calc_forces(surf::TwoDSurf)

    # First term in eqn (2.30) Ramesh et al. in coefficient form
    cnc = 2*pi*(surf.kinem.u*cos(surf.kinem.alpha)/surf.uref + surf.kinem.hdot*sin(surf.kinem.alpha)/surf.uref)*(surf.a0[1] + surf.aterm[1]/2.)

    # Second term in eqn (2.30) Ramesh et al. in coefficient form
    cnnc = 2*pi*(3*surf.c*surf.a0dot[1]/(4*surf.uref) + surf.c*surf.adot[1]/(4*surf.uref) + surf.c*surf.adot[2]/(8*surf.uref))

    # Suction force given in eqn (2.31) Ramesh et al.
    cs = 2*pi*surf.a0[1]*surf.a0[1]

    #The components of normal force and moment from induced velocities are calulcated in dimensional units and nondimensionalized later
    nonl=0
    nonl_m=0
    for ib = 1:surf.ndiv-1
        nonl = nonl + (surf.uind[ib]*cos(surf.kinem.alpha) - surf.wind[ib]*sin(surf.kinem.alpha))*surf.bv[ib].s
        nonl_m = nonl_m + (surf.uind[ib]*cos(surf.kinem.alpha) - surf.wind[ib]*sin(surf.kinem.alpha))*surf.x[ib]*surf.bv[ib].s
    end
    nonl = nonl*2. /(surf.uref*surf.uref*surf.c)
    nonl_m = nonl_m*2. /(surf.uref*surf.uref*surf.c*surf.c)

    # Normal force coefficient
    cn = cnc + cnnc + nonl

    # Lift and drag coefficients
    cl = cn*cos(surf.kinem.alpha) + cs*sin(surf.kinem.alpha)
    cd = cn*sin(surf.kinem.alpha)-cs*cos(surf.kinem.alpha)

    #Pitching moment is clockwise or nose up positive
    cm = cn*surf.pvt - 2*pi*((surf.kinem.u*cos(surf.kinem.alpha)/surf.uref + surf.kinem.hdot*sin(surf.kinem.alpha)/surf.uref)*(surf.a0[1]/4. + surf.aterm[1]/4. - surf.aterm[2]/8.) + (surf.c/surf.uref)*(7. *surf.a0dot[1]/16. + 3. *surf.adot[1]/16. + surf.adot[2]/16. - surf.adot[3]/64.)) - nonl_m
    return cl, cd, cm
end

function writeStamp(dirname::String, t::Float64, surf::TwoDSurf, curfield::TwoDFlowField)
    dirvec = readdir()
    if dirname in dirvec
        rm(dirname, recursive=true)
    end
    mkdir(dirname)
    cd(dirname)

    f = open("timeKinem", "w")
    Serialization.serialize(f, ["#time \t", "alpha (deg) \t", "h/c \t", "u/uref \t"])
    DelimitedFiles.writedlm(f, [t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u]')
    close(f)

    f = open("FourierCoeffs", "w")
    Serialization.serialize(f, ["#Fourier coeffs (0-n) \t", "d/dt (Fourier coeffs) \n"])
    matfour = zeros(surf.naterm+1, 2)
    matfour[:,1] = [surf.a0[1];surf.aterm[:]]
    matfour[1:4,2] = [surf.a0dot[1];surf.adot[:]]
    DelimitedFiles.writedlm(f, matfour)
    close(f)

    # f = open("forces", "w")
    # write(f, ["#cl \t", "cd \t", "cm \t", "Gamma \t", "cn \t", "cs \t", "cnc \t",
    #           "cnnc \t", "nonl \t", "cm_n \t", "cm_pvt \t", "nonl_m \n"])
    # cl, cd, cm, gamma, cn, cs, cnc, cncc, nonl, cm_n, cm_pvt, nonl_m = calc_forces_more(surf)
    # writedlm(f, [cl, cd, cm, gamma, cn, cs, cnc, cnnc, nonl, cm_n, cm_pvt, nonl_m])
    # close(f)

    f = open("tev", "w")
    Serialization.serialize(f, ["#strength \t", "x-position \t", "z-position \n"])
    tevmat = zeros(length(curfield.tev), 3)
    for i = 1:length(curfield.tev)
        tevmat[i,:] = [curfield.tev[i].s curfield.tev[i].x curfield.tev[i].z]
    end
    DelimitedFiles.writedlm(f, tevmat)
    close(f)

    f = open("lev", "w")
    Serialization.serialize(f, ["#strength \t", "x-position \t", "z-position \n"])
    levmat = zeros(length(curfield.lev), 3)
    for i = 1:length(curfield.lev)
        levmat[i,:] = [curfield.lev[i].s curfield.lev[i].x curfield.lev[i].z]
    end
    DelimitedFiles.writedlm(f, levmat)
    close(f)

    f = open("boundv", "w")
    Serialization.serialize(f, ["#strength \t", "x-position \t", "z-position \n"])
    bvmat = zeros(length(surf.bv), 3)
    for i = 1:length(surf.bv)
        bvmat[i,:] = [surf.bv[i].s surf.bv[i].x surf.bv[i].z]
    end
    DelimitedFiles.writedlm(f, bvmat)
    close(f)
    cd("..")
end

function writeStamp(dirname::String, t::Float64, surf::Vector{TwoDSurf}, curfield::TwoDFlowField)
    dirvec = readdir()
    if dirname in dirvec
        rm(dirname, recursive=true)
    end
    mkdir(dirname)
    cd(dirname)

    nsurf = length(surf)
    for i = 1:nsurf
        f = open("timeKinem-$i", "w")
        Serialization.Serialization.serialize(f, ["#time \t", "alpha (deg) \t", "h/c \t", "u/uref \t"])
        DelimitedFiles.writedlm(f, [t, surf[i].kinem.alpha, surf[i].kinem.h, surf[i].kinem.u]')
        close(f)

        f = open("FourierCoeffs-$i", "w")
        Serialization.serialize(f, ["#Fourier coeffs (0-n) \t", "d/dt (Fourier coeffs) \n"])
        matfour = zeros(surf[i].naterm+1, 2)
        matfour[:,1] = [surf[i].a0[1];surf[i].aterm[:]]
        matfour[1:4,2] = [surf[i].a0dot[1];surf[i].adot[:]]
        DelimitedFiles.writedlm(f, matfour)
        close(f)

        # f = open("forces", "w")
        # write(f, ["#cl \t", "cd \t", "cm \t", "Gamma \t", "cn \t", "cs \t", "cnc \t",
        #           "cnnc \t", "nonl \t", "cm_n \t", "cm_pvt \t", "nonl_m \n"])
        # cl, cd, cm, gamma, cn, cs, cnc, cncc, nonl, cm_n, cm_pvt, nonl_m = calc_forces_more(surf)
        # writedlm(f, [cl, cd, cm, gamma, cn, cs, cnc, cnnc, nonl, cm_n, cm_pvt, nonl_m])
        # close(f)
    end

    f = open("tev", "w")
    Serialization.serialize(f, ["#strength \t", "x-position \t", "z-position \n"])
    tevmat = zeros(length(curfield.tev), 3)
    for i = 1:length(curfield.tev)
        tevmat[i,:] = [curfield.tev[i].s curfield.tev[i].x curfield.tev[i].z]
    end
    DelimitedFiles.writedlm(f, tevmat)
    close(f)

    f = open("lev", "w")
    Serialization.serialize(f, ["#strength \t", "x-position \t", "z-position \n"])
    levmat = zeros(length(curfield.lev), 3)
    levcount = 0
    for i = 1:length(curfield.lev)
        if curfield.lev[i].vc != 0.
            levcount += 1
            levmat[levcount,:] = [curfield.lev[i].s curfield.lev[i].x curfield.lev[i].z]
        end
    end

    levmat = levmat[1:levcount,:]
    DelimitedFiles.writedlm(f, levmat)
    close(f)

    for is = 1:nsurf
        f = open("boundv-$is", "w")
        Serialization.serialize(f, ["#strength \t", "x-position \t", "z-position \n"])
        bvmat = zeros(length(surf[is].bv), 3)
        for i = 1:length(surf[is].bv)
            bvmat[i,:] = [surf[is].bv[i].s surf[is].bv[i].x surf[is].bv[i].z]
        end
        DelimitedFiles.writedlm(f, bvmat)
        close(f)
    end

    cd("..")
end


function calc_forces(surf::TwoDSurfThick, int_wax_prev :: Vector{Float64}, dt :: Float64)


    l_x = zeros(surf.ndiv)
    t_z = zeros(surf.ndiv)
    l_z = zeros(surf.ndiv)
    t_x = zeros(surf.ndiv)
    ws_x = zeros(surf.ndiv)
    ws_z = zeros(surf.ndiv)
    wa_x = zeros(surf.ndiv)
    wa_z = zeros(surf.ndiv)

    #All terms involving l_x, l_z (having the singularity) are integrated in theta to eliminate the singularity

    for i = 1:surf.ndiv
        l_x[i] = surf.a0[1]*(1+cos(surf.theta[i]))
        l_z[i] = -surf.a0[1]*sin(surf.theta[i])

        ws_x[i] = 0.5*(surf.uind_u[i] + surf.uind_l[i])
        wa_x[i] = 0.5*(surf.uind_u[i] - surf.uind_l[i])

        ws_z[i] = 0.5*(surf.wind_u[i] - surf.wind_l[i])
        wa_z[i] = 0.5*(surf.wind_u[i] + surf.wind_l[i])

        for ia = 1:surf.naterm
            l_x[i] += surf.aterm[ia]*sin(ia*surf.theta[i])*sin(surf.theta[i])
            t_x[i] -= surf.bterm[ia]*cos(ia*surf.theta[i])
            l_z[i] += surf.aterm[ia]*cos(ia*surf.theta[i])*sin(surf.theta[i])
            t_z[i] += surf.bterm[ia]*sin(ia*surf.theta[i])
        end
        l_x[i] = l_x[i]*surf.uref*surf.c/2.
        t_x[i] = t_x[i]*surf.uref
        l_z[i] = l_z[i]*surf.uref*surf.c/2.
        t_z[i] = t_z[i]*surf.uref
    end

    int_lxtx = simpleTrapz(l_x.*t_x, surf.theta)
    int_lxwsx = simpleTrapz(l_x.*ws_x, surf.theta)
    int_txwax = simpleTrapz(t_x.*wa_x, surf.x)
    int_wsxwax = simpleTrapz(ws_x.*wa_x, surf.x)

    int_lztz = simpleTrapz(l_z.*t_z, surf.theta)
    int_lzwsz = simpleTrapz(l_z.*ws_z, surf.theta)
    int_tzwaz = simpleTrapz(t_z.*wa_z, surf.x)
    int_wszwaz = simpleTrapz(ws_z.*wa_z, surf.x)

    cnc1 = (4. /(surf.uref*surf.uref*surf.c))*(int_lxtx + int_lxwsx + int_txwax +
    int_wsxwax + int_lztz + int_lzwsz + int_tzwaz + int_wszwaz)

    int_wax_kinem = simpleTrapz(wa_x.*(surf.kinem.u*cos(surf.kinem.alpha) +
    surf.kinem.hdot*sin(surf.kinem.alpha) .- surf.kinem.alphadot*surf.cam), surf.x)
    int_lx_etac = simpleTrapz(l_x*surf.kinem.alphadot.*surf.cam, surf.theta)

    cnc2 = 2*pi*(surf.kinem.u*cos(surf.kinem.alpha)/surf.uref +
    surf.kinem.hdot*sin(surf.kinem.alpha)/surf.uref)*(surf.a0[1] + surf.aterm[1]/2.) +
    (4. /(surf.uref*surf.uref*surf.c))*(int_wax_kinem - int_lx_etac)

    int_tx_etat = simpleTrapz((t_x + wa_x)*surf.kinem.alphadot.*surf.thick, surf.x)

    int_wsz_kinem = simpleTrapz(ws_z.*(surf.kinem.u*sin(surf.kinem.alpha)/surf.uref -
    surf.kinem.hdot*cos(surf.kinem.alpha) .+ surf.kinem.alphadot*(surf.x .- surf.pvt*surf.c)), surf.x)

    cnc3 = pi*surf.bterm[1]*(surf.kinem.u*sin(surf.kinem.alpha)/surf.uref -
    surf.kinem.hdot*cos(surf.kinem.alpha)/surf.uref -
    surf.kinem.alphadot*surf.c*(0.5 - surf.pvt)/surf.uref) -
    surf.kinem.alphadot*surf.c*pi*surf.bterm[2]/(4*surf.uref) +
    (4. /(surf.uref*surf.uref*surf.c))*(int_wsz_kinem - int_tx_etat)

    int_wax = zeros(surf.ndiv)
    d_int_wax = zeros(surf.ndiv)

    for i = 1:surf.ndiv
        int_wax[i] = simpleTrapz(wa_x[1:i], surf.x[1:i])

        d_int_wax[i] = (int_wax[i] - int_wax_prev[i])/dt
    end

    cnnc = (2. *pi*surf.c/surf.uref)*(3. *surf.a0dot[1]/4. + surf.adot[1]/4. + surf.adot[2]/8.) +
    (4. /(surf.uref*surf.uref*surf.c))*simpleTrapz(d_int_wax, surf.x)

    # Suction force given in eqn (2.31) Ramesh et al.
    cs = 2*pi*surf.a0[1]*surf.a0[1]

    # Normal force coefficient
    cn = cnc1 + cnc2 + cnc3 + cnnc

    # Lift and drag coefficients
    cl = cn*cos(surf.kinem.alpha) + cs*sin(surf.kinem.alpha)
    cd = cn*sin(surf.kinem.alpha)-cs*cos(surf.kinem.alpha)

    #Pitching moment is clockwise or nose up positive

    return cnc1, cnc2, cnc3, cnnc, cn, cs, cl, cd
end

function calc_q_cp(surf::TwoDSurfThick)

    u1c = zeros(surf.ndiv)
    u1t = zeros(surf.ndiv)
    w1t = zeros(surf.ndiv)
    w1c = zeros(surf.ndiv)
    for ib = 2:surf.ndiv-1
        u1c[ib] = surf.a0[1]*(1 + cos(surf.theta[ib]))/sin(surf.theta[ib])
        w1t[ib] = surf.b0[1]*(1 + cos(surf.theta[ib]))/sin(surf.theta[ib])
        w1t[ib] += surf.b0[2]*(1 - cos(surf.theta[ib]))/sin(surf.theta[ib])
        w1c[ib] = -surf.a0[1]
        u1t[ib] = surf.b0[1] - surf.b0[2]
        for ia = 1:surf.naterm
            u1c[ib] += surf.aterm[ia]*sin(ia*surf.theta[ib])
            u1t[ib] -= surf.bterm[ia]*cos(ia*surf.theta[ib])
            w1t[ib] += surf.bterm[ia]*sin(ia*surf.theta[ib])
            w1c[ib] += surf.aterm[ia]*cos(ia*surf.theta[ib])
        end
        u1c[ib] = u1c[ib]*surf.uref
        u1t[ib] = u1t[ib]*surf.uref
        w1t[ib] = w1t[ib]*surf.uref
        w1c[ib] = w1c[ib]*surf.uref
    end

    w1t[1] = 2*w1t[2] - w1t[3]
    w1t[surf.ndiv] = 2*w1t[surf.ndiv-1] - w1t[surf.ndiv-2]
    u1t[1] = 2*u1t[2] - u1t[3]
    u1t[surf.ndiv] = 2*u1t[surf.ndiv-1] - u1t[surf.ndiv-2]
    w1c[1] = 2*w1c[2] - w1c[3]
    w1c[surf.ndiv] = 2*w1c[surf.ndiv-1] - w1c[surf.ndiv-2]
    u1c[1] = 2*u1c[2] - u1c[3]
    u1c[surf.ndiv] = 2*u1c[surf.ndiv-1] - u1c[surf.ndiv-2]

    lamb = mean(surf.cam_slope[1:5])
    rho = surf.rho
    qu = zeros(surf.ndiv)
    ql = zeros(surf.ndiv)
    quc = zeros(surf.ndiv)
    qlc = zeros(surf.ndiv)

    for i = 1:surf.ndiv
        u_u = u1c[i] + u1t[i] + surf.uind_u[i] + surf.kinem.u*cos(surf.kinem.alpha) +
        surf.kinem.hdot*sin(surf.kinem.alpha) - surf.kinem.alphadot*(surf.cam[i] + surf.thick[i])

        w_u = w1c[i] + w1t[i] + surf.wind_u[i] + surf.kinem.u*sin(surf.kinem.alpha) -
        surf.kinem.hdot*cos(surf.kinem.alpha) + surf.kinem.alphadot*(surf.x[i] - surf.pvt*surf.c)

        vu[i] = sqrt(u_u^2 + w_u^2)

        u_l = -u1c[i] + u1t[i] + surf.uind_l[i] + surf.kinem.u*cos(surf.kinem.alpha) +
        surf.kinem.hdot*sin(surf.kinem.alpha) - surf.kinem.alphadot*(surf.cam[i] - surf.thick[i])

        w_l = w1c[i] - w1t[i] + surf.wind_l[i] + surf.kinem.u*sin(surf.kinem.alpha) -
        surf.kinem.hdot*cos(surf.kinem.alpha) + surf.kinem.alphadot*(surf.x[i] - surf.pvt*surf.c)

        vl[i] = sqrt(u_l^2 + w_l^2)

        vu_van[i] = sqrt((surf.x[i]+lamb*sqrt(2*rho*surf.x[i]))/(surf.x[i] + lamb*sqrt(2*rho*surf.x[i]) + rho/2.))*(vu[i]/surf.uref+rho/(4*surf.x[i]))
        vl_van[i] = sqrt((surf.x[i]-lamb*sqrt(2*rho*surf.x[i]))/(surf.x[i] - lamb*sqrt(2*rho*surf.x[i]) + rho/2.))*(vl[i]/surf.uref+rho/(4*surf.x[i]))

    end
    return vu, vl, vu_van, vl_van
end

function writeStamp(dirname::String, t::Float64, surf::TwoDSurfThick, curfield::TwoDFlowField)
    dirvec = readdir()
    if dirname in dirvec
        rm(dirname, recursive=true)
    end
    mkdir(dirname)
    cd(dirname)

    f = open("timeKinem", "w")
    Serialization.serialize(f, ["#time \t", "alpha (deg) \t", "h/c \t", "u/uref \t"])
    DelimitedFiles.writedlm(f, [t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u]')
    close(f)

    f = open("FourierCoeffsA", "w")
    Serialization.serialize(f, ["#Fourier coeffs (0-n) \t", "d/dt (Fourier coeffs) \n"])
    mat = zeros(surf.naterm+1, 2)
    mat[:,1] = [surf.a0[1];surf.aterm[:]]
    mat[1:4,2] = [surf.a0dot[1];surf.adot[:]]
    DelimitedFiles.writedlm(f, mat)
    close(f)

    f = open("FourierCoeffsB", "w")
    Serialization.serialize(f, ["#Fourier coeffs (0-n) \t", "d/dt (Fourier coeffs) \n"])
    mat = zeros(surf.naterm, 1)
    mat[:,1] = surf.bterm[:]
    DelimitedFiles.writedlm(f, mat)
    close(f)

    # f = open("forces", "w")
    # write(f, ["#cl \t", "cd \t", "cm \t", "Gamma \t", "cn \t", "cs \t", "cnc \t",
    #           "cnnc \t", "nonl \t", "cm_n \t", "cm_pvt \t", "nonl_m \n"])
    # cl, cd, cm, gamma, cn, cs, cnc, cncc, nonl, cm_n, cm_pvt, nonl_m = calc_forces_more(surf)
    # writedlm(f, [cl, cd, cm, gamma, cn, cs, cnc, cnnc, nonl, cm_n, cm_pvt, nonl_m])
    # close(f)

    f = open("tev", "w")
    Serialization.serialize(f, ["#strength \t", "x-position \t", "z-position \n"])
    tevmat = zeros(length(curfield.tev), 3)
    for i = 1:length(curfield.tev)
        tevmat[i,:] = [curfield.tev[i].s curfield.tev[i].x curfield.tev[i].z]
    end
    DelimitedFiles.writedlm(f, tevmat)
    close(f)

    f = open("lev", "w")
    Serialization.serialize(f, ["#strength \t", "x-position \t", "z-position \n"])
    levmat = zeros(length(curfield.lev), 3)
    for i = 1:length(curfield.lev)
        levmat[i,:] = [curfield.lev[i].s curfield.lev[i].x curfield.lev[i].z]
    end
    DelimitedFiles.writedlm(f, levmat)
    close(f)

    f = open("boundv", "w")
    Serialization.serialize(f, ["#strength \t", "x-position \t", "z-position \n"])
    bvmat = zeros(length(surf.bv), 3)
    for i = 1:length(surf.bv)
        bvmat[i,:] = [surf.bv[i].s surf.bv[i].x surf.bv[i].z]
    end
    DelimitedFiles.writedlm(f, bvmat)
    close(f)

    f = open("src", "w")
    Serialization.serialize(f, ["#strength \t", "x-position \t", "z-position \n"])
    bvmat = zeros(length(surf.src), 3)
    for i = 1:length(surf.src)
        bvmat[i,:] = [surf.src[i].s surf.src[i].x surf.src[i].z]
    end
    DelimitedFiles.writedlm(f, bvmat)
    close(f)
    
    cd("..")
end

function getEndCycle(mat, k)

    T = pi/k
    
    end_cycle = mat[end,1]
    
    #Find number of cycles
    ncyc = 0
    for i = 1:1000
        t_l = real(i)*T
        if t_l > end_cycle
            break
        end
        ncyc = ncyc + 1
    end
    
    start_t = real(ncyc-1)*T
    end_t = real(ncyc)*T
    start_ind = argmin(abs.(mat[:,1] .- start_t))
    end_ind = argmin(abs.(mat[:,1] .- end_t))

    nlast = end_ind - start_ind + 1
    
    newmat = zeros(nlast, 8)
    newmat[:,1] = (mat[start_ind:end_ind,1] .- start_t)/T
    for i = 2:7
        newmat[:,i] = mat[start_ind:end_ind,i]
    end

    f = open("resultsSummaryEndCycle", "w")
    Serialization.serialize(f, ["#t/T \t", "alpha (deg) \t", "h/c \t", "u/uref \t", "A0 \t", "Cl \t", "Cd \t", "Cm \n"])
    DelimitedFiles.writedlm(f, newmat)
    close(f)

    return newmat
end
