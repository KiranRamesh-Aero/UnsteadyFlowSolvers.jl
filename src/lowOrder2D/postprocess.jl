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
    nonl = nonl*2./(surf.uref*surf.uref*surf.c)
    nonl_m = nonl_m*2./(surf.uref*surf.uref*surf.c*surf.c)

    # Normal force coefficient
    cn = cnc + cnnc + nonl

    # Lift and drag coefficients
    cl = cn*cos(surf.kinem.alpha) + cs*sin(surf.kinem.alpha)
    cd = cn*sin(surf.kinem.alpha)-cs*cos(surf.kinem.alpha)

    #Pitching moment is clockwise or nose up positive
    cm = cn*surf.pvt - 2*pi*((surf.kinem.u*cos(surf.kinem.alpha)/surf.uref + surf.kinem.hdot*sin(surf.kinem.alpha)/surf.uref)*(surf.a0[1]/4. + surf.aterm[1]/4. - surf.aterm[2]/8.) + (surf.c/surf.uref)*(7.*surf.a0dot[1]/16. + 3.*surf.adot[1]/16. + surf.adot[2]/16. - surf.adot[3]/64.)) - nonl_m
    return cl, cd, cm
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

    cnc1 = (4./(surf.uref*surf.uref*surf.c))*(int_lxtx + int_lxwsx + int_txwax
    + int_wsxwax + int_lztz + int_lzwsz + int_tzwaz + int_wszwaz)

    int_wax_kinem = simpleTrapz(wa_x.*(surf.kinem.u*cos(surf.kinem.alpha)
    + surf.kinem.hdot*sin(surf.kinem.alpha) - surf.kinem.alphadot*surf.cam), surf.x)
    int_lx_etac = simpleTrapz(l_x*surf.kinem.alphadot.*surf.cam, surf.theta)

    cnc2 = 2*pi*(surf.kinem.u*cos(surf.kinem.alpha)/surf.uref
    + surf.kinem.hdot*sin(surf.kinem.alpha)/surf.uref)*(surf.a0[1] + surf.aterm[1]/2.)
    + (4./(surf.uref*surf.uref*surf.c))*(int_wax_kinem - int_lx_etac)

    int_tx_etat = simpleTrapz((t_x + wa_x)*surf.kinem.alphadot.*surf.thick, surf.x)

    int_wsz_kinem = simpleTrapz(ws_z.*(surf.kinem.u*sin(surf.kinem.alpha)/surf.uref
    - surf.kinem.hdot*cos(surf.kinem.alpha) + surf.kinem.alphadot*(surf.x - surf.pvt*surf.c)), surf.x)

    cnc3 = pi*surf.bterm[1]*(surf.kinem.u*sin(surf.kinem.alpha)/surf.uref -
    surf.kinem.hdot*cos(surf.kinem.alpha)/surf.uref
    -surf.kinem.alphadot*surf.c*(0.5 - surf.pvt)/surf.uref)
    - surf.kinem.alphadot*surf.c*pi*surf.bterm[2]/(4*surf.uref)
    + (4./(surf.uref*surf.uref*surf.c))*(int_wsz_kinem - int_tx_etat)

    int_wax = zeros(surf.ndiv)
    d_int_wax = zeros(surf.ndiv)

    for i = 1:surf.ndiv
        int_wax[i] = simpleTrapz(wa_x[1:i], surf.x[1:i])

        d_int_wax[i] = (int_wax[i] - int_wax_prev[i])/dt
    end

    cnnc = (2*pi*surf.c/surf.uref)*(3*surf.a0dot[1]/4 + surf.adot[1]/4 + surf.adot[2]/8)
    + (4./(surf.uref*surf.uref*surf.c))*simpleTrapz(d_int_wax, surf.x)

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


function calc_forces(surf::TwoDSurf2DOF)

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
    nonl = nonl*2./(surf.uref*surf.uref*surf.c)
    nonl_m = nonl_m*2./(surf.uref*surf.uref*surf.c*surf.c)

    # Normal force coefficient
    cn = cnc + cnnc + nonl

    # Lift and drag coefficients
    cl = cn*cos(surf.kinem.alpha) + cs*sin(surf.kinem.alpha)
    cd = cn*sin(surf.kinem.alpha)-cs*cos(surf.kinem.alpha)

    #Pitching moment is clockwise or nose up positive
    cm = cn*surf.pvt - 2*pi*((surf.kinem.u*cos(surf.kinem.alpha)/surf.uref + surf.kinem.hdot*sin(surf.kinem.alpha)/surf.uref)*(surf.a0[1]/4. + surf.aterm[1]/4. - surf.aterm[2]/8.) + (surf.c/surf.uref)*(7.*surf.a0dot[1]/16. + 3.*surf.adot[1]/16. + surf.adot[2]/16. - surf.adot[3]/64.)) - nonl_m
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
    write(f, ["#time \t", "alpha (deg) \t", "h/c \t", "u/uref \t"])
    writedlm(f, [t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u]')
    close(f)

    f = open("FourierCoeffs", "w")
    write(f, ["#Fourier coeffs (0-n) \t", "d/dt (Fourier coeffs) \n"])
    matfour = zeros(surf.naterm+1, 2)
    matfour[:,1] = [surf.a0[1];surf.aterm[:]]
    matfour[1:4,2] = [surf.a0dot[1];surf.adot[:]]
    writedlm(f, matfour)
    close(f)

    # f = open("forces", "w")
    # write(f, ["#cl \t", "cd \t", "cm \t", "Gamma \t", "cn \t", "cs \t", "cnc \t",
    #           "cnnc \t", "nonl \t", "cm_n \t", "cm_pvt \t", "nonl_m \n"])
    # cl, cd, cm, gamma, cn, cs, cnc, cncc, nonl, cm_n, cm_pvt, nonl_m = calc_forces_more(surf)
    # writedlm(f, [cl, cd, cm, gamma, cn, cs, cnc, cnnc, nonl, cm_n, cm_pvt, nonl_m])
    # close(f)

    f = open("tev", "w")
    write(f, ["#strength \t", "x-position \t", "z-position \n"])
    tevmat = zeros(length(curfield.tev), 3)
    for i = 1:length(curfield.tev)
        tevmat[i,:] = [curfield.tev[i].s curfield.tev[i].x curfield.tev[i].z]
    end
    writedlm(f, tevmat)
    close(f)

    f = open("lev", "w")
    write(f, ["#strength \t", "x-position \t", "z-position \n"])
    levmat = zeros(length(curfield.lev), 3)
    for i = 1:length(curfield.lev)
        levmat[i,:] = [curfield.lev[i].s curfield.lev[i].x curfield.lev[i].z]
    end
    writedlm(f, levmat)
    close(f)

    f = open("boundv", "w")
    write(f, ["#strength \t", "x-position \t", "z-position \n"])
    bvmat = zeros(length(surf.bv), 3)
    for i = 1:length(surf.bv)
        bvmat[i,:] = [surf.bv[i].s surf.bv[i].x surf.bv[i].z]
    end
    writedlm(f, bvmat)
    close(f)
    cd("..")
end

function writeStamp(dirname::String, t::Float64, surf::TwoDSurf2DOF, curfield::TwoDFlowField)
    dirvec = readdir()
    if dirname in dirvec
        rm(dirname, recursive=true)
    end
    mkdir(dirname)
    cd(dirname)

    f = open("timeKinem", "w")
    write(f, ["#time \t", "alpha (deg) \t", "h/c \t", "u/uref \t"])
    writedlm(f, [t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u]')
    close(f)

    f = open("FourierCoeffs", "w")
    write(f, ["#Fourier coeffs (0-n) \t", "d/dt (Fourier coeffs) \n"])
    matfour = zeros(surf.naterm+1, 2)
    matfour[:,1] = [surf.a0[1];surf.aterm[:]]
    matfour[1:4,2] = [surf.a0dot[1];surf.adot[:]]
    writedlm(f, matfour)
    close(f)

    # f = open("forces", "w")
    # write(f, ["#cl \t", "cd \t", "cm \t", "Gamma \t", "cn \t", "cs \t", "cnc \t",
    #           "cnnc \t", "nonl \t", "cm_n \t", "cm_pvt \t", "nonl_m \n"])
    # cl, cd, cm, gamma, cn, cs, cnc, cncc, nonl, cm_n, cm_pvt, nonl_m = calc_forces_more(surf)
    # writedlm(f, [cl, cd, cm, gamma, cn, cs, cnc, cnnc, nonl, cm_n, cm_pvt, nonl_m])
    # close(f)

    f = open("tev", "w")
    write(f, ["#strength \t", "x-position \t", "z-position \n"])
    tevmat = zeros(length(curfield.tev), 3)
    for i = 1:length(curfield.tev)
        tevmat[i,:] = [curfield.tev[i].s curfield.tev[i].x curfield.tev[i].z]
    end
    writedlm(f, tevmat)
    close(f)

    f = open("lev", "w")
    write(f, ["#strength \t", "x-position \t", "z-position \n"])
    levmat = zeros(length(curfield.lev), 3)
    for i = 1:length(curfield.lev)
        levmat[i,:] = [curfield.lev[i].s curfield.lev[i].x curfield.lev[i].z]
    end
    writedlm(f, levmat)
    close(f)

    f = open("boundv", "w")
    write(f, ["#strength \t", "x-position \t", "z-position \n"])
    bvmat = zeros(length(surf.bv), 3)
    for i = 1:length(surf.bv)
        bvmat[i,:] = [surf.bv[i].s surf.bv[i].x surf.bv[i].z]
    end
    writedlm(f, bvmat)
    close(f)
    cd("..")
end
