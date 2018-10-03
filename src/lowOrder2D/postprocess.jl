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
    Serialization.Serialization.serialize(f, ["#time \t", "alpha (deg) \t", "h/c \t", "u/uref \t"])
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

function writeStamp(dirname::String, t::Float64, surf::TwoDSurf2DOF, curfield::TwoDFlowField)
    dirvec = readdir()
    if dirname in dirvec
        rm(dirname, recursive=true)
    end
    mkdir(dirname)
    cd(dirname)

    f = open("timeKinem", "w")
    write(f, ["#time \t", "alpha (deg) \t", "h/c \t", "u/uref \t"])
    DelimitedFiles.writedlm(f, [t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u]')
    close(f)

    f = open("FourierCoeffs", "w")
    write(f, ["#Fourier coeffs (0-n) \t", "d/dt (Fourier coeffs) \n"])
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
    write(f, ["#strength \t", "x-position \t", "z-position \n"])
    tevmat = zeros(length(curfield.tev), 3)
    for i = 1:length(curfield.tev)
        tevmat[i,:] = [curfield.tev[i].s curfield.tev[i].x curfield.tev[i].z]
    end
    DelimitedFiles.writedlm(f, tevmat)
    close(f)

    f = open("lev", "w")
    write(f, ["#strength \t", "x-position \t", "z-position \n"])
    levmat = zeros(length(curfield.lev), 3)
    for i = 1:length(curfield.lev)
        levmat[i,:] = [curfield.lev[i].s curfield.lev[i].x curfield.lev[i].z]
    end
    DelimitedFiles.writedlm(f, levmat)
    close(f)

    f = open("boundv", "w")
    write(f, ["#strength \t", "x-position \t", "z-position \n"])
    bvmat = zeros(length(surf.bv), 3)
    for i = 1:length(surf.bv)
        bvmat[i,:] = [surf.bv[i].s surf.bv[i].x surf.bv[i].z]
    end
    DelimitedFiles.writedlm(f, bvmat)
    close(f)
    cd("..")
end
