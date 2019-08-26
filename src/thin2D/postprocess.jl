function calc_forces(surf::TwoDSurf, vels::Vector{Float64})

    # First term in eqn (2.30) Ramesh et al. in coefficient form
    cnc = 2*pi*((surf.kinem.u + vels[1])*cos(surf.kinem.alpha)/surf.uref + (surf.kinem.hdot - vels[2])*sin(surf.kinem.alpha)/surf.uref)*(surf.a0[1] + surf.aterm[1]/2.)

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
    cm = cn*surf.pvt - 2*pi*(((surf.kinem.u + vels[1])*cos(surf.kinem.alpha)/surf.uref + (surf.kinem.hdot - vels[2])*sin(surf.kinem.alpha)/surf.uref)*(surf.a0[1]/4. + surf.aterm[1]/4. - surf.aterm[2]/8.) + (surf.c/surf.uref)*(7. *surf.a0dot[1]/16. + 3. *surf.adot[1]/16. + surf.adot[2]/16. - surf.adot[3]/64.)) - nonl_m
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
    writedlm(f, [t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u]')
    close(f)

    f = open("FourierCoeffs", "w")
    Serialization.serialize(f, ["#Fourier coeffs (0-n) \t", "d/dt (Fourier coeffs) \n"])
    matfour = zeros(surf.naterm+1, 2)
    matfour[:,1] = [surf.a0[1];surf.aterm[:]]
    matfour[:,2] = [surf.a0dot[1];surf.adot[:]]
    writedlm(f, matfour)
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
    writedlm(f, tevmat)
    close(f)

    f = open("lev", "w")
    Serialization.serialize(f, ["#strength \t", "x-position \t", "z-position \n"])
    levmat = zeros(length(curfield.lev), 3)
    for i = 1:length(curfield.lev)
        levmat[i,:] = [curfield.lev[i].s curfield.lev[i].x curfield.lev[i].z]
    end
    writedlm(f, levmat)
    close(f)

    f = open("boundv", "w")
    Serialization.serialize(f, ["#strength \t", "x-position \t", "z-position \n"])
    bvmat = zeros(length(surf.bv), 3)
    for i = 1:length(surf.bv)
        bvmat[i,:] = [surf.bv[i].s surf.bv[i].x surf.bv[i].z]
    end
    writedlm(f, bvmat)
    close(f)

    f = open("delcp_edgevel", "w")
    Serialization.serialize(f, ["#x \t", "delcp \t", "delcp_inner \t", "delcp_outer \t", "qu \t", "ql \n"])
    mat = zeros(length(surf.x), 6)
    mat[:,1] = surf.x[:]
    delcp, delcp_in, delcp_out = calc_delcp(surf, [curfield.u[1]; curfield.w[1]])
    qu, ql = calc_edgeVel(surf, [curfield.u[1]; curfield.w[1]])
    mat[:,2] = delcp[:]
    mat[:,3] = delcp_in[:]
    mat[:,4] = delcp_out[:]
    mat[:,5] = qu[:]
    mat[:,6] = ql[:]
    writedlm(f, mat)
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
        matfour[:,2] = [surf[i].a0dot[1];surf[i].adot[:]]
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

function calc_delcp(surf::TwoDSurf, vels::Vector{Float64})

    p_in = zeros(surf.ndiv)
    p_out = zeros(surf.ndiv)
    gam = zeros(surf.ndiv)
    gamint = zeros(surf.ndiv)
    p_com = zeros(surf.ndiv)
   
      
    for i = 1:surf.ndiv
        udash = surf.kinem.u*cos(surf.kinem.alpha) + surf.kinem.hdot*sin(surf.kinem.alpha) + surf.uind[i]*cos(surf.kinem.alpha) - surf.wind[i]*sin(surf.kinem.alpha)
        
        p_in[i] = 2*(4*surf.uref*sqrt(surf.c)*surf.a0[1]*sqrt(surf.x[i])*udash/(surf.rho*surf.c + 2*surf.x[i]) + 4*surf.uref*sqrt(surf.c)*surf.a0dot[1]*sqrt(surf.x[i]))
    end
    for i = 2:surf.ndiv
        udash = surf.kinem.u*cos(surf.kinem.alpha) + surf.kinem.hdot*sin(surf.kinem.alpha) + surf.uind[i]*cos(surf.kinem.alpha) - surf.wind[i]*sin(surf.kinem.alpha)
        gam[i] = surf.a0[1]*cot(surf.theta[i]/2)
        gamint[i] = surf.a0dot[1]*(surf.theta[i] + sin(surf.theta[i]))  
        gammod = zeros(surf.ndiv)
        
        for n = 1:surf.naterm
            gam[i] += surf.aterm[n]*sin(n*surf.theta[i])
            gammod[i] += surf.aterm[n]*sin(n*surf.theta[i])
            if n == 1
                gamint[i] += surf.adot[1]*(surf.theta[i]/2 - sin(2*surf.theta[i])/4)
            else
                gamint[i] += surf.adot[n]/2*(sin((n-1)*surf.theta[i])/(n-1) - sin((n+1)*surf.theta[i])/(n+1))
            end
        end
        gam[i] = gam[i]*surf.uref
        gammod[i] = gammod[i]*surf.uref
        gamint[i] = gamint[i]*surf.uref*surf.c
        
        p_out[i] = 2*(2*udash*gam[i] + gamint[i]) 
        p_com[i] = 2*(2*udash*surf.uref*surf.a0[1]*(cos(surf.theta[i]/2)-1)/sin(surf.theta[i]/2) + 2*udash*gammod[i] + gamint[i] + 4*surf.uref*udash*surf.c*surf.a0[1]*sin(surf.theta[i]/2)/(surf.rho*surf.c + 2*surf.c*(sin(surf.theta[i]/2))^2))
    end
    
    return p_com, p_in, p_out
end

function calc_edgeVel(surf::TwoDSurf, vels::Vector{Float64})

    gammod = zeros(surf.ndiv)
    q_com_u = zeros(surf.ndiv)
    q_com_l = zeros(surf.ndiv)

    q_com_u[1] = surf.uref*surf.a0[1]/sqrt(0.5*surf.rho)
    q_com_l[1] = surf.uref*surf.a0[1]/sqrt(0.5*surf.rho)
    
    for i = 2:surf.ndiv
        udash = (surf.kinem.u + vels[1])*cos(surf.kinem.alpha) + (surf.kinem.hdot - vels[2])*sin(surf.kinem.alpha) + surf.uind[i]*cos(surf.kinem.alpha) - surf.wind[i]*sin(surf.kinem.alpha)
        
        for n = 1:surf.naterm
            gammod[i] += surf.aterm[n]*sin(n*surf.theta[i])
        end
        gammod[i] = gammod[i]*surf.uref
                
        q_com_u[i] = surf.uref*surf.a0[1]*(cos(surf.theta[i]/2) - 1)/sin(surf.theta[i]/2) + surf.uref*gammod[i] + (sqrt(surf.c)*surf.uref*surf.a0[1] + sqrt(surf.x[i])*udash)/sqrt(surf.x[i] + surf.rho*surf.c/2)
        q_com_l[i] = -surf.uref*surf.a0[1]*(cos(surf.theta[i]/2) - 1)/sin(surf.theta[i]/2) - surf.uref*gammod[i] + (-sqrt(surf.c)*surf.uref*surf.a0[1] + sqrt(surf.x[i])*udash)/sqrt(surf.x[i] + surf.rho*surf.c/2)
    end

    return q_com_u, q_com_l
end

function calc_stag(surf :: TwoDSurf)
    stag = (surf.kinem.u/(surf.kinem.u*cos(surf.kinem.alpha) + surf.kinem.hdot*sin(surf.kinem.alpha)
                   + surf.uind[1]*cos(surf.kinem.alpha) - surf.wind[1]*sin(surf.kinem.alpha))*surf.a0[1])^2
    return stag
end
