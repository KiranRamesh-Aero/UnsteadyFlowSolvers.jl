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
    DelimitedFiles.writedlm(f, [t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u]')
    close(f)

    f = open("FourierCoeffs", "w")
    Serialization.serialize(f, ["#Fourier coeffs (0-n) \t", "d/dt (Fourier coeffs) \n"])
    matfour = zeros(surf.naterm+1, 2)
    matfour[:,1] = [surf.a0[1];surf.aterm[:]]
    matfour[:,2] = [surf.a0dot[1];surf.adot[:]]
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

    f = open("delcp_edgevel", "w")
    Serialization.serialize(f, ["#x \t", "delcp \t", "delcp_inner \t", "delcp_outer \t", "qu \t", "ql \n"])
    mat = zeros(length(surf.x), 6)
    mat[:,1] = surf.x[:]
    delcp, delcp_in, delcp_out = UNSflow.calc_delcp(surf, [curfield.u[1]; curfield.w[1]])
    qu, ql = UNSflow.calc_edgeVel(surf, [curfield.u[1]; curfield.w[1]])
    mat[:,2] = delcp[:]
    mat[:,3] = delcp_in[:]
    mat[:,4] = delcp_out[:]
    mat[:,5] = qu[:]
    mat[:,6] = ql[:]
    DelimitedFiles.writedlm(f, mat)
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


function calc_forces(surf::TwoDSurfThick, int_wax_prev::Vector{Float64}, int_c_prev::Vector{Float64}, int_t_prev::Vector{Float64}, dt :: Float64)
    
    l_x = zeros(surf.ndiv)
    t_z = zeros(surf.ndiv)
    l_z = zeros(surf.ndiv)
    t_x = zeros(surf.ndiv)
    ws_x = zeros(surf.ndiv)
    ws_z = zeros(surf.ndiv)
    wa_x = zeros(surf.ndiv)
    wa_z = zeros(surf.ndiv)

    #All terms involving l_x, l_z (having the singularity) are integrated in theta to eliminate the singularity

    for i = 2:surf.ndiv
        l_x[i] = surf.a0[1]*(1+cos(surf.theta[i]))/sin(surf.theta[i])
        l_z[i] = -surf.a0[1]
        
        ws_x[i] = 0.5*(surf.uind_u[i] + surf.uind_l[i])
        wa_x[i] = 0.5*(surf.uind_u[i] - surf.uind_l[i])

        ws_z[i] = 0.5*(surf.wind_u[i] - surf.wind_l[i])
        wa_z[i] = 0.5*(surf.wind_u[i] + surf.wind_l[i])

        for ia = 1:surf.naterm
            l_x[i] += surf.aterm[ia]*sin(ia*surf.theta[i])
            t_x[i] -= surf.bterm[ia]*cos(ia*surf.theta[i])
            l_z[i] += surf.aterm[ia]*cos(ia*surf.theta[i])
            t_z[i] += surf.bterm[ia]*sin(ia*surf.theta[i])
        end
        l_x[i] = l_x[i]*surf.uref
        t_x[i] = t_x[i]*surf.uref
        l_z[i] = l_z[i]*surf.uref
        t_z[i] = t_z[i]*surf.uref
    end

    nonl1 = 2*simpleTrapz((ws_x .- surf.kinem.alphadot*surf.cam .+ t_x).*l_x, surf.x)
    nonl2 = 2*simpleTrapz(wa_x.*(surf.kinem.u*cos(surf.kinem.alpha) + surf.kinem.hdot*sin(surf.kinem.alpha) .- surf.kinem.alphadot*surf.cam .+ t_x .+ ws_x), surf.x)

    nonl3 = 2*simpleTrapz((wa_x .- l_z).*t_z, surf.x)

    nonl4 = 2*simpleTrapz(ws_z.*(surf.kinem.u*sin(surf.kinem.alpha) - surf.kinem.hdot*cos(surf.kinem.alpha) .+ surf.kinem.alphadot.*(surf.x .- surf.pvt*surf.c) .+ wa_z .+ l_z), surf.x)

    cnc_nonl = (nonl1 + nonl2 + nonl3 + nonl4)/(0.5*surf.uref^2*surf.c)

    cnc = 2*pi*(surf.kinem.u*cos(surf.kinem.alpha)/surf.uref + surf.kinem.hdot*sin(surf.kinem.alpha)/surf.uref)*(surf.a0[1] + 0.5*surf.aterm[1]) + pi*surf.bterm[1]*(surf.kinem.u*sin(surf.kinem.alpha)/surf.uref - surf.kinem.hdot*cos(surf.kinem.alpha)/surf.uref + surf.kinem.alphadot*surf.c*(0.5-surf.pvt)/surf.uref) - pi*surf.bterm[2]*surf.kinem.alphadot*surf.c/(4*surf.uref) + cnc_nonl
    

    int_wax = zeros(surf.ndiv)
    d_int_wax = zeros(surf.ndiv)
    int_c = zeros(surf.ndiv)
    d_int_c = zeros(surf.ndiv)
    int_t = zeros(surf.ndiv)
    d_int_t = zeros(surf.ndiv)
    
    for i = 1:surf.ndiv
        int_wax[i] = simpleTrapz(wa_x[1:i], surf.x[1:i])
        int_c[i] = simpleTrapz(surf.cam_slope[1:i].*(ws_z[1:i] + t_z[1:i]), surf.x[1:i])
        int_t[i] = simpleTrapz(surf.thick_slope[1:i].*(surf.kinem.u*sin(surf.kinem.alpha) - surf.kinem.hdot*cos(surf.kinem.alpha) .+ surf.kinem.alphadot.*(surf.x[1:i] .- surf.pvt*surf.c) .+ wa_x[1:i] .+ l_z[1:i]), surf.x[1:i])
        
        d_int_wax[i] = (int_wax[i] - int_wax_prev[i])/dt
        d_int_c[i] = (int_c[i] = int_c_prev[i])/dt
        d_int_t[i] = (int_t[i] = int_t_prev[i])/dt
    end

    cnnc = (2. *pi*surf.c/surf.uref)*(3. *surf.a0dot[1]/4. + surf.adot[1]/4. + surf.adot[2]/8.) +
        (4. /(surf.uref*surf.uref*surf.c))*simpleTrapz(d_int_wax, surf.x) + 
        (4. /(surf.uref*surf.uref*surf.c))*simpleTrapz(d_int_c, surf.x) +
        (4. /(surf.uref*surf.uref*surf.c))*simpleTrapz(d_int_t, surf.x)
    
    # Suction force given in eqn (2.31) Ramesh et al.
    #cs = 2*pi*surf.a0[1]*surf.a0[1]
    cs = 4*pi*surf.a0[1]*surf.a0[1]/surf.rho[1]
    
    # Normal force coefficient
    cn = cnc + cnnc

    # Lift and drag coefficients
    cl = cn*cos(surf.kinem.alpha) + cs*sin(surf.kinem.alpha)
    cd = cn*sin(surf.kinem.alpha)-cs*cos(surf.kinem.alpha)

    #Pitching moment is clockwise or nose up positive

    
    # int_lxtx = simpleTrapz(l_x.*t_x, surf.theta)
    # int_lxwsx = simpleTrapz(l_x.*ws_x, surf.theta)
    # int_txwax = simpleTrapz(t_x.*wa_x, surf.x)
    # int_wsxwax = simpleTrapz(ws_x.*wa_x, surf.x)

    # int_lztz = simpleTrapz(l_z.*t_z, surf.theta)
    # int_lzwsz = simpleTrapz(l_z.*ws_z, surf.theta)
    # int_tzwaz = simpleTrapz(t_z.*wa_z, surf.x)
    # int_wszwaz = simpleTrapz(ws_z.*wa_z, surf.x)




    
    # cnc1 = (4. /(surf.uref*surf.uref*surf.c))*(int_lxtx + int_lxwsx + int_txwax +
    #                                            int_wsxwax + int_lztz + int_lzwsz + int_tzwaz + int_wszwaz)
    
    # int_wax_kinem = simpleTrapz(wa_x.*(surf.kinem.u*cos(surf.kinem.alpha) +
    #                                    surf.kinem.hdot*sin(surf.kinem.alpha) .- surf.kinem.alphadot*surf.cam), surf.x)
    # int_lx_etac = simpleTrapz(l_x*surf.kinem.alphadot.*surf.cam, surf.theta)
    
    # cnc2 = 2*pi*(surf.kinem.u*cos(surf.kinem.alpha)/surf.uref +
    #              surf.kinem.hdot*sin(surf.kinem.alpha)/surf.uref)*(surf.a0[1] + surf.aterm[1]/2.) +
    # (4. /(surf.uref*surf.uref*surf.c))*(int_wax_kinem - int_lx_etac)
    
    # int_tx_etat = simpleTrapz((t_x + wa_x)*surf.kinem.alphadot.*surf.thick, surf.x)
    
    # int_wsz_kinem = simpleTrapz(ws_z.*(surf.kinem.u*sin(surf.kinem.alpha)/surf.uref -
    #                                    surf.kinem.hdot*cos(surf.kinem.alpha) .+ surf.kinem.alphadot*(surf.x .- surf.pvt*surf.c)), surf.x)
    
    # cnc3 = pi*surf.bterm[1]*(surf.kinem.u*sin(surf.kinem.alpha)/surf.uref -
    #                          surf.kinem.hdot*cos(surf.kinem.alpha)/surf.uref -
    #                          surf.kinem.alphadot*surf.c*(0.5 - surf.pvt)/surf.uref) -
    #                          surf.kinem.alphadot*surf.c*pi*surf.bterm[2]/(4*surf.uref) +
    # (4. /(surf.uref*surf.uref*surf.c))*(int_wsz_kinem - int_tx_etat)
    
    # int_wax = zeros(surf.ndiv)
    # d_int_wax = zeros(surf.ndiv)
    
    # for i = 1:surf.ndiv
    #     int_wax[i] = simpleTrapz(wa_x[1:i], surf.x[1:i])
        
    #     d_int_wax[i] = (int_wax[i] - int_wax_prev[i])/dt
    # end
    
    # cnnc = (2. *pi*surf.c/surf.uref)*(3. *surf.a0dot[1]/4. + surf.adot[1]/4. + surf.adot[2]/8.) +
    # (4. /(surf.uref*surf.uref*surf.c))*simpleTrapz(d_int_wax, surf.x)
    
    # # Suction force given in eqn (2.31) Ramesh et al.
    # cs = 2*pi*surf.a0[1]*surf.a0[1]
    
    # # Normal force coefficient
    # cn = cnc1 + cnc2 + cnc3 + cnnc

    # # Lift and drag coefficients
    # cl = cn*cos(surf.kinem.alpha) + cs*sin(surf.kinem.alpha)
    # cd = cn*sin(surf.kinem.alpha)-cs*cos(surf.kinem.alpha)

    # #Pitching moment is clockwise or nose up positive

return cnc, cnnc, cn, cs, cl, cd, int_wax, int_c, int_t
end

function calc_delcp(surf::TwoDSurf, vels::Vector{Float64})

    p_in = zeros(surf.ndiv)
    p_out = zeros(surf.ndiv)
    gam = zeros(surf.ndiv)
    gamint = zeros(surf.ndiv)
    p_com = zeros(surf.ndiv)
    gammod = zeros(surf.ndiv)
      
    for i = 1:surf.ndiv
        udash = surf.kinem.u*cos(surf.kinem.alpha) + surf.kinem.hdot*sin(surf.kinem.alpha) + surf.uind[i]*cos(surf.kinem.alpha) - surf.wind[i]*sin(surf.kinem.alpha)
        
        p_in[i] = 8*surf.uref*sqrt(surf.c)*surf.a0[1]*sqrt(surf.x[i])*udash/(surf.rho*surf.c + 2*surf.x[i]) + 4*surf.uref*sqrt(surf.c)*surf.a0dot[1]*sqrt(surf.x[i])
    end
    for i = 2:surf.ndiv
        udash = surf.kinem.u*cos(surf.kinem.alpha) + surf.kinem.hdot*sin(surf.kinem.alpha) + surf.uind[i]*cos(surf.kinem.alpha) - surf.wind[i]*sin(surf.kinem.alpha)
        gam[i] = surf.a0[1]*cot(surf.theta[i]/2)
        gamint[i] = surf.a0dot[1]*(surf.theta[i] + sin(surf.theta[i]))  
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
        
        p_out[i] = 2*udash*gam[i] + gamint[i] 
        p_com[i] = 2*udash*surf.uref*surf.a0[1]*(cos(surf.theta[i]/2)-1)/sin(surf.theta[i]/2) + 2*udash*gammod[i] + gamint[i] + 8*surf.uref*udash*surf.c*surf.a0[1]*sin(surf.theta[i]/2)/(surf.rho*surf.c + 2*surf.c*(sin(surf.theta[i]/2))^2)
    end
    
    return p_com, p_in, p_out
end

function calc_delcp(surf::TwoDSurfThick, vels::Vector{Float64})

    p_in = zeros(surf.ndiv)
    p_out = zeros(surf.ndiv)
    gam = zeros(surf.ndiv)
    gamint = zeros(surf.ndiv)
    p_com = zeros(surf.ndiv)
    gammod = zeros(surf.ndiv)

    srcval = zeros(surf.ndiv)
    src_z = zeros(surf.ndiv)
    
    for i = 2:surf.ndiv

        wtx = 0.5*(surf.uind_u[i]*cos(surf.kinem.alpha) - surf.wind_u[i]*sin(surf.kinem.alpha) +
                   surf.uind_l[i]*cos(surf.kinem.alpha) - surf.wind_l[i]*sin(surf.kinem.alpha))
        
        for n = 1:surf.naterm
            srcval[i] -= surf.bterm[n]*cos(n*surf.theta[i])
            src_z[i] += surf.bterm[n]*sin(n*surf.theta[i])
        end
        srcval[i] = srcval[i]*surf.uref
        src_z[i] = src_z[i]*surf.uref
        
        udash = (surf.kinem.u + vels[1])*cos(surf.kinem.alpha) + (surf.kinem.hdot - vels[2])*sin(surf.kinem.alpha) + wtx + srcval[i]
        
        p_in[i] = 8*surf.uref*sqrt(surf.c)*surf.a0[1]*sqrt(surf.x[i])*udash/(surf.rho*surf.c + 2*surf.x[i]) + 4*surf.uref*sqrt(surf.c)*surf.a0dot[1]*sqrt(surf.x[i])
        
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
        
        wlz = 0.5*(surf.wind_u[i]*cos(surf.kinem.alpha) + surf.uind_u[i]*sin(surf.kinem.alpha) +
                   surf.wind_l[i]*cos(surf.kinem.alpha) + surf.uind_l[i]*sin(surf.kinem.alpha))
        
        wtz = 0.5*(surf.wind_u[i]*cos(surf.kinem.alpha) + surf.uind_u[i]*sin(surf.kinem.alpha) -
                   surf.wind_l[i]*cos(surf.kinem.alpha) - surf.uind_l[i]*sin(surf.kinem.alpha))
        
        wtx = 0.5*(surf.uind_u[i]*cos(surf.kinem.alpha) - surf.wind_u[i]*sin(surf.kinem.alpha) +
                   surf.uind_l[i]*cos(surf.kinem.alpha) - surf.wind_l[i]*sin(surf.kinem.alpha))

        wlx = 0.5*(surf.uind_u[i]*cos(surf.kinem.alpha) - surf.wind_u[i]*sin(surf.kinem.alpha) -
                   surf.uind_l[i]*cos(surf.kinem.alpha) + surf.wind_l[i]*sin(surf.kinem.alpha))

        t1 = 2*wlx*udash
        t2 = 2*((surf.kinem.u + vels[1])*sin(surf.kinem.alpha) - (surf.kinem.hdot - vels[2])*cos(surf.kinem.alpha) + surf.kinem.alphadot*(surf.x[i] - surf.pvt*surf.c) + wlz + gam[i])*(src_z[i] + wtz) 
        
        #Omit additional app mass terms
        p_out[i] = 2*udash*gam[i] + gamint[i] + t1 + t2

        p_com[i] = 2*udash*surf.uref*surf.a0[1]*(cos(surf.theta[i]/2)-1)/sin(surf.theta[i]/2) + 2*udash*gammod[i] + gamint[i] + t1 + t2 + 8*surf.uref*udash*surf.c*surf.a0[1]*sin(surf.theta[i]/2)/(surf.rho*surf.c + 2*surf.c*(sin(surf.theta[i]/2))^2)
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

function calc_edgeVel(surf::TwoDSurfThick, vels::Vector{Float64})

    gammod = zeros(surf.ndiv)
    q_com_u = zeros(surf.ndiv)
    q_com_l = zeros(surf.ndiv)
    srcval = zeros(surf.ndiv)

    for n = 1:surf.naterm
        srcval[1] -= surf.bterm[n]*cos(n*surf.theta[1])/surf.uref
    end
    q_com_u[1] = surf.uref*surf.a0[1]/sqrt(0.5*surf.rho)
    q_com_l[1] = surf.uref*surf.a0[1]/sqrt(0.5*surf.rho)
    
    for i = 2:surf.ndiv

        wtx = 0.5*(surf.uind_u[i]*cos(surf.kinem.alpha) - surf.wind_u[i]*sin(surf.kinem.alpha) +
                   surf.uind_l[i]*cos(surf.kinem.alpha) - surf.wind_l[i]*sin(surf.kinem.alpha))

        for n = 1:surf.naterm
            gammod[i] += surf.aterm[n]*sin(n*surf.theta[i])
            srcval[i] -= surf.bterm[n]*cos(n*surf.theta[i])
        end
        gammod[i] = gammod[i]*surf.uref
        srcval[i] = srcval[i]*surf.uref
        
        udash = (surf.kinem.u + vels[1])*cos(surf.kinem.alpha) + (surf.kinem.hdot - vels[2])*sin(surf.kinem.alpha) + wtx + srcval[i]
             
        q_com_u[i] = surf.uref*surf.a0[1]*(cos(surf.theta[i]/2) - 1)/sin(surf.theta[i]/2) + surf.uref*gammod[i] + (sqrt(surf.c)*surf.uref*surf.a0[1] + sqrt(surf.x[i])*udash)/sqrt(surf.x[i] + surf.rho*surf.c/2)
        q_com_l[i] = -surf.uref*surf.a0[1]*(cos(surf.theta[i]/2) - 1)/sin(surf.theta[i]/2) - surf.uref*gammod[i] + (-sqrt(surf.c)*surf.uref*surf.a0[1] + sqrt(surf.x[i])*udash)/sqrt(surf.x[i] + surf.rho*surf.c/2)
    end

    return q_com_u, q_com_l
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
    mat[:,2] = [surf.a0dot[1];surf.adot[:]]
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
