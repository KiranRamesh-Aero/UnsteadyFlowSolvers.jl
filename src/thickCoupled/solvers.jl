function IBL_shape_attached(Re, surf::TwoDSurfThick, curfield::TwoDFlowField, nsteps::Int64 = 300, dtstar::Float64 = 0.015, startflag = 0, writeflag = 0, writeInterval = 1000., delvort = delNone(); maxwrite = 50, nround=6)
    
    # If a restart directory is provided, read in the simulation data
    if startflag == 0
        mat = zeros(0, 13)
        t = 0.
    elseif startflag == 1
        dirvec = readdir()
        dirresults = map(x->(v = tryparse(Float64,x); isnull(v) ? 0.0 : get(v)),dirvec)
        latestTime = maximum(dirresults)
        mat = readdlm("resultsSummary")
        t = mat[end,1]
    else
        throw("invalid start flag, should be 0 or 1")
    end
    mat = mat'

    dt = dtstar*surf.c/surf.uref

    # if writeflag is on, determine the timesteps to write at
    if writeflag == 1
        writeArray = Int64[]
        tTot = nsteps*dt
        for i = 1:maxwrite
            tcur = writeInterval*real(i)
            if t > tTot
                break
            else
                push!(writeArray, Int(round(tcur/dt)))
            end
        end
    end

    vcore = 1.3*dt*surf.c

    int_wax = zeros(surf.ndiv)
    int_c = zeros(surf.ndiv)
    int_t = zeros(surf.ndiv)

    #At this time coded for symmetric situations only
    qu = zeros(surf.ndiv)
    ql = zeros(surf.ndiv)

    spos = zeros(2*surf.ndiv-1)
    sc = zeros(2*surf.ndiv-2)
    dsc = zeros(2*surf.ndiv-2)
    qc = zeros(2*surf.ndiv-2)
    qc_prev = zeros(2*surf.ndiv-2)
    qcx = zeros(2*surf.ndiv-2)
    qct = zeros(2*surf.ndiv-2)
    su = zeros(surf.ndiv)
    
    thick_orig = zeros(surf.ndiv)
    thick_orig[:] = surf.thick[:]
    thick_slope_orig = zeros(surf.ndiv)
    thick_slope_orig[:] = surf.thick_slope[:]

    cam_orig = zeros(surf.ndiv)
    cam_orig[:] = surf.cam[:]
    cam_slope_orig = zeros(surf.ndiv)
    cam_slope_orig[:] = surf.cam_slope[:]

    
    #Initialise boundary layer
    del, _, E, _ = initDelE(2*surf.ndiv-2)

    spos[1] = 0.
    dsdx = zeros(surf.ndiv)
    for i = 2:surf.ndiv
        dsdx[i] = sqrt(1 + (surf.cam_slope[i] - surf.thick_slope[i])^2)
    end
    dsdx[1] = dsdx[2]
    su[1] = 0.
    for i = surf.ndiv-1:-1:1
        su[surf.ndiv+1-i] = simpleTrapz(dsdx[surf.ndiv:-1:i], surf.x[surf.ndiv:-1:i])
    end
    spos[1:surf.ndiv] = -su[:]
    for i = 2:surf.ndiv
        dsdx[i] = sqrt(1 + (surf.cam_slope[i] + surf.thick_slope[i])^2)
    end
    dsdx[1] = dsdx[2]
    su[1] = 0.
    for i = 2:surf.ndiv
        su[i] = simpleTrapz(dsdx[1:i], surf.x[1:i])
    end
    spos[surf.ndiv+1:2*surf.ndiv-1] = spos[surf.ndiv] .+ su[2:surf.ndiv]
    for i = 1:2*surf.ndiv-2
        sc[i] = (spos[i] + spos[i+1])/2
    end

    dspos = diff(spos)
    
    #All BL parameters go from trailing edge on lower surface to leading edge to trailing edge on uper surface. 

    
    phi_u = zeros(surf.ndiv)
    phi_l = zeros(surf.ndiv)

    x_w = collect(surf.c*1.01:surf.c*0.01:surf.c*3)
    nw = length(x_w)
    wfn = zeros(nw)

    bound_circ = 0.
    
    tevstr = zeros(100)
    tevdist = zeros(100)
    restev = zeros(100,2)
    restev_prev = zeros(2,2)
    
    phi_u_temp = zeros(surf.ndiv)
    phi_l_temp = zeros(surf.ndiv)

    for istep = 1:nsteps

        #Udpate current time
        t = t + dt

        #Update kinematic parameters
        update_kinem(surf, t)

        #Update flow field parameters if any
        update_externalvel(curfield, t)

        #Update bound vortex positions
        update_boundpos(surf, dt)

        #Update incduced velocities on airfoil
        update_indbound(surf, curfield)

        ntev = length(curfield.tev)
        if ntev == 0
            xloc_tev = surf.bnd_x_chord[surf.ndiv] + 0.5*surf.kinem.u*dt*cos(surf.kinem.alpha)
            zloc_tev = surf.bnd_z_chord[surf.ndiv] - 0.5*surf.kinem.u*dt*sin(surf.kinem.alpha)
        else
            xloc_tev = surf.bnd_x_chord[surf.ndiv] + (1. /3.)*(curfield.tev[ntev].x - surf.bnd_x_chord[surf.ndiv])
            zloc_tev = surf.bnd_z_chord[surf.ndiv] + (1. /3.)*(curfield.tev[ntev].z - surf.bnd_z_chord[surf.ndiv])
        end
        
        #Outer iteration to solve coupled problem
        iter = 0
        res = 1
        del_iter = zeros(2*surf.ndiv-2)
        del_prev = zeros(2*surf.ndiv-2)
        soln = zeros(2*surf.naterm+1)
        E_iter = zeros(2*surf.ndiv-2)

        resTol = 1e-5
        iterMax = 20

        qc_prev[:] = qc[:]       
        
        #while (iter < iterMax)
        while res > resTol
            
            iter += 1

            function tev_iter!(FC, J, x)

                #x vector = [aterm; bterm; ate; tevstr]
                
                #FC = zeros(2*surf.naterm+2)
                #           J = zeros(2*surf.naterm+2, 2*surf.naterm+2)

                aterm = x[1:surf.naterm]
                bterm = x[surf.naterm+1:2*surf.naterm]
                ate = x[2*surf.naterm+1]
                tevstr = x[2*surf.naterm+2]
                
                dummyvort = TwoDVort(xloc_tev, zloc_tev, tevstr, vcore, 0., 0.)
                
                uu, wu = ind_vel([dummyvort], surf.bnd_x_u, surf.bnd_z_u)
                ul, wl = ind_vel([dummyvort], surf.bnd_x_l, surf.bnd_z_l)
                
                wlz = 0.5*((surf.wind_u .+ wu).*cos(surf.kinem.alpha) .+ (surf.uind_u .+ uu).*sin(surf.kinem.alpha) .+
                           (surf.wind_l .+ wl)*cos(surf.kinem.alpha) .+ (surf.uind_l .+ ul)*sin(surf.kinem.alpha))
                
                wtz = 0.5*((surf.wind_u .+ wu).*cos(surf.kinem.alpha) .+ (surf.uind_u .+ uu).*sin(surf.kinem.alpha) .-
                           (surf.wind_l .+ wl).*cos(surf.kinem.alpha) .- (surf.uind_l .+ ul).*sin(surf.kinem.alpha))
                
                wtx = 0.5*((surf.uind_u .+ uu).*cos(surf.kinem.alpha) .- (surf.wind_u .+ wu).*sin(surf.kinem.alpha) .+
                           (surf.uind_l .+ ul).*cos(surf.kinem.alpha) .- (surf.wind_l .+ wl).*sin(surf.kinem.alpha))
                
                wlx = 0.5*((surf.uind_u .+ uu).*cos(surf.kinem.alpha) .- (surf.wind_u .+ wu).*sin(surf.kinem.alpha) .-
                           (surf.uind_l .+ ul).*cos(surf.kinem.alpha) .+ (surf.wind_l .+ wl).*sin(surf.kinem.alpha))

                dummyvort = TwoDVort(xloc_tev, zloc_tev, 1., vcore, 0., 0.)
                uu, wu = ind_vel([dummyvort], surf.bnd_x_u, surf.bnd_z_u)
                ul, wl = ind_vel([dummyvort], surf.bnd_x_l, surf.bnd_z_l)

                
                wlz_t = 0.5*(wu.*cos(surf.kinem.alpha) .+ uu.*sin(surf.kinem.alpha) .+
                             wl*cos(surf.kinem.alpha) .+ ul*sin(surf.kinem.alpha))
                wtz_t = 0.5*(wu.*cos(surf.kinem.alpha) .+ uu.*sin(surf.kinem.alpha) .-
                             wl.*cos(surf.kinem.alpha) .- ul.*sin(surf.kinem.alpha))
                wtx_t = 0.5*(uu.*cos(surf.kinem.alpha) .- wu.*sin(surf.kinem.alpha) .+
                             ul.*cos(surf.kinem.alpha) .- wl.*sin(surf.kinem.alpha))
                wlx_t = 0.5*(uu.*cos(surf.kinem.alpha) .- wu.*sin(surf.kinem.alpha) .-
                             ul.*cos(surf.kinem.alpha) .+ wl.*sin(surf.kinem.alpha))
                
                rng = 1:surf.naterm

                Lx = zeros(surf.ndiv); Lz = zeros(surf.ndiv)
                Tx = zeros(surf.ndiv); Tz = zeros(surf.ndiv)
                phi_u_integ = zeros(surf.ndiv)
                phi_l_integ = zeros(surf.ndiv)
                
                i = surf.ndiv-1
                vref_x_u = (surf.kinem.u + curfield.u[1])*cos(surf.kinem.alpha) + (surf.kinem.hdot - curfield.w[1])*sin(surf.kinem.alpha) - surf.kinem.alphadot*(surf.cam[i] + surf.thick[i])
                vref_x_l = (surf.kinem.u + curfield.u[1])*cos(surf.kinem.alpha) + (surf.kinem.hdot - curfield.w[1])*sin(surf.kinem.alpha) - surf.kinem.alphadot*(surf.cam[i] - surf.thick[i])
                vref_z = (surf.kinem.u + curfield.u[1])*sin(surf.kinem.alpha) - (surf.kinem.hdot - curfield.w[1])*cos(surf.kinem.alpha) + surf.kinem.alphadot*(surf.x[i] - surf.pvt*surf.c)
                
                for i = 1:surf.ndiv
                    Lx[i] = surf.uref*(sum(aterm[rng]'*sin.(rng*surf.theta[i])) + ate*tan(surf.theta[i]/2))
                    Lz[i] = surf.uref*(sum(aterm[rng]'*cos.(rng*surf.theta[i])) + ate)
                    Tz[i] = surf.uref*sum(bterm[rng]'*sin.(rng*surf.theta[i]))
                    Tx[i] = -surf.uref*sum(bterm[rng]'*cos.(rng*surf.theta[i]))
                    if i ==1
                        phi_u_integ[i] = Lz[i] + Tz[i] + wtz[i] + wlz[i]
                        phi_l_integ[i] = -(Lz[i] - Tz[i] - wtz[i] +wlz[i] )
                    else
                        phi_u_integ[i] = (Lx[i] + Tx[i] + wtx[i] + wlx[i]) + (surf.cam_slope[i] + surf.thick_slope[i])*(Lz[i] + Tz[i] + wtz[i] + wlz[i])
                        phi_l_integ[i] = (-Lx[i] + Tx[i] + wtx[i] - wlx[i]) + (surf.cam_slope[i] - surf.thick_slope[i])*(Lz[i] - Tz[i] - wtz[i] + wlz[i])
                    end    
                end
                
                phi_u_bc = 0; phi_l_bc = 0
                for i = 2:surf.ndiv-1
                    phi_u_bc += 0.5*(phi_u_integ[i]/sqrt(1. + (surf.thick_slope[i] + surf.cam_slope[i])^2) + phi_u_integ[i-1]/sqrt(1. + (surf.thick_slope[i-1] + surf.cam_slope[i-1])^2))*sqrt((surf.x[i] - surf.x[i-1])^2 + (surf.cam[i] + surf.thick[i] - surf.cam[i-1] - surf.thick[i-1])^2)
                    phi_l_bc += 0.5*(phi_l_integ[i]/sqrt(1. + (-surf.thick_slope[i] + surf.cam_slope[i])^2) + phi_l_integ[i-1]/sqrt(1. + (-surf.thick_slope[i-1] + surf.cam_slope[i-1])^2))*sqrt((surf.x[i] - surf.x[i-1])^2 + (surf.cam[i] - surf.thick[i] - surf.cam[i-1] + surf.thick[i-1])^2)
                end

                if !(FC == nothing)
                    for i = 2:surf.ndiv-1
                        rhs_l = -(surf.kinem.u + curfield.u[1])*sin(surf.kinem.alpha) - surf.kinem.alphadot*(surf.x[i] - surf.pvt*surf.c) + (surf.kinem.hdot - curfield.w[1])*cos(surf.kinem.alpha) - wlz[i] + surf.cam_slope[i]*((surf.kinem.u + curfield.u[1])*cos(surf.kinem.alpha) + (surf.kinem.hdot - curfield.w[1])*sin(surf.kinem.alpha) + wtx[i] - surf.kinem.alphadot*surf.cam[i]) + surf.thick_slope[i]*(wlx[i] - surf.kinem.alphadot*surf.thick[i])
                        rhs_nonl = surf.cam_slope[i]*(wlx[i] - surf.kinem.alphadot*surf.thick[i]) + surf.thick_slope[i]*((surf.kinem.u + curfield.u[1])*cos(surf.kinem.alpha) + (surf.kinem.hdot - curfield.w[1])*sin(surf.kinem.alpha) + wtx[i] - surf.kinem.alphadot*surf.cam[i]) - wtz[i]
                        
                        FC[i-1] = Lz[i] - surf.cam_slope[i]*Tx[i] - surf.thick_slope[i]*Lx[i] - rhs_l
                        FC[surf.ndiv+i-3] = Tz[i] - surf.cam_slope[i]*Lx[i] - surf.thick_slope[i]*Tx[i] - rhs_nonl
                        
                    end
                    
                    #Kutta condition
                    i = surf.ndiv-1
                    qu = sqrt((vref_x_u + Lx[i] + Tx[i] + wlx[i] + wtx[i])^2 + (vref_z + Lz[i] + Tz[i] + wlz[i] + wtz[i])^2)
                    ql = sqrt((vref_x_l - Lx[i] + Tx[i] - wlx[i] + wtx[i])^2 + (vref_z + Lz[i] - Tz[i] + wlz[i] - wtz[i])^2)
                    FC[2*surf.ndiv-3] = tevstr - 0.5*(qu^2 - ql^2)*dt

                    #Kelvin condition
                    bc = phi_u_bc - phi_l_bc
                    FC[2*surf.ndiv-2] = bc - bound_circ + tevstr
                end

                if !(J == nothing)
                    for i = 2:surf.ndiv-1
                        for n = 1:surf.naterm
                            J[i-1,n] = cos(n*surf.theta[i]) - surf.thick_slope[i]*sin(n*surf.theta[i]) 
                            J[i-1,n+surf.naterm] = surf.cam_slope[i]*cos(n*surf.theta[i])
                            J[surf.ndiv+i-3,n] = -surf.cam_slope[i]*sin(n*surf.theta[i])
                            J[surf.ndiv+i-3,surf.naterm+n] = sin(n*surf.theta[i]) + surf.thick_slope[i]*cos(n*surf.theta[i])            
                        end
                        J[i-1,2*surf.naterm+1] = 1. - surf.thick_slope[i]*tan(surf.theta[i]/2)
                        J[surf.ndiv+i-3,2*surf.naterm+1] = -surf.cam_slope[i]*tan(surf.theta[i]/2)
                        J[i-1,2+2*surf.naterm] = -surf.cam_slope[i]*wtx_t[i] - surf.thick_slope[i]*wlx_t[i] + wlz_t[i]
                        J[surf.ndiv+i-3,2+2*surf.naterm] = -surf.cam_slope[i]*wlx_t[i] - surf.thick_slope[i]*wtx_t[i] + wtz_t[i]
                    end
                    
                    #Kutta condition
                    i = surf.ndiv-1
                    for n = 1:surf.naterm
                        J[2*surf.ndiv-3,n] = -2*dt*(0.5*(vref_x_u + vref_x_l) + Tx[i] + wtx[i])*sin(n*surf.theta[i]) - 2*dt*(Tz[i] + wtz[i])*cos(n*surf.theta[i])
                        J[2*surf.ndiv-3,n+surf.naterm] = 2*dt*(0.5*(vref_x_u - vref_x_l) + Lx[i] + wlx[i])*cos(n*surf.theta[i]) - 2*dt*(vref_z + Lz[i] + wlz[i])*sin(n*surf.theta[i])
                    end
                    J[2*surf.ndiv-3,2*surf.naterm+1] = -2*dt*(0.5*(vref_x_u + vref_x_l) + Tx[i] + wtx[i])*tan(surf.theta[i]/2) - 2*dt*(Tz[i] + wtz[i])
                    J[2*surf.ndiv-3,2*surf.naterm+2] = 1. - 2*dt*(0.5*(vref_x_u + vref_x_l) + Tx[i] + wtx[i])*wlx_t[i] - 2*dt*(0.5*(vref_x_u - vref_x_l) + Lx[i] + wlx[i])*wtx_t[i] - 2*dt*(vref_z + Lz[i] + wlz[i])*wtz_t[i] -2*dt*(Tz[i] + wtz[i])*wlz_t[i]
                    
                    
                    J[2*surf.ndiv-2,:] .= 0.
                    #Kelvin condition
                    i = 1
                    ate_u = tan(surf.theta[i]/2) + (surf.cam_slope[i] + surf.thick_slope[i])
                    ate_l = -tan(surf.theta[i]/2) + (surf.cam_slope[i] - surf.thick_slope[i])
                    tev_u = wtx_t[i] + wlx_t[i] + (surf.cam_slope[i] + surf.thick_slope[i])*(wlz_t[i] + wtz_t[i])
                    tev_l = wtx_t[i] - wlx_t[i] + (surf.cam_slope[i] - surf.thick_slope[i])*(wlz_t[i] - wtz_t[i])
                    den_u = sqrt(1. + (surf.cam_slope[i] + surf.thick_slope[i])^2)
                    den_l = sqrt(1. + (surf.cam_slope[i] - surf.thick_slope[i])^2)
                    ate_u_int = 0; ate_l_int = 0; tev_u_int = 0; tev_l_int = 0

                    for i = 2:surf.ndiv-1
                        
                        ate_u_p = ate_u
                        tev_u_p = tev_u
                        ate_l_p = ate_l
                        tev_l_p = tev_l
                        den_u_p = den_u
                        den_l_p = den_l
                        
                        ate_u = tan(surf.theta[i]/2) + (surf.cam_slope[i] + surf.thick_slope[i])
                        ate_l = -tan(surf.theta[i]/2) + (surf.cam_slope[i] - surf.thick_slope[i])
                        tev_u = wtx_t[i] + wlx_t[i] + (surf.cam_slope[i] + surf.thick_slope[i])*(wlz_t[i] + wtz_t[i])
                        tev_l = wtx_t[i] - wlx_t[i] + (surf.cam_slope[i] - surf.thick_slope[i])*(wlz_t[i] - wtz_t[i])
                        den_u = sqrt(1. + (surf.cam_slope[i] + surf.thick_slope[i])^2)
                        den_l = sqrt(1. + (surf.cam_slope[i] - surf.thick_slope[i])^2)
                        
                        ds_u = sqrt((surf.x[i] - surf.x[i-1])^2 + (surf.cam[i] + surf.thick[i] - surf.cam[i-1] - surf.thick[i-1])^2)
                        ds_l = sqrt((surf.x[i] - surf.x[i-1])^2 + (surf.cam[i] - surf.thick[i] - surf.cam[i-1] + surf.thick[i-1])^2)
                        J[2*surf.ndiv-2,2*surf.naterm+1] += 0.5*(ate_u/den_u + ate_u_p/den_u_p)*ds_u 
                        J[2*surf.ndiv-2,2*surf.naterm+1] -= 0.5*(ate_l/den_l + ate_l_p/den_l_p)*ds_l
                        
                        J[2*surf.ndiv-2,2*surf.naterm+2] += 0.5*(tev_u/den_u + tev_u_p/den_u_p)*ds_u
                        J[2*surf.ndiv-2,2*surf.naterm+2] -= 0.5*(tev_l/den_l + tev_l_p/den_l_p)*ds_l
                    end                
                    J[2*surf.ndiv-2,2*surf.naterm+2] += 1.
                    
                    for n = 1:surf.naterm
                        i = 1
                        an_u = sin(n*surf.theta[i]) + (surf.cam_slope[i] + surf.thick_slope[i])*cos(n*surf.theta[i])
                        an_l = -sin(n*surf.theta[i]) + (surf.cam_slope[i] - surf.thick_slope[i])*cos(n*surf.theta[i])
                        bn_u = -cos(n*surf.theta[i]) + (surf.cam_slope[i] + surf.thick_slope[i])*sin(n*surf.theta[i])
                        bn_l = -cos(n*surf.theta[i]) - (surf.cam_slope[i] - surf.thick_slope[i])*sin(n*surf.theta[i])
                        den_u = sqrt(1. + (surf.cam_slope[i] + surf.thick_slope[i])^2)
                        den_l = sqrt(1. + (surf.cam_slope[i] - surf.thick_slope[i])^2)
                        
                        an_u_int = 0; an_l_int = 0; bn_u_int = 0; bn_l_int = 0;
                        
                        for i = 2:surf.ndiv-1

                            an_u_p = an_u
                            bn_u_p = bn_u
                            an_l_p = an_l
                            bn_l_p = bn_l
                            den_u_p = den_u
                            den_l_p = den_l
                            
                            an_u = sin(n*surf.theta[i]) + (surf.cam_slope[i] + surf.thick_slope[i])*cos(n*surf.theta[i])
                            an_l = -sin(n*surf.theta[i]) + (surf.cam_slope[i] - surf.thick_slope[i])*cos(n*surf.theta[i])
                            bn_u = -cos(n*surf.theta[i]) + (surf.cam_slope[i] + surf.thick_slope[i])*sin(n*surf.theta[i])
                            bn_l = -cos(n*surf.theta[i]) - (surf.cam_slope[i] - surf.thick_slope[i])*sin(n*surf.theta[i])
                            den_u = sqrt(1. + (surf.cam_slope[i] + surf.thick_slope[i])^2)
                            den_l = sqrt(1. + (surf.cam_slope[i] - surf.thick_slope[i])^2)

                            ds_u = sqrt((surf.x[i] - surf.x[i-1])^2 + (surf.cam[i] + surf.thick[i] - surf.cam[i-1] - surf.thick[i-1])^2)
                            ds_l = sqrt((surf.x[i] - surf.x[i-1])^2 + (surf.cam[i] - surf.thick[i] - surf.cam[i-1] + surf.thick[i-1])^2)
                            J[2*surf.ndiv-2,n] += 0.5*(an_u/den_u + an_u_p/den_u_p)*ds_u
                            J[2*surf.ndiv-2,n+surf.naterm] += 0.5*(bn_u/den_u + bn_u_p/den_u_p)*ds_u
                            J[2*surf.ndiv-2,n] -= 0.5*(an_l/den_l + an_l_p/den_l_p)*ds_l
                            J[2*surf.ndiv-2,n+surf.naterm] -= 0.5*(bn_l/den_l + bn_l_p/den_l_p)*ds_l
                        end
                    end
                end
                
            end
                        
            xstart = [surf.aterm; surf.bterm; surf.ate[1]; -0.001]
            
            soln = nlsolve(only_fj!(tev_iter!), xstart, xtol = 1e-6, method = :newton)
            soln = soln.zero
            
            #assign the solution
            surf.aterm[:] = soln[1:surf.naterm]
            surf.bterm[:] = soln[surf.naterm+1:2*surf.naterm]
            surf.ate[1] = soln[2*surf.naterm+1]
            tevstr = soln[2*surf.naterm+2]
            temptev = TwoDVort(xloc_tev, zloc_tev, tevstr, vcore, 0., 0.)
            
            a_wfn = log10((qc[end]*del_iter[end] + qc[1]*del_iter[1])/sqrt(Re))
            wfn[1] = 1/sqrt(Re)*((qc[end]*del_iter[end] - qc[end-1]*del_iter[end-1])/(surf.x[end] - surf.x[end-1]) + (qc[1]*del_iter[1] - qc[2]*del_iter[2])/(surf.x[1] - surf.x[2]))
            for i = 2:nw
                wfn[i] = (10^(a_wfn - 3.2*(x_w[i] - surf.c)) - 10^(a_wfn - 3.2*(x_w[i-1] - surf.c)))/(x_w[i] - x_w[i-1])
            end
            
            #Source strengths in wake and induced velocity
            uind_src = zeros(surf.ndiv)
            for i = 1:surf.ndiv
                for iw = 1:nw-1
                    str = 0.5(wfn[iw] + wfn[iw+1])
                    xloc = 0.5*(x_w[iw] + x_w[iw+1])
                    uind_src[i] += 1/(2*pi)*str/(xloc - surf.x[i])*(x_w[iw+1] - x_w[iw])
                end
            end
            uind_src[:] .= 0.
            surf.uind_u[:] += uind_src[:]
            surf.uind_l[:] += uind_src[:]

            ind_new_u_u, ind_new_w_u = ind_vel([temptev], surf.bnd_x_u, surf.bnd_z_u)
            ind_new_u_l, ind_new_w_l = ind_vel([temptev], surf.bnd_x_l, surf.bnd_z_l)
            
            surf.uind_u[:] += ind_new_u_u[:]
            surf.wind_u[:] += ind_new_w_u[:]
            surf.uind_l[:] += ind_new_u_l[:]
            surf.wind_l[:] += ind_new_w_l[:]

            qu, ql = calc_edgeVel(surf, [curfield.u[1], curfield.w[1]])
            
            surf.uind_u[:] -= uind_src[:]
            surf.uind_l[:] -= uind_src[:]
            
            surf.uind_u[:] -= ind_new_u_u[:]
            surf.wind_u[:] -= ind_new_w_u[:]
            surf.uind_l[:] -= ind_new_u_l[:]
            surf.wind_l[:] -= ind_new_w_l[:]
            
            #smoothScaledEnd!(surf.x, qu, 10)

            #Full vector of qc with positive going from LE to TE
            for i = 1:surf.ndiv-1
                qc[i] = (ql[surf.ndiv-i] + ql[surf.ndiv+1-i])/2
                qc[surf.ndiv-1+i] = (qu[i] + qu[i+1])/2
            end
            
            qcx[1:surf.ndiv-1] = -diff(qc[1:surf.ndiv])./diff(sc[1:surf.ndiv])
            #qcx[1] = qcx[2]
            qcx[surf.ndiv+1:2*surf.ndiv-2] = diff(qc[surf.ndiv:2*surf.ndiv-2])./diff(sc[surf.ndiv:2*surf.ndiv-2])
            qcx[surf.ndiv] = qcx[surf.ndiv+1]
            #smoothEnd!(qucx, 10)
            
            if istep == 1
                qct[:] .= 0.
            else
                qct[:] = (qc[:] - qc_prev[:])/dt
            end

            qc_prev[:] = qc[:]

            
            stag = find_stag(surf, qu, ql)

            #Find the closest point to stagnation
            xs_ind = argmin(abs.(stag .- surf.x))

            #Split into two problems 
            if qu[1] > 0.
                #Stag on lower surface
                scu = [sc[surf.ndiv+1-xs_ind:surf.ndiv-1]; sc[surf.ndiv:2*surf.ndiv-2]] .- sc[surf.ndiv+1-xs_ind]
                scl = sc[surf.ndiv-xs_ind] .- reverse(sc[1:surf.ndiv-xs_ind]) 
                qcu = [-qc[surf.ndiv+1-xs_ind:surf.ndiv-1]; qc[surf.ndiv:2*surf.ndiv-2]]
                qcl = reverse(qc[1:surf.ndiv-xs_ind]) 
                qcux = [-qcx[surf.ndiv+1-xs_ind:surf.ndiv-1]; qcx[surf.ndiv:2*surf.ndiv-2]]
                qclx = reverse(qcx[1:surf.ndiv-xs_ind]) 
                qcut = [-qct[surf.ndiv+1-xs_ind:surf.ndiv-1]; qct[surf.ndiv:2*surf.ndiv-2]]
                qclt = reverse(qct[1:surf.ndiv-xs_ind]) 
                delu = [del[surf.ndiv+1-xs_ind:surf.ndiv-1]; del[surf.ndiv:2*surf.ndiv-2]]
                dell = reverse(del[1:surf.ndiv-xs_ind]) 
                Eu = [E[surf.ndiv+1-xs_ind:surf.ndiv-1]; E[surf.ndiv:2*surf.ndiv-2]]
                El = reverse(E[1:surf.ndiv-xs_ind])
                dsposu = [dspos[surf.ndiv+1-xs_ind:surf.ndiv-1]; dspos[surf.ndiv:2*surf.ndiv-2]]
                dsposl = reverse(dspos[1:surf.ndiv-xs_ind])
            else
                scu = sc[surf.ndiv-1+xs_ind:2*surf.ndiv-2] .- sc[surf.ndiv-1+xs_ind]
                scl = -([reverse(sc[surf.ndiv:surf.ndiv-2+xs_ind]); reverse(sc[1:surf.ndiv-1])] .- sc[surf.ndiv-2+xs_ind])
                qcu = qc[surf.ndiv-1+xs_ind:2*surf.ndiv-2]
                qcl = [-reverse(qc[surf.ndiv:surf.ndiv-2+xs_ind]); reverse(qc[1:surf.ndiv-1])]
                qcux = qcx[surf.ndiv-1+xs_ind:2*surf.ndiv-2]
                qclx = [-reverse(qcx[surf.ndiv:surf.ndiv-2+xs_ind]); reverse(qcx[1:surf.ndiv-1])]
                qcut = qct[surf.ndiv-1+xs_ind:2*surf.ndiv-2]
                qclt = [-reverse(qct[surf.ndiv:surf.ndiv-2+xs_ind]); reverse(qct[1:surf.ndiv-1])]
                delu = del[surf.ndiv-1+xs_ind:2*surf.ndiv-2]
                dell = [reverse(del[surf.ndiv:surf.ndiv-2+xs_ind]); reverse(del[1:surf.ndiv-1])]
                Eu = E[surf.ndiv-1+xs_ind:2*surf.ndiv-2]
                El = [reverse(E[surf.ndiv:surf.ndiv-2+xs_ind]); reverse(E[1:surf.ndiv-1])]
                dsposu = dspos[surf.ndiv-1+xs_ind:2*surf.ndiv-2]
                dsposl = [reverse(dspos[surf.ndiv:surf.ndiv-2+xs_ind]); reverse(dspos[1:surf.ndiv-1])]
            end

            qcux = reverse(smoothScaledEnd!(reverse(scu), reverse(qcux), 5))
            qclx = reverse(smoothScaledEnd!(reverse(scl), reverse(qclx), 5))
            smoothScaledEnd!(scu, qcux, 10)
            smoothScaledEnd!(scl, qclx, 10)
            
            
            #Solve the FV problem at cell centres
            plot(scu, qcux)
            plot(scl, qclx)
            
            wu = [delu delu.*(Eu .+ 1)]
            wl = [dell dell.*(El .+ 1)]
            
            wsolnu, i_sep = FVMIBLgridvar(wu, qcu, qcut, qcux, dsposu, t-dt, t)
            wsolnl, i_sep = FVMIBLgridvar(wl, qcl, qclt, qclx, dsposl, t-dt, t) 

            #Assign the solution
            delu = wsolnu[:,1]
            Eu = wsolnu[:,2]./wsolnu[:,1] .- 1.
            dell = wsolnl[:,1]
            El = wsolnl[:,2]./wsolnl[:,1] .- 1.

            #smoothScaledEnd!(scu, delu, 5)
            #smoothScaledEnd!(scu, Eu, 5)
            #smoothScaledEnd!(scl, dell, 5)
            #smoothScaledEnd!(scl, El, 5)
            
            
            del_prev[:] = del_iter[:]

            #Arrange solution
            
            del_iter[:] = [reverse(dell);delu] 
            E_iter[:] = [reverse(El);Eu] 
            
            #smoothScaledEnd!(sc, del_iter, 5)
            if mod(istep,20) == 0
                figure(1)
                plot(scu, delu)
                plot(scl, dell)
                figure(2)
                plot(scu, Eu)
                plot(scl, El)
            end
            
            #error("here")
            
            #Find suitable naca coefficients to fit the modified airfoil

            newthick = zeros(surf.ndiv)
            newcam = zeros(surf.ndiv)
            for i = 2:surf.ndiv-1
                th_u = 0.5*(qc[surf.ndiv-2+i]*del_iter[surf.ndiv-2+i] + qc[surf.ndiv-1+i]*del_iter[surf.ndiv-1+i])/sqrt(1. + (cam_slope_orig[i] + thick_slope_orig[i])^2)
                th_l = 0.5*(qc[surf.ndiv+1-i]*del_iter[surf.ndiv+1-i] + qc[surf.ndiv-i]*del_iter[surf.ndiv-i])/sqrt(1. + (cam_slope_orig[i] - thick_slope_orig[i])^2)
                newthick[i] = thick_orig[i] + 1/sqrt(Re)*0.5*(th_u + th_l)
                newcam[i] = cam_orig[i] + 1/sqrt(Re)*0.5*(th_u - th_l)
            end
            newthick[1] = thick_orig[1]
            newcam[1] = cam_orig[1]
            
            i = surf.ndiv
            th_u = qc[2*surf.ndiv-2]*del_iter[2*surf.ndiv-2]/sqrt(1. + (cam_slope_orig[i] + thick_slope_orig[i])^2)
            th_l = qc[1]*del_iter[1]/sqrt(1. + (cam_slope_orig[i] - thick_slope_orig[i])^2)    
            newthick[i] = thick_orig[i] + 1/sqrt(Re)*0.5(th_u + th_l)
            newcam[i] = cam_orig[i] + 1/sqrt(Re)*0.5(th_u - th_l)
            
            bstart = [-0.1260; -0.3516; 0.2843; -0.1015]
            thickCoef = find_nacaThickCoef(surf, newthick, bstart)
            bstart = zeros(4)
            camCoef = find_nacaCamCoef(surf, newcam, bstart)
            
            th = parse(Int, surf.coord_file[7:8])/100.
            b1 = 0.2969
            bt = [b1; thickCoef]            
            @. nacath(x) = 5*th*(bt[1]*sqrt(x) + bt[2]*x + bt[3]*x^2 + bt[4]*x^3 + bt[5]*x^4)

            bc = camCoef
            @. nacacam(x) = bc[1]*x + bc[2]*x^2 + bc[3]*x^3 + bc[4]*x^4

            #Find new shape of airfoil
            for i = 1:surf.ndiv
                surf.thick[i] = nacath(surf.x[i])
                surf.thick_slope[i] = ForwardDiff.derivative(nacath, surf.x[i])
                surf.cam[i] = nacacam(surf.x[i])
                surf.cam_slope[i] = ForwardDiff.derivative(nacacam, surf.x[i])
            end
            surf.thick_slope[1] = thick_slope_orig[1]
            #plot(surf.x, surf.cam .+ surf.thick)
            #plot(surf.x, surf.cam .- surf.thick)

            #Check for convergence
            res =  sum(abs.(del_prev .- del_iter))
            println(iter, "   ", res)

            #if iter == iterMax
            if res <= resTol
                println("converged")
                del[:] = del_iter[:]
                E[:] = E_iter[:]
                push!(curfield.tev, TwoDVort(xloc_tev, zloc_tev, tevstr, vcore, 0., 0.))
            end

            # if iter == 3 && mod(istep,10) == 0
            #     figure(1)
            #     plot(surf.x, qu)
            #     figure(2)
            #     plot(surf.x, surf.thick)
            #     axis("equal")
            #     figure(3)
            #     plot(surf.x[2:end], delu)
            # end
            
        end
        

        #Update induced velocities to include effect of last shed vortex
        #add_indbound_lasttev(surf, curfield)
        update_indbound(surf, curfield)

        #Calculate adot
        update_atermdot(surf, dt)
        
        #Set previous values of aterm to be used for derivatives in next time step
        surf.ateprev[1] = surf.ate[1]
        for ia = 1:3
            surf.aprev[ia] = surf.aterm[ia]
        end
        
        #Calculate bound vortex strengths
        update_bv_src(surf)
        
        #Wake rollup
        wakeroll(surf, curfield, dt)

        qu, ql, phi_u, phi_l, cpu, cpl = calc_edgeVel_cp(surf, [curfield.u[1]; curfield.w[1]], phi_u, phi_l, dt)
        
        #println(cpu[end], " ", cpl[end], 0.5*(qu[end]^2 - ql[end]^2)*dt, " ", tevstr)
        
        #Force calculation
        cn, cs, cl, cd, cm = calc_forces(surf, cpu, cpl)
        
        bound_circ = phi_u[end-1] - phi_l[end-1]
        
        # write flow details if required
        if writeflag == 1
            if istep in writeArray
                dirname = "$(round(t,sigdigits=nround))"
                writeStamp(dirname, t, surf, curfield, qu, ql, cpu, cpl)
            end
        end
        vle = qu[1]
        
        stag = find_stag(surf, qu, ql)
        
        mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, vle,
                        cl, cd, cm, cn, cs, bound_circ, stag/surf.c, cpu[1]])

    end
    
    mat = mat'
    
    f = open("resultsSummary", "w")
    Serialization.serialize(f, ["#time \t", "alpha (rad) \t", "h/c \t", "u/uref \t", "Vle \t", "Cl \t", "Cd \t", "Cm \t", "Cn \t", "Cs \t", "bc \t", "xs \n"])
    writedlm(f, mat)
    close(f)
    
    mat, surf, curfield


    
    # # write flow details if required
    # if writeflag == 1
    #     if istep in writeArray
    #         dirname = "$(round(t,sigdigits=nround))"
    #         writeStamp(dirname, t, surf, curfield, qu, ql, cpu, cpl, suc, delu, Eu, thick_orig, quc, qucx, quct)
    #     end
    # end
    
end



function IBL_shape_attached_dblwk(Re, surf::TwoDSurfThick, curfield::TwoDFlowField, nsteps::Int64 = 300, dtstar::Float64 = 0.015, startflag = 0, writeflag = 0, writeInterval = 1000., delvort = delNone(); maxwrite = 50, nround=6)
    
    # If a restart directory is provided, read in the simulation data
    if startflag == 0
        mat = zeros(0, 13)
        t = 0.
    elseif startflag == 1
        dirvec = readdir()
        dirresults = map(x->(v = tryparse(Float64,x); isnull(v) ? 0.0 : get(v)),dirvec)
        latestTime = maximum(dirresults)
        mat = readdlm("resultsSummary")
        t = mat[end,1]
    else
        throw("invalid start flag, should be 0 or 1")
    end
    mat = mat'

    dt = dtstar*surf.c/surf.uref

    # if writeflag is on, determine the timesteps to write at
    if writeflag == 1
        writeArray = Int64[]
        tTot = nsteps*dt
        for i = 1:maxwrite
            tcur = writeInterval*real(i)
            if t > tTot
                break
            else
                push!(writeArray, Int(round(tcur/dt)))
            end
        end
    end

    vcore = 1.3*dt*surf.c

    int_wax = zeros(surf.ndiv)
    int_c = zeros(surf.ndiv)
    int_t = zeros(surf.ndiv)

    #At this time coded for symmetric situations only
    qu = zeros(surf.ndiv)
    ql = zeros(surf.ndiv)

    spos = zeros(2*surf.ndiv-1)
    sc = zeros(2*surf.ndiv-2)
    dsc = zeros(2*surf.ndiv-2)
    qc = zeros(2*surf.ndiv-2)
    qc_prev = zeros(2*surf.ndiv-2)
    qcx = zeros(2*surf.ndiv-2)
    qct = zeros(2*surf.ndiv-2)
    su = zeros(surf.ndiv)
    
    thick_orig = zeros(surf.ndiv)
    thick_orig[:] = surf.thick[:]
    thick_slope_orig = zeros(surf.ndiv)
    thick_slope_orig[:] = surf.thick_slope[:]

    cam_orig = zeros(surf.ndiv)
    cam_orig[:] = surf.cam[:]
    cam_slope_orig = zeros(surf.ndiv)
    cam_slope_orig[:] = surf.cam_slope[:]

    
    #Initialise boundary layer
    del, _, E, _ = initDelE(2*surf.ndiv-2)

    spos[1] = 0.
    dsdx = zeros(surf.ndiv)
    for i = 2:surf.ndiv
        dsdx[i] = sqrt(1 + (surf.cam_slope[i] - surf.thick_slope[i])^2)
    end
    dsdx[1] = dsdx[2]
    su[1] = 0.
    for i = surf.ndiv-1:-1:1
        su[surf.ndiv+1-i] = simpleTrapz(dsdx[surf.ndiv:-1:i], surf.x[surf.ndiv:-1:i])
    end
    spos[1:surf.ndiv] = -su[:]
    for i = 2:surf.ndiv
        dsdx[i] = sqrt(1 + (surf.cam_slope[i] + surf.thick_slope[i])^2)
    end
    dsdx[1] = dsdx[2]
    su[1] = 0.
    for i = 2:surf.ndiv
        su[i] = simpleTrapz(dsdx[1:i], surf.x[1:i])
    end
    spos[surf.ndiv+1:2*surf.ndiv-1] = spos[surf.ndiv] .+ su[2:surf.ndiv]
    for i = 1:2*surf.ndiv-2
        sc[i] = (spos[i] + spos[i+1])/2
    end

    dspos = diff(spos)
    
    #All BL parameters go from trailing edge on lower surface to leading edge to trailing edge on uper surface. 

    
    phi_u = zeros(surf.ndiv)
    phi_l = zeros(surf.ndiv)

    x_w = collect(surf.c*1.01:surf.c*0.01:surf.c*3)
    nw = length(x_w)
    wfn = zeros(nw)

    bound_circ = 0.
    
    tevstr = zeros(100)
    tevdist = zeros(100)
    restev = zeros(100,2)
    restev_prev = zeros(2,2)
    
    phi_u_temp = zeros(surf.ndiv)
    phi_l_temp = zeros(surf.ndiv)

    for istep = 1:nsteps

        #Udpate current time
        t = t + dt

        #Update kinematic parameters
        update_kinem(surf, t)

        #Update flow field parameters if any
        update_externalvel(curfield, t)

        #Update bound vortex positions
        update_boundpos(surf, dt)

        #Update incduced velocities on airfoil
        update_indbound(surf, curfield)

        ntev = length(curfield.tev)
        zu = surf.bnd_z_u[surf.ndiv] + del[2*surf.ndiv-2]/sqrt(Re)
        zl = surf.bnd_z_l[surf.ndiv] - del[1]/sqrt(Re)
        
        if ntev == 0
            xloc_tev = surf.bnd_x_chord[surf.ndiv] + 0.5*surf.kinem.u*dt*cos(surf.kinem.alpha)
            zloc_tev = zl - 0.5*surf.kinem.u*dt*sin(surf.kinem.alpha)
            xloc_lev = surf.bnd_x_chord[surf.ndiv] + 0.5*surf.kinem.u*dt*cos(surf.kinem.alpha)
            zloc_lev = zu - 0.5*surf.kinem.u*dt*sin(surf.kinem.alpha)
        else
            xloc_tev = surf.bnd_x_chord[surf.ndiv] + (1. /3.)*(curfield.tev[ntev].x - surf.bnd_x_chord[surf.ndiv])
            zloc_tev = zl + (1. /3.)*(curfield.tev[ntev].z - zl)
            xloc_lev = surf.bnd_x_chord[surf.ndiv] + (1. /3.)*(curfield.lev[ntev].x - surf.bnd_x_chord[surf.ndiv])
            zloc_lev = zu + (1. /3.)*(curfield.tev[ntev].z - zu)
        end
        
        #Outer iteration to solve coupled problem
        iter = 0
        res = 1
        del_iter = zeros(2*surf.ndiv-2)
        del_prev = zeros(2*surf.ndiv-2)
        soln = zeros(2*surf.naterm+1)
        E_iter = zeros(2*surf.ndiv-2)

        resTol = 1e-5
        iterMax = 20

        qc_prev[:] = qc[:]       

        #while (iter < iterMax)
        while res > resTol
            
            iter += 1

            function tev_iter!(FC, J, x)

                #x vector = [aterm; bterm; ate; tevstr]
                
                #FC = zeros(2*surf.naterm+2)
                #           J = zeros(2*surf.naterm+2, 2*surf.naterm+2)

                aterm = x[1:surf.naterm]
                bterm = x[surf.naterm+1:2*surf.naterm]
                ate = x[2*surf.naterm+1]
                tevstr = x[2*surf.naterm+2]
                levstr = x[2*surf.naterm+3]
                
                dummyvort1 = TwoDVort(xloc_tev, zloc_tev, tevstr, vcore, 0., 0.)
                dummyvort2 = TwoDVort(xloc_lev, zloc_lev, levstr, vcore, 0., 0.)
                
                uu, wu = ind_vel([dummyvort1; dummyvort2], surf.bnd_x_u, surf.bnd_z_u)
                ul, wl = ind_vel([dummyvort1; dummyvort2], surf.bnd_x_l, surf.bnd_z_l)
                
                wlz = 0.5*((surf.wind_u .+ wu).*cos(surf.kinem.alpha) .+ (surf.uind_u .+ uu).*sin(surf.kinem.alpha) .+
                           (surf.wind_l .+ wl)*cos(surf.kinem.alpha) .+ (surf.uind_l .+ ul)*sin(surf.kinem.alpha))
                
                wtz = 0.5*((surf.wind_u .+ wu).*cos(surf.kinem.alpha) .+ (surf.uind_u .+ uu).*sin(surf.kinem.alpha) .-
                           (surf.wind_l .+ wl).*cos(surf.kinem.alpha) .- (surf.uind_l .+ ul).*sin(surf.kinem.alpha))
                
                wtx = 0.5*((surf.uind_u .+ uu).*cos(surf.kinem.alpha) .- (surf.wind_u .+ wu).*sin(surf.kinem.alpha) .+
                           (surf.uind_l .+ ul).*cos(surf.kinem.alpha) .- (surf.wind_l .+ wl).*sin(surf.kinem.alpha))
                
                wlx = 0.5*((surf.uind_u .+ uu).*cos(surf.kinem.alpha) .- (surf.wind_u .+ wu).*sin(surf.kinem.alpha) .-
                           (surf.uind_l .+ ul).*cos(surf.kinem.alpha) .+ (surf.wind_l .+ wl).*sin(surf.kinem.alpha))

                dummyvort1 = TwoDVort(xloc_tev, zloc_tev, 1., vcore, 0., 0.)
                dummyvort2 = TwoDVort(xloc_lev, zloc_lev, 1., vcore, 0., 0.)
                
                uu, wu = ind_vel([dummyvort1], surf.bnd_x_u, surf.bnd_z_u)
                ul, wl = ind_vel([dummyvort1], surf.bnd_x_l, surf.bnd_z_l)
                
                
                wlz_t = 0.5*(wu.*cos(surf.kinem.alpha) .+ uu.*sin(surf.kinem.alpha) .+
                             wl*cos(surf.kinem.alpha) .+ ul*sin(surf.kinem.alpha))
                wtz_t = 0.5*(wu.*cos(surf.kinem.alpha) .+ uu.*sin(surf.kinem.alpha) .-
                             wl.*cos(surf.kinem.alpha) .- ul.*sin(surf.kinem.alpha))
                wtx_t = 0.5*(uu.*cos(surf.kinem.alpha) .- wu.*sin(surf.kinem.alpha) .+
                             ul.*cos(surf.kinem.alpha) .- wl.*sin(surf.kinem.alpha))
                wlx_t = 0.5*(uu.*cos(surf.kinem.alpha) .- wu.*sin(surf.kinem.alpha) .-
                             ul.*cos(surf.kinem.alpha) .+ wl.*sin(surf.kinem.alpha))
                uu, wu = ind_vel([dummyvort2], surf.bnd_x_u, surf.bnd_z_u)
                ul, wl = ind_vel([dummyvort2], surf.bnd_x_l, surf.bnd_z_l)

                
                wlz_l = 0.5*(wu.*cos(surf.kinem.alpha) .+ uu.*sin(surf.kinem.alpha) .+
                             wl*cos(surf.kinem.alpha) .+ ul*sin(surf.kinem.alpha))
                wtz_l = 0.5*(wu.*cos(surf.kinem.alpha) .+ uu.*sin(surf.kinem.alpha) .-
                             wl.*cos(surf.kinem.alpha) .- ul.*sin(surf.kinem.alpha))
                wtx_l = 0.5*(uu.*cos(surf.kinem.alpha) .- wu.*sin(surf.kinem.alpha) .+
                             ul.*cos(surf.kinem.alpha) .- wl.*sin(surf.kinem.alpha))
                wlx_l = 0.5*(uu.*cos(surf.kinem.alpha) .- wu.*sin(surf.kinem.alpha) .-
                             ul.*cos(surf.kinem.alpha) .+ wl.*sin(surf.kinem.alpha))
                
                rng = 1:surf.naterm

                Lx = zeros(surf.ndiv); Lz = zeros(surf.ndiv)
                Tx = zeros(surf.ndiv); Tz = zeros(surf.ndiv)
                phi_u_integ = zeros(surf.ndiv)
                phi_l_integ = zeros(surf.ndiv)
                
                i = surf.ndiv-1
                vref_x_u = (surf.kinem.u + curfield.u[1])*cos(surf.kinem.alpha) + (surf.kinem.hdot - curfield.w[1])*sin(surf.kinem.alpha) - surf.kinem.alphadot*(surf.cam[i] + surf.thick[i])
                vref_x_l = (surf.kinem.u + curfield.u[1])*cos(surf.kinem.alpha) + (surf.kinem.hdot - curfield.w[1])*sin(surf.kinem.alpha) - surf.kinem.alphadot*(surf.cam[i] - surf.thick[i])
                vref_z = (surf.kinem.u + curfield.u[1])*sin(surf.kinem.alpha) - (surf.kinem.hdot - curfield.w[1])*cos(surf.kinem.alpha) + surf.kinem.alphadot*(surf.x[i] - surf.pvt*surf.c)
                
                for i = 1:surf.ndiv
                    Lx[i] = surf.uref*(sum(aterm[rng]'*sin.(rng*surf.theta[i])) + ate*tan(surf.theta[i]/2))
                    Lz[i] = surf.uref*(sum(aterm[rng]'*cos.(rng*surf.theta[i])) + ate)
                    Tz[i] = surf.uref*sum(bterm[rng]'*sin.(rng*surf.theta[i]))
                    Tx[i] = -surf.uref*sum(bterm[rng]'*cos.(rng*surf.theta[i]))
                    if i ==1
                        phi_u_integ[i] = Lz[i] + Tz[i] + wtz[i] + wlz[i]
                        phi_l_integ[i] = -(Lz[i] - Tz[i] - wtz[i] +wlz[i] )
                    else
                        phi_u_integ[i] = (Lx[i] + Tx[i] + wtx[i] + wlx[i]) + (surf.cam_slope[i] + surf.thick_slope[i])*(Lz[i] + Tz[i] + wtz[i] + wlz[i])
                        phi_l_integ[i] = (-Lx[i] + Tx[i] + wtx[i] - wlx[i]) + (surf.cam_slope[i] - surf.thick_slope[i])*(Lz[i] - Tz[i] - wtz[i] + wlz[i])
                    end    
                end
                
                phi_u_bc = 0; phi_l_bc = 0
                for i = 2:surf.ndiv-1
                    phi_u_bc += 0.5*(phi_u_integ[i]/sqrt(1. + (surf.thick_slope[i] + surf.cam_slope[i])^2) + phi_u_integ[i-1]/sqrt(1. + (surf.thick_slope[i-1] + surf.cam_slope[i-1])^2))*sqrt((surf.x[i] - surf.x[i-1])^2 + (surf.cam[i] + surf.thick[i] - surf.cam[i-1] - surf.thick[i-1])^2)
                    phi_l_bc += 0.5*(phi_l_integ[i]/sqrt(1. + (-surf.thick_slope[i] + surf.cam_slope[i])^2) + phi_l_integ[i-1]/sqrt(1. + (-surf.thick_slope[i-1] + surf.cam_slope[i-1])^2))*sqrt((surf.x[i] - surf.x[i-1])^2 + (surf.cam[i] - surf.thick[i] - surf.cam[i-1] + surf.thick[i-1])^2)
                end

                if !(FC == nothing)
                    for i = 2:surf.ndiv-1
                        rhs_l = -(surf.kinem.u + curfield.u[1])*sin(surf.kinem.alpha) - surf.kinem.alphadot*(surf.x[i] - surf.pvt*surf.c) + (surf.kinem.hdot - curfield.w[1])*cos(surf.kinem.alpha) - wlz[i] + surf.cam_slope[i]*((surf.kinem.u + curfield.u[1])*cos(surf.kinem.alpha) + (surf.kinem.hdot - curfield.w[1])*sin(surf.kinem.alpha) + wtx[i] - surf.kinem.alphadot*surf.cam[i]) + surf.thick_slope[i]*(wlx[i] - surf.kinem.alphadot*surf.thick[i])
                        rhs_nonl = surf.cam_slope[i]*(wlx[i] - surf.kinem.alphadot*surf.thick[i]) + surf.thick_slope[i]*((surf.kinem.u + curfield.u[1])*cos(surf.kinem.alpha) + (surf.kinem.hdot - curfield.w[1])*sin(surf.kinem.alpha) + wtx[i] - surf.kinem.alphadot*surf.cam[i]) - wtz[i]
                        
                        FC[i-1] = Lz[i] - surf.cam_slope[i]*Tx[i] - surf.thick_slope[i]*Lx[i] - rhs_l
                        FC[surf.ndiv+i-3] = Tz[i] - surf.cam_slope[i]*Lx[i] - surf.thick_slope[i]*Tx[i] - rhs_nonl
                        
                    end
                    
                    #Kutta condition
                    i = surf.ndiv-1
                    qu = sqrt((vref_x_u + Lx[i] + Tx[i] + wlx[i] + wtx[i])^2 + (vref_z + Lz[i] + Tz[i] + wlz[i] + wtz[i])^2)
                    ql = sqrt((vref_x_l - Lx[i] + Tx[i] - wlx[i] + wtx[i])^2 + (vref_z + Lz[i] - Tz[i] + wlz[i] - wtz[i])^2)
                    FC[2*surf.ndiv-3] = levstr + tevstr - 0.5*(qu^2 - ql^2)*dt

                    #Kelvin condition
                    bc = phi_u_bc - phi_l_bc
                    FC[2*surf.ndiv-2] = bc - bound_circ + tevstr + levstr
                    
                    #LEV condition
                    FC[2*surf.ndiv-1] = levstr - 0.5*qu^2*dt
                end

                if !(J == nothing)
                    for i = 2:surf.ndiv-1
                        for n = 1:surf.naterm
                            J[i-1,n] = cos(n*surf.theta[i]) - surf.thick_slope[i]*sin(n*surf.theta[i]) 
                            J[i-1,n+surf.naterm] = surf.cam_slope[i]*cos(n*surf.theta[i])
                            J[surf.ndiv+i-3,n] = -surf.cam_slope[i]*sin(n*surf.theta[i])
                            J[surf.ndiv+i-3,surf.naterm+n] = sin(n*surf.theta[i]) + surf.thick_slope[i]*cos(n*surf.theta[i])            
                        end
                        J[i-1,2*surf.naterm+1] = 1. - surf.thick_slope[i]*tan(surf.theta[i]/2)
                        J[surf.ndiv+i-3,2*surf.naterm+1] = -surf.cam_slope[i]*tan(surf.theta[i]/2)
                        J[i-1,2+2*surf.naterm] = -surf.cam_slope[i]*wtx_t[i] - surf.thick_slope[i]*wlx_t[i] + wlz_t[i]
                        J[i-1,3+2*surf.naterm] = -surf.cam_slope[i]*wtx_l[i] - surf.thick_slope[i]*wlx_l[i] + wlz_l[i]
                        J[surf.ndiv+i-3,2+2*surf.naterm] = -surf.cam_slope[i]*wlx_t[i] - surf.thick_slope[i]*wtx_t[i] + wtz_t[i]
                        J[surf.ndiv+i-3,3+2*surf.naterm] = -surf.cam_slope[i]*wlx_l[i] - surf.thick_slope[i]*wtx_l[i] + wtz_l[i]
                    end
                    
                    #Kutta condition
                    i = surf.ndiv-1
                    for n = 1:surf.naterm
                        J[2*surf.ndiv-3,n] = -2*dt*(0.5*(vref_x_u + vref_x_l) + Tx[i] + wtx[i])*sin(n*surf.theta[i]) - 2*dt*(Tz[i] + wtz[i])*cos(n*surf.theta[i])
                        J[2*surf.ndiv-3,n+surf.naterm] = 2*dt*(0.5*(vref_x_u - vref_x_l) + Lx[i] + wlx[i])*cos(n*surf.theta[i]) - 2*dt*(vref_z + Lz[i] + wlz[i])*sin(n*surf.theta[i])
                    end
                    J[2*surf.ndiv-3,2*surf.naterm+1] = -2*dt*(0.5*(vref_x_u + vref_x_l) + Tx[i] + wtx[i])*tan(surf.theta[i]/2) - 2*dt*(Tz[i] + wtz[i])
                    J[2*surf.ndiv-3,2*surf.naterm+2] = 1. - 2*dt*(0.5*(vref_x_u + vref_x_l) + Tx[i] + wtx[i])*wlx_t[i] - 2*dt*(0.5*(vref_x_u - vref_x_l) + Lx[i] + wlx[i])*wtx_t[i] - 2*dt*(vref_z + Lz[i] + wlz[i])*wtz_t[i] -2*dt*(Tz[i] + wtz[i])*wlz_t[i]
                    J[2*surf.ndiv-3,2*surf.naterm+3] = 1. - 2*dt*(0.5*(vref_x_u + vref_x_l) + Tx[i] + wtx[i])*wlx_l[i] - 2*dt*(0.5*(vref_x_u - vref_x_l) + Lx[i] + wlx[i])*wtx_l[i] - 2*dt*(vref_z + Lz[i] + wlz[i])*wtz_l[i] -2*dt*(Tz[i] + wtz[i])*wlz_l[i]
                    utotu = vref_x_u + Lx[i] + Tx[i] + wlx[i] + wtx[i]
                    wtotu = vref_z + Lz[i] + Tz[i] + wlz[i] + wtz[i]
                    
                    for n = 1:surf.naterm
                        J[2*surf.ndiv-1,n] = -(utotu*sin(n*surf.theta[i]) + wtotu*cos(n*surf.theta[i]))*dt
                        J[2*surf.ndiv-1,n+surf.naterm] = -(-utotu*cos(n*surf.theta[i]) + wtotu*sin(n*surf.theta[i]))*dt 
                    end
                    J[2*surf.ndiv-1,2*surf.naterm+1] = -(utotu*tan(surf.theta[i]/2) + wtotu)*dt
                    J[2*surf.ndiv-1,2*surf.naterm+2] = -(utotu*(wlx_t[i] + wtx_t[i]) + wtotu*(wlz_t[i] + wtz_t[i]))*dt
                    J[2*surf.ndiv-1,2*surf.naterm+3] = 1. -(utotu*(wlx_l[i] + wtx_l[i]) + wtotu*(wlz_l[i] + wtz_l[i]))*dt
                    
                    J[2*surf.ndiv-2,:] .= 0.
                    #Kelvin condition
                    i = 1
                    ate_u = tan(surf.theta[i]/2) + (surf.cam_slope[i] + surf.thick_slope[i])
                    ate_l = -tan(surf.theta[i]/2) + (surf.cam_slope[i] - surf.thick_slope[i])
                    tev_u = wtx_t[i] + wlx_t[i] + (surf.cam_slope[i] + surf.thick_slope[i])*(wlz_t[i] + wtz_t[i])
                    tev_l = wtx_t[i] - wlx_t[i] + (surf.cam_slope[i] - surf.thick_slope[i])*(wlz_t[i] - wtz_t[i])
                    lev_u = wtx_l[i] + wlx_l[i] + (surf.cam_slope[i] + surf.thick_slope[i])*(wlz_l[i] + wtz_l[i])
                    lev_l = wtx_l[i] - wlx_l[i] + (surf.cam_slope[i] - surf.thick_slope[i])*(wlz_l[i] - wtz_l[i])
                    
                    den_u = sqrt(1. + (surf.cam_slope[i] + surf.thick_slope[i])^2)
                    den_l = sqrt(1. + (surf.cam_slope[i] - surf.thick_slope[i])^2)
                    ate_u_int = 0; ate_l_int = 0; tev_u_int = 0; tev_l_int = 0

                    for i = 2:surf.ndiv-1
                        
                        ate_u_p = ate_u
                        tev_u_p = tev_u
                        lev_u_p = lev_u
                        ate_l_p = ate_l
                        tev_l_p = tev_l
                        lev_l_p = lev_l
                        den_u_p = den_u
                        den_l_p = den_l
                        
                        ate_u = tan(surf.theta[i]/2) + (surf.cam_slope[i] + surf.thick_slope[i])
                        ate_l = -tan(surf.theta[i]/2) + (surf.cam_slope[i] - surf.thick_slope[i])
                        tev_u = wtx_t[i] + wlx_t[i] + (surf.cam_slope[i] + surf.thick_slope[i])*(wlz_t[i] + wtz_t[i])
                        tev_l = wtx_t[i] - wlx_t[i] + (surf.cam_slope[i] - surf.thick_slope[i])*(wlz_t[i] - wtz_t[i])
                        lev_u = wtx_l[i] + wlx_l[i] + (surf.cam_slope[i] + surf.thick_slope[i])*(wlz_l[i] + wtz_l[i])
                        lev_l = wtx_l[i] - wlx_l[i] + (surf.cam_slope[i] - surf.thick_slope[i])*(wlz_l[i] - wtz_l[i])
                        den_u = sqrt(1. + (surf.cam_slope[i] + surf.thick_slope[i])^2)
                        den_l = sqrt(1. + (surf.cam_slope[i] - surf.thick_slope[i])^2)
                        
                        ds_u = sqrt((surf.x[i] - surf.x[i-1])^2 + (surf.cam[i] + surf.thick[i] - surf.cam[i-1] - surf.thick[i-1])^2)
                        ds_l = sqrt((surf.x[i] - surf.x[i-1])^2 + (surf.cam[i] - surf.thick[i] - surf.cam[i-1] + surf.thick[i-1])^2)
                        J[2*surf.ndiv-2,2*surf.naterm+1] += 0.5*(ate_u/den_u + ate_u_p/den_u_p)*ds_u 
                        J[2*surf.ndiv-2,2*surf.naterm+1] -= 0.5*(ate_l/den_l + ate_l_p/den_l_p)*ds_l
                        
                        J[2*surf.ndiv-2,2*surf.naterm+2] += 0.5*(tev_u/den_u + tev_u_p/den_u_p)*ds_u
                        J[2*surf.ndiv-2,2*surf.naterm+2] -= 0.5*(tev_l/den_l + tev_l_p/den_l_p)*ds_l
                        J[2*surf.ndiv-2,2*surf.naterm+3] += 0.5*(lev_u/den_u + lev_u_p/den_u_p)*ds_u
                        J[2*surf.ndiv-2,2*surf.naterm+3] -= 0.5*(lev_l/den_l + lev_l_p/den_l_p)*ds_l
                    end
                    
                    J[2*surf.ndiv-2,2*surf.naterm+2] += 1.
                    J[2*surf.ndiv-2,2*surf.naterm+3] += 1.
                    
                    
                    for n = 1:surf.naterm
                        i = 1
                        an_u = sin(n*surf.theta[i]) + (surf.cam_slope[i] + surf.thick_slope[i])*cos(n*surf.theta[i])
                        an_l = -sin(n*surf.theta[i]) + (surf.cam_slope[i] - surf.thick_slope[i])*cos(n*surf.theta[i])
                        bn_u = -cos(n*surf.theta[i]) + (surf.cam_slope[i] + surf.thick_slope[i])*sin(n*surf.theta[i])
                        bn_l = -cos(n*surf.theta[i]) - (surf.cam_slope[i] - surf.thick_slope[i])*sin(n*surf.theta[i])
                        den_u = sqrt(1. + (surf.cam_slope[i] + surf.thick_slope[i])^2)
                        den_l = sqrt(1. + (surf.cam_slope[i] - surf.thick_slope[i])^2)
                        
                        an_u_int = 0; an_l_int = 0; bn_u_int = 0; bn_l_int = 0;
                        
                        for i = 2:surf.ndiv-1

                            an_u_p = an_u
                            bn_u_p = bn_u
                            an_l_p = an_l
                            bn_l_p = bn_l
                            den_u_p = den_u
                            den_l_p = den_l
                            
                            an_u = sin(n*surf.theta[i]) + (surf.cam_slope[i] + surf.thick_slope[i])*cos(n*surf.theta[i])
                            an_l = -sin(n*surf.theta[i]) + (surf.cam_slope[i] - surf.thick_slope[i])*cos(n*surf.theta[i])
                            bn_u = -cos(n*surf.theta[i]) + (surf.cam_slope[i] + surf.thick_slope[i])*sin(n*surf.theta[i])
                            bn_l = -cos(n*surf.theta[i]) - (surf.cam_slope[i] - surf.thick_slope[i])*sin(n*surf.theta[i])
                            den_u = sqrt(1. + (surf.cam_slope[i] + surf.thick_slope[i])^2)
                            den_l = sqrt(1. + (surf.cam_slope[i] - surf.thick_slope[i])^2)

                            ds_u = sqrt((surf.x[i] - surf.x[i-1])^2 + (surf.cam[i] + surf.thick[i] - surf.cam[i-1] - surf.thick[i-1])^2)
                            ds_l = sqrt((surf.x[i] - surf.x[i-1])^2 + (surf.cam[i] - surf.thick[i] - surf.cam[i-1] + surf.thick[i-1])^2)
                            J[2*surf.ndiv-2,n] += 0.5*(an_u/den_u + an_u_p/den_u_p)*ds_u
                            J[2*surf.ndiv-2,n+surf.naterm] += 0.5*(bn_u/den_u + bn_u_p/den_u_p)*ds_u
                            J[2*surf.ndiv-2,n] -= 0.5*(an_l/den_l + an_l_p/den_l_p)*ds_l
                            J[2*surf.ndiv-2,n+surf.naterm] -= 0.5*(bn_l/den_l + bn_l_p/den_l_p)*ds_l
                            
                        end
                    end
                end
                    
            end
            
            xstart = [surf.aterm; surf.bterm; surf.ate[1]; -0.001; 0.001]
            
            soln = nlsolve(only_fj!(tev_iter!), xstart, xtol = 1e-6, method = :newton)
            soln = soln.zero
            
            #assign the solution
            surf.aterm[:] = soln[1:surf.naterm]
            surf.bterm[:] = soln[surf.naterm+1:2*surf.naterm]
            surf.ate[1] = soln[2*surf.naterm+1]
            tevstr = soln[2*surf.naterm+2]
            levstr = soln[2*surf.naterm+3]
            
            temptev = TwoDVort(xloc_tev, zloc_tev, tevstr, vcore, 0., 0.)
            templev = TwoDVort(xloc_lev, zloc_lev, levstr, vcore, 0., 0.)
            
            # a_wfn = log10((qc[end]*del_iter[end] + qc[1]*del_iter[1])/sqrt(Re))
            # wfn[1] = 1/sqrt(Re)*((qc[end]*del_iter[end] - qc[end-1]*del_iter[end-1])/(surf.x[end] - surf.x[end-1]) + (qc[1]*del_iter[1] - qc[2]*del_iter[2])/(surf.x[1] - surf.x[2]))
            # for i = 2:nw
            #     wfn[i] = (10^(a_wfn - 3.2*(x_w[i] - surf.c)) - 10^(a_wfn - 3.2*(x_w[i-1] - surf.c)))/(x_w[i] - x_w[i-1])
            # end
            
            # #Source strengths in wake and induced velocity
            # uind_src = zeros(surf.ndiv)
            # for i = 1:surf.ndiv
            #     for iw = 1:nw-1
            #         str = 0.5(wfn[iw] + wfn[iw+1])
            #         xloc = 0.5*(x_w[iw] + x_w[iw+1])
            #         uind_src[i] += 1/(2*pi)*str/(xloc - surf.x[i])*(x_w[iw+1] - x_w[iw])
            #     end
            # end
            # uind_src[:] .= 0.
            # surf.uind_u[:] += uind_src[:]
            # surf.uind_l[:] += uind_src[:]

            ind_new_u_u, ind_new_w_u = ind_vel([temptev; templev], surf.bnd_x_u, surf.bnd_z_u)
            ind_new_u_l, ind_new_w_l = ind_vel([temptev; templev], surf.bnd_x_l, surf.bnd_z_l)
            
            surf.uind_u[:] += ind_new_u_u[:]
            surf.wind_u[:] += ind_new_w_u[:]
            surf.uind_l[:] += ind_new_u_l[:]
            surf.wind_l[:] += ind_new_w_l[:]

            qu, ql = calc_edgeVel(surf, [curfield.u[1], curfield.w[1]])
            
            #surf.uind_u[:] -= uind_src[:]
            #surf.uind_l[:] -= uind_src[:]
            
            surf.uind_u[:] -= ind_new_u_u[:]
            surf.wind_u[:] -= ind_new_w_u[:]
            surf.uind_l[:] -= ind_new_u_l[:]
            surf.wind_l[:] -= ind_new_w_l[:]
            
            #smoothScaledEnd!(surf.x, qu, 10)

            #Full vector of qc with positive going from LE to TE
            for i = 1:surf.ndiv-1
                qc[i] = (ql[surf.ndiv-i] + ql[surf.ndiv+1-i])/2
                qc[surf.ndiv-1+i] = (qu[i] + qu[i+1])/2
            end
            
            qcx[1:surf.ndiv-1] = -diff(qc[1:surf.ndiv])./diff(sc[1:surf.ndiv])
            #qcx[1] = qcx[2]
            qcx[surf.ndiv+1:2*surf.ndiv-2] = diff(qc[surf.ndiv:2*surf.ndiv-2])./diff(sc[surf.ndiv:2*surf.ndiv-2])
            qcx[surf.ndiv] = qcx[surf.ndiv+1]
            #smoothEnd!(qucx, 10)
            
            if istep == 1
                qct[:] .= 0.
            else
                qct[:] = (qc[:] - qc_prev[:])/dt
            end

            qc_prev[:] = qc[:]

            
            stag = find_stag(surf, qu, ql)

            #Find the closest point to stagnation
            xs_ind = argmin(abs.(stag .- surf.x))

            #Split into two problems 
            if qu[1] > 0.
                #Stag on lower surface
                scu = [sc[surf.ndiv+1-xs_ind:surf.ndiv-1]; sc[surf.ndiv:2*surf.ndiv-2]] .- sc[surf.ndiv+1-xs_ind]
                scl = sc[surf.ndiv-xs_ind] .- reverse(sc[1:surf.ndiv-xs_ind]) 
                qcu = [-qc[surf.ndiv+1-xs_ind:surf.ndiv-1]; qc[surf.ndiv:2*surf.ndiv-2]]
                qcl = reverse(qc[1:surf.ndiv-xs_ind]) 
                qcux = [-qcx[surf.ndiv+1-xs_ind:surf.ndiv-1]; qcx[surf.ndiv:2*surf.ndiv-2]]
                qclx = reverse(qcx[1:surf.ndiv-xs_ind]) 
                qcut = [-qct[surf.ndiv+1-xs_ind:surf.ndiv-1]; qct[surf.ndiv:2*surf.ndiv-2]]
                qclt = reverse(qct[1:surf.ndiv-xs_ind]) 
                delu = [del[surf.ndiv+1-xs_ind:surf.ndiv-1]; del[surf.ndiv:2*surf.ndiv-2]]
                dell = reverse(del[1:surf.ndiv-xs_ind]) 
                Eu = [E[surf.ndiv+1-xs_ind:surf.ndiv-1]; E[surf.ndiv:2*surf.ndiv-2]]
                El = reverse(E[1:surf.ndiv-xs_ind])
                dsposu = [dspos[surf.ndiv+1-xs_ind:surf.ndiv-1]; dspos[surf.ndiv:2*surf.ndiv-2]]
                dsposl = reverse(dspos[1:surf.ndiv-xs_ind])
            else
                scu = sc[surf.ndiv-1+xs_ind:2*surf.ndiv-2] .- sc[surf.ndiv-1+xs_ind]
                scl = -([reverse(sc[surf.ndiv:surf.ndiv-2+xs_ind]); reverse(sc[1:surf.ndiv-1])] .- sc[surf.ndiv-2+xs_ind])
                qcu = qc[surf.ndiv-1+xs_ind:2*surf.ndiv-2]
                qcl = [-reverse(qc[surf.ndiv:surf.ndiv-2+xs_ind]); reverse(qc[1:surf.ndiv-1])]
                qcux = qcx[surf.ndiv-1+xs_ind:2*surf.ndiv-2]
                qclx = [-reverse(qcx[surf.ndiv:surf.ndiv-2+xs_ind]); reverse(qcx[1:surf.ndiv-1])]
                qcut = qct[surf.ndiv-1+xs_ind:2*surf.ndiv-2]
                qclt = [-reverse(qct[surf.ndiv:surf.ndiv-2+xs_ind]); reverse(qct[1:surf.ndiv-1])]
                delu = del[surf.ndiv-1+xs_ind:2*surf.ndiv-2]
                dell = [reverse(del[surf.ndiv:surf.ndiv-2+xs_ind]); reverse(del[1:surf.ndiv-1])]
                Eu = E[surf.ndiv-1+xs_ind:2*surf.ndiv-2]
                El = [reverse(E[surf.ndiv:surf.ndiv-2+xs_ind]); reverse(E[1:surf.ndiv-1])]
                dsposu = dspos[surf.ndiv-1+xs_ind:2*surf.ndiv-2]
                dsposl = [reverse(dspos[surf.ndiv:surf.ndiv-2+xs_ind]); reverse(dspos[1:surf.ndiv-1])]
            end

            qcux = reverse(smoothScaledEnd!(reverse(scu), reverse(qcux), 5))
            qclx = reverse(smoothScaledEnd!(reverse(scl), reverse(qclx), 5))
            smoothScaledEnd!(scu, qcux, 10)
            smoothScaledEnd!(scl, qclx, 10)

            #plot(scu, qcux)
            #plot(scl, qclx)
            
            
            #Solve the FV problem at cell centres
            #plot(sc)
            #plot(sc, qcx)
            
            wu = [delu delu.*(Eu .+ 1)]
            wl = [dell dell.*(El .+ 1)]
            
            wsolnu, i_sep = FVMIBLgridvar(wu, qcu, qcut, qcux, dsposu, t-dt, t)
            wsolnl, i_sep = FVMIBLgridvar(wl, qcl, qclt, qclx, dsposl, t-dt, t) 

            #Assign the solution
            delu = wsolnu[:,1]
            Eu = wsolnu[:,2]./wsolnu[:,1] .- 1.
            dell = wsolnl[:,1]
            El = wsolnl[:,2]./wsolnl[:,1] .- 1.

            #smoothScaledEnd!(scu, delu, 5)
            #smoothScaledEnd!(scu, Eu, 5)
            #smoothScaledEnd!(scl, dell, 5)
            #smoothScaledEnd!(scl, El, 5)
            
            
            del_prev[:] = del_iter[:]

            #Arrange solution
            
            del_iter[:] = [reverse(dell);delu] 
            E_iter[:] = [reverse(El);Eu] 
            
            #smoothScaledEnd!(sc, del_iter, 5)
            if mod(istep,20) == 0
                figure(1)
                plot(scu, delu)
                plot(scl, dell)
                figure(2)
                plot(scu, Eu)
                plot(scl, El)
            end
            
            #error("here")
            
            #Find suitable naca coefficients to fit the modified airfoil

            newthick = zeros(surf.ndiv)
            newcam = zeros(surf.ndiv)
            for i = 2:surf.ndiv-1
                th_u = 0.5*(qc[surf.ndiv-2+i]*del_iter[surf.ndiv-2+i] + qc[surf.ndiv-1+i]*del_iter[surf.ndiv-1+i])/sqrt(1. + (cam_slope_orig[i] + thick_slope_orig[i])^2)
                th_l = 0.5*(qc[surf.ndiv+1-i]*del_iter[surf.ndiv+1-i] + qc[surf.ndiv-i]*del_iter[surf.ndiv-i])/sqrt(1. + (cam_slope_orig[i] - thick_slope_orig[i])^2)
                newthick[i] = thick_orig[i] + 1/sqrt(Re)*0.5*(th_u + th_l)
                newcam[i] = cam_orig[i] + 1/sqrt(Re)*0.5*(th_u - th_l)
            end
            newthick[1] = thick_orig[1]
            newcam[1] = cam_orig[1]
            
            i = surf.ndiv
            th_u = qc[2*surf.ndiv-2]*del_iter[2*surf.ndiv-2]/sqrt(1. + (cam_slope_orig[i] + thick_slope_orig[i])^2)
            th_l = qc[1]*del_iter[1]/sqrt(1. + (cam_slope_orig[i] - thick_slope_orig[i])^2)    
            newthick[i] = thick_orig[i] + 1/sqrt(Re)*0.5(th_u + th_l)
            newcam[i] = cam_orig[i] + 1/sqrt(Re)*0.5(th_u - th_l)
            
            bstart = [-0.1260; -0.3516; 0.2843; -0.1015]
            thickCoef = find_nacaThickCoef(surf, newthick, bstart)
            bstart = zeros(4)
            camCoef = find_nacaCamCoef(surf, newcam, bstart)
            
            th = parse(Int, surf.coord_file[7:8])/100.
            b1 = 0.2969
            bt = [b1; thickCoef]            
            @. nacath(x) = 5*th*(bt[1]*sqrt(x) + bt[2]*x + bt[3]*x^2 + bt[4]*x^3 + bt[5]*x^4)

            bc = camCoef
            @. nacacam(x) = bc[1]*x + bc[2]*x^2 + bc[3]*x^3 + bc[4]*x^4

            #Find new shape of airfoil
            for i = 1:surf.ndiv
                surf.thick[i] = nacath(surf.x[i])
                surf.thick_slope[i] = ForwardDiff.derivative(nacath, surf.x[i])
                surf.cam[i] = nacacam(surf.x[i])
                surf.cam_slope[i] = ForwardDiff.derivative(nacacam, surf.x[i])
            end
            surf.thick_slope[1] = thick_slope_orig[1]
            #plot(surf.x, surf.cam .+ surf.thick)
            #plot(surf.x, surf.cam .- surf.thick)

            #Check for convergence
            res =  sum(abs.(del_prev .- del_iter))
            println(iter, "   ", res)

            #if iter == iterMax
            if res <= resTol
                println("converged")
                del[:] = del_iter[:]
                E[:] = E_iter[:]
                push!(curfield.tev, TwoDVort(xloc_tev, zloc_tev, tevstr, vcore, 0., 0.))
                push!(curfield.lev, TwoDVort(xloc_lev, zloc_lev, levstr, vcore, 0., 0.))
            end

            # if iter == 3 && mod(istep,10) == 0
            #     figure(1)
            #     plot(surf.x, qu)
            #     figure(2)
            #     plot(surf.x, surf.thick)
            #     axis("equal")
            #     figure(3)
            #     plot(surf.x[2:end], delu)
            # end
            
        end
        

        #Update induced velocities to include effect of last shed vortex
        #add_indbound_lasttev(surf, curfield)
        update_indbound(surf, curfield)

        #Calculate adot
        update_atermdot(surf, dt)
        
        #Set previous values of aterm to be used for derivatives in next time step
        surf.ateprev[1] = surf.ate[1]
        for ia = 1:3
            surf.aprev[ia] = surf.aterm[ia]
        end
        
        #Calculate bound vortex strengths
        update_bv_src(surf)
        
        #Wake rollup
        wakeroll(surf, curfield, dt)

        qu, ql, phi_u, phi_l, cpu, cpl = calc_edgeVel_cp(surf, [curfield.u[1]; curfield.w[1]], phi_u, phi_l, dt)
        
        #println(cpu[end], " ", cpl[end], 0.5*(qu[end]^2 - ql[end]^2)*dt, " ", tevstr)
        
        #Force calculation
        cn, cs, cl, cd, cm = calc_forces(surf, cpu, cpl)
        
        bound_circ = phi_u[end-1] - phi_l[end-1]
        
        # write flow details if required
        if writeflag == 1
            if istep in writeArray
                dirname = "$(round(t,sigdigits=nround))"
                writeStamp(dirname, t, surf, curfield, qu, ql, cpu, cpl)
            end
        end
        vle = qu[1]
        
        stag = find_stag(surf, qu, ql)
        
        mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, vle,
                        cl, cd, cm, cn, cs, bound_circ, stag/surf.c, cpu[1]])

    end
    
    mat = mat'
    
    f = open("resultsSummary", "w")
    Serialization.serialize(f, ["#time \t", "alpha (rad) \t", "h/c \t", "u/uref \t", "Vle \t", "Cl \t", "Cd \t", "Cm \t", "Cn \t", "Cs \t", "bc \t", "xs \n"])
    writedlm(f, mat)
    close(f)
    
    mat, surf, curfield


    
    # # write flow details if required
    # if writeflag == 1
    #     if istep in writeArray
    #         dirname = "$(round(t,sigdigits=nround))"
    #         writeStamp(dirname, t, surf, curfield, qu, ql, cpu, cpl, suc, delu, Eu, thick_orig, quc, qucx, quct)
    #     end
    # end
    
end


