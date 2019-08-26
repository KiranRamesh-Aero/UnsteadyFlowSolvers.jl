function IBL_simul_iter_shed(surf::TwoDSurfThickBL, curfield::TwoDFlowField, ncell::Int64, nsteps::Int64 = 300, dtstar::Float64 = 0.015, startflag = 0, writeflag = 0, writeInterval = 1000., delvort = delNone(); maxwrite = 50, nround=6)

    
    
    int_wax = zeros(surf.ndiv)
    int_c = zeros(surf.ndiv)
    int_t = zeros(surf.ndiv)

    RHStransp = zeros(surf.ndiv*2-2)

    quc = zeros(surf.ndiv-1)
    suc = zeros(surf.ndiv-1)
    dsu = zeros(surf.ndiv-1)
    quxc = zeros(surf.ndiv-1)
    qutc = zeros(surf.ndiv-1)
    quc_prev = zeros(surf.ndiv-1)

    for i = 1:surf.ndiv-1
        suc[i] = (surf.su[i] + surf.su[i+1])/2
    end
    dsu[:] = diff(surf.su)

    wtu = zeros(surf.ndiv)
    wtl = zeros(surf.ndiv)
    
    t = 0.
    dt = dtstar

    vcore = 1.3*dt
    
    mat = zeros(0, 12)
    
    Re = 1000000

    q_u = zeros(surf.ndiv)
    q_l = zeros(surf.ndiv)

    i_xsl = surf.ndiv
    i_xsu = surf.ndiv

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

        #Solve the coupled problem
        if istep == 1
           soln = surf.LHS \ surf.RHS
            surf.aterm[:] = soln[1:surf.naterm]
            surf.bterm[:] = soln[surf.naterm+1:2*surf.naterm]
            surf.qu[:], surf.ql[:] = calc_edgeVel(surf, [curfield.u[1], curfield.w[1]])
            #smoothEnd!(surf.qu, 4)
            for i = 1:surf.ndiv-1
                quc_prev[i] = (surf.qu[i] + surf.qu[i+1])/2
            end
            
        else
            quc_prev[:] = quc[:]
        end
        
        if istep == 1
            swu = zeros(2)
            quw = zeros(2)
            swu[1] = surf.su[end]
            swu[2] = surf.su[end] + surf.uref*dt
            R = [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)]
            u_af = surf.kinem.u + curfield.u[1]
            w_af = -(surf.kinem.hdot - curfield.w[1])
            quw[1] = surf.qu[end]
            quw_x, quw_z = R*[u_af; w_af]
            quw[2] = sqrt(quw_x^2 + quw_z^2)
            delw, _, Ew, _ = initDelE(1)
        else
            n_w = minimum(length(curfield.lev), 50)
            swu = zeros(2+n_w)
            quw = zeros(2+n_w)
            swu[1] = surf.su[end]
            swu[2] = surf.su[end] + surf.uref*dt
            
            alpha = surf.kinem.alpha
            R = [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)]
            u_af = surf.kinem.u + curfield.u[1]
            w_af = -(surf.kinem.hdot - curfield.w[1])
            quw_x, quw_z = R*[u_af; w_af]
            quw[2] = sqrt(quw_x^2 + quw_z^2)
            
            for i = 1:n_w
                swu[i+2] = swu[i+1] + sqrt((curfield.lev[end-i].x - curfield.lev[end-i+1].x)^2 + (curfield.lev[end-i].z - curfield.lev[end-i+1].z)^2)
                quw_x, quw_z = R*[(curfield.lev[end-i].vx + u_af); (curfield.lev[end-i].vz + w_af)]
                quw[i+2] = sqrt(quw_x^2 + quw_z^2)
            end
            push!(delw, surf.delu[end])
            push!(Ew, surf.Eu[end])
        end
        dswu = zeros(length(swu)-1)
        swuc = zeros(length(swu)-1)
        quwc = zeros(length(quw)-1)
        
        for i = 1:length(swuc)
            swuc[i] = 0.5*(swu[i] + swu[i+1])
            dswu[i] = swu[i+1] - swu[i]
            quwc[i] = 0.5*(quw[i] + quw[i+1])
        end
        quwxc = diff(quwc)./dswu
        
        suvec = [surf.su; swu]
        dsuvec = [dsu; dswu]
        quvec = [quc; quwc]
        quxvec = [quxc; quwxc]
        qutvec = [qutc; zeros(length(quwc))]
        
        iter = 0
        res = 1
        del_iter = zeros(surf.ndiv-1 + length(delw))
        del_prev = zeros(surf.ndiv-1 + length(delw))
        soln = zeros(2*surf.naterm+1)
        E_iter = zeros(surf.ndiv-1 + length(Ew))
        RHS_v = zeros(2*surf.ndiv-2)


        
        while (res > 1e-5)
            
            iter += 1
            
            tev_us_str = 0.5*surf.qu[i_xsu]^2*dt
            
            xloc_u, zloc_u = calc_vloc_u(surf, dt, i_xsu)
            
            temptev1 = TwoDVort(xloc_u, zloc_u, tev_us_str, vcore, 0., 0.)
            
            RHS_v[2*surf.ndiv-3] = -100*tev_us_str/(surf.uref*surf.c)

            ind1_new_u_u, ind1_new_w_u = ind_vel([temptev1], surf.bnd_x_u, surf.bnd_z_u)
            ind1_new_u_l, ind1_new_w_l = ind_vel([temptev1], surf.bnd_x_l, surf.bnd_z_l)
            
            surf.uind_u[:] += ind1_new_u_u[:]
            surf.wind_u[:] += ind1_new_w_u[:]
            surf.uind_l[:] += ind1_new_u_l[:]
            surf.wind_l[:] += ind1_new_w_l[:]

            xloc_l, zloc_l = calc_vloc_l(surf, dt, i_xsu)
            
            update_thickLHS2V(surf, curfield, dt, vcore, xloc_l, zloc_l) 
            update_thickRHS(surf, curfield)
            
            soln[:] = surf.LHS[1:surf.ndiv*2-2, 1:surf.naterm*2+1] \ (surf.RHS[1:surf.ndiv*2-2] + RHStransp[:] + RHS_v[:])
            
            temptev2 = TwoDVort(xloc_tev, zloc_tev-0.01, soln[2*surf.naterm+1], vcore, 0., 0.)
            
            ind2_new_u_u, ind2_new_w_u = ind_vel([temptev2], surf.bnd_x_u, surf.bnd_z_u)
            ind2_new_u_l, ind2_new_w_l = ind_vel([temptev2], surf.bnd_x_l, surf.bnd_z_l)
            
            surf.uind_u[:] += ind2_new_u_u[:]
            surf.wind_u[:] += ind2_new_w_u[:]
            surf.uind_l[:] += ind2_new_u_l[:]
            surf.wind_l[:] += ind2_new_w_l[:]
            
            #Calculated edge velocity
            
            for i = 1:surf.ndiv
                l_x = 0; l_z = 0; t_x = 0; t_z = 0;
                for n = 1:surf.naterm
                    l_x += soln[n]*sin(n*surf.theta[i])
                    l_z += soln[n]*cos(n*surf.theta[i])
                    t_x -= soln[n+surf.naterm]*cos(n*surf.theta[i])
                    t_z += soln[n+surf.naterm]*sin(n*surf.theta[i])
                end
                l_x *= surf.uref
                l_z *= surf.uref
                t_x *= surf.uref
                t_z *= surf.uref
                
                w_x_u = surf.uind_u[i]*cos(surf.kinem.alpha) - surf.wind_u[i]*sin(surf.kinem.alpha)
                w_z_u = surf.uind_u[i]*sin(surf.kinem.alpha) + surf.wind_u[i]*cos(surf.kinem.alpha)
                w_x_l = surf.uind_l[i]*cos(surf.kinem.alpha) - surf.wind_l[i]*sin(surf.kinem.alpha)
                w_z_l = surf.uind_l[i]*sin(surf.kinem.alpha) + surf.wind_l[i]*cos(surf.kinem.alpha)
                
                vels = [curfield.u[1], curfield.w[1]]
                
                utot_u = (surf.kinem.u + vels[1])*cos(surf.kinem.alpha) + (surf.kinem.hdot - vels[2])*sin(surf.kinem.alpha) - surf.kinem.alphadot*(surf.cam[i] + surf.thick[i]) + w_x_u + l_x + t_x
                wtot_u = (surf.kinem.u + vels[1])*sin(surf.kinem.alpha) + (surf.kinem.hdot - vels[2])*cos(surf.kinem.alpha) + surf.kinem.alphadot*(surf.x[i] - surf.pvt*surf.c) + w_z_u + l_z + t_z
                utot_l = (surf.kinem.u + vels[1])*cos(surf.kinem.alpha) + (surf.kinem.hdot - vels[2])*sin(surf.kinem.alpha) - surf.kinem.alphadot*(surf.cam[i] - surf.thick[i]) + w_x_l - l_x + t_x
                wtot_l = (surf.kinem.u + vels[1])*sin(surf.kinem.alpha) + (surf.kinem.hdot - vels[2])*cos(surf.kinem.alpha) + surf.kinem.alphadot*(surf.x[i] - surf.pvt*surf.c) + w_z_l + l_z - t_z
                
                q_u[i] = (1. /sqrt(1. + (surf.cam_slope[i] + surf.thick_slope[i])^2))*(utot_u + (surf.cam_slope[i] + surf.thick_slope[i])*wtot_u)
                q_l[i] = (1. /sqrt(1. + (surf.cam_slope[i] - surf.thick_slope[i])^2))*(utot_l + (surf.cam_slope[i] - surf.thick_slope[i])*wtot_l)
            end
            
            #Remove the infulence of last TEV
            surf.uind_u[:] -= ind1_new_u_u[:]
            surf.wind_u[:] -= ind1_new_w_u[:]
            surf.uind_l[:] -= ind1_new_u_l[:]
            surf.wind_l[:] -= ind1_new_w_l[:]

            surf.uind_u[:] -= ind2_new_u_u[:]
            surf.wind_u[:] -= ind2_new_w_u[:]
            surf.uind_l[:] -= ind2_new_u_l[:]
            surf.wind_l[:] -= ind2_new_w_l[:]
                  
            #Solve the FV problem at cell centres
            for i = 1:surf.ndiv-1
                quc[i] = (q_u[i] + q_u[i+1])/2
            end
            
            quxc[2:end] = diff(quc)./diff(suc)
            quxc[1] = quxc[2]
            qutc[:] = (quc[:] - quc_prev[:])/dt
            
            quvec = [quc; quwc]
            quxvec = [quxc; quwxc]
            qutvec = [qutc; zeros(length(quwc))]
            
            deltot = [surf.delu; delw]
            Etot = [surf.Eu; Ew]
            w = [deltot deltot.*(Etot .+ 1)]
            
            wsoln, i_xsu = FVMIBLgridvar(w, quc, qutc, quxc, dsu, t, t+dt) 
            if i_xsu == 0
                i_xsu = surf.ndiv
            end
            
            del_prev[:] = del_iter[:]
            
            del_iter[:] = wsoln[:,1]
            E_iter[:] = wsoln[:,2]./wsoln[:,1] .- 1.
            
            #smoothEdges!(del_iter, 20)
            #smoothEdges!(E_iter, 20)

            #del_iter[end-10:end] = ones(11)*del_iter[end-11]

            wtu[2:end] = (1/sqrt(Re))*diff(q_u.*del_iter[1:surf.ndiv])./diff(surf.x)
            wtu[1] = wtu[2]

            #smoothEdges!(wtu, 20)
            
            wtl[:] = -wtu[:]

            RHStransp[:] .= 0.
            
            #Add transpiration velocity to RHS
            for i = 2:surf.ndiv-1
                RHStransp[i-1] = 0.5*(wtu[i] + wtl[i])
                RHStransp[surf.ndiv+i-3] = 0.5*(wtu[i] - wtl[i])
            end
            
            res =  sum(abs.(del_prev .- del_iter))

            println(iter, "   ", res)
            
        end

        
        #Assign the solution
        surf.aterm[:] = soln[1:surf.naterm]
        surf.bterm[:] = soln[surf.naterm+1:2*surf.naterm]
        tevstr = soln[2*surf.naterm+1]
        push!(curfield.tev, TwoDVort(xloc_tev, zloc_tev-0.01, tevstr, vcore, 0., 0.))
        push!(curfield.tev, TwoDVort(xloc_tev, zloc_tev+0.01, 0.5*surf.qu[end]^2*dt, vcore, 0., 0.))
        surf.delu[:] = del_iter[:]
        surf.Eu[:] = E_iter[:]

        surf.qu[:], surf.ql[:] = calc_edgeVel(surf, [curfield.u[1], curfield.w[1]])

        smoothEnd!(surf.qu, 4)

        #Calculate adot
        update_atermdot(surf, dt)
        
         #Set previous values of aterm to be used for derivatives in next time step
         surf.a0prev[1] = surf.a0[1]
         for ia = 1:3
             surf.aprev[ia] = surf.aterm[ia]
         end

         #Calculate bound vortex strengths
         update_bv_src(surf)

         #Add effect of transpiration to sources and vortices

         #Wake rollup
         wakeroll(surf, curfield, dt)

         #Force calculation
        cnc, cnnc, cn, cs, cl, cd, int_wax, int_c, int_t = calc_forces(surf, int_wax, int_c, int_t, dt)

        vle = 0.
        stag = 0.
        
        #mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, vle,
    #                    cl, cd, cnc, cnnc, cn, cs, stag])
        
        println("here")
    end
    
    #mat = mat'

    #f = open("resultsSummary", "w")
    #Serialization.serialize(f, ["#time \t", "alpha (deg) \t", "h/c \t", "u/uref \t", "A0 \t", "Cl \t", "Cd \t", "Cm \n"])
    #writedlm(f, mat)
    #close(f)
    
    return mat, surf, curfield

end


function IBL_simul(surf::TwoDSurfThickBL, curfield::TwoDFlowField, ncell::Int64, nsteps::Int64 = 300, dtstar::Float64 = 0.015, startflag = 0, writeflag = 0, writeInterval = 1000., delvort = delNone(); maxwrite = 50, nround=6)

    vcore = 0.02*surf.c
    
    int_wax = zeros(surf.ndiv)
    int_c = zeros(surf.ndiv)
    int_t = zeros(surf.ndiv)

    RHStransp = zeros(surf.ndiv*2-2)

    quc = zeros(surf.ndiv-1)
    suc = zeros(surf.ndiv-1)
    dsu = zeros(surf.ndiv-1)
    quxc = zeros(surf.ndiv-1)
    qutc = zeros(surf.ndiv-1)
    quc_prev = zeros(surf.ndiv-1)

    for i = 1:surf.ndiv-1
        suc[i] = (surf.su[i] + surf.su[i+1])/2
    end
    dsu[:] = diff(surf.su)

    wtu = zeros(surf.ndiv)
    wtl = zeros(surf.ndiv)
    
    t = 0.
    dt = dtstar

    mat = zeros(0, 12)
    
    Re = 10000
    
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

        #Set up the matrix problem
        surf, xloc_tev, zloc_tev = update_thickLHS(surf, curfield, dt, vcore)

        #Construct RHS vector
        update_thickRHS(surf, curfield)

        #Solve the coupled problem
        if istep == 1
            soln = surf.LHS \ surf.RHS
            surf.aterm[:] = soln[1:surf.naterm]
            surf.bterm[:] = soln[surf.naterm+1:2*surf.naterm]
                                            
            surf.qu[:], surf.ql[:] = calc_edgeVel(surf, [curfield.u[1], curfield.w[1]])
            
            for i = 1:surf.ndiv-1
                quc_prev[i] = (surf.qu[i] + surf.qu[i+1])/2
            end
        else
            quc_prev[:] = quc[:]
        end
        
        function interaction_simul(x)
            
            #x vector = [aterm; bterm; tevstr; delu, Eu]

            FC = zeros(4*surf.ndiv-4)
            
            FC[1:surf.ndiv*2-2] = surf.LHS[1:surf.ndiv*2-2, 1:surf.naterm*2+1]*x[1:surf.naterm*2+1] - (surf.RHS[1:surf.ndiv*2-2] + RHStransp[:])
           
                       
            #Add influence of last shed vortex to calculate Ue
            
            temptev = TwoDVort(xloc_tev, zloc_tev, x[2*surf.naterm+1], vcore, 0., 0.)
            
            ind_new_u_u, ind_new_w_u = ind_vel([temptev], surf.bnd_x_u, surf.bnd_z_u)
            ind_new_u_l, ind_new_w_l = ind_vel([temptev], surf.bnd_x_l, surf.bnd_z_l)
            
            surf.uind_u[:] += ind_new_u_u[:]
            surf.wind_u[:] += ind_new_w_u[:]
            surf.uind_l[:] += ind_new_u_l[:]
            surf.wind_l[:] += ind_new_w_l[:]
            
            #Calculated edge velocity
            q_u = zeros(surf.ndiv)
            q_l = zeros(surf.ndiv)
            
            for i = 1:surf.ndiv
                
                l_x = 0; l_z = 0; t_x = 0; t_z = 0;
                for n = 1:surf.naterm
                    l_x += x[n]*sin(n*surf.theta[i])
                    l_z += x[n]*cos(n*surf.theta[i])
                    t_x -= x[n+surf.naterm]*cos(n*surf.theta[i])
                    t_z += x[n+surf.naterm]*sin(n*surf.theta[i])
                end
                l_x *= surf.uref
                l_z *= surf.uref
                t_x *= surf.uref
                t_z *= surf.uref
                
                w_x_u = surf.uind_u[i]*cos(surf.kinem.alpha) - surf.wind_u[i]*sin(surf.kinem.alpha)
                w_z_u = surf.uind_u[i]*sin(surf.kinem.alpha) + surf.wind_u[i]*cos(surf.kinem.alpha)
                w_x_l = surf.uind_l[i]*cos(surf.kinem.alpha) - surf.wind_l[i]*sin(surf.kinem.alpha)
                w_z_l = surf.uind_l[i]*sin(surf.kinem.alpha) + surf.wind_l[i]*cos(surf.kinem.alpha)

                vels = [curfield.u[1], curfield.w[1]]
                
                utot_u = (surf.kinem.u + vels[1])*cos(surf.kinem.alpha) + (surf.kinem.hdot - vels[2])*sin(surf.kinem.alpha) - surf.kinem.alphadot*(surf.cam[i] + surf.thick[i]) + w_x_u + l_x + t_x
                wtot_u = (surf.kinem.u + vels[1])*sin(surf.kinem.alpha) + (surf.kinem.hdot - vels[2])*cos(surf.kinem.alpha) + surf.kinem.alphadot*(surf.x[i] - surf.pvt*surf.c) + w_z_u + l_z + t_z
                utot_l = (surf.kinem.u + vels[1])*cos(surf.kinem.alpha) + (surf.kinem.hdot - vels[2])*sin(surf.kinem.alpha) - surf.kinem.alphadot*(surf.cam[i] - surf.thick[i]) + w_x_l - l_x + t_x
                wtot_l = (surf.kinem.u + vels[1])*sin(surf.kinem.alpha) + (surf.kinem.hdot - vels[2])*cos(surf.kinem.alpha) + surf.kinem.alphadot*(surf.x[i] - surf.pvt*surf.c) + w_z_l + l_z - t_z
                
                q_u[i] = (1. /sqrt(1. + (surf.cam_slope[i] + surf.thick_slope[i])^2))*(utot_u + (surf.cam_slope[i] + surf.thick_slope[i])*wtot_u)
                q_l[i] = (1. /sqrt(1. + (surf.cam_slope[i] - surf.thick_slope[i])^2))*(utot_l + (surf.cam_slope[i] - surf.thick_slope[i])*wtot_l)
            end
            #q_u, q_l = calc_edgeVel(surf, [curfield.u[1], curfield.w[1]])

            #Remove the infulence of last TEV
            surf.uind_u[:] -= ind_new_u_u[:]
            surf.wind_u[:] -= ind_new_w_u[:]
            surf.uind_l[:] -= ind_new_u_l[:]
            surf.wind_l[:] -= ind_new_w_l[:]

            
            #Solve the FV problem at cell centres
            for i = 1:surf.ndiv-1
                quc[i] = (q_u[i] + q_u[i+1])/2
            end
            
            quxc[2:end] = diff(quc)./diff(suc)
            quxc[1] = quxc[2]
            
            qutc[:] = (quc[:] - quc_prev[:])/dt
            
            #qutc[:] .=  0.

            w = [surf.delu surf.delu[:].*(surf.Eu[:] .+ 1)]
            
            wupd = FVMIBLgridvar(w, quc, qutc, quxc, dsu, t, t+dt)
            
            x_w = [x[2*surf.naterm+2:2*surf.naterm+surf.ndiv] x[2*surf.naterm+2:2*surf.naterm+surf.ndiv].*(x[2*surf.naterm+surf.ndiv+1:2*surf.naterm+2*surf.ndiv-1] .+ 1)]
            res_bl = x_w - wupd
            
            FC[2*surf.ndiv-1:3*surf.ndiv-3] = res_bl[:,1]
            FC[3*surf.ndiv-2:4*surf.ndiv-4] = res_bl[:,2]

            delvec = [x[2*surf.naterm+2:2*surf.naterm+surf.ndiv]; x[2*surf.naterm+surf.ndiv]]
            
            wtu[2:end] = (1/sqrt(Re))*diff(q_u.*delvec)./diff(surf.x)
            wtu[1] = wtu[2]
            wtl[:] = -wtu[:]

            RHStransp[:] .= 0.
            
            #Add transpiration velocity to RHS
            for i = 2:surf.ndiv-1
                RHStransp[i-1] = 0.5*(wtu[i] + wtl[i])
                RHStransp[surf.ndiv+i-3] = 0.5*(wtu[i] - wtl[i])
            end
            
            x[2*surf.naterm+2*surf.ndiv-1:4*surf.ndiv-4] .= 0.
            #return sum(FC.^2)
            return FC
        end

        xstart = [surf.aterm; surf.bterm; -0.01; surf.delu; surf.Eu; zeros(5)]
        #x0 = [surf.aterm; surf.bterm; -0.001; surf.delu; surf.Eu]
        
        #x1 = [surf.aterm .+ 1e-6.*rand(surf.naterm); surf.bterm .+ 1e-6.*rand(surf.naterm); -0.0011; surf.delu .+ 1e-6.*rand(surf.ndiv-1);surf.Eu .- 1e-6.*rand(surf.ndiv-1)]
        
        soln = nlsolve(interaction_simul, xstart, xtol = 1e-6, show_trace=true, iterations=3, method=:newton)
        #optimize!(LeastSquaresProblem(x=soln, f! = interaction_simul!, output_length=4*surf.ndiv-4, autodiff=:central))

        #soln = optimize!(LeastSquaresProblem(x = x0, f! = interaction_simul!, output_length = length(x0)))
        #soln = LeastSquaresOptim.optimize(interaction_simul, x0)
                        
                
        soln = soln.zero
        
        #Assign the solution
        # surf.aterm[:] = soln.zero[1:surf.naterm]
        # surf.bterm[:] = soln.zero[surf.naterm+1:2*surf.naterm]
        # tevstr = soln.zero[2*surf.naterm+1]
        # push!(curfield.tev, TwoDVort(xloc_tev, zloc_tev, tevstr, vcore, 0., 0.))
        # surf.delu[:] = soln.zero[2*surf.naterm+2:2*surf.naterm+surf.ndiv]
        # surf.Eu[:] = soln.zero[2*surf.naterm+surf.ndiv+1:2*surf.naterm+2*surf.ndiv-1]

        surf.aterm[:] = soln[1:surf.naterm]
        surf.bterm[:] = soln[surf.naterm+1:2*surf.naterm]
        tevstr = soln[2*surf.naterm+1]
        push!(curfield.tev, TwoDVort(xloc_tev, zloc_tev, tevstr, vcore, 0., 0.))
        surf.delu[:] = soln[2*surf.naterm+2:2*surf.naterm+surf.ndiv]
        surf.Eu[:] = soln[2*surf.naterm+surf.ndiv+1:2*surf.naterm+2*surf.ndiv-1]

        surf.qu[:], surf.ql[:] = calc_edgeVel(surf, [curfield.u[1], curfield.w[1]])

        smoothEdges!(surf.delu, 10)
        smoothEdges!(surf.Eu, 10)
        smoothEnd!(surf.qu, 10)
        
        #Calculate adot
        update_atermdot(surf, dt)
        
         #Set previous values of aterm to be used for derivatives in next time step
         surf.a0prev[1] = surf.a0[1]
         for ia = 1:3
             surf.aprev[ia] = surf.aterm[ia]
         end

         #Calculate bound vortex strengths
         update_bv_src(surf)

         #Add effect of transpiration to sources and vortices

         #Wake rollup
         wakeroll(surf, curfield, dt)

         #Force calculation
        cnc, cnnc, cn, cs, cl, cd, int_wax, int_c, int_t = calc_forces(surf, int_wax, int_c, int_t, dt)

        vle = 0.
        stag = 0.
        
        #mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, vle,
    #                    cl, cd, cnc, cnnc, cn, cs, stag])
        
        println("here")
    end
    
    #mat = mat'

    #f = open("resultsSummary", "w")
    #Serialization.serialize(f, ["#time \t", "alpha (deg) \t", "h/c \t", "u/uref \t", "A0 \t", "Cl \t", "Cd \t", "Cm \n"])
    #writedlm(f, mat)
    #close(f)
    
    return mat, surf, curfield

end


function IBL_simul_iter(surf::TwoDSurfThickBL, curfield::TwoDFlowField, ncell::Int64, nsteps::Int64 = 300, dtstar::Float64 = 0.015, startflag = 0, writeflag = 0, writeInterval = 1000., delvort = delNone(); maxwrite = 50, nround=6)

    vcore = 0.02*surf.c
    
    int_wax = zeros(surf.ndiv)
    int_c = zeros(surf.ndiv)
    int_t = zeros(surf.ndiv)

    RHStransp = zeros(surf.ndiv*2-2)

    quc = zeros(surf.ndiv-1)
    suc = zeros(surf.ndiv-1)
    dsu = zeros(surf.ndiv-1)
    quxc = zeros(surf.ndiv-1)
    qutc = zeros(surf.ndiv-1)
    quc_prev = zeros(surf.ndiv-1)

    for i = 1:surf.ndiv-1
        suc[i] = (surf.su[i] + surf.su[i+1])/2
    end
    dsu[:] = diff(surf.su)

    wtu = zeros(surf.ndiv)
    wtl = zeros(surf.ndiv)
    
    t = 0.
    dt = dtstar

    mat = zeros(0, 12)
    
    Re = 10000
    
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

        #Set up the matrix problem
        surf, xloc_tev, zloc_tev = update_thickLHS(surf, curfield, dt, vcore)

        #Construct RHS vector
        update_thickRHS(surf, curfield)

        #Solve the coupled problem
        if istep == 1
           soln = surf.LHS \ surf.RHS
            surf.aterm[:] = soln[1:surf.naterm]
            surf.bterm[:] = soln[surf.naterm+1:2*surf.naterm]
            
            surf.qu[:], surf.ql[:] = calc_edgeVel(surf, [curfield.u[1], curfield.w[1]])
            smoothEnd!(surf.qu, 4)
            for i = 1:surf.ndiv-1
                quc_prev[i] = (surf.qu[i] + surf.qu[i+1])/2
            end
            
        else
            quc_prev[:] = quc[:]
        end

        iter = 0
        res = 1
        del_iter = zeros(surf.ndiv-1)
        del_prev = zeros(surf.ndiv-1)
        soln = zeros(2*surf.naterm+1)
        E_iter = zeros(surf.ndiv-1)
        while (res > 1e-5)

            iter += 1
            
            soln[:] = surf.LHS[1:surf.ndiv*2-2, 1:surf.naterm*2+1] \ (surf.RHS[1:surf.ndiv*2-2] + RHStransp[:])
            
            temptev = TwoDVort(xloc_tev, zloc_tev, soln[2*surf.naterm+1], vcore, 0., 0.)
            
            ind_new_u_u, ind_new_w_u = ind_vel([temptev], surf.bnd_x_u, surf.bnd_z_u)
            ind_new_u_l, ind_new_w_l = ind_vel([temptev], surf.bnd_x_l, surf.bnd_z_l)
            
            surf.uind_u[:] += ind_new_u_u[:]
            surf.wind_u[:] += ind_new_w_u[:]
            surf.uind_l[:] += ind_new_u_l[:]
            surf.wind_l[:] += ind_new_w_l[:]
            
            #Calculated edge velocity
            q_u = zeros(surf.ndiv)
            q_l = zeros(surf.ndiv)
            
            for i = 1:surf.ndiv
                l_x = 0; l_z = 0; t_x = 0; t_z = 0;
                for n = 1:surf.naterm
                    l_x += soln[n]*sin(n*surf.theta[i])
                    l_z += soln[n]*cos(n*surf.theta[i])
                    t_x -= soln[n+surf.naterm]*cos(n*surf.theta[i])
                    t_z += soln[n+surf.naterm]*sin(n*surf.theta[i])
                end
                l_x *= surf.uref
                l_z *= surf.uref
                t_x *= surf.uref
                t_z *= surf.uref
                
                w_x_u = surf.uind_u[i]*cos(surf.kinem.alpha) - surf.wind_u[i]*sin(surf.kinem.alpha)
                w_z_u = surf.uind_u[i]*sin(surf.kinem.alpha) + surf.wind_u[i]*cos(surf.kinem.alpha)
                w_x_l = surf.uind_l[i]*cos(surf.kinem.alpha) - surf.wind_l[i]*sin(surf.kinem.alpha)
                w_z_l = surf.uind_l[i]*sin(surf.kinem.alpha) + surf.wind_l[i]*cos(surf.kinem.alpha)
                
                vels = [curfield.u[1], curfield.w[1]]
                
                utot_u = (surf.kinem.u + vels[1])*cos(surf.kinem.alpha) + (surf.kinem.hdot - vels[2])*sin(surf.kinem.alpha) - surf.kinem.alphadot*(surf.cam[i] + surf.thick[i]) + w_x_u + l_x + t_x
                wtot_u = (surf.kinem.u + vels[1])*sin(surf.kinem.alpha) + (surf.kinem.hdot - vels[2])*cos(surf.kinem.alpha) + surf.kinem.alphadot*(surf.x[i] - surf.pvt*surf.c) + w_z_u + l_z + t_z
                utot_l = (surf.kinem.u + vels[1])*cos(surf.kinem.alpha) + (surf.kinem.hdot - vels[2])*sin(surf.kinem.alpha) - surf.kinem.alphadot*(surf.cam[i] - surf.thick[i]) + w_x_l - l_x + t_x
                wtot_l = (surf.kinem.u + vels[1])*sin(surf.kinem.alpha) + (surf.kinem.hdot - vels[2])*cos(surf.kinem.alpha) + surf.kinem.alphadot*(surf.x[i] - surf.pvt*surf.c) + w_z_l + l_z - t_z
                
                q_u[i] = (1. /sqrt(1. + (surf.cam_slope[i] + surf.thick_slope[i])^2))*(utot_u + (surf.cam_slope[i] + surf.thick_slope[i])*wtot_u)
                q_l[i] = (1. /sqrt(1. + (surf.cam_slope[i] - surf.thick_slope[i])^2))*(utot_l + (surf.cam_slope[i] - surf.thick_slope[i])*wtot_l)
            end

            #Remove the infulence of last TEV
            surf.uind_u[:] -= ind_new_u_u[:]
            surf.wind_u[:] -= ind_new_w_u[:]
            surf.uind_l[:] -= ind_new_u_l[:]
            surf.wind_l[:] -= ind_new_w_l[:]

            #Solve the FV problem at cell centres
            for i = 1:surf.ndiv-1
                quc[i] = (q_u[i] + q_u[i+1])/2
            end

            smoothEnd!(quc, 4)
            
            quxc[2:end] = diff(quc)./diff(suc)
            quxc[1] = quxc[2]

                       
            qutc[:] = (quc[:] - quc_prev[:])/dt
            
            #    qutc[:] .= 0.
            
            w = [surf.delu surf.delu.*(surf.Eu .+ 1)]

            wsoln = FVMIBLgridvar(w, quc, qutc, quxc, dsu, t, t+dt) 
            
            del_prev[:] = del_iter[:]
            
            del_iter[:] = wsoln[:,1]
            E_iter[:] = wsoln[:,2]./wsoln[:,1] .- 1.
            
            #smoothEdges!(del_iter, 20)
            #smoothEdges!(E_iter, 20)

            #del_iter[end-10:end] = ones(11)*del_iter[end-11]
            
            wtu[2:end] = (1/sqrt(Re))*diff(q_u.*[del_iter; del_iter[end]])./diff(surf.x)
            wtu[1] = wtu[2]

            #smoothEdges!(wtu, 20)
            
            wtl[:] = -wtu[:]

            RHStransp[:] .= 0.
            
            #Add transpiration velocity to RHS
            for i = 2:surf.ndiv-1
                RHStransp[i-1] = 0.5*(wtu[i] + wtl[i])
                RHStransp[surf.ndiv+i-3] = 0.5*(wtu[i] - wtl[i])
            end
            
            res =  sum(abs.(del_prev .- del_iter))

            println(iter, "   ", res)
            
        end

        
        #Assign the solution
        surf.aterm[:] = soln[1:surf.naterm]
        surf.bterm[:] = soln[surf.naterm+1:2*surf.naterm]
        tevstr = soln[2*surf.naterm+1]
        push!(curfield.tev, TwoDVort(xloc_tev, zloc_tev, tevstr, vcore, 0., 0.))
        surf.delu[:] = del_iter[:]
        surf.Eu[:] = E_iter[:]

        surf.qu[:], surf.ql[:] = calc_edgeVel(surf, [curfield.u[1], curfield.w[1]])

        smoothEnd!(surf.qu, 4)

        #Calculate adot
        update_atermdot(surf, dt)
        
         #Set previous values of aterm to be used for derivatives in next time step
         surf.a0prev[1] = surf.a0[1]
         for ia = 1:3
             surf.aprev[ia] = surf.aterm[ia]
         end

         #Calculate bound vortex strengths
         update_bv_src(surf)

         #Add effect of transpiration to sources and vortices

         #Wake rollup
         wakeroll(surf, curfield, dt)

         #Force calculation
        cnc, cnnc, cn, cs, cl, cd, int_wax, int_c, int_t = calc_forces(surf, int_wax, int_c, int_t, dt)

        vle = 0.
        stag = 0.
        
        #mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, vle,
    #                    cl, cd, cnc, cnnc, cn, cs, stag])
        
        println("here")
    end
    
    #mat = mat'

    #f = open("resultsSummary", "w")
    #Serialization.serialize(f, ["#time \t", "alpha (deg) \t", "h/c \t", "u/uref \t", "A0 \t", "Cl \t", "Cd \t", "Cm \n"])
    #writedlm(f, mat)
    #close(f)
    
    return mat, surf, curfield

end



    # function transpTogether(surf::TwoDSurfThickBL, curfield::TwoDFlowField, ncell::Int64, nsteps::Int64 = 300, dtstar::Float64 = 0.015, startflag = 0, writeflag = 0, writeInterval = 1000., delvort = delNone(); maxwrite = 50, nround=6)

#     # If a restart directory is provided, read in the simulation data
#     if startflag == 0
#         mat = zeros(0, 12)
#         t = 0.

#     elseif startflag == 1
#         dirvec = readdir()
#         dirresults = map(x->(v = tryparse(Float64,x); typeof(v) == Nothing ? 0.0 : v),dirvec)
#         latestTime = maximum(dirresults)
#         mat = DelimitedFiles.readdlm("resultsSummary")
#         t = mat[end,1]
#     else
#         throw("invalid start flag, should be 0 or 1")
#     end
#     mat = mat'

#     dt = dtstar*surf.c/surf.uref

#     Re = 1000000


#     # if writeflag is on, determine the timesteps to write at
#     if writeflag == 1
#         writeArray = Int64[]
#         tTot = nsteps*dt
#         for i = 1:maxwrite
#             tcur = writeInterval*real(i)
#             if t > tTot
#                 break
#             else
#                 push!(writeArray, Int(round(tcur/dt)))
#             end
#         end
#     end

#     vcore = 0.02*surf.c

#     int_wax = zeros(surf.ndiv)
#     int_c = zeros(surf.ndiv)
#     int_t = zeros(surf.ndiv)

#     quprevc = zeros(surf.ndiv-1)

#     for istep = 1:nsteps

#         t = t + dt

#         #Update kinematic parameters
#         update_kinem(surf, t)

#         #Update flow field parameters if any
#         update_externalvel(curfield, t)

#         #Update bound vortex positions
#         update_boundpos(surf, dt)

#         #Update incduced velocities on airfoil
#         update_indbound(surf, curfield)
    
#         #Set up the matrix problem
#         surf, xloc_tev, zloc_tev = update_thickLHS(surf, curfield, dt, vcore)

#         #Construct RHS vector
#         update_thickRHS(surf, curfield)

#         res = 1.

#         RHStransp = zeros(surf.ndiv*2-2)
        
#         iter_delu = zeros(surf.ndiv-1)
#         iter_Eu = zeros(surf.ndiv-1)
#         wtu = zeros(surf.ndiv)
#         wtl = zeros(surf.ndiv)

#         quc = zeros(surf.ndiv-1)
        
        
#         #Iterate for viscous solution and interaction

#         iter = 0

#         while (iter < 2 || res > 5e-3)
            
#             iter += 1

#             println("res")

#             if iter > 1
#                 pop!(curfield.tev)
#             end

#             #Update induced velocities to include effect of last shed vortex
#             update_indbound(surf, curfield)

#             soln = surf.LHS[1:surf.ndiv*2-2, 1:surf.naterm*2+1] \ (surf.RHS[1:surf.ndiv*2-2] + RHStransp[:])

#             res = sqrt(sum((soln[1:end-1] .- [surf.aterm; surf.bterm]).^2))

#             #Assign the solution
#             for i = 1:surf.naterm
#                 surf.aterm[i] = soln[i]
#                 surf.bterm[i] = soln[i+surf.naterm]
#             end
#             tevstr = soln[2*surf.naterm+1]*surf.uref*surf.c

#             push!(curfield.tev, TwoDVort(xloc_tev, zloc_tev, tevstr, vcore, 0., 0.))
#             #Update induced velocities to include effect of last shed vortex
#             update_indbound(surf, curfield)

#             #Dont update qu during iteration - just the inviscid value.
#             #This basically means no strong coupling.
    #             #if iter == 1
#                 surf.qu[:], surf.ql[:] = calc_edgeVel(surf, [curfield.u[1], curfield.w[1]])
#             #end

#             quxc = zeros(surf.ndiv-1)
#             qutc = zeros(surf.ndiv-1)
#             dsu = zeros(surf.ndiv-1)
#             suc = zeros(surf.ndiv-1)
            
    #                 #Solve the FV problem at cell centres
    #             for i = 1:surf.ndiv-1
#                 quc[i] = (surf.qu[i] + surf.qu[i+1])/2
#                 suc[i] = (surf.su[i] + surf.su[i+1])/2
#             end
    
#             quxc[2:end] = diff(quc)./diff(suc)
#             quxc[1] = 2*quxc[2] - quxc[3]
            
#             smoothEdges!(quxc, 5)
            
#             if istep == 1
#                 quprevc[:] = quc[:]
#             end

#             qutc[:] .= (quc[:] .- quprevc[:])./dt

#             dsu = diff(surf.su)
            
#             w0 = [surf.delu surf.delu.*(surf.Eu .+ 1)]

#             w0 = FVMIBLgridvar(w0, quc, qutc, quxc, dsu, t-dt, t)
            
#             iter_delu[:] = w0[:,1]
#             iter_Eu[:] = w0[:,2]./w0[:,1] .- 1.

#             #println(iter_delu)
            
#             #plot(surf.x, iter_delu)
            
#             smoothEdges!(iter_delu, 5)
            
#             wtu[2:end] = (1/sqrt(Re))*diff(surf.qu.*[iter_delu; iter_delu[end]])./diff(surf.x)
#             wtu[1] = wtu[2]

#             # smoothEdges!(wtu, 10)

#              wtl[:] = -wtu[:]
            
#             if istep ==80 && iter == 2
#                 figure(1)
#                 plot(suc, iter_delu)
#                 figure(2)
#                 plot(suc, quxc)
#                 figure(3)
#                 plot(surf.x, wtu)
#                 error("here")
#             end

#             RHStransp[:] .= 0.

#             #Add transpiration velocity to RHS
#             for i = 2:surf.ndiv-1
#                 RHStransp[i-1] = 0.5*(wtu[i] + wtl[i])
#                 RHStransp[surf.ndiv+i-3] = 0.5*(wtu[i] - wtl[i])
#             end

#             println(istep, "   ", res)
#         end

#         quprevc[:] = quc[:]

#         #Assign bl
#         surf.delu[:] = iter_delu[:]
#         #surf.delu[end] = surf.delu[end-1]
#         #surf.Eu[end] = surf.Eu[end-1]

#         surf.Eu[:] = iter_Eu[:]

#         #Calculate adot
#         update_atermdot(surf, dt)

#         #Set previous values of aterm to be used for derivatives in next time step
#         surf.a0prev[1] = surf.a0[1]
#         for ia = 1:3
#             surf.aprev[ia] = surf.aterm[ia]
#         end

#         #Calculate bound vortex strengths
#         update_bv_src(surf)

#         #Add effect of transpiration to sources and vortices

#         #Wake rollup
#         wakeroll(surf, curfield, dt)

#         #Force calculation
#         cnc, cnnc, cn, cs, cl, cd, int_wax, int_c, int_t = calc_forces(surf, int_wax, int_c, int_t, dt)

#         vle = surf.qu[1]

#         if vle > 0.
#             qspl = Spline1D(surf.x, surf.ql)
#             stag = try
#                 roots(qspl, maxn=1)[1]
# OA            catch
#                 0.
#             end
#         else
#             qspl = Spline1D(surf.x, surf.qu)
#             stag = try
#                 roots(qspl, maxn=1)[1]
#             catch
#                 0.
#             end
#         end

#         mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, vle,
#         cl, cd, cnc, cnnc, cn, cs, stag])

#         println("here")
#     end

#     mat = mat'



#     f = open("resultsSummary", "w")
#     Serialization.serialize(f, ["#time \t", "alpha (deg) \t", "h/c \t", "u/uref \t", "A0 \t", "Cl \t", "Cd \t", "Cm \n"])
#     writedlm(f, mat)
#     close(f)

#     return mat, surf, curfield

# end


# function IBLshed(surf::TwoDSurfThickBL, curfield::TwoDFlowField, ncell::Int64, nsteps::Int64 = 300, dtstar::Float64 = 0.015, startflag = 0, writeflag = 0, writeInterval = 1000., delvort = delNone(); maxwrite = 50, nround=6)

#     # If a restart directory is provided, read in the simulation data
#     if startflag == 0
#         mat = zeros(0, 12)
#         t = 0.

#     elseif startflag == 1
#         dirvec = readdir()
#         dirresults = map(x->(v = tryparse(Float64,x); typeof(v) == Nothing ? 0.0 : v),dirvec)
#         latestTime = maximum(dirresults)
#         mat = DelimitedFiles.readdlm("resultsSummary")
#         t = mat[end,1]
#     else
#         throw("invalid start flag, should be 0 or 1")
#     end
#     mat = mat'

#     dt = dtstar*surf.c/surf.uref

#     Re = 10000


#     # if writeflag is on, determine the timesteps to write at
#     if writeflag == 1
#         writeArray = Int64[]
#         tTot = nsteps*dt
#         for i = 1:maxwrite
#             tcur = writeInterval*real(i)
#             if t > tTot
#                 break
#             else
#                 push!(writeArray, Int(round(tcur/dt)))
#             end
#         end
#     end

#     vcore = 0.02*surf.c

#     int_wax = zeros(surf.ndiv)
#     int_c = zeros(surf.ndiv)
#     int_t = zeros(surf.ndiv)

#     quprevc = zeros(surf.ndiv-1)

#     i_xsu = surf.ndiv-1
#     i_xsl = surf.ndiv-1
    
#     for istep = 1:nsteps

#         t = t + dt

#         #Update kinematic parameters
#         update_kinem(surf, t)

#         #Update flow field parameters if any
#         update_externalvel(curfield, t)

#         #Update bound vortex positions
#         update_boundpos(surf, dt)

#         res = 1.

#         RHStransp = zeros(surf.ndiv*2-2)
        
#         iter_delu = zeros(surf.ndiv-1)
#         iter_Eu = zeros(surf.ndiv-1)
#         quc = zeros(surf.ndiv-1)
        
        
#         #Iterate for viscous solution and interaction

#         iter = 0

#         while (iter < 2 || res > 5e-3)
            
#             iter += 1
            
#             #println("res")
            
#             if iter > 1
#                 pop!(curfield.tev)
#                 pop!(curfield.lev)
#             end
            
#             nlev = length(curfield.lev)
#             ntev = length(curfield.tev)
            
#             println(i_xsu, "    ", surf.qu[i_xsu])

#             vx = surf.qu[i_xsu]/sqrt(1. + (surf.thick_slope[i_xsu] + surf.cam_slope[i_xsu])^2)
#             vz = surf.qu[i_xsu]*(surf.thick_slope[i_xsu] + surf.cam_slope[i_xsu])/sqrt(1. + (surf.thick_slope[i_xsu] + surf.cam_slope[i_xsu])^2)

#             alpha = surf.kinem.alpha
#             R = [cos(alpha) -sin(alpha); sin(alpha) cos(alpha)]
                        
#             vx1, vz1 = R*[vx; 0]
#             vx2, vz2 = R*[0; vz]
#             vx = vx1 + vx2
#             vz = vz1 + vz2
            
#             #if surf.levflag[1] == 0
#                 xloc_lev = surf.bnd_x_u[i_xsu] + 0.5*vx*dt
#                 zloc_lev = surf.bnd_z_u[i_xsu] + 0.5*vz*dt + 0.01
#             #else
#             #    xloc_lev = surf.bnd_x_u[i_xsu] + (1. /3.)*(curfield.lev[nlev].x - surf.bnd_x_u[i_xsu])
#  #               zloc_lev = surf.bnd_z_u[i_xsu] + (1. /3.)*(curfield.lev[nlev].z - surf.bnd_z_u[i_xsu])
#   #          end
#             levstr = 0.5*surf.qu[i_xsu]^2*dt
#             push!(curfield.lev, TwoDVort(xloc_lev, zloc_lev, levstr, vcore, 0., 0.))

#             #if i_xsu != 139
#             println(surf.levflag[1], "      ", xloc_lev, "    ", zloc_lev)
#             #    error("here")
#             #end

#             #Update induced velocities to include effect of last shed vortex
#             update_indbound(surf, curfield)

#             surf, xloc_tev, zloc_tev = update_thickLHS2V(surf, curfield, dt, vcore, i_xsl)
#             update_thickRHS(surf, curfield)
            
#             soln = surf.LHS[1:surf.ndiv*2-2, 1:surf.naterm*2+1] \ (surf.RHS[1:surf.ndiv*2-2] + RHStransp[:])

#             res = sqrt(sum((soln[1:end-1] .- [surf.aterm; surf.bterm]).^2))

#             #Assign the solution
#             for i = 1:surf.naterm
#                 surf.aterm[i] = soln[i]
#                 surf.bterm[i] = soln[i+surf.naterm]
#             end
#             tevstr = soln[2*surf.naterm+1]*surf.uref*surf.c
            
#             push!(curfield.tev, TwoDVort(xloc_tev, zloc_tev, tevstr, vcore, 0., 0.))
            
#             #Update induced velocities to include effect of last shed vortex
#             update_indbound(surf, curfield)
            
#             #Dont update qu during iteration - just the inviscid value.
            
#             surf.qu[:], surf.ql[:] = calc_edgeVel(surf, [curfield.u[1], curfield.w[1]])
#             #end
            
#             quxc = zeros(surf.ndiv-1)
#             qutc = zeros(surf.ndiv-1)
#             dsu = zeros(surf.ndiv-1)
#             suc = zeros(surf.ndiv-1)
            
#                 #Solve the FV problem at cell centres
#             for i = 1:surf.ndiv-1
#                 quc[i] = (surf.qu[i] + surf.qu[i+1])/2
#                 suc[i] = (surf.su[i] + surf.su[i+1])/2
#             end
    
#             quxc[2:end] = diff(quc)./diff(suc)
#             quxc[1] = 2*quxc[2] - quxc[3]
            
#             smoothEdges!(quxc, 5)
            
#             if istep == 1
#                 quprevc[:] = quc[:]
#             end

#             qutc[:] .= (quc[:] .- quprevc[:])./dt

#             dsu = diff(surf.su)
            
#             w0 = [surf.delu surf.delu.*(surf.Eu .+ 1)]

#             w0, i_xsu = FVMIBLgridvar(w0, quc, qutc, quxc, dsu, t-dt, t)
#             i_xsl = i_xsu

#             #println(i_xsu)
            
#             iter_delu[:] = w0[:,1]
#             iter_Eu[:] = w0[:,2]./w0[:,1] .- 1.

#             #println(iter_delu)
            
#             #plot(surf.x, iter_delu)
            
#             smoothEdges!(iter_delu, 5)
            
#             surf.wtu[2:end] = (1/sqrt(Re))*diff(surf.qu.*[iter_delu; iter_delu[end]])./diff(surf.x)
#             surf.wtu[1] = surf.wtu[2]

#             smoothEdges!(surf.wtu, 10)

#             surf.wtl[:] = -surf.wtu[:]
            
#             if istep ==70 && iter == 5
#                 figure(1)
#                 plot(suc, iter_delu)
#                 figure(2)
#                 plot(suc, quxc)
#                 figure(3)
#                 plot(surf.x, surf.wtu)
#                 error("here")
#             end

#             RHStransp[:] .= 0.

#             #Add transpiration velocity to RHS
#             for i = 2:surf.ndiv-1
#                 RHStransp[i-1] = 0.5*(surf.wtu[i] + surf.wtl[i])
#                 RHStransp[surf.ndiv+i-3] = 0.5*(surf.wtu[i] - surf.wtl[i])
#             end

#             println(istep, "   ", res)
#         end

#         if i_xsu != 139
#             surf.levflag[1] = 1
#         end
        
#         quprevc[:] = quc[:]
        
#         #Assign bl
#         surf.delu[:] = iter_delu[:]
#         #surf.delu[end] = surf.delu[end-1]
#         #surf.Eu[end] = surf.Eu[end-1]

#         surf.Eu[:] = iter_Eu[:]

#         #Calculate adot
#         update_atermdot(surf, dt)

#         #Set previous values of aterm to be used for derivatives in next time step
#         surf.a0prev[1] = surf.a0[1]
#         for ia = 1:3
#             surf.aprev[ia] = surf.aterm[ia]
#         end

#         #Calculate bound vortex strengths
#         update_bv_src(surf)

#         #Add effect of transpiration to sources and vortices

#         #Wake rollup
#         wakeroll(surf, curfield, dt)

#         #Force calculation
#         cnc, cnnc, cn, cs, cl, cd, int_wax, int_c, int_t = calc_forces(surf, int_wax, int_c, int_t, dt)
        
#         # write flow details if required
#         if writeflag == 1
#             if istep in writeArray
#                 dirname = "$(round(t,sigdigits=nround))"
#                 writeStamp(dirname, t, surf, curfield)
#             end
#         end
                
#         vle = surf.qu[1]

#         if vle > 0.
#             qspl = Spline1D(surf.x, surf.ql)
#             stag = try
#                 roots(qspl, maxn=1)[1]
# OA            catch
#                 0.
#             end
#         else
#             qspl = Spline1D(surf.x, surf.qu)
#             stag = try
#                 roots(qspl, maxn=1)[1]
#             catch
#                 0.
#             end
#         end

#         mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, vle,
#         cl, cd, cnc, cnnc, cn, cs, stag])

#         println("here")
#     end

#     mat = mat'



#     f = open("resultsSummary", "w")
#     Serialization.serialize(f, ["#time \t", "alpha (deg) \t", "h/c \t", "u/uref \t", "A0 \t", "Cl \t", "Cd \t", "Cm \n"])
#     writedlm(f, mat)
#     close(f)

#     return mat, surf, curfield

# end

# function IBLsheddisp(surf::TwoDSurfThickBL, curfield::TwoDFlowField, ncell::Int64, nsteps::Int64 = 300, dtstar::Float64 = 0.015, startflag = 0, writeflag = 0, writeInterval = 1000., delvort = delNone(); maxwrite = 50, nround=6)

#     # If a restart directory is provided, read in the simulation data
#     if startflag == 0
#         mat = zeros(0, 12)
#         t = 0.

#     elseif startflag == 1
#         dirvec = readdir()
#         dirresults = map(x->(v = tryparse(Float64,x); typeof(v) == Nothing ? 0.0 : v),dirvec)
#         latestTime = maximum(dirresults)
#         mat = DelimitedFiles.readdlm("resultsSummary")
#         t = mat[end,1]
#     else
#         throw("invalid start flag, should be 0 or 1")
#     end
#     mat = mat'

#     dt = dtstar*surf.c/surf.uref

#     Re = 1000000

#     thick_orig = zeros(surf.ndiv)
#     cam_orig = zeros(surf.ndiv)
#     thick_slope_orig = zeros(surf.ndiv)
#     cam_slope_orig = zeros(surf.ndiv)
    
#     thick_orig[:] = surf.thick[:]
#     cam_orig[:] = surf.cam[:]
#     thick_slope_orig[:] = surf.thick_slope[:]
#     cam_slope_orig[:] = surf.cam_slope[:]
    
#     # if writeflag is on, determine the timesteps to write at
#     if writeflag == 1
#         writeArray = Int64[]
#         tTot = nsteps*dt
#         for i = 1:maxwrite
#             tcur = writeInterval*real(i)
#             if t > tTot
#                 break
#             else
#                 push!(writeArray, Int(round(tcur/dt)))
#             end
#         end
#     end

#     vcore = 0.02*surf.c

#     int_wax = zeros(surf.ndiv)
#     int_c = zeros(surf.ndiv)
#     int_t = zeros(surf.ndiv)

#     quprevc = zeros(surf.ndiv-1)

#     i_xsu = surf.ndiv-1
#     i_xsl = surf.ndiv-1
    
#     for istep = 1:nsteps

#         t = t + dt

#         #Update kinematic parameters
#         update_kinem(surf, t)

#         #Update flow field parameters if any
#         update_externalvel(curfield, t)

#         #Update bound vortex positions
#         update_boundpos(surf, dt)

#         res = 1.

#         RHStransp = zeros(surf.ndiv*2-2)
        
#         iter_delu = zeros(surf.ndiv)
#         iter_Eu = zeros(surf.ndiv-1)
#         quc = zeros(surf.ndiv-1)
        
        
#         #Iterate for viscous solution and interaction

#         iter = 0

#         while (iter < 2 || res > 1e-2)
            
#             iter += 1
            
#             #println("res")
            
#             if iter > 1
#                 pop!(curfield.tev)
#                 pop!(curfield.lev)
#             end
            
#             nlev = length(curfield.lev)
#             ntev = length(curfield.tev)
            
#             #println(i_xsu, "    ", surf.qu[i_xsu])
            
#             vx = surf.qu[i_xsu]/sqrt(1. + (surf.thick_slope[i_xsu] + surf.cam_slope[i_xsu])^2)
#             vz = surf.qu[i_xsu]*(surf.thick_slope[i_xsu] + surf.cam_slope[i_xsu])/sqrt(1. + (surf.thick_slope[i_xsu] + surf.cam_slope[i_xsu])^2)

#             alpha = surf.kinem.alpha
#             R = [cos(alpha) -sin(alpha); sin(alpha) cos(alpha)]
                        
#             vx1, vz1 = R*[vx; 0]
#             vx2, vz2 = R*[0; vz]
#             vx = vx1 + vx2
#             vz = vz1 + vz2
            
#             #if surf.levflag[1] == 0
#             xloc_lev = surf.bnd_x_u[i_xsu] + 0.5*vx*dt
#             zloc_lev = surf.bnd_z_u[i_xsu] + 0.5*vz*dt + 0.001
#             #else
#             #    xloc_lev = surf.bnd_x_u[i_xsu] + (1. /3.)*(curfield.lev[nlev].x - surf.bnd_x_u[i_xsu])
#  #               zloc_lev = surf.bnd_z_u[i_xsu] + (1. /3.)*(curfield.lev[nlev].z - surf.bnd_z_u[i_xsu])
#   #          end
#             levstr = 0.15*surf.qu[i_xsu]^2*dt
#             push!(curfield.lev, TwoDVort(xloc_lev, zloc_lev, levstr, vcore, 0., 0.))

#             #if i_xsu != 139
#             #println(surf.levflag[1], "      ", xloc_lev, "    ", zloc_lev)
#             #    error("here")
#             #end

#             #Update induced velocities to include effect of last shed vortex
#             update_indbound(surf, curfield)

#             surf, xloc_tev, zloc_tev = update_thickLHS2V(surf, curfield, dt, vcore, i_xsl)
#             update_thickRHS(surf, curfield)
            
#             soln = surf.LHS[1:surf.ndiv*2-2, 1:surf.naterm*2+1] \ (surf.RHS[1:surf.ndiv*2-2] + RHStransp[:])

#             res = sqrt(sum((soln[1:end-1] .- [surf.aterm; surf.bterm]).^2))

#             #Assign the solution
#             for i = 1:surf.naterm
#                 surf.aterm[i] = soln[i]
#                 surf.bterm[i] = soln[i+surf.naterm]
#             end
#             tevstr = soln[2*surf.naterm+1]*surf.uref*surf.c
            
#             push!(curfield.tev, TwoDVort(xloc_tev, zloc_tev, tevstr, vcore, 0., 0.))
            
#             #Update induced velocities to include effect of last shed vortex
#             update_indbound(surf, curfield)
            
#             #Dont update qu during iteration - just the inviscid value.
            
#             surf.qu[:], surf.ql[:] = calc_edgeVel(surf, [curfield.u[1], curfield.w[1]])
#             #end
            
#             quxc = zeros(surf.ndiv-1)
#             qutc = zeros(surf.ndiv-1)
#             dsu = zeros(surf.ndiv-1)
#             suc = zeros(surf.ndiv-1)
            
#             #Solve the FV problem at cell centres
#             for i = 1:surf.ndiv-1
#                 quc[i] = (surf.qu[i] + surf.qu[i+1])/2
#                 suc[i] = (surf.su[i] + surf.su[i+1])/2
#             end
    
#             quxc[2:end] = diff(quc)./diff(suc)
#             quxc[1] = 2*quxc[2] - quxc[3]
            
#             smoothEdges!(quxc, 5)
            
#             if istep == 1
#                 quprevc[:] = quc[:]
#             end

#             qutc[:] .= (quc[:] .- quprevc[:])./dt

#             dsu = diff(surf.su)
            
#             w0 = [surf.delu surf.delu.*(surf.Eu .+ 1)]
            
#             w0, i_xsu = FVMIBLgridvar(w0, quc, qutc, quxc, dsu, t-dt, t)
#             i_xsl = i_xsu

#             #println(i_xsu)
            
#             iter_delu[:] = [w0[:,1];w0[end,1]]
#             iter_Eu[:] = w0[:,2]./w0[:,1] .- 1.

#             #println(iter_delu)
            
#             #plot(surf.x, iter_delu)
            
#             smoothEdges!(iter_delu, 5)

#             #Update airfoil shape
#             for i = 2:surf.ndiv
#                 th_u = 1/sqrt(Re)*iter_delu[i]/sqrt(1. + (cam_slope_orig[i] + thick_slope_orig[i])^2)
#                 th_l = -1/sqrt(Re)*iter_delu[i]/sqrt(1. + (cam_slope_orig[i] - thick_slope_orig[i])^2)
#                 surf.thick[i] = thick_orig[i] + (th_u - th_l)/2
#                 surf.cam[i] = cam_orig[i] + (th_u + th_l)/2
#             end
#             #Update slopes
#             thick_spl = Spline1D(surf.x, surf.thick)
#             surf.thick_slope[:] = derivative(thick_spl, surf.x)
#             cam_spl = Spline1D(surf.x, surf.cam)
#             surf.cam_slope[:] = derivative(cam_spl, surf.x)
            
#             if istep == 10 && iter == 1
#                 figure(1)
#                 plot(suc, iter_delu[1:end-1])
#                 figure(2)
#                 plot(suc, quxc)
#                 figure(3)
#                 plot(surf.x, surf.thick+surf.cam)
#                 plot(surf.x, thick_orig+cam_orig)
#                 error("here")
#             end

#             println(istep, "   ", res)
#         end

#         if i_xsu != 139
#             surf.levflag[1] = 1
#         end
        
#         quprevc[:] = quc[:]
        
#         #Assign bl
#         surf.delu[:] = iter_delu[1:end-1]
#         #surf.delu[end] = surf.delu[end-1]
#         #surf.Eu[end] = surf.Eu[end-1]

#         surf.Eu[:] = iter_Eu[:]

#         #Calculate adot
#         update_atermdot(surf, dt)

#         #Set previous values of aterm to be used for derivatives in next time step
#         surf.a0prev[1] = surf.a0[1]
#         for ia = 1:3
#             surf.aprev[ia] = surf.aterm[ia]
#         end

#         #Calculate bound vortex strengths
#         update_bv_src(surf)

#         #Add effect of transpiration to sources and vortices

#         #Wake rollup
#         wakeroll(surf, curfield, dt)

#         #Force calculation
#         cnc, cnnc, cn, cs, cl, cd, int_wax, int_c, int_t = calc_forces(surf, int_wax, int_c, int_t, dt)
        
#         # write flow details if required
#         if writeflag == 1
#             if istep in writeArray
#                 dirname = "$(round(t,sigdigits=nround))"
#                 writeStamp(dirname, t, surf, curfield)
#             end
#         end
                
#         vle = surf.qu[1]

#         if vle > 0.
#             qspl = Spline1D(surf.x, surf.ql)
#             stag = try
#                 roots(qspl, maxn=1)[1]
# OA            catch
#                 0.
#             end
#         else
#             qspl = Spline1D(surf.x, surf.qu)
#             stag = try
#                 roots(qspl, maxn=1)[1]
#             catch
#                 0.
#             end
#         end

#         mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, vle,
#         cl, cd, cnc, cnnc, cn, cs, stag])

#         println("here")
#     end

#     mat = mat'



#     f = open("resultsSummary", "w")
#     Serialization.serialize(f, ["#time \t", "alpha (deg) \t", "h/c \t", "u/uref \t", "A0 \t", "Cl \t", "Cd \t", "Cm \n"])
#     writedlm(f, mat)
#     close(f)

#     return mat, surf, curfield

# end


