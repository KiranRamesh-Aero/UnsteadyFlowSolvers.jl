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


function IBLshed(surf::TwoDSurfThickBL, curfield::TwoDFlowField, ncell::Int64, nsteps::Int64 = 300, dtstar::Float64 = 0.015, startflag = 0, writeflag = 0, writeInterval = 1000., delvort = delNone(); maxwrite = 50, nround=6)

    # If a restart directory is provided, read in the simulation data
    if startflag == 0
        mat = zeros(0, 12)
        t = 0.

    elseif startflag == 1
        dirvec = readdir()
        dirresults = map(x->(v = tryparse(Float64,x); typeof(v) == Nothing ? 0.0 : v),dirvec)
        latestTime = maximum(dirresults)
        mat = DelimitedFiles.readdlm("resultsSummary")
        t = mat[end,1]
    else
        throw("invalid start flag, should be 0 or 1")
    end
    mat = mat'

    dt = dtstar*surf.c/surf.uref

    Re = 1000000


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

    vcore = 0.02*surf.c

    int_wax = zeros(surf.ndiv)
    int_c = zeros(surf.ndiv)
    int_t = zeros(surf.ndiv)

    quprevc = zeros(surf.ndiv-1)

    i_xsu = surf.ndiv-1
    i_xsl = surf.ndiv-1
    
    for istep = 1:nsteps

        t = t + dt

        #Update kinematic parameters
        update_kinem(surf, t)

        #Update flow field parameters if any
        update_externalvel(curfield, t)

        #Update bound vortex positions
        update_boundpos(surf, dt)

        res = 1.

        RHStransp = zeros(surf.ndiv*2-2)
        
        iter_delu = zeros(surf.ndiv-1)
        iter_Eu = zeros(surf.ndiv-1)
        quc = zeros(surf.ndiv-1)
        
        
        #Iterate for viscous solution and interaction

        iter = 0

        while (iter < 2 || res > 5e-3)
            
            iter += 1
            
            #println("res")
            
            if iter > 1
                pop!(curfield.tev)
                pop!(curfield.lev)
            end
            
            nlev = length(curfield.lev)
            ntev = length(curfield.tev)
            
            println(i_xsu, "    ", surf.qu[i_xsu])

            vx = surf.qu[i_xsu]/sqrt(1. + (surf.thick_slope[i_xsu] + surf.cam_slope[i_xsu])^2)
            vz = surf.qu[i_xsu]*(surf.thick_slope[i_xsu] + surf.cam_slope[i_xsu])/sqrt(1. + (surf.thick_slope[i_xsu] + surf.cam_slope[i_xsu])^2)

            alpha = surf.kinem.alpha
            R = [cos(alpha) -sin(alpha); sin(alpha) cos(alpha)]
                        
            vx1, vz1 = R*[vx; 0]
            vx2, vz2 = R*[0; vz]
            vx = vx1 + vx2
            vz = vz1 + vz2
            
            #if surf.levflag[1] == 0
                xloc_lev = surf.bnd_x_u[i_xsu] + 0.5*vx*dt
                zloc_lev = surf.bnd_z_u[i_xsu] + 0.5*vz*dt + 0.01
            #else
            #    xloc_lev = surf.bnd_x_u[i_xsu] + (1. /3.)*(curfield.lev[nlev].x - surf.bnd_x_u[i_xsu])
 #               zloc_lev = surf.bnd_z_u[i_xsu] + (1. /3.)*(curfield.lev[nlev].z - surf.bnd_z_u[i_xsu])
  #          end
            levstr = 0.5*surf.qu[i_xsu]^2*dt
            push!(curfield.lev, TwoDVort(xloc_lev, zloc_lev, levstr, vcore, 0., 0.))

            #if i_xsu != 139
            println(surf.levflag[1], "      ", xloc_lev, "    ", zloc_lev)
            #    error("here")
            #end

            #Update induced velocities to include effect of last shed vortex
            update_indbound(surf, curfield)

            surf, xloc_tev, zloc_tev = update_thickLHS2V(surf, curfield, dt, vcore, i_xsl)
            update_thickRHS(surf, curfield)
            
            soln = surf.LHS[1:surf.ndiv*2-2, 1:surf.naterm*2+1] \ (surf.RHS[1:surf.ndiv*2-2] + RHStransp[:])

            res = sqrt(sum((soln[1:end-1] .- [surf.aterm; surf.bterm]).^2))

            #Assign the solution
            for i = 1:surf.naterm
                surf.aterm[i] = soln[i]
                surf.bterm[i] = soln[i+surf.naterm]
            end
            tevstr = soln[2*surf.naterm+1]*surf.uref*surf.c
            
            push!(curfield.tev, TwoDVort(xloc_tev, zloc_tev, tevstr, vcore, 0., 0.))
            
            #Update induced velocities to include effect of last shed vortex
            update_indbound(surf, curfield)
            
            #Dont update qu during iteration - just the inviscid value.
            
            surf.qu[:], surf.ql[:] = calc_edgeVel(surf, [curfield.u[1], curfield.w[1]])
            #end
            
            quxc = zeros(surf.ndiv-1)
            qutc = zeros(surf.ndiv-1)
            dsu = zeros(surf.ndiv-1)
            suc = zeros(surf.ndiv-1)
            
                #Solve the FV problem at cell centres
            for i = 1:surf.ndiv-1
                quc[i] = (surf.qu[i] + surf.qu[i+1])/2
                suc[i] = (surf.su[i] + surf.su[i+1])/2
            end
    
            quxc[2:end] = diff(quc)./diff(suc)
            quxc[1] = 2*quxc[2] - quxc[3]
            
            smoothEdges!(quxc, 5)
            
            if istep == 1
                quprevc[:] = quc[:]
            end

            qutc[:] .= (quc[:] .- quprevc[:])./dt

            dsu = diff(surf.su)
            
            w0 = [surf.delu surf.delu.*(surf.Eu .+ 1)]

            w0, i_xsu = FVMIBLgridvar(w0, quc, qutc, quxc, dsu, t-dt, t)
            i_xsl = i_xsu

            #println(i_xsu)
            
            iter_delu[:] = w0[:,1]
            iter_Eu[:] = w0[:,2]./w0[:,1] .- 1.

            #println(iter_delu)
            
            #plot(surf.x, iter_delu)
            
            smoothEdges!(iter_delu, 5)
            
            surf.wtu[2:end] = (1/sqrt(Re))*diff(surf.qu.*[iter_delu; iter_delu[end]])./diff(surf.x)
            surf.wtu[1] = surf.wtu[2]

            smoothEdges!(surf.wtu, 10)

            surf.wtl[:] = -surf.wtu[:]
            
            if istep ==70 && iter == 5
                figure(1)
                plot(suc, iter_delu)
                figure(2)
                plot(suc, quxc)
                figure(3)
                plot(surf.x, surf.wtu)
                error("here")
            end

            RHStransp[:] .= 0.

            #Add transpiration velocity to RHS
            for i = 2:surf.ndiv-1
                RHStransp[i-1] = 0.5*(surf.wtu[i] + surf.wtl[i])
                RHStransp[surf.ndiv+i-3] = 0.5*(surf.wtu[i] - surf.wtl[i])
            end

            println(istep, "   ", res)
        end

        if i_xsu != 139
            surf.levflag[1] = 1
        end
        
        quprevc[:] = quc[:]
        
        #Assign bl
        surf.delu[:] = iter_delu[:]
        #surf.delu[end] = surf.delu[end-1]
        #surf.Eu[end] = surf.Eu[end-1]

        surf.Eu[:] = iter_Eu[:]

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
        
        # write flow details if required
        if writeflag == 1
            if istep in writeArray
                dirname = "$(round(t,sigdigits=nround))"
                writeStamp(dirname, t, surf, curfield)
            end
        end
                
        vle = surf.qu[1]

        if vle > 0.
            qspl = Spline1D(surf.x, surf.ql)
            stag = try
                roots(qspl, maxn=1)[1]
OA            catch
                0.
            end
        else
            qspl = Spline1D(surf.x, surf.qu)
            stag = try
                roots(qspl, maxn=1)[1]
            catch
                0.
            end
        end

        mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, vle,
        cl, cd, cnc, cnnc, cn, cs, stag])

        println("here")
    end

    mat = mat'



    f = open("resultsSummary", "w")
    Serialization.serialize(f, ["#time \t", "alpha (deg) \t", "h/c \t", "u/uref \t", "A0 \t", "Cl \t", "Cd \t", "Cm \n"])
    writedlm(f, mat)
    close(f)

    return mat, surf, curfield

end

