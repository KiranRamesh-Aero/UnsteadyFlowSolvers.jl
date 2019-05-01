# function transpCoupled(surf::TwoDSurfThickBL, curfield::TwoDFlowField, ncell::Int64, nsteps::Int64 = 300, dtstar::Float64 = 0.015, startflag = 0, writeflag = 0, writeInterval = 1000., delvort = delNone(); maxwrite = 50, nround=6)

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

#         #Solve inviscid problem
#         invsoln = surf.LHS[1:surf.ndiv*2-2, 1:surf.naterm*2+1] \ surf.RHS[1:surf.ndiv*2-2]



#         #Assign the solution
#         for i = 1:surf.naterm
#             surf.aterm[i] = invsoln[i]
#             surf.bterm[i] = invsoln[i+surf.naterm]
#         end
#         tevstr = invsoln[2*surf.naterm+1]*surf.uref*surf.c
#         push!(curfield.tev, TwoDVort(xloc_tev, zloc_tev, tevstr, vcore, 0., 0.))

#         avisc = zeros(surf.naterm)
#         bvisc = zeros(surf.naterm)

#         #Update induced velocities to include effect of last shed vortex
#         update_indbound(surf, curfield)

#         res = 1.

#         qux = zeros(surf.ndiv)
#         qlx = zeros(surf.ndiv)
#         qut = zeros(surf.ndiv)
#         qlt = zeros(surf.ndiv)

#         iter_delu = zeros(surf.ndiv)
#         iter_dell = zeros(surf.ndiv)
#         iter_Eu = zeros(surf.ndiv)
#         iter_El = zeros(surf.ndiv)
#         wtu = zeros(surf.ndiv)
#         wtl = zeros(surf.ndiv)

#         #Iterate for viscous solution and interaction

#         iter = 0
#         while iter < 1

#             iter += 1

#             surf.aterm[:] = invsoln[1:surf.naterm] .+ avisc[:]
#             surf.bterm[:] = invsoln[surf.naterm+1:2*surf.naterm] .+ bvisc[:]

#             surf.qu[:], surf.ql[:]= calc_edgeVel(surf, [curfield.u[1], curfield.w[1]])

#             qux[2:end] = diff(surf.qu)./diff(surf.theta)
#             qux[1] = 2*qux[2] - qux[3]
#             qlx[2:end] = diff(surf.ql)./diff(surf.theta)
#             qlx[1] = 2*qlx[2] - qlx[3]

#             if iter == 3
#                 error("now")
#             end

#             if istep == 1
#                 surf.quprev[:] = surf.qu[:]
#                 surf.qlprev[:] = surf.ql[:]
#             end

#             qut[:] = (surf.qu[:] .- surf.quprev[:])./dt
#             qlt[:] = (surf.ql[:] .- surf.qlprev[:])./dt

#             w0 = [surf.delu surf.delu.*(surf.Eu .+ 1)]
#             w, j1 ,j2 = FVMIBLorig(w0, surf.qu, qut, qux, surf.theta, dt)

#             iter_delu = w[:,1]
#             iter_Eu = (w[:,2]./w[:,1]) .- 1.0
#             iter_dell[:] = surf.delu[:]
#             iter_El[:] = surf.Eu[:]

#             plot(surf.x, surf.delu)

#             wtu[2:end] = (1/100)*diff(surf.qu.*iter_delu)./diff(surf.x)
#             wtu[1] = wtu[2]
#             wtl[:] = -wtu[:]

#             figure()
#             plot(surf.x, wtu)
#             error("here")


#             #wtu[wtu .> 0.1] .= 0.1
#             #wtu[wtu .< -0.1] .= -0.1

#             RHStransp = zeros(surf.ndiv*2-2)

#             #Add transpiration velocity to RHS
#             for i = 2:surf.ndiv-1
#                 RHStransp[i-1] = 0.5*(sqrt(1 + (surf.cam_slope[i] + surf.thick_slope[i])^2)*wtu[i]
#                                        + sqrt(1 + (surf.cam_slope[i] - surf.thick_slope[i])^2)*wtl[i])

#                 RHStransp[surf.ndiv+i-3] = 0.5*(sqrt(1 + (surf.cam_slope[i] + surf.thick_slope[i])^2)*wtu[i]
#                                                       - sqrt(1 + (surf.cam_slope[i] - surf.thick_slope[i])^2)*wtl[i])
#             end

#             plot(surf.x[2:end-1], RHStransp[surf.ndiv-1:2*surf.ndiv-4])
#             println(RHStransp)

#             #Solve viscous problem
#             viscsoln = surf.LHS[1:surf.ndiv*2-2, 1:surf.naterm*2] \ RHStransp[1:surf.ndiv*2-2]

#             res = sqrt(sum((viscsoln[:] .- [avisc;bvisc]).^2))

#             println(res)

#             avisc[:] = viscsoln[1:surf.naterm]
#             bvisc[:] = viscsoln[surf.naterm+1:end]

#             #            println(bvisc)

#         end

#         #Assign bl
#         surf.delu = iter_delu[:]
#         surf.dell = iter_dell[:]
#         surf.Eu[:] = iter_Eu[:]
#         surf.El[:] = iter_El[:]

#         #Calculate adot
#         update_atermdot(surf, dt)

#         #Set previous values of aterm to be used for derivatives in next time step
#         surf.a0prev[1] = surf.a0[1]
#         for ia = 1:3
#             surf.aprev[ia] = surf.aterm[ia]
#         end

#         surf.quprev[:] = surf.qu[:]
#         surf.qlprev[:] = surf.ql[:]

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
#     DelimitedFiles.writedlm(f, mat)
#     close(f)

#     return mat, surf, curfield

# end




function transpTogether(surf::TwoDSurfThickBL, curfield::TwoDFlowField, ncell::Int64, nsteps::Int64 = 300, dtstar::Float64 = 0.015, startflag = 0, writeflag = 0, writeInterval = 1000., delvort = delNone(); maxwrite = 50, nround=6)

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

    for istep = 1:nsteps

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

        res = 1.

        RHStransp = zeros(surf.ndiv*2-2)
        
        iter_delu = zeros(surf.ndiv-1)
        iter_Eu = zeros(surf.ndiv-1)
        wtu = zeros(surf.ndiv)
        wtl = zeros(surf.ndiv)

        quc = zeros(surf.ndiv-1)
        
        
        #Iterate for viscous solution and interaction

        iter = 0

        while (iter < 2 || res > 5e-3)
            
            iter += 1

            println("res")

            if iter > 1
                pop!(curfield.tev)
            end

            #Update induced velocities to include effect of last shed vortex
            update_indbound(surf, curfield)

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
            #This basically means no strong coupling.
            #if iter == 1
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

            w0 = FVMIBLgridvar(w0, quc, qutc, quxc, dsu, t-dt, t)
            
            iter_delu[:] = w0[:,1]
            iter_Eu[:] = w0[:,2]./w0[:,1] .- 1.

            #println(iter_delu)
            
            #plot(surf.x, iter_delu)
            
            smoothEdges!(iter_delu, 5)
            
            wtu[2:end] = (1/sqrt(Re))*diff(surf.qu.*[iter_delu; iter_delu[end]])./diff(surf.x)
            wtu[1] = wtu[2]

            # smoothEdges!(wtu, 10)

             wtl[:] = -wtu[:]
            
            if istep ==50 && iter == 2
                figure(1)
                plot(suc, iter_delu)
                figure(2)
                plot(suc, quxc)
                figure(3)
                plot(surf.x, wtu)
                error("here")
            end

            RHStransp[:] .= 0.

            #Add transpiration velocity to RHS
            for i = 2:surf.ndiv-1
                RHStransp[i-1] = 0.5*(wtu[i] + wtl[i])
                RHStransp[surf.ndiv+i-3] = 0.5*(wtu[i] - wtl[i])
            end

            println(istep, "   ", res)
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
            
            if nlev == 0
                xloc_lev = surf.bnd_x_u[i_xsu] + 0.5*vx*dt
                zloc_lev = surf.bnd_z_u[i_xsu] + 0.5*vz*dt
            else
                xloc_lev = surf.bnd_x_u[i_xsu] + (1. /3.)*(curfield.lev[nlev].x - surf.bnd_x_u[i_xsu])
                zloc_lev = surf.bnd_z_u[i_xsu] + (1. /3.)*(curfield.lev[nlev].z - surf.bnd_z_u[i_xsu])
            end
            levstr = 0.5*surf.qu[i_xsu]^2*dt
            push!(curfield.lev, TwoDVort(xloc_lev, zloc_lev, levstr, vcore, 0., 0.))
            
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
            
            if istep ==40 && iter == 2
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



# function transpTogetherWake(surf::TwoDSurfThickBL, curfield::TwoDFlowField, ncell::Int64, nsteps::Int64 = 300, dtstar::Float64 = 0.015, startflag = 0, writeflag = 0, writeInterval = 1000., delvort = delNone(); maxwrite = 50, nround=6)

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
        
#         iter_delu = zeros(2*surf.ndiv-2)
#         iter_Eu = zeros(2*surf.ndiv-2)
#         wtu = zeros(surf.ndiv)
#         wtl = zeros(surf.ndiv)

#         quc = zeros(2*surf.ndiv-2)
        
        
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

#             #Form wake approximation
#             sw = zeros(surf.ndiv)
#             quw = zeros(surf.ndiv)
#             delw = zeros(surf.ndiv)
            
#             for i = 1:surf.ndiv
#                 sw[i] = surf.su[end] + 2*surf.c*(i-1)/(surf.ndiv-1)
#             end
#             for i = 1:surf.ndiv
#                 a = log10(surf.delu[end])
#                 delw[i] = 10^(a - 3.2*(sw[i] - surf.su[end]))
#                 if surf.qu[end] < 1
#                     a = log10(1. - surf.qu[end])
#                     quw[i] = -10^(a - 3.2*(sw[i] - surf.su[end])) + 1
#                 else
#                     a = log10(surf.qu[end] - 1.)
#                     quw[i] = 10^(a - 3.2*(sw[i] - surf.su[end])) + 1
#                 end
#             end
            
#             quxc = zeros(2*surf.ndiv-2)
#             qutc = zeros(2*surf.ndiv-2)
#             dsu = zeros(2*surf.ndiv-2)
#             suc = zeros(2*surf.ndiv-2)
            
#             #Solve the FV problem at cell centres
#             for i = 1:surf.ndiv-1
#                 quc[i] = (surf.qu[i] + surf.qu[i+1])/2
#                 quc[i+surf.ndiv-1] = (quw[i] + quw[i+1])/2
#                 suc[i] = (surf.su[i] + surf.su[i+1])/2
#                 suc[i+surf.ndiv-1] = (sw[i] + sw[i+1])/2
#             end
    
#             quxc[2:end] = diff(quc)./diff(suc)
#             quxc[1] = 2*quxc[2] - quxc[3]
            
#             #smoothEdges!(quxc, 5)
            
#             if istep == 1
#                 quprevc[:] = quc[:]
#             end

#             qutc[:] .= (quc[:] .- quprevc[:])./dt

#             dsu = [diff(surf.su); diff(sw)]
            
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
            
#             if istep ==50 && iter == 2
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
