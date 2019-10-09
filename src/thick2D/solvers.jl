# function lautat(surf::TwoDSurfThick, curfield::TwoDFlowField, nsteps::Int64 = 500,
#                 dtstar::Float64 = 0.015, startflag = 0, writeflag = 0, writeInterval = 1000.,
#                 delvort = delNone(); maxwrite = 100, nround=6)

#     # If a restart directory is provided, read in the simulation data
#     if startflag == 0
#         mat = zeros(0, 12)
#         t = 0.
#     elseif startflag == 1
#         dirvec = readdir()
#         dirresults = map(x->(v = tryparse(Float64,x); isnull(v) ? 0.0 : get(v)),dirvec)
#         latestTime = maximum(dirresults)
#         mat = readdlm("resultsSummary")
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

#     phi_u = zeros(surf.ndiv)
#     phi_l = zeros(surf.ndiv)

#     uind_up = zeros(surf.ndiv)
#     wind_up = zeros(surf.ndiv)
#     uind_lp = zeros(surf.ndiv)
#     wind_lp = zeros(surf.ndiv)

#     bound_circ = 0.
    
#     for istep = 1:nsteps

#         #Udpate current time
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
#         update_thickRHS(surf, curfield, bound_circ[1])
#         #surf.RHS[2*surf.ndiv-3] = bound_circ
#         #Now solve the matrix problem
#         #soln = surf.LHS[[1:surf.ndiv*2-3;2*surf.ndiv-1], 1:surf.naterm*2+2] \ surf.RHS[[1:surf.ndiv*2-3; 2*surf.ndiv-1]]
#         soln = surf.LHS[1:surf.ndiv*2-3, 1:surf.naterm*2+1] \ surf.RHS[1:surf.ndiv*2-3]
        
#         #Assign the solution
#         for i = 1:surf.naterm
#             surf.aterm[i] = soln[i]
#             surf.bterm[i] = soln[i+surf.naterm]
#         end
#         tevstr = soln[2*surf.naterm+1]
#         push!(curfield.tev, TwoDVort(xloc_tev, zloc_tev, tevstr, vcore, 0., 0.))
        
#         #Calculate adot
#         update_atermdot(surf, dt)

#         #Set previous values of aterm to be used for derivatives in next time step
#         surf.a0prev[1] = surf.a0[1]
#         for ia = 1:3
#             surf.aprev[ia] = surf.aterm[ia]
#         end

#         #Update the bound circulation value
#         bound_circ = surf.LHS[2*surf.ndiv-3]*[surf.aterm;surf.bterm;0.]
        
#         #Update induced velocities to include effect of last shed vortex
#         update_indbound(surf, curfield)

#         uind_up[:] = surf.uind_u[:]
#         wind_up[:] = surf.wind_u[:]
#         uind_lp[:] = surf.uind_l[:]
#         wind_lp[:] = surf.wind_l[:]
        
#         #Calculate bound vortex strengths
#         update_bv_src(surf)

#         #Wake rollup
#         wakeroll(surf, curfield, dt)

#         #Force calculation
#         cnc, cnnc, cn, cs, cl, cd, int_wax, int_c, int_t = calc_forces(surf, int_wax, int_c, int_t, dt)
#         #println(phi_u)
#         qu, ql, phi_u, phi_l, cpu, cpl = calc_edgeVel_cp(surf, [curfield.u[1]; curfield.w[1]], phi_u, phi_l, dt)
#         #println(phi_l)
#         # write flow details if required
#         if writeflag == 1
#             if istep in writeArray
#                 dirname = "$(round(t,sigdigits=nround))"
#                 writeStamp(dirname, t, surf, curfield, qu, ql, cpu, cpl)
#             end
#         end

#         #plot(surf.x, qu)
        
        
#         #LE velocity and stagnation point location
#         #vle = (surf.kinem.u + curfield.u[1])*sin(surf.kinem.alpha) + (curfield.w[1] - surf.kinem.hdot)*cos(surf.kinem.alpha) - surf.kinem.alphadot*surf.pvt*surf.c + sum(surf.aterm) + surf.wind_u[1]  
        

#         vle = qu[1]

#         # if vle > 0.
#         #     qspl = Spline1D(surf.x, ql)
#         #     stag = try
#         #         roots(qspl, maxn=1)[1]
#         #     catch
#         #         0.
#         #     end
#         # else
#         #     qspl = Spline1D(surf.x, qu)
#         #     stag = try
#         #         roots(qspl, maxn=1)[1]
#         #     catch
#         #         0.
#         #     end
#         # end

#          if vle > 0.
#              stag = surf.x[argmin(ql)]
#          else
#              stag = surf.x[argmin(qu)]
#          end
#              # else
#         #     qspl = Spline1D(surf.x, qu)
#         #     stag = try
#         #         roots(qspl, maxn=1)[1]
#         #     catch
#         #         0.
#         #     end
#         # end


        
#         mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, vle,
#                         cl, cd, cnc, cnnc, cn, bound_circ[1], stag])

#     end
    
#     mat = mat'
    
#     f = open("resultsSummary", "w")
#     Serialization.serialize(f, ["#time \t", "alpha (rad) \t", "h/c \t", "u/uref \t", "A0 \t", "Cl \t", "Cd \t", "Cm \n"])
#     writedlm(f, mat)
#     close(f)

#     mat, surf, curfield
# end

function lautat(surf::TwoDSurfThick, curfield::TwoDFlowField, nsteps::Int64 = 500,
                dtstar::Float64 = 0.015, startflag = 0, writeflag = 0, writeInterval = 1000.,
                delvort = delNone(); maxwrite = 100, nround=6)

    # If a restart directory is provided, read in the simulation data
    if startflag == 0
        mat = zeros(0, 12)
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

    vcore = 0.02*surf.c

    int_wax = zeros(surf.ndiv)
    int_c = zeros(surf.ndiv)
    int_t = zeros(surf.ndiv)

    phi_u = zeros(surf.ndiv)
    phi_l = zeros(surf.ndiv)

    uind_up = zeros(surf.ndiv)
    wind_up = zeros(surf.ndiv)
    uind_lp = zeros(surf.ndiv)
    wind_lp = zeros(surf.ndiv)

    bound_circ = 0.

    tevstr = zeros(100)
    restev = zeros(100)
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
            xloc_tev = surf.bnd_x_chord[surf.ndiv] + 0.5*surf.kinem.u*dt
            zloc_tev = surf.bnd_z_chord[surf.ndiv]
        else
            xloc_tev = surf.bnd_x_chord[surf.ndiv] + (1. /3.)*(curfield.tev[ntev].x - surf.bnd_x_chord[surf.ndiv])
            zloc_tev = surf.bnd_z_chord[surf.ndiv] + (1. /3.)*(curfield.tev[ntev].z - surf.bnd_z_chord[surf.ndiv])
        end
        #Iteratively solve for TE strength
        
        #No need to update LHS
        tevstr[1] = -0.1
        
        push!(curfield.tev, TwoDVort(xloc_tev, zloc_tev, tevstr[1], vcore, 0., 0.))
        update_indbound(surf, curfield)
        #Update RHS vector
        update_thickRHS(surf, curfield)
        soln = surf.LHS[1:surf.ndiv*2-4, 1:surf.naterm*2] \ surf.RHS[1:surf.ndiv*2-4]
        
        #Assign the solution
        for i = 1:surf.naterm
            surf.aterm[i] = soln[i]
            surf.bterm[i] = soln[i+surf.naterm]
        end
        
        qu, ql, phi_u_temp, phi_l_temp, cpu, cpl = calc_edgeVel_cp(surf, [curfield.u[1]; curfield.w[1]], phi_u, phi_l, dt)
        restev[1] = (phi_u_temp[surf.ndiv] - phi_l_temp[surf.ndiv]) - (phi_u[surf.ndiv] - phi_l[surf.ndiv]) + tevstr[1]
        
        pop!(curfield.tev)
        tevstr[2] = 0.1
        #-((phi_u_temp[surf.ndiv] - phi_l_temp[surf.ndiv]) - (phi_u[surf.ndiv] - phi_l[surf.ndiv]))
        push!(curfield.tev, TwoDVort(xloc_tev, zloc_tev, tevstr[2], vcore, 0., 0.))
        update_indbound(surf, curfield)
        #Update RHS vector
        update_thickRHS(surf, curfield)
        soln = surf.LHS[1:surf.ndiv*2-4, 1:surf.naterm*2] \ surf.RHS[1:surf.ndiv*2-4]
        #Assign the solution
        for i = 1:surf.naterm
            surf.aterm[i] = soln[i]
            surf.bterm[i] = soln[i+surf.naterm]
        end
        
        qu, ql, phi_u_temp, phi_l_temp, cpu, cpl = calc_edgeVel_cp(surf, [curfield.u[1]; curfield.w[1]], phi_u, phi_l, dt)
        restev[2] = (phi_u_temp[surf.ndiv] - phi_l_temp[surf.ndiv]) - (phi_u[surf.ndiv] - phi_l[surf.ndiv]) + tevstr[2]
        
        iter = 2
        
        while restev[iter] > 1e-6
            iter += 1
            pop!(curfield.tev)
            tevstr[iter] = tevstr[iter-1] - restev[iter-1]*(tevstr[iter-1] - tevstr[iter-2])/(restev[iter-1] - restev[iter-2])
            push!(curfield.tev, TwoDVort(xloc_tev, zloc_tev, tevstr[iter], vcore, 0., 0.))

            update_indbound(surf, curfield)
            #Update RHS vector
            update_thickRHS(surf, curfield)
            soln = surf.LHS[1:surf.ndiv*2-4, 1:surf.naterm*2] \ surf.RHS[1:surf.ndiv*2-4]
            
            #Assign the solution
            for i = 1:surf.naterm
                surf.aterm[i] = soln[i]
                surf.bterm[i] = soln[i+surf.naterm]
            end
            
            qu, ql, phi_u_temp, phi_l_temp, cpu, cpl = calc_edgeVel_cp(surf, [curfield.u[1]; curfield.w[1]], phi_u, phi_l, dt)

            restev[iter] = (phi_u_temp[surf.ndiv] - phi_l_temp[surf.ndiv]) - (phi_u[surf.ndiv] - phi_l[surf.ndiv]) + tevstr[iter]
            
            #println(istep, " " , iter, " ", tevstr[iter], " " , restev[iter], " ",  phi_u[end], " ", phi_l[end], " ", phi_u_temp[end], phi_l_temp[end])

        end
                
        
        #println(istep, " " , iter, " ", tevstr[iter], " ", restev[iter])
        #error("here")
        
        #Calculate adot
        update_atermdot(surf, dt)

        #Set previous values of aterm to be used for derivatives in next time step
        surf.a0prev[1] = surf.a0[1]
        for ia = 1:3
            surf.aprev[ia] = surf.aterm[ia]
        end

        #Update the bound circulation value
        bound_circ = surf.LHS[2*surf.ndiv-3]*[surf.aterm;surf.bterm;0.]
        
        #Update induced velocities to include effect of last shed vortex
        update_indbound(surf, curfield)

        uind_up[:] = surf.uind_u[:]
        wind_up[:] = surf.wind_u[:]
        uind_lp[:] = surf.uind_l[:]
        wind_lp[:] = surf.wind_l[:]
        
        #Calculate bound vortex strengths
        update_bv_src(surf)

        #Wake rollup
        wakeroll(surf, curfield, dt)

        #Force calculation
        cnc, cnnc, cn, cs, cl, cd, int_wax, int_c, int_t = calc_forces(surf, int_wax, int_c, int_t, dt)
        #println(phi_u)
        qu, ql, phi_u, phi_l, cpu, cpl = calc_edgeVel_cp(surf, [curfield.u[1]; curfield.w[1]], phi_u, phi_l, dt)

        #Calcualte forces from pressure
        ds = sqrt((surf.x[2] - surf.x[1])^2 + (surf.cam[2] + surf.thick[2] - surf.cam[1] - surf.thick[1])^2)
        vert = 1. /sqrt(1. + (surf.cam_slope[2] + surf.thick_slope[2])^2)
        hor = -(surf.cam_slope[2] + surf.thick_slope[2])/sqrt(1. + (surf.cam_slope[2] + surf.thick_slope[2])^2)
        cs = 0.5*(cpu[2]*hor + cpu[1])*ds
        cn = -0.5*(cpu[2]*vert)*ds
        ds = sqrt((surf.x[2] - surf.x[1])^2 + (surf.cam[2] - surf.thick[2] - surf.cam[1] + surf.thick[1])^2)
        vert = -1. /sqrt(1. + (surf.cam_slope[2] - surf.thick_slope[2])^2)
        hor = (surf.cam_slope[2] - surf.thick_slope[2])/sqrt(1. + (surf.cam_slope[2] - surf.thick_slope[2])^2)
        cs += 0.5*(cpl[2]*hor + cpl[1])*ds
        cn -= 0.5*(cpl[2]*vert)*ds

        for i = 3:surf.ndiv
            ds = sqrt((surf.x[i] - surf.x[i-1])^2 + (surf.cam[i] + surf.thick[i] - surf.cam[i-1] - surf.thick[i-1])^2)
            vert = 1. /sqrt(1. + (surf.cam_slope[i] + surf.thick_slope[i])^2)
            vert_p = 1. /sqrt(1. + (surf.cam_slope[i-1] + surf.thick_slope[i-1])^2)
            hor = -(surf.cam_slope[i] + surf.thick_slope[i])/sqrt(1. + (surf.cam_slope[i] + surf.thick_slope[i])^2)
            hor_p = -(surf.cam_slope[i-1] + surf.thick_slope[i-1])/sqrt(1. + (surf.cam_slope[i-1] + surf.thick_slope[i-1])^2)
            
            cs += 0.5*(cpu[i]*hor + cpu[i-1]*hor_p)*ds
            cn -= 0.5*(cpu[i]*vert + cpu[i-1]*vert_p)*ds

            ds = sqrt((surf.x[i] - surf.x[i-1])^2 + (surf.cam[i] - surf.thick[i] - surf.cam[i-1] + surf.thick[i-1])^2)
            vert = -1. /sqrt(1. + (surf.cam_slope[i] - surf.thick_slope[i])^2)
            vert_p = -1. /sqrt(1. + (surf.cam_slope[i-1] - surf.thick_slope[i-1])^2)
            hor = (surf.cam_slope[i] - surf.thick_slope[i])/sqrt(1. + (surf.cam_slope[i] + surf.thick_slope[i])^2)
            hor_p = (surf.cam_slope[i-1] - surf.thick_slope[i-1])/sqrt(1. + (surf.cam_slope[i-1] - surf.thick_slope[i-1])^2)
            
            cs += 0.5*(cpl[i]*hor + cpl[i-1]*hor_p)*ds
            cn -= 0.5*(cpl[i]*vert + cpl[i-1]*vert_p)*ds
        end
        cl = cn*cos(surf.kinem.alpha) + cs*sin(surf.kinem.alpha)
        cd = cn*sin(surf.kinem.alpha) - cs*cos(surf.kinem.alpha)
        
        #println(phi_l)
        # write flow details if required
        if writeflag == 1
            if istep in writeArray
                dirname = "$(round(t,sigdigits=nround))"
                writeStamp(dirname, t, surf, curfield, qu, ql, cpu, cpl)
            end
        end

        #plot(surf.x, qu)
        
        
        #LE velocity and stagnation point location
        #vle = (surf.kinem.u + curfield.u[1])*sin(surf.kinem.alpha) + (curfield.w[1] - surf.kinem.hdot)*cos(surf.kinem.alpha) - surf.kinem.alphadot*surf.pvt*surf.c + sum(surf.aterm) + surf.wind_u[1]  
        

        vle = qu[1]

        # if vle > 0.
        #     qspl = Spline1D(surf.x, ql)
        #     stag = try
        #         roots(qspl, maxn=1)[1]
        #     catch
        #         0.
        #     end
        # else
        #     qspl = Spline1D(surf.x, qu)
        #     stag = try
        #         roots(qspl, maxn=1)[1]
        #     catch
        #         0.
        #     end
        # end

         if vle > 0.
             stag = surf.x[argmin(ql)]
         else
             stag = surf.x[argmin(qu)]
         end
             # else
        #     qspl = Spline1D(surf.x, qu)
        #     stag = try
        #         roots(qspl, maxn=1)[1]
        #     catch
        #         0.
        #     end
        # end


        
        mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, vle,
                        cl, cd, cnc, cs, cn, phi_u[end]-phi_l[end], stag])

    end
    
    mat = mat'
    
    f = open("resultsSummary", "w")
    Serialization.serialize(f, ["#time \t", "alpha (rad) \t", "h/c \t", "u/uref \t", "A0 \t", "Cl \t", "Cd \t", "Cm \n"])
    writedlm(f, mat)
    close(f)

    mat, surf, curfield
end


function wkg_lautat(surf::TwoDSurfThick, curfield::TwoDFlowField, nsteps::Int64 = 500,
                dtstar::Float64 = 0.015, startflag = 0, writeflag = 0, writeInterval = 1000.,
                delvort = delNone(); maxwrite = 100, nround=6)

    # If a restart directory is provided, read in the simulation data
    if startflag == 0
        mat = zeros(0, 12)
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

    vcore = 0.02*surf.c

    int_wax = zeros(surf.ndiv)
    int_c = zeros(surf.ndiv)
    int_t = zeros(surf.ndiv)

    phi_u = zeros(surf.ndiv)
    phi_l = zeros(surf.ndiv)
    
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
        surf, xloc_tev, zloc_tev = wkg_update_thickLHS(surf, curfield, dt, vcore)

        #Construct RHS vector
        wkg_update_thickRHS(surf, curfield)
        
        #Now solve the matrix problem
        #soln = surf.LHS[[1:surf.ndiv*2-3;2*surf.ndiv-1], 1:surf.naterm*2+2] \ surf.RHS[[1:surf.ndiv*2-3; 2*surf.ndiv-1]]
        soln = surf.LHS[1:surf.ndiv*2-3, 1:surf.naterm*2+1] \ surf.RHS[1:surf.ndiv*2-3]
        
        #Assign the solution
        for i = 1:surf.naterm
            surf.aterm[i] = soln[i]
            surf.bterm[i] = soln[i+surf.naterm]
        end
        tevstr = soln[2*surf.naterm+1]*surf.uref*surf.c
        push!(curfield.tev, TwoDVort(xloc_tev, zloc_tev, tevstr, vcore, 0., 0.))
        
        #Calculate adot
        update_atermdot(surf, dt)

        #Set previous values of aterm to be used for derivatives in next time step
        surf.a0prev[1] = surf.a0[1]
        for ia = 1:3
            surf.aprev[ia] = surf.aterm[ia]
        end

        #Update induced velocities to include effect of last shed vortex
        update_indbound(surf, curfield)

        #Calculate bound vortex strengths
        update_bv_src(surf)

        #Wake rollup
        wakeroll(surf, curfield, dt)

        #Force calculation
        cnc, cnnc, cn, cs, cl, cd, int_wax, int_c, int_t = calc_forces(surf, int_wax, int_c, int_t, dt)
        #println(phi_u)
        qu, ql, phi_u, phi_l, cpu, cpl = calc_edgeVel_cp(surf, [curfield.u[1]; curfield.w[1]], phi_u, phi_l, dt)
        #println(phi_l)
        # write flow details if required
        if writeflag == 1
            if istep in writeArray
                dirname = "$(round(t,sigdigits=nround))"
                writeStamp(dirname, t, surf, curfield, qu, ql, cpu, cpl)
            end
        end

        #LE velocity and stagnation point location
        #vle = (surf.kinem.u + curfield.u[1])*sin(surf.kinem.alpha) + (curfield.w[1] - surf.kinem.hdot)*cos(surf.kinem.alpha) - surf.kinem.alphadot*surf.pvt*surf.c + sum(surf.aterm) + surf.wind_u[1]  
        

        vle = qu[1]

        if vle > 0.
            qspl = Spline1D(surf.x, ql)
            stag = try
                roots(qspl, maxn=1)[1]
            catch
                0.
            end
        else
            qspl = Spline1D(surf.x, qu)
            stag = try
                roots(qspl, maxn=1)[1]
            catch
                0.
            end
        end

        mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, vle,
                        cl, cd, cnc, cnnc, cn, cs, stag])

    end
    
    mat = mat'
    
    f = open("resultsSummary", "w")
    Serialization.serialize(f, ["#time \t", "alpha (rad) \t", "h/c \t", "u/uref \t", "A0 \t", "Cl \t", "Cd \t", "Cm \n"])
    writedlm(f, mat)
    close(f)

    mat, surf, curfield
end





# function lautat_kutta(surf::TwoDSurfThick, curfield::TwoDFlowField, nsteps::Int64 = 500,
#                 dtstar::Float64 = 0.015, startflag = 0, writeflag = 0, writeInterval = 1000.,
#                 delvort = delNone(); maxwrite = 100, nround=6)

#     # If a restart directory is provided, read in the simulation data
#     if startflag == 0
#         mat = zeros(0, 12)
#         t = 0.
#     elseif startflag == 1
#         dirvec = readdir()
#         dirresults = map(x->(v = tryparse(Float64,x); isnull(v) ? 0.0 : get(v)),dirvec)
#         latestTime = maximum(dirresults)
#         mat = readdlm("resultsSummary")
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

#     phi_u = zeros(surf.ndiv)
#     phi_l = zeros(surf.ndiv)

#     RHSvort = zeros(2*surf.ndiv-3)

#     qu = zeros(surf.ndiv)
#     ql = zeros(surf.ndiv)

#     tol_iter = 1e-6
    
#     for istep = 1:nsteps

#         #Udpate current time
#         t = t + dt

#         #Update kinematic parameters
#         update_kinem(surf, t)

#         #Update flow field parameters if any
#         update_externalvel(curfield, t)

#         #Update bound vortex positions
#         update_boundpos(surf, dt)

#         #Update incduced velocities on airfoil
#         update_indbound(surf, curfield)

#         res = 1
#         nlev = length(curfield.lev)
#         a_iter = zeros(surf.naterm)
#         b_iter = zeros(surf.naterm)
#         a_iter[:] = surf.aterm[:]
#         b_iter[:] = surf.bterm[:]

#         #Iterate for solution
#         function soln_iter()
#             #Add vortex on upper surface
#             if nlev == 0
#                 xloc_lev = surf.bnd_x_u[surf.ndiv] + 0.5*surf.kinem.u*dt
#                 zloc_lev =  surf.bnd_z_u[surf.ndiv]
#             else
#                 xloc_lev = surf.bnd_x_u[surf.ndiv] + (1. /3.)*(curfield.lev[nlev].x - surf.bnd_x_u[surf.ndiv])
#                 zloc_lev = surf.bnd_z_u[surf.ndiv]+(1. /3.)*(curfield.lev[nlev].z - surf.bnd_z_u[surf.ndiv])
#             end
            
#             qu, ql = calc_edgeVel_newFourier(surf, [curfield.u[1]; curfield.w[1]], a_iter, b_iter)
#             levstr = 0.5*qu[end]^2*dt
#             u_vort = TwoDVort(xloc_lev, zloc_lev, levstr, vcore, 0., 0.)
            
#             ind_new_u_u, ind_new_w_u = ind_vel([u_vort], surf.bnd_x_u, surf.bnd_z_u)
#             ind_new_u_l, ind_new_w_l = ind_vel([u_vort], surf.bnd_x_l, surf.bnd_z_l)
            
#             surf.uind_u[:] += ind_new_u_u[:]
#             surf.wind_u[:] += ind_new_w_u[:]
#             surf.uind_l[:] += ind_new_u_l[:]
#             surf.wind_l[:] += ind_new_w_l[:]
            
#             #Set up the matrix problem
#             surf, xloc_tev, zloc_tev = update_thickLHS_kutta(surf, curfield, dt, vcore)
            
#             #Construct RHS vector
#             update_thickRHS(surf, curfield)
            
#         RHSvort[2*surf.ndiv-3] = -100*levstr/(surf.uref*surf.c)

#         surf.uind_u[:] -= ind_new_u_u[:]
#         surf.wind_u[:] -= ind_new_w_u[:]
#         surf.uind_l[:] -= ind_new_u_l[:]
#         surf.wind_l[:] -= ind_new_w_l[:]
   
        
#         #Now solve the matrix problem
#         #soln = surf.LHS[[1:surf.ndiv*2-3;2*surf.ndiv-1], 1:surf.naterm*2+2] \ surf.RHS[[1:surf.ndiv*2-3; 2*surf.ndiv-1]]
#         soln = surf.LHS[1:surf.ndiv*2-3, 1:surf.naterm*2+1] \ (surf.RHS[1:surf.ndiv*2-3] + RHSvort[:])
        
#         #Assign the solution
#         for i = 1:surf.naterm
#             surf.aterm[i] = soln[i]
#             surf.bterm[i] = soln[i+surf.naterm]
#         end
#         tevstr = soln[2*surf.naterm+1]*surf.uref*surf.c
#         push!(curfield.tev, TwoDVort(xloc_tev, zloc_tev, tevstr, vcore, 0., 0.))
#         push!(curfield.lev, u_vort)
        
#         #Calculate adot
#         update_atermdot(surf, dt)

#         #Set previous values of aterm to be used for derivatives in next time step
#         surf.a0prev[1] = surf.a0[1]
#         for ia = 1:3
#             surf.aprev[ia] = surf.aterm[ia]
#         end

#         #Update induced velocities to include effect of last shed vortex
#         update_indbound(surf, curfield)

#         #Calculate bound vortex strengths
#         update_bv_src(surf)

#         #Wake rollup
#         wakeroll(surf, curfield, dt)

#         #Force calculation
#         cnc, cnnc, cn, cs, cl, cd, int_wax, int_c, int_t = calc_forces(surf, int_wax, int_c, int_t, dt)
#         #println(phi_u)
#         qu, ql, phi_u, phi_l, cpu, cpl = calc_edgeVel_cp(surf, [curfield.u[1]; curfield.w[1]], phi_u, phi_l, dt)
#         #println(phi_l)
#         # write flow details if required
#         if writeflag == 1
#             if istep in writeArray
#                 dirname = "$(round(t,sigdigits=nround))"
#                 writeStamp(dirname, t, surf, curfield, qu, ql, cpu, cpl)
#             end
#         end

#         #LE velocity and stagnation point location
#         #vle = (surf.kinem.u + curfield.u[1])*sin(surf.kinem.alpha) + (curfield.w[1] - surf.kinem.hdot)*cos(surf.kinem.alpha) - surf.kinem.alphadot*surf.pvt*surf.c + sum(surf.aterm) + surf.wind_u[1]  
        

#         vle = qu[1]

#         if vle > 0.
#             qspl = Spline1D(surf.x, ql)
#             stag = try
#                 roots(qspl, maxn=1)[1]
#             catch
#                 0.
#             end
#         else
#             qspl = Spline1D(surf.x, qu)
#             stag = try
#                 roots(qspl, maxn=1)[1]
#             catch
#                 0.
#             end
#         end

#         mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, vle,
#                         cl, cd, cnc, cnnc, cn, cs, stag])

#     end
    
#     mat = mat'
    
#     f = open("resultsSummary", "w")
#     Serialization.serialize(f, ["#time \t", "alpha (rad) \t", "h/c \t", "u/uref \t", "A0 \t", "Cl \t", "Cd \t", "Cm \n"])
#     writedlm(f, mat)
#     close(f)

#     mat, surf, curfield
# end

# function ldvm(surf::TwoDSurfThick, curfield::TwoDFlowField, nsteps::Int64 = 500,
#                 dtstar::Float64 = 0.015, startflag = 0, writeflag = 0, writeInterval = 1000.,
#                 delvort = delNone(); maxwrite = 100, nround=6)

#     # If a restart directory is provided, read in the simulation data
#     if startflag == 0
#         mat = zeros(0, 11)
#         t = 0.
#     elseif startflag == 1
#         dirvec = readdir()
#         dirresults = map(x->(v = tryparse(Float64,x); isnull(v) ? 0.0 : get(v)),dirvec)
#         latestTime = maximum(dirresults)
#         mat = readdlm("resultsSummary")
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

#         #Udpate current time
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
        
#         #Now solve the matrix problem
#         soln = surf.LHS[1:surf.ndiv*2-3, 1:surf.naterm*2+2] \ surf.RHS[1:surf.ndiv*2-3]
        
#         #Assign the solution
#         surf.a0[1] = soln[1]
#         for i = 1:surf.naterm
#             surf.aterm[i] = soln[i+1]
#             surf.bterm[i] = soln[i+surf.naterm+1]
#         end

#         #Calculate adot
#         surf.a0dot[1] = (surf.a0[1] - surf.a0prev[1])/dt
#         for ia = 1:surf.naterm
#             surf.adot[ia] = (surf.aterm[ia]-surf.aprev[ia])/dt
#         end

#         #Check if LEV shedding is true
#         lesp = sqrt(2. /surf.rho)*surf.a0[1]

#         if abs(lesp) > surf.lespcrit[1]

#             qu, ql = UNSflow.calc_edgeVel(surf, [0.; 0.])
#             qshed = maximum(qu)
#             levstr = dt*surf.uref^2*qshed^2/15

#             if surf.levflag[1] == 0
#                 le_vel_x = sqrt(2. /surf.rho)*surf.uref*surf.a0[1]*sin(surf.kinem.alpha)
#                 le_vel_z = sqrt(2. /surf.rho)*surf.uref*surf.a0[1]*cos(surf.kinem.alpha)
#                 xloc_lev = surf.bnd_x_u[1] + 0.5*le_vel_x*dt
#                 zloc_lev = surf.bnd_z_u[1] + 0.5*le_vel_z*dt
#             else
#                 xloc_lev = surf.bnd_x_u[1]+(1. /3.)*(curfield.lev[end].x - surf.bnd_x_u[1])
#                 zloc_lev = surf.bnd_z_u[1]+(1. /3.)*(curfield.lev[end].z - surf.bnd_z_u[1])
#             end

#             push!(curfield.lev, TwoDVort(xloc_lev, zloc_lev, levstr, vcore, 0., 0.))

#             #Update incduced velocities on airfoil
#             update_indbound(surf, curfield)
            
#             #Set up the matrix problem
#             surf, xloc_tev, zloc_tev = update_thickLHS(surf, curfield, dt, vcore)
            
#             #Construct RHS vector
#             update_thickRHS(surf, curfield)
            
#             #Now solve the matrix problem
#             soln = surf.LHS[1:surf.ndiv*2-3, 1:surf.naterm*2+2] \ surf.RHS[1:surf.ndiv*2-3]

#             #Assign the solution
#             surf.a0[1] = soln[1]
#             for i = 1:surf.naterm
#                 surf.aterm[i] = soln[i+1]
#                 surf.bterm[i] = soln[i+surf.naterm+1]
#             end
            
#             tevstr = soln[2*surf.naterm+2]*surf.uref*surf.c
#             push!(curfield.tev, TwoDVort(xloc_tev, zloc_tev, tevstr, vcore, 0., 0.))
               
#             surf.levflag[1] = 1
#         else
#             tevstr = soln[2*surf.naterm+2]*surf.uref*surf.c
#             push!(curfield.tev, TwoDVort(xloc_tev, zloc_tev, tevstr, vcore, 0., 0.))
#             surf.levflag[1] = 0
#         end
        
#         #Set previous values of aterm to be used for derivatives in next time step
#         surf.a0prev[1] = surf.a0[1]
#         for ia = 1:surf.naterm
#             surf.aprev[ia] = surf.aterm[ia]
#         end

#         #Update induced velocities to include effect of last shed vortex
#         update_indbound(surf, curfield)

#         #Calculate bound vortex strengths
#         update_bv_src(surf)

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

#         mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, surf.a0[1],
#                         cl, cd, cnc, cnnc, cn, cs])

#     end

#     mat = mat'

#     f = open("resultsSummary", "w")
#     Serialization.serialize(f, ["#time \t", "alpha (rad) \t", "h/c \t", "u/uref \t", "A0 \t", "Cl \t", "Cd \t", "Cm \n"])
#     DelimitedFiles.writedlm(f, mat)
#     close(f)

#     mat, surf, curfield
# end
