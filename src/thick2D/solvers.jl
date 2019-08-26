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
