function QSLLT_lautat(surf :: ThreeDSurfSimple, field :: ThreeDFieldSimple, nsteps :: Int64 = 500, dtstar :: Float64 = 0.015, startflag = 0, writeflag = 0, writeInterval = 1000., delvort = delNone(); maxwrite = 50, nround = 6)

    # If a restart directory is provided, read in the simulation data
    if startflag == 0
        mat = zeros(0, 4+surf.nspan*8)
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

    # if writeflag is on, determine the timesteps to write at
    if writeflag == 1
        writeArray = Int64[]
        tTot = nsteps*dtstar
        for i = 1:maxwrite
            tcur = writeInterval*real(i)
            if t > tTot
                break
            else
                push!(writeArray, Int(round(tcur/dtstar)))
            end
        end
    end
    
    
    dt = dtstar*surf.cref/surf.uref

    cl = zeros(surf.nspan)
    cd = zeros(surf.nspan)
    cm = zeros(surf.nspan)

    #Initialise flowfield
    for is = 1:surf.nspan
        push!(field.f2d, TwoDFlowField())
    end
    
    for istep = 1:nsteps
        #Udpate current time
        t = t + dt
        
        for is = 1:surf.nspan
            #Update external flowfield
            update_externalvel(field.f2d[is], t)
            
            #Update kinematic parameters
            update_kinem(surf.s2d[is], t)

            #Update bound vortex positions
            update_boundpos(surf.s2d[is], dt)

            #Add a TEV with dummy strength
            place_tev(surf.s2d[is], field.f2d[is], dt)

            update_indbound(surf.s2d[is], field.f2d[is])
        end

        kelv = KelvinConditionQSLLT(surf, field)

        #Solve for TEV strength to satisfy Kelvin condition
        
        soln = nlsolve(not_in_place(kelv), 0.0001*ones(surf.nspan))
        
        for is = 1:surf.nspan
            field.f2d[is].tev[end].s = soln.zero[is]
        end
        
        calc_a2toan3d(surf)
        
        for is = 1:surf.nspan
            #Update 3D effect on A0 and A0dot
            surf.s2d[is].a0[1] += surf.a03d[is]
                       
            #Update rest of Fourier terms
            update_a2toan(surf.s2d[is])

            #Update 3D effect on An
            for ia = 1:surf.naterm
                surf.s2d[is].aterm[ia] += surf.aterm3d[ia,is]
            end
            
            #Update derivatives
            update_adot(surf.s2d[is],dt)
            
            #Set previous values of aterm to be used for derivatives in next time step
            surf.s2d[is].a0prev[1] = surf.s2d[is].a0[1]
            for ia = 1:3
                surf.s2d[is].aprev[ia] = surf.s2d[is].aterm[ia]
            end
            
            #Calculate bound vortex strengths
            update_bv(surf.s2d[is])
        end

        #Wake rollup
        wakeroll(surf, field, dt)

        cl3d, cd3d, cm3d, cl, cd, cm = calc_forces(surf, field, dt)

        # write flow details if required
        if writeflag == 1
            if istep in writeArray
                dirname = "$(round(t, digits=nround))"
                writeStamp(dirname, t, surf, field)
            end
        end
        
        matvect = [t, cl3d, cd3d, cm3d]
        for is = 1:surf.nspan
            matvect = [matvect; [surf.s2d[is].kinem.alpha, surf.s2d[is].kinem.h,
                                 surf.s2d[is].kinem.u, surf.s2d[is].a0[1], cl[is]*surf.s2d[is].c/surf.cref, cd[is]*surf.s2d[is].c/surf.cref, cm[is]*surf.s2d[is].c*surf.s2d[is].c/(surf.cref*surf.cref), surf.a03d[is]]]
        end
        mat = hcat(mat,matvect)
    end
    mat = mat'
    
    f = open("resultsSummary", "w")
    Serialization.serialize(f, ["#time \t", "CL3D \t", "CD3D \t", "CM3D \t", "alpha -1 (deg) \t", "h/c -1 \t", "u/uref -1 \t", "A0 -1 \t", "Cl -1 \t", "Cd -2\t", "Cm -2 \t", "A03D-1 \t", "alpha -2 ...\n"])
    DelimitedFiles.writedlm(f, mat)
    close(f)

    mat, surf, field, cl, cd, cm

end


function QSLLT_ldvm(surf :: ThreeDSurfSimple, field :: ThreeDFieldSimple, nsteps :: Int64 = 500, dtstar :: Float64 = 0.015, startflag = 0, writeflag = 0, writeInterval = 1000., delvort = delNone(); maxwrite = 50, nround = 6)

    # If a restart directory is provided, read in the simulation data
    if startflag == 0
        mat = zeros(0, 4+surf.nspan*8)
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

    # if writeflag is on, determine the timesteps to write at
    if writeflag == 1
        writeArray = Int64[]
        tTot = nsteps*dtstar
        for i = 1:maxwrite
            tcur = writeInterval*real(i)
            if t > tTot
                break
            else
                push!(writeArray, Int(round(tcur/dtstar)))
            end
        end
    end
    
    
    dt = dtstar*surf.cref/surf.uref

    cl = zeros(surf.nspan)
    cd = zeros(surf.nspan)
    cm = zeros(surf.nspan)
    lesp = zeros(surf.nspan)
    
    #Initialise flowfield
    for is = 1:surf.nspan
        push!(field.f2d, TwoDFlowField())
    end
    
    for istep = 1:nsteps
        #Udpate current time
        t = t + dt
        
        for is = 1:surf.nspan
            #Update external flowfield
            update_externalvel(field.f2d[is], t)
            
            #Update kinematic parameters
            update_kinem(surf.s2d[is], t)

            #Update bound vortex positions
            update_boundpos(surf.s2d[is], dt)

            #Add a TEV with dummy strength
            place_tev(surf.s2d[is], field.f2d[is], dt)

            update_indbound(surf.s2d[is], field.f2d[is])
        end

        kelv = KelvinConditionQSLLT(surf, field)

        #Solve for TEV strength to satisfy Kelvin condition
        
        soln = nlsolve(not_in_place(kelv), 0.0001*ones(surf.nspan))
        
        for is = 1:surf.nspan
            field.f2d[is].tev[end].s = soln.zero[is]
        end
        
        for is = 1:surf.nspan
            lesp[is] = surf.s2d[is].a0[1]
        end

        shed_ind = Int[]

        #Check for LESP condition
        for is = 1:surf.nspan
            if abs(lesp[is]) > surf.s2d[is].lespcrit[1]
                push!(shed_ind, is)
            end
        end
        nshed = length(shed_ind)

        #println(map(q->q.a0[1], surf.s2d))
        #println(shed_ind)
        #error("h")
        
        if nshed > 0
            #2D iteration at section if LESP_crit is exceeded 
            #Add LEVs with dummy strength
            for is = 1:surf.nspan
                if is in shed_ind
                    place_lev(surf.s2d[is], field.f2d[is], dt)
                end
            end
            
            kelvkutta = KelvinKuttaQSLLT(surf, field, shed_ind)
            soln = nlsolve(not_in_place(kelvkutta), [-0.01*ones(surf.nspan); 0.01*ones(nshed)])

            levcount = 0
            for is = 1:surf.nspan
                if is in shed_ind
                    field.f2d[is].tev[end].s = soln.zero[is]
                    levcount +=1
                    field.f2d[is].lev[end].s = soln.zero[levcount+surf.nspan]
                    surf.s2d[is].levflag[1] = 1
                else
                    field.f2d[is].tev[end].s = soln.zero[is]
                    surf.s2d[is].levflag[1] = 0
                end
            end
        else
            for is = 1:surf.nspan
                surf.s2d[is].levflag[1] = 0
            end
        end

        calc_a2toan3d(surf)
            
        for is = 1:surf.nspan
            #Update 3D effect on A0 and A0dot
            surf.s2d[is].a0[1] += surf.a03d[is]
            
            #Update rest of Fourier terms
            update_a2toan(surf.s2d[is])
            
            #Update 3D effect on An
            for ia = 1:surf.naterm
                surf.s2d[is].aterm[ia] += surf.aterm3d[ia,is]
            end

            #Update derivatives
            update_adot(surf.s2d[is],dt)
                
            #Set previous values of aterm to be used for derivatives in next time step
           
            surf.s2d[is].a0prev[1] = surf.s2d[is].a0[1]
            for ia = 1:3
                surf.s2d[is].aprev[ia] = surf.s2d[is].aterm[ia]
            end
            
            #Calculate bound vortex strengths
            update_bv(surf.s2d[is])
            
            # Delete or merge vortices if required
            controlVortCount(delvort, surf.s2d[is].bnd_x[Int(round(surf.ndiv/2))], surf.s2d[is].bnd_z[Int(round(surf.s2d[is].ndiv/2))], field.f2d[is])
            
            #Wake rollup
            wakeroll(surf.s2d[is], field.f2d[is], dt)
        end
        
        cl3d, cd3d, cm3d, cl, cd, cm = calc_forces(surf, field, dt)
        
        # write flow details if required
        if writeflag == 1
            if istep in writeArray
                dirname = "$(round(t, digits=nround))"
                writeStamp(dirname, t, surf, field)
            end
        end
        
        matvect = [t, cl3d, cd3d, cm3d]
        for is = 1:surf.nspan
            matvect = [matvect; [surf.s2d[is].kinem.alpha, surf.s2d[is].kinem.h,
                                 surf.s2d[is].kinem.u, surf.s2d[is].a0[1], cl[is]*surf.s2d[is].c/surf.cref, cd[is]*surf.s2d[is].c/surf.cref, cm[is]*surf.s2d[is].c*surf.s2d[is].c/(surf.cref*surf.cref), surf.a03d[is]]]
        end
        mat = hcat(mat,matvect)
    end
    mat = mat'
    
    f = open("resultsSummary", "w")
    Serialization.serialize(f, ["#time \t", "CL3D \t", "CD3D \t", "CM3D \t", "alpha -1 (deg) \t", "h/c -1 \t", "u/uref -1 \t", "A0 -1 \t", "Cl -1 \t", "Cd -2\t", "Cm -2 \t", "A03D-1 \t", "alpha -2 ...\n"])
    DelimitedFiles.writedlm(f, mat)
    close(f)

    mat, surf, field, cl, cd, cm

end

# function QSLLT_ldvm(surf :: ThreeDSurfSimple, field :: ThreeDFieldSimple, nsteps :: Int64 = 500, dtstar :: Float64 = 0.015, startflag = 0, writeflag = 0, writeInterval = 1000., delvort = delNone(); maxwrite = 50, nround = 6)

#     # If a restart directory is provided, read in the simulation data
#     if startflag == 0
#         mat = zeros(0, 4+surf.nspan*8)
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

#     # if writeflag is on, determine the timesteps to write at
#     if writeflag == 1
#         writeArray = Int64[]
#         tTot = nsteps*dtstar
#         for i = 1:maxwrite
#             tcur = writeInterval*real(i)
#             if t > tTot
#                 break
#             else
#                 push!(writeArray, Int(round(tcur/dtstar)))
#             end
#         end
#     end
    
    
#     dt = dtstar*surf.cref/surf.uref

#     cl = zeros(surf.nspan)
#     cd = zeros(surf.nspan)
#     cm = zeros(surf.nspan)
#     lesp = zeros(surf.nspan)
#     levstr = zeros(surf.nspan)
    
#     #Initialise flowfield
#     for is = 1:surf.nspan
#         push!(field.f2d, TwoDFlowField())
#     end
    
#     for istep = 1:nsteps
#         #Udpate current time
#         t = t + dt
        
#         for is = 1:surf.nspan
#             #Update external flowfield
#             update_externalvel(field.f2d[is], t)
            
#             #Update kinematic parameters
#             update_kinem(surf.s2d[is], t)

#             #Update bound vortex positions
#             update_boundpos(surf.s2d[is], dt)

#             #Add a TEV with dummy strength
#             place_tev(surf.s2d[is], field.f2d[is], dt)

#             update_indbound(surf.s2d[is], field.f2d[is])
#         end

#         kelv = KelvinConditionQSLLT(surf, field)

#         #Solve for TEV strength to satisfy Kelvin condition
        
#         soln = nlsolve(not_in_place(kelv), 0.0001*ones(surf.nspan))
        
#         for is = 1:surf.nspan
#             field.f2d[is].tev[end].s = soln.zero[is]
#         end
        
#         calc_a2toan3d(surf)
        
#         for is = 1:surf.nspan
#             #Update 3D effect on A0 and A0dot
#             surf.s2d[is].a0[1] += surf.a03d[is]
                       
#             #Update rest of Fourier terms
#             update_a2toan(surf.s2d[is])

#             #Update 3D effect on An
#             for ia = 1:surf.naterm
#                 surf.s2d[is].aterm[ia] += surf.aterm3d[ia,is]
#             end
            
#             #Update derivatives
#             update_adot(surf.s2d[is],dt)

#             lesp[is] = surf.s2d[is].a0[1]
#         end

#         shed_ind = Int[]

#         #Check for LESP condition
#         for is = 1:surf.nspan
#             if abs(lesp[is]) > surf.s2d[is].lespcrit[1]
#                 push!(shed_ind, is)
#             end
#         end
#         nshed = length(shed_ind)

#         #println(map(q->q.a0[1], surf.s2d))
#         #println(shed_ind)
#         #error("h")
        
#         if nshed > 0
#             #2D iteration at section if LESP_crit is exceeded 
#             #Add LEVs with dummy strength
#             for is = 1:surf.nspan
#                 if is in shed_ind
#                     place_lev(surf.s2d[is], field.f2d[is], dt)
#                 else
#                     push!(field.f2d[is].lev, TwoDVort(0., 0., 0., 0., 0., 0.))
#                 end
#             end

#             #Start global iteration
#             iter = 0
#             tev_previter = zeros(surf.nspan)
#             while true
#                 iter += 1

                
#                 for is = 1:surf.nspan
#                     if is in shed_ind
#                         levstr[is] -= field.f2d[is].lev[end].s  
#                         kelvkutta = KelvinKuttaQSLLT_sep(surf, field, is)
#                         soln = nlsolve(not_in_place(kelvkutta), [-0.01; 0.01])
#                         field.f2d[is].tev[end].s = soln.zero[1]
#                         field.f2d[is].lev[end].s = soln.zero[2]
#                         levstr[is] += field.f2d[is].lev[end].s  
#                     else
#                         kelv = KelvinConditionQSLLT_sep(surf, field, is)
#                         soln = nlsolve(not_in_place(kelv), [-0.01;])
#                         field.f2d[is].tev[end].s = soln.zero[1]
#                         levstr[is] = 0.
#                     end
#                 end

#                 #Update cumulative LEV strength
#                 calc_a0a13d_wlev2(surf, levstr)
                
#                 if sum(abs.(map(q->q.tev[end].s, field.f2d) - tev_previter)) < 1e-4
#                     break
#                 end
#                 tev_previter = map(q->q.tev[end].s, field.f2d)
#             end
#             println(iter)
            
#             calc_a2toan3d(surf)
            
#             for is = 1:surf.nspan
#                 #Update 3D effect on A0 and A0dot
#                 surf.s2d[is].a0[1] += surf.a03d[is]
                       
#                 #Update rest of Fourier terms
#                 update_a2toan(surf.s2d[is])
                
#                 #Update 3D effect on An
#                     for ia = 1:surf.naterm
#                     surf.s2d[is].aterm[ia] += surf.aterm3d[ia,is]
#                 end
#             end
            
#         else
#             for is = 1:surf.nspan
#                 surf.s2d[is].levflag[1] = 0
#             end
#         end
        
#         #Set previous values of aterm to be used for derivatives in next time step
#         for is = 1:surf.nspan
#             surf.s2d[is].a0prev[1] = surf.s2d[is].a0[1]
#             for ia = 1:3
#                 surf.s2d[is].aprev[ia] = surf.s2d[is].aterm[ia]
#             end
            
#             #Calculate bound vortex strengths
#             update_bv(surf.s2d[is])
            
#             # Delete or merge vortices if required
#             controlVortCount(delvort, surf.s2d[is].bnd_x[Int(round(surf.ndiv/2))], surf.s2d[is].bnd_z[Int(round(surf.s2d[is].ndiv/2))], field.f2d[is])
            
#             #Wake rollup
#             wakeroll(surf.s2d[is], field.f2d[is], dt)
#         end
        
#         cl3d, cd3d, cm3d, cl, cd, cm = calc_forces(surf)
        
#         # write flow details if required
#         if writeflag == 1
#             if istep in writeArray
#                 dirname = "$(round(t, digits=nround))"
#                 writeStamp(dirname, t, surf, field)
#             end
#         end
        
#         matvect = [t, cl3d, cd3d, cm3d]
#         for is = 1:surf.nspan
#             matvect = [matvect; [surf.s2d[is].kinem.alpha, surf.s2d[is].kinem.h,
#                                  surf.s2d[is].kinem.u, surf.s2d[is].a0[1], cl[is]*surf.s2d[is].c/surf.cref, cd[is]*surf.s2d[is].c/surf.cref, cm[is]*surf.s2d[is].c*surf.s2d[is].c/(surf.cref*surf.cref), surf.a03d[is]]]
#         end
#         mat = hcat(mat,matvect)
#     end
#     mat = mat'
    
#     f = open("resultsSummary", "w")
#     Serialization.serialize(f, ["#time \t", "CL3D \t", "CD3D \t", "CM3D \t", "alpha -1 (deg) \t", "h/c -1 \t", "u/uref -1 \t", "A0 -1 \t", "Cl -1 \t", "Cd -2\t", "Cm -2 \t", "A03D-1 \t", "alpha -2 ...\n"])
#     DelimitedFiles.writedlm(f, mat)
#     close(f)

#     mat, surf, field, cl, cd, cm

# end
