function lautat(surf::TwoDSurf, curfield::TwoDFlowField, nsteps::Int64 = 500; dtstar::Float64 = 0.015, delvort = DelVortDef(0, 0, 0.), mat = Array(Float64, 0, 8), kelv_enf = 0., writefile = "Nil")

    if (size(mat,1) > 0)
        t = mat[end,1]
    else
        t = 0.
    end

    mat = mat'

    dt = dtstar*surf.c/surf.uref

    #If required, write an archive file
    if writefile != "Nil"
        jldopen(writefile, "w") do file
            write(file, "dt",  dt)
            write(file, "dtstar",  dt)
            write(file, "Insurf",  surf)
            write(file, "Infield", curfield)
            write(file, "delvort", delvort)
            write(file, "nsteps", nsteps)
            cl, cd, cm, gamma, cn, cs, cnc, cncc, nonl, cm_n, cm_pvt, nonl_m = calc_forces_more(surf)
            g = g_create(file, "Init")
            g = write_stamp(surf, curfield, t, kelv_enf, g)
        end
    end
    
    
    #Intialise flowfield
    for istep = 1:nsteps
        #Udpate current time
        t = t + dt

        #Update kinematic parameters
        update_kinem(surf, t)

        #Update bound vortex positions
        update_boundpos(surf, dt)

        #Add a TEV with dummy strength
        place_tev(surf,curfield,dt)

        kelv = KelvinCondition(surf,curfield)
        #Solve for TEV strength to satisfy Kelvin condition
        #curfield.tev[length(curfield.tev)].s = secant_method(kelv, 0., -0.01)
        soln = nlsolve(not_in_place(kelv), [-0.01])
        curfield.tev[length(curfield.tev)].s = soln.zero[1]

        #Update adot
        update_a2a3adot(surf,dt)

        #Check for LEV and shed if yes
        #Set previous values of aterm to be used for derivatives in next time step
        surf.a0prev[1] = surf.a0[1]
        for ia = 1:3
            surf.aprev[ia] = surf.aterm[ia]
        end

        #Update rest of Fourier terms
        #update_a2toan(surf)

        #Calculate bound vortex strengths
        #update_bv(surf)

        #Remove vortices that are far away from airfoil
        if (delvort.flag == 1)
            if length(curfield.tev) > delvort.limit
                if (sqrt((curfield.tev[1].x- surf.bnd_x[div(surf.ndiv,2)])^2 + (curfield.tev[1].z- surf.bnd_z[div(surf.ndiv,2)])^2) > delvort.dist*surf.c)
                    kelv_enf = kelv_enf + curfield.tev[1].s
                    for i = 1:length(curfield.tev)-1
                        curfield.tev[i] = curfield.tev[i+1]
                    end
                    pop!(curfield.tev)
                end
            end
            if length(curfield.lev) > delvort.limit
                if (sqrt((curfield.lev[1].x- surf.bnd_x[div(surf.ndiv,2)])^2 + (curfield.lev[1].z- surf.bnd_z[div(surf.ndiv,2)])^2) > delvort.dist*surf.c)
                    kelv_enf = kelv_enf + curfield.lev[1].s
                    for i = 1:length(curfield.lev)-1
                        curfield.lev[i] = curfield.lev[i+1]
                    end
                    pop!(curfield.lev)
                end
            end
        end


        #wakeroll(surf, curfield)
        if writefile == "Nil"
            cl, cd, cm = calc_forces(surf)
        else
            cl, cd, cm, gamma, cn, cs, cnc, cncc, nonl, cm_n, cm_pvt, nonl_m = calc_forces_more(surf)
        end
        
        mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, surf.a0[1], cl, cd, cm])
        
        #If required, write an archive file
        if writefile != "Nil"
            jldopen(writefile, "r+") do file
                g = g_create(file, "t$istep")
                g = write_stamp(surf, curfield, t, kelv_enf, g)
            end
        end
        
    end
    mat = mat'

    mat, surf, curfield, kelv_enf

end

function lautat(surf::Vector{TwoDSurf}, curfield::TwoDFlowFieldMultSurf, nsteps::Int64 = 500; dtstar::Float64 = 0.015)

    
    #Deleting vortices is not currently supported
    #Resuming simulations is not currently supproted
    #Writefile is not currently supported
    #Move these to arguments when these are started

    nsurf = length(surf)
    mat = Array(Float64, nsteps, 8, nsurf)
    
    t = 0.
    
    dt = 100
    for i = 1:length(surf)
        dt = min(dt, dtstar*surf[i].c/surf[i].uref)
    end
    #Calculate tstep based on frequency as well.
    
    cl = zeros(nsurf)
    cd = zeros(nsurf)
    cm = zeros(nsurf)
    
    #Intialise flowfield
    for istep = 1:nsteps
        #Udpate current time
        t = t + dt
        
        #Update kinematic parameters
        update_kinem(surf, t)
        
        #Update flow field parameters if any
        update_externalvel(curfield, t)
        
        #Update bound vortex positions
        update_boundpos(surf, dt)
        
        #Add a TEV with dummy strength
        place_tev(surf,curfield,dt)
        
        kelv = KelvinConditionMultSurf(surf, curfield)
        
        #Solve for TEV strength to satisfy Kelvin condition
        soln = nlsolve(not_in_place(kelv), -0.01*ones(nsurf))
        
        for i = 1:nsurf
            curfield.tev[length(curfield.tev)][i].s = soln.zero[i]
        end
        
        #Calculate adot
        update_a2a3adot(surf, dt)
        
        for i = 1:nsurf
            surf[i].a0prev[1] = surf[i].a0[1]
            for ia = 1:3
                surf[i].aprev[ia] = surf[i].aterm[ia]
            end
        end
        
        #Update rest of Fourier terms
        update_a2toan(surf)
        
        #Calculate bound vortex strengths
        update_bv(surf)
        
        wakeroll(surf, curfield, dt)
        
        cl, cd, cm = calc_forces(surf)
        
        for i = 1:nsurf
            mat[istep,1,i] = t
            mat[istep,2,i] = surf[i].kinem.alpha
            mat[istep,3,i] = surf[i].kinem.h
            mat[istep,4,i] = surf[i].kinem.u
            mat[istep,5,i] = surf[i].a0[1]
            mat[istep,6,i] = cl[i]
            mat[istep,7,i] = cd[i]
            mat[istep,8,i] = cm[i]
        end
    end
    
    mat, surf, curfield
end

function lautat_wakeroll(surf::TwoDSurf, curfield::TwoDFlowField, nsteps::Int64 = 500, dtstar::Float64 = 0.015, delvort = DelVortDef(0, 0, 0.), mat = Array(Float64, 0, 8), kelv_enf = 0.)

    if (size(mat,1) > 0)
        t = mat[end,1]
    else
        t = 0.
    end

    mat = mat'

    dt = dtstar*surf.c/surf.uref

    #Intialise flowfield
    for istep = 1:nsteps
        #Udpate current time
        t = t + dt

        #Update kinematic parameters
        update_kinem(surf, t)

        #Update bound vortex positions
        update_boundpos(surf, dt)

        #Add a TEV with dummy strength
        place_tev(surf,curfield,dt)

        kelv = KelvinCondition(surf,curfield)
        #Solve for TEV strength to satisfy Kelvin condition
        #curfield.tev[length(curfield.tev)].s = secant_method(kelv, 0., -0.01)
        soln = nlsolve(not_in_place(kelv), [-0.01])
        curfield.tev[length(curfield.tev)].s = soln.zero[1]

        #Update adot
        update_a2a3adot(surf,dt)

        #Check for LEV and shed if yes
        #Set previous values of aterm to be used for derivatives in next time step
        surf.a0prev[1] = surf.a0[1]
        for ia = 1:3
            surf.aprev[ia] = surf.aterm[ia]
        end

        #Update rest of Fourier terms
        update_a2toan(surf)

        #Calculate bound vortex strengths
        update_bv(surf)

        #Remove vortices that are far away from airfoil
        if (delvort.flag == 1)
            if length(curfield.tev) > delvort.limit
                if (sqrt((curfield.tev[1].x- surf.bnd_x[div(surf.ndiv,2)])^2 + (curfield.tev[1].z- surf.bnd_z[div(surf.ndiv,2)])^2) > delvort.dist*surf.c)
                    kelv_enf = kelv_enf + curfield.tev[1].s
                    for i = 1:length(curfield.tev)-1
                        curfield.tev[i] = curfield.tev[i+1]
                    end
                    pop!(curfield.tev)
                end
            end
            if length(curfield.lev) > delvort.limit
                if (sqrt((curfield.lev[1].x- surf.bnd_x[div(surf.ndiv,2)])^2 + (curfield.lev[1].z- surf.bnd_z[div(surf.ndiv,2)])^2) > delvort.dist*surf.c)
                    kelv_enf = kelv_enf + curfield.lev[1].s
                    for i = 1:length(curfield.lev)-1
                        curfield.lev[i] = curfield.lev[i+1]
                    end
                    pop!(curfield.lev)
                end
            end
        end


        wakeroll(surf, curfield, dt)

        cl, cd, cm = calc_forces(surf)
        mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, surf.a0[1], cl, cd, cm])
    end
    mat = mat'
    mat, surf, curfield, kelv_enf

end

function lautat(surf::TwoDSurfLV, curfield::TwoDFlowField, nsteps::Int64 = 500, dtstar::Float64 = 0.015, delvort = DelVortDef(0, 0, 0.), mat = Array(Float64, 0, 8), kelv_enf = 0.)
    #Lumped vortex method
    if (size(mat,1) > 0)
        t = mat[end,1]
    else
        t = 0.
    end

    mat = mat'

    dt = dtstar*surf.c/surf.uref

    #temp variables used in simualtion
    uw = zeros(2)
    uwx = zeros(surf.npanel)
    uwz = zeros(surf.npanel)
    utx = zeros(surf.npanel)
    utz = zeros(surf.npanel)

    cp = zeros(surf.npanel)
    sum_gam = zeros(surf.npanel)
    sum_gam_prev = zeros(surf.npanel)
    
    #Intialise flowfield
    for istep = 1:nsteps
        #Udpate current time
        t = t + dt

        #Update kinematic parameters
        update_kinem(surf, t)

        #Update intertial coordinates - freestream is horizontal
        #velocity and plunge velocity is vertical velocity.
        
        surf.X0[1] = surf.X0[1] - surf.kinem.u*dt
        surf.X0[2] = surf.X0[2] + surf.kinem.hdot*dt
                
        # Update Global to Local Transformation
        surf.tlg[1,1] = cos(surf.kinem.alpha)
        surf.tlg[1,2] = -sin(surf.kinem.alpha)
        surf.tlg[2,1] = sin(surf.kinem.alpha)
        surf.tlg[2,2] = cos(surf.kinem.alpha)
        
        # Update Local to Global Transformation
        surf.tgl[1,1] = cos(surf.kinem.alpha)
        surf.tgl[1,2] = sin(surf.kinem.alpha)
        surf.tgl[2,1] = -sin(surf.kinem.alpha)
        surf.tgl[2,2] = cos(surf.kinem.alpha)
                        
        #Add a TEV with dummy strength
        place_tev(surf,curfield,dt)
        
        # Update Collocation Points and Vortex Points to Inertial Reference
        for i=1:surf.npanel
            (surf.lv[i].xc_I, surf.lv[i].zc_I) = surf.tgl*[surf.lv[i].xc; surf.lv[i].zc] + [surf.X0[1];surf.X0[2]]
            (surf.lv[i].xv_I, surf.lv[i].zv_I) = surf.tgl*[surf.lv[i].xv; surf.lv[i].zv] + [surf.X0[1];surf.X0[2]]
        end

        # Normal Velocity is a combination of self-induced velocity, kinematic velocity and wake induced velocity.
        
        #Update IC - relative position of shed TEV changes with Ansari's 1/3 law
        j = surf.npanel + 1
        for i = 1:surf.npanel
            (xloc, zloc) = surf.tlg*[curfield.tev[length(curfield.tev)].x; curfield.tev[length(curfield.tev)].z] - [surf.X0[1]; surf.X0[2]]
            
            dummytev = TwoDVort(xloc, zloc, 1.0, 0.02*surf.c, 0.0, 0.0)
            ui, wi = ind_vel([dummytev;], surf.lv[i].xc, surf.lv[i].zc)
            uw[1] = ui[1]
            uw[2] = wi[1]
            surf.IC[i,j] = dot(uw,[surf.lv[i].nx; surf.lv[i].nz])
        end
        
        # Wake Induced Velocity and kinematics induced velocity are know quantities and can be summed together to form the RHS
        # Establish RHS

        #Update induced velocities on airfoil - except last wake vortex
        update_indbound(surf, curfield)
        
        # Kinematics induced Velocity
        for i = 1:surf.npanel
            (utx[i], utz[i]) = surf.tlg*[+surf.kinem.u; -surf.kinem.hdot] + [-surf.kinem.alphadot*surf.lv[i].zc; surf.kinem.alphadot*(surf.lv[i].xc - surf.pvt)]
        end
        
        for i = 1:surf.npanel
            (uwx[i], uwz[i]) = [utx[i]; utz[i]] + [surf.lv[i].uind; surf.lv[i].wind]
        end

        #Right-Hand Side
        for i = 1:surf.npanel
            surf.rhs[i] = dot(-[uwx[i]; uwz[i]], [surf.lv[i].nx; surf.lv[i].nz])
        end
        surf.rhs[surf.npanel+1] = surf.gamma_prev[1]
        
        circsol = surf.IC \ surf.rhs
        surf.gamma_prev[1] = surf.gamma[1]
        surf.gamma[1] = 0
        for i = 1:surf.npanel
            surf.lv[i].s = circsol[i]
            surf.gamma[1] = surf.gamma[1] + circsol[i]
        end
        
        curfield.tev[length(curfield.tev)].s = circsol[surf.npanel+1]

        for i = 1:surf.npanel
            sum_gam_prev[i] = sum_gam[i]
        end
        
        sum_gam[1] = surf.lv[1].s
        for i = 2:surf.npanel
            sum_gam[i] = sum_gam[i-1] + surf.lv[i].s
        end

        # Force calculations
        cl = 0
        cm = 0
        cd = 0
        for i = 1:surf.npanel
            cp[i] = (1./(surf.uref*surf.uref))*(dot([utx[i] + uwx[i]; utz[i] + uwz[i]], [surf.lv[i].tx; surf.lv[i].tz])*surf.lv[i].s/surf.lv[i].l + (sum_gam[i]-sum_gam_prev[i])/dt)
            cl = cl + cp[i]*surf.lv[i].l*surf.lv[i].nz/surf.c
            cm = cm - cp[i]*surf.lv[i].l*surf.lv[i].nz*(surf.lv[i].xv - surf.pvt*surf.c)/(surf.c*surf.c)
            cd = cd + uwz[i]*surf.lv[i].s + surf.lv[i].l*surf.lv[i].nx*(sum_gam[i] - sum_gam_prev[i])/dt
        end
        cd = 2*cd/(surf.uref*surf.uref*surf.c)
        mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, surf.gamma[1], cl, cd, cm])
    end
    mat = mat'
    mat, surf, curfield, kelv_enf
end    
    
function lautat_wakeroll(surf::TwoDSurfLV, curfield::TwoDFlowField, nsteps::Int64 = 500, dtstar::Float64 = 0.015, delvort = DelVortDef(0, 0, 0.), mat = Array(Float64, 0, 8), kelv_enf = 0.)
    #Lumped vortex method
    if (size(mat,1) > 0)
        t = mat[end,1]
    else
        t = 0.
    end

    mat = mat'

    dt = dtstar*surf.c/surf.uref

    #temp variables used in simualtion
    uw = zeros(2)
    uwx = zeros(surf.npanel)
    uwz = zeros(surf.npanel)
    utx = zeros(surf.npanel)
    utz = zeros(surf.npanel)

    cp = zeros(surf.npanel)
    sum_gam = zeros(surf.npanel)
    sum_gam_prev = zeros(surf.npanel)
    
    #Intialise flowfield
    for istep = 1:nsteps
        #Udpate current time
        t = t + dt

        #Update kinematic parameters
        update_kinem(surf, t)

        #Update intertial coordinates - freestream is horizontal
        #velocity and plunge velocity is vertical velocity.
        
        surf.X0[1] = surf.X0[1] - surf.kinem.u*dt
        surf.X0[2] = surf.X0[2] + surf.kinem.hdot*dt
                
        # Update Global to Local Transformation
        surf.tlg[1,1] = cos(surf.kinem.alpha)
        surf.tlg[1,2] = -sin(surf.kinem.alpha)
        surf.tlg[2,1] = sin(surf.kinem.alpha)
        surf.tlg[2,2] = cos(surf.kinem.alpha)
        
        # Update Local to Global Transformation
        surf.tgl[1,1] = cos(surf.kinem.alpha)
        surf.tgl[1,2] = sin(surf.kinem.alpha)
        surf.tgl[2,1] = -sin(surf.kinem.alpha)
        surf.tgl[2,2] = cos(surf.kinem.alpha)
                        
        #Add a TEV with dummy strength
        place_tev(surf,curfield,dt)
        
        # Update Collocation Points and Vortex Points to Inertial Reference
        for i=1:surf.npanel
            (surf.lv[i].xc_I, surf.lv[i].zc_I) = surf.tgl*[surf.lv[i].xc; surf.lv[i].zc] + [surf.X0[1];surf.X0[2]]
            (surf.lv[i].xv_I, surf.lv[i].zv_I) = surf.tgl*[surf.lv[i].xv; surf.lv[i].zv] + [surf.X0[1];surf.X0[2]]
        end

        # Normal Velocity is a combination of self-induced velocity, kinematic velocity and wake induced velocity.
        
        #Update IC - relative position of shed TEV changes with Ansari's 1/3 law
        j = surf.npanel + 1
        for i = 1:surf.npanel
            (xloc, zloc) = surf.tlg*[curfield.tev[length(curfield.tev)].x; curfield.tev[length(curfield.tev)].z] - [surf.X0[1]; surf.X0[2]]
            
            dummytev = TwoDVort(xloc, zloc, 1.0, 0.02*surf.c, 0.0, 0.0)
            ui, wi = ind_vel([dummytev;], surf.lv[i].xc, surf.lv[i].zc)
            uw[1] = ui[1]
            uw[2] = wi[1]
            surf.IC[i,j] = dot(uw,[surf.lv[i].nx; surf.lv[i].nz])
        end
        
        # Wake Induced Velocity and kinematics induced velocity are know quantities and can be summed together to form the RHS
        # Establish RHS

        #Update induced velocities on airfoil - except last wake vortex
        update_indbound(surf, curfield)
        
        # Kinematics induced Velocity
        for i = 1:surf.npanel
            (utx[i], utz[i]) = surf.tlg*[+surf.kinem.u; -surf.kinem.hdot] + [-surf.kinem.alphadot*surf.lv[i].zc; surf.kinem.alphadot*(surf.lv[i].xc - surf.pvt)]
        end
        
        for i = 1:surf.npanel
            (uwx[i], uwz[i]) = [utx[i]; utz[i]] + [surf.lv[i].uind; surf.lv[i].wind]
        end

        #Right-Hand Side
        for i = 1:surf.npanel
            surf.rhs[i] = dot(-[uwx[i]; uwz[i]], [surf.lv[i].nx; surf.lv[i].nz])
        end
        surf.rhs[surf.npanel+1] = surf.gamma_prev[1]
        
        circsol = surf.IC \ surf.rhs
        surf.gamma_prev[1] = surf.gamma[1]
        surf.gamma[1] = 0
        for i = 1:surf.npanel
            surf.lv[i].s = circsol[i]
            surf.gamma[1] = surf.gamma[1] + circsol[i]
        end
        
        curfield.tev[length(curfield.tev)].s = circsol[surf.npanel+1]

        for i = 1:surf.npanel
            sum_gam_prev[i] = sum_gam[i]
        end
        
        sum_gam[1] = surf.lv[1].s
        for i = 2:surf.npanel
            sum_gam[i] = sum_gam[i-1] + surf.lv[i].s
        end

        #Wake rollup
        wakeroll(surf, curfield, dt)
        
        
        # Force calculations
        cl = 0
        cm = 0
        cd = 0
        for i = 1:surf.npanel
            cp[i] = (1./(surf.uref*surf.uref))*(dot([utx[i] + uwx[i]; utz[i] + uwz[i]], [surf.lv[i].tx; surf.lv[i].tz])*surf.lv[i].s/surf.lv[i].l + (sum_gam[i]-sum_gam_prev[i])/dt)
            cl = cl + cp[i]*surf.lv[i].l*surf.lv[i].nz/surf.c
            cm = cm - cp[i]*surf.lv[i].l*surf.lv[i].nz*(surf.lv[i].xv - surf.pvt*surf.c)/(surf.c*surf.c)
            cd = cd + uwz[i]*surf.lv[i].s + surf.lv[i].l*surf.lv[i].nx*(sum_gam[i] - sum_gam_prev[i])/dt
        end
        cd = 2*cd/(surf.uref*surf.uref*surf.c)
        mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, surf.gamma[1], cl, cd, cm])
    end
    mat = mat'
    mat, surf, curfield, kelv_enf
end    
    


function lautat_wakeroll_more(surf::TwoDSurf, curfield::TwoDFlowField, nsteps::Int64 = 500, dtstar::Float64 = 0.015, delvort = DelVortDef(0, 0, 0.), mat = Array(Float64, 0, 11), kelv_enf = 0.)

    if (size(mat,1) > 0)
        t = mat[end,1]
    else
        t = 0.
    end

    mat = mat'

    dt = dtstar*surf.c/surf.uref

    #Intialise flowfield
    for istep = 1:nsteps
        #Udpate current time
        t = t + dt

        #Update kinematic parameters
        update_kinem(surf, t)

        #Update bound vortex positions
        update_boundpos(surf, dt)

        #Add a TEV with dummy strength
        place_tev(surf,curfield,dt)

        kelv = KelvinCondition(surf,curfield)
        #Solve for TEV strength to satisfy Kelvin condition
        #curfield.tev[length(curfield.tev)].s = secant_method(kelv, 0., -0.01)
        soln = nlsolve(not_in_place(kelv), [-0.01])
        curfield.tev[length(curfield.tev)].s = soln.zero[1]

        #Update adot
        update_a2a3adot(surf,dt)

        #Check for LEV and shed if yes
        #Set previous values of aterm to be used for derivatives in next time step
        surf.a0prev[1] = surf.a0[1]
        for ia = 1:3
            surf.aprev[ia] = surf.aterm[ia]
        end

        #Update rest of Fourier terms
        update_a2toan(surf)

        #Calculate bound vortex strengths
        update_bv(surf)

        #Remove vortices that are far away from airfoil
        if (delvort.flag == 1)
            if length(curfield.tev) > delvort.limit
                if (sqrt((curfield.tev[1].x- surf.bnd_x[div(surf.ndiv,2)])^2 + (curfield.tev[1].z- surf.bnd_z[div(surf.ndiv,2)])^2) > delvort.dist*surf.c)
                    kelv_enf = kelv_enf + curfield.tev[1].s
                    for i = 1:length(curfield.tev)-1
                        curfield.tev[i] = curfield.tev[i+1]
                    end
                    pop!(curfield.tev)
                end
            end
            if length(curfield.lev) > delvort.limit
                if (sqrt((curfield.lev[1].x- surf.bnd_x[div(surf.ndiv,2)])^2 + (curfield.lev[1].z- surf.bnd_z[div(surf.ndiv,2)])^2) > delvort.dist*surf.c)
                    kelv_enf = kelv_enf + curfield.lev[1].s
                    for i = 1:length(curfield.lev)-1
                        curfield.lev[i] = curfield.lev[i+1]
                    end
                    pop!(curfield.lev)
                end
            end
        end


        wakeroll(surf, curfield, dt)

        cl, cd, cm, bc, cn, cs = calc_forces_more(surf)
        mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, surf.a0[1], cl, cd, cm, bc, cn, cs])
    end
    mat = mat'
    mat, surf, curfield, kelv_enf

end


function theodorsen(theo::TheoDef)
    # Inputs:
    # h_amp = h_amp/c, plunge amplitude, positive up
    # alpha_amp = pitch amplitude, positive LE up
    # phi = phase angle by which pitch leads plunge
    # pvt = pvt/c, non-dimensional chordwise rotation point (0 to 1)
    # alpha_mean = mean angle of attack of chord line
    # alpha_zl = zero-lift angle of attack of airfoil
    # reduced frequency

    # Motion definitions
    # h = h_amp*exp(i*wt)
    # alpha = alpha_mean + alpha_amp*exp(i*(wt + phi))

    #define a
    a = (theo.pvt-0.5)/0.5

    wt = [0:2*pi/360:2*pi;]

    #Theodorsen function
    C = besselh(1,2,theo.k)./(besselh(1,2,theo.k) + im*besselh(0,2,theo.k))

    # steady-state Cl
    Cl_ss = 2*pi*(theo.alpha_mean - theo.alpha_zl)

    # plunge contribution
    Cl_h = 2*pi*theo.k^2*theo.h_amp*exp(im*wt) - im*4*pi*theo.k*C*theo.h_amp*exp(im*wt)

    # pitch contribution
    Cl_alpha = (im*pi*theo.k + pi*theo.k^2*a)*theo.alpha_amp*exp(im*(wt+theo.phi)) + (1 + im*theo.k*(0.5-a))*2*pi*C*theo.alpha_amp*exp(im*(wt+theo.phi))

    # total contributions
    Cl_tot = Cl_ss + Cl_h + Cl_alpha

    return wt/(2*pi), Cl_h, Cl_alpha, Cl_tot

end

function theodorsen(theo::TheoDefwFlap)
    # Inputs:
    # h_amp = h_amp/c, plunge amplitude, positive up
    # alpha_amp = pitch amplitude, positive LE up
    # beta_amp = flap amplitude, positive LE up
    # phi = phase angle by which pitch leads plunge
    # psi = phase angle by which flap rotation leads plunge
    # pvt = pvt/c, non-dimensional chordwise rotation point (0 to 1)
    # alpha_mean = mean angle of attack of chord line
    # alpha_zl = zero-lift angle of attack of airfoil
    # xf = xf/c, non-dimensional flap location

    # Motion definitions
    # h = h_amp*exp(i*wt)
    # alpha = alpha_mean + alpha_amp*exp(i*(wt + phi))
    # beta =  beta_amp*exp(i*(wt + psi))

    #define a and c
    a = (theo.pvt-0.5)/0.5
    c = (theo.xf-0.5)/0.5

    #Define the required coefficients
    T1 = -(2+c*c)*sqrt(1-c*c)/3 + c*acos(c)
    T4 = c*sqrt(1-c*c) - acos(c)
    T11 = (2-c)*sqrt(1-c*c) + (1-2*c)*acos(c)
    T12 = (2+c)*sqrt(1-c*c) + (1-2*c)*acos(c)
    T2 = T4*(T11+T12)
    T3 = -(1-c*c)*(5*c*c+4)/8 + c*(7+2*c*c)*sqrt(1-c*c)*acos(c)/4 - (1/8+c*c)*acos(c)*acos(c)
    T5 = -(1-c*c)+2*c*sqrt(1-c*c)*acos(c)-acos(c)*acos(c)
    T6 = T2
    T7 = c*(7+2*c*c)*sqrt(1-c*c)/8 - (1/8+c*c)*acos(c)
    T8 = -(1+2*c*c)*sqrt(1-c*c)/3 + c*acos(c)
    T9 = ((1-c*c)^(3/2)/3 + a*T4)/2
    T10 = sqrt(1-c*c) + acos(c)
    T13 = -0.5*(T7+(c-a)*T1)
    T14 = 1/16 + a*c/2
    T15 = T4 + T10
    T16 = T1 - T8 -(c-a)*T4 + T11/2
    T17 = -2*T9 - T1 + (a-1/2)*T4
    T18 = T5 - T4*T10
    T19 = T4*T11
    T20 = T10 - 2*sqrt(1-c*c)

    wt = [0:2*pi/360:2*pi;]

    #Theodorsen function
    C = besselh(1,2,theo.k)./(besselh(1,2,theo.k) + im*besselh(0,2,theo.k))

    # steady-state Cl
    Cl_ss = 2*pi*(theo.alpha_mean - theo.alpha_zl)

    #Derived in a Julia notebook (Available in the UNSflow package)
    # Lift
    Cl_h = 2*pi*theo.h_amp*theo.k*(-2*im*C+theo.k)*exp(im*wt)
    Cl_alpha = -pi*theo.alpha_amp*(2*C*(im*theo.k*(a-0.5)-1) - theo.k*(a*theo.k+im))*exp(im*(theo.phi+wt))
    Cl_beta = theo.beta_amp*(C*(im*theo.k*T11+2*T10) + theo.k*(theo.k*T1-im*T4))*exp(im*(wt+theo.psi))
    # -------------------------------------
    # Pitching moment
    Cmal_h = -pi*theo.h_amp*theo.k*(2*im*C*(a + 0.5) - a*theo.k)*exp(im*wt)
    Cmal_al = -0.5*pi*theo.alpha_amp*(2*C*(a + 0.5)*(im*theo.k*(a + 1) - 1) - theo.k*(theo.k*(a^2 + 0.125) + im*(a - 0.5)))*exp(im*(theo.phi + wt))
    Cmal_be = 0.5*theo.beta_amp*(C*(a+0.5)*(im*theo.k*T11* + 2*T10) + 2*T13*theo.k^2 - im*theo.k*T16 - T15)*exp(im*(wt + theo.psi))
    # -------------------------------------
    # Hinge moment
    Cmbe_h = theo.h_amp*theo.k*(im*C*T12 + theo.k*T1)*exp(im*wt)
    Cmbe_al = 0.5*theo.alpha_amp*(C*T12*(im*theo.k*(a - 0.5) - 1) + theo.k*(2*theo.k*T13 - im*T17))*exp(im*(theo.phi + wt))
    Cmbe_be = -(theo.beta_amp/(4*pi))*(C*T12*(im*theo.k*T11 + 2*T10) + 2*T3*theo.k^2 - im*theo.k*T19 +2*T18)*exp(im*(wt + theo.psi))

    # total contributions
    Cl_tot = Cl_ss + Cl_h + Cl_alpha + Cl_beta
    Cmal_tot = Cmal_h + Cmal_al + Cmal_be
    Cmbe_tot = Cmbe_h + Cmbe_al + Cmbe_be

#   return wt/(2*pi), Cl_ss, Cl_h, Cl_alpha, Cl_beta, Cl_tot
#   return wt/(2*pi), Cmal_h, Cmal_al, Cmal_be, Cmal_tot
#   return wt/(2*pi), Cmbe_h, Cmbe_al, Cmbe_be, Cmbe_tot
    return wt/(2*pi), Cl_tot, Cmal_tot, Cmbe_tot
end

function ldvm_varU(surf::TwoDSurf, curfield::TwoDFlowField, recalib::Vector{Float64}, lespcalib::Vector{Float64}, re_ref::Float64, nsteps::Int64 = 500, dtstar::Float64 = 0.015; delvort = DelVortDef(0, 0, 0.), mat = Array(Float64, 0, 10), kelv_enf = 0., writefile = "Nil")

    lespspl = Spline1D(recalib, lespcalib)

    if (size(mat,1) > 0)
        t = mat[end,1]
    else
        t = 0.
    end
    
    #mat = zeros(nsteps,11)
    mat = mat'
    
    dt = dtstar*surf.c/surf.uref
    
    #If required, write an archive file
    if writefile != "Nil"
        jldopen(writefile, "w") do file
            write(file, "dt",  dt)
            write(file, "dtstar",  dt)
            write(file, "Insurf",  surf)
            write(file, "Infield", curfield)
            write(file, "delvort", delvort)
            write(file, "nsteps", nsteps)
            cl, cd, cm, gamma, cn, cs, cnc, cncc, nonl, cm_n, cm_pvt, nonl_m = calc_forces_more(surf)
            g = g_create(file, "Init")
            g = write_stamp(surf, curfield, t, kelv_enf, g)
        end
    end
    
    
    #Intialise flowfield
    for istep = 1:nsteps
        #Udpate current time
        t = t + dt
        
        #Update kinematic parameters
        update_kinem(surf, t)
        
        
        #Update flow field parameters if any
        update_externalvel(curfield, t)
        
        #Update bound vortex positions
        update_boundpos(surf, dt)
        
        #Add a TEV with dummy strength
        place_tev(surf,curfield,dt)
        
        kelv = KelvinCondition(surf,curfield)
        #Solve for TEV strength to satisfy Kelvin condition
        #curfield.tev[length(curfield.tev)].s = secant_method(kelv, 0., -0.01)
        soln = nlsolve(not_in_place(kelv), [-0.01])
        curfield.tev[length(curfield.tev)].s = soln.zero[1]
        
        #Check for LESP condition
        #Update values with converged value of shed tev
        #Update incduced velocities on airfoil
        update_indbound(surf, curfield)
        
        #Calculate downwash
        update_downwash(surf, [curfield.u[1],curfield.w[1]])
        
        #Calculate first two fourier coefficients
        update_a0anda1(surf)
        
        
        le_vel_x = surf.kinem.u*cos(surf.kinem.alpha) + surf.kinem.hdot*sin(surf.kinem.alpha) + surf.uind[1]
        le_vel_z = surf.kinem.u*sin(surf.kinem.alpha) - surf.kinem.alphadot*surf.pvt*surf.c - surf.kinem.hdot*cos(surf.kinem.alpha) + surf.wind[1]
        vmag = sqrt(le_vel_x*le_vel_x+le_vel_z*le_vel_z)
        re_le = re_ref*vmag
        
        surf.lespcrit[1] = Dierckx.evaluate(lespspl, re_le)
        
        lesp = surf.a0[1]
        
        #Update adot
        update_a2a3adot(surf,dt)
        
        #2D iteration if LESP_crit is exceeded
        if (abs(lesp)>surf.lespcrit[1])
            #Remove the previous tev
            pop!(curfield.tev)
            #Add a TEV with dummy strength
            place_tev(surf,curfield,dt)
            
            #Add a LEV with dummy strength
            place_lev(surf,curfield,dt)
            
            kelvkutta = KelvinKutta(surf,curfield)
            #Solve for TEV and LEV strengths to satisfy Kelvin condition and Kutta condition at leading edge
            
            soln = nlsolve(not_in_place(kelvkutta), [-0.01; 0.01])
            (curfield.tev[length(curfield.tev)].s, curfield.lev[length(curfield.lev)].s) = soln.zero[1], soln.zero[2]
            
            surf.levflag[1] = 1
        else
            surf.levflag[1] = 0
        end
        
        
        #Update rest of Fourier terms
        update_a2toan(surf)
        
        #Set previous values of aterm to be used for derivatives in next time step
        surf.a0prev[1] = surf.a0[1]
        for ia = 1:3
            surf.aprev[ia] = surf.aterm[ia]
        end
        
        #Calculate bound vortex strengths
        update_bv(surf)
        
        #Remove vortices that are far away from airfoil
        if (delvort.flag == 1)
            if length(curfield.tev) > delvort.limit
                if (sqrt((curfield.tev[1].x- surf.bnd_x[div(surf.ndiv,2)])^2 + (curfield.tev[1].z- surf.bnd_z[div(surf.ndiv,2)])^2) > delvort.dist*surf.c)
                    kelv_enf = kelv_enf + curfield.tev[1].s
                    for i = 1:length(curfield.tev)-1
                        curfield.tev[i] = curfield.tev[i+1]
                    end
                    pop!(curfield.tev)
                end
            end
            if length(curfield.lev) > delvort.limit
                if (sqrt((curfield.lev[1].x- surf.bnd_x[div(surf.ndiv,2)])^2 + (curfield.lev[1].z- surf.bnd_z[div(surf.ndiv,2)])^2) > delvort.dist*surf.c)
                    kelv_enf = kelv_enf + curfield.lev[1].s
                    for i = 1:length(curfield.lev)-1
                        curfield.lev[i] = curfield.lev[i+1]
                    end
                    pop!(curfield.lev)
                end
            end
        end
        wakeroll(surf, curfield, dt)
        
        #cl, cd, cm, cn, cs = calc_forces(surf)
        
        if writefile == "Nil"
            cl, cd, cm = calc_forces(surf)
        else
            cl, cd, cm, gamma, cn, cs, cnc, cncc, nonl, cm_n, cm_pvt, nonl_m = calc_forces_more(surf)
        end
        
        mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, surf.a0[1], cl, cd, cm, surf.lespcrit[1], re_le])
        
        #If required, write an archive file
        if writefile != "Nil"
            jldopen(writefile, "r+") do file
                g = g_create(file, "t$istep")
                g = write_stamp(surf, curfield, t, kelv_enf, g)
            end
        end
    end

mat = mat'
mat, surf, curfield, kelv_enf

end

        
function ldvm(surf::TwoDSurf, curfield::TwoDFlowField, nsteps::Int64 = 500, dtstar::Float64 = 0.015; delvort = DelVortDef(0, 0, 0.), mat = Array(Float64, 0, 8), kelv_enf = 0., writefile = "Nil")

    if (size(mat,1) > 0)
        t = mat[end,1]
    else
        t = 0.
    end

    #mat = zeros(nsteps,11)
    mat = mat'

    dt = dtstar*surf.c/surf.uref
   
    #If required, write an archive file
    if writefile != "Nil"
        jldopen(writefile, "w") do file
            write(file, "dt",  dt)
            write(file, "dtstar",  dt)
            write(file, "Insurf",  surf)
            write(file, "Infield", curfield)
            write(file, "delvort", delvort)
            write(file, "nsteps", nsteps)
            cl, cd, cm, gamma, cn, cs, cnc, cncc, nonl, cm_n, cm_pvt, nonl_m = calc_forces_more(surf)
            g = g_create(file, "Init")
            g = write_stamp(surf, curfield, t, kelv_enf, g)
        end
    end
    
    
    #Intialise flowfield
    for istep = 1:nsteps
        #Udpate current time
        t = t + dt

        #Update kinematic parameters
        update_kinem(surf, t)

        #Update flow field parameters if any
        update_externalvel(curfield, t)

        #Update bound vortex positions
        update_boundpos(surf, dt)

        #Add a TEV with dummy strength
        place_tev(surf,curfield,dt)

        kelv = KelvinCondition(surf,curfield)
        #Solve for TEV strength to satisfy Kelvin condition
        #curfield.tev[length(curfield.tev)].s = secant_method(kelv, 0., -0.01)
        soln = nlsolve(not_in_place(kelv), [-0.01])
        curfield.tev[length(curfield.tev)].s = soln.zero[1]

        #Check for LESP condition
        #Update values with converged value of shed tev
        #Update incduced velocities on airfoil
        update_indbound(surf, curfield)

        #Calculate downwash
        update_downwash(surf, [curfield.u[1],curfield.w[1]])

        #Calculate first two fourier coefficients
        update_a0anda1(surf)

        lesp = surf.a0[1]
        
        #Update adot
        update_a2a3adot(surf,dt)

        #2D iteration if LESP_crit is exceeded
        if (abs(lesp)>surf.lespcrit[1])
            #Remove the previous tev
            pop!(curfield.tev)
            #Add a TEV with dummy strength
            place_tev(surf,curfield,dt)

            #Add a LEV with dummy strength
            place_lev(surf,curfield,dt)

            kelvkutta = KelvinKutta(surf,curfield)
            #Solve for TEV and LEV strengths to satisfy Kelvin condition and Kutta condition at leading edge

            soln = nlsolve(not_in_place(kelvkutta), [-0.01; 0.01])
            (curfield.tev[length(curfield.tev)].s, curfield.lev[length(curfield.lev)].s) = soln.zero[1], soln.zero[2]

            surf.levflag[1] = 1
        else
            surf.levflag[1] = 0
        end


        #Update rest of Fourier terms
        update_a2toan(surf)

        #Set previous values of aterm to be used for derivatives in next time step
        surf.a0prev[1] = surf.a0[1]
        for ia = 1:3
            surf.aprev[ia] = surf.aterm[ia]
        end

        #Calculate bound vortex strengths
        update_bv(surf)

        #Remove vortices that are far away from airfoil
        if (delvort.flag == 1)
            if length(curfield.tev) > delvort.limit
                if (sqrt((curfield.tev[1].x- surf.bnd_x[div(surf.ndiv,2)])^2 + (curfield.tev[1].z- surf.bnd_z[div(surf.ndiv,2)])^2) > delvort.dist*surf.c)
                    kelv_enf = kelv_enf + curfield.tev[1].s
                    for i = 1:length(curfield.tev)-1
                        curfield.tev[i] = curfield.tev[i+1]
                    end
                    pop!(curfield.tev)
                end
            end
            if length(curfield.lev) > delvort.limit
                if (sqrt((curfield.lev[1].x- surf.bnd_x[div(surf.ndiv,2)])^2 + (curfield.lev[1].z- surf.bnd_z[div(surf.ndiv,2)])^2) > delvort.dist*surf.c)
                    kelv_enf = kelv_enf + curfield.lev[1].s
                    for i = 1:length(curfield.lev)-1
                        curfield.lev[i] = curfield.lev[i+1]
                    end
                    pop!(curfield.lev)
                end
            end
        end
        wakeroll(surf, curfield, dt)

        #cl, cd, cm, cn, cs = calc_forces(surf)

        if writefile == "Nil"
            cl, cd, cm = calc_forces(surf)
        else
            cl, cd, cm, gamma, cn, cs, cnc, cncc, nonl, cm_n, cm_pvt, nonl_m = calc_forces_more(surf)
        end

        mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, surf.a0[1], cl, cd, cm])

        #If required, write an archive file
        if writefile != "Nil"
            jldopen(writefile, "r+") do file
                g = g_create(file, "t$istep")
                g = write_stamp(surf, curfield, t, kelv_enf, g)
            end
        end
    end

    mat = mat'
    mat, surf, curfield, kelv_enf
end

function ldvm_klb(surf::TwoDSurf, curfield::TwoDFlowField, nsteps::Int64 = 500, dtstar::Float64 = 0.015; delvort = DelVortDef(0, 0, 0.), mat = Array(Float64, 0, 10), kelv_enf = 0., writefile = "Nil", tf = 3.0, alpha1 = 7*pi/180, S1 = 3.0*pi/180, S2 = 2.3*pi/180)
    
    if (size(mat,1) > 0)
        t = mat[end,1]
    else
        t = 0.
    end
    
    #mat = zeros(nsteps,11)
    mat = mat'
    
    dt = dtstar*surf.c/surf.uref
    
    #If required, write an archive file
    if writefile != "Nil"
        jldopen(writefile, "w") do file
            write(file, "dt",  dt)
            write(file, "dtstar",  dt)
            write(file, "Insurf",  surf)
            write(file, "Infield", curfield)
            write(file, "delvort", delvort)
            write(file, "nsteps", nsteps)
            cl, cd, cm, gamma, cn, cs, cnc, cncc, nonl, cm_n, cm_pvt, nonl_m = calc_forces_more(surf)
            g = g_create(file, "Init")
            g = write_stamp(surf, curfield, t, kelv_enf, g)
        end
    end
    
    dfn = 0.
    dfn_pr = 0.
    fsep = 1.
    f0sep = 1.
    f0sep_pr = 1.
    
    for istep = 1:nsteps
        #Udpate current time
        t = t + dt
        
        #Update kinematic parameters
        update_kinem(surf, t)
        
        #Update flow field parameters if any
        update_externalvel(curfield, t)

        #Update bound vortex positions
        update_boundpos(surf, dt)
        
        #Add a TEV with dummy strength
        place_tev(surf,curfield,dt)
        
        kelv = KelvinCondition(surf,curfield)
        #Solve for TEV strength to satisfy Kelvin condition
        #curfield.tev[length(curfield.tev)].s = secant_method(kelv, 0., -0.01)
        soln = nlsolve(not_in_place(kelv), [-0.01])
        curfield.tev[length(curfield.tev)].s = soln.zero[1]
        
        #Check for LESP condition
        #Update values with converged value of shed tev
        #Update incduced velocities on airfoil
        update_indbound(surf, curfield)
        
        #Calculate downwash
        update_downwash(surf, [curfield.u[1],curfield.w[1]])
        
        #Calculate first two fourier coefficients
        update_a0anda1(surf)
        
        #Calculate fsep
        f0sep_pr = f0sep
        if abs(surf.a0[1]) < alpha1
            f0sep = 1 - 0.3*exp((abs(surf.a0[1]) - alpha1)/S1)
        else
            f0sep = 0.04 + 0.66*exp((alpha1 - abs(surf.a0[1]))/S2)
        end

        dfn_pr = dfn
        dfn = dfn_pr*exp(-2*dt/tf)+(f0sep - f0sep_pr)*exp(-dt/tf)
        fsep = f0sep - dfn
              
        lesp = surf.a0[1]
        
        #Update adot
        update_a2a3adot(surf,dt)

        #2D iteration if LESP_crit is exceeded
        if (abs(lesp)>surf.lespcrit[1])
            #Remove the previous tev
            pop!(curfield.tev)
            #Add a TEV with dummy strength
            place_tev(surf,curfield,dt)

            #Add a LEV with dummy strength
            place_lev(surf,curfield,dt)

            kelvkutta = KelvinKutta(surf,curfield)
            #Solve for TEV and LEV strengths to satisfy Kelvin condition and Kutta condition at leading edge

            soln = nlsolve(not_in_place(kelvkutta), [-0.01; 0.01])
            (curfield.tev[length(curfield.tev)].s, curfield.lev[length(curfield.lev)].s) = soln.zero[1], soln.zero[2]

            surf.levflag[1] = 1
        else
            surf.levflag[1] = 0
        end


        #Update rest of Fourier terms
        update_a2toan(surf)

        #Set previous values of aterm to be used for derivatives in next time step
        surf.a0prev[1] = surf.a0[1]
        for ia = 1:3
            surf.aprev[ia] = surf.aterm[ia]
        end

        #Calculate bound vortex strengths
        update_bv(surf)

        #Remove vortices that are far away from airfoil
        if (delvort.flag == 1)
            if length(curfield.tev) > delvort.limit
                if (sqrt((curfield.tev[1].x- surf.bnd_x[div(surf.ndiv,2)])^2 + (curfield.tev[1].z- surf.bnd_z[div(surf.ndiv,2)])^2) > delvort.dist*surf.c)
                    kelv_enf = kelv_enf + curfield.tev[1].s
                    for i = 1:length(curfield.tev)-1
                        curfield.tev[i] = curfield.tev[i+1]
                    end
                    pop!(curfield.tev)
                end
            end
            if length(curfield.lev) > delvort.limit
                if (sqrt((curfield.lev[1].x- surf.bnd_x[div(surf.ndiv,2)])^2 + (curfield.lev[1].z- surf.bnd_z[div(surf.ndiv,2)])^2) > delvort.dist*surf.c)
                    kelv_enf = kelv_enf + curfield.lev[1].s
                    for i = 1:length(curfield.lev)-1
                        curfield.lev[i] = curfield.lev[i+1]
                    end
                    pop!(curfield.lev)
                end
            end
        end
        wakeroll(surf, curfield, dt)

        #cl, cd, cm, cn, cs = calc_forces(surf)

        #if writefile == "Nil"
        cl, cd, cm = calc_forces(surf, fsep)
        #else
            #This is currently not in use - will need to be changed
         #   cl, cd, cm, gamma, cn, cs, cnc, cncc, nonl, cm_n, cm_pvt, nonl_m = calc_forces_more(surf, fsep)
        #end
        
        mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, surf.a0[1], cl, cd, cm, f0sep, fsep])
        
        #If required, write an archive file
        if writefile != "Nil"
            jldopen(writefile, "r+") do file
                g = g_create(file, "t$istep")
                g = write_stamp(surf, curfield, t, kelv_enf, g)
            end
        end
    end

mat = mat'
mat, surf, curfield, kelv_enf
end


function ldvm(surf::Vector{TwoDSurf}, curfield::TwoDFlowFieldMultSurf, nsteps::Int64 = 500; dtstar::Float64 = 0.015)

    #Deleting vortices is not currently supported
    #Resuming simulations is not currently supproted
    #Writefile is not currently supported
    #Move these to arguments when these are started
    
    #Some functions have been given shed_ind as an optional parameter with multiple dispatch - to be called for lautat and ldvm. This is confusing programming practice and should be changed. 
    
    
    nsurf = length(surf)
    mat = Array(Float64, nsteps, 8, nsurf)
    
    t = 0.
    
    dt = 100
    for i = 1:length(surf)
        dt = min(dt, dtstar*surf[i].c/surf[i].uref)
    end
    #Calculate tstep based on frequency as well.
    
    cl = zeros(nsurf)
    cd = zeros(nsurf)
    cm = zeros(nsurf)
    lesp = zeros(nsurf)
    
    eps = 1e-6
    
    #shed_ind should be a vector of vectors with rows corresponding to lev count and columns corresponding to surfaces.
    shed_ind = Vector{Int}[]
    
    #Intialise flowfield
    for istep = 1:nsteps
        #Udpate current time
        t = t + dt
        
        #Update kinematic parameters
        update_kinem(surf, t)
        
        #Update flow field parameters if any
        update_externalvel(curfield, t)
        
        #Update bound vortex positions
        update_boundpos(surf, dt)
        
        iter = 0
        while true
            iter += 1
            if iter > 1
                tev_prev = map(q->q.s, curfield.tev[end])
                pop!(curfield.tev)
            end
            
            #Add a TEV with dummy strength
            place_tev(surf, curfield, dt)
            
            for i = 1:nsurf
                kelv = KelvinConditionMultSurfSep(surf, curfield, shed_ind, i)
                soln = nlsolve(not_in_place(kelv), [-0.01])
                curfield.tev[end][i].s = soln.zero[1]
            end
            
            if iter > 1
                teverr = mean(abs(map(q->q.s, curfield.tev[end]) - tev_prev))
                #println(iter, teverr, map(q->q.s, curfield.tev[end]), tev_prev)
                if (teverr < eps)
                    #TEV conditions alone actually implies the others and would be sufficient
                    break
                end
                if (iter > 10)
                    println("1D convergence failed")
                    break
                end
            end
        end
        
        # #Check for LESP condition
        # #Update values with converged value of shed tev
        # #Update incduced velocities on airfoil
        update_indbound(surf, curfield, shed_ind)
        
        # #Calculate downwash
        update_downwash(surf, [curfield.u[1],curfield.w[1]])
        
        # #Calculate first two fourier coefficients
        update_a0anda1(surf)
        
        update_a2a3adot(surf, dt) 
        
        #2D iteration if LESP_crit is exceeded
        shedv = Int[]
        for i = 1:nsurf
            if (abs(surf[i].a0[1])>surf[i].lespcrit[1])
                push!(shedv, i)
            end
        end
        nshed = length(shedv)
        
        if nshed > 0
            iter = 0
            while true
                iter += 1
                
                # Remove the previous lev and tev if this is a further iteration
                if iter > 1
                    tev_prev = map(q->q.s, curfield.tev[end])
                    lev_prev = map(q->q.s, curfield.lev[end])
                    
                    pop!(curfield.tev)
                    pop!(curfield.lev)
                    pop!(shed_ind)
                else
                    pop!(curfield.tev)
                end
                push!(shed_ind, shedv)
                
                #Add a TEV with dummy strength
                place_tev(surf, curfield, dt)
            
                #Add a LEV with dummy strength
                place_lev(surf, curfield, dt, shed_ind)
                
                for i = 1:nsurf
                    if i in shed_ind[end]
                        kelvkutta = KelvinKuttaMultSurfSep(surf, curfield, shed_ind, i)
                        #Solve for TEV and LEV strengths to satisfy Kelvin condition and Kutta condition at leading edge
                    
                        soln = nlsolve(not_in_place(kelvkutta), [-0.01; 0.01])
                        
                        #Assign the solution
                        curfield.tev[end][i].s = soln.zero[1]
                        curfield.lev[end][i].s = soln.zero[2]
                    else
                        kelv= KelvinConditionMultSurfSep(surf, curfield, shed_ind, i)
                        soln = nlsolve(not_in_place(kelv), [-0.01])
                        
                        #Assign the solution
                        curfield.tev[end][i].s = soln.zero[1]
                    end
                end
                
                #check for convergence
                #Shed indices must be the same
                shedv = Int[]
                for i = 1:nsurf
                    if (abs(surf[i].a0[1])>surf[i].lespcrit[1] || i in shed_ind[end])
                        push!(shedv, i)
                    end
                end
                nshed = length(shedv)
                
                # TEV and LEV iterations should be converged 
                
                if iter > 1 
                    teverr = mean(abs(map(q->q.s, curfield.tev[end]) - tev_prev))
                    leverr = mean(abs(map(q->q.s, curfield.lev[end]) - lev_prev))
                    if (teverr < eps && leverr < eps)
                        #TEV conditions alone actually implies the others and would be sufficient
                        break
                    end
                    
                    if (iter > 10)
                        println("2D convergence failed")
                        break
                    end
                    
                end
            end  
            
            #Set the levflag parameter after the iteration is over
            for i = 1:nsurf
                if i in shed_ind[end]
                    surf[i].levflag[1] = 1
                else
                    surf[i].levflag[1] = 0
                end
            end
        end

#Update rest of Fourier terms
update_a2toan(surf)

for i = 1:nsurf
    surf[i].a0prev[1] = surf[i].a0[1]
    for ia = 1:3
        surf[i].aprev[ia] = surf[i].aterm[ia]
    end
end

#Calculate bound vortex strengths
update_bv(surf)

wakeroll(surf, curfield, dt, shed_ind)

cl, cd, cm = calc_forces(surf)

for i = 1:nsurf
    mat[istep,1,i] = t
    mat[istep,2,i] = surf[i].kinem.alpha
    mat[istep,3,i] = surf[i].kinem.h
    mat[istep,4,i] = surf[i].kinem.u
    mat[istep,5,i] = surf[i].a0[1]
    mat[istep,6,i] = cl[i]
    mat[istep,7,i] = cd[i]
    mat[istep,8,i] = cm[i]
end
end

mat, surf, curfield
end

    
function ldvm_more(surf::TwoDSurf, curfield::TwoDFlowField, nsteps::Int64 = 500; dtstar::Float64 = 0.015, delvort = DelVortDef(0, 0, 0.), mat = Array(Float64, 0, 11), kelv_enf = 0.)

    if (size(mat,1) > 0)
        t = mat[end,1]
    else
        t = 0.
    end

    #mat = zeros(nsteps,11)
    mat = mat'

    dt = dtstar*surf.c/surf.uref
    #t = 0.

    #Intialise flowfield
    for istep = 1:nsteps
        #Udpate current time
        t = t + dt

        #Update kinematic parameters
        update_kinem(surf, t)

        #Update flow field parameters if any
        update_externalvel(curfield, t)

        #Update bound vortex positions
        update_boundpos(surf, dt)

        #Add a TEV with dummy strength
        place_tev(surf,curfield,dt)

        kelv = KelvinCondition(surf,curfield)
        #Solve for TEV strength to satisfy Kelvin condition
        #curfield.tev[length(curfield.tev)].s = secant_method(kelv, 0., -0.01)
        soln = nlsolve(not_in_place(kelv), [-0.01])
        curfield.tev[length(curfield.tev)].s = soln.zero[1]

        #Check for LESP condition
        #Update values with converged value of shed tev
        #Update incduced velocities on airfoil
        update_indbound(kelv.surf, kelv.field)

        #Calculate downwash
        update_downwash(kelv.surf, [curfield.u[1],curfield.w[1]])

        #Calculate first two fourier coefficients
        update_a0anda1(kelv.surf)

        lesp = surf.a0[1]

        #Update adot
        update_a2a3adot(surf,dt)

        #2D iteration if LESP_crit is exceeded
        if (abs(lesp)>surf.lespcrit[1])
            #Remove the previous tev
            pop!(curfield.tev)
            #Add a TEV with dummy strength
            place_tev(surf,curfield,dt)

            #Add a LEV with dummy strength
            place_lev(surf,curfield,dt)

            kelvkutta = KelvinKutta(surf,curfield)
            #Solve for TEV and LEV strengths to satisfy Kelvin condition and Kutta condition at leading edge

            soln = nlsolve(not_in_place(kelvkutta), [-0.01; 0.01])
            (curfield.tev[length(curfield.tev)].s, curfield.lev[length(curfield.lev)].s) = soln.zero[1], soln.zero[2]

            surf.levflag[1] = 1
        else
            surf.levflag[1] = 0
        end


        #Update rest of Fourier terms
        update_a2toan(surf)

        #Set previous values of aterm to be used for derivatives in next time step
        surf.a0prev[1] = surf.a0[1]
        for ia = 1:3
            surf.aprev[ia] = surf.aterm[ia]
        end

        #Calculate bound vortex strengths
        update_bv(surf)

        #Remove vortices that are far away from airfoil
        if (delvort.flag == 1)
            if length(curfield.tev) > delvort.limit
                if (sqrt((curfield.tev[1].x- surf.bnd_x[div(surf.ndiv,2)])^2 + (curfield.tev[1].z- surf.bnd_z[div(surf.ndiv,2)])^2) > delvort.dist*surf.c)
                    kelv_enf = kelv_enf + curfield.tev[1].s
                    for i = 1:length(curfield.tev)-1
                        curfield.tev[i] = curfield.tev[i+1]
                    end
                    pop!(curfield.tev)
                end
            end
            if length(curfield.lev) > delvort.limit
                if (sqrt((curfield.lev[1].x- surf.bnd_x[div(surf.ndiv,2)])^2 + (curfield.lev[1].z- surf.bnd_z[div(surf.ndiv,2)])^2) > delvort.dist*surf.c)
                    kelv_enf = kelv_enf + curfield.lev[1].s
                    for i = 1:length(curfield.lev)-1
                        curfield.lev[i] = curfield.lev[i+1]
                    end
                    pop!(curfield.lev)
                end
            end
        end
        wakeroll(surf, curfield, dt)

        #cl, cd, cm, cn, cs = calc_forces(surf)
        cl, cd, cm, bc, cn, cs = calc_forces_more(surf)
        #bnd_circ = (surf.a0[1] + surf.aterm[1]/2.)
        #mat[istep,:] = [t surf.kinem.alpha surf.kinem.h surf.kinem.u surf.a0[1] cl cd cm bnd_circ cn cs]
        mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, surf.a0[1], cl, cd, cm, bc, cn, cs])
    end

    mat = mat'
    mat, surf, curfield, kelv_enf
end

function ldvm_E(surf::TwoDSurf, curfield::TwoDFlowField, nsteps::Int64 = 500, dtstar::Float64 = 0.015, delvort = DelVortDef(0, 0, 0.), mat = Array(Float64, 0, 8), kelv_enf = 0.)

    if (size(mat,1) > 0)
        t = mat[end,1]
    else
        t = 0.
    end

    #mat = zeros(nsteps,11)
    mat = mat'

    dt = dtstar*surf.c/surf.uref
    #t = 0.

    #Intialise flowfield
    for istep = 1:nsteps
        #Udpate current time
        t = t + dt

        #Update kinematic parameters
        update_kinem(surf, t)

        #Update flow field parameters if any
        update_externalvel(curfield, t)

        #Update bound vortex positions
        update_boundpos(surf, dt)

        #Add a TEV with dummy strength
        place_tev(surf,curfield,dt)

        kelv = KelvinCondition(surf,curfield)
        #Solve for TEV strength to satisfy Kelvin condition
        #curfield.tev[length(curfield.tev)].s = secant_method(kelv, 0., -0.01)
        soln = nlsolve(not_in_place(kelv), [-0.01])
        curfield.tev[length(curfield.tev)].s = soln.zero[1]

        #Check for LESP condition
        #Update values with converged value of shed tev
        #Update incduced velocities on airfoil
        update_indbound(kelv.surf, kelv.field)

        #Calculate downwash
        update_downwash(kelv.surf, [curfield.u[1],curfield.w[1]])

        #Calculate first two fourier coefficients
        update_a0anda1(kelv.surf)

        lesp = surf.a0[1]

        #Update adot - No need for this, done after 2D NR
        #update_a2a3adot(surf,dt)

        #2D iteration if LESP_crit is exceeded
        if (abs(lesp)>surf.lespcrit[1])
            #Remove the previous tev
            pop!(curfield.tev)
            #Add a TEV with dummy strength
            place_tev(surf,curfield,dt)

            #Add a LEV with dummy strength
            place_lev(surf,curfield,dt)

            kelvkutta = KelvinKutta(surf,curfield)
            #Solve for TEV and LEV strengths to satisfy Kelvin condition and Kutta condition at leading edge

            soln = nlsolve(not_in_place(kelvkutta), [-0.01; 0.01])
            (curfield.tev[length(curfield.tev)].s, curfield.lev[length(curfield.lev)].s) = soln.zero[1], soln.zero[2]

            surf.levflag[1] = 1
        else
            surf.levflag[1] = 0
        end


        #Update rest of Fourier terms
        update_a2toan(surf)

        #Update derivatives of Fourier coefficients
        update_adot(surf,dt)

        #Set previous values of aterm to be used for derivatives in next time step
        surf.a0prev[1] = surf.a0[1]
        for ia = 1:3
            surf.aprev[ia] = surf.aterm[ia]
        end

        #Calculate bound vortex strengths
        update_bv(surf)

        #Remove vortices that are far away from airfoil
        if (delvort.flag == 1)
            if length(curfield.tev) > delvort.limit
                if (sqrt((curfield.tev[1].x- surf.bnd_x[div(surf.ndiv,2)])^2 + (curfield.tev[1].z- surf.bnd_z[div(surf.ndiv,2)])^2) > delvort.dist*surf.c)
                    kelv_enf = kelv_enf + curfield.tev[1].s
                    for i = 1:length(curfield.tev)-1
                        curfield.tev[i] = curfield.tev[i+1]
                    end
                    pop!(curfield.tev)
                end
            end
            if length(curfield.lev) > delvort.limit
                if (sqrt((curfield.lev[1].x- surf.bnd_x[div(surf.ndiv,2)])^2 + (curfield.lev[1].z- surf.bnd_z[div(surf.ndiv,2)])^2) > delvort.dist*surf.c)
                    kelv_enf = kelv_enf + curfield.lev[1].s
                    for i = 1:length(curfield.lev)-1
                        curfield.lev[i] = curfield.lev[i+1]
                    end
                    pop!(curfield.lev)
                end
            end
        end
        wakeroll(surf, curfield, dt)

        if (surf.levflag[1] == 1)
            cl, cd, cm = calc_forces_E(surf,curfield.lev[length(curfield.lev)].s, dt)
        else
            cl, cd, cm = calc_forces(surf)
        end

        #bnd_circ = (surf.a0[1] + surf.aterm[1]/2.)
        #mat[istep,:] = [t surf.kinem.alpha surf.kinem.h surf.kinem.u surf.a0[1] cl cd cm bnd_circ cn cs]
        mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, surf.a0[1], cl, cd, cm])
    end

    mat = mat'
    mat, surf, curfield, kelv_enf
end

function ldvm_E_more(surf::TwoDSurf, curfield::TwoDFlowField, nsteps::Int64 = 500, dtstar::Float64 = 0.015, delvort = DelVortDef(0, 0, 0.), mat = Array(Float64, 0, 11), kelv_enf = 0.)

    if (size(mat,1) > 0)
        t = mat[end,1]
    else
        t = 0.
    end

    #mat = zeros(nsteps,11)
    mat = mat'

    dt = dtstar*surf.c/surf.uref
    #t = 0.

    #Intialise flowfield
    for istep = 1:nsteps
        #Udpate current time
        t = t + dt

        #Update kinematic parameters
        update_kinem(surf, t)

        #Update flow field parameters if any
        update_externalvel(curfield, t)

        #Update bound vortex positions
        update_boundpos(surf, dt)

        #Add a TEV with dummy strength
        place_tev(surf,curfield,dt)

        kelv = KelvinCondition(surf,curfield)
        #Solve for TEV strength to satisfy Kelvin condition
        #curfield.tev[length(curfield.tev)].s = secant_method(kelv, 0., -0.01)
        soln = nlsolve(not_in_place(kelv), [-0.01])
        curfield.tev[length(curfield.tev)].s = soln.zero[1]

        #Check for LESP condition
        #Update values with converged value of shed tev
        #Update incduced velocities on airfoil
        update_indbound(kelv.surf, kelv.field)

        #Calculate downwash
        update_downwash(kelv.surf, [curfield.u[1],curfield.w[1]])

        #Calculate first two fourier coefficients
        update_a0anda1(kelv.surf)

        lesp = surf.a0[1]

        #Update adot - No need for this, done after 2D NR
        #update_a2a3adot(surf,dt)

        #2D iteration if LESP_crit is exceeded
        if (abs(lesp)>surf.lespcrit[1])
            #Remove the previous tev
            pop!(curfield.tev)
            #Add a TEV with dummy strength
            place_tev(surf,curfield,dt)

            #Add a LEV with dummy strength
            place_lev(surf,curfield,dt)

            kelvkutta = KelvinKutta(surf,curfield)
            #Solve for TEV and LEV strengths to satisfy Kelvin condition and Kutta condition at leading edge

            soln = nlsolve(not_in_place(kelvkutta), [-0.01; 0.01])
            (curfield.tev[length(curfield.tev)].s, curfield.lev[length(curfield.lev)].s) = soln.zero[1], soln.zero[2]

            surf.levflag[1] = 1
        else
            surf.levflag[1] = 0
        end


        #Update rest of Fourier terms
        update_a2toan(surf)

        #Update derivatives of Fourier coefficients
        update_adot(surf,dt)

        #Set previous values of aterm to be used for derivatives in next time step
        surf.a0prev[1] = surf.a0[1]
        for ia = 1:3
            surf.aprev[ia] = surf.aterm[ia]
        end

        #Calculate bound vortex strengths
        update_bv(surf)

        #Remove vortices that are far away from airfoil
        if (delvort.flag == 1)
            if length(curfield.tev) > delvort.limit
                if (sqrt((curfield.tev[1].x- surf.bnd_x[div(surf.ndiv,2)])^2 + (curfield.tev[1].z- surf.bnd_z[div(surf.ndiv,2)])^2) > delvort.dist*surf.c)
                    kelv_enf = kelv_enf + curfield.tev[1].s
                    for i = 1:length(curfield.tev)-1
                        curfield.tev[i] = curfield.tev[i+1]
                    end
                    pop!(curfield.tev)
                end
            end
            if length(curfield.lev) > delvort.limit
                if (sqrt((curfield.lev[1].x- surf.bnd_x[div(surf.ndiv,2)])^2 + (curfield.lev[1].z- surf.bnd_z[div(surf.ndiv,2)])^2) > delvort.dist*surf.c)
                    kelv_enf = kelv_enf + curfield.lev[1].s
                    for i = 1:length(curfield.lev)-1
                        curfield.lev[i] = curfield.lev[i+1]
                    end
                    pop!(curfield.lev)
                end
            end
        end
        wakeroll(surf, curfield, dt)

        if (surf.levflag[1] == 1)
            cl, cd, cm, bc, cn, cs = calc_forces_E_more(surf,curfield.lev[length(curfield.lev)].s, dt)
        else
            cl, cd, cm, bc, cn, cs = calc_forces_more(surf)
        end

        #bnd_circ = (surf.a0[1] + surf.aterm[1]/2.)
        #mat[istep,:] = [t surf.kinem.alpha surf.kinem.h surf.kinem.u surf.a0[1] cl cd cm bnd_circ cn cs]
        mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, surf.a0[1], cl, cd, cm, bc, cn, cs])
    end

    mat = mat'
    mat, surf, curfield, kelv_enf
end


function ldvm(surf::TwoDSurfwFlap, curfield::TwoDFlowField, nsteps::Int64 = 500, dtstar::Float64 = 0.015, delvort = DelVort(0, 0, 0.), mat = Array(Float64, 0, 9), kelv_enf = 0.)

    if (size(mat,1) > 0)
        t = mat[end,1]
    else
        t = 0.
    end

    mat = mat'

    dt = dtstar*surf.c/surf.uref

    #Intialise flowfield
    for istep = 1:nsteps
        #Udpate current time
        t = t + dt

        #Update kinematic parameters
        update_kinem(surf, t)

        #Update deformation
        update_deform(surf, t)

        #Update bound vortex positions
        update_boundpos(surf, dt)

        #Add a TEV with dummy strength
        place_tev(surf,curfield,dt)

        kelv = KelvinConditionwFlap(surf,curfield)
        #Solve for TEV strength to satisfy Kelvin condition
        #curfield.tev[length(curfield.tev)].s = secant_method(kelv, 0., -0.01)
        soln = nlsolve(not_in_place(kelv), [-0.01])
        curfield.tev[length(curfield.tev)].s = soln.zero[1]

        #Check for LESP condition
        #Update values with converged value of shed tev
        #Update incduced velocities on airfoil
        update_indbound(kelv.surf, kelv.field)

        #Calculate downwash
        update_downwash(kelv.surf)

        #Calculate first two fourier coefficients
        update_a0anda1(kelv.surf)

        lesp = surf.a0[1]

        #Update adot
        update_a2a3adot(surf,dt)

        #2D iteration if LESP_crit is exceeded
        if (abs(lesp)>surf.lespcrit[1])
            #Remove the previous tev
            pop!(curfield.tev)
            #Add a TEV with dummy strength
            place_tev(surf,curfield,dt)

            #Add a LEV with dummy strength
            place_lev(surf,curfield,dt)

            kelvkutta = KelvinKuttawFlap(surf,curfield)
            #Solve for TEV and LEV strengths to satisfy Kelvin condition and Kutta condition at leading edge

            soln = nlsolve(not_in_place(kelvkutta), [-0.01; 0.01])
            (curfield.tev[length(curfield.tev)].s, curfield.lev[length(curfield.lev)].s) = soln.zero[1], soln.zero[2]

            surf.levflag[1] = 1
        else
            surf.levflag[1] = 0
        end


        #Update rest of Fourier terms
        update_a2toan(surf)

        #Set previous values of aterm to be used for derivatives in next time step
        surf.a0prev[1] = surf.a0[1]
        for ia = 1:3
            surf.aprev[ia] = surf.aterm[ia]
        end

        #Calculate bound vortex strengths
        update_bv(surf)

        #Remove vortices that are far away from airfoil
        if (delvort.flag == 1)
            if length(curfield.tev) > delvort.limit
                if (sqrt((curfield.tev[1].x- surf.bnd_x[div(surf.ndiv,2)])^2 + (curfield.tev[1].z- surf.bnd_z[div(surf.ndiv,2)])^2) > delvort.dist*surf.c)
                    kelv_enf = kelv_enf + curfield.tev[1].s
                    for i = 1:length(curfield.tev)-1
                        curfield.tev[i] = curfield.tev[i+1]
                    end
                    pop!(curfield.tev)
                end
            end
            if length(curfield.lev) > delvort.limit
                if (sqrt((curfield.lev[1].x- surf.bnd_x[div(surf.ndiv,2)])^2 + (curfield.lev[1].z- surf.bnd_z[div(surf.ndiv,2)])^2) > delvort.dist*surf.c)
                    kelv_enf = kelv_enf + curfield.lev[1].s
                    for i = 1:length(curfield.lev)-1
                        curfield.lev[i] = curfield.lev[i+1]
                    end
                    pop!(curfield.lev)
                end
            end
        end

        wakeroll(surf, curfield, dt)

        cl, cd, cm, cm_be = calc_forces(surf, dt)

        mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, surf.a0[1], cl, cd, cm, cm_be])

    end

    mat = mat'
    mat, surf, curfield, kelv_enf
    #Plot flowfield viz
#    figure(0)
#    view_vorts(surf, curfield)

end

function ldvm(surf::TwoDSurf_2DOF, curfield::TwoDFlowField, nsteps::Int64 = 500, dtstar::Float64 = 0.015, delvort = DelVortDef(0, 0, 0.), mat = Array(Float64, 0, 8), kelv_enf = 0.)

    if (size(mat,1) > 0)
        t = mat[end,1]
    else
        t = 0.
    end

    #mat = zeros(nsteps,8)
    mat = mat'

    dt = dtstar*surf.c/surf.uref
    #t = 0.
    #kelv_enf = 0

    cl = 0.
    cm = 0.
    
    #Intialise flowfield
    for istep = 1:nsteps
        #Udpate current time
        t = t + dt

        #Update kinematic parameters (based on 2DOF response)
        #if (t > dt) # Allow initial condition
        update_kinem(surf, dt, cl, cm)
        #end
        #Update bound vortex positions
        update_boundpos(surf, dt)

        #Add a TEV with dummy strength
        place_tev(surf,curfield,dt)

        kelv = KelvinCondition2DOF(surf,curfield,kelv_enf)
        #Solve for TEV strength to satisfy Kelvin condition
        #curfield.tev[length(curfield.tev)].s = secant_method(kelv, 0., -0.01)
        soln = nlsolve(not_in_place(kelv), [-0.01])
        curfield.tev[length(curfield.tev)].s = soln.zero[1]

        #Check for LESP condition
        #Update values with converged value of shed tev
        #Update incduced velocities on airfoil
        update_indbound(kelv.surf, kelv.field)

        #Calculate downwash
        update_downwash(kelv.surf)

        #Calculate first two fourier coefficients
        update_a0anda1(kelv.surf)

        lesp = surf.a0[1]

        #Update adot
        update_a2a3adot(surf,dt)

        #2D iteration if LESP_crit is exceeded
        if (abs(lesp)>surf.lespcrit[1])
            #Remove the previous tev
            pop!(curfield.tev)
            #Add a TEV with dummy strength
            place_tev(surf,curfield,dt)

            #Add a LEV with dummy strength
            place_lev(surf,curfield,dt)

            kelvkutta = KelvinKutta2DOF(surf,curfield,kelv_enf)
            #Solve for TEV and LEV strengths to satisfy Kelvin condition and Kutta condition at leading edge

            soln = nlsolve(not_in_place(kelvkutta), [-0.01; 0.01])
            (curfield.tev[length(curfield.tev)].s, curfield.lev[length(curfield.lev)].s) = soln.zero[1], soln.zero[2]

            surf.levflag[1] = 1
        else
            surf.levflag[1] = 0
        end


        #Update rest of Fourier terms
        update_a2toan(surf)

        #Set previous values of aterm to be used for derivatives in next time step
        surf.a0prev[1] = surf.a0[1]
        for ia = 1:3
            surf.aprev[ia] = surf.aterm[ia]
        end

        #Calculate bound vortex strengths
        update_bv(surf)

        wakeroll(surf, curfield, dt)

        #Remove vortices that are far away from airfoil
        if (delvort.flag == 1)
            if length(curfield.tev) > delvort.limit
                if (sqrt((curfield.tev[1].x- surf.bnd_x[div(surf.ndiv,2)])^2 + (curfield.tev[1].z- surf.bnd_z[div(surf.ndiv,2)])^2) > delvort.dist*surf.c)
                    kelv_enf = kelv_enf + curfield.tev[1].s
                    shift!(curfield.tev)
                end
            end
            if length(curfield.lev) > delvort.limit
                if (sqrt((curfield.lev[1].x- surf.bnd_x[div(surf.ndiv,2)])^2 + (curfield.lev[1].z- surf.bnd_z[div(surf.ndiv,2)])^2) > delvort.dist*surf.c)
                    kelv_enf = kelv_enf + curfield.lev[1].s
                    shift!(curfield.lev)
                end
            end
        end

        

        
        #Calculate forces
        cl, cd, cm = calc_forces(surf)

        #Update other kinematic terms in KinemPar2DOF
        #update_kinem2DOF(surf)
        
        #Using the force data, update - hddot and alphaddot
        #calc_struct2DOF(surf, cl, cm)

        #Write out nondimensional quantities
        mat = hcat(mat,[t*surf.uref/surf.c, surf.kinem.alpha, surf.kinem.h/surf.c, surf.kinem.u, surf.a0[1], cl, cd, cm])
        #mat[istep,:] = [t surf.kinem.alpha surf.kinem.h surf.kinem.u surf.a0[1] cl cd cm]
    end
    mat = mat'
    mat, surf, curfield, kelv_enf
    #Plot flowfield viz
#    figure(0)
#    view_vorts(surf, curfield)

end

function ldvm_klb(surf::TwoDSurf_2DOF, curfield::TwoDFlowField, nsteps::Int64 = 500, dtstar::Float64 = 0.015, delvort = DelVortDef(0, 0, 0.), mat = Array(Float64, 0, 10), kelv_enf = 0., tf = 3.0, alpha1 = 7*pi/180, S1 = 3.0*pi/180, S2 = 2.3*pi/180)
    
    if (size(mat,1) > 0)
        t = mat[end,1]
    else
        t = 0.
    end

    #mat = zeros(nsteps,8)
    mat = mat'

    dt = dtstar*surf.c/surf.uref

    cl = 0.
    cm = 0.

    dfn = 0.
    dfn_pr = 0.
    fsep = 1.
    f0sep = 1.
    f0sep_pr = 1.
    
    #Intialise flowfield
    for istep = 1:nsteps
        #Udpate current time
        t = t + dt

        #Update kinematic parameters (based on 2DOF response)
        #if (t > dt) # Allow initial condition
        update_kinem(surf, dt, cl, cm)
        #end
        #Update bound vortex positions
        update_boundpos(surf, dt)

        #Add a TEV with dummy strength
        place_tev(surf,curfield,dt)

        kelv = KelvinCondition2DOF(surf,curfield,kelv_enf)
        #Solve for TEV strength to satisfy Kelvin condition
        #curfield.tev[length(curfield.tev)].s = secant_method(kelv, 0., -0.01)
        soln = nlsolve(not_in_place(kelv), [-0.01])
        curfield.tev[length(curfield.tev)].s = soln.zero[1]

        #Check for LESP condition
        #Update values with converged value of shed tev
        #Update incduced velocities on airfoil
        update_indbound(kelv.surf, kelv.field)

        #Calculate downwash
        update_downwash(kelv.surf)

        #Calculate first two fourier coefficients
        update_a0anda1(kelv.surf)

        #Calculate fsep
        f0sep_pr = f0sep
        if abs(surf.a0[1]) < alpha1
            f0sep = 1 - 0.3*exp((abs(surf.a0[1]) - alpha1)/S1)
        else
            f0sep = 0.04 + 0.66*exp((alpha1 - abs(surf.a0[1]))/S2)
        end

        dfn_pr = dfn
        dfn = dfn_pr*exp(-2*dt/tf)+(f0sep - f0sep_pr)*exp(-dt/tf)
        fsep = f0sep - dfn
        
        lesp = surf.a0[1]

        #Update adot
        update_a2a3adot(surf,dt)

        #2D iteration if LESP_crit is exceeded
        if (abs(lesp)>surf.lespcrit[1])
            #Remove the previous tev
            pop!(curfield.tev)
            #Add a TEV with dummy strength
            place_tev(surf,curfield,dt)

            #Add a LEV with dummy strength
            place_lev(surf,curfield,dt)

            kelvkutta = KelvinKutta2DOF(surf,curfield,kelv_enf)
            #Solve for TEV and LEV strengths to satisfy Kelvin condition and Kutta condition at leading edge

            soln = nlsolve(not_in_place(kelvkutta), [-0.01; 0.01])
            (curfield.tev[length(curfield.tev)].s, curfield.lev[length(curfield.lev)].s) = soln.zero[1], soln.zero[2]

            surf.levflag[1] = 1
        else
            surf.levflag[1] = 0
        end


        #Update rest of Fourier terms
        update_a2toan(surf)

        #Set previous values of aterm to be used for derivatives in next time step
        surf.a0prev[1] = surf.a0[1]
        for ia = 1:3
            surf.aprev[ia] = surf.aterm[ia]
        end

        #Calculate bound vortex strengths
        update_bv(surf)

        wakeroll(surf, curfield, dt)

        #Remove vortices that are far away from airfoil
        if (delvort.flag == 1)
            if length(curfield.tev) > delvort.limit
                if (sqrt((curfield.tev[1].x- surf.bnd_x[div(surf.ndiv,2)])^2 + (curfield.tev[1].z- surf.bnd_z[div(surf.ndiv,2)])^2) > delvort.dist*surf.c)
                    kelv_enf = kelv_enf + curfield.tev[1].s
                    shift!(curfield.tev)
                end
            end
            if length(curfield.lev) > delvort.limit
                if (sqrt((curfield.lev[1].x- surf.bnd_x[div(surf.ndiv,2)])^2 + (curfield.lev[1].z- surf.bnd_z[div(surf.ndiv,2)])^2) > delvort.dist*surf.c)
                    kelv_enf = kelv_enf + curfield.lev[1].s
                    shift!(curfield.lev)
                end
            end
        end

        

        
        #Calculate forces
        cl, cd, cm = calc_forces(surf, fsep)

        #Update other kinematic terms in KinemPar2DOF
        #update_kinem2DOF(surf)
        
        #Using the force data, update - hddot and alphaddot
        #calc_struct2DOF(surf, cl, cm)

        #Write out nondimensional quantities
        mat = hcat(mat,[t*surf.uref/surf.c, surf.kinem.alpha, surf.kinem.h/surf.c, surf.kinem.u, surf.a0[1], cl, cd, cm, f0sep, fsep])
        #mat[istep,:] = [t surf.kinem.alpha surf.kinem.h surf.kinem.u surf.a0[1] cl cd cm]
    end
    mat = mat'
    mat, surf, curfield, kelv_enf
    #Plot flowfield viz
    #    figure(0)
    #    view_vorts(surf, curfield)

end

function ldvm_E(surf::TwoDSurf_2DOF, curfield::TwoDFlowField, nsteps::Int64 = 500, dtstar::Float64 = 0.015, delvort = DelVortDef(0, 0, 0.), mat = Array(Float64, 0, 8), kelv_enf = 0.)

    if (size(mat,1) > 0)
        t = mat[end,1]
    else
        t = 0.
    end

    #mat = zeros(nsteps,8)
    mat = mat'

    dt = dtstar*surf.c/surf.uref
    #t = 0.
    #kelv_enf = 0

    #Intialise flowfield
    for istep = 1:nsteps
        #Udpate current time
        t = t + dt

        #Update kinematic parameters (based on 2DOF response)
        if (t > dt) # Allow initial condition
            update_kinem(surf, dt)
        end
        #Update bound vortex positions
        update_boundpos(surf, dt)

        #Add a TEV with dummy strength
        place_tev(surf,curfield,dt)

        kelv = KelvinCondition2DOF(surf,curfield,kelv_enf)
        #Solve for TEV strength to satisfy Kelvin condition
        #curfield.tev[length(curfield.tev)].s = secant_method(kelv, 0., -0.01)
        soln = nlsolve(not_in_place(kelv), [-0.01])
        curfield.tev[length(curfield.tev)].s = soln.zero[1]

        #Check for LESP condition
        #Update values with converged value of shed tev
        #Update incduced velocities on airfoil
        update_indbound(kelv.surf, kelv.field)

        #Calculate downwash
        update_downwash(kelv.surf)

        #Calculate first two fourier coefficients
        update_a0anda1(kelv.surf)

        lesp = surf.a0[1]


        #2D iteration if LESP_crit is exceeded
        if (abs(lesp)>surf.lespcrit[1])
            #Remove the previous tev
            pop!(curfield.tev)
            #Add a TEV with dummy strength
            place_tev(surf,curfield,dt)

            #Add a LEV with dummy strength
            place_lev(surf,curfield,dt)

            kelvkutta = KelvinKutta2DOF(surf,curfield,kelv_enf)
            #Solve for TEV and LEV strengths to satisfy Kelvin condition and Kutta condition at leading edge

            soln = nlsolve(not_in_place(kelvkutta), [-0.01; 0.01])
            (curfield.tev[length(curfield.tev)].s, curfield.lev[length(curfield.lev)].s) = soln.zero[1], soln.zero[2]

            surf.levflag[1] = 1
        else
            surf.levflag[1] = 0
        end


        #Update rest of Fourier terms
        update_a2toan(surf)

        #Update derivatives of Fourier coefficients
        update_adot(surf,dt)

        #Set previous values of aterm to be used for derivatives in next time step
        surf.a0prev[1] = surf.a0[1]
        for ia = 1:3
            surf.aprev[ia] = surf.aterm[ia]
        end

        #Calculate bound vortex strengths
        update_bv(surf)

        #Remove vortices that are far away from airfoil
        if (delvort.flag == 1)
            if length(curfield.tev) > delvort.limit
                if (sqrt((curfield.tev[1].x- surf.bnd_x[div(surf.ndiv,2)])^2 + (curfield.tev[1].z- surf.bnd_z[div(surf.ndiv,2)])^2) > delvort.dist*surf.c)
                    kelv_enf = kelv_enf + curfield.tev[1].s
                    for i = 1:length(curfield.tev)-1
                        curfield.tev[i] = curfield.tev[i+1]
                    end
                    pop!(curfield.tev)
                end
            end
            if length(curfield.lev) > delvort.limit
                if (sqrt((curfield.lev[1].x- surf.bnd_x[div(surf.ndiv,2)])^2 + (curfield.lev[1].z- surf.bnd_z[div(surf.ndiv,2)])^2) > delvort.dist*surf.c)
                    kelv_enf = kelv_enf + curfield.lev[1].s
                    for i = 1:length(curfield.lev)-1
                        curfield.lev[i] = curfield.lev[i+1]
                    end
                    pop!(curfield.lev)
                end
            end
        end

        wakeroll(surf, curfield, dt)

        #Update kinematic terms in KinemPar2DOF
        update_kinem2DOF(surf)

        #Calculate forces
        if (surf.levflag[1] == 1)
            cl, cd, cm = calc_forces_E(surf,curfield.lev[length(curfield.lev)].s, dt)
        else
            cl, cd, cm = calc_forces(surf)
        end

        #Using the force data, update - hddot and alphaddot
        calc_struct2DOF(surf, cl, cm)
        mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, surf.a0[1], cl, cd, cm])
        #mat[istep,:] = [t surf.kinem.alpha surf.kinem.h surf.kinem.u surf.a0[1] cl cd cm]
    end
    mat = mat'
    mat, surf, curfield, kelv_enf
    #Plot flowfield viz
#    figure(0)
#    view_vorts(surf, curfield)

end

function drone_trajectory_problem(surf::TwoDFreeSurf, curfield::TwoDFlowField, nsteps::Int64 = 500, dtstar::Float64 = 0.015, cf::Float64 = 0)
    mat = zeros(nsteps,10)
    ind_fr = 0
    fr_freq = 1000

    dtstar = 0.015
    dt = dtstar*surf.c/surf.uref
    t = 0.
    kelv_enf = 0

    #Intialise flowfield
    for istep = 1:nsteps

        #Dynamically determine dt
        #If the shooting vorrtex is in the vicinity of airfoil, use a smaller time step
        if (sqrt((curfield.extv[1].x-surf.bv[35].x)^2 + (curfield.extv[1].z-surf.bv[35].z)^2) < 10*surf.c)
            dt = minimum([(0.015*0.2*3)/maximum([abs(surf.kinem.alphadot) abs(surf.kinem.hdot) abs(surf.kinem.u) abs(curfield.extv[1].vx)/2.]) 0.015]) #this is dimensional value - c/u from K cancels with dtstar
        else
            dt = minimum([(0.015*0.2*3)/maximum([abs(surf.kinem.alphadot) abs(surf.kinem.hdot) abs(surf.kinem.u)]) 0.015]) #this is dimensional value - c/u from K cancels with dtstar
        end

        #Udpate current time
        t = t + dt

        #Update kinematic parameters (based on 2DOF response)
        if (istep > 1) # Allow initial condition
            update_kinem(surf, dt)
        end
        println(istep," ", t," ",length(curfield.tev)," ", length(curfield.lev))

        #Check if ground is breached (for drone crash problem)
        if (surf.kinem.h < 0)
            mat, surf, curfield
            break
        end
        #Update bound vortex positions
        update_boundpos(surf, dt)

        #Add a TEV with dummy strength
        place_tev(surf,curfield,dt)

        kelv = KelvinCondition2DFree(surf,curfield, kelv_enf)
        #Solve for TEV strength to satisfy Kelvin condition
        #curfield.tev[length(curfield.tev)].s = secant_method(kelv, 0., -0.01)
        soln = nlsolve(not_in_place(kelv), [-0.01])
        curfield.tev[length(curfield.tev)].s = soln.zero[1]

        #Check for LESP condition
        #Update values with converged value of shed tev
        #Update incduced velocities on airfoil
        update_indbound(kelv.surf, kelv.field)

        #Calculate downwash
        update_downwash(kelv.surf)

        #Calculate first two fourier coefficients
        update_a0anda1(kelv.surf)

        if surf.kinem.u == 0
            lesp = surf.a0[1]
        else
            lesp = surf.a0[1]/surf.kinem.u
        end

        #Update adot
        update_a2a3adot(surf,dt)


        #2D iteration if LESP_crit is exceeded
        if (abs(lesp)>surf.lespcrit[1])
            #Remove the previous tev
            pop!(curfield.tev)
            #Add a TEV with dummy strength
            place_tev(surf,curfield,dt)

            #Add a LEV with dummy strength
            place_lev(surf,curfield,dt)

            kelvkutta = KelvinKutta2DFree(surf,curfield, kelv_enf)
            #Solve for TEV and LEV strengths to satisfy Kelvin condition and Kutta condition at leading edge

            soln = nlsolve(not_in_place(kelvkutta), [-0.01; 0.01])
            (curfield.tev[length(curfield.tev)].s, curfield.lev[length(curfield.lev)].s) = soln.zero[1], soln.zero[2]

            surf.levflag[1] = 1
        else
            surf.levflag[1] = 0
        end


        #Update rest of Fourier terms
        update_a2toan(surf)

        #Set previous values of aterm to be used for derivatives in next time step
        surf.a0prev[1] = surf.a0[1]
        for ia = 1:3
            surf.aprev[ia] = surf.aterm[ia]
        end

        #Remove vortices that are far away from airfoil
        if length(curfield.tev) > 100
            if (sqrt((curfield.tev[1].x- surf.bnd_x[div(surf.ndiv,2)])^2 + (curfield.tev[1].z- surf.bnd_z[div(surf.ndiv,2)])^2) > 10*surf.c)
                kelv_enf = kelv_enf + curfield.tev[1].s
                for i = 1:length(curfield.tev)-1
                    curfield.tev[i] = curfield.tev[i+1]
                end
                pop!(curfield.tev)
            end
        end
        if length(curfield.lev) > 100
            if (sqrt((curfield.lev[1].x- surf.bnd_x[div(surf.ndiv,2)])^2 + (curfield.lev[1].z- surf.bnd_z[div(surf.ndiv,2)])^2) > 10*surf.c)
                kelv_enf = kelv_enf + curfield.lev[1].s
                for i = 1:length(curfield.lev)-1
                    curfield.lev[i] = curfield.lev[i+1]
                end
                pop!(curfield.lev)
            end
        end

        #Calculate bound vortex strengths
        update_bv(surf)

        wakeroll(surf, curfield, dt)

        #Update kinematic terms in KinemPar2DOF
        update_kinem2DOF(surf)

        #Calculate forces
        cl, cd, cm = calc_forces(surf)

        #Using the force data, update - hddot and alphaddot
        calc_moveFree(surf, cl, cd, cm, cf)

        mat[istep,:] = [t surf.kinem.alpha surf.kinem.h surf.kinem.u surf.a0[1] cl cd cm surf.bnd_x[div(surf.ndiv,2)] surf.bnd_z[div(surf.ndiv,2)]]
        if (rem(istep,fr_freq) == 0)
            ind_fr += 1
            outfile = open("dump/field.$istep", "w")
            flowmat = zeros(1+length(curfield.tev)+length(curfield.lev)+length(curfield.extv)+length(surf.bv),4)
            flowmat[1,:] = [length(curfield.tev), length(curfield.lev), length(curfield.extv), length(surf.bv)]
            flowmat[2:end,1] = [map(q->q.s, curfield.tev); map(q->q.s, curfield.lev); map(q->q.s, curfield.extv); map(q->q.s, surf.bv)]
            flowmat[2:end,2] = [map(q->q.x, curfield.tev); map(q->q.x, curfield.lev); map(q->q.x, curfield.extv); map(q->q.x, surf.bv)]
            flowmat[2:end,3] = [map(q->q.z, curfield.tev); map(q->q.z, curfield.lev); map(q->q.z, curfield.extv); map(q->q.z, surf.bv)]
            writedlm(outfile, flowmat)
            close(outfile)
        end
    end

    mat, surf, curfield
    #Plot flowfield viz
#    figure(0)
#    view_vorts(surf, curfield)

end

function lesp_design_max(h_amp::Float64, alphadef::MotionDef)
  hdef = EldUpIntDef(h_amp,alphadef.K*h_amp/alphadef.amp,alphadef.a)
  udef = ConstDef(1.)
  full_kinem = KinemDef(alphadef, hdef, udef)
  lespcrit = [20;]
  pvt = 0.25

  surf = TwoDSurf(1., 1., "sd7003_fine.dat", pvt, 70, 35, "Prescribed", full_kinem,lespcrit)

  curfield = TwoDFlowField()

  nsteps =round(Int,3.5/0.015)+1

  data, surf, curfield = ldvm(surf, curfield, nsteps)

  return abs(maximum(data[:,5]) - 0.25)
end

function design_solve(alphadef::MotionDef)
  iter_h = zeros(10)
  ld = zeros(10)
  iter_h[1] = 0.
  iter_h[2] = 0.1
  ld[1] = lesp_design_max(iter_h[1],alphadef)
  iter_max = 11
  iter = 1
  eps = 1e-08

  while (ld[iter] > eps)
    if (iter > iter_max)
      error("Iteration has failed")
    end
    iter = iter + 1
    ld[iter] = lesp_design_max(iter_h[iter],alphadef)
    dld = (ld[iter] - ld[iter-1])/(iter_h[iter] - iter_h[iter-1])
    iter_h[iter+1] = iter_h[iter] - ld[iter]/dld
  end
  return iter_h[iter]
end

# function QScorrect_lautat(surf :: ThreeDSurf, field :: ThreeDFlowField, nsteps :: Int64, dtstar :: Float64)

#     psi = zeros(surf.nspan)
#     dpsi = pi/surf.nspan

#     for i = 1:surf.nspan
#         psi[i] = (real(i)-0.5)*dpsi
#     end

#     mat = Array{Float64, 2}[]

#     if surf.kindef.vartype == "Constant"

#         #This is just the 2D solution corrected to 3D - kinematics at all strips is the same

#         kinem2d = KinemDef(surf.kindef.alpha, surf.kindef.h, surf.kindef.u)

#         surf2d = TwoDSurf(surf.patchdata[1].coord_file, surf.patchdata[1].pvt, kinem2d, [surf.patchdata[1].lc;])

#         #If 3D flow field is defined with disturbances or external vortices, these should be transferred to the 2D flowfield
#         curfield2d = TwoDFlowField()

#         mat2d, surf2d, curfield2d = lautat_wakeroll_more(surf2d, curfield2d, nsteps, dtstar)

#         #Fill the 2D solution into the 3D data structures

#         for i = 1:surf.nspan
#             push!(mat, mat2d)
#             surf.cam[i,:] = surf2d.cam
#             surf.cam_slope[i,:] = surf2d.cam_slope
#             surf.kinem[i,:] = surf2d.kinem
#             surf.bnd_x[i,:] = surf2d.bnd_x
#             surf.bnd_z[i,:] = surf2d.bnd_z
#             surf.uind[i,:] = surf2d.uind
#             surf.vind[i,:] = 0
#             surf.wind[i,:] = surf2d.wind
#             surf.downwash[i,:] = surf2d.downwash
#             surf.a0[i] = surf2d.a0[1]
#             surf.aterm[i,:] = surf2d.aterm
#             surf.a0dot[i] = surf2d.a0dot[1]
#             surf.adot[i,:] = surf2d.adot
#             surf.a0prev[i] = surf2d.a0prev[1]
#             surf.aprev[i,:] = surf2d.aprev
#             surf.levflag[i] = surf2d.levflag[1]

#             #Convert the bound vortices into 3d vortices and add them
#             for j = 1:length(surf2d.bv)
#                 surf.bv[i,j] = ThreeDVort([surf2d.bv[j].x; surf.yle[i]; surf2d.bv[j].z], [0.; surf2d.bv[j].s; 0.], surf2d.bv[j].vc, surf2d.bv[j].vx, 0., surf2d.bv[j].vz)
#             end
#         end

#         for j = 1:surf.nspan
#             for i = 1:length(curfield2d.lev)
#                 push!(field.lev, ThreeDVort([curfield2d.lev[i].x; surf.yle[j]; curfield2d.lev[i].z], [0.; curfield2d.lev[i].s; 0.], curfield2d.lev[i].vc, curfield2d.lev[i].vx, 0., curfield2d.lev[i].vz))
#             end
#             for i = 1:length(curfield2d.tev)
#                 push!(field.tev, ThreeDVort([curfield2d.tev[i].x; surf.yle[j]; curfield2d.tev[i].z], [0.; curfield2d.tev[i].s; 0.], curfield2d.tev[i].vc, curfield2d.tev[i].vx, 0., curfield2d.tev[i].vz))
#             end
#             for i = 1:length(curfield2d.extv)
#                 push!(field.extv, ThreeDVort([curfield2d.extv[i].x; surf.yle[j]; curfield2d.extv[i].z], [0.; curfield2d.extv[i].s; 0.], curfield2d.extv[i].vc, curfield2d.extv[i].vx, 0., curfield2d.extv[i].vz))
#             end
#         end

#         AR = surf.bref/surf.cref

#         lhs = zeros(surf.nspan,surf.nbterm)
#         rhs = zeros(surf.nspan)
#         bcoeff = zeros(nsteps,surf.nbterm)

#         cnc_finite = zeros(nsteps)
#         cnnc_finite = zeros(nsteps)

#         # There is no apparent mass correction in this method

#         for i = 1:nsteps
#             for j = 1:surf.nspan
#                 for n = 1:surf.nbterm
#                     lhs[j,n] = sin(n*psi[j])*(sin(psi[j]) + (n*pi/(2*AR)))
#                 end
#                 rhs[j] = pi*sin(psi[j])*mat[j][i,9]/(2*AR)
#             end

#             bcoeff[i,:] = \(lhs, rhs)
#         end

#         a03d = zeros(nsteps,surf.nspan)
#         cd_ind = zeros(nsteps)

#         for i = 1:nsteps
#             cd_ind[i] = 0
#             for n = 1:surf.nbterm
#                 cd_ind[i] = cd_ind[i] + real(n)*bcoeff[i,n]^2
#             end
#             cd_ind[i] = cd_ind[i]*pi*AR
#             for j = 1:surf.nspan
#                 a03d[i,j] = 0

#                 for n = 1:surf.nbterm
#                     a03d[i,j] = a03d[i,j] - real(n)*bcoeff[i,n]*sin(n*psi[j])/sin(psi[j])
#                 end
#             end
#         end

#         cn = zeros(nsteps)
#         cs = zeros(nsteps)
#         cl = zeros(nsteps)
#         cd = zeros(nsteps)
#         cn3d = zeros(surf.nspan)
#         cs3d = zeros(surf.nspan)
#         cl3d = zeros(surf.nspan)
#         cd3d = zeros(surf.nspan)

#         for i = 1:nsteps
#             cn[i] = 0
#             cs[i] = 0
#             update_kinem(surf2d, mat[1][i,1])

#             for j = 1:surf.nspan
#                 cn3d[j] = mat[j][i,10] + (2*pi/surf.uref)*(mat[j][i,4]*cos(mat[j][i,2]) + surf2d.kinem.hdot*sin(mat[j][i,2]))*a03d[i,j]
#                 cs3d[j] = mat[j][i,11] + 2*pi*a03d[i,j]^2
#                 cl3d[j] = cn3d[j]*cos(mat[j][i,2]) + cs3d[j]*sin(mat[j][i,2])
#                 cd3d[j] = cn3d[j]*sin(mat[j][i,2]) - cs3d[j]*cos(mat[j][i,2])
#                 cn[i] = cn[i] + cn3d[j]*sin(psi[j])*dpsi/2
#                 cs[i] = cs[i] + cs3d[j]*sin(psi[j])*dpsi/2
#                 cl[i] = cl[i] + cl3d[j]*sin(psi[j])*dpsi/2
#                 cd[i] = cd[i] + cd3d[j]*sin(psi[j])*dpsi/2
#             end
#         end
#         return cl, cd, cd_ind, surf, field, mat2d, a03d

#     end


# end

# function QScorrect_ldvm(surf :: ThreeDSurf, field :: ThreeDFlowField, nsteps :: Int64, dtstar :: Float64)

#     mat = Array{Float64, 2}[]

#     if surf.kindef.vartype == "Constant"

#         # Kinematics at all strips is the same

#         kinem2d = KinemDef(surf.kindef.alpha, surf.kindef.h, surf.kindef.u)

#         surf2d = TwoDSurf(surf.patchdata[1].coord_file, surf.patchdata[1].pvt, kinem2d, [surf.patchdata[1].lc;])

#         #If 3D flow field is defined with disturbances or external vortices, these should be transferred to the 2D flowfield
#         curfield2d = TwoDFlowField()

#         mat2d, surf2d, curfield2d = lautat_wakeroll_more(surf2d, curfield2d, nsteps, dtstar)

#         #Fill the 2D solution into the 3D data structures

#         for i = 1:surf.nspan
#             push!(mat, mat2d)
#             surf.cam[i,:] = surf2d.cam
#             surf.cam_slope[i,:] = surf2d.cam_slope
#             surf.kinem[i,:] = surf2d.kinem
#             surf.bnd_x[i,:] = surf2d.bnd_x
#             surf.bnd_z[i,:] = surf2d.bnd_z
#             surf.uind[i,:] = surf2d.uind
#             surf.vind[i,:] = 0
#             surf.wind[i,:] = surf2d.wind
#             surf.downwash[i,:] = surf2d.downwash
#             surf.a0[i] = surf2d.a0[1]
#             surf.aterm[i,:] = surf2d.aterm
#             surf.a0dot[i] = surf2d.a0dot[1]
#             surf.adot[i,:] = surf2d.adot
#             surf.a0prev[i] = surf2d.a0prev[1]
#             surf.aprev[i,:] = surf2d.aprev
#             surf.levflag[i] = surf2d.levflag[1]

#             #Convert the bound vortices into 3d vortices and add them
#             for j = 1:length(surf2d.bv)
#                 surf.bv[i,j] = ThreeDVort([surf2d.bv[j].x; surf.yle[i]; surf2d.bv[j].z], [0.; surf2d.bv[j].s; 0.], surf2d.bv[j].vc, surf2d.bv[j].vx, 0., surf2d.bv[j].vz)
#             end
#         end

#         for j = 1:surf.nspan
#             for i = 1:length(curfield2d.lev)
#                 push!(field.lev, ThreeDVort([curfield2d.lev[i].x; surf.yle[j]; curfield2d.lev[i].z], [0.; curfield2d.lev[i].s; 0.], curfield2d.lev[i].vc, curfield2d.lev[i].vx, 0., curfield2d.lev[i].vz))
#             end
#             for i = 1:length(curfield2d.tev)
#                 push!(field.tev, ThreeDVort([curfield2d.tev[i].x; surf.yle[j]; curfield2d.tev[i].z], [0.; curfield2d.tev[i].s; 0.], curfield2d.tev[i].vc, curfield2d.tev[i].vx, 0., curfield2d.tev[i].vz))
#             end
#             for i = 1:length(curfield2d.extv)
#                 push!(field.extv, ThreeDVort([curfield2d.extv[i].x; surf.yle[j]; curfield2d.extv[i].z], [0.; curfield2d.extv[i].s; 0.], curfield2d.extv[i].vc, curfield2d.extv[i].vx, 0., curfield2d.extv[i].vz))
#             end
#         end

#         AR = surf.bref/surf.cref

#         lhs = zeros(surf.nspan,surf.nbterm)
#         rhs = zeros(surf.nspan)
#         bcoeff = zeros(nsteps,surf.nbterm)

#         cnc_finite = zeros(nsteps)
#         cnnc_finite = zeros(nsteps)

#         # There is no apparent mass correction in this method

#         for i = 1:nsteps
#             for j = 1:surf.nspan
#                 for n = 1:surf.nbterm
#                     lhs[j,n] = sin(n*surf.psi[j])*(sin(surf.psi[j]) + (n*pi/(2*AR)))
#                 end
#                 rhs[j] = pi*sin(surf.psi[j])*mat[j][i,9]/(2*AR)
#             end

#             bcoeff[i,:] = \(lhs, rhs)
#         end

#         a03d = zeros(nsteps,surf.nspan)
#         cd_ind = zeros(nsteps)

#         for i = 1:nsteps
#             cd_ind[i] = 0
#             for n = 1:surf.nbterm
#                 cd_ind[i] = cd_ind[i] + real(n)*bcoeff[i,n]^2
#             end
#             cd_ind[i] = cd_ind[i]*pi*AR
#             for j = 1:surf.nspan
#                 a03d[i,j] = 0

#                 for n = 1:surf.nbterm
#                     a03d[i,j] = a03d[i,j] - real(n)*bcoeff[i,n]*sin(n*surf.psi[j])/sin(surf.psi[j])
#                 end
#             end
#         end

#         cn = zeros(nsteps)
#         cs = zeros(nsteps)
#         cl = zeros(nsteps)
#         cd = zeros(nsteps)
#         cn3d = zeros(surf.nspan)
#         cs3d = zeros(surf.nspan)
#         cl3d = zeros(surf.nspan)
#         cd3d = zeros(surf.nspan)

#         for i = 1:nsteps
#             cn[i] = 0
#             cs[i] = 0
#             update_kinem(surf2d, mat[1][i,1])

#             for j = 1:surf.nspan
#                 cn3d[j] = mat[j][i,10] + (2*pi/surf.uref)*(mat[j][i,4]*cos(mat[j][i,2]) + surf2d.kinem.hdot*sin(mat[j][i,2]))*a03d[i,j]
#                 cs3d[j] = mat[j][i,11] + 2*pi*a03d[i,j]^2
#                 cl3d[j] = cn3d[j]*cos(mat[j][i,2]) + cs3d[j]*sin(mat[j][i,2])
#                 cd3d[j] = cn3d[j]*sin(mat[j][i,2]) - cs3d[j]*cos(mat[j][i,2])
#             end
#             for j = 1:surf.nspan-1
#                 cn[i] = cn[i] + 0.5*(cn3d[j] + cn3d[j+1])*sin(0.5*(surf.psi[j] + surf.psi[j+1]))*(surf.psi[j+1] - surf.psi[j])/2
#                 cs[i] = cs[i] + 0.5*(cs3d[j] + cs3d[j+1])*sin(0.5*(surf.psi[j] + surf.psi[j+1]))*(surf.psi[j+1] - surf.psi[j])/2
#                 cl[i] = cl[i] + 0.5*(cl3d[j] + cl3d[j+1])*sin(0.5*(surf.psi[j] + surf.psi[j+1]))*(surf.psi[j+1] - surf.psi[j])/2
#                 cd[i] = cd[i] + 0.5*(cd3d[j] + cd3d[j+1])*sin(0.5*(surf.psi[j] + surf.psi[j+1]))*(surf.psi[j+1] - surf.psi[j])/2
#             end
#         end
#         return cl, cd, cd_ind, surf, field, mat2d, a03d

#     end

# end

function QSLLT_lautat(surf :: ThreeDSurfSimple, field :: ThreeDFieldSimple, nsteps :: Int64, dtstar :: Float64)

    mat = Array(Float64, 0, 4)

    mat = mat'

    dt = dtstar*surf.cref/surf.uref

    t = 0.

    cl = zeros(surf.nspan)
    cd = zeros(surf.nspan)
    cm = zeros(surf.nspan)

    lhs = zeros(surf.nspan, surf.nspan)
    rhs = zeros(surf.nspan)
    bcoeff = zeros(surf.nspan)

    for istep = 1:nsteps
        #Udpate current time
        t = t + dt

        for i = 1:surf.nspan
            #Define the flow field
            push!(field.f2d, TwoDFlowField())
            #Update kinematic parameters
            update_kinem(surf.s2d[i], t)

            #Update flow field parameters if any
            update_externalvel(field.f2d[i], t)

            #Update bound vortex positions
            update_boundpos(surf.s2d[i], dt)

            #Add a TEV with dummy strength
            place_tev(surf.s2d[i], field.f2d[i], dt)
        end

        kelv = KelvinConditionLLTldvm(surf, field)

        #Solve for TEV strength to satisfy Kelvin condition
        
        soln = nlsolve(not_in_place(kelv), [-0.01*ones(surf.nspan); zeros(surf.nspan)])
        
        for i = 1:surf.nspan
            field.f2d[i].tev[length(field.f2d[i].tev)].s = soln.zero[i]

            #Update incduced velocities on airfoil
            update_indbound(surf.s2d[i], field.f2d[i])
            #Calculate downwash
            update_downwash(surf.s2d[i], [field.f2d[i].u[1],field.f2d[i].w[1]])

            #Calculate first two fourier coefficients
            update_a0anda1(surf.s2d[i])
            surf.bc[i] = surf.s2d[i].a0[1] + 0.5*surf.s2d[i].aterm[1]
        end

        for i = 1:surf.nspan
            surf.a03d[i] = 0
            for n = 1:surf.nspan
                nn = 2*n - 1
                surf.a03d[i] = surf.a03d[i] - real(nn)*soln.zero[n+surf.nspan]*sin(nn*surf.psi[i])/sin(surf.psi[i])
            end
        end


        for i = 1:surf.nspan
            #Update 3D effect on A0
            surf.s2d[i].a0[1] = surf.s2d[i].a0[1] + surf.a03d[i]

            #Update rest of Fourier terms
            update_a2toan(surf.s2d[i])

            #Update derivatives of Fourier coefficients
            update_adot(surf.s2d[i],dt)

            #Set previous values of aterm to be used for derivatives in next time step
            surf.s2d[i].a0prev[1] = surf.s2d[i].a0[1]
            for ia = 1:3
                surf.s2d[i].aprev[ia] = surf.s2d[i].aterm[ia]
            end

            #Calculate bound vortex strengths
            update_bv(surf.s2d[i])

            # #Remove vortices that are far away from airfoil
            # if (delvort.flag == 1)
            #     if length(field2d[i].tev) > delvort.limit
            #         if (sqrt((field2d[i].tev[1].x- surf2d[i].bnd_x[div(surf2d[i].ndiv,2)])^2 + (field2d[i].tev[1].z- surf2d[i].bnd_z[div(surf2d[i].ndiv,2)])^2) > delvort.dist*surf2d[i].c)
            #             kelv_enf = kelv_enf + field2d[i].tev[1].s
            #             for i = 1:length(field2d[i].tev)-1
            #                 field2d[i].tev[i] = field2d[i].tev[i+1]
            #             end
            #             pop!(field2d[i].tev)
            #         end
            #     end
            #     if length(field2d[i].lev) > delvort.limit
            #         if (sqrt((field2d[i].lev[1].x- surf2d[i].bnd_x[div(surf2d[i].ndiv,2)])^2 + (field2d[i].lev[1].z- surf2d[i].bnd_z[div(surf2d[i].ndiv,2)])^2) > delvort.dist*surf2d[i].c)
            #             kelv_enf = kelv_enf + field2d[i].lev[1].s
            #             for i = 1:length(field2d[i].lev)-1
            #                 field2d[i].lev[i] = field2d[i].lev[i+1]
            #             end
            #         pop!(field2d[i].lev)
            #         end
            #     end
            # end
            wakeroll(surf.s2d[i], field.f2d[i], dt)
            
            if (surf.s2d[i].levflag[1] == 1)
                cl[i], cd[i], cm[i] = calc_forces_E(surf.s2d[i], field.f2d[i].lev[length(field.f2d[i].lev)].s, dt)
            else
                cl[i], cd[i], cm[i] = calc_forces(surf.s2d[i])
            end

        end

        cl3d = 0
        cd3d = 0
        cm3d = 0

        for i = 1:surf.nspan-1
            cl3d = cl3d + 0.5*(cl[i] + cl[i+1])*sin(0.5*(surf.psi[i] + surf.psi[i+1]))*(surf.psi[i+1] - surf.psi[i])/2
            cd3d = cd3d + 0.5*(cd[i] + cd[i+1])*sin(0.5*(surf.psi[i] + surf.psi[i+1]))*(surf.psi[i+1] - surf.psi[i])/2
            cm3d = cm3d + 0.5*(cm[i] + cm[i+1])*sin(0.5*(surf.psi[i] + surf.psi[i+1]))*(surf.psi[i+1] - surf.psi[i])/2
        end

        mat = hcat(mat, [t, cl3d, cd3d, cm3d])
    end
    mat = mat'
    mat, surf, field

end

function QSLLT_lautat2(surf :: ThreeDSurfSimple, curfield :: ThreeDFieldStrip, nsteps :: Int64, dtstar :: Float64 = 0.015; delvort = DelVortDef(0, 0, 0.), mat = Array(Float64, 0, 4), kelv_enf = zeros(surf.nspan), writefile = "Nil", writeind :: Vector{Int} = Int[])
    
    #Alternate way of solving the newton iteration. For lautat both
    #types of iteration work, but for ldvm only the type used in this
    #code will work.

    if (size(mat,1) > 0)
        t = mat[end,1]
    else
        t = 0.
    end

    mat = mat'

    dt = dtstar*surf.cref/surf.uref

    eps = 1e-6
    
    cl = zeros(surf.nspan)
    cd = zeros(surf.nspan)
    cm = zeros(surf.nspan)
    a03dprev = zeros(surf.nspan)
    
    shed_ind = Vector{Int}[]
    
    #If required, write an archive file
    if writefile != "Nil"
        jldopen(writefile, "w") do file
            write(file, "dt",  dt)
            write(file, "dtstar",  dt)
            write(file, "Insurf",  surf)
            write(file, "Infield", curfield)
            write(file, "delvort", delvort)
            write(file, "nsteps", nsteps)
            write(file, "writeind", writeind)
            g = g_create(file, "Init")
            g = write_stamp(surf, curfield, t, kelv_enf, g)
        end
    end
    
    for istep = 1:nsteps
        #Udpate current time
        t = t + dt
        
        #Update kinematic parameters
        update_kinem(surf.s2d, t)
        
        #Update bound vortex positions
        update_boundpos(surf.s2d, dt)

        #Update IC if required
        
        #Update flow field parameters if any
        update_externalvel(curfield, t)

        place_tev(surf.s2d, curfield, dt)
        
        iter = 0
        while true
            iter += 1
            if iter > 1
                a03dprev[:] = surf.a03d[:]
                tev_prev = map(q->q.s, curfield.tev[end])
            end
            
            for i = 1:surf.nspan
                kelv = KelvinConditionLLTldvmSep(surf, curfield, shed_ind, i, kelv_enf)
                soln = nlsolve(not_in_place(kelv), [-0.01])
                curfield.tev[end][i].s = soln.zero[1]
            end

            calc_a03d(surf)
            
            if iter > 1
                teverr = mean(abs(map(q->q.s, curfield.tev[end]) - tev_prev))
                a03derr = mean(abs(surf.a03d - a03dprev))
                
                if (teverr < eps && a03derr < eps)
                    break
                end
                if (iter > 20)
                    println("1D convergence failed")
                    println(teverr, a03derr)
                    break
                end
            end
        end
        
        #Update rest of Fourier terms
        update_a2toan(surf.s2d)
        
        #Update derivatives of Fourier coefficients
        update_adot(surf.s2d,dt)

        #Update a03ddot
        for i = 1:surf.nspan
            surf.a03ddot[i] = (surf.a03d[i] - surf.a03dprev[i])/dt
        end
        
        #Set previous values of aterm to be used for derivatives in next time step
        for i = 1:surf.nspan
            surf.s2d[i].a0prev[1] = surf.s2d[i].a0[1]
            for ia = 1:3
                surf.s2d[i].aprev[ia] = surf.s2d[i].aterm[ia]
            end
            surf.a03dprev[i] = surf.a03d[i]
        end
        
        #For the convenience of not changing 2D codes, a03d is simply
        #added to a0 as it has the same influence. For proper
        #accounting however, it is subtracted again so that the terms
        #are kept separate

        for i = 1:surf.nspan
            #Update 3D effect on A0
            surf.s2d[i].a0[1] = surf.s2d[i].a0[1] + surf.a03d[i]
            surf.s2d[i].a0dot[1] = surf.s2d[i].a0dot[1] + surf.a03ddot[i]
        end
                
        #Calculate bound vortex strengths
        update_bv(surf.s2d)
        
        #Remove vortices that are far away from airfoil
        if (delvort.flag == 1)
            if length(curfield.tev) > delvort.limit
                #check at the first spanwise location and if satisfied delete the whole row
                if (sqrt((curfield.tev[1][1].x -  surf.s2d[1].bnd_x[div(surf.ndiv,2)])^2 + (curfield.tev[1][1].z - surf.s2d[1].bnd_z[div(surf.ndiv,2)])^2) > delvort.dist*surf.cref)
                    for i = 1:surf.nspan
                        kelv_enf[i] = kelv_enf[i] + curfield.tev[1][i].s
                        for j = 1:length(curfield.tev) - 1
                            curfield.tev[j][i] = curfield.tev[j+1][i]
                        end
                    end
                    pop!(curfield.tev)
                end
            end

            if length(curfield.lev) > delvort.limit
                if (sqrt((curfield.lev[1][1].x - surf.s2d[1].bnd_x[div(surf.ndiv,2)])^2 + (curfield.lev[1][1].z - surf.s2d[1].bnd_z[div(surf.ndiv,2)])^2) > delvort.dist*surf.cref)
                    for i = 1:surf.nspan
                        kelv_enf[i] = kelv_enf[i] + curfield.lev[1][i].s

                        for i = 1:length(curfield.lev)-1
                            curfield.lev[j][i] = curfield.lev[j+1][i]
                        end
                    end
                    pop!(curfield.lev)
                end
            end
        end

        wakeroll(surf.s2d, curfield, dt, shed_ind)

        #bunched a0 with a03d is used here too. Remove after
        cl, cd, cm = calc_forces(surf.s2d)
        
        cl3d = 0
        cd3d = 0
        cm3d = 0
        
        for i = 1:surf.nspan-1
            cl3d = cl3d + 0.5*(cl[i] + cl[i+1])*sin(0.5*(surf.psi[i] + surf.psi[i+1]))*(surf.psi[i+1] - surf.psi[i])
            cd3d = cd3d + 0.5*(cd[i] + cd[i+1])*sin(0.5*(surf.psi[i] + surf.psi[i+1]))*(surf.psi[i+1] - surf.psi[i])
            cm3d = cm3d + 0.5*(cm[i] + cm[i+1])*sin(0.5*(surf.psi[i] + surf.psi[i+1]))*(surf.psi[i+1] - surf.psi[i])
            
        end
        
        mat = hcat(mat, [t, cl3d, cd3d, cm3d])
        
        #If required, write an archive file
        if writefile != "Nil"
            if istep in writeind
                jldopen(writefile, "r+") do file
                    g = g_create(file, "t$istep")
                    g = write_stamp(surf, curfield, t, kelv_enf, g)
                end
            end
        end

for i = 1:surf.nspan
    #Update 3D effect on A0
    surf.s2d[i].a0[1] = surf.s2d[i].a0[1] - surf.a03d[i]
    surf.s2d[i].a0dot[1] = surf.s2d[i].a0dot[1] - surf.a03ddot[i]
end

end
mat = mat'
mat, surf, curfield, kelv_enf

end

function QSWeiss_lautat(surf :: ThreeDSurfWeiss, curfield :: ThreeDFieldStrip, nsteps :: Int64, dtstar :: Float64 = 0.015; delvort = DelVortDef(0, 0, 0.), mat = Array(Float64, 0, 4), kelv_enf = zeros(surf.nlat), writefile = "Nil", writeind :: Vector{Int} = Int[])
    
    #Alternate way of solving the newton iteration. For lautat both
    #types of iteration work, but for ldvm only the type used in this
    #code will work.

    if (size(mat,1) > 0)
        t = mat[end,1]
    else
        t = 0.
    end

    mat = mat'

    dt = dtstar*surf.cref/surf.uref

    eps = 1e-6
    
    cl = zeros(surf.nlat)
    cd = zeros(surf.nlat)
    cm = zeros(surf.nlat)
    a03dprev = zeros(surf.nlat)
    
    shed_ind = Vector{Int}[]
    
    #If required, write an archive file
    if writefile != "Nil"
        jldopen(writefile, "w") do file
            write(file, "dt",  dt)
            write(file, "dtstar",  dt)
            write(file, "Insurf",  surf)
            write(file, "Infield", curfield)
            write(file, "delvort", delvort)
            write(file, "nsteps", nsteps)
            write(file, "writeind", writeind)
            g = g_create(file, "Init")
            g = write_stamp(surf, curfield, t, kelv_enf, g)
        end
    end

    #The initial value of plunge (const through chord) used to update IC
    if typeof(surf.kindef.h) == BendingDef 
        hs_prev = zeros(surf.nlat)
        he_prev = zeros(surf.nlat)
        for i = 1:surf.nlat
            if i == 1
                hs_prev[i] = surf.s2d[i].kinem.h - 0.5*(surf.s2d[i+1].kinem.h - surf.s2d[i].kinem.h)
            else
                hs_prev[i] = surf.s2d[i].kinem.h - 0.5*(surf.s2d[i].kinem.h - surf.s2d[i-1].kinem.h)
            end
            if i == nlat
                he_prev[i] = surf.s2d[i].kinem.h + 0.5*(surf.s2d[i].kinem.h - surf.s2d[i-1].kinem.h)
            else
                he_prev[i] = surf.s2d[i].kinem.h + 0.5*(surf.s2d[i+1].kinem.h - surf.s2d[i].kinem.h)
            end
        end       
    end
    
    for istep = 1:nsteps
        #Udpate current time
        t = t + dt
        
        #Update kinematic parameters
        update_kinem(surf.s2d, t)
        
        #Update bound vortex positions
        update_boundpos(surf.s2d, dt)
        
        #Update IC if required
        if typeof(surf.kindef.h) == BendingDef 
            update_IC(surf, hs_prev, he_prev)

            for i = 1:surf.nlat
                if i == 1
                    hs_prev[i] = surf.s2d[i].kinem.h - 0.5*(surf.s2d[i+1].kinem.h - surf.s2d[i].kinem.h)
                else
                    hs_prev[i] = surf.s2d[i].kinem.h - 0.5*(surf.s2d[i].kinem.h - surf.s2d[i-1].kinem.h)
                end
                if i == nlat
                    he_prev[i] = surf.s2d[i].kinem.h + 0.5*(surf.s2d[i].kinem.h - surf.s2d[i-1].kinem.h)
                else
                    he_prev[i] = surf.s2d[i].kinem.h + 0.5*(surf.s2d[i+1].kinem.h - surf.s2d[i].kinem.h)
                end
            end
        end
        
        #Update flow field parameters if any
        update_externalvel(curfield, t)
        
        place_tev(surf.s2d, curfield, dt)
        
        iter = 0
        while true
            iter += 1
            if iter > 1
                a03dprev[:] = surf.a03d[:]
                tev_prev = map(q->q.s, curfield.tev[end])
            end
            
            for i = 1:surf.nlat
                kelv = KelvinConditionWeiss(surf, curfield, shed_ind, i, kelv_enf)
                soln = nlsolve(not_in_place(kelv), [-0.01])
                curfield.tev[end][i].s = soln.zero[1]
            end
            
            surf.a03d[:] = surf.IC*surf.gamma/pi
            surf.a13d[:] = -surf.IC*surf.gamma/pi
            
            if iter > 1
                teverr = mean(abs(map(q->q.s, curfield.tev[end]) - tev_prev))
                a03derr = mean(abs(surf.a03d - a03dprev))
                
                if (teverr < eps && a03derr < eps)
                    break
                end
                if (iter > 20)
                    println("1D convergence failed")
                    println(teverr, a03derr)
                    break
                end
            end
        end
        
        #Update rest of Fourier terms
        update_a2toan(surf.s2d)
        
        #Update derivatives of Fourier coefficients
        update_adot(surf.s2d,dt)
        
        #Update a03ddot
        for i = 1:surf.nlat
            surf.a03ddot[i] = (surf.a03d[i] - surf.a03dprev[i])/dt
        end
        
        #Set previous values of aterm to be used for derivatives in next time step
        for i = 1:surf.nlat
            surf.s2d[i].a0prev[1] = surf.s2d[i].a0[1]
            for ia = 1:3
                surf.s2d[i].aprev[ia] = surf.s2d[i].aterm[ia]
            end
            surf.a03dprev[i] = surf.a03d[i]
        end
        
        #For the convenience of not changing 2D codes, a03d is simply
        #added to a0 as it has the same influence. For proper
        #accounting however, it is subtracted again so that the terms
        #are kept separate
        
        for i = 1:surf.nlat
            #Update 3D effect on A0
            surf.s2d[i].a0[1] = surf.s2d[i].a0[1] + surf.a03d[i]
            surf.s2d[i].a0dot[1] = surf.s2d[i].a0dot[1] + surf.a03ddot[i]
        end
        
        #Calculate bound vortex strengths
        update_bv(surf.s2d)
        
        #Remove vortices that are far away from airfoil
        if (delvort.flag == 1)
            if length(curfield.tev) > delvort.limit
                #check at the first spanwise location and if satisfied delete the whole row
                if (sqrt((curfield.tev[1][1].x -  surf.s2d[1].bnd_x[div(surf.ndiv,2)])^2 + (curfield.tev[1][1].z - surf.s2d[1].bnd_z[div(surf.ndiv,2)])^2) > delvort.dist*surf.cref)
                    for i = 1:surf.nlat
                        kelv_enf[i] = kelv_enf[i] + curfield.tev[1][i].s
                        for j = 1:length(curfield.tev) - 1
                            curfield.tev[j][i] = curfield.tev[j+1][i]
                        end
                    end
                    pop!(curfield.tev)
                end
            end
            
            if length(curfield.lev) > delvort.limit
                if (sqrt((curfield.lev[1][1].x - surf.s2d[1].bnd_x[div(surf.ndiv,2)])^2 + (curfield.lev[1][1].z - surf.s2d[1].bnd_z[div(surf.ndiv,2)])^2) > delvort.dist*surf.cref)
                    for i = 1:surf.nlat
                        kelv_enf[i] = kelv_enf[i] + curfield.lev[1][i].s
                        
                        for i = 1:length(curfield.lev)-1
                            curfield.lev[j][i] = curfield.lev[j+1][i]
                        end
                    end
                    pop!(curfield.lev)
                end
            end
        end
        
        wakeroll(surf.s2d, curfield, dt, shed_ind)
        
        #bunched a0 with a03d is used here too. Remove after
        cl, cd, cm = calc_forces(surf.s2d)
        
        cl3d = sum(cl.*surf.chord.*surf.ds)/surf.sref
cd3d = sum(cd.*surf.chord.*surf.ds)/surf.sref
cm3d = sum(cm.*surf.chord.*surf.ds)/surf.sref



mat = hcat(mat, [t, cl3d, cd3d, cm3d])

#If required, write an archive file
if writefile != "Nil"
    if istep in writeind
        jldopen(writefile, "r+") do file
            g = g_create(file, "t$istep")
            g = write_stamp(surf, curfield, t, kelv_enf, g)
        end
    end
end

for i = 1:surf.nlat
    #Update 3D effect on A0
    surf.s2d[i].a0[1] = surf.s2d[i].a0[1] - surf.a03d[i]
    surf.s2d[i].a0dot[1] = surf.s2d[i].a0dot[1] - surf.a03ddot[i]
end

end
mat = mat'
mat, surf, curfield, kelv_enf

end

    
    
    # function QSLLT_ldvm(surf :: ThreeDSurf, field :: ThreeDFieldSimple, nsteps :: Int64, dtstar :: Float64)

#     mat = Array(Float64, 0, 4)

#     mat = mat'

#     surf2d = TwoDSurf[]
#     field2d = TwoDFlowField[]
#     kinem2d = KinemDef[]

#     dt = dtstar*surf.cref/surf.uref

#     t = 0.

#     AR = surf.bref/surf.cref

#     bc = zeros(surf.nspan)
#     a03d = zeros(surf.nspan)
#     cl = zeros(surf.nspan)
#     cd = zeros(surf.nspan)
#     cm = zeros(surf.nspan)

#     lhs = zeros(surf.nspan, surf.nspan)
#     rhs = zeros(surf.nspan)
#     bcoeff = zeros(surf.nspan)

#     if surf.kindef.vartype == "Constant"

#         for i = 1:surf.nspan
#             # Kinematics at all strips is the same

#             push!(kinem2d, KinemDef(surf.kindef.alpha, surf.kindef.h, surf.kindef.u))
#             push!(surf2d, TwoDSurf(surf.patchdata[1].coord_file, surf.patchdata[1].pvt, kinem2d[i], [surf.patchdata[1].lc;]))
#             #If 3D flow field is defined with disturbances or external vortices, these should be transferred to the 2D flowfield
#             push!(field2d, TwoDFlowField())
#         end
#     end

#     for istep = 1:nsteps
#         #Udpate current time
#         t = t + dt

#         for i = 1:surf.nspan
#             #Update kinematic parameters
#             update_kinem(surf2d[i], t)

#             #Update flow field parameters if any
#             update_externalvel(field2d[i], t)

#             #Update bound vortex positions
#             update_boundpos(surf2d[i], dt)

#             #Add a TEV with dummy strength
#             place_tev(surf2d[i], field2d[i], dt)
#         end

#         kelv = KelvinConditionLLTldvm(surf, surf2d, field2d)

#         #Solve for TEV strength to satisfy Kelvin condition

#         soln = nlsolve(not_in_place(kelv), [-0.01*ones(surf.nspan); zeros(surf.nspan)])

#         for i = 1:surf.nspan
#             field2d[i].tev[length(field2d[i].tev)].s = soln.zero[i]

#             #Update incduced velocities on airfoil
#             update_indbound(surf2d[i], field2d[i])

#             #Calculate downwash
#             update_downwash(surf2d[i], [field2d[i].u[1],field2d[i].w[1]])

#             #Calculate first two fourier coefficients
#             update_a0anda1(surf2d[i])
#         end

#         for i = 1:surf.nspan
#             a03d[i] = 0
#             for n = 1:surf.nspan
#                 nn = 2*n - 1
#                 a03d[i] = a03d[i] - real(nn)*soln.zero[n+surf.nspan]*sin(nn*surf.psi[i])/sin(surf.psi[i])
#             end
#         end

#         nshed = Int(0)
#         for i = 1:surf.nspan
#             #Update 3D effect on A0
#             surf2d[i].a0[1] = surf2d[i].a0[1] + a03d[i]

#             #2D iteration if LESP_crit is exceeded
#             if abs(surf2d[i].a0[1]) > surf2d[i].lespcrit[1]
#                 #Remove the previous tev
#                 pop!(field2d[i].tev)
#                 #Add a TEV with dummy strength
#                 place_tev(surf2d[i],field2d[i],dt)

#                 #Add a LEV with dummy strength
#                 place_lev(surf2d[i],field2d[i],dt)
#                 surf2d[i].levflag[1] = 1
#                 nshed += 1
#             else
#                 surf2d[i].levflag[1] = 0
#             end
#         end
#         if nshed > 0
#             kelvkutta = KelvinKuttaLLTldvm(surf,surf2d,field2d, nshed)

#             #Solve for TEV and LEV strengths to satisfy Kelvin condition and Kutta condition at leading edge
#             soln = nlsolve(not_in_place(kelvkutta), [-0.01*ones(surf.nspan); 0.01*ones(nshed); zeros(surf.nspan)])

#             cntr = surf.nspan + 1

#             for i = 1:surf.nspan
#                 field2d[i].tev[length(field2d[i].tev)].s = soln.zero[i]
#             end
#             for i = 1:surf.nspan
#                 if surf2d[i].levflag[1] == 1
#                     field2d[i].lev[length(field2d[i].lev)].s = soln.zero[cntr]
#                     cntr += 1
#                 end
#             end

#             for i = 1:surf.nspan
#                 #Update incduced velocities on airfoil
#                 update_indbound(surf2d[i], field2d[i])

#                 #Calculate downwash
#                 update_downwash(surf2d[i], [field2d[i].u[1],field2d[i].w[1]])

#                 #Calculate first two fourier coefficients
#                 update_a0anda1(surf2d[i])
#             end

#             for i = 1:surf.nspan
#                 a03d[i] = 0
#                 for n = 1:surf.nspan
#                     nn = 2*n - 1
#                     a03d[i] = a03d[i] - real(nn)*soln.zero[n+surf.nspan+nshed]*sin(nn*surf.psi[i])/sin(surf.psi[i])
#                 end
#             end

#             for i = 1:surf.nspan
#                 #Update 3D effect on A0
#                 surf2d[i].a0[1] = surf2d[i].a0[1] + a03d[i]
#             end
#         end
#         for i = 1:surf.nspan
#             #Update rest of Fourier terms
#             update_a2toan(surf2d[i])

#             #Update derivatives of Fourier coefficients
#             update_adot(surf2d[i],dt)

#             #Set previous values of aterm to be used for derivatives in next time step
#             surf2d[i].a0prev[1] = surf2d[i].a0[1]
#             for ia = 1:3
#                 surf2d[i].aprev[ia] = surf2d[i].aterm[ia]
#             end

#             #Calculate bound vortex strengths
#             update_bv(surf2d[i])

#             # #Remove vortices that are far away from airfoil
#             # if (delvort.flag == 1)
#             #     if length(field2d[i].tev) > delvort.limit
#             #         if (sqrt((field2d[i].tev[1].x- surf2d[i].bnd_x[div(surf2d[i].ndiv,2)])^2 + (field2d[i].tev[1].z- surf2d[i].bnd_z[div(surf2d[i].ndiv,2)])^2) > delvort.dist*surf2d[i].c)
#             #             kelv_enf = kelv_enf + field2d[i].tev[1].s
#             #             for i = 1:length(field2d[i].tev)-1
#             #                 field2d[i].tev[i] = field2d[i].tev[i+1]
#             #             end
#             #             pop!(field2d[i].tev)
#             #         end
#             #     end
#             #     if length(field2d[i].lev) > delvort.limit
#             #         if (sqrt((field2d[i].lev[1].x- surf2d[i].bnd_x[div(surf2d[i].ndiv,2)])^2 + (field2d[i].lev[1].z- surf2d[i].bnd_z[div(surf2d[i].ndiv,2)])^2) > delvort.dist*surf2d[i].c)
#             #             kelv_enf = kelv_enf + field2d[i].lev[1].s
#             #             for i = 1:length(field2d[i].lev)-1
#             #                 field2d[i].lev[i] = field2d[i].lev[i+1]
#             #             end
#             #         pop!(field2d[i].lev)
#             #         end
#             #     end
#             # end
#             wakeroll(surf2d[i], field2d[i], dt)

#             if (surf2d[i].levflag[1] == 1)
#                 cl[i], cd[i], cm[i] = calc_forces_E(surf2d[i], field2d[i].lev[length(field2d[i].lev)].s, dt)
#             else
#                 cl[i], cd[i], cm[i] = calc_forces(surf2d[i])
#             end

#         end

#         cl3d = 0
#         cd3d = 0
#         cm3d = 0

#         for i = 1:surf.nspan-1
#             cl3d = cl3d + 0.5*(cl[i] + cl[i+1])*sin(0.5*(surf.psi[i] + surf.psi[i+1]))*(surf.psi[i+1] - surf.psi[i])/2
#             cd3d = cd3d + 0.5*(cd[i] + cd[i+1])*sin(0.5*(surf.psi[i] + surf.psi[i+1]))*(surf.psi[i+1] - surf.psi[i])/2
#             cm3d = cm3d + 0.5*(cm[i] + cm[i+1])*sin(0.5*(surf.psi[i] + surf.psi[i+1]))*(surf.psi[i+1] - surf.psi[i])/2
#         end

#         mat = hcat(mat, [t, cl3d, cd3d, cm3d])
#     end
#     mat = mat'
#     mat, surf2d, field2d

# end



