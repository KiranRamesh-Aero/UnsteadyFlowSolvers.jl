#=
LVE.jl
    Written by Matthew White
    5/21/2019
=#
function LVE(surf::TwoDSurf,curfield::TwoDFlowField,nsteps::Int64 = 500,dtstar::Float64 = 0.015,frameCnt = 20,view = "square",LEV_Stop = 1000)
    ## Initialize variables
    mat = zeros(8 ,0) # loads output matrix
    prevCirc = zeros(surf.ndiv-1)
    circChange = zeros(surf.ndiv-1)
    vor_loc = zeros(surf.ndiv-1,2) # (x1,z1 ; x2,z2 ...) Inertial frame
    coll_loc = zeros(surf.ndiv-1,2) # (x1,z1 ; x2,z2 ...)
    global RHS = zeros(surf.ndiv)
    global relVel = zeros(surf.ndiv-1,2)
    IFR_vor = 0 .* vor_loc # Global frame
    IFR_field = TwoDFlowField()
    IFR_surf = TwoDSurf(surf.coord_file,surf.pvt,surf.kindef,surf.lespcrit; ndiv = surf.ndiv, camberType = "linear")
    global dp = zeros(surf.ndiv-1,1) # change in panel pressure
    aniStep = round((nsteps-1) / frameCnt) # Frames to capture animation on
    aniStep < 1 ? aniStep = 1 : aniStep
    frames = 0
    prevNow = 0
    LEVhist = hcat([],[])
    TEVhist = hcat([],[])


    #   length of each panel
    ds = sqrt.(( surf.x[1:end-1]-surf.x[2:end]).^2 + (surf.cam[1:end-1]-surf.cam[2:end]).^2)
    #   Surf cam_slope does not give the correct panel locations
    cam_slope = asind.((surf.cam[2:end]-surf.cam[1:end-1])./ds) # [deg]

    # Normal vectors for each panel
    n = hcat(-sind.(cam_slope),cosd.(cam_slope))
    # Tangential vectors for each panel
    tau = hcat(cosd.(cam_slope),sind.(cam_slope))

    ## Vortex Locations
    # Located at 25% of each panel
    #   Vortex loactions (at 25% of panel length)
    vor_loc[:,1] = surf.x[1:end-1] + .25 .* ds .* cosd.(cam_slope)
    vor_loc[:,2] = surf.cam[1:end-1] + .25 .* ds .* sind.(cam_slope)
    #   Collocation points (at 75% panel chordwise)
    coll_loc[:,1] = surf.x[1:end-1] + .75 .* ds .* cosd.(cam_slope)
    coll_loc[:,2] = surf.cam[1:end-1] + .75 .* ds .* sind.(cam_slope)

    ## Critical gamma (LEV)
    x = 1 - 2*ds[1] / surf.c
    gamma_crit = surf.lespcrit / 1.13 * ( surf.kinem.u * surf.c * ( acos(x) + sin(acos(x) )))
    gamma_crit = gamma_crit[1]

    ## Refresh vortex values
    refresh_vortex(surf,vor_loc)

    ## Time Stepping Loop
    t = -dtstar

    for i = 0:nsteps
        t += dtstar

        # Update Kinematics
        update_kinem(surf, t)

        # Inertial reference frame
        IFR_vor[:,1],IFR_vor[:,2] = IFR(surf,vor_loc[:,1],vor_loc[:,2],t)
        IFR_surf = refresh_vortex(IFR_surf,IFR_vor)

        ## TEV setup
        place_tev2(IFR_surf,IFR_field,dtstar,t)

        ## Influence Coeffcients
        x_w, z_w = BFR(surf, IFR_field.tev[end].x, IFR_field.tev[end].z, t)
        global a = influence_coeff(surf,curfield,coll_loc,n,dtstar,x_w,z_w)

        ## RHS Vector
        RHS[end] = sum(prevCirc) # previous total circulation
        # u_w, w_w
        if length(curfield.tev) > 0 # Trailing edge vortexs exist
            surf.uind[1:end-1], surf.wind[1:end-1] = ind_vel([curfield.tev; curfield.lev],coll_loc[:,1],coll_loc[:,2])
        end
        for j = 1:surf.ndiv-1
            alpha = surf.kinem.alpha
            # relVel = [U(t) ; W(t)]
            global relVel[j,:] = [cos(alpha) -sin(alpha) ; sin(alpha) cos(alpha)] * [surf.kinem.u ; -surf.kinem.hdot] + [-surf.kinem.alphadot*coll_loc[j,2] ; surf.kinem.alphadot * (coll_loc[j,1] - surf.pvt)]

            # RHS = -[U_t + u_w , W_t + w_w] dot n
            global RHS[j] = -((relVel[j,1] + surf.uind[j])*n[j,1] + (relVel[j,2] + surf.wind[j])*n[j,2])
        end

        ## Vortex Strenghths a*circ = RHS
        global circ = a\RHS

        ## Circulation changes before each panel
        for j = 1:surf.ndiv-1
            gamma1 = sum(prevCirc[1:j-1])
            gamma2 = sum(circ[1:j-1])
            circChange[j] = (gamma2-gamma1) ./ dtstar
        end

        ## LEV Ejection
        #   Test for LEV criteria
        if abs(circ[1]) > gamma_crit
            if surf.levflag[1] == 0 && t < LEV_Stop # if lev is first in batch
                gamma_crit_use = abs(circ[1])/circ[1] * gamma_crit # = gamma_crit w/ sign of circ[1]
                global t_start = t
                #println("LEV start batch:")
            else # if lev sheds but is not first in batch
                m_slp = -surf.lespcrit / (LEV_Stop - t_start)
                c_slp = surf.lespcrit - (m_slp*t_start)
                LESP_crit_use = (m_slp*t) + c_slp
                x = 1 - 2*ds[1] / surf.c
                gamma_crit_use = abs(circ[1])/circ[1] * LESP_crit_use / 1.13 * ( surf.kinem.u * surf.c * ( acos(x) + sin(acos(x) ))) # abs(circ[1])/circ[1] keeps sign of original LEV
            end
            #println("LEV Ejection: Step $i, time $t")
            gamma_crit_use = gamma_crit_use[1]

            # Place LEV
            place_lev2(surf,IFR_field,dtstar,t)

            surf.levflag[1] = 1 # Set LEV flag to track constant LEV streams

            # Convert coords to BFR
            x_lev, z_lev = BFR(surf, IFR_field.lev[end].x, IFR_field.lev[end].z, t)
            global a,a1 = mod_influence_coeff(surf,curfield,coll_loc,n,dtstar,x_w,z_w,x_lev,z_lev)

            # Recalculate RHS
            # u_w, w_w
            if length(curfield.tev) > 0 # Trailing edge vortexs exist
                surf.uind[1:end-1], surf.wind[1:end-1] = ind_vel([curfield.tev; curfield.lev],coll_loc[:,1],coll_loc[:,2])
            end
            for j = 1:surf.ndiv-1
                # RHS = -[U_t + u_w , W_t + w_w] dot n
                global RHS[j] = -((relVel[j,1] + surf.uind[j])*n[j,1] + (relVel[j,2] + surf.wind[j])*n[j,2]) - gamma_crit_use*a1[j]
            end
            RHS[end] -= gamma_crit_use # first RHS - gamma_crit_use*a1[end] (which is always 1)

            global circ = a\RHS # Solve new circluations

            IFR_field.lev[end].s = circ[end] #store LEV strength

            circ[2:end] = circ[1:end-1] # shift circs
            circ[1] = gamma_crit_use

            # LEV history data [t,s]
            stamp = hcat(t,IFR_field.lev[end].s)
            LEVhist = vcat(LEVhist,stamp)
        else
            surf.levflag[1] = 0
        end

        prevCirc = circ[1:end-1]

        for j = 1:surf.ndiv-1 # Set bv circs
            surf.bv[j].s = circ[j]
            IFR_surf.bv[j].s = circ[j]
        end
        if length(IFR_field.tev) > 1 # Set TEV circ
            IFR_field.tev[end].s = circ[end]
        end

        # LEV history data [t,s]
        stamp = hcat(t,IFR_field.tev[end].s)
        TEVhist = vcat(TEVhist,stamp)

        # Position mesh
        if i % aniStep == 0. && Int(i / aniStep) < frameCnt

            temp = Int(i / aniStep) + 1
            now = Dates.format(Dates.now(), "HH:MM")
            if prevNow == now
                println("$temp / $frameCnt")
            else
                println("$temp / $frameCnt,  $now")
            end
            prevNow = now
            # Velocities at mesh points
            vorts = [IFR_field.tev;IFR_field.lev]
            global mesh = meshgrid(surf,vorts,.25,t,100,view)
            vorts = [vorts ; IFR_surf.bv]
            for j = 1:size(vorts,1)
                u,w = ind_vel(vorts[j],mesh.x,mesh.z)
                global mesh.uMat = mesh.uMat .+ u
                global mesh.wMat = mesh.wMat .+ w
                global mesh.velMag = sqrt.(mesh.uMat.^2 + mesh.wMat.^2)
                global mesh.t = t
                global mesh.circ = circ
            end
            if frames == 0
                frames = [mesh;]
            else
                push!(frames,mesh)
            end
        end
        if i > 0
            ## Pressures (Change in pressure for each panel)
            for j = 1:length(dp)
                global dp[j] = surf.rho * ( ( (relVel[j,1] + surf.uind[j])*tau[j,1] + (relVel[j,2] + surf.wind[j])*tau[j,2] ) * circ[j]/ds[j] + circChange[j])
            end
            ## Loads
            # Cn normal
            ndem = .5 * surf.rho * surf.kinem.u^2 * surf.c # non-dimensionalization constant
            cn = sum( dp[j]*ds[j] for j = 1:length(dp) ) / ndem
            # LESP
            x = 1 - 2*ds[1] / surf.c
            surf.a0[1] = 1.13 * circ[1] / ( surf.kinem.u * surf.c * ( acos(x) + sin(acos(x) )))
            # Cs suction
            cs = 2*pi*surf.a0[1]^2
            # Cl lift
            alpha = surf.kinem.alpha
            cl = cn*cos(alpha) + cs*sin(alpha)
            # Cd drag
            cd = cn*sin(alpha) - cs*cos(alpha)
            # Cm moment
            ndem = ndem * surf.c # change ndem to moment ndem
            cm = -sum( dp[j]*ds[j]*(vor_loc[j,1]-surf.pvt) for j = 1:length(dp) ) / ndem

            mat = hcat(mat,[t, surf.kinem.alpha*180/pi, surf.kinem.h, surf.kinem.u, surf.a0[1], cl, cd, cm]) # Pressures and loads

        end
        ## Wake Rollup
        IFR_field = wakeroll(IFR_surf,IFR_field,dtstar)

        # Convert IFR_field values to curfield (IFR -> BFR)
        push!(curfield.tev,TwoDVort(0,0,IFR_field.tev[end].s,0.02*surf.c,0.,0.))
        if size(IFR_field.lev,1) > 0 && surf.levflag[1] == 1
            push!(curfield.lev,TwoDVort(0,0,IFR_field.lev[end].s,0.02*surf.c,0.,0.))
        end
        vor_BFR(surf,t+dtstar,IFR_field,curfield) # Calculate IFR positions for next iteration

    end
    mat = mat'
    test1 = TEVhist
    test2 = LEVhist
    return frames,IFR_field,mat,test1,test2
end
