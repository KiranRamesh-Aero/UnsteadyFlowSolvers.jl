#=
LVE.jl
    Written by Matthew White
    5/21/2019
=#
function LVE(surf::TwoDSurf,curfield::TwoDFlowField,nsteps::Int64 = 500,dtstar::Float64 = 0.015,frameCnt = 20,view = "square")
    ## Initialize variables
    mat = zeros(8 ,0) # loads output matrix
    prevCirc = 0
    vor_loc = zeros(surf.ndiv-1,2) # (x1,z1 ; x2,z2 ...) Inertial frame
    coll_loc = zeros(surf.ndiv-1,2) # (x1,z1 ; x2,z2 ...)
    global RHS = zeros(surf.ndiv)
    IFR_vor = 0 .* vor_loc # Global frame
    IFR_field = TwoDFlowField()
    dp = zeros(surf.ndiv-1,1) # change in panel pressure
    aniStep = round((nsteps-1) / frameCnt) # Frames to capture animation on
    aniStep < 1 ? aniStep = 1 : aniStep
    frames = 0
    prevNow = 0

    #   length of each panel
    ds = sqrt.(( surf.x[1:end-1]-surf.x[2:end]).^2 + (surf.cam[1:end-1]-surf.cam[2:end]).^2)
    #   Surf cam_slope does not give the correct panel locations
    cam_slope = asind.((surf.cam[2:end]-surf.cam[1:end-1])./ds)

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
    vcore = 1.3 * dtstar * surf.c # From Dr. Narisipur's paper

    ## Time Stepping Loop
    t = 0

    for i = 0:nsteps
        t += dtstar

        # Update Kinematics
        update_kinem(surf, t)

        ## Refresh vortex values
        refresh_vortex(surf,vor_loc)

        ## Influence Coeffcients
        global a = influence_coeff(surf,curfield,coll_loc,n)

        ## RHS Vector
        RHS[end] = prevCirc # previous total circulation
        for j = 1:surf.ndiv-1
            theta = surf.kinem.alpha
            # relVel = [U_t ; W_t]
            global relVel = [cos(theta) -sin(theta) ; sin(theta) cos(theta)] * [surf.kinem.u ; -surf.kinem.hdot] + [-surf.kinem.alphadot*coll_loc[j,2] ; surf.kinem.alphadot * (coll_loc[j,1] - surf.pvt)]
            u_w = 0
            w_w = 0
            if length(curfield.tev) > 0 # Trailing edge vortexs exist
                u_w,w_w = ind_vel(curfield.tev,coll_loc[j,1],coll_loc[j,2])
                global u_w = u_w[1]
                global w_w = w_w[1]
            end
            # RHS = -[U_t + u_w , W_t + w_w] dot n
            global RHS[j] = -(relVel[1] + u_w)*n[j,1] - (relVel[2] + w_w)*n[j,2]
        end

        ## Vortex Strenghths a*circ = RHS
        global circ = a\RHS
        circChange = sum(circ[1:end-1]) - prevCirc #bound vorticy circ change
        prevCirc = sum(circ[1:end-1])
        for j = 1:surf.ndiv-1 # Set bv circs
            surf.bv[j].s = circ[j]
        end
        # Inertial reference frame
        IFR_vor[:,1],IFR_vor[:,2] = IFR(surf,vor_loc[:,1],vor_loc[:,2],t)
        IFR_surf = refresh_vortex(surf,IFR_vor)

        ## TEV setup
        place_tev2(IFR_surf,IFR_field,dtstar)
        if length(IFR_field.tev) == 1
            IFR_field.tev[1].s = 0
        else
            IFR_field.tev[end].s = circ[end]
        end

        # Position mesh
        if i % aniStep == 0.

            temp = Int(i / aniStep)
            now = Dates.format(Dates.now(), "HH:MM")
            if prevNow == now
                println("$temp / $frameCnt")
            else
                println("$temp / $frameCnt,  $now")
            end
            prevNow = now
            # Velocities at mesh points
            vorts = IFR_field.tev[:]
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
                dp[j] = surf.rho * ( ( (relVel[1] + u_w)*tau[j,1] + (relVel[2] + w_w)*tau[j,2] ) * circ[j]/ds[j] + circChange )
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

            mat = hcat(mat,[t, surf.kinem.alpha*180/pi, surf.kinem.h, surf.kinem.u, surf.a0[1], cl, cd, cm])
        end
        ## Wake Rollup
        IFR_field = wakeroll(IFR_surf,IFR_field,dtstar)

        # Convert IFR_field values to curfield (IFR to BFR)
        push!(curfield.tev,TwoDVort(0,0,IFR_field.tev[end].s,0.02*surf.c,0.,0.))
        vor_BFR(surf,dtstar,IFR_field,curfield)
    end
    mat = mat'
    test = ds
    return frames,IFR_field,mat,test
end
