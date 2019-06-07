#=
LVE.jl
    Written by Matthew White
    5/21/2019
=#
function LVE(surf::TwoDSurf,curfield::TwoDFlowField,nsteps::Int64 = 500,dtstar::Float64 = 0.015,frameCnt = 20)
    ## Initialize variables
    prevCirc = 0
    prevWake = 0
    vor_loc = zeros(surf.ndiv-1,2) # (x1,z1 ; x2,z2 ...) Inertial frame
    coll_loc = zeros(surf.ndiv-1,2) # (x1,z1 ; x2,z2 ...)
    global a = zeros(surf.ndiv,surf.ndiv)
    a[end,:] .= 1. # for wake portion
    global RHS = zeros(surf.ndiv)
    glo_vor = 0 .* vor_loc # Global frame
    glo_field = TwoDFlowField()
    aniStep = round(nsteps / frameCnt) # Frames to capture animation on
    aniStep < 1 ? aniStep = 1 : aniStep
    frames = 0

    # Normal vectors for each panel
    n = hcat(sin.(surf.cam_slope),cos.(surf.cam_slope))
    # Tangential vectors for each panel
    tau = hcat(cos.(surf.cam_slope),-sin.(surf.cam_slope))

    ## Vortex Locations
    # Located at 25% of each panel
    #   length of each panel
    ds = sqrt.(( surf.x[1:end-1]-surf.x[2:end]).^2 + (surf.cam[1:end-1]-surf.cam[2:end]).^2)
    #   Surf cam_slope does not give the correct panel locations
    cam_slope = asind.((surf.cam[2:end]-surf.cam[1:end-1])./ds)
    #   Vortex loactions (at 25% of panel length)
    vor_loc[:,1] = surf.x[1:end-1] + .25 .* ds .* cosd.(cam_slope)
    vor_loc[:,2] = surf.cam[1:end-1] + .25 .* ds .* sind.(cam_slope)
    #   Collocation points (at 75% panel chordwise)
    coll_loc[:,1] = surf.x[1:end-1] + .75 .* ds .* cosd.(cam_slope)
    coll_loc[:,2] = surf.cam[1:end-1] + .75 .* ds .* sind.(cam_slope)
    vcore = 1.3 * dtstar * surf.c # From Dr. Narisipur's paper

    ## Time Stepping Loop
    t = 0

    for i = 1:nsteps
        t += dtstar

        # Update Kinematics
        update_kinem(surf, t)

        ## Refresh vortex values
        refresh_vortex(surf,vor_loc,true)

        ## Influence Coeffcients
        for j = 1:surf.ndiv-1, k = 1:surf.ndiv-1
            # j is test location, k is vortex source
            u,w = ind_vel(surf.bv[k],coll_loc[j,1],coll_loc[j,2])
            global a[j,k] = u*n[j,1] + w*n[j,2]
        end

        ## RHS Vector
        RHS[end] = prevCirc # previous total circulation
        for j = 1:surf.ndiv-1
            theta = surf.kinem.alpha
            # relVel = [U_t ; W_t]
            relVel = [cosd(theta) sind(theta) ; -sind(theta) cosd(theta)] * [-surf.kinem.u ; -surf.kinem.hdot]
            #U = -surf.kinem.u
            #W = -surf.kinem.hdot
            u_w = 0
            w_w = 0
            if length(curfield.tev) > 0 # Trailing edge vortexs exist
                u_w,w_w = ind_vel(curfield.tev,coll_loc[j,1],coll_loc[j,2])
                u_w = u_w[1]
                w_w = w_w[1]
            end
            # RHS = -[U_t + u_w , W_t + w_w] dot n
            global RHS[j] = -(relVel[1] + u_w)*n[j,1] - (relVel[2] + w_w)*n[j,2]
        end

        ## Vortex Strenghths a*circ = RHS
        global circ = a\RHS
        prevCirc = sum(circ[1:end-1])
        for j = 1:surf.ndiv-1
            surf.bv[j].s = circ[j]
        end
        ## TEV setup
        if i > 1 # Create first vortex after intial circulation calculation
            place_tev2(surf,curfield,dtstar)
            curfield.tev[end].s = circ[end] - prevWake
            prevWake = sum(map(q -> q.s, curfield.tev))
        end

        ## Velocity Components
        # Velocities at collocation points
        surf.uind .= surf.kinem.u
        surf.wind .= surf.kinem.hdot
        vorts = [surf.bv ; curfield.tev[:]]
        for j = 1:surf.ndiv-1, k = 1:size(vorts,1)
            u,w = ind_vel(vorts[k],coll_loc[j,1],coll_loc[j,2])
            surf.uind[j] += u
            surf.wind[j] += w
        end

        # Gloabal reference frame
        glo_vor[:,1],glo_vor[:,2] = globalFrame(surf,vor_loc[:,1],vor_loc[:,2],t)
        refresh_vortex(surf,glo_vor,false)
        if size(curfield.tev,1) > 0
            tev_X,tev_Z = globalFrame(surf,curfield.tev[end].x,curfield.tev[end].z,t)
            push!(glo_field.tev,TwoDVort(tev_X,tev_Z,curfield.tev[end].s,0.02*surf.c,0.,0.))
        end
        # position mesh
        if i % aniStep == 0.
            view = "square"
            global mesh = meshgrid(surf,.25,t,100,view)
            # Velocities at mesh points
            vorts = [surf.bv ; glo_field.tev[:]]
            for j = 1:size(vorts,1)
                u,w = ind_vel(vorts[j],mesh.x,mesh.z)
                global mesh.uMat = mesh.uMat .+ u
                global mesh.wMat = mesh.wMat .+ w
                global mesh.velMag = sqrt.(mesh.uMat.^2 + mesh.wMat.^2)
            end
            if frames == 0
                frames = [mesh;]
            else
                push!(frames,mesh)
            end
            #meshPlot(mesh,200,:kdc) Not working
        end
        ## Pressures

        ## Loads


        ## Wake Rollup

        global newsurf = surf
    end
    test = 0
    return newsurf,frames,curfield,glo_field,test
end
