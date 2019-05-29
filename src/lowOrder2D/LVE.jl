#=
LVE.jl
    Written by Matthew White
    5/21/2019
=#
function LVE(surf::TwoDSurf,curfield::TwoDFlowField,nsteps::Int64 = 500,dtstar::Float64 = 0.015)
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
    vor_loc = zeros(surf.ndiv-1,2) # (x1,z1 ; x2,z2 ...)
    vor_loc[:,1] = surf.x[1:end-1] + .25 .* ds .* cosd.(cam_slope)
    vor_loc[:,2] = surf.cam[1:end-1] + .25 .* ds .* sind.(cam_slope)
    #   Collocation points (at 75% panel chordwise)
    coll_loc = zeros(surf.ndiv-1,2) # (x1,z1 ; x2,z2 ...)
    coll_loc[:,1] = surf.x[1:end-1] + .75 .* ds .* cosd.(cam_slope)
    coll_loc[:,2] = surf.cam[1:end-1] + .75 .* ds .* sind.(cam_slope)
    vcore = 1.3 * dtstar * surf.c # From Dr. Narisipur's paper
    # initialize influence matrix
    global a = zeros(surf.ndiv,surf.ndiv)
    a[end,:] .= 1. # for wake portion

    ## Time Stepping Loop
    t = 0

    for i = 1:1#nsteps
        t += dtstar

        # Update Kinematics
        update_kinem(surf, t)
        ## Gloabal reference frame
        X,Z = globalFrame(surf,surf.x,surf.cam,t)
        vorX,vorZ = globalFrame(surf,vor_loc[:,1],vor_loc[:,2],t)
        ## Refresh vortex values
        refresh_vortex(surf,vor_loc,true)

        ## Influence Coeffcients
        for j = 1:surf.ndiv-1, k = 1:surf.ndiv-1
            # j is test location, k is vortex source
            u,w = ind_vel(surf.bv[k],coll_loc[j,1],coll_loc[j,2])
            global a[j,k] = u*n[j,1] + w*n[j,2]
        end

        ## RHS Vector
        global RHS = zeros(surf.ndiv)
        RHS[end] = 0 # previous total circulation
        for j = 1:surf.ndiv-1
            U = surf.kinem.u
            W = surf.kinem.hdot
            if length(curfield.tev) == 0 # No trailing edge vortexs
                u_w = 0
                w_w = 0
            else
                u_w,w_w = ind_vel(curfield.tev,coll_loc[j,1],coll_loc[j,2])
            end
            # RHS = -[U+u_w , W+w_w] dot n
            global RHS[j] = -(U + u_w)*n[j,1] - (W + w_w)*n[j,2]
        end

        ## Vortex Strenghths a*circ = RHS
        global circ = a\RHS

        ## Velocity Components

        ## Pressures

        ## Loads


        ## Wake Rollup
        global newsurf = surf
    end
    test = circ
    newsurf,test
    end
