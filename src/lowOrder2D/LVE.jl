#=
LVE.jl
    Written by Matthew White
    5/21/2019
=#
function LVE(surf::TwoDSurf,curfield::TwoDFlowField,nsteps::Int64 = 500,dtstar::Float64 = 0.015)
    # Normal vectors for each panel
    n = hcat(-sin.(surf.cam_slope),cos.(surf.cam_slope))
    # Tangential vectors for each panel
    tau = hcat(cos.(surf.cam_slope),sin.(surf.cam_slope))

    ## Vortex Locations
    # Located at 25% of each panel
    #   length of each panel
    ds = sqrt.(( surf.x[1:end-1]-surf.x[2:end]).^2 + (surf.cam[1:end-1]-surf.cam[2:end]).^2)
    #   Surf cam_slope does not give the correct panel locations
    cam_slope = asind.((surf.cam[2:end]-surf.cam[1:end-1])./ds)
    #   Vortex loactions
    vor_loc = zeros(surf.ndiv-1,2) # (x1,z1 ; x2,z2 ...)
    vor_loc[:,1] = surf.x[1:end-1] + .25 .* ds .* cosd.(cam_slope)
    vor_loc[:,2] = surf.cam[1:end-1] + .25 .* ds .* sind.(cam_slope)
    vcore = 1.3 * dtstar # From Dr. Narisipur's paper

    ## Time Stepping Loop
    t = 0

    for i = 1:1#nsteps
        t += dtstar

        # Update Kinematics
        update_kinem(surf, t)
        ## Inertial Frame (stationary) currently wrong
        #X = (surf.x - pvt).*cosd(alphadef(t)) + hdef(t)*sind(alphadef(t)) + udef(t)*t

        ## Influence Coeffcients
        a = ones(length(vor_loc))
        for j = 1:length(vor_loc), k = 1:length(vor_loc)
            # j is test location, k is vortex source
            #a[j,k]
        end
        ## RHS Vector


        ## Vortex Strenghths


        ## Velocity Components

        ## Pressures

        ## Loads


        ## Wake Rollup

    end
    test = vor_loc
    surf,test
    end
