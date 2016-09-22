function lautat(surf::TwoDSurf, curfield::TwoDFlowField, nsteps::Int64)
    outfile = open("results.dat", "w")

    dtstar = 0.015
    dt = dtstar*surf.c/surf.uref
    t = 0.

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

        #wakeroll(surf, curfield)

        cl, cd, cm = calc_forces(surf)
        write(outfile, join((t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, surf.a0[1], cl, cd, cm)," "), "\n")

    end
    close(outfile)

    #Plot flowfield viz and A0 history
    figure(0)
    view_vorts(surf, curfield)

end

function lautat_wakeroll(surf::TwoDSurf, curfield::TwoDFlowField, nsteps::Int64)
    outfile = open("results.dat", "w")

    dtstar = 0.015
    dt = dtstar*surf.c/surf.uref
    nsteps = 500
    t = 0.

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

        wakeroll(surf, curfield, dt)

        cl, cd, cm = calc_forces(surf)
        write(outfile, join((t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, surf.a0[1], cl, cd, cm)," "), "\n")

    end
    close(outfile)

    #Plot flowfield viz
    figure(0)
    view_vorts(surf, curfield)

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


function ldvm(surf::TwoDSurf, curfield::TwoDFlowField, nsteps::Int64 = 500, dtstar::Float64 = 0.015, delvort = DelVortDef(0, 0, 0.))
    mat = zeros(nsteps,11)

    dt = dtstar*surf.c/surf.uref
    t = 0.

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

        cl, cd, cm, cn, cs = calc_forces(surf)
        bnd_circ = surf.uref*surf.c*pi*(surf.a0[1] + surf.aterm[1]/2.)
        mat[istep,:] = [t surf.kinem.alpha surf.kinem.h surf.kinem.u surf.a0[1] cl cd cm bnd_circ cn cs]
    end



    mat, surf, curfield
    #Plot flowfield viz
#    figure(0)
#    view_vorts(surf, curfield)

end


function ldvm(surf::TwoDSurfwFlap, curfield::TwoDFlowField, nsteps::Int64 = 500, dtstar::Float64 = 0.015)
    mat = zeros(nsteps,9)

    dt = dtstar*surf.c/surf.uref
    t = 0.

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

        wakeroll(surf, curfield, dt)

        cl, cd, cm, cm_be = calc_forces(surf, dt)

        mat[istep,:] = [t surf.kinem.alpha surf.kinem.h surf.kinem.u surf.a0[1] cl cd cm cm_be]

    end

    return mat, surf, curfield
    #Plot flowfield viz
#    figure(0)
#    view_vorts(surf, curfield)

end

function ldvm(surf::TwoDSurf_2DOF, curfield::TwoDFlowField, nsteps::Int64 = 500, dtstar::Float64 = 0.015, delvort = DelVortDef(0, 0, 0.))
    mat = zeros(nsteps,8)

    dt = dtstar*surf.c/surf.uref
    t = 0.
    kelv_enf = 0

    #Intialise flowfield
    for istep = 1:nsteps
        #Udpate current time
        t = t + dt

        #Update kinematic parameters (based on 2DOF response)
        if (istep > 1) # Allow initial condition
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
        cl, cd, cm = calc_forces(surf)

        #Using the force data, update - hddot and alphaddot
        calc_struct2DOF(surf, cl, cm)
        mat[istep,:] = [t surf.kinem.alpha surf.kinem.h surf.kinem.u surf.a0[1] cl cd cm]
    end

    mat, surf, curfield
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
