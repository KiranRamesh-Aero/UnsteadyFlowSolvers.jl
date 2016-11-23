function lautat(surf::TwoDSurf, curfield::TwoDFlowField, nsteps::Int64 = 500, dtstar::Float64 = 0.015, delvort = DelVortDef(0, 0, 0.), mat = Array(Float64, 0, 8), kelv_enf = 0.)

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

        cl, cd, cm = calc_forces(surf)

        mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, surf.a0[1], cl, cd, cm])
    end
    mat = mat'

    mat, surf, curfield, kelv_enf

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


function ldvm(surf::TwoDSurf, curfield::TwoDFlowField, nsteps::Int64 = 500, dtstar::Float64 = 0.015, delvort = DelVortDef(0, 0, 0.), mat = Array(Float64, 0, 8), kelv_enf = 0.)

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
        cl, cd, cm = calc_forces(surf)
        #bnd_circ = (surf.a0[1] + surf.aterm[1]/2.)
        #mat[istep,:] = [t surf.kinem.alpha surf.kinem.h surf.kinem.u surf.a0[1] cl cd cm bnd_circ cn cs]
        mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, surf.a0[1], cl, cd, cm])
    end

    mat = mat'
    mat, surf, curfield, kelv_enf
end

function ldvm_more(surf::TwoDSurf, curfield::TwoDFlowField, nsteps::Int64 = 500, dtstar::Float64 = 0.015, delvort = DelVortDef(0, 0, 0.), mat = Array(Float64, 0, 11), kelv_enf = 0.)

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
        mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, surf.a0[1], cl, cd, cm])
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

function QSLLT_lautat(surf :: ThreeDSurf, field :: ThreeDFlowField, nsteps :: Int64, dtstar :: Float64)

    psi = zeros(surf.nspan)
    dpsi = pi/surf.nspan
    
    for i = 1:surf.nspan
        psi[i] = (real(i)-0.5)*dpsi
    end

    mat = Array{Float64, 2}[] 
    
    if surf.kindef.vartype == "Constant"
        
        #This is just the 2D solution corrected to 3D - kinematics at all strips is the same

        kinem2d = KinemDef(surf.kindef.alpha, surf.kindef.h, surf.kindef.u)
        
        surf2d = TwoDSurf(surf.patchdata[1].coord_file, surf.patchdata[1].pvt, kinem2d, [surf.patchdata[1].lc;])

        #If 3D flow field is defined with disturbances or external vortices, these should be transferred to the 2D flowfield
        curfield2d = TwoDFlowField()
        
        mat2d, surf2d, curfield2d = lautat_wakeroll_more(surf2d, curfield2d, nsteps, dtstar)

        #Fill the 2D solution into the 3D data structures

        for i = 1:surf.nspan
            push!(mat, mat2d)
            surf.cam[i,:] = surf2d.cam
            surf.cam_slope[i,:] = surf2d.cam_slope
            surf.kinem[i,:] = surf2d.kinem
            surf.bnd_x[i,:] = surf2d.bnd_x
            surf.bnd_z[i,:] = surf2d.bnd_z
            surf.uind[i,:] = surf2d.uind
            surf.vind[i,:] = 0
            surf.wind[i,:] = surf2d.wind
            surf.downwash[i,:] = surf2d.downwash
            surf.a0[i] = surf2d.a0[1]
            surf.aterm[i,:] = surf2d.aterm
            surf.a0dot[i] = surf2d.a0dot[1]
            surf.adot[i,:] = surf2d.adot
            surf.a0prev[i] = surf2d.a0prev[1]
            surf.aprev[i,:] = surf2d.aprev
            surf.levflag[i] = surf2d.levflag[1]

            #Convert the bound vortices into 3d vortices and add them
            for j = 1:length(surf2d.bv)
                surf.bv[i,j] = ThreeDVort([surf2d.bv[j].x; surf.yle[i]; surf2d.bv[j].z], [0.; surf2d.bv[j].s; 0.], surf2d.bv[j].vc, surf2d.bv[j].vx, 0., surf2d.bv[j].vz)
            end
        end
        
        for j = 1:surf.nspan
            for i = 1:length(curfield2d.lev)
                push!(field.lev, ThreeDVort([curfield2d.lev[i].x; surf.yle[j]; curfield2d.lev[i].z], [0.; curfield2d.lev[i].s; 0.], curfield2d.lev[i].vc, curfield2d.lev[i].vx, 0., curfield2d.lev[i].vz))
            end
            for i = 1:length(curfield2d.tev)
                push!(field.tev, ThreeDVort([curfield2d.tev[i].x; surf.yle[j]; curfield2d.tev[i].z], [0.; curfield2d.tev[i].s; 0.], curfield2d.tev[i].vc, curfield2d.tev[i].vx, 0., curfield2d.tev[i].vz))
            end
            for i = 1:length(curfield2d.extv)
                push!(field.extv, ThreeDVort([curfield2d.extv[i].x; surf.yle[j]; curfield2d.extv[i].z], [0.; curfield2d.extv[i].s; 0.], curfield2d.extv[i].vc, curfield2d.extv[i].vx, 0., curfield2d.extv[i].vz))
            end
        end

        AR = surf.bref/surf.cref
        
        lhs = zeros(surf.nspan,surf.nbterm)
        rhs = zeros(surf.nspan)
        bcoeff = zeros(nsteps,surf.nbterm)
        
        cnc_finite = zeros(nsteps)
        cnnc_finite = zeros(nsteps)

        # There is no apparent mass correction in this method

        for i = 1:nsteps
            for j = 1:surf.nspan
                for n = 1:surf.nbterm
                    lhs[j,n] = sin(n*psi[j])*(sin(psi[j]) + (n*pi/(2*AR)))
                end
                rhs[j] = pi*sin(psi[j])*mat[j][i,9]/(2*AR)
            end
            
            bcoeff[i,:] = \(lhs, rhs)
        end
        
        a03d = zeros(nsteps,surf.nspan)
        cd_ind = zeros(nsteps)

        for i = 1:nsteps
            cd_ind[i] = 0
            for n = 1:surf.nbterm
                cd_ind[i] = cd_ind[i] + real(n)*bcoeff[i,n]^2
            end
            cd_ind[i] = cd_ind[i]*pi*AR
            for j = 1:surf.nspan
                a03d[i,j] = 0

                for n = 1:surf.nbterm
                    a03d[i,j] = a03d[i,j] - real(n)*bcoeff[i,n]*sin(n*psi[j])/sin(psi[j])
                end
            end
        end

        cn = zeros(nsteps)
        cs = zeros(nsteps)
        cl = zeros(nsteps)
        cd = zeros(nsteps)
        cn3d = zeros(surf.nspan)
        cs3d = zeros(surf.nspan)
        cl3d = zeros(surf.nspan)
        cd3d = zeros(surf.nspan)
        
        for i = 1:nsteps
            cn[i] = 0
            cs[i] = 0
            update_kinem(surf2d, mat[1][i,1])
                        
            for j = 1:surf.nspan
                cn3d[j] = mat[j][i,10] + (2*pi/surf.uref)*(mat[j][i,4]*cos(mat[j][i,2]) + surf2d.kinem.hdot*sin(mat[j][i,2]))*a03d[i,j]
                cs3d[j] = mat[j][i,11] + 2*pi*a03d[i,j]^2
                cl3d[j] = cn3d[j]*cos(mat[j][i,2]) + cs3d[j]*sin(mat[j][i,2])
                cd3d[j] = cn3d[j]*sin(mat[j][i,2]) - cs3d[j]*cos(mat[j][i,2]) 
                cn[i] = cn[i] + cn3d[j]*sin(psi[j])*dpsi/2
                cs[i] = cs[i] + cs3d[j]*sin(psi[j])*dpsi/2
                cl[i] = cl[i] + cl3d[j]*sin(psi[j])*dpsi/2
                cd[i] = cd[i] + cd3d[j]*sin(psi[j])*dpsi/2
            end
        end
        return cl, cd, cd_ind, surf, field, mat2d, a03d
        
    end

    
end

function QSLLT_ldvm(surf :: ThreeDSurf, field :: ThreeDFlowField, nsteps :: Int64, dtstar :: Float64)

    mat = Array{Float64, 2}[] 
    
    if surf.kindef.vartype == "Constant"
        
        # Kinematics at all strips is the same

        kinem2d = KinemDef(surf.kindef.alpha, surf.kindef.h, surf.kindef.u)
        
        surf2d = TwoDSurf(surf.patchdata[1].coord_file, surf.patchdata[1].pvt, kinem2d, [surf.patchdata[1].lc;])

        #If 3D flow field is defined with disturbances or external vortices, these should be transferred to the 2D flowfield
        curfield2d = TwoDFlowField()
        
        mat2d, surf2d, curfield2d = lautat_wakeroll_more(surf2d, curfield2d, nsteps, dtstar)

        #Fill the 2D solution into the 3D data structures

        for i = 1:surf.nspan
            push!(mat, mat2d)
            surf.cam[i,:] = surf2d.cam
            surf.cam_slope[i,:] = surf2d.cam_slope
            surf.kinem[i,:] = surf2d.kinem
            surf.bnd_x[i,:] = surf2d.bnd_x
            surf.bnd_z[i,:] = surf2d.bnd_z
            surf.uind[i,:] = surf2d.uind
            surf.vind[i,:] = 0
            surf.wind[i,:] = surf2d.wind
            surf.downwash[i,:] = surf2d.downwash
            surf.a0[i] = surf2d.a0[1]
            surf.aterm[i,:] = surf2d.aterm
            surf.a0dot[i] = surf2d.a0dot[1]
            surf.adot[i,:] = surf2d.adot
            surf.a0prev[i] = surf2d.a0prev[1]
            surf.aprev[i,:] = surf2d.aprev
            surf.levflag[i] = surf2d.levflag[1]

            #Convert the bound vortices into 3d vortices and add them
            for j = 1:length(surf2d.bv)
                surf.bv[i,j] = ThreeDVort([surf2d.bv[j].x; surf.yle[i]; surf2d.bv[j].z], [0.; surf2d.bv[j].s; 0.], surf2d.bv[j].vc, surf2d.bv[j].vx, 0., surf2d.bv[j].vz)
            end
        end
        
        for j = 1:surf.nspan
            for i = 1:length(curfield2d.lev)
                push!(field.lev, ThreeDVort([curfield2d.lev[i].x; surf.yle[j]; curfield2d.lev[i].z], [0.; curfield2d.lev[i].s; 0.], curfield2d.lev[i].vc, curfield2d.lev[i].vx, 0., curfield2d.lev[i].vz))
            end
            for i = 1:length(curfield2d.tev)
                push!(field.tev, ThreeDVort([curfield2d.tev[i].x; surf.yle[j]; curfield2d.tev[i].z], [0.; curfield2d.tev[i].s; 0.], curfield2d.tev[i].vc, curfield2d.tev[i].vx, 0., curfield2d.tev[i].vz))
            end
            for i = 1:length(curfield2d.extv)
                push!(field.extv, ThreeDVort([curfield2d.extv[i].x; surf.yle[j]; curfield2d.extv[i].z], [0.; curfield2d.extv[i].s; 0.], curfield2d.extv[i].vc, curfield2d.extv[i].vx, 0., curfield2d.extv[i].vz))
            end
        end

        AR = surf.bref/surf.cref
        
        lhs = zeros(surf.nspan,surf.nbterm)
        rhs = zeros(surf.nspan)
        bcoeff = zeros(nsteps,surf.nbterm)
        
        cnc_finite = zeros(nsteps)
        cnnc_finite = zeros(nsteps)

        # There is no apparent mass correction in this method

        for i = 1:nsteps
            for j = 1:surf.nspan
                for n = 1:surf.nbterm
                    lhs[j,n] = sin(n*surf.psi[j])*(sin(surf.psi[j]) + (n*pi/(2*AR)))
                end
                rhs[j] = pi*sin(surf.psi[j])*mat[j][i,9]/(2*AR)
            end
            
            bcoeff[i,:] = \(lhs, rhs)
        end
        
        a03d = zeros(nsteps,surf.nspan)
        cd_ind = zeros(nsteps)

        for i = 1:nsteps
            cd_ind[i] = 0
            for n = 1:surf.nbterm
                cd_ind[i] = cd_ind[i] + real(n)*bcoeff[i,n]^2
            end
            cd_ind[i] = cd_ind[i]*pi*AR
            for j = 1:surf.nspan
                a03d[i,j] = 0

                for n = 1:surf.nbterm
                    a03d[i,j] = a03d[i,j] - real(n)*bcoeff[i,n]*sin(n*surf.psi[j])/sin(surf.psi[j])
                end
            end
        end

        cn = zeros(nsteps)
        cs = zeros(nsteps)
        cl = zeros(nsteps)
        cd = zeros(nsteps)
        cn3d = zeros(surf.nspan)
        cs3d = zeros(surf.nspan)
        cl3d = zeros(surf.nspan)
        cd3d = zeros(surf.nspan)
        
        for i = 1:nsteps
            cn[i] = 0
            cs[i] = 0
            update_kinem(surf2d, mat[1][i,1])
                        
            for j = 1:surf.nspan
                cn3d[j] = mat[j][i,10] + (2*pi/surf.uref)*(mat[j][i,4]*cos(mat[j][i,2]) + surf2d.kinem.hdot*sin(mat[j][i,2]))*a03d[i,j]
                cs3d[j] = mat[j][i,11] + 2*pi*a03d[i,j]^2
                cl3d[j] = cn3d[j]*cos(mat[j][i,2]) + cs3d[j]*sin(mat[j][i,2])
                cd3d[j] = cn3d[j]*sin(mat[j][i,2]) - cs3d[j]*cos(mat[j][i,2]) 
            end
            for j = 1:surf.nspan-1
                cn[i] = cn[i] + 0.5*(cn3d[j] + cn3d[j+1])*sin(0.5*(surf.psi[j] + surf.psi[j+1]))*(surf.psi[j+1] - surf.psi[j])/2
                cs[i] = cs[i] + 0.5*(cs3d[j] + cs3d[j+1])*sin(0.5*(surf.psi[j] + surf.psi[j+1]))*(surf.psi[j+1] - surf.psi[j])/2
                cl[i] = cl[i] + 0.5*(cl3d[j] + cl3d[j+1])*sin(0.5*(surf.psi[j] + surf.psi[j+1]))*(surf.psi[j+1] - surf.psi[j])/2
                cd[i] = cd[i] + 0.5*(cd3d[j] + cd3d[j+1])*sin(0.5*(surf.psi[j] + surf.psi[j+1]))*(surf.psi[j+1] - surf.psi[j])/2
            end
        end
        return cl, cd, cd_ind, surf, field, mat2d, a03d
        
    end
    
end

function LLT_ldvm(surf :: ThreeDSurf, field :: ThreeDFlowField, nsteps :: Int64, dtstar :: Float64)
    
    mat = Array(Float64, 0, 4)

    mat = mat'
    
    surf2d = TwoDSurf[]
    field2d = TwoDFlowField[]
    kinem2d = KinemDef[]

    dt = dtstar*surf.cref/surf.uref

    t = 0.

    AR = surf.bref/surf.cref

    bc = zeros(surf.nspan)
    a03d = zeros(surf.nspan)
    cl = zeros(surf.nspan)
    cd = zeros(surf.nspan)
    cm = zeros(surf.nspan)

    lhs = zeros(surf.nspan, surf.nbterm)
    rhs = zeros(surf.nspan)
    bcoeff = zeros(surf.nbterm)
    
    if surf.kindef.vartype == "Constant"

        for i = 1:surf.nspan
            # Kinematics at all strips is the same
            
            push!(kinem2d, KinemDef(surf.kindef.alpha, surf.kindef.h, surf.kindef.u))
            push!(surf2d, TwoDSurf(surf.patchdata[1].coord_file, surf.patchdata[1].pvt, kinem2d[i], [surf.patchdata[1].lc;]))
            #If 3D flow field is defined with disturbances or external vortices, these should be transferred to the 2D flowfield
            push!(field2d, TwoDFlowField())
        end
    end

    for istep = 1:nsteps
        #Udpate current time
        t = t + dt

        for i = 1:surf.nspan
            #Update kinematic parameters
            update_kinem(surf2d[i], t)
        
            #Update flow field parameters if any
            update_externalvel(field2d[i], t)

            #Update bound vortex positions
            update_boundpos(surf2d[i], dt)

            #Add a TEV with dummy strength
            place_tev(surf2d[i], field2d[i], dt)
        end
        
        kelv = KelvinConditionLLTldvm(surf, surf2d, field2d)
        
        #Solve for TEV strength to satisfy Kelvin condition
        
        soln = nlsolve(not_in_place(kelv), -0.01*ones(surf.nspan))
        
        for i = 1:surf.nspan
            field2d[i].tev[length(field2d[i].tev)].s = soln.zero[i]

            #Update incduced velocities on airfoil
            update_indbound(surf2d[i], field2d[i])
            
            #Calculate downwash
            update_downwash(surf2d[i], [field2d[i].u[1],field2d[i].w[1]])
            
            #Calculate first two fourier coefficients
            update_a0anda1(surf2d[i])

            bc[i] = surf2d[i].a0[1] + 0.5*surf2d[i].aterm[1]
        end

        for i = 1:surf.nspan
            for n = 1:surf.nbterm
                lhs[i,n] = sin(n*surf.psi[i])*(sin(surf.psi[i]) + (n*pi/(2*AR)))
            end
            rhs[i] = pi*sin(surf.psi[i])*bc[i]/(2*AR)
        end
        bcoeff[:] = \(lhs, rhs)   

        for i = 1:surf.nspan
            a03d[i] = 0
            for n = 1:surf.nbterm
            a03d[i] = a03d[i] - real(n)*bcoeff[n]*sin(n*surf.psi[i])/sin(surf.psi[i])
            end
        end

        for i = 1:surf.nspan
            #Update 3D effect on A0
            surf2d[i].a0[1] = surf2d[i].a0[1] + a03d[i]
                
            #Update rest of Fourier terms
            update_a2toan(surf2d[i])
            
            #Update derivatives of Fourier coefficients
            update_adot(surf2d[i],dt)
            
            #Set previous values of aterm to be used for derivatives in next time step
            surf2d[i].a0prev[1] = surf2d[i].a0[1]
            for ia = 1:3
                surf2d[i].aprev[ia] = surf2d[i].aterm[ia]
            end
            
            #Calculate bound vortex strengths
            update_bv(surf2d[i])
            
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
            wakeroll(surf2d[i], field2d[i], dt)

            if (surf.levflag[i] == 1) 
                cl[i], cd[i], cm[i] = calc_forces_E(surf2d[i], field2d[i].lev[length(field2d[i].lev)].s, dt)
            else
                cl[i], cd[i], cm[i] = calc_forces(surf2d[i])
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
    mat, surf2d, field2d
    
end
        
        
           
