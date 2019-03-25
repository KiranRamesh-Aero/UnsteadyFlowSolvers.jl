function update_atermdot(surf::TwoDSurfThick,dt)
    surf.a0dot[1] = (surf.a0[1] - surf.a0prev[1])/dt
    for ia = 1:surf.naterm
        surf.adot[ia] = (surf.aterm[ia]-surf.aprev[ia])/dt
    end
    return surf
end

function update_indbound(surf::TwoDSurfThick, curfield::TwoDFlowField)
    surf.uind_u[1:surf.ndiv], surf.wind_u[1:surf.ndiv] = ind_vel([curfield.tev;
                                                                  curfield.lev; curfield.extv], surf.bnd_x_u, surf.bnd_z_u)

    surf.uind_l[1:surf.ndiv], surf.wind_l[1:surf.ndiv] = ind_vel([curfield.tev;
                                                                  curfield.lev; curfield.extv], surf.bnd_x_l, surf.bnd_z_l)

    return surf
end

function update_kinem(surf::TwoDSurfThick, t)

    # Pitch kinematics
    if (typeof(surf.kindef.alpha) == EldUpDef)
        surf.kinem.alpha = surf.kindef.alpha(t)
        surf.kinem.alphadot = ForwardDiff.derivative(surf.kindef.alpha,t)*surf.uref/surf.c
    elseif (typeof(surf.kindef.alpha) == EldUptstartDef)
        surf.kinem.alpha = surf.kindef.alpha(t)
        surf.kinem.alphadot = ForwardDiff.derivative(surf.kindef.alpha,t)*surf.uref/surf.c
    elseif (typeof(surf.kindef.alpha) == EldRampReturnDef)
        surf.kinem.alpha = surf.kindef.alpha(t)
        surf.kinem.alphadot = ForwardDiff.derivative(surf.kindef.alpha,t)*surf.uref/surf.c
    elseif (typeof(surf.kindef.alpha) == ConstDef)
        surf.kinem.alpha = surf.kindef.alpha(t)
        surf.kinem.alphadot = 0.
    elseif (typeof(surf.kindef.alpha) == SinDef)
        surf.kinem.alpha = surf.kindef.alpha(t)
        surf.kinem.alphadot = ForwardDiff.derivative(surf.kindef.alpha,t)*surf.uref/surf.c
    elseif (typeof(surf.kindef.alpha) == CosDef)
        surf.kinem.alpha = surf.kindef.alpha(t)
        surf.kinem.alphadot = ForwardDiff.derivative(surf.kindef.alpha,t)*surf.uref/surf.c
    elseif (typeof(surf.kindef.alpha) == VAWTalphaDef)
        surf.kinem.alpha = surf.kindef.alpha(t)
        surf.kinem.alphadot = ForwardDiff.derivative(surf.kindef.alpha,t)*surf.uref/surf.c
    end
    # ---------------------------------------------------------------------------------------------

    # Plunge kinematics
    if (typeof(surf.kindef.h) == EldUpDef)
        surf.kinem.h = surf.kindef.h(t)*surf.c
        surf.kinem.hdot = ForwardDiff.derivative(surf.kindef.h,t)*surf.uref
    elseif (typeof(surf.kindef.h) == EldUptstartDef)
        surf.kinem.h = surf.kindef.h(t)*surf.c
        surf.kinem.hdot = ForwardDiff.derivative(surf.kindef.h,t)*surf.uref
    elseif (typeof(surf.kindef.h) == EldUpIntDef)
        surf.kinem.h = surf.kindef.h(t)*surf.c
        surf.kinem.hdot = ForwardDiff.derivative(surf.kindef.h,t)*surf.uref
    elseif (typeof(surf.kindef.h) == EldUpInttstartDef)
        surf.kinem.h = surf.kindef.h(t)*surf.c
        surf.kinem.hdot = ForwardDiff.derivative(surf.kindef.h,t)*surf.uref
    elseif (typeof(surf.kindef.h) == EldRampReturnDef)
        surf.kinem.h = surf.kindef.h(t)*surf.c
        surf.kinem.hdot = ForwardDiff.derivative(surf.kindef.h,t)*surf.uref
    elseif (typeof(surf.kindef.h) == ConstDef)
        surf.kinem.h = surf.kindef.h(t)*surf.c
        surf.kinem.hdot = 0.
    elseif (typeof(surf.kindef.h) == SinDef)
        surf.kinem.h = surf.kindef.h(t)*surf.c
        surf.kinem.hdot = ForwardDiff.derivative(surf.kindef.h,t)*surf.uref
    elseif (typeof(surf.kindef.h) == CosDef)
        surf.kinem.h = surf.kindef.h(t)*surf.c
        surf.kinem.hdot =  ForwardDiff.derivative(surf.kindef.h,t)*surf.uref
    elseif (typeof(surf.kindef.h) == VAWThDef)
        surf.kinem.h = surf.kindef.h(t)*surf.c
        surf.kinem.hdot = ForwardDiff.derivative(surf.kindef.h,t)*surf.uref
    end
    # ---------------------------------------------------------------------------------------------

    # Forward velocity
    if (typeof(surf.kindef.u) == EldUpDef)
        surf.kinem.u = surf.kindef.u(t)*surf.uref
        surf.kinem.udot = ForwardDiff.derivative(surf.kindef.u,t)*surf.uref*surf.uref/surf.c
    elseif (typeof(surf.kindef.u) == EldRampReturnDef)
        surf.kinem.u, surf.kinem.udot = surf.kindef.u(t)
        surf.kinem.u = surf.kinem.u*surf.uref
        surf.kinem.udot = surf.kinem.udot*surf.uref*surf.uref/surf.c
    elseif (typeof(surf.kindef.u) == ConstDef)
        surf.kinem.u = surf.kindef.u(t)*surf.uref
        surf.kinem.udot = 0.
    elseif (typeof(surf.kindef.u) == SinDef)
        surf.kinem.u = surf.kindef.u(t)*surf.uref
        surf.kinem.udot = ForwardDiff.derivative(surf.kindef.u,t)*surf.uref*surf.uref/surf.c
    elseif (typeof(surf.kindef.u) == CosDef)
        surf.kinem.u = surf.kindef.u(t)*surf.uref
        surf.kinem.udot = ForwardDiff.derivative(surf.kindef.u,t)*surf.uref*surf.uref/surf.c
    elseif (typeof(surf.kindef.u) == LinearDef)
        surf.kinem.u = surf.kindef.u(t)*surf.uref
        surf.kinem.udot = ForwardDiff.derivative(surf.kindef.u,t)*surf.uref*surf.uref/surf.c
    elseif (typeof(surf.kindef.u) == VAWTuDef)
        surf.kinem.u = surf.kindef.u(t)*surf.uref
        surf.kinem.udot = ForwardDiff.derivative(surf.kindef.u,t)*surf.uref*surf.uref/surf.c
    end
    # ---------------------------------------------------------------------------------------------
    return surf
end

function wakeroll(surf::TwoDSurfThick, curfield::TwoDFlowField, dt)

    nlev = length(curfield.lev)
    ntev = length(curfield.tev)
    nextv = length(curfield.extv)

    #Clean induced velocities
    for i = 1:ntev
        curfield.tev[i].vx = 0
        curfield.tev[i].vz = 0
    end

    for i = 1:nlev
        curfield.lev[i].vx = 0
        curfield.lev[i].vz = 0
    end

    for i = 1:nextv
        curfield.extv[i].vx = 0
        curfield.extv[i].vz = 0
    end

    #Velocities induced by free vortices on each other
    mutual_ind([curfield.tev; curfield.lev; curfield.extv])

    #Add the influence of velocities induced by bound vortices
    utemp = zeros(ntev + nlev + nextv)
    wtemp = zeros(ntev + nlev + nextv)
    utemp, wtemp = ind_vel(surf.bv, [map(q -> q.x, curfield.tev); map(q -> q.x, curfield.lev); map(q -> q.x, curfield.extv)], [map(q -> q.z, curfield.tev); map(q -> q.z, curfield.lev); map(q -> q.z, curfield.extv) ])

    for i = 1:ntev
        curfield.tev[i].vx += utemp[i]
        curfield.tev[i].vz += wtemp[i]
    end
    for i = ntev+1:ntev+nlev
        curfield.lev[i-ntev].vx += utemp[i]
        curfield.lev[i-ntev].vz += wtemp[i]
    end
    for i = ntev+nlev+1:ntev+nlev+nextv
        curfield.extv[i-ntev-nlev].vx += utemp[i]
        curfield.extv[i-ntev-nlev].vz += wtemp[i]
    end

    #Add effect of sources
    utemp, wtemp = ind_vel_src(surf.src, [map(q -> q.x, curfield.tev);
                                          map(q -> q.x, curfield.lev); map(q -> q.x, curfield.extv)],
                               [map(q -> q.z, curfield.tev); map(q -> q.z, curfield.lev);
                                map(q -> q.z, curfield.extv)])

    for i = 1:ntev
        curfield.tev[i].vx += utemp[i]
        curfield.tev[i].vz += wtemp[i]
    end
    for i = ntev+1:ntev+nlev
        curfield.lev[i-ntev].vx += utemp[i]
        curfield.lev[i-ntev].vz += wtemp[i]
    end
    for i = ntev+nlev+1:ntev+nlev+nextv
        curfield.extv[i-ntev-nlev].vx += utemp[i]
        curfield.extv[i-ntev-nlev].vz += wtemp[i]
    end

    #Convect free vortices with their induced velocities
    for i = 1:ntev
        curfield.tev[i].x += dt*curfield.tev[i].vx
        curfield.tev[i].z += dt*curfield.tev[i].vz
    end
    for i = 1:nlev
        curfield.lev[i].x += dt*curfield.lev[i].vx
        curfield.lev[i].z += dt*curfield.lev[i].vz
    end
    for i = 1:nextv
        curfield.extv[i].x += dt*curfield.extv[i].vx
        curfield.extv[i].z += dt*curfield.extv[i].vz
    end

    return curfield
end

function update_boundpos(surf::TwoDSurfThick, dt::Float64)
    for i = 1:surf.ndiv
        surf.bnd_x_u[i] = surf.bnd_x_u[i] + dt*((surf.pvt*surf.c - surf.x[i])*
                                                sin(surf.kinem.alpha)*surf.kinem.alphadot - surf.kinem.u + (surf.cam[i] +
                                                                                                            surf.thick[i])*cos(surf.kinem.alpha)*surf.kinem.alphadot)

        surf.bnd_z_u[i] = surf.bnd_z_u[i] + dt*(surf.kinem.hdot + (surf.pvt*surf.c -
                                                                   surf.x[i])*cos(surf.kinem.alpha)*surf.kinem.alphadot - (surf.cam[i] +
                                                                                                                           surf.thick[i])*sin(surf.kinem.alpha)*surf.kinem.alphadot)

        surf.bnd_x_l[i] = surf.bnd_x_l[i] + dt*((surf.pvt*surf.c - surf.x[i])*
                                                sin(surf.kinem.alpha)*surf.kinem.alphadot - surf.kinem.u + (surf.cam[i] -
                                                                                                            surf.thick[i])*cos(surf.kinem.alpha)*surf.kinem.alphadot)

        surf.bnd_z_l[i] = surf.bnd_z_l[i] + dt*(surf.kinem.hdot + (surf.pvt*
                                                                   surf.c - surf.x[i])*cos(surf.kinem.alpha)*surf.kinem.alphadot - (surf.cam[i] -
                                                                                                                                    surf.thick[i])*sin(surf.kinem.alpha)*surf.kinem.alphadot)

        surf.bnd_x_chord[i] = surf.bnd_x_chord[i] + dt*((surf.pvt*surf.c - surf.x[i])*
                                                        sin(surf.kinem.alpha)*surf.kinem.alphadot - surf.kinem.u)
        surf.bnd_z_chord[i] = surf.bnd_z_chord[i] + dt*(surf.kinem.hdot + (surf.pvt*surf.c -
                                                                           surf.x[i])*cos(surf.kinem.alpha)*surf.kinem.alphadot)
    end

    return surf
end

function update_thickLHS(surf::TwoDSurfThick, curfield::TwoDFlowField, dt::Float64, vcore::Float64)
    #Calculate the missing column in LHS that depends on last shed vortex location

    ntev = length(curfield.tev)

    if ntev == 0
        xloc_tev = surf.bnd_x_chord[surf.ndiv] + 0.5*surf.kinem.u*dt
        zloc_tev = surf.bnd_z_chord[surf.ndiv]
    else
        xloc_tev = surf.bnd_x_chord[surf.ndiv] + (1. /3.)*(curfield.tev[ntev].x - surf.bnd_x_chord[surf.ndiv])
        zloc_tev = surf.bnd_z_chord[surf.ndiv] + (1. /3.)*(curfield.tev[ntev].z - surf.bnd_z_chord[surf.ndiv])
    end

    dummyvort = TwoDVort(xloc_tev, zloc_tev, 1., vcore, 0., 0.)

    uu, wu = ind_vel([dummyvort], surf.bnd_x_u, surf.bnd_z_u)
    ul, wl = ind_vel([dummyvort], surf.bnd_x_l, surf.bnd_z_l)

    for i = 2:surf.ndiv-1


        #Sweep all rows (corresponding to ndiv) for lifting equation

        #Sweep columns for aterms
        for n = 1:surf.naterm
           surf.LHS[i-1,n] = cos(n*surf.theta[i]) - surf.thick_slope[i]*sin(n*surf.theta[i])
        end

        #Sweep columns for bterm
        for n = 1:surf.naterm
           surf.LHS[i-1,n+surf.naterm] = surf.cam_slope[i]*cos(n*surf.theta[i])
        end

        #TEV term must be updated in the loop after its location is known
        #Sweep all rows (corresponding to ndiv) for nonlifting equation

        for n = 1:surf.naterm
           surf.LHS[surf.ndiv+i-3,n]  = -surf.cam_slope[i]*sin(n*surf.theta[i])
        end
        for n = 1:surf.naterm
           surf.LHS[surf.ndiv+i-3,surf.naterm+n] = sin(n*surf.theta[i]) + surf.thick_slope[i]*cos(n*surf.theta[i])
        end


        wlz = 0.5*(wu[i]*cos(surf.kinem.alpha) + uu[i]*sin(surf.kinem.alpha) +
                   wl[i]*cos(surf.kinem.alpha) + ul[i]*sin(surf.kinem.alpha))

        wtz = 0.5*(wu[i]*cos(surf.kinem.alpha) + uu[i]*sin(surf.kinem.alpha) -
                   wl[i]*cos(surf.kinem.alpha) - ul[i]*sin(surf.kinem.alpha))

        wtx = 0.5*(uu[i]*cos(surf.kinem.alpha) - wu[i]*sin(surf.kinem.alpha) +
                   ul[i]*cos(surf.kinem.alpha) - wl[i]*sin(surf.kinem.alpha))

        wlx = 0.5*(uu[i]*cos(surf.kinem.alpha) - wu[i]*sin(surf.kinem.alpha) -
                   ul[i]*cos(surf.kinem.alpha) + wl[i]*sin(surf.kinem.alpha))

        surf.LHS[i-1,1+2*surf.naterm] = -surf.cam_slope[i]*wtx -
            surf.thick_slope[i]*wlx + wlz

        surf.LHS[surf.ndiv+i-3,1+2*surf.naterm] = -surf.cam_slope[i]*wlx -
            surf.thick_slope[i]*wtx + wtz

    end

    #surf.LHS[2*surf.ndiv-1,2+2*surf.naterm] = -wu[1]*cos(surf.kinem.alpha) -
    #uu[1]*sin(surf.kinem.alpha)


    return surf, xloc_tev, zloc_tev
end

# function update_thickLHSLEV(surf::TwoDSurfThick, curfield::TwoDFlowField, dt::Float64, vcore::Float64)
#     #Calculate the missing column in LHS that depends on last shed vortex location

#     ntev = length(curfield.tev)
#     nlev = length(curfield.lev)

#     if ntev == 0
#         xloc_tev = surf.bnd_x_chord[surf.ndiv] + 0.5*surf.kinem.u*dt
#         zloc_tev = surf.bnd_z_chord[surf.ndiv]
#     else
#         xloc_tev = surf.bnd_x_chord[surf.ndiv] + (1. /3.)*(curfield.tev[ntev].x - surf.bnd_x_chord[surf.ndiv])
#         zloc_tev = surf.bnd_z_chord[surf.ndiv] + (1. /3.)*(curfield.tev[ntev].z - surf.bnd_z_chord[surf.ndiv])
#     end

#     if surf.levflag[1] == 0

#         if surf.a0[1] >= 0.
#             le_vel_x = sqrt(2. /surf.rho)*surf.uref*surf.a0[1]*sin(surf.kinem.alpha)
#             le_vel_z = sqrt(2. /surf.rho)*surf.uref*surf.a0[1]*cos(surf.kinem.alpha)
#             xloc_lev = surf.bnd_x_u[1] + 0.5*le_vel_x*dt
#             zloc_lev = surf.bnd_z_u[1] + 0.5*le_vel_z*dt
#         else
#             le_vel_x = sqrt(2. /surf.rho)*surf.uref*surf.a0[1]*sin(surf.kinem.alpha)
#             le_vel_z = sqrt(2. /surf.rho)*surf.uref*surf.a0[1]*cos(surf.kinem.alpha)
#             xloc_lev = surf.bnd_x_l[1] + 0.5*le_vel_x*dt
#             zloc_lev = surf.bnd_z_l[1] + 0.5*le_vel_z*dt
#         end
#     else
#         if surf.a0[1] >= 0.
#             xloc_lev = surf.bnd_x_u[1]+(1. /3.)*(curfield.lev[nlev].x - surf.bnd_x_u[1])
#             zloc_lev = surf.bnd_z_u[1]+(1. /3.)*(curfield.lev[nlev].z - surf.bnd_z_u[1])
#         else
#             xloc_lev = surf.bnd_x_l[1]+(1. /3.)*(curfield.lev[nlev].x - surf.bnd_x_l[1])
#             zloc_lev = surf.bnd_z_l[1]+(1. /3.)*(curfield.lev[nlev].z - surf.bnd_z_l[1])
#         end
#     end

#     dummyvort = TwoDVort(xloc_tev, zloc_tev, 1., vcore, 0., 0.)

#     uu, wu = ind_vel([dummyvort], surf.bnd_x_u, surf.bnd_z_u)
#     ul, wl = ind_vel([dummyvort], surf.bnd_x_l, surf.bnd_z_l)

#     for i = 2:surf.ndiv-1
#         wlz = 0.5*(wu[i]*cos(surf.kinem.alpha) + uu[i]*sin(surf.kinem.alpha) +
#                    wl[i]*cos(surf.kinem.alpha) + ul[i]*sin(surf.kinem.alpha))

#         wtz = 0.5*(wu[i]*cos(surf.kinem.alpha) + uu[i]*sin(surf.kinem.alpha) -
#                    wl[i]*cos(surf.kinem.alpha) - ul[i]*sin(surf.kinem.alpha))

#         wtx = 0.5*(uu[i]*cos(surf.kinem.alpha) - wu[i]*sin(surf.kinem.alpha) +
#                    ul[i]*cos(surf.kinem.alpha) - wl[i]*sin(surf.kinem.alpha))

#         wlx = 0.5*(uu[i]*cos(surf.kinem.alpha) - wu[i]*sin(surf.kinem.alpha) -
#                    ul[i]*cos(surf.kinem.alpha) + wl[i]*sin(surf.kinem.alpha))

#         surf.LHS[i-1,2+2*surf.naterm] = -surf.cam_slope[i]*wtx -
#             surf.thick_slope[i]*wlx + wlz

#         surf.LHS[surf.ndiv+i-3,2+2*surf.naterm] = -surf.cam_slope[i]*wlx -
#             surf.thick_slope[i]*wtx + wtz

#     end

#     dummyvort = TwoDVort(xloc_lev, zloc_lev, 1., vcore, 0., 0.)

#     uu, wu = ind_vel([dummyvort], surf.bnd_x_u, surf.bnd_z_u)
#     ul, wl = ind_vel([dummyvort], surf.bnd_x_l, surf.bnd_z_l)

#     for i = 2:surf.ndiv-1
#         wlz = 0.5*(wu[i]*cos(surf.kinem.alpha) + uu[i]*sin(surf.kinem.alpha) +
#                    wl[i]*cos(surf.kinem.alpha) + ul[i]*sin(surf.kinem.alpha))

#         wtz = 0.5*(wu[i]*cos(surf.kinem.alpha) + uu[i]*sin(surf.kinem.alpha) -
#                    wl[i]*cos(surf.kinem.alpha) - ul[i]*sin(surf.kinem.alpha))

#         wtx = 0.5*(uu[i]*cos(surf.kinem.alpha) - wu[i]*sin(surf.kinem.alpha) +
#                    ul[i]*cos(surf.kinem.alpha) - wl[i]*sin(surf.kinem.alpha))

#         wlx = 0.5*(uu[i]*cos(surf.kinem.alpha) - wu[i]*sin(surf.kinem.alpha) -
#                    ul[i]*cos(surf.kinem.alpha) + wl[i]*sin(surf.kinem.alpha))

#         surf.LHS[i-1,3+2*surf.naterm] = -surf.cam_slope[i]*wtx -
#             surf.thick_slope[i]*wlx + wlz

#         surf.LHS[surf.ndiv+i-3,3+2*surf.naterm] = -surf.cam_slope[i]*wlx -
#             surf.thick_slope[i]*wtx + wtz
#     end

#     return surf, xloc_tev, zloc_tev, xloc_lev, zloc_lev
# end

function update_thickRHS(surf::TwoDSurfThick, curfield::TwoDFlowField)
    for i = 2:surf.ndiv-1
        #RHS for lifting equation

        wlz = 0.5*(surf.wind_u[i]*cos(surf.kinem.alpha) + surf.uind_u[i]*sin(surf.kinem.alpha) +
                   surf.wind_l[i]*cos(surf.kinem.alpha) + surf.uind_l[i]*sin(surf.kinem.alpha))

        wtz = 0.5*(surf.wind_u[i]*cos(surf.kinem.alpha) + surf.uind_u[i]*sin(surf.kinem.alpha) -
                   surf.wind_l[i]*cos(surf.kinem.alpha) - surf.uind_l[i]*sin(surf.kinem.alpha))

        wtx = 0.5*(surf.uind_u[i]*cos(surf.kinem.alpha) - surf.wind_u[i]*sin(surf.kinem.alpha) +
                   surf.uind_l[i]*cos(surf.kinem.alpha) - surf.wind_l[i]*sin(surf.kinem.alpha))

        wlx = 0.5*(surf.uind_u[i]*cos(surf.kinem.alpha) - surf.wind_u[i]*sin(surf.kinem.alpha) -
                   surf.uind_l[i]*cos(surf.kinem.alpha) + surf.wind_l[i]*sin(surf.kinem.alpha))

        surf.RHS[i-1] = -(surf.kinem.u + curfield.u[1])*sin(surf.kinem.alpha)/surf.uref -
            surf.kinem.alphadot*(surf.x[i] - surf.pvt*surf.c)/surf.uref +
            (surf.kinem.hdot - curfield.w[1])*cos(surf.kinem.alpha)/surf.uref - wlz/surf.uref +
            surf.cam_slope[i]*((surf.kinem.u + curfield.u[1])*cos(surf.kinem.alpha) +
                               (surf.kinem.hdot - curfield.w[1])*sin(surf.kinem.alpha) + wtx -
                               surf.kinem.alphadot*surf.cam[i]) +
                               surf.thick_slope[i]*(wlx - surf.kinem.alphadot*surf.thick[i])

        surf.RHS[surf.ndiv+i-3] = surf.cam_slope[i]*(wlx - surf.kinem.alphadot*surf.thick[i]) +
            surf.thick_slope[i]*((surf.kinem.u + curfield.u[1])*cos(surf.kinem.alpha) +
                                 (surf.kinem.hdot - curfield.w[1])*sin(surf.kinem.alpha) + wtx -
                                 surf.kinem.alphadot*surf.cam[i]) - wtz
    end

    #RHS for Kelvin condition (negative strength of all previously shed vortices)
    surf.RHS[2*surf.ndiv-3] = -100*(sum(map(q->q.s, curfield.tev)) + sum(map(q->q.s, curfield.lev)))/(surf.uref*surf.c)

    #Kutta condition
    surf.RHS[2*surf.ndiv-2] = 0.
    #surf.RHS[2*surf.ndiv-3] = pi*(surf.a0prev[1] + 0.5*surf.aprev[1])

        #RHS for LE velocity condition
    #surf.RHS[2*surf.ndiv-1] = sin(surf.kinem.alpha) - surf.kinem.hdot*cos(surf.kinem.alpha)/surf.uref -
    #surf.kinem.alphadot*surf.pvt*surf.c/surf.uref +
    #   (surf.wind_u[1]*cos(surf.kinem.alpha) + surf.uind_u[1]*sin(surf.kinem.alpha))/surf.uref


   # surf.RHS[2*surf.ndiv-2] = -(surf.kinem.u + curfield.u[1])*cos(surf.kinem.alpha) - (surf.kinem.hdot - curfield.w[1])*sin(surf.kinem.alpha)

    #surf.RHS[2*surf.ndiv-1] = -(surf.kinem.u + curfield.u[1])*cos(surf.kinem.alpha) - (surf.kinem.hdot - curfield.w[1])*sin(surf.kinem.alpha)

end

# function update_thickRHSLEV(surf::TwoDSurfThick, curfield::TwoDFlowField)
#     for i = 2:surf.ndiv-1
#         #RHS for lifting equation

#         wlz = 0.5*(surf.wind_u[i]*cos(surf.kinem.alpha) + surf.uind_u[i]*sin(surf.kinem.alpha) +
#                    surf.wind_l[i]*cos(surf.kinem.alpha) + surf.uind_l[i]*sin(surf.kinem.alpha))

#         wtz = 0.5*(surf.wind_u[i]*cos(surf.kinem.alpha) + surf.uind_u[i]*sin(surf.kinem.alpha) -
#                    surf.wind_l[i]*cos(surf.kinem.alpha) - surf.uind_l[i]*sin(surf.kinem.alpha))

#         wtx = 0.5*(surf.uind_u[i]*cos(surf.kinem.alpha) - surf.wind_u[i]*sin(surf.kinem.alpha) +
#                    surf.uind_l[i]*cos(surf.kinem.alpha) - surf.wind_l[i]*sin(surf.kinem.alpha))

#         wlx = 0.5*(surf.uind_u[i]*cos(surf.kinem.alpha) - surf.wind_u[i]*sin(surf.kinem.alpha) -
#                    surf.uind_l[i]*cos(surf.kinem.alpha) + surf.wind_l[i]*sin(surf.kinem.alpha))

#         surf.RHS[i-1] = -(surf.kinem.u + curfield.u[1])*sin(surf.kinem.alpha)/surf.uref -
#             surf.kinem.alphadot*(surf.x[i] - surf.pvt*surf.c)/surf.uref +
#             (surf.kinem.hdot - curfield.w[1])*cos(surf.kinem.alpha)/surf.uref - wlz/surf.uref +
#             surf.cam_slope[i]*((surf.kinem.u + curfield.u[1])*cos(surf.kinem.alpha) +
#                                (surf.kinem.hdot - curfield.w[1])*sin(surf.kinem.alpha) + wtx -
#                                surf.kinem.alphadot*surf.cam[i]) +
#                                surf.thick_slope[i]*(wlx - surf.kinem.alphadot*surf.thick[i])

#         surf.RHS[surf.ndiv+i-3] = surf.cam_slope[i]*(wlx - surf.kinem.alphadot*surf.thick[i]) +
#             surf.thick_slope[i]*((surf.kinem.u + curfield.u[1])*cos(surf.kinem.alpha) +
#                                  (surf.kinem.hdot - curfield.w[1])*sin(surf.kinem.alpha) + wtx -
#                                  surf.kinem.alphadot*surf.cam[i]) - wtz
#     end

#     #RHS for Kelvin condition (negative strength of all previously shed vortices)
#     surf.RHS[2*surf.ndiv-3] = pi*(surf.a0prev[1] + 0.5*surf.aprev[1])

#     # #RHS for LE Kutta condition
#     # if surf.a0[1] > 0.
#     #     surf.RHS[2*surf.ndiv-2] = 1000*sqrt(surf.rho/ 2. )*surf.lespcrit[1]
#     # else
#     #     surf.RHS[2*surf.ndiv-2] = -1000*sqrt(surf.rho/ 2. )*surf.lespcrit[1]
#     # end

#     # #RHS for LE velocity condition
#     # surf.RHS[2*surf.ndiv-1] = sin(surf.kinem.alpha) - surf.kinem.hdot*cos(surf.kinem.alpha)/surf.uref -
#     #     surf.kinem.alphadot*surf.pvt*surf.c/surf.uref +
#     #     (surf.wind_u[1]*cos(surf.kinem.alpha) + surf.uind_u[1]*sin(surf.kinem.alpha))/surf.uref

# end


#Update position and strengths of bound vortices and sources
function update_bv_src(surf::TwoDSurfThick)
    gamma = zeros(surf.ndiv)
    src = zeros(surf.ndiv)
    for ib = 1:surf.ndiv
        gamma[ib] = surf.a0[1]*(1 + cos(surf.theta[ib]))
        for ia = 1:surf.naterm
            gamma[ib] += surf.aterm[ia]*sin(ia*surf.theta[ib])*sin(surf.theta[ib])
            src[ib] += surf.bterm[ia]*sin(ia*surf.theta[ib])*sin(surf.theta[ib])
        end
        gamma[ib] = gamma[ib]*surf.uref*surf.c
        src[ib] = src[ib]*surf.uref*surf.c
    end

    for ib = 2:surf.ndiv
        surf.bv[ib-1].s = (gamma[ib]+gamma[ib-1])*(surf.theta[2]-surf.theta[1])/2.
        surf.bv[ib-1].x = (surf.bnd_x_chord[ib] + surf.bnd_x_chord[ib-1])/2.
        surf.bv[ib-1].z = (surf.bnd_z_chord[ib] + surf.bnd_z_chord[ib-1])/2.
        surf.src[ib-1].s = (src[ib] + src[ib-1])*(surf.theta[2]-surf.theta[1])/2.
        surf.src[ib-1].x = surf.bv[ib-1].x
        surf.src[ib-1].z = surf.bv[ib-1].z
    end
end

# Function to calculate induced velocities by a set of sources at a target location
function ind_vel_src(src::Vector{TwoDSource},t_x,t_z)
    #'s' stands for source and 't' for target
    uind = zeros(length(t_x))
    wind = zeros(length(t_x))

    for itr = 1:length(t_x)
	for isr = 1:length(src)
            xdist = src[isr].x - t_x[itr]
            zdist = src[isr].z - t_z[itr]
            distsq = xdist*xdist + zdist*zdist
            uind[itr] = uind[itr] - src[isr].s*xdist/(2*pi*distsq)
            wind[itr] = wind[itr] - src[isr].s*zdist/(2*pi*distsq)
        end
    end
    return uind, wind
end
