function update_atermdot(surf::TwoDSurfThick,dt)
    surf.atedot[1] = (surf.ate[1] - surf.ateprev[1])/dt
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

function add_indbound_lasttev(surf::TwoDSurfThick, curfield::TwoDFlowField)
    uind_u, wind_u = ind_vel([curfield.tev[end]], surf.bnd_x_u, surf.bnd_z_u)
    uind_l, wind_l = ind_vel([curfield.tev[end]], surf.bnd_x_l, surf.bnd_z_l)
    
    surf.uind_u += uind_u
    surf.wind_u += wind_u
    surf.uind_l += uind_l
    surf.wind_l += wind_l
    
    return surf
end

function add_indbound_lastlev(surf::TwoDSurfThick, curfield::TwoDFlowField)
    uind_u, wind_u = ind_vel([curfield.lev[end]], surf.bnd_x_u, surf.bnd_z_u)
    uind_l, wind_l = ind_vel([curfield.lev[end]], surf.bnd_x_l, surf.bnd_z_l)
    
    surf.uind_u += uind_u
    surf.wind_u += wind_u
    surf.uind_l += uind_l
    surf.wind_l += wind_l
    
    return surf
end

function minus_indbound_lasttev(surf::TwoDSurfThick, curfield::TwoDFlowField)
    uind_u, wind_u = ind_vel([curfield.tev[end]], surf.bnd_x_u, surf.bnd_z_u)
    uind_l, wind_l = ind_vel([curfield.tev[end]], surf.bnd_x_l, surf.bnd_z_l)
    
    surf.uind_u -= uind_u
    surf.wind_u -= wind_u
    surf.uind_l -= uind_l
    surf.wind_l -= wind_l
    
    return surf
end

function minus_indbound_lastlev(surf::TwoDSurfThick, curfield::TwoDFlowField)
    uind_u, wind_u = ind_vel([curfield.lev[end]], surf.bnd_x_u, surf.bnd_z_u)
    uind_l, wind_l = ind_vel([curfield.lev[end]], surf.bnd_x_l, surf.bnd_z_l)
    
    surf.uind_u -= uind_u
    surf.wind_u -= wind_u
    surf.uind_l -= uind_l
    surf.wind_l -= wind_l
    
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
    elseif (typeof(surf.kindef.alpha) == EldDownDef)
        surf.kinem.alpha = surf.kindef.alpha(t)
        surf.kinem.alphadot = ForwardDiff.derivative(surf.kindef.alpha,t)*surf.uref/surf.c
    elseif (typeof(surf.kindef.alpha) == EldUpDownDef)
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
    elseif (typeof(surf.kindef.h) == EldDownDef)
        surf.kinem.h = surf.kindef.h(t)*surf.c
        surf.kinem.hdot = ForwardDiff.derivative(surf.kindef.h,t)*surf.uref
    elseif (typeof(surf.kindef.h) == EldUpDownDef)
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
    elseif (typeof(surf.kindef.u) == EldDownDef)
        surf.kinem.u, surf.kinem.udot = surf.kindef.u(t)
        surf.kinem.u = surf.kinem.u*surf.uref
        surf.kinem.udot = surf.kinem.udot*surf.uref*surf.uref/surf.c
    elseif (typeof(surf.kindef.u) == EldUpDownDef)
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


#Update position and strengths of bound vortices and sources
function update_bv_src(surf::TwoDSurfThick)
    gamma = zeros(surf.ndiv)
    src = zeros(surf.ndiv)
    for ib = 1:surf.ndiv
        gamma[ib] = surf.ate[1]*(1. - cos(surf.theta[ib]))
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



function update_LHSRHS(surf::TwoDSurfThick, curfield::TwoDFlowField, dt::Float64, vcore::Float64, bc_prev::Float64)
    #Calculate the missing column in LHS that depends on last shed vortex location

    ntev = length(curfield.tev)

    if ntev == 0
        xloc_tev = surf.bnd_x_chord[surf.ndiv] + 0.5*surf.kinem.u*dt*cos(surf.kinem.alpha)
        zloc_tev = surf.bnd_z_chord[surf.ndiv] - 0.5*surf.kinem.u*dt*cos(surf.kinem.alpha)
    else
        xloc_tev = surf.bnd_x_chord[surf.ndiv] + (1. /3.)*(curfield.tev[ntev].x - surf.bnd_x_chord[surf.ndiv])
        zloc_tev = surf.bnd_z_chord[surf.ndiv] + (1. /3.)*(curfield.tev[ntev].z - surf.bnd_z_chord[surf.ndiv])
    end

    dummyvort = TwoDVort(xloc_tev, zloc_tev, 1., vcore, 0., 0.)

    uu, wu = ind_vel([dummyvort], surf.bnd_x_u, surf.bnd_z_u)
    ul, wl = ind_vel([dummyvort], surf.bnd_x_l, surf.bnd_z_l)

    wlz_t = 0.5*(wu.*cos(surf.kinem.alpha) .+ uu.*sin(surf.kinem.alpha) .+ wl.*cos(surf.kinem.alpha) .+ ul.*sin(surf.kinem.alpha))
    wtz_t = 0.5*(wu.*cos(surf.kinem.alpha) .+ uu.*sin(surf.kinem.alpha) .- wl.*cos(surf.kinem.alpha) .- ul.*sin(surf.kinem.alpha))
    wtx_t = 0.5*(uu.*cos(surf.kinem.alpha) .- wu.*sin(surf.kinem.alpha) .+ ul.*cos(surf.kinem.alpha) .- wl.*sin(surf.kinem.alpha))
    wlx_t = 0.5*(uu.*cos(surf.kinem.alpha) .- wu.*sin(surf.kinem.alpha) .- ul.*cos(surf.kinem.alpha) .+ wl.*sin(surf.kinem.alpha))
    
    wlz = 0.5*(surf.wind_u.*cos(surf.kinem.alpha) .+ surf.uind_u.*sin(surf.kinem.alpha) .+
               surf.wind_l*cos(surf.kinem.alpha) .+ surf.uind_l*sin(surf.kinem.alpha))
    
    wtz = 0.5*(surf.wind_u.*cos(surf.kinem.alpha) .+ surf.uind_u.*sin(surf.kinem.alpha) .-
               surf.wind_l.*cos(surf.kinem.alpha) .- surf.uind_l.*sin(surf.kinem.alpha))
    
    wtx = 0.5*(surf.uind_u.*cos(surf.kinem.alpha) .- surf.wind_u.*sin(surf.kinem.alpha) .+
               surf.uind_l.*cos(surf.kinem.alpha) .- surf.wind_l.*sin(surf.kinem.alpha))
    
    wlx = 0.5*(surf.uind_u.*cos(surf.kinem.alpha) .- surf.wind_u.*sin(surf.kinem.alpha) .-
               surf.uind_l.*cos(surf.kinem.alpha) .+ surf.wind_l.*sin(surf.kinem.alpha))

    
    sw1 = 0; sw2 = 0; sw3 = 0; sw4 = 0;
    sw1_t = 0; sw2_t = 0; sw3_t = 0; sw4_t = 0;
    
    for i = 2:surf.ndiv-1
        surf.LHS[i-1,1+2*surf.naterm] = -surf.cam_slope[i]*wtx_t[i] -
            surf.thick_slope[i]*wlx_t[i] + wlz_t[i]
        
        surf.LHS[surf.ndiv+i-3,1+2*surf.naterm] = -surf.cam_slope[i]*wlx_t[i] -
            surf.thick_slope[i]*wtx_t[i] + wtz_t[i]
        
        surf.RHS[i-1] = -(surf.kinem.u + curfield.u[1])*sin(surf.kinem.alpha)/surf.uref -
            surf.kinem.alphadot*(surf.x[i] - surf.pvt*surf.c)/surf.uref +
            (surf.kinem.hdot - curfield.w[1])*cos(surf.kinem.alpha)/surf.uref - wlz[i]/surf.uref +
            surf.cam_slope[i]*((surf.kinem.u + curfield.u[1])*cos(surf.kinem.alpha) +
                               (surf.kinem.hdot - curfield.w[1])*sin(surf.kinem.alpha) + wtx[i] -
                               surf.kinem.alphadot*surf.cam[i]) +
                               surf.thick_slope[i]*(wlx[i] - surf.kinem.alphadot*surf.thick[i])
        
        surf.RHS[surf.ndiv+i-3] = surf.cam_slope[i]*(wlx[i] - surf.kinem.alphadot*surf.thick[i]) +
            surf.thick_slope[i]*((surf.kinem.u + curfield.u[1])*cos(surf.kinem.alpha) +
                                 (surf.kinem.hdot - curfield.w[1])*sin(surf.kinem.alpha) + wtx[i] -
                                 surf.kinem.alphadot*surf.cam[i]) - wtz[i]
    end
    
    for i = 2:surf.ndiv-1
        ds = sqrt((surf.x[i] - surf.x[i-1])^2 + (surf.cam[i] + surf.thick[i] - surf.cam[i-1] - surf.thick[i-1])^2)
        den_u = sqrt(1. + (surf.cam_slope[i] + surf.thick_slope[i])^2)
        den_u_p = sqrt(1. + (surf.cam_slope[i-1] + surf.thick_slope[i-1])^2)
        den_l = sqrt(1. + (surf.cam_slope[i] - surf.thick_slope[i])^2)
        den_l_p = sqrt(1. + (surf.cam_slope[i-1] - surf.thick_slope[i-1])^2)
        
        sw1_t += 0.5*ds*((wtx_t[i] + wlx_t[i])/den_u + (wtx_t[i-1] + wlx_t[i-1])/den_u_p)
        sw2_t += 0.5*ds*((surf.cam_slope[i] + surf.thick_slope[i])*(wtz_t[i] + wlz_t[i])/den_u + (surf.cam_slope[i-1] + surf.thick_slope[i-1])*(wtz_t[i-1] + wlz_t[i-1])/den_u_p)
        sw1 += 0.5*ds*((wtx[i] + wlx[i])/den_u + (wtx[i-1] + wlx[i-1])/den_u_p)
        sw2 += 0.5*ds*((surf.cam_slope[i] + surf.thick_slope[i])*(wtz[i] + wlz[i])/den_u + (surf.cam_slope[i-1] + surf.thick_slope[i-1])*(wtz[i-1] + wlz[i-1])/den_u_p)
        
        ds = sqrt((surf.x[i] - surf.x[i-1])^2 + (surf.cam[i] - surf.thick[i] - surf.cam[i-1] + surf.thick[i-1])^2)
        sw3_t += 0.5*ds*((wtx_t[i] - wlx_t[i])/den_l + (wtx_t[i-1] - wlx_t[i-1])/den_l_p)
        sw4_t += 0.5*ds*((surf.cam_slope[i] - surf.thick_slope[i])*(wlz_t[i] - wtz_t[i])/den_l + (surf.cam_slope[i-1] - surf.thick_slope[i-1])*(wlz_t[i-1] - wtz_t[i-1])/den_l_p)
        
        
        sw3 += 0.5*ds*((wtx[i] - wlx[i])/den_l + (wtx[i-1] - wlx[i-1])/den_l_p)
        sw4 += 0.5*ds*((surf.cam_slope[i] - surf.thick_slope[i])*(wlz[i] - wtz[i])/den_l + (surf.cam_slope[i-1] - surf.thick_slope[i-1])*(wlz[i-1] - wtz[i-1])/den_l_p)
    end
    
    surf.LHS[2*surf.ndiv-3, 2*surf.naterm+1] = (1. + sw1_t + sw2_t - sw3_t - sw4_t)
    surf.RHS[2*surf.ndiv-3] = (bc_prev - sw1 - sw2 + sw3 + sw4)
    
    #Kutta condition
    i = surf.ndiv-1
    vref_x_u = (surf.kinem.u + curfield.u[1])*cos(surf.kinem.alpha) + (surf.kinem.hdot - curfield.w[1])*sin(surf.kinem.alpha) - surf.kinem.alphadot*(surf.cam[i] + surf.thick[i])
    vref_x_l = (surf.kinem.u + curfield.u[1])*cos(surf.kinem.alpha) + (surf.kinem.hdot - curfield.w[1])*sin(surf.kinem.alpha) - surf.kinem.alphadot*(surf.cam[i] - surf.thick[i])
    vref_z = (surf.kinem.u + curfield.u[1])*sin(surf.kinem.alpha) - (surf.kinem.hdot - curfield.w[1])*cos(surf.kinem.alpha) + surf.kinem.alphadot*(surf.x[i] - surf.pvt*surf.c)
    tx = -surf.bterm'*cos.(collect(1:surf.naterm)*surf.theta[i])
    tz = surf.bterm'*sin.(collect(1:surf.naterm)*surf.theta[i])
    lx = surf.aterm'*sin.(collect(1:surf.naterm)*surf.theta[i]) + surf.ate[1]*tan(surf.theta[i]/2)
    lz = surf.aterm'*cos.(collect(1:surf.naterm)*surf.theta[i]) + surf.ate[1]
    
    for n = 1:surf.naterm
        surf.LHS[2*surf.ndiv-2,n] = (-dt*(vref_x_u + vref_x_l + 2*tx + 2*wtx[i])*sin(n*surf.theta[i]))
        #surf.LHS[2*surf.ndiv-2,n+surf.naterm] = -2*dt*(vref_z + lz + wlz[i])*sin(n*surf.theta[i])
    end
    surf.LHS[2*surf.ndiv-2,2*surf.naterm+1] = (1. - dt*(vref_x_l + vref_x_u + 2*tx + 2*wtx[i])*wlx_t[i]) #- 2*dt*(vref_z + lz + wlz[i])*wtz_t[i]
    surf.LHS[2*surf.ndiv-2,2*surf.naterm+2] = (-dt*(vref_x_l + vref_x_u + 2*tx + 2*wtx[i])*tan(surf.theta[i]/2))
    surf.RHS[2*surf.ndiv-2] = (dt*(vref_x_u + vref_x_l + 2*tx + 2*wtx[i])*wlx[i]) #+ 2*dt*(vref_z + lz + wlz[i])*wtz[i]
    
    return surf, xloc_tev, zloc_tev
end

function update_LHSRHS_kutta(surf::TwoDSurfThick, curfield::TwoDFlowField, dt::Float64, tevstr::Float64, aorig::Array{Float64}, borig::Array{Float64})
    i = surf.ndiv-1

    wlz = 0.5*(surf.wind_u[i].*cos(surf.kinem.alpha) .+ surf.uind_u[i].*sin(surf.kinem.alpha) .+
               surf.wind_l[i]*cos(surf.kinem.alpha) .+ surf.uind_l[i]*sin(surf.kinem.alpha))
    
    wtz = 0.5*(surf.wind_u[i].*cos(surf.kinem.alpha) .+ surf.uind_u[i].*sin(surf.kinem.alpha) .-
               surf.wind_l[i].*cos(surf.kinem.alpha) .- surf.uind_l[i].*sin(surf.kinem.alpha))
    
    wtx = 0.5*(surf.uind_u[i].*cos(surf.kinem.alpha) .- surf.wind_u[i].*sin(surf.kinem.alpha) .+
               surf.uind_l[i].*cos(surf.kinem.alpha) .- surf.wind_l[i].*sin(surf.kinem.alpha))
    
    wlx = 0.5*(surf.uind_u[i].*cos(surf.kinem.alpha) .- surf.wind_u[i].*sin(surf.kinem.alpha) .-
               surf.uind_l[i].*cos(surf.kinem.alpha) .+ surf.wind_l[i].*sin(surf.kinem.alpha))
    
    
    vref_x_u = (surf.kinem.u + curfield.u[1])*cos(surf.kinem.alpha) + (surf.kinem.hdot - curfield.w[1])*sin(surf.kinem.alpha) - surf.kinem.alphadot*(surf.cam[i] + surf.thick[i])
    vref_x_l = (surf.kinem.u + curfield.u[1])*cos(surf.kinem.alpha) + (surf.kinem.hdot - curfield.w[1])*sin(surf.kinem.alpha) - surf.kinem.alphadot*(surf.cam[i] - surf.thick[i])
    vref_z = (surf.kinem.u + curfield.u[1])*sin(surf.kinem.alpha) - (surf.kinem.hdot - curfield.w[1])*cos(surf.kinem.alpha) + surf.kinem.alphadot*(surf.x[i] - surf.pvt*surf.c)
    tx = -borig'*cos.(collect(1:surf.naterm)*surf.theta[i])
    tz = borig'*sin.(collect(1:surf.naterm)*surf.theta[i])
    lx = aorig'*sin.(collect(1:surf.naterm)*surf.theta[i]) + surf.ate[1]*tan(surf.theta[i]/2)
    lz = aorig'*cos.(collect(1:surf.naterm)*surf.theta[i]) + surf.ate[1]
    
    for n = 1:surf.naterm
        surf.LHS[2*surf.ndiv-3,n] = (-dt*(vref_x_u + vref_x_l + 2*tx + 2*wtx)*sin(n*surf.theta[i]))
        #surf.LHS[2*surf.ndiv-3,n+surf.naterm] = -2*dt*(vref_z + lz + wlz)*sin(n*surf.theta[i])
    end
    surf.LHS[2*surf.ndiv-3,2*surf.naterm+1] = (-dt*(vref_x_l + vref_x_u + 2*tx + 2*wtx)*tan(surf.theta[i]/2))
    surf.RHS[2*surf.ndiv-3] = (dt*(vref_x_u + vref_x_l + 2*tx + 2*wtx)*wlx) - tevstr #+ 2*dt*(vref_z + lz + wlz)*wtz
end

function update_RHS(surf::TwoDSurfThick, curfield::TwoDFlowField)
     
    wlz = 0.5*(surf.wind_u.*cos(surf.kinem.alpha) .+ surf.uind_u.*sin(surf.kinem.alpha) .+
               surf.wind_l*cos(surf.kinem.alpha) .+ surf.uind_l*sin(surf.kinem.alpha))
    
    wtz = 0.5*(surf.wind_u.*cos(surf.kinem.alpha) .+ surf.uind_u.*sin(surf.kinem.alpha) .-
               surf.wind_l.*cos(surf.kinem.alpha) .- surf.uind_l.*sin(surf.kinem.alpha))
    
    wtx = 0.5*(surf.uind_u.*cos(surf.kinem.alpha) .- surf.wind_u.*sin(surf.kinem.alpha) .+
               surf.uind_l.*cos(surf.kinem.alpha) .- surf.wind_l.*sin(surf.kinem.alpha))
    
    wlx = 0.5*(surf.uind_u.*cos(surf.kinem.alpha) .- surf.wind_u.*sin(surf.kinem.alpha) .-
               surf.uind_l.*cos(surf.kinem.alpha) .+ surf.wind_l.*sin(surf.kinem.alpha))

    for i = 2:surf.ndiv-1
        surf.RHS[i-1] = -(surf.kinem.u + curfield.u[1])*sin(surf.kinem.alpha)/surf.uref -
            surf.kinem.alphadot*(surf.x[i] - surf.pvt*surf.c)/surf.uref +
            (surf.kinem.hdot - curfield.w[1])*cos(surf.kinem.alpha)/surf.uref - wlz[i]/surf.uref +
            surf.cam_slope[i]*((surf.kinem.u + curfield.u[1])*cos(surf.kinem.alpha) +
                               (surf.kinem.hdot - curfield.w[1])*sin(surf.kinem.alpha) + wtx[i] -
                               surf.kinem.alphadot*surf.cam[i]) +
                               surf.thick_slope[i]*(wlx[i] - surf.kinem.alphadot*surf.thick[i])
        
        surf.RHS[surf.ndiv+i-3] = surf.cam_slope[i]*(wlx[i] - surf.kinem.alphadot*surf.thick[i]) +
            surf.thick_slope[i]*((surf.kinem.u + curfield.u[1])*cos(surf.kinem.alpha) +
                                 (surf.kinem.hdot - curfield.w[1])*sin(surf.kinem.alpha) + wtx[i] -
                                 surf.kinem.alphadot*surf.cam[i]) - wtz[i]
    end

    return surf
end

function calc_phi(surf::TwoDSurfThick)

    uphi_u = 0
    wphi_u = 0
    uphi_l = 0
    wphi_l = 0

    phi_u = 0.
    phi_l = 0.
    
    for i = 1:surf.ndiv-1

        l_x = 0; l_z = 0; t_x = 0; t_z = 0;
        for n = 1:surf.naterm
            l_x += surf.aterm[n]*sin(n*surf.theta[i])
            l_z += surf.aterm[n]*cos(n*surf.theta[i])
            t_x -= surf.bterm[n]*cos(n*surf.theta[i])
            t_z += surf.bterm[n]*sin(n*surf.theta[i])
        end
        l_x += surf.ate[1]*tan(surf.theta[i]/2)
        l_z += surf.ate[1]

        l_x *= surf.uref
        l_z *= surf.uref
        t_x *= surf.uref
        t_z *= surf.uref

        w_x_u = surf.uind_u[i]*cos(surf.kinem.alpha) - surf.wind_u[i]*sin(surf.kinem.alpha)
        w_z_u = surf.uind_u[i]*sin(surf.kinem.alpha) + surf.wind_u[i]*cos(surf.kinem.alpha)
        w_x_l = surf.uind_l[i]*cos(surf.kinem.alpha) - surf.wind_l[i]*sin(surf.kinem.alpha)
        w_z_l = surf.uind_l[i]*sin(surf.kinem.alpha) + surf.wind_l[i]*cos(surf.kinem.alpha)

        uphi_u_p = uphi_u
        wphi_u_p = wphi_u
        uphi_l_p = uphi_l
        wphi_l_p = wphi_l
        
        uphi_u = w_x_u + l_x + t_x
        wphi_u = w_z_u + l_z + t_z
        uphi_l = w_x_l - l_x + t_x
        wphi_l = w_z_l + l_z - t_z
        
        if i != 1
            
            ds = sqrt((surf.x[i] - surf.x[i-1])^2 + (surf.cam[i] + surf.thick[i] - surf.cam[i-1] - surf.thick[i-1])^2)
            if i ==2
                val1 = 0.5*(uphi_u/sqrt(1. + (surf.thick_slope[i] + surf.cam_slope[i])^2))  
                val2 = 0.5*(wphi_u*(surf.thick_slope[i] + surf.cam_slope[i])/sqrt(1. + (surf.thick_slope[i] + surf.cam_slope[i])^2) + wphi_u_p)  
            else
                val1 = 0.5*(uphi_u/sqrt(1. + (surf.thick_slope[i] + surf.cam_slope[i])^2) + uphi_u_p/sqrt(1. + (surf.thick_slope[i-1] + surf.cam_slope[i-1])^2))  
                val2 = 0.5*(wphi_u*(surf.thick_slope[i] + surf.cam_slope[i])/sqrt(1. + (surf.thick_slope[i] + surf.cam_slope[i])^2) + wphi_u_p*(surf.thick_slope[i-1] + surf.cam_slope[i-1])/sqrt(1. + (surf.thick_slope[i-1] + surf.cam_slope[i-1])^2))  
            end
            phi_u += val1*ds + val2*ds

            ds = sqrt((surf.x[i] - surf.x[i-1])^2 + (surf.cam[i] - surf.thick[i] - surf.cam[i-1] + surf.thick[i-1])^2)
            if i == 2
                val1 = 0.5*(uphi_l/sqrt(1. + (-surf.thick_slope[i] + surf.cam_slope[i])^2))
                val2 = 0.5*(wphi_l*(-surf.thick_slope[i] + surf.cam_slope[i])/sqrt(1. + (-surf.thick_slope[i] + surf.cam_slope[i])^2) + wphi_l_p)  
            else
                val1 = 0.5*(uphi_l/sqrt(1. + (-surf.thick_slope[i] + surf.cam_slope[i])^2) + uphi_l_p/sqrt(1. + (-surf.thick_slope[i-1] + surf.cam_slope[i-1])^2))  
                val2 = 0.5*(wphi_l*(-surf.thick_slope[i] + surf.cam_slope[i])/sqrt(1. + (-surf.thick_slope[i] + surf.cam_slope[i])^2) + wphi_l_p*(-surf.thick_slope[i-1] + surf.cam_slope[i-1])/sqrt(1. + (-surf.thick_slope[i-1] + surf.cam_slope[i-1])^2))  
            end
            phi_l += val1*ds + val2*ds
             
        end
        
    end
    
    return phi_u, phi_l
end

function find_stag(surf::TwoDSurfThick, qu::Array{Float64}, ql::Array{Float64})

    if qu[1] > 0
        #stagnation on lower surface
        i_s = 1
        for i = 1:surf.ndiv-1
            if ql[i]*ql[i+1] < 0
                #stagnation lies between these points
                i_s = i
                break
            end
        end
        xs = (surf.x[i_s]*ql[i_s+1] - surf.x[i_s+1]*ql[i_s])/(ql[i_s+1] - ql[i_s])
    else
        i_s = 1
        for i = 1:surf.ndiv-1
            if qu[i]*qu[i+1] < 0
                #stagnation lies between these points
                i_s = i
                break
            end
        end
        xs = (surf.x[i_s]*qu[i_s+1] - surf.x[i_s+1]*qu[i_s])/(qu[i_s+1] - qu[i_s])
    end

    return xs
    
end

# Places a trailing edge vortex
function place_tev(surf::TwoDSurfThick,field::TwoDFlowField,dt)
    ntev = length(field.tev)
    if ntev == 0
        xloc = surf.bnd_x_chord[surf.ndiv] + 0.5*surf.kinem.u*dt*cos(surf.kinem.alpha)
        zloc = surf.bnd_z_chord[surf.ndiv] - 0.5*surf.kinem.u*dt*sin(surf.kinem.alpha)
    else
        xloc = surf.bnd_x_chord[surf.ndiv]+(1. /3.)*(field.tev[ntev].x - surf.bnd_x_chord[surf.ndiv])

        zloc = surf.bnd_z_chord[surf.ndiv]+(1. /3.)*(field.tev[ntev].z - surf.bnd_z_chord[surf.ndiv])
    end
    push!(field.tev,TwoDVort(xloc,zloc,0.,0.02*surf.c,0.,0.))
    return field
end




# function update_RHS(surf::TwoDSurfThick, curfield::TwoDFlowField)

#     wlz = 0.5*(surf.wind_u.*cos(surf.kinem.alpha) .+ surf.uind_u.*sin(surf.kinem.alpha) .+
#                surf.wind_l*cos(surf.kinem.alpha) .+ surf.uind_l*sin(surf.kinem.alpha))

#     wtz = 0.5*(surf.wind_u.*cos(surf.kinem.alpha) .+ surf.uind_u.*sin(surf.kinem.alpha) .-
#                    surf.wind_l.*cos(surf.kinem.alpha) .- surf.uind_l.*sin(surf.kinem.alpha))

#     wtx = 0.5*(surf.uind_u.*cos(surf.kinem.alpha) .- surf.wind_u.*sin(surf.kinem.alpha) .+
#                surf.uind_l.*cos(surf.kinem.alpha) .- surf.wind_l.*sin(surf.kinem.alpha))

#     wlx = 0.5*(surf.uind_u.*cos(surf.kinem.alpha) .- surf.wind_u.*sin(surf.kinem.alpha) .-
#                surf.uind_l.*cos(surf.kinem.alpha) .+ surf.wind_l.*sin(surf.kinem.alpha))
    
#     for i = 2:surf.ndiv-1
#         #RHS for lifting equation
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
    
#     #surf.RHS[2*surf.ndiv-3] = bound_circ
# end

# function update_thickRHS(surf::TwoDSurfThick, curfield::TwoDFlowField)
#     for i = 1:surf.ndiv-1
#         #RHS for lifting equation

#         wlz = 0.5*(surf.wind_u[i]*cos(surf.kinem.alpha) + surf.uind_u[i]*sin(surf.kinem.alpha) +
#                    surf.wind_l[i]*cos(surf.kinem.alpha) + surf.uind_l[i]*sin(surf.kinem.alpha))

#         wtz = 0.5*(surf.wind_u[i]*cos(surf.kinem.alpha) + surf.uind_u[i]*sin(surf.kinem.alpha) -
#                    surf.wind_l[i]*cos(surf.kinem.alpha) - surf.uind_l[i]*sin(surf.kinem.alpha))

#         wtx = 0.5*(surf.uind_u[i]*cos(surf.kinem.alpha) - surf.wind_u[i]*sin(surf.kinem.alpha) +
#                    surf.uind_l[i]*cos(surf.kinem.alpha) - surf.wind_l[i]*sin(surf.kinem.alpha))

#         wlx = 0.5*(surf.uind_u[i]*cos(surf.kinem.alpha) - surf.wind_u[i]*sin(surf.kinem.alpha) -
#                    surf.uind_l[i]*cos(surf.kinem.alpha) + surf.wind_l[i]*sin(surf.kinem.alpha))

#         surf.RHS[i] = -(surf.kinem.u + curfield.u[1])*sin(surf.kinem.alpha)/surf.uref -
#             surf.kinem.alphadot*(surf.x[i] - surf.pvt*surf.c)/surf.uref +
#             (surf.kinem.hdot - curfield.w[1])*cos(surf.kinem.alpha)/surf.uref - wlz/surf.uref +
#             surf.cam_slope[i]*((surf.kinem.u + curfield.u[1])*cos(surf.kinem.alpha) +
#                                (surf.kinem.hdot - curfield.w[1])*sin(surf.kinem.alpha) + wtx -
#                                surf.kinem.alphadot*surf.cam[i]) +
#                                surf.thick_slope[i]*(wlx - surf.kinem.alphadot*surf.thick[i])
#         surf.RHS[surf.ndiv+i-1] = surf.cam_slope[i]*(wlx - surf.kinem.alphadot*surf.thick[i]) +
#             surf.thick_slope[i]*((surf.kinem.u + curfield.u[1])*cos(surf.kinem.alpha) +
#                                  (surf.kinem.hdot - curfield.w[1])*sin(surf.kinem.alpha) + wtx -
#                                  surf.kinem.alphadot*surf.cam[i]) - wtz
#     end

#     #RHS for Kelvin condition (negative strength of all previously shed vortices)
#     surf.RHS[2*surf.ndiv-1] = -100*(sum(map(q->q.s, curfield.tev)) + sum(map(q->q.s, curfield.lev)))/(surf.uref*surf.c)
    
#     wzu = surf.wind_u[surf.ndiv]*cos(surf.kinem.alpha) + surf.uind_u[surf.ndiv]*sin(surf.kinem.alpha)
#     wzl = surf.wind_l[surf.ndiv]*cos(surf.kinem.alpha) + surf.uind_l[surf.ndiv]*sin(surf.kinem.alpha)
#     wxu = surf.uind_u[surf.ndiv]*cos(surf.kinem.alpha) - surf.wind_u[surf.ndiv]*sin(surf.kinem.alpha)
#     wxl = surf.uind_l[surf.ndiv]*cos(surf.kinem.alpha) - surf.wind_l[surf.ndiv]*sin(surf.kinem.alpha)

#     #Kutta condition
#     surf.RHS[2*surf.ndiv] = ((surf.kinem.u + curfield.u[1])*cos(surf.kinem.alpha) + (surf.kinem.hdot - curfield.w[1])*sin(surf.kinem.alpha) + surf.kinem.alphadot*(surf.cam[end] - surf.thick[end]) + wxl + (surf.cam_slope[end] - surf.thick_slope[end])*((surf.kinem.u + curfield.u[1])*sin(surf.kinem.alpha) - (surf.kinem.hdot - curfield.w[1])*cos(surf.kinem.alpha) + surf.kinem.alphadot*(surf.x[surf.ndiv] - surf.pvt*surf.c) + wzl))/sqrt(1. + (surf.cam_slope[surf.ndiv] - surf.thick_slope[surf.ndiv])^2) - ((surf.kinem.u + curfield.u[1])*cos(surf.kinem.alpha) + (surf.kinem.hdot - curfield.w[1])*sin(surf.kinem.alpha) + surf.kinem.alphadot*(surf.cam[end] + surf.thick[end]) + wxu + (surf.cam_slope[end] + surf.thick_slope[end])*((surf.kinem.u + curfield.u[1])*sin(surf.kinem.alpha) - (surf.kinem.hdot - curfield.w[1])*cos(surf.kinem.alpha) + surf.kinem.alphadot*(surf.x[surf.ndiv] - surf.pvt*surf.c) + wzu))/sqrt(1. + (surf.cam_slope[surf.ndiv] + surf.thick_slope[surf.ndiv])^2) 

#     i = surf.ndiv
#     wlz = 0.5*(surf.wind_u[i]*cos(surf.kinem.alpha) + surf.uind_u[i]*sin(surf.kinem.alpha) +
#                        surf.wind_l[i]*cos(surf.kinem.alpha) + surf.uind_l[i]*sin(surf.kinem.alpha))
    
#     wtz = 0.5*(surf.wind_u[i]*cos(surf.kinem.alpha) + surf.uind_u[i]*sin(surf.kinem.alpha) -
#                    surf.wind_l[i]*cos(surf.kinem.alpha) - surf.uind_l[i]*sin(surf.kinem.alpha))
    
#     wtx = 0.5*(surf.uind_u[i]*cos(surf.kinem.alpha) - surf.wind_u[i]*sin(surf.kinem.alpha) +
#                    surf.uind_l[i]*cos(surf.kinem.alpha) - surf.wind_l[i]*sin(surf.kinem.alpha))
    
#     wlx = 0.5*(surf.uind_u[i]*cos(surf.kinem.alpha) - surf.wind_u[i]*sin(surf.kinem.alpha) -
#                    surf.uind_l[i]*cos(surf.kinem.alpha) + surf.wind_l[i]*sin(surf.kinem.alpha))


#     surf.RHS[2*surf.ndiv+1] = surf.cam_slope[i]*(wlx - surf.kinem.alphadot*surf.thick[i]) +
#         surf.thick_slope[i]*((surf.kinem.u + curfield.u[1])*cos(surf.kinem.alpha) +
#                              (surf.kinem.hdot - curfield.w[1])*sin(surf.kinem.alpha) + wtx -
#                              surf.kinem.alphadot*surf.cam[i]) - wtz
    

    
#     #surf.RHS[2*surf.ndiv-3] = pi*(surf.a0prev[1] + 0.5*surf.aprev[1])

#         #RHS for LE velocity condition
#     #surf.RHS[2*surf.ndiv-1] = sin(surf.kinem.alpha) - surf.kinem.hdot*cos(surf.kinem.alpha)/surf.uref -
#     #surf.kinem.alphadot*surf.pvt*surf.c/surf.uref +
#     #   (surf.wind_u[1]*cos(surf.kinem.alpha) + surf.uind_u[1]*sin(surf.kinem.alpha))/surf.uref


#    # surf.RHS[2*surf.ndiv-2] = -(surf.kinem.u + curfield.u[1])*cos(surf.kinem.alpha) - (surf.kinem.hdot - curfield.w[1])*sin(surf.kinem.alpha)

#     #surf.RHS[2*surf.ndiv-1] = -(surf.kinem.u + curfield.u[1])*cos(surf.kinem.alpha) - (surf.kinem.hdot - curfield.w[1])*sin(surf.kinem.alpha)

# end

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


# function wkg_update_thickLHS(surf::TwoDSurfThick, curfield::TwoDFlowField, dt::Float64, vcore::Float64)
#     #Calculate the missing column in LHS that depends on last shed vortex location

#     ntev = length(curfield.tev)

#     if ntev == 0
#         xloc_tev = surf.bnd_x_chord[surf.ndiv] + 0.5*surf.kinem.u*dt
#         zloc_tev = surf.bnd_z_chord[surf.ndiv]
#     else
#         xloc_tev = surf.bnd_x_chord[surf.ndiv] + (1. /3.)*(curfield.tev[ntev].x - surf.bnd_x_chord[surf.ndiv])
#         zloc_tev = surf.bnd_z_chord[surf.ndiv] + (1. /3.)*(curfield.tev[ntev].z - surf.bnd_z_chord[surf.ndiv])
#     end

#     dummyvort = TwoDVort(xloc_tev, zloc_tev, 1., vcore, 0., 0.)

#     uu, wu = ind_vel([dummyvort], surf.bnd_x_u, surf.bnd_z_u)
#     ul, wl = ind_vel([dummyvort], surf.bnd_x_l, surf.bnd_z_l)

#     for i = 2:surf.ndiv-1


#         #Sweep all rows (corresponding to ndiv) for lifting equation

#         #Sweep columns for aterms
#         for n = 1:surf.naterm
#            surf.LHS[i-1,n] = cos(n*surf.theta[i]) - surf.thick_slope[i]*sin(n*surf.theta[i])
#         end

#         #Sweep columns for bterm
#         for n = 1:surf.naterm
#            surf.LHS[i-1,n+surf.naterm] = surf.cam_slope[i]*cos(n*surf.theta[i])
#         end

#         #TEV term must be updated in the loop after its location is known
#         #Sweep all rows (corresponding to ndiv) for nonlifting equation

#         for n = 1:surf.naterm
#            surf.LHS[surf.ndiv+i-3,n]  = -surf.cam_slope[i]*sin(n*surf.theta[i])
#         end
#         for n = 1:surf.naterm
#            surf.LHS[surf.ndiv+i-3,surf.naterm+n] = sin(n*surf.theta[i]) + surf.thick_slope[i]*cos(n*surf.theta[i])
#         end


#         wlz = 0.5*(wu[i]*cos(surf.kinem.alpha) + uu[i]*sin(surf.kinem.alpha) +
#                    wl[i]*cos(surf.kinem.alpha) + ul[i]*sin(surf.kinem.alpha))

#         wtz = 0.5*(wu[i]*cos(surf.kinem.alpha) + uu[i]*sin(surf.kinem.alpha) -
#                    wl[i]*cos(surf.kinem.alpha) - ul[i]*sin(surf.kinem.alpha))

#         wtx = 0.5*(uu[i]*cos(surf.kinem.alpha) - wu[i]*sin(surf.kinem.alpha) +
#                    ul[i]*cos(surf.kinem.alpha) - wl[i]*sin(surf.kinem.alpha))

#         wlx = 0.5*(uu[i]*cos(surf.kinem.alpha) - wu[i]*sin(surf.kinem.alpha) -
#                    ul[i]*cos(surf.kinem.alpha) + wl[i]*sin(surf.kinem.alpha))

#         surf.LHS[i-1,1+2*surf.naterm] = -surf.cam_slope[i]*wtx -
#             surf.thick_slope[i]*wlx + wlz

#         surf.LHS[surf.ndiv+i-3,1+2*surf.naterm] = -surf.cam_slope[i]*wlx -
#             surf.thick_slope[i]*wtx + wtz

#     end

#     #surf.LHS[2*surf.ndiv-1,2+2*surf.naterm] = -wu[1]*cos(surf.kinem.alpha) -
#     #uu[1]*sin(surf.kinem.alpha)


#     return surf, xloc_tev, zloc_tev
# end

# function wkg_update_thickRHS(surf::TwoDSurfThick, curfield::TwoDFlowField)
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
#     surf.RHS[2*surf.ndiv-3] = -100*(sum(map(q->q.s, curfield.tev)) + sum(map(q->q.s, curfield.lev)))/(surf.uref*surf.c)

#     #Kutta condition
#     surf.RHS[2*surf.ndiv-2] = 0.
#     #surf.RHS[2*surf.ndiv-3] = pi*(surf.a0prev[1] + 0.5*surf.aprev[1])

#         #RHS for LE velocity condition
#     #surf.RHS[2*surf.ndiv-1] = sin(surf.kinem.alpha) - surf.kinem.hdot*cos(surf.kinem.alpha)/surf.uref -
#     #surf.kinem.alphadot*surf.pvt*surf.c/surf.uref +
#     #   (surf.wind_u[1]*cos(surf.kinem.alpha) + surf.uind_u[1]*sin(surf.kinem.alpha))/surf.uref


#    # surf.RHS[2*surf.ndiv-2] = -(surf.kinem.u + curfield.u[1])*cos(surf.kinem.alpha) - (surf.kinem.hdot - curfield.w[1])*sin(surf.kinem.alpha)

#     #surf.RHS[2*surf.ndiv-1] = -(surf.kinem.u + curfield.u[1])*cos(surf.kinem.alpha) - (surf.kinem.hdot - curfield.w[1])*sin(surf.kinem.alpha)

# end







# function update_thickLHS(surf::TwoDSurfThick, curfield::TwoDFlowField, dt::Float64, vcore::Float64)
#     #Calculate the missing column in LHS that depends on last shed vortex location

#     ntev = length(curfield.tev)

#     if ntev == 0
#         xloc_tev = surf.bnd_x_chord[surf.ndiv] + 0.5*surf.kinem.u*dt
#         zloc_tev = surf.bnd_z_chord[surf.ndiv]
#     else
#         xloc_tev = surf.bnd_x_chord[surf.ndiv] + (1. /3.)*(curfield.tev[ntev].x - surf.bnd_x_chord[surf.ndiv])
#         zloc_tev = surf.bnd_z_chord[surf.ndiv] + (1. /3.)*(curfield.tev[ntev].z - surf.bnd_z_chord[surf.ndiv])
#     end

#     dummyvort = TwoDVort(xloc_tev, zloc_tev, 1., vcore, 0., 0.)

#     uu, wu = ind_vel([dummyvort], surf.bnd_x_u, surf.bnd_z_u)
#     ul, wl = ind_vel([dummyvort], surf.bnd_x_l, surf.bnd_z_l)

#     surf.LHS[surf.ndiv*2-3,:] .= 0.

#     for i = 2:surf.ndiv-1

#         #Sweep all rows (corresponding to ndiv) for lifting equation

#         #Sweep columns for aterms
#         for n = 1:surf.naterm
#             surf.LHS[i-1,n] = cos(n*surf.theta[i]) - surf.thick_slope[i]*sin(n*surf.theta[i])
#         end
        
#         #Sweep columns for bterm
#         for n = 1:surf.naterm
#            surf.LHS[i-1,n+surf.naterm] = surf.cam_slope[i]*cos(n*surf.theta[i])
#         end

#         #TEV term must be updated in the loop after its location is known
#         #Sweep all rows (corresponding to ndiv) for nonlifting equation

#         for n = 1:surf.naterm
#            surf.LHS[surf.ndiv+i-3,n]  = -surf.cam_slope[i]*sin(n*surf.theta[i])
#         end
#         for n = 1:surf.naterm
#             surf.LHS[surf.ndiv+i-3,surf.naterm+n] = sin(n*surf.theta[i]) + surf.thick_slope[i]*cos(n*surf.theta[i])
#         end
        
        
#         wlz = 0.5*(wu[i]*cos(surf.kinem.alpha) + uu[i]*sin(surf.kinem.alpha) +
#                    wl[i]*cos(surf.kinem.alpha) + ul[i]*sin(surf.kinem.alpha))
        
#         wtz = 0.5*(wu[i]*cos(surf.kinem.alpha) + uu[i]*sin(surf.kinem.alpha) -
#                    wl[i]*cos(surf.kinem.alpha) - ul[i]*sin(surf.kinem.alpha))

#         wtx = 0.5*(uu[i]*cos(surf.kinem.alpha) - wu[i]*sin(surf.kinem.alpha) +
#                    ul[i]*cos(surf.kinem.alpha) - wl[i]*sin(surf.kinem.alpha))

#         wlx = 0.5*(uu[i]*cos(surf.kinem.alpha) - wu[i]*sin(surf.kinem.alpha) -
#                    ul[i]*cos(surf.kinem.alpha) + wl[i]*sin(surf.kinem.alpha))

#         surf.LHS[i-1,1+2*surf.naterm] = -surf.cam_slope[i]*wtx -
#             surf.thick_slope[i]*wlx + wlz

#         surf.LHS[surf.ndiv+i-3,1+2*surf.naterm] = -surf.cam_slope[i]*wlx -
#             surf.thick_slope[i]*wtx + wtz
#     end
    
#     for n = 1:surf.naterm
#         s1 = 0. ; s2 = 0. ; s3 = 0. ; s4 = 0. ; s5 = 0. ; s6 = 0.
#         for i = 2:surf.ndiv
#             #dx is ds
#             dx = sqrt((surf.x[i] - surf.x[i-1])^2 + (surf.cam[i] + surf.thick[i] - surf.cam[i-1] - surf.thick[i-1])^2)
#             den_pl = sqrt(1 + (surf.thick_slope[i] + surf.cam_slope[i])^2)
#             den_mi = sqrt(1 + (-surf.thick_slope[i] + surf.cam_slope[i])^2)
#             den_pl_p = sqrt(1 + (surf.thick_slope[i-1] + surf.cam_slope[i-1])^2)
#             den_mi_p = sqrt(1 + (-surf.thick_slope[i-1] + surf.cam_slope[i-1])^2)
            
#             s1 += 0.5*dx*(sin(n*surf.theta[i])/den_pl + sin(n*surf.theta[i])/den_mi + sin(n*surf.theta[i-1])/den_pl_p + sin(n*surf.theta[i-1])/den_mi_p)
#             s2 += 0.5*dx*(surf.cam_slope[i]*cos(n*surf.theta[i])/den_pl - surf.cam_slope[i]*cos(n*surf.theta[i])/den_mi + surf.cam_slope[i-1]*cos(n*surf.theta[i-1])/den_pl_p - surf.cam_slope[i-1]*cos(n*surf.theta[i-1])/den_mi_p)
#             s3 += 0.5*dx*(surf.thick_slope[i]*cos(n*surf.theta[i])/den_pl + surf.thick_slope[i]*cos(n*surf.theta[i])/den_mi + surf.thick_slope[i-1]*cos(n*surf.theta[i-1])/den_pl_p + surf.thick_slope[i-1]*cos(n*surf.theta[i-1])/den_mi_p)

#             dx = sqrt((surf.x[i] - surf.x[i-1])^2 + (surf.cam[i] - surf.thick[i] - surf.cam[i-1] + surf.thick[i-1])^2)

#             s4 -= 0.5*dx*(cos(n*surf.theta[i])/den_pl - cos(n*surf.theta[i])/den_mi + cos(n*surf.theta[i-1])/den_pl_p - cos(n*surf.theta[i-1])/den_mi_p)
#             s5 += 0.5*dx*(surf.cam_slope[i]*sin(n*surf.theta[i])/den_pl + surf.cam_slope[i]*sin(n*surf.theta[i])/den_mi + surf.cam_slope[i-1]*sin(n*surf.theta[i-1])/den_pl_p + surf.cam_slope[i-1]*sin(n*surf.theta[i-1])/den_mi_p)
#             s6 += 0.5*dx*(surf.thick_slope[i]*sin(n*surf.theta[i])/den_pl - surf.thick_slope[i]*sin(n*surf.theta[i])/den_mi + surf.thick_slope[i-1]*sin(n*surf.theta[i-1])/den_pl_p - surf.thick_slope[i-1]*sin(n*surf.theta[i-1])/den_mi_p)
#         end
#         surf.LHS[2*surf.ndiv-3, n] = surf.uref*(s1 + s2 + s3)
#         surf.LHS[2*surf.ndiv-3, n+surf.naterm] = surf.uref*(s4 + s5 + s6)
#         surf.LHS[2*surf.ndiv-3, 2*surf.naterm+1] = 1.
#     end 
        
#     return surf, xloc_tev, zloc_tev
# end

    
    # function update_thickLHS(surf::TwoDSurfThick, curfield::TwoDFlowField, dt::Float64, vcore::Float64)
    #     #Calculate the missing column in LHS that depends on last shed vortex location
    
    #     ntev = length(curfield.tev)
    
    #     if ntev == 0
    #         xloc_tev = surf.bnd_x_chord[surf.ndiv] + 0.5*surf.kinem.u*dt
#         zloc_tev = surf.bnd_z_chord[surf.ndiv]
    #     else
    #         xloc_tev = surf.bnd_x_chord[surf.ndiv] + (1. /3.)*(curfield.tev[ntev].x - surf.bnd_x_chord[surf.ndiv])
    #         zloc_tev = surf.bnd_z_chord[surf.ndiv] + (1. /3.)*(curfield.tev[ntev].z - surf.bnd_z_chord[surf.ndiv])
#     end

#     dummyvort = TwoDVort(xloc_tev, zloc_tev, 1., vcore, 0., 0.)

#     uu, wu = ind_vel([dummyvort], surf.bnd_x_u, surf.bnd_z_u)
#     ul, wl = ind_vel([dummyvort], surf.bnd_x_l, surf.bnd_z_l)

#     for i = 1:surf.ndiv-1

#         #Sweep all rows (corresponding to ndiv) for lifting equation

#         #Sweep columns for aterms
#         for n = 1:surf.naterm
#             surf.LHS[i,n] = cos(n*surf.theta[i]) - surf.thick_slope[i]*sin(n*surf.theta[i])
#         end

#         #Sweep columns for bterm
#         for n = 1:surf.naterm
#             surf.LHS[i,n+surf.naterm] = surf.cam_slope[i]*cos(n*surf.theta[i])
#         end
        
#         #TEV term must be updated in the loop after its location is known
#         #Sweep all rows (corresponding to ndiv) for nonlifting equation
        
#         for n = 1:surf.naterm
#             surf.LHS[surf.ndiv+i-1,n]  = -surf.cam_slope[i]*sin(n*surf.theta[i])
#         end
        
#         for n = 1:surf.naterm
#             surf.LHS[surf.ndiv+i-1,surf.naterm+n] = sin(n*surf.theta[i]) + surf.thick_slope[i]*cos(n*surf.theta[i])
#         end

        
#         wlz = 0.5*(wu[i]*cos(surf.kinem.alpha) + uu[i]*sin(surf.kinem.alpha) +
#                    wl[i]*cos(surf.kinem.alpha) + ul[i]*sin(surf.kinem.alpha))
        
#         wtz = 0.5*(wu[i]*cos(surf.kinem.alpha) + uu[i]*sin(surf.kinem.alpha) -
#                    wl[i]*cos(surf.kinem.alpha) - ul[i]*sin(surf.kinem.alpha))
        
#         wtx = 0.5*(uu[i]*cos(surf.kinem.alpha) - wu[i]*sin(surf.kinem.alpha) +
#                    ul[i]*cos(surf.kinem.alpha) - wl[i]*sin(surf.kinem.alpha))
        
#         wlx = 0.5*(uu[i]*cos(surf.kinem.alpha) - wu[i]*sin(surf.kinem.alpha) -
#                    ul[i]*cos(surf.kinem.alpha) + wl[i]*sin(surf.kinem.alpha))
        
#         surf.LHS[i,1+2*surf.naterm] = -surf.cam_slope[i]*wtx -
#             surf.thick_slope[i]*wlx + wlz
        
#         surf.LHS[surf.ndiv+i-1,1+2*surf.naterm] = -surf.cam_slope[i]*wlx -
#             surf.thick_slope[i]*wtx + wtz        
#     end
    
#     surf.LHS[2*surf.ndiv,:] .= 0.
    
#     for n = 1:surf.naterm
#         surf.LHS[2*surf.ndiv, n] = (surf.cam_slope[surf.ndiv] + surf.thick_slope[surf.ndiv])*(-1)^n/sqrt(1. + (surf.cam_slope[surf.ndiv] + surf.thick_slope[surf.ndiv])^2) - (surf.cam_slope[surf.ndiv] - surf.thick_slope[surf.ndiv])*(-1)^n/sqrt(1. + (surf.cam_slope[surf.ndiv] - surf.thick_slope[surf.ndiv])^2)
#         surf.LHS[2*surf.ndiv, n+surf.naterm] = (-1)^n/sqrt(1. + (surf.cam_slope[surf.ndiv] - surf.thick_slope[surf.ndiv])^2) - (-1)^n/sqrt(1. + (surf.cam_slope[surf.ndiv] + surf.thick_slope[surf.ndiv])^2)
#     end
    
#     wzu = wu[surf.ndiv]*cos(surf.kinem.alpha) + uu[surf.ndiv]*sin(surf.kinem.alpha)
#     wzl = wl[surf.ndiv]*cos(surf.kinem.alpha) + ul[surf.ndiv]*sin(surf.kinem.alpha)     
#     wxu = uu[surf.ndiv]*cos(surf.kinem.alpha) - wu[surf.ndiv]*sin(surf.kinem.alpha)
#     wxl = ul[surf.ndiv]*cos(surf.kinem.alpha) - wl[surf.ndiv]*sin(surf.kinem.alpha)
    
#     #Add BCs at TE and LE
    
#     surf.LHS[2*surf.ndiv, 2*surf.naterm+1] = (wxu + (surf.cam_slope[surf.ndiv] + surf.thick_slope[surf.ndiv])*wzu)/sqrt(1. + (surf.cam_slope[surf.ndiv] + surf.thick_slope[surf.ndiv])^2) - (wxl + (surf.cam_slope[surf.ndiv] - surf.thick_slope[surf.ndiv])*wzl)/sqrt(1. + (surf.cam_slope[surf.ndiv] - surf.thick_slope[surf.ndiv])^2)  
    
#     # #     #surf.LHS[2*surf.ndiv-1,2+2*surf.naterm] = -wu[1]*cos(surf.kinem.alpha) -
# # #     #uu[1]*sin(surf.kinem.alpha)

#     #Nonlifting equations at the trailing edge
#     i = surf.ndiv
#     for n = 1:surf.naterm
#         surf.LHS[2*surf.ndiv+1,n]  = -surf.cam_slope[i]*sin(n*surf.theta[i])
#     end
    
#     for n = 1:surf.naterm
#         surf.LHS[2*surf.ndiv+1,surf.naterm+n] = sin(n*surf.theta[i]) + surf.thick_slope[i]*cos(n*surf.theta[i])
#     end
    
#     wlz = 0.5*(wu[i]*cos(surf.kinem.alpha) + uu[i]*sin(surf.kinem.alpha) +
#                        wl[i]*cos(surf.kinem.alpha) + ul[i]*sin(surf.kinem.alpha))
    
#     wtz = 0.5*(wu[i]*cos(surf.kinem.alpha) + uu[i]*sin(surf.kinem.alpha) -
#                    wl[i]*cos(surf.kinem.alpha) - ul[i]*sin(surf.kinem.alpha))
    
#     wtx = 0.5*(uu[i]*cos(surf.kinem.alpha) - wu[i]*sin(surf.kinem.alpha) +
#                    ul[i]*cos(surf.kinem.alpha) - wl[i]*sin(surf.kinem.alpha))
    
#     wlx = 0.5*(uu[i]*cos(surf.kinem.alpha) - wu[i]*sin(surf.kinem.alpha) -
#                    ul[i]*cos(surf.kinem.alpha) + wl[i]*sin(surf.kinem.alpha))
    
#     surf.LHS[2*surf.ndiv+1,1+2*surf.naterm] = -surf.cam_slope[i]*wlx -
#             surf.thick_slope[i]*wtx + wtz        

#     return surf, xloc_tev, zloc_tev
# end


# function update_thickLHS_kutta(surf::TwoDSurfThick, curfield::TwoDFlowField, dt::Float64, vcore::Float64)
#     #Calculate the missing column in LHS that depends on last shed vortex location

#     ntev = length(curfield.tev)

#     if ntev == 0
#         xloc_tev = surf.bnd_x_l[surf.ndiv] + 0.5*surf.kinem.u*dt
#         zloc_tev = surf.bnd_z_l[surf.ndiv]
#     else
#         xloc_tev = surf.bnd_x_l[surf.ndiv] + (1. /3.)*(curfield.tev[ntev].x - surf.bnd_x_l[surf.ndiv])
#         zloc_tev = surf.bnd_z_l[surf.ndiv] + (1. /3.)*(curfield.tev[ntev].z - surf.bnd_z_l[surf.ndiv])
#     end

#     dummyvort = TwoDVort(xloc_tev, zloc_tev, 1., vcore, 0., 0.)
    
#     uu, wu = ind_vel([dummyvort], surf.bnd_x_u, surf.bnd_z_u)
#     ul, wl = ind_vel([dummyvort], surf.bnd_x_l, surf.bnd_z_l)

#     for i = 2:surf.ndiv-1


#         #Sweep all rows (corresponding to ndiv) for lifting equation

#         #Sweep columns for aterms
#         for n = 1:surf.naterm
#            surf.LHS[i-1,n] = cos(n*surf.theta[i]) - surf.thick_slope[i]*sin(n*surf.theta[i])
#         end

#         #Sweep columns for bterm
#         for n = 1:surf.naterm
#            surf.LHS[i-1,n+surf.naterm] = surf.cam_slope[i]*cos(n*surf.theta[i])
#         end

#         #TEV term must be updated in the loop after its location is known
#         #Sweep all rows (corresponding to ndiv) for nonlifting equation

#         for n = 1:surf.naterm
#            surf.LHS[surf.ndiv+i-3,n]  = -surf.cam_slope[i]*sin(n*surf.theta[i])
#         end
#         for n = 1:surf.naterm
#            surf.LHS[surf.ndiv+i-3,surf.naterm+n] = sin(n*surf.theta[i]) + surf.thick_slope[i]*cos(n*surf.theta[i])
#         end


#         wlz = 0.5*(wu[i]*cos(surf.kinem.alpha) + uu[i]*sin(surf.kinem.alpha) +
#                    wl[i]*cos(surf.kinem.alpha) + ul[i]*sin(surf.kinem.alpha))

#         wtz = 0.5*(wu[i]*cos(surf.kinem.alpha) + uu[i]*sin(surf.kinem.alpha) -
#                    wl[i]*cos(surf.kinem.alpha) - ul[i]*sin(surf.kinem.alpha))

#         wtx = 0.5*(uu[i]*cos(surf.kinem.alpha) - wu[i]*sin(surf.kinem.alpha) +
#                    ul[i]*cos(surf.kinem.alpha) - wl[i]*sin(surf.kinem.alpha))

#         wlx = 0.5*(uu[i]*cos(surf.kinem.alpha) - wu[i]*sin(surf.kinem.alpha) -
#                    ul[i]*cos(surf.kinem.alpha) + wl[i]*sin(surf.kinem.alpha))

#         surf.LHS[i-1,1+2*surf.naterm] = -surf.cam_slope[i]*wtx -
#             surf.thick_slope[i]*wlx + wlz

#         surf.LHS[surf.ndiv+i-3,1+2*surf.naterm] = -surf.cam_slope[i]*wlx -
#             surf.thick_slope[i]*wtx + wtz

#     end

#     #surf.LHS[2*surf.ndiv-1,2+2*surf.naterm] = -wu[1]*cos(surf.kinem.alpha) -
#     #uu[1]*sin(surf.kinem.alpha)


#     return surf, xloc_tev, zloc_tev
# end

# function update_thickLHS(surf::TwoDSurfThick, curfield::TwoDFlowField, dt::Float64, vcore::Float64, xloc_tev::Float64, zloc_tev::Float64)
#     #Calculate the missing column in LHS that depends on last shed vortex location

#     ntev = length(curfield.tev)

#     dummyvort = TwoDVort(xloc_tev, zloc_tev, 1., vcore, 0., 0.)

#     uu, wu = ind_vel([dummyvort], surf.bnd_x_u, surf.bnd_z_u)
#     ul, wl = ind_vel([dummyvort], surf.bnd_x_l, surf.bnd_z_l)

#     for i = 2:surf.ndiv-1


#         #Sweep all rows (corresponding to ndiv) for lifting equation

#         #Sweep columns for aterms
#         for n = 1:surf.naterm
#            surf.LHS[i-1,n] = cos(n*surf.theta[i]) - surf.thick_slope[i]*sin(n*surf.theta[i])
#         end

#         #Sweep columns for bterm
#         for n = 1:surf.naterm
#            surf.LHS[i-1,n+surf.naterm] = surf.cam_slope[i]*cos(n*surf.theta[i])
#         end

#         #TEV term must be updated in the loop after its location is known
#         #Sweep all rows (corresponding to ndiv) for nonlifting equation

#         for n = 1:surf.naterm
#            surf.LHS[surf.ndiv+i-3,n]  = -surf.cam_slope[i]*sin(n*surf.theta[i])
#         end
#         for n = 1:surf.naterm
#            surf.LHS[surf.ndiv+i-3,surf.naterm+n] = sin(n*surf.theta[i]) + surf.thick_slope[i]*cos(n*surf.theta[i])
#         end


#         wlz = 0.5*(wu[i]*cos(surf.kinem.alpha) + uu[i]*sin(surf.kinem.alpha) +
#                    wl[i]*cos(surf.kinem.alpha) + ul[i]*sin(surf.kinem.alpha))

#         wtz = 0.5*(wu[i]*cos(surf.kinem.alpha) + uu[i]*sin(surf.kinem.alpha) -
#                    wl[i]*cos(surf.kinem.alpha) - ul[i]*sin(surf.kinem.alpha))

#         wtx = 0.5*(uu[i]*cos(surf.kinem.alpha) - wu[i]*sin(surf.kinem.alpha) +
#                    ul[i]*cos(surf.kinem.alpha) - wl[i]*sin(surf.kinem.alpha))

#         wlx = 0.5*(uu[i]*cos(surf.kinem.alpha) - wu[i]*sin(surf.kinem.alpha) -
#                    ul[i]*cos(surf.kinem.alpha) + wl[i]*sin(surf.kinem.alpha))

#         surf.LHS[i-1,1+2*surf.naterm] = -surf.cam_slope[i]*wtx -
#             surf.thick_slope[i]*wlx + wlz

#         surf.LHS[surf.ndiv+i-3,1+2*surf.naterm] = -surf.cam_slope[i]*wlx -
#             surf.thick_slope[i]*wtx + wtz

#     end

#     for n = 1:surf.naterm
#         s1 = 0. ; s2 = 0. ; s3 = 0. ; s4 = 0. ; s5 = 0. ; s6 = 0.
#         for i = 2:surf.ndiv
#             dx = surf.x[i] - surf.x[i-1]
#             den_pl = sqrt(1 + (surf.thick_slope[i] + surf.cam_slope[i])^2)
#             den_mi = sqrt(1 + (-surf.thick_slope[i] + surf.cam_slope[i])^2)
#             den_pl_p = sqrt(1 + (surf.thick_slope[i-1] + surf.cam_slope[i-1])^2)
#             den_mi_p = sqrt(1 + (-surf.thick_slope[i-1] + surf.cam_slope[i-1])^2)
            
#             s1 += 0.5*dx*(sin(n*surf.theta[i])/den_pl + sin(n*surf.theta[i])/den_mi + sin(n*surf.theta[i-1])/den_pl_p + sin(n*surf.theta[i-1])/den_mi_p)
#             s2 += 0.5*dx*(surf.cam_slope[i]*cos(n*surf.theta[i])/den_pl - surf.cam_slope[i]*cos(n*surf.theta[i])/den_mi + surf.cam_slope[i-1]*cos(n*surf.theta[i-1])/den_pl_p - surf.cam_slope[i-1]*cos(n*surf.theta[i-1])/den_mi_p)
#             s3 += 0.5*dx*(surf.thick_slope[i]*cos(n*surf.theta[i])/den_pl + surf.thick_slope[i]*cos(n*surf.theta[i])/den_mi + surf.thick_slope[i-1]*cos(n*surf.theta[i-1])/den_pl_p + surf.thick_slope[i-1]*cos(n*surf.theta[i-1])/den_mi_p)
#             s4 -= 0.5*dx*(cos(n*surf.theta[i])/den_pl - cos(n*surf.theta[i])/den_mi + cos(n*surf.theta[i-1])/den_pl_p - cos(n*surf.theta[i-1])/den_mi_p)
#             s5 += 0.5*dx*(surf.cam_slope[i]*sin(n*surf.theta[i])/den_pl + surf.cam_slope[i]*sin(n*surf.theta[i])/den_mi + surf.cam_slope[i-1]*sin(n*surf.theta[i-1])/den_pl_p + surf.cam_slope[i-1]*sin(n*surf.theta[i-1])/den_mi_p)
#             s6 += 0.5*dx*(surf.thick_slope[i]*sin(n*surf.theta[i])/den_pl - surf.thick_slope[i]*sin(n*surf.theta[i])/den_mi + surf.thick_slope[i-1]*sin(n*surf.theta[i-1])/den_pl_p - surf.thick_slope[i-1]*sin(n*surf.theta[i-1])/den_mi_p)
#         end
#         surf.LHS[2*surf.ndiv-3, n] = surf.uref*(s1 + s2 + s3)
#         surf.LHS[2*surf.ndiv-3, n+surf.naterm] = surf.uref*(s4 + s5 + s6)
#         surf.LHS[2*surf.ndiv-3, 2*surf.naterm+1] = 1.
#     end 
#     return surf
# end
    


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
