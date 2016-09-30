#Function for calculating a simple 3D correction based on quasi-steady LLT

function simple_LLT(mat::Array{Float64,2}, surf::TwoDSurf, s::Float64, n_bterm::Int64, n_span::Int64)
    # s is the semi-span
    nsteps = length(mat[:,1])

    bnd_circ2d = mat[:,9]
    t = mat[:,1]
    alpha = mat[:,2]
    h = mat[:,3]
    hdot = zeros(nsteps)
    hdot[1] = 0
    hdot[2:end] = diff(h)./diff(t)
    u = mat[:,4]


    lhs = zeros(n_span,n_bterm)
    rhs = zeros(n_span)
    b_coeff = zeros(n_bterm)
    bcoeff_prev = zeros(n_bterm)
    bdot = zeros(n_bterm)
    sp_gam = zeros(nsteps,n_span)
    psi = zeros(n_span)
    dt = mat[2,1]-mat[1,1]
    cnc_f = zeros(nsteps)
    cnnc_f = zeros(nsteps)

    for i = 1:n_span
        psi[i] = (real(i)/(n_span+1))*pi
    end


    for i = 1:nsteps
        for j = 1:n_span
            for n = 1:n_bterm
                lhs[j,n] = sin(n*psi[j])*(sin(psi[j]) + (n*surf.c*pi/(4*s)))
            end
            rhs[j] = surf.c*pi*sin(psi[j])*bnd_circ2d[i]/(4*s)
        end
        bcoeff_prev = b_coeff
        b_coeff = \(lhs, rhs)
        bdot = (b_coeff - bcoeff_prev)/dt

        # for j = 1:n_span
        #     sp_gam[i,j] = 0
        #     for n = 1:n_bterm
        #         sp_gam[i,j] = sp_gam[i,j] + 4*s*u*b_coeff[n]*sin(n*psi[j])
        #     end
        # end

        sum_bcoeff = 0
        for n = 1:n_bterm
            if rem(n,2) != 0
                sum_bcoeff = sum_bcoeff + b_coeff[n]
            end
        end
        cnc_f[i] = -2*pi*(u[i]*cos(alpha[i])/surf.uref + hdot[i]*sin(alpha[i])/surf.uref)*(sum_bcoeff)
        sum_bdot = 0
        for n = 1:n_bterm
            if rem(n,2) != 0
                sum_bdot = sum_bdot + bdot[n]
            end
        end
        cnnc_f[i] = -(2*pi*surf.c/(surf.uref))*(3*sum_bdot/4)
    end
    return cnc_f, cnnc_f
end

# ---------------------------------------------------------------------------------------------
# Function for calculating a_2 and a_3 fourier coefficients
function update_a2a3adot(surf::TwoDSurf,dt)
    for ia = 2:3
        surf.aterm[ia] = trapz(surf.downwash.*cos(ia*surf.theta),surf.theta)
        surf.aterm[ia] = 2.*surf.aterm[ia]/(surf.uref*pi)
    end
    surf.a0dot[1] = (surf.a0[1] - surf.a0prev[1])/dt
    for ia = 1:3
        surf.adot[ia] = (surf.aterm[ia]-surf.aprev[ia])/dt
    end
    return surf
end

function update_a2a3adot(surf::TwoDSurf_2DOF,dt)
    for ia = 2:3
        surf.aterm[ia] = trapz(surf.downwash.*cos(ia*surf.theta),surf.theta)
        surf.aterm[ia] = 2.*surf.aterm[ia]/(surf.uref*pi)
    end
    surf.a0dot[1] = (surf.a0[1] - surf.a0prev[1])/dt
    for ia = 1:3
        surf.adot[ia] = (surf.aterm[ia]-surf.aprev[ia])/dt
    end
    return surf
end

function update_a2a3adot(surf::TwoDFreeSurf,dt)
    for ia = 2:3
        surf.aterm[ia] = trapz(surf.downwash.*cos(ia*surf.theta),surf.theta)
        surf.aterm[ia] = 2.*surf.aterm[ia]/(surf.uref*pi)
    end
    surf.a0dot[1] = (surf.a0[1] - surf.a0prev[1])/dt
    for ia = 1:3
        surf.adot[ia] = (surf.aterm[ia]-surf.aprev[ia])/dt
    end
    return surf
end

function update_a2a3adot(surf::TwoDSurfwFlap,dt)
    for ia = 2:3
        surf.aterm[ia] = trapz(surf.downwash.*cos(ia*surf.theta),surf.theta)
        surf.aterm[ia] = 2.*surf.aterm[ia]/(surf.uref*pi)
    end
    surf.a0dot[1] = (surf.a0[1] - surf.a0prev[1])/dt
    for ia = 1:3
        surf.adot[ia] = (surf.aterm[ia]-surf.aprev[ia])/dt
    end
    return surf
end

# ---------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------
# Function for updating the induced velocities u, w: eqns (2.12), (2.13) in Ramesh et al.
function update_indbound(surf::TwoDSurf, curfield::TwoDFlowField)
    surf.uind[1:surf.ndiv], surf.wind[1:surf.ndiv] = ind_vel([curfield.tev; curfield.lev], surf.bnd_x, surf.bnd_z)
    return surf
end

function update_indbound(surf::TwoDSurfwFlap, curfield::TwoDFlowField)
    surf.uind[1:surf.ndiv], surf.wind[1:surf.ndiv] = ind_vel([curfield.tev; curfield.lev], surf.bnd_x, surf.bnd_z)
    return surf
end

function update_indbound(surf::TwoDSurf_2DOF, curfield::TwoDFlowField)
    surf.uind[1:surf.ndiv], surf.wind[1:surf.ndiv] = ind_vel([curfield.tev; curfield.lev], surf.bnd_x, surf.bnd_z)
    return surf
end

function update_indbound(surf::TwoDFreeSurf, curfield::TwoDFlowField)
    surf.uind[1:surf.ndiv], surf.wind[1:surf.ndiv] = ind_vel([curfield.tev; curfield.lev; curfield.extv], surf.bnd_x, surf.bnd_z)
    return surf
end

# ---------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------
# Function for updating the downwash W: eqn (2.5) in Ramesh et al.
function update_downwash(surf::TwoDSurf, vels::Vector{Float64})
    for ib = 1:surf.ndiv
        surf.downwash[ib] = -(surf.kinem.u + vels[1])*sin(surf.kinem.alpha) - surf.uind[ib]*sin(surf.kinem.alpha) + (surf.kinem.hdot - vels[2])*cos(surf.kinem.alpha) - surf.wind[ib]*cos(surf.kinem.alpha) - surf.kinem.alphadot*(surf.x[ib] - surf.pvt*surf.c) + surf.cam_slope[ib]*(surf.uind[ib]*cos(surf.kinem.alpha) + (surf.kinem.u + vels[1])*cos(surf.kinem.alpha) + (surf.kinem.hdot - vels[2])*sin(surf.kinem.alpha) - surf.wind[ib]*sin(surf.kinem.alpha))
    end
    return surf
end
# ---------------------------------------------------------------------------------------------

function update_downwash(surf::TwoDSurf_2DOF)
    for ib = 1:surf.ndiv
        surf.downwash[ib] = -surf.kinem.u*sin(surf.kinem.alpha) - surf.uind[ib]*sin(surf.kinem.alpha) + surf.kinem.hdot*cos(surf.kinem.alpha) - surf.wind[ib]*cos(surf.kinem.alpha) - surf.kinem.alphadot*(surf.x[ib] - surf.pvt*surf.c) + surf.cam_slope[ib]*(surf.uind[ib]*cos(surf.kinem.alpha) + surf.kinem.u*cos(surf.kinem.alpha) + surf.kinem.hdot*sin(surf.kinem.alpha) - surf.wind[ib]*sin(surf.kinem.alpha))
    end
    return surf
end

function update_downwash(surf::TwoDFreeSurf)
    for ib = 1:surf.ndiv
        surf.downwash[ib] = -surf.kinem.u*sin(surf.kinem.alpha) - surf.uind[ib]*sin(surf.kinem.alpha) + surf.kinem.hdot*cos(surf.kinem.alpha) - surf.wind[ib]*cos(surf.kinem.alpha) - surf.kinem.alphadot*(surf.x[ib] - surf.pvt*surf.c) + surf.cam_slope[ib]*(surf.uind[ib]*cos(surf.kinem.alpha) + surf.kinem.u*cos(surf.kinem.alpha) + surf.kinem.hdot*sin(surf.kinem.alpha) - surf.wind[ib]*sin(surf.kinem.alpha))
    end
    return surf
end
# ---------------------------------------------------------------------------------------------

# Function for updating the downwash W: eqn (2.5) in Ramesh et al.
# Additional terms: alphadot*eta and time derivative of eta to account for deformations
# Added by Laura Merchant 2016
function update_downwash(surf::TwoDSurfwFlap)
    for ib = 1:surf.ndiv
        surf.downwash[ib] = -surf.kinem.u*sin(surf.kinem.alpha) - surf.uind[ib]*sin(surf.kinem.alpha) + surf.kinem.hdot*cos(surf.kinem.alpha) - surf.wind[ib]*cos(surf.kinem.alpha) - surf.kinem.alphadot*(surf.x[ib] - surf.pvt*surf.c) + surf.cam_tder[ib] + surf.cam_slope[ib]*(surf.uind[ib]*cos(surf.kinem.alpha) + surf.kinem.u*cos(surf.kinem.alpha) - surf.kinem.alphadot*surf.cam[ib] + surf.kinem.hdot*sin(surf.kinem.alpha) - surf.wind[ib]*sin(surf.kinem.alpha))
    end
    return surf
end
# ---------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------
# Function for a_0 and a_1 fourier coefficients
function update_a0anda1(surf::TwoDSurf)
    surf.a0[1] = trapz(surf.downwash,surf.theta)
    surf.aterm[1] = trapz(surf.downwash.*cos(surf.theta),surf.theta)
    surf.a0[1] = -surf.a0[1]/(surf.uref*pi)
    surf.aterm[1] = 2.*surf.aterm[1]/(surf.uref*pi)
    return surf
end

function update_a0anda1(surf::TwoDSurf_2DOF)
    surf.a0[1] = trapz(surf.downwash,surf.theta)
    surf.aterm[1] = trapz(surf.downwash.*cos(surf.theta),surf.theta)
    surf.a0[1] = -surf.a0[1]/(surf.uref*pi)
    surf.aterm[1] = 2.*surf.aterm[1]/(surf.uref*pi)
    return surf
end

function update_a0anda1(surf::TwoDFreeSurf)
    surf.a0[1] = trapz(surf.downwash,surf.theta)
    surf.aterm[1] = trapz(surf.downwash.*cos(surf.theta),surf.theta)
    surf.a0[1] = -surf.a0[1]/(surf.uref*pi)
    surf.aterm[1] = 2.*surf.aterm[1]/(surf.uref*pi)
    return surf
end

function update_a0anda1(surf::TwoDSurfwFlap)
    surf.a0[1] = trapz(surf.downwash,surf.theta)
    surf.aterm[1] = trapz(surf.downwash.*cos(surf.theta),surf.theta)
    surf.a0[1] = -surf.a0[1]/(surf.uref*pi)
    surf.aterm[1] = 2.*surf.aterm[1]/(surf.uref*pi)
    return surf
end
# ---------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------
# Function for calculating the fourier coefficients a_2 upwards to a_n
function update_a2toan(surf::TwoDSurf)
    for ia = 2:surf.naterm
        surf.aterm[ia] = trapz(surf.downwash.*cos(ia*surf.theta),surf.theta)
        surf.aterm[ia] = 2.*surf.aterm[ia]/(surf.uref*pi)
    end
    return surf
end

function update_a2toan(surf::TwoDSurf_2DOF)
    for ia = 2:surf.naterm
        surf.aterm[ia] = trapz(surf.downwash.*cos(ia*surf.theta),surf.theta)
        surf.aterm[ia] = 2.*surf.aterm[ia]/(surf.uref*pi)
    end
    return surf
end

function update_a2toan(surf::TwoDFreeSurf)
    for ia = 2:surf.naterm
        surf.aterm[ia] = trapz(surf.downwash.*cos(ia*surf.theta),surf.theta)
        surf.aterm[ia] = 2.*surf.aterm[ia]/(surf.uref*pi)
    end
    return surf
end

function update_a2toan(surf::TwoDSurfwFlap)
    for ia = 2:surf.naterm
        surf.aterm[ia] = trapz(surf.downwash.*cos(ia*surf.theta),surf.theta)
        surf.aterm[ia] = 2.*surf.aterm[ia]/(surf.uref*pi)
    end
    return surf
end
# ---------------------------------------------------------------------------------------------

#Function to update the external flowfield
function update_externalvel(curfield::TwoDFlowField, t)
    if (typeof(curfield.velX) == CosDef)
        curfield.u[1] = curfield.velX(t)
        curfield.w[1] = curfield.velZ(t)
    elseif (typeof(curfield.velX) == SinDef)
        curfield.u[1] = curfield.velX(t)
        curfield.w[1] = curfield.velZ(t)
    elseif (typeof(curfield.velX) == ConstDef)
        curfield.u[1] = curfield.velX(t)
        curfield.w[1] = curfield.velZ(t)
    end
end

# ---------------------------------------------------------------------------------------------
# Function updating the dimensional kinematic parameters
function update_kinem(surf::TwoDSurf, t)

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
    end
    # ---------------------------------------------------------------------------------------------
    return surf
end
# END kinem function
# ---------------------------------------------------------------------------------------------

function update_kinem(surf::TwoDSurf_2DOF, dt)
    #Kinematics are updated according to the 2DOF response
    surf.kinem.alphadot = surf.kinem.alphadot_pr + (dt/12.)*(23*surf.kinem.alphaddot_pr - 16*surf.kinem.alphaddot_pr2 + 5*surf.kinem.alphaddot_pr3)
    surf.kinem.hdot = surf.kinem.hdot_pr + (dt/12)*(23*surf.kinem.hddot_pr-16*surf.kinem.hddot_pr2 + 5*surf.kinem.hddot_pr3)
    surf.kinem.alpha = surf.kinem.alpha_pr + (dt/12)*(23*surf.kinem.alphadot_pr-16*surf.kinem.alphadot_pr2 + 5*surf.kinem.alphadot_pr3)
    surf.kinem.h = surf.kinem.h_pr + (dt/12)*(23*surf.kinem.hdot_pr - 16*surf.kinem.hdot_pr2 + 5*surf.kinem.hdot_pr3)
end


function update_kinem(surf::TwoDFreeSurf, dt)
    #Kinematics are updated according to the 2DOF response
    surf.kinem.alphadot = surf.kinem.alphadot_pr + (dt/12.)*(23*surf.kinem.alphaddot_pr - 16*surf.kinem.alphaddot_pr2 + 5*surf.kinem.alphaddot_pr3)
    surf.kinem.hdot = surf.kinem.hdot_pr + (dt/12)*(23*surf.kinem.hddot_pr-16*surf.kinem.hddot_pr2 + 5*surf.kinem.hddot_pr3)
    surf.kinem.alpha = surf.kinem.alpha_pr + (dt/12)*(23*surf.kinem.alphadot_pr-16*surf.kinem.alphadot_pr2 + 5*surf.kinem.alphadot_pr3)
    surf.kinem.h = surf.kinem.h_pr + (dt/12)*(23*surf.kinem.hdot_pr - 16*surf.kinem.hdot_pr2 + 5*surf.kinem.hdot_pr3)
    surf.kinem.u = surf.kinem.u_pr + (dt/12)*(23*surf.kinem.udot_pr - 16*surf.kinem.udot_pr2 + 5*surf.kinem.udot_pr3)
end



# ---------------------------------------------------------------------------------------------
# Function updating the kinematic parameters
function update_kinem(surf::TwoDSurfwFlap, t)
# Added by Laura Merchant 2016
    # ---------------------------------------------------------------------------------------------
    # Pitch kinematics
    if (typeof(surf.kindef.alpha) == EldUpDef)
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
    end

    # ---------------------------------------------------------------------------------------------
    # Plunge kinematics
    if (typeof(surf.kindef.h) == EldUpDef)
        surf.kinem.h = surf.kindef.h(t)*surf.c
        surf.kinem.hdot = ForwardDiff.derivative(surf.kindef.h,t)*surf.uref
    elseif (typeof(surf.kindef.h) == EldUpIntDef)
        surf.kinem.h = surf.kindef.h(t)*surf.c
        surf.kinem.hdot = EldUpDef(30,0.2,0.8)(t)*surf.kindef.h.amp*surf.uref
        #surf.kinem.hdot = ForwardDiff.derivative(surf.kindef.h,t)*surf.uref
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
      surf.kinem.hdot = ForwardDiff.derivative(surf.kindef.h,t)*surf.uref
    end

    # ---------------------------------------------------------------------------------------------
    # Velocity
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
    end

    # ---------------------------------------------------------------------------------------------
    # Deformation kinematics
    if (typeof(surf.kindef.n) == CosDef)
        surf.kinem.n = surf.kindef.n(t)
        surf.kinem.ndot = ForwardDiff.derivative(surf.kindef.n,t)*surf.uref/surf.c
    elseif (typeof(surf.kindef.n) == EldRampReturnDef)
        surf.kinem.n = surf.kindef.n(t)
        surf.kinem.ndot = ForwardDiff.derivative(surf.kindef.n,t)*surf.uref/surf.c
    elseif (typeof(surf.kindef.n) == ConstDef)
        surf.kinem.n = surf.kindef.n(t)
        surf.kinem.ndot = 0.
    elseif (typeof(surf.kindef.n) == SinDef)
        surf.kinem.n = surf.kindef.n(t)
        surf.kinem.ndot = ForwardDiff.derivative(surf.kindef.n,t)*surf.uref/surf.c
    end

    return surf
end
# END kinemwFlap function
# ---------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------
# Updates the bound vorticity distribution: eqn (2.1) in Ramesh et al. (2013)
# determines the strength of bound vortices
# determines the x, z components of the bound vortices
function update_bv(surf::TwoDSurf)
    gamma = zeros(surf.ndiv)
    for ib = 1:surf.ndiv
        gamma[ib] = (surf.a0[1]*(1 + cos(surf.theta[ib])))
        for ia = 1:surf.naterm
            gamma[ib] = gamma[ib] + surf.aterm[ia]*sin(ia*surf.theta[ib])*sin(surf.theta[ib])
        end
        gamma[ib] = gamma[ib]*surf.uref*surf.c
    end



    for ib = 2:surf.ndiv
        surf.bv[ib-1].s = (gamma[ib]+gamma[ib-1])*(surf.theta[2]-surf.theta[1])/2.
        surf.bv[ib-1].x = (surf.bnd_x[ib] + surf.bnd_x[ib-1])/2.
        surf.bv[ib-1].z = (surf.bnd_z[ib] + surf.bnd_z[ib-1])/2.
    end
end

function update_bv(surf::TwoDSurf_2DOF)
    gamma = zeros(surf.ndiv)
    for ib = 1:surf.ndiv
        gamma[ib] = (surf.a0[1]*(1 + cos(surf.theta[ib])))
        for ia = 1:surf.naterm
            gamma[ib] = gamma[ib] + surf.aterm[ia]*sin(ia*surf.theta[ib])*sin(surf.theta[ib])
        end
        gamma[ib] = gamma[ib]*surf.uref*surf.c
    end

    for ib = 2:surf.ndiv
        surf.bv[ib-1].s = (gamma[ib]+gamma[ib-1])*(surf.theta[2]-surf.theta[1])/2.
        surf.bv[ib-1].x = (surf.bnd_x[ib] + surf.bnd_x[ib-1])/2.
        surf.bv[ib-1].z = (surf.bnd_z[ib] + surf.bnd_z[ib-1])/2.
    end
end

function update_bv(surf::TwoDFreeSurf)
    gamma = zeros(surf.ndiv)
    for ib = 1:surf.ndiv
        gamma[ib] = (surf.a0[1]*(1 + cos(surf.theta[ib])))
        for ia = 1:surf.naterm
            gamma[ib] = gamma[ib] + surf.aterm[ia]*sin(ia*surf.theta[ib])*sin(surf.theta[ib])
        end
        gamma[ib] = gamma[ib]*surf.uref*surf.c
    end

    for ib = 2:surf.ndiv
        surf.bv[ib-1].s = (gamma[ib]+gamma[ib-1])*(surf.theta[2]-surf.theta[1])/2.
        surf.bv[ib-1].x = (surf.bnd_x[ib] + surf.bnd_x[ib-1])/2.
        surf.bv[ib-1].z = (surf.bnd_z[ib] + surf.bnd_z[ib-1])/2.
    end
end

# ---------------------------------------------------------------------------------------------
# Update KinemPar2DOF for 2DOF simulations
function update_kinem2DOF(surf::TwoDSurf_2DOF)
    surf.kinem.alpha_pr = surf.kinem.alpha
    surf.kinem.alpha_pr2 = surf.kinem.alpha_pr
    surf.kinem.alpha_pr3 = surf.kinem.alpha_pr2
    surf.kinem.h_pr = surf.kinem.h
    surf.kinem.h_pr2 = surf.kinem.h_pr
    surf.kinem.h_pr3 = surf.kinem.h_pr2
    surf.kinem.alphadot_pr = surf.kinem.alphadot
    surf.kinem.alphadot_pr2 = surf.kinem.alphadot_pr
    surf.kinem.alphadot_pr3 = surf.kinem.alphadot_pr2
    surf.kinem.hdot_pr = surf.kinem.hdot
    surf.kinem.hdot_pr2 = surf.kinem.hdot_pr
    surf.kinem.hdot_pr3 = surf.kinem.hdot_pr2
    surf.kinem.alphaddot_pr = surf.kinem.alphaddot
    surf.kinem.alphaddot_pr2 = surf.kinem.alphaddot_pr
    surf.kinem.alphaddot_pr3 = surf.kinem.alphaddot_pr2
    surf.kinem.hddot_pr = surf.kinem.hddot
    surf.kinem.hddot_pr2 = surf.kinem.hddot_pr
    surf.kinem.hddot_pr3 = surf.kinem.hddot_pr2
    return surf
end
# ---------------------------------------------------------------------------------------------
# Update KinemPar2DOF for Free simulations
function update_kinem2DOF(surf::TwoDFreeSurf)
    surf.kinem.alpha_pr = surf.kinem.alpha
    surf.kinem.h_pr = surf.kinem.h
    surf.kinem.u_pr = surf.kinem.u
    surf.kinem.alphadot_pr = surf.kinem.alphadot
    surf.kinem.alphadot_pr2 = surf.kinem.alphadot_pr
    surf.kinem.alphadot_pr3 = surf.kinem.alphadot_pr2
    surf.kinem.udot_pr = surf.kinem.udot
    surf.kinem.udot_pr2 = surf.kinem.udot_pr
    surf.kinem.udot_pr3 = surf.kinem.udot_pr2
    surf.kinem.hdot_pr = surf.kinem.hdot
    surf.kinem.hdot_pr2 = surf.kinem.hdot_pr
    surf.kinem.hdot_pr3 = surf.kinem.hdot_pr2
    surf.kinem.alphaddot_pr = surf.kinem.alphaddot
    surf.kinem.alphaddot_pr2 = surf.kinem.alphaddot_pr
    surf.kinem.alphaddot_pr3 = surf.kinem.alphaddot_pr2
    surf.kinem.hddot_pr = surf.kinem.hddot
    surf.kinem.hddot_pr2 = surf.kinem.hddot_pr
    surf.kinem.hddot_pr3 = surf.kinem.hddot_pr2
    return surf
end
# ---------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------
# Calculate hddot and alphaddot from forces based on 2DOF structural system
function calc_struct2DOF(surf::TwoDSurf_2DOF, cl::Float64, cm::Float64)
    m11 = 2./surf.c
    m12 = -surf.strpar.x_alpha*cos(surf.kinem.alpha)
    m21 = -2.*surf.strpar.x_alpha*cos(surf.kinem.alpha)/surf.c
    m22 = surf.strpar.r_alpha*surf.strpar.r_alpha

    R1 = 4*surf.strpar.kappa*surf.uref*surf.uref*cl/(pi*surf.c*surf.c) - 2*surf.strpar.w_h*surf.strpar.w_h*(surf.strpar.cubic_h_1*surf.kinem.h + surf.strpar.cubic_h_3*surf.kinem.h^3)/surf.c - surf.strpar.x_alpha*sin(surf.kinem.alpha)*surf.kinem.alphadot*surf.kinem.alphadot

    R2 = 8*surf.strpar.kappa*surf.uref*surf.uref*cm/(pi*surf.c*surf.c) - surf.strpar.w_alpha*surf.strpar.w_alpha*surf.strpar.r_alpha*surf.strpar.r_alpha*(surf.strpar.cubic_alpha_1*surf.kinem.alpha + surf.strpar.cubic_alpha_3*surf.kinem.alpha^3)

    surf.kinem.hddot = (1/(m11*m22-m21*m12))*(m22*R1-m12*R2)
    surf.kinem.alphaddot = (1/(m11*m22-m21*m12))*(-m21*R1+m11*R2)
    return surf
end

function calc_moveFree(surf::TwoDFreeSurf, cl::Float64, cd :: Float64,
cm :: Float64, cf :: Float64)
    accl_g = 9.8

    surf.kinem.udot = -2*surf.strpar.kappa*surf.uref*surf.uref*(cd+cf*sin(surf.kinem.alpha))/(pi*surf.c)
    surf.kinem.hddot = 2*surf.strpar.kappa*surf.uref*surf.uref*(cl+cf*cos(surf.kinem.alpha))/(pi*surf.c) - accl_g
    surf.kinem.alphaddot = 8*surf.strpar.kappa*surf.uref*surf.uref*cm/(pi*surf.c*surf.c*surf.strpar.r_g*surf.strpar.r_g)
end


function update_bv(surf::TwoDSurfwFlap)
    gamma = zeros(surf.ndiv)
    for ib = 1:surf.ndiv
        gamma[ib] = (surf.a0[1]*(1 + cos(surf.theta[ib])))
        for ia = 1:surf.naterm
            gamma[ib] = gamma[ib] + surf.aterm[ia]*sin(ia*surf.theta[ib])*sin(surf.theta[ib])
        end
        gamma[ib] = gamma[ib]*surf.uref*surf.c
    end

    for ib = 1:surf.ndiv-1
       surf.bv_prev[ib].s = surf.bv[ib].s
    end

    for ib = 2:surf.ndiv
        surf.bv[ib-1].s = (gamma[ib]+gamma[ib-1])*(surf.theta[2]-surf.theta[1])/2.
        surf.bv[ib-1].x = (surf.bnd_x[ib] + surf.bnd_x[ib-1])/2.
        surf.bv[ib-1].z = (surf.bnd_z[ib] + surf.bnd_z[ib-1])/2.
    end
    return surf
end
# ---------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------
# Numerical integration method: Trapezoidal
function trapz{T<:Real}(y::Vector{T}, x::Vector{T})
    local len = length(y)
    if (len != length(x))
        error("Vectors must be of same length")
    end
    r = 0.0
    for i in 2:len
        r += (x[i] - x[i-1]) * (y[i] + y[i-1])
    end
    r/2.0
end
# ---------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------
# Aerofoil camber calculation from coordinate file
function camber_calc(x::Vector,airfoil::ASCIIString)
    #Determine camber and camber slope on airfoil from airfoil input file

    ndiv = length(x);
    c = x[ndiv];

    cam = Array(Float64,ndiv)
    cam_slope = Array(Float64,ndiv)

    in_air = readdlm(airfoil);
    xcoord = in_air[:,1];
    ycoord = in_air[:,2];
    ncoord = length(xcoord);
    xcoord_sum = Array(Float64,ncoord);
    xcoord_sum[1] = 0;
    for i = 1:ncoord-1
        xcoord_sum[i+1] = xcoord_sum[i] + abs(xcoord[i+1]-xcoord[i]);
    end
    y_spl = Spline1D(xcoord_sum,ycoord);
    y_ans = Array(Float64,2*ndiv);

    for i=1:ndiv
        y_ans[i] = evaluate(y_spl,x[i]/c);
    end
    for i=ndiv+1:2*ndiv
        y_ans[i] = evaluate(y_spl,(x[ndiv]/c) + (x[i-ndiv]/c));
    end
    cam[1:ndiv] = [(y_ans[i] + y_ans[(2*ndiv) + 1 - i])*c/2 for i = ndiv:-1:1];
    cam[1] = 0;
    cam_spl = Spline1D(x,cam);
    cam_slope[1:ndiv] = Dierckx.derivative(cam_spl,x);
    return cam, cam_slope
end
# ---------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------
# Updating the deformation of the aerofoil: camber, spatial and time derivatives (from flap)
# Added by Laura Merchant 2016
function update_deform(surf::TwoDSurfwFlap, t)
    ndiv = length(surf.x);

    for i=1:ndiv
        if surf.x[i] < surf.x_b[1]
            surf.cam[i] = surf.cam_af[i];
        else
            surf.cam[i] = surf.cam_af[i] + (surf.x_b[1]-surf.x[i])*tan(surf.kinem.n);
            surf.cam_tder[i] = (surf.x_b[1]-surf.x[i])*surf.kinem.ndot*sec(surf.kinem.n)*sec(surf.kinem.n);
        end
    end

    cam_spl = Spline1D(surf.x,surf.cam);
    surf.cam_slope[1:surf.ndiv] = Dierckx.derivative(cam_spl,surf.x);
    return surf
end
# ---------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------
# Function for calculating the wake rollup
function wakeroll(surf::TwoDSurf, curfield::TwoDFlowField, dt)
    #Clean induced velocities
    for i = 1:length(curfield.tev)
        curfield.tev[i].vx = 0
        curfield.tev[i].vz = 0
    end

    for i = 1:length(curfield.lev)
        curfield.lev[i].vx = 0
        curfield.lev[i].vz = 0
    end

    #Velocities induced by free vortices on each other
    mutual_ind([curfield.tev; curfield.lev])

    #Add the influence of velocities induced by bound vortices
    utemp = zeros(length(curfield.tev)+length(curfield.lev))
    wtemp = zeros(length(curfield.tev)+length(curfield.lev))
    utemp, wtemp = ind_vel(surf.bv, [map(q -> q.x, curfield.tev); map(q -> q.x, curfield.lev)], [map(q -> q.z, curfield.tev); map(q -> q.z, curfield.lev)])

    for i = 1:length(curfield.tev)
        curfield.tev[i].vx += utemp[i]
        curfield.tev[i].vz += wtemp[i]
    end
    for i = length(curfield.tev)+1:length(utemp)
        curfield.lev[i-length(curfield.tev)].vx += utemp[i]
        curfield.lev[i-length(curfield.tev)].vz += wtemp[i]
    end

    #Convect free vortices with their induced velocities
    for i = 1:length(curfield.tev)
        curfield.tev[i].x += dt*curfield.tev[i].vx
        curfield.tev[i].z += dt*curfield.tev[i].vz
    end
    for i = 1:length(curfield.lev)
            curfield.lev[i].x += dt*curfield.lev[i].vx
        curfield.lev[i].z += dt*curfield.lev[i].vz
    end
    return curfield
end

function wakeroll(surf::TwoDSurf_2DOF, curfield::TwoDFlowField, dt)
    #Clean induced velocities
    for i = 1:length(curfield.tev)
        curfield.tev[i].vx = 0
        curfield.tev[i].vz = 0
    end

    for i = 1:length(curfield.lev)
        curfield.lev[i].vx = 0
        curfield.lev[i].vz = 0
    end

    #Velocities induced by free vortices on each other
    mutual_ind([curfield.tev; curfield.lev])

    #Add the influence of velocities induced by bound vortices
    utemp = zeros(length(curfield.tev)+length(curfield.lev))
    wtemp = zeros(length(curfield.tev)+length(curfield.lev))
    utemp, wtemp = ind_vel(surf.bv, [map(q -> q.x, curfield.tev); map(q -> q.x, curfield.lev)], [map(q -> q.z, curfield.tev); map(q -> q.z, curfield.lev)])

    for i = 1:length(curfield.tev)
        curfield.tev[i].vx += utemp[i]
        curfield.tev[i].vz += wtemp[i]
    end
    for i = length(curfield.tev)+1:length(utemp)
        curfield.lev[i-length(curfield.tev)].vx += utemp[i]
        curfield.lev[i-length(curfield.tev)].vz += wtemp[i]
    end

    #Convect free vortices with their induced velocities
    for i = 1:length(curfield.tev)
        curfield.tev[i].x += dt*curfield.tev[i].vx
        curfield.tev[i].z += dt*curfield.tev[i].vz
    end
    for i = 1:length(curfield.lev)
            curfield.lev[i].x += dt*curfield.lev[i].vx
        curfield.lev[i].z += dt*curfield.lev[i].vz
    end
    return curfield
end

function wakeroll(surf::TwoDFreeSurf, curfield::TwoDFlowField, dt)
    #Clean induced velocities
    for i = 1:length(curfield.tev)
        curfield.tev[i].vx = 0
        curfield.tev[i].vz = 0
    end

    for i = 1:length(curfield.lev)
        curfield.lev[i].vx = 0
        curfield.lev[i].vz = 0
    end

    #Assume external vortex is not affected by flow vortices
    #Velocities induced by free vortices on each other
    mutual_ind([curfield.tev; curfield.lev])

    #Add the influence of velocities induced by bound vortices
    utemp = zeros(length(curfield.tev)+length(curfield.lev))
    wtemp = zeros(length(curfield.tev)+length(curfield.lev))
    utemp, wtemp = ind_vel(surf.bv, [map(q -> q.x, curfield.tev); map(q -> q.x, curfield.lev)], [map(q -> q.z, curfield.tev); map(q -> q.z, curfield.lev)])

    for i = 1:length(curfield.tev)
        curfield.tev[i].vx += utemp[i]
        curfield.tev[i].vz += wtemp[i]
    end
    for i = length(curfield.tev)+1:length(utemp)
        curfield.lev[i-length(curfield.tev)].vx += utemp[i]
        curfield.lev[i-length(curfield.tev)].vz += wtemp[i]
    end

    #Convect free vortices with their induced velocities
    for i = 1:length(curfield.tev)
        curfield.tev[i].x += dt*curfield.tev[i].vx
        curfield.tev[i].z += dt*curfield.tev[i].vz
    end
    for i = 1:length(curfield.lev)
        curfield.lev[i].x += dt*curfield.lev[i].vx
        curfield.lev[i].z += dt*curfield.lev[i].vz
    end
    for i = 1:length(curfield.extv)
        curfield.extv[i].x += dt*curfield.extv[i].vx
        curfield.extv[i].z += dt*curfield.extv[i].vz
    end

    return curfield
end

function wakeroll(surf::TwoDSurfwFlap, curfield::TwoDFlowField, dt)
    #Clean induced velocities
    for i = 1:length(curfield.tev)
        curfield.tev[i].vx = 0
        curfield.tev[i].vz = 0
    end

    for i = 1:length(curfield.lev)
        curfield.lev[i].vx = 0
        curfield.lev[i].vz = 0
    end

    #Velocities induced by free vortices on each other
    mutual_ind([curfield.tev; curfield.lev])

    #Add the influence of velocities induced by bound vortices
    utemp = zeros(length(curfield.tev)+length(curfield.lev))
    wtemp = zeros(length(curfield.tev)+length(curfield.lev))
    utemp, wtemp = ind_vel(surf.bv, [map(q -> q.x, curfield.tev); map(q -> q.x, curfield.lev)], [map(q -> q.z, curfield.tev); map(q -> q.z, curfield.lev)])

    for i = 1:length(curfield.tev)
        curfield.tev[i].vx += utemp[i]
        curfield.tev[i].vz += wtemp[i]
    end
    for i = length(curfield.tev)+1:length(utemp)
        curfield.lev[i-length(curfield.tev)].vx += utemp[i]
        curfield.lev[i-length(curfield.tev)].vz += wtemp[i]
    end

    #Convect free vortices with their induced velocities
    for i = 1:length(curfield.tev)
        curfield.tev[i].x += dt*curfield.tev[i].vx
        curfield.tev[i].z += dt*curfield.tev[i].vz
    end
    for i = 1:length(curfield.lev)
            curfield.lev[i].x += dt*curfield.lev[i].vx
        curfield.lev[i].z += dt*curfield.lev[i].vz
    end
    return curfield
end
# ---------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------
# Places a trailing edge vortex
function place_tev(surf::TwoDSurf,field::TwoDFlowField,dt)
    ntev = length(field.tev)
    if ntev == 0
        xloc = surf.bnd_x[surf.ndiv] + 0.5*surf.kinem.u*dt
        zloc = surf.bnd_z[surf.ndiv]
        else
        xloc = surf.bnd_x[surf.ndiv]+(1./3.)*(field.tev[ntev].x - surf.bnd_x[surf.ndiv])

        zloc = surf.bnd_z[surf.ndiv]+(1./3.)*(field.tev[ntev].z - surf.bnd_z[surf.ndiv])
    end
    push!(field.tev,TwoDVort(xloc,zloc,0.,0.02*surf.c,0.,0.))
    return field
end

function place_tev(surf::TwoDSurf_2DOF,field::TwoDFlowField,dt)
    ntev = length(field.tev)
    if ntev == 0
        xloc = surf.bnd_x[surf.ndiv] + 0.5*surf.kinem.u*dt
        zloc = surf.bnd_z[surf.ndiv]
        else
        xloc = surf.bnd_x[surf.ndiv]+(1./3.)*(field.tev[ntev].x - surf.bnd_x[surf.ndiv])

        zloc = surf.bnd_z[surf.ndiv]+(1./3.)*(field.tev[ntev].z - surf.bnd_z[surf.ndiv])
    end
    push!(field.tev,TwoDVort(xloc,zloc,0.,0.02*surf.c,0.,0.))
    return field
end

function place_tev(surf::TwoDFreeSurf,field::TwoDFlowField,dt)
    ntev = length(field.tev)
    if ntev == 0
        xloc = surf.bnd_x[surf.ndiv] + 0.5*surf.kinem.u*dt
        zloc = surf.bnd_z[surf.ndiv]
        else
        xloc = surf.bnd_x[surf.ndiv]+(1./3.)*(field.tev[ntev].x - surf.bnd_x[surf.ndiv])

        zloc = surf.bnd_z[surf.ndiv]+(1./3.)*(field.tev[ntev].z - surf.bnd_z[surf.ndiv])
    end
    push!(field.tev,TwoDVort(xloc,zloc,0.,0.02*surf.c,0.,0.))
    return field
end

function place_tev(surf::TwoDSurfwFlap,field::TwoDFlowField,dt)
    ntev = length(field.tev)
    if ntev == 0
        xloc = surf.bnd_x[surf.ndiv] + 0.5*surf.kinem.u*dt
        zloc = surf.bnd_z[surf.ndiv]
        else
        xloc = surf.bnd_x[surf.ndiv]+(1./3.)*(field.tev[ntev].x - surf.bnd_x[surf.ndiv])

        zloc = surf.bnd_z[surf.ndiv]+(1./3.)*(field.tev[ntev].z - surf.bnd_z[surf.ndiv])
    end
    push!(field.tev,TwoDVort(xloc,zloc,0.,0.02*surf.c,0.,0.))
    return field
end
# ---------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------
# Places a leading edge vortex
function place_lev(surf::TwoDSurf,field::TwoDFlowField,dt)
    nlev = length(field.lev)

    le_vel_x = surf.kinem.u - surf.kinem.alphadot*sin(surf.kinem.alpha)*surf.pvt*surf.c + surf.uind[1]
    le_vel_z = -surf.kinem.alphadot*cos(surf.kinem.alpha)*surf.pvt*surf.c- surf.kinem.hdot + surf.wind[1]

    if (surf.levflag[1] == 0) then
        xloc = surf.bnd_x[1] + 0.5*le_vel_x*dt
        zloc = surf.bnd_z[1] + 0.5*le_vel_z*dt
    else
        xloc = surf.bnd_x[1]+(1./3.)*(field.lev[nlev].x - surf.bnd_x[1])
        zloc = surf.bnd_z[1]+(1./3.)*(field.lev[nlev].z - surf.bnd_z[1])
    end

    push!(field.lev,TwoDVort(xloc,zloc,0.,0.02*surf.c,0.,0.))

    return field
end

function place_lev(surf::TwoDSurf_2DOF,field::TwoDFlowField,dt)
    nlev = length(field.lev)

    le_vel_x = surf.kinem.u - surf.kinem.alphadot*sin(surf.kinem.alpha)*surf.pvt*surf.c + surf.uind[1]
    le_vel_z = -surf.kinem.alphadot*cos(surf.kinem.alpha)*surf.pvt*surf.c- surf.kinem.hdot + surf.wind[1]

    if (surf.levflag[1] == 0) then
        xloc = surf.bnd_x[1] + 0.5*le_vel_x*dt
        zloc = surf.bnd_z[1] + 0.5*le_vel_z*dt
    else
        xloc = surf.bnd_x[1]+(1./3.)*(field.lev[nlev].x - surf.bnd_x[1])
        zloc = surf.bnd_z[1]+(1./3.)*(field.lev[nlev].z - surf.bnd_z[1])
    end

    push!(field.lev,TwoDVort(xloc,zloc,0.,0.02*surf.c,0.,0.))

    return field
end

function place_lev(surf::TwoDFreeSurf,field::TwoDFlowField,dt)
    nlev = length(field.lev)

    le_vel_x = surf.kinem.u - surf.kinem.alphadot*sin(surf.kinem.alpha)*surf.pvt*surf.c + surf.uind[1]
    le_vel_z = -surf.kinem.alphadot*cos(surf.kinem.alpha)*surf.pvt*surf.c- surf.kinem.hdot + surf.wind[1]

    if (surf.levflag[1] == 0) then
        xloc = surf.bnd_x[1] + 0.5*le_vel_x*dt
        zloc = surf.bnd_z[1] + 0.5*le_vel_z*dt
    else
        xloc = surf.bnd_x[1]+(1./3.)*(field.lev[nlev].x - surf.bnd_x[1])
        zloc = surf.bnd_z[1]+(1./3.)*(field.lev[nlev].z - surf.bnd_z[1])
    end

    push!(field.lev,TwoDVort(xloc,zloc,0.,0.02*surf.c,0.,0.))

    return field
end

function place_lev(surf::TwoDSurfwFlap,field::TwoDFlowField,dt)
    nlev = length(field.lev)

    le_vel_x = surf.kinem.u - surf.kinem.alphadot*sin(surf.kinem.alpha)*surf.pvt*surf.c + surf.uind[1]
    le_vel_z = -surf.kinem.alphadot*cos(surf.kinem.alpha)*surf.pvt*surf.c- surf.kinem.hdot + surf.wind[1]

    if (surf.levflag[1] == 0) then
        xloc = surf.bnd_x[1] + 0.5*le_vel_x*dt
        zloc = surf.bnd_z[1] + 0.5*le_vel_z*dt
    else
        xloc = surf.bnd_x[1]+(1./3.)*(field.lev[nlev].x - surf.bnd_x[1])
        zloc = surf.bnd_z[1]+(1./3.)*(field.lev[nlev].z - surf.bnd_z[1])
    end

    push!(field.lev,TwoDVort(xloc,zloc,0.,0.02*surf.c,0.,0.))

    return field
end
# ---------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------
# Function determining the effects of interacting vorticies - velocities induced on each other - classical n-body problem
function mutual_ind(vorts::Vector{TwoDVort})
    for i = 1:length(vorts)
        for j = i+1:length(vorts)
            dx = vorts[i].x - vorts[j].x
            dz = vorts[i].z - vorts[j].z
            #source- tar
            dsq = dx*dx + dz*dz

            magitr = 1./(2*pi*sqrt(vorts[j].vc*vorts[j].vc*vorts[j].vc*vorts[j].vc + dsq*dsq))
            magjtr = 1./(2*pi*sqrt(vorts[i].vc*vorts[i].vc*vorts[i].vc*vorts[i].vc + dsq*dsq))

            vorts[j].vx -= dz * vorts[i].s * magjtr
            vorts[j].vz += dx * vorts[i].s * magjtr

            vorts[i].vx += dz * vorts[j].s * magitr
            vorts[i].vz -= dx * vorts[j].s * magitr
        end
    end
    return vorts
end


# ---------------------------------------------------------------------------------------------
# Function for updating the positions of the bound vortices
function update_boundpos(surf::TwoDSurf, dt::Float64)
    for i = 1:surf.ndiv
        surf.bnd_x[i] = surf.bnd_x[i] + dt*((surf.pvt*surf.c - surf.x[i])*sin(surf.kinem.alpha)*surf.kinem.alphadot - surf.kinem.u + surf.cam[i]*cos(surf.kinem.alpha)*surf.kinem.alphadot)
        surf.bnd_z[i] = surf.bnd_z[i] + dt*(surf.kinem.hdot + (surf.pvt*surf.c - surf.x[i])*cos(surf.kinem.alpha)*surf.kinem.alphadot - surf.cam[i]*sin(surf.kinem.alpha)*surf.kinem.alphadot)
    end
    return surf
end

function update_boundpos(surf::TwoDSurf_2DOF, dt::Float64)
    for i = 1:surf.ndiv
        surf.bnd_x[i] = surf.bnd_x[i] + dt*((surf.pvt*surf.c - surf.x[i])*sin(surf.kinem.alpha)*surf.kinem.alphadot - surf.kinem.u + surf.cam[i]*cos(surf.kinem.alpha)*surf.kinem.alphadot)
        surf.bnd_z[i] = surf.bnd_z[i] + dt*(surf.kinem.hdot + (surf.pvt*surf.c - surf.x[i])*cos(surf.kinem.alpha)*surf.kinem.alphadot - surf.cam[i]*sin(surf.kinem.alpha)*surf.kinem.alphadot)
    end
    return surf
end

function update_boundpos(surf::TwoDFreeSurf, dt::Float64)
    for i = 1:surf.ndiv
        surf.bnd_x[i] = surf.bnd_x[i] + dt*((surf.pvt*surf.c - surf.x[i])*sin(surf.kinem.alpha)*surf.kinem.alphadot - surf.kinem.u + surf.cam[i]*cos(surf.kinem.alpha)*surf.kinem.alphadot)
        surf.bnd_z[i] = surf.bnd_z[i] + dt*(surf.kinem.hdot + (surf.pvt*surf.c - surf.x[i])*cos(surf.kinem.alpha)*surf.kinem.alphadot - surf.cam[i]*sin(surf.kinem.alpha)*surf.kinem.alphadot)
    end
    return surf
end

function update_boundpos(surf::TwoDSurfwFlap, dt::Float64)
    for i = 1:surf.ndiv
        surf.bnd_x[i] = surf.bnd_x[i] + dt*((surf.pvt*surf.c - surf.x[i])*sin(surf.kinem.alpha)*surf.kinem.alphadot - surf.kinem.u + surf.cam[i]*cos(surf.kinem.alpha)*surf.kinem.alphadot)
        surf.bnd_z[i] = surf.bnd_z[i] + dt*(surf.kinem.hdot + (surf.pvt*surf.c - surf.x[i])*cos(surf.kinem.alpha)*surf.kinem.alphadot - surf.cam[i]*sin(surf.kinem.alpha)*surf.kinem.alphadot)
    end
    return surf
end
# ---------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------
# Function to calculate induced velocities by a set of votices at a target location
function ind_vel(src::Vector{TwoDVort},t_x,t_z)
    #'s' stands for source and 't' for target
    uind = zeros(length(t_x))
    wind = zeros(length(t_x))

    for itr = 1:length(t_x)
        for isr = 1:length(src)
            xdist = src[isr].x - t_x[itr]
            zdist = src[isr].z - t_z[itr]
            distsq = xdist*xdist + zdist*zdist
            uind[itr] = uind[itr] - src[isr].s*zdist/(2*pi*sqrt(src[isr].vc*src[isr].vc*src[isr].vc*src[isr].vc + distsq*distsq))
            wind[itr] = wind[itr] + src[isr].s*xdist/(2*pi*sqrt(src[isr].vc*src[isr].vc*src[isr].vc*src[isr].vc + distsq*distsq))
        end
    end
    return uind, wind
end
# ---------------------------------------------------------------------------------------------
