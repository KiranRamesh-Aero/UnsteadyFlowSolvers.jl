#Function for estimating a problem's time step

#Simple linear interpolation function
function interp(x1 ::Float64, x2 :: Float64, y1 :: Float64, y2 :: Float64, x::Float64)
    y = y1 + (y2 - y1)*(x - x1)/(x2 - x1)
    return y
end

    
function find_tstep(kin:: Array{CosDef})
    dtstar = 0.015
    for i = 1:length(kin)
        dt_tmp = 0.015*0.2/(kin[i].k*kin[i].amp)
        dtstar = minimum([dt_tmp dtstar])
    end
    return dtstar
end

function find_tstep(kin:: Array{SinDef})
    dtstar = 0.015
    for i = 1:length(kin)
        dt_tmp = 0.015*0.2/(kin[i].k*kin[i].amp)
        dtstar = minimum([dt_tmp dtstar])
    end
    return dtstar
end

function find_tstep(kin:: EldUpDef)
    dtstar = 0.015
    dtstar = minimum([0.015*0.2/kin.K 0.015])
    return dtstar
end                        

function find_tstep(kin:: EldUptstartDef)
    dtstar = 0.015
    dtstar = minimum([0.015*0.2/kin.K 0.015])
    return dtstar
end

function find_tstep(kin:: EldRampReturnDef)
    dtstar = 0.015
    dtstar = minimum([0.015*0.2/kin.K 0.015])
    return dtstar
end

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

function update_adot(surf::TwoDSurf,dt)
    surf.a0dot[1] = (surf.a0[1] - surf.a0prev[1])/dt
    for ia = 1:3
        surf.adot[ia] = (surf.aterm[ia]-surf.aprev[ia])/dt
    end
    return surf
end

function update_adot(surf::TwoDSurf_2DOF,dt)
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

function update_indbound(surf::TwoDSurfLV, curfield::TwoDFlowField)
    for i = 1:surf.npanel
        ui, wi = ind_vel([curfield.tev[1:length(curfield.tev)-1]; curfield.lev], surf.lv[i].xc_I, surf.lv[i].zc_I)
        surf.lv[i].uind_I = ui[1]
        surf.lv[i].wind_I = wi[1]
        (surf.lv[i].uind, surf.lv[i].wind) = surf.tlg*[surf.lv[i].uind_I; surf.lv[i].wind_I] 
    end
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


function update_kinem(surf::TwoDSurfLV, t)

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

function update_kinem(surf::TwoDSurf_2DOF, dt, cl, cm)
    #Update previous terms 
    surf.kinem.alpha_pr = surf.kinem.alpha
    surf.kinem.h_pr = surf.kinem.h
    
    surf.kinem.alphadot_pr3 = surf.kinem.alphadot_pr2
    surf.kinem.alphadot_pr2 = surf.kinem.alphadot_pr
    surf.kinem.alphadot_pr = surf.kinem.alphadot
    
    surf.kinem.hdot_pr3 = surf.kinem.hdot_pr2
    surf.kinem.hdot_pr2 = surf.kinem.hdot_pr
    surf.kinem.hdot_pr = surf.kinem.hdot
    
    surf.kinem.alphaddot_pr3 = surf.kinem.alphaddot_pr2
    surf.kinem.alphaddot_pr2 = surf.kinem.alphaddot_pr
        
    surf.kinem.hddot_pr3 = surf.kinem.hddot_pr2
    surf.kinem.hddot_pr2 = surf.kinem.hddot_pr

    # Calculate hddot and alphaddot from forces based on 2DOF structural system
    m11 = 2./surf.c
    m12 = -surf.strpar.x_alpha*cos(surf.kinem.alpha)
    m21 = -2.*surf.strpar.x_alpha*cos(surf.kinem.alpha)/surf.c
    m22 = surf.strpar.r_alpha*surf.strpar.r_alpha
    
    R1 = 4*surf.strpar.kappa*surf.uref*surf.uref*cl/(pi*surf.c*surf.c) - 2*surf.strpar.w_h*surf.strpar.w_h*(surf.strpar.cubic_h_1*surf.kinem.h + surf.strpar.cubic_h_3*surf.kinem.h^3)/surf.c - surf.strpar.x_alpha*sin(surf.kinem.alpha)*surf.kinem.alphadot*surf.kinem.alphadot
    
    R2 = 8*surf.strpar.kappa*surf.uref*surf.uref*cm/(pi*surf.c*surf.c) - surf.strpar.w_alpha*surf.strpar.w_alpha*surf.strpar.r_alpha*surf.strpar.r_alpha*(surf.strpar.cubic_alpha_1*surf.kinem.alpha + surf.strpar.cubic_alpha_3*surf.kinem.alpha^3)
    
    surf.kinem.hddot_pr = (1/(m11*m22-m21*m12))*(m22*R1-m12*R2)
    surf.kinem.alphaddot_pr = (1/(m11*m22-m21*m12))*(-m21*R1+m11*R2)

    
    #Kinematics are updated according to the 2DOF response
    surf.kinem.alphadot = surf.kinem.alphadot_pr + (dt/12.)*(23*surf.kinem.alphaddot_pr - 16*surf.kinem.alphaddot_pr2 + 5*surf.kinem.alphaddot_pr3)
    surf.kinem.hdot = surf.kinem.hdot_pr + (dt/12)*(23*surf.kinem.hddot_pr-16*surf.kinem.hddot_pr2 + 5*surf.kinem.hddot_pr3)
    surf.kinem.alpha = surf.kinem.alpha_pr + (dt/12)*(23*surf.kinem.alphadot_pr-16*surf.kinem.alphadot_pr2 + 5*surf.kinem.alphadot_pr3)
    surf.kinem.h = surf.kinem.h_pr + (dt/12)*(23*surf.kinem.hdot_pr - 16*surf.kinem.hdot_pr2 + 5*surf.kinem.hdot_pr3)
    return surf
end


function update_kinem(surf::TwoDFreeSurf, dt)
    #Kinematics are updated according to the 2DOF response

    

    surf.kinem.alphadot = surf.kinem.alphadot_pr + (dt/12.)*(23*surf.kinem.alphaddot_pr - 16*surf.kinem.alphaddot_pr2 + 5*surf.kinem.alphaddot_pr3)
    surf.kinem.hdot = surf.kinem.hdot_pr + (dt/12)*(23*surf.kinem.hddot_pr-16*surf.kinem.hddot_pr2 + 5*surf.kinem.hddot_pr3)
    surf.kinem.alpha = surf.kinem.alpha_pr + (dt/12)*(23*surf.kinem.alphadot_pr-16*surf.kinem.alphadot_pr2 + 5*surf.kinem.alphadot_pr3)
    surf.kinem.h = surf.kinem.h_pr + (dt/12)*(23*surf.kinem.hdot_pr - 16*surf.kinem.hdot_pr2 + 5*surf.kinem.hdot_pr3)
    surf.kinem.u = surf.kinem.u_pr + (dt/12)*(23*surf.kinem.udot_pr - 16*surf.kinem.udot_pr2 + 5*surf.kinem.udot_pr3)
    return surf
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
# function update_kinem2DOF(surf::TwoDSurf_2DOF)
#     surf.kinem.alpha_pr = surf.kinem.alpha
#     surf.kinem.h_pr = surf.kinem.h

#     surf.kinem.alphadot_pr3 = surf.kinem.alphadot_pr2
#     surf.kinem.alphadot_pr2 = surf.kinem.alphadot_pr
#     surf.kinem.alphadot_pr = surf.kinem.alphadot
   
#     surf.kinem.hdot_pr3 = surf.kinem.hdot_pr2
#     surf.kinem.hdot_pr2 = surf.kinem.hdot_pr
#     surf.kinem.hdot_pr = surf.kinem.hdot

#     surf.kinem.alphaddot_pr3 = surf.kinem.alphaddot_pr2
#     surf.kinem.alphaddot_pr2 = surf.kinem.alphaddot_pr
#     surf.kinem.alphaddot_pr = surf.kinem.alphaddot

#     surf.kinem.hddot_pr3 = surf.kinem.hddot_pr2
#     surf.kinem.hddot_pr2 = surf.kinem.hddot_pr
#     surf.kinem.hddot_pr = surf.kinem.hddot
 

#     return surf
# end
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
# function calc_struct2DOF(surf::TwoDSurf_2DOF, cl::Float64, cm::Float64)
#     m11 = 2./surf.c
#     m12 = -surf.strpar.x_alpha*cos(surf.kinem.alpha)
#     m21 = -2.*surf.strpar.x_alpha*cos(surf.kinem.alpha)/surf.c
#     m22 = surf.strpar.r_alpha*surf.strpar.r_alpha

#     R1 = 4*surf.strpar.kappa*surf.uref*surf.uref*cl/(pi*surf.c*surf.c) - 2*surf.strpar.w_h*surf.strpar.w_h*(surf.strpar.cubic_h_1*surf.kinem.h + surf.strpar.cubic_h_3*surf.kinem.h^3)/surf.c - surf.strpar.x_alpha*sin(surf.kinem.alpha)*surf.kinem.alphadot*surf.kinem.alphadot

#     R2 = 8*surf.strpar.kappa*surf.uref*surf.uref*cm/(pi*surf.c*surf.c) - surf.strpar.w_alpha*surf.strpar.w_alpha*surf.strpar.r_alpha*surf.strpar.r_alpha*(surf.strpar.cubic_alpha_1*surf.kinem.alpha + surf.strpar.cubic_alpha_3*surf.kinem.alpha^3)

#     surf.kinem.hddot_pr = (1/(m11*m22-m21*m12))*(m22*R1-m12*R2)
#     surf.kinem.alphaddot_pr = (1/(m11*m22-m21*m12))*(-m21*R1+m11*R2)
#     return surf
# end

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
function camber_calc(x::Vector,airfoil::String)
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

function wakeroll(surf::TwoDSurfLV, curfield::TwoDFlowField, dt)
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

    bv = TwoDVort[]
    for i = 1:surf.npanel
        push!(bv, TwoDVort(surf.lv[i].xv,surf.lv[i].zv,surf.lv[i].s,0.02*surf.c,0,0))
    end
    
    utemp, wtemp = ind_vel(bv, [map(q -> q.x, curfield.tev); map(q -> q.x, curfield.lev)], [map(q -> q.z, curfield.tev); map(q -> q.z, curfield.lev)])

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

function place_tev(surf::TwoDSurfLV,field::TwoDFlowField,dt)
    ntev = length(field.tev)
    if ntev == 0
        xloc = surf.x[surf.npanel+1] + surf.tevloc*surf.kinem.u*dt
        zloc = surf.cam[surf.npanel+1]
    else
        (xteI, zteI) = surf.tlg*[surf.x[surf.npanel+1]; surf.cam[surf.npanel+1]] + [surf.X0[1]; surf.X0[2]]
        xloc = xteI +(1./3.)*(field.tev[ntev].x - xteI)
        zloc = zteI + (1./3.)*(field.tev[ntev].z - zteI)
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

    if (surf.levflag[1] == 0)
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

    if (surf.levflag[1] == 0)
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

function fgfromrhoGauss(rho::Float64)
    f = sqrt(2./pi)*exp(-rho^2/2.)
    g = erf(rho/sqrt(2.)) - rho*f
    return f, g
end

function mutual_ind(vorts::Vector{ThreeDVort})
      for i = 1:length(vorts)
          for j = i+1:length(vorts)
              dx = vorts[i].x[1] - vorts[j].x[1]
              dy = vorts[i].x[2] - vorts[j].x[2]
              dz = vorts[i].x[3] - vorts[j].x[3]
              #source - tar

              r = sqrt(dx*dx + dz*dz + dy*dy)
              rhoi = r/vorts[i].vc
              rhoj = r/vorts[i].vc

              f, g = fgfromrhoGauss(rhoi)
            
              vorts[j].v += -g*cross([dx; dy; dz],vorts[i].s)/(4.*pi*r^3)
                
              vorts[j].ds += -(-g*cross(vorts[j].s, vorts[i].s)/rhoi^3 + 
                ((3*g/rhoi^3 - f)*dot(vorts[j].s, [dx;dy;dz])*cross([dx;dy;dz],vorts[i].s))/r^2)/(4*pi*vorts[i].vc^3)  
                
            
            f, g = fgfromrhoGauss(rhoj)
            vorts[i].v += -g*cross([dx; dy; dz],vorts[j].s)/(4.*pi*r^3)
                
            vorts[i].ds += -(-g*cross(vorts[i].s, vorts[j].s)/rhoj^3 + 
            ((3*g/rhoj^3 - f)*dot(vorts[i].s, [dx;dy;dz])*cross([dx;dy;dz],vorts[j].s))/r^2)/(4*pi*vorts[j].vc^3)  
                

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
#Routines for integral BL calculation (Matsushita et al.)

function FfromE(E::Float64)
      F = 4.8274*E^4 - 5.9816*E^3 + 4.0274*E^2 + 0.23247*E + 0.15174
end

function BfromE(E::Float64)
    if E < -0.0616
        B = -225.86*E^3 - 3016.6*E^2 - 208.68*E - 17.915
    elseif E > -0.0395
        B = 131.9*E^3 - 167.32*E^2 + 76.642*E - 11.068
    else
        B = 0.5*(-225.86*E^3 - 3016.6*E^2 - 208.68*E - 17.915 + 131.9*E^3 - 167.32*E^2 + 76.642*E - 11.068)
    end
end

function SfromE(E::Float64)
    if E < -0.0582
        S = 451.55*E^3 + 2010.*E^2 + 138.96*E + 11.296
    elseif E > -0.042
        S = -96.739*E^3 + 117.74*E^2 - 46.432*E + 6.8074
    else
        S = 0.5*(451.55*E^3 + 2010.*E^2 + 138.96*E + 11.296 - 96.739*E^3 + 117.74*E^2 - 46.432*E + 6.8074)
    end
end

function dfdefromE(E::Float64)
    dfde = 4*4.8274*E^3 - 3*5.9816*E^2 + 2*4.0274*E + 0.23247
end

function calcdt(dx::Float64, cfl::Float64, lamb::Array{Float64,2})
    dt = 10000
    for ic = 2:size(lamb,2)-1
        dti = cfl*(dx/(abs(lamb[1,ic]) + abs(lamb[2,ic])))
        if dti < dt
            dt=dti
        end
    end
    return dt
end

function calcEigen(ue::Vector{Float64}, E::Vector{Float64}, F::Vector{Float64}, dfde::Vector{Float64})
    ncell = length(E) - 2
    lamb = zeros(2,ncell+2)
    for ic = 1:ncell+2
        aq = 1.
        bq = -ue[ic]*(dfde[ic] - 1.)
        cq = ue[ic]*ue[ic]*(E[ic]*dfde[ic] - F[ic])
        lamb[1,ic] = (-bq + sqrt(bq*bq - 4*aq*cq))/(2*aq)
        lamb[2,ic] = (-bq - sqrt(bq*bq - 4*aq*cq))/(2*aq)

        #Always have lamb1 > lamb2
        if lamb[2,ic] > lamb[1,ic]
            temp  = lamb[2,ic]
            lamb[2,ic] = lamb[1,ic]
            lamb[1,ic] = temp
        end
    end
    return lamb
end


function calc_flux(lamb::Array{Float64,2}, ue::Vector{Float64}, E::Vector{Float64}, del::Vector{Float64}, F::Vector{Float64})
    ncell = length(ue) - 2
    flux = zeros(2,2,ncell+2)
    Apos = zeros(2,2)
    Aneg = zeros(2,2)

    for ic = 1:ncell+2
        if lamb[1,ic] >= 0. && lamb[2,ic] >= 0.
            flux[1,1,ic] = ue[ic]*E[ic]*del[ic]
            flux[1,2,ic] = ue[ic]*F[ic]*del[ic]
        elseif lamb[1,ic] < 0. && lamb[2,ic] < 0.
            flux[1,:,ic] = 0.
            else
            Apos[1,1] = (ue[ic]*lamb[1,ic]/(lamb[1,ic] - lamb[2,ic]))*(-1. - lamb[2,ic]/ue[ic])
            Apos[1,2] = ue[ic]*lamb[1,ic]/(lamb[1,ic] - lamb[2,ic])
            Apos[2,1] = -(ue[ic]*lamb[1,ic]/(lamb[1,ic] - lamb[2,ic]))*(1 + lamb[1,ic]/ue[ic])*(1 + lamb[2,ic]/ue[ic])
            Apos[2,2] = (ue[ic]*lamb[1,ic]/(lamb[1,ic] - lamb[2,ic]))*(1 + lamb[1,ic]/ue[ic])
            flux[1,1,ic] =  Apos[1,1]*del[ic] + Apos[1,2]*(E[ic] + 1.)*del[ic]
            flux[1,2,ic] = Apos[2,1]*del[ic] + Apos[2,2]*(E[ic] + 1.)*del[ic]
        end
    end
    for ic = 1:ncell+2
        if lamb[1,ic] >= 0. && lamb[2,ic] >= 0.
            flux[2,:,ic] = 0.
         elseif lamb[1,ic] < 0. && lamb[2,ic] < 0.
            flux[2,1,ic] = ue[ic]*E[ic]*del[ic]
            flux[2,2,ic] = ue[ic]*F[ic]*del[ic]
         else
            Aneg[1,1] = (ue[ic]*lamb[2,ic]/(lamb[1,ic] - lamb[2,ic]))*(1. + lamb[1,ic]/ue[ic])
            Aneg[1,2] = -ue[ic]*lamb[2,ic]/(lamb[1,ic] - lamb[2,ic])
            Aneg[2,1] = (ue[ic]*lamb[2,ic]/(lamb[1,ic] - lamb[2,ic]))*(1. + lamb[1,ic]/ue[ic])*(1. + lamb[2,ic]/ue[ic])
            Aneg[2,2] = (ue[ic]*lamb[2,ic]/(lamb[1,ic] - lamb[2,ic]))*(-1. - lamb[2,ic]/ue[ic])
            flux[2,1,ic] = Aneg[1,1]*del[ic] + Aneg[1,2]*(E[ic] + 1.)*del[ic]
            flux[2,2,ic] = Aneg[2,1]*del[ic] + Aneg[2,2]*(E[ic] + 1.)*del[ic]
         end
    end
    return flux
end


function IBLm_cyl(ntime = 1000, cfl = 0.5, ncell = 64, ttime = 1.4)

    sepflag = Int(0)
    dx = pi/(real(ncell) + 1.)
    t = 0.005

    x = zeros(ncell+2)
    ue = zeros(ncell+2)
    uex = zeros(ncell+2)
    uet = zeros(ncell+2)

    #Set sources
    for i = 1:ncell+2
        x[i] = real(i-1)*dx
        ue[i] = 2*sin(x[i])
        uex[i] = 2*cos(x[i])
    end

    uet[1:ncell+2] = 0

    E = zeros(ncell+2)
    B = zeros(ncell+2)
    del = zeros(ncell+2)
    F = zeros(ncell+2)
    dfde = zeros(ncell+2)
    S = zeros(ncell+2)
    unk = zeros(2,ncell+2)

    lamb = zeros(2,ncell+2)
    unkh = zeros(2,ncell+2)
    rhs = zeros(2,ncell+2)
    flux = zeros(2,2,ncell+2)
    crit = zeros(ncell+2)
    Apos = zeros(2,2)
    Aneg = zeros(2,2)

    #Set initial conditions
    E[:] = 0.4142
    for ic = 1:ncell+2
        B[ic] = BfromE(E[ic])
        F[ic] = FfromE(E[ic])
    end
    del = sqrt(B*t)

    unk[1,:] = del
    unk[2,:] = del.*(E + 1.)

    #Main loop over time steps
    for i = 1:ntime
        #i = 1
        unk[:,1] = 2*unk[:,2] - unk[:,3]
        unk[:,ncell+2] = 2*unk[:,ncell+1] - unk[:,ncell]
        #Calculate derived quantities
        for ic = 1:ncell+2
            del[ic] = unk[1,ic]
            E[ic] = (unk[2,ic]./del[ic]) - 1.
            F[ic] = FfromE(E[ic])
            B[ic] = BfromE(E[ic])
            S[ic] = SfromE(E[ic])
            dfde[ic] = dfdefromE(E[ic])
        end

        #Compute eigenvalues
        lamb = calcEigen(ue, E, F, dfde)

        #Compute timestep
        dt = calcdt(dx, cfl, lamb)

        #Compute fluxes
        flux = calc_flux(lamb, ue, E, del, F)

        #compute rhs
        for ic = 2:ncell+1
        rhs[1,ic] = B[ic]/(2*del[ic]) - del[ic]*uet[ic]/ue[ic] - (E[ic] + 1.)*del[ic]*uex[ic]
            rhs[2,ic] = S[ic]/del[ic] - 2*E[ic]*del[ic]*uet[ic]/ue[ic] - 2*F[ic]*del[ic]*uex[ic]
        end

        for ic = 2:ncell+1
            unkh[1,ic] = unk[1,ic] - (dt/dx)*(flux[1,1,ic] - flux[1,1,ic-1] + flux[2,1,ic+1]
                                              - flux[2,1,ic]) + dt*rhs[1,ic]
            unkh[2,ic] = unk[2,ic] - (dt/dx)*(flux[1,2,ic] - flux[1,2,ic-1] + flux[2,2,ic+1]
                                              - flux[2,2,ic]) + dt*rhs[2,ic]
        end
        ic = ncell+1
        #unkh[1,ic] = unk[1,ic] - (dt/dx)*(flux[1,1,ic] - flux[1,1,ic-1]) + dt*rhs[1,ic]
        #unkh[2,ic] = unk[2,ic] - (dt/dx)*(flux[1,2,ic] - flux[1,2,ic-1]) + dt*rhs[2,ic]

        #Update ghost cells
        unkh[:,1] = 2*unkh[:,2] - unkh[:,3]
        unkh[:,ncell+2] = 2*unkh[:,ncell+1] - unkh[:,ncell]

        #Calculate derived quantities
        for ic = 1:ncell+2
            del[ic] = unkh[1,ic]
            E[ic] = (unkh[2,ic]./del[ic]) - 1.
            F[ic] = FfromE(E[ic])
            B[ic] = BfromE(E[ic])
            S[ic] = SfromE(E[ic])
            dfde[ic] = dfdefromE(E[ic])
        end

        #Compute eigenvalues
        lamb = calcEigen(ue, E, F, dfde)

        #Compute fluxes
        fluxhalf = calc_flux(lamb, ue, E, del, F)

        #compute rhs
        for ic = 2:ncell+1
            rhs[1,ic] = B[ic]/(2*del[ic]) - del[ic]*uet[ic]/ue[ic] - (E[ic] + 1.)*del[ic]*uex[ic]
            rhs[2,ic] = S[ic]/del[ic] - 2*E[ic]*del[ic]*uet[ic]/ue[ic] - 2*F[ic]*del[ic]*uex[ic]
        end

        ic = 2

        unk[1,ic] = 0.5*(unk[1,ic] + unkh[1,ic]) -
        (0.5*dt/dx)*(flux[1,1,ic] - flux[1,1,ic-1] - flux[2,1,ic] +
        2*flux[2,1,ic+1] - flux[2,1,ic+2] + fluxhalf[1,1,ic] -
        fluxhalf[1,1,ic-1] + fluxhalf[2,1,ic+1] - fluxhalf[2,1,ic]) +
        0.5*dt*rhs[1,ic]

        unk[2,ic] = 0.5*(unk[2,ic] + unkh[2,ic]) -
        (0.5*dt/dx)*(flux[1,2,ic] - flux[1,2,ic-1] - flux[2,2,ic] +
        2*flux[2,2,ic+1] - flux[2,2,ic+2] + fluxhalf[1,2,ic] -
        fluxhalf[1,2,ic-1] + fluxhalf[2,2,ic+1] - fluxhalf[2,2,ic]) +
        0.5*dt*rhs[2,ic]

        for ic = 3:ncell

            unk[1,ic] = 0.5*(unk[1,ic] + unkh[1,ic]) -
            (0.5*dt/dx)*(flux[1,1,ic] - 2*flux[1,1,ic-1] +
            flux[1,1,ic-2] - flux[2,1,ic] + 2*flux[2,1,ic+1] -
            flux[2,1,ic+2] + fluxhalf[1,1,ic] -fluxhalf[1,1,ic-1] +
            fluxhalf[2,1,ic+1] - fluxhalf[2,1,ic]) + 0.5*dt*rhs[1,ic]

            unk[2,ic] = 0.5*(unk[2,ic] + unkh[2,ic]) -
            (0.5*dt/dx)*(flux[1,2,ic] - 2*flux[1,2,ic-1] +
            flux[1,2,ic-2] - flux[2,2,ic] + 2*flux[2,2,ic+1] -
            flux[2,2,ic+2] + fluxhalf[1,2,ic] - fluxhalf[1,2,ic-1] +
            fluxhalf[2,2,ic+1] - fluxhalf[2,2,ic]) + 0.5*dt*rhs[2,ic]

        end
        ic = ncell+1

        unk[1,ic] = 0.5*(unk[1,ic] + unkh[1,ic]) -
        (0.5*dt/dx)*(flux[1,1,ic] - 2*flux[1,1,ic-1] + flux[1,1,ic-2]
        +fluxhalf[1,1,ic] - fluxhalf[1,1,ic-1]) + 0.5*dt*rhs[1,ic]

        unk[2,ic] = 0.5*(unk[2,ic] + unkh[2,ic]) -
        (0.5*dt/dx)*(flux[1,2,ic] - 2*flux[1,2,ic-1] + flux[1,2,ic-2]
        +fluxhalf[1,2,ic] - fluxhalf[1,2,ic-1]) + 0.5*dt*rhs[2,ic]

        t = t + dt
        if t > ttime
            break
        end

        for ic = 2:ncell+1
            crit[ic] = abs((del[ic+1] - del[ic])/(del[ic] - del[ic-1]))
            if abs(crit[ic]) > 10.
                sepflag = 1
                println(t, " ",x[ic]," ", ue[ic]," ", crit[ic])
                break
            end
        end
        if t > ttime || sepflag == 1
            break
        end
    end
t, x, del, E
end

