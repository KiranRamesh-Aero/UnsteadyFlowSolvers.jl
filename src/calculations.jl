#Function for estimating a problem's time step

#Simple linear interpolation function
function interp(x1 ::Float64, x2 :: Float64, y1 :: Float64, y2 :: Float64, x::Float64)
    y = y1 + (y2 - y1)*(x - x1)/(x2 - x1)
    return y
end

    
function find_tstep(kin:: Array{CosDef})
    dtstar = 15
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

function find_tstep(kin :: BendingDef)
    dtstar = 15
    amp = evaluate(kin.spl, kin.spl.t[end])
    dt_tmp = 0.015*0.2/(kin.k*amp)
    dtstar = minimum([dt_tmp dtstar])
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
function update_a2a3adot(surf::Vector{TwoDSurf},dt)
    for i = 1:length(surf)
        for ia = 2:3
            surf[i].aterm[ia] = trapz(surf[i].downwash.*cos(ia*surf[i].theta),surf[i].theta)
            surf[i].aterm[ia] = 2.*surf[i].aterm[ia]/(surf[i].uref*pi)
        end
        surf[i].a0dot[1] = (surf[i].a0[1] - surf[i].a0prev[1])/dt
        for ia = 1:3
            surf[i].adot[ia] = (surf[i].aterm[ia]-surf[i].aprev[ia])/dt
        end
    end
    return surf
end

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

function update_adot(surf::Vector{TwoDSurf},dt)
    for i = 1:length(surf)
        surf[i].a0dot[1] = (surf[i].a0[1] - surf[i].a0prev[1])/dt
        for ia = 1:3
            surf[i].adot[ia] = (surf[i].aterm[ia]-surf[i].aprev[ia])/dt
        end
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
    surf.uind[1:surf.ndiv], surf.wind[1:surf.ndiv] = ind_vel([curfield.tev; curfield.lev; curfield.extv], surf.bnd_x, surf.bnd_z)
    return surf
end

function update_indbound(surf::Vector{TwoDSurf}, curfield::TwoDFlowFieldMultSurf, shed_ind::Vector{Vector{Int}})

    
    for i = 1:curfield.nsurf
        uind = zeros(surf[i].ndiv)
        wind = zeros(surf[i].ndiv)
   
        for j = 1:surf[i].ndiv
            surf[i].uind[j] = 0.
            surf[i].wind[j] = 0.
        end
        
        for n = 1:length(curfield.tev)
            uind, wind = ind_vel(curfield.tev[n], surf[i].bnd_x, surf[i].bnd_z)
            for j = 1:surf[i].ndiv
                surf[i].uind[j] += uind[j]
                surf[i].wind[j] += wind[j]
            end
        end
        
        for n = 1:length(curfield.lev)
            tempv = TwoDVort[]
            for ii = 1:curfield.nsurf
                if ii in shed_ind[n]
                    push!(tempv,curfield.lev[n][ii])
                end
            end       
            uind, wind = ind_vel(tempv, surf[i].bnd_x, surf[i].bnd_z)
            for j = 1:surf[i].ndiv
                surf[i].uind[j] += uind[j]
                surf[i].wind[j] += wind[j]
            end
        end

        for n = 1:curfield.nsurf
            if n != i
                #Calculate bv for this surface
                #update_bv(surf[n])
                uind, wind = ind_vel(surf[n].bv, surf[i].bnd_x, surf[i].bnd_z)
                for j = 1:surf[i].ndiv
                    surf[i].uind[j] += uind[j]
                    surf[i].wind[j] += wind[j]
                end
            end
        end

    end
    return surf
end

function update_indbound(surf::Vector{TwoDSurf}, curfield::ThreeDFieldStrip, shed_ind::Vector{Vector{Int}})

    #The vortices associated with each spanwise location dont affect each other. Can be added
    for i = 1:curfield.nspan
        uind = zeros(surf[i].ndiv)
        wind = zeros(surf[i].ndiv)
        
        for j = 1:surf[i].ndiv
            surf[i].uind[j] = 0.
            surf[i].wind[j] = 0.
        end
        
        for n = 1:length(curfield.tev)
            uind, wind = ind_vel([curfield.tev[n][i]], surf[i].bnd_x, surf[i].bnd_z)
            for j = 1:surf[i].ndiv
                surf[i].uind[j] += uind[j]
                surf[i].wind[j] += wind[j]
            end
        end
        
        for n = 1:length(curfield.lev)
            if i in shed_ind[n]
                uind, wind = ind_vel([curfield.lev[n][i]], surf[i].bnd_x, surf[i].bnd_z)
                for j = 1:surf[i].ndiv
                    surf[i].uind[j] += uind[j]
                    surf[i].wind[j] += wind[j]
                end
            end       
        end

        # for n = 1:curfield.nsurf
        #     if n != i
        #         #Calculate bv for this surface
        #         #update_bv(surf[n])
        #         uind, wind = ind_vel(surf[n].bv, surf[i].bnd_x, surf[i].bnd_z)
        #         for j = 1:surf[i].ndiv
        #             surf[i].uind[j] += uind[j]
        #             surf[i].wind[j] += wind[j]
        #         end
        #     end
        # end

    end
    return surf
end


function update_indbound(surf::Vector{TwoDSurf}, curfield::TwoDFlowFieldMultSurf)
    
    for i = 1:curfield.nsurf
        uind = zeros(surf[i].ndiv)
        wind = zeros(surf[i].ndiv)
   
        for j = 1:surf[i].ndiv
            surf[i].uind[j] = 0.
            surf[i].wind[j] = 0.
        end
        
        for n = 1:length(curfield.tev)
            uind, wind = ind_vel(curfield.tev[n], surf[i].bnd_x, surf[i].bnd_z)
            for j = 1:surf[i].ndiv
                surf[i].uind[j] += uind[j]
                surf[i].wind[j] += wind[j]
            end
        end

        #Add influence of other surfaces
        for n = 1:curfield.nsurf
            if n != i
                #Calculate bv for this surface
                #update_bv(surf[n])
                uind, wind = ind_vel(surf[n].bv, surf[i].bnd_x, surf[i].bnd_z)
                for j = 1:surf[i].ndiv
                    surf[i].uind[j] += uind[j]
                    surf[i].wind[j] += wind[j]
                end
            end
        end
    end            

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
function update_downwash(surf::Vector{TwoDSurf}, vels::Vector{Float64})
    for i = 1:length(surf)
        for ib = 1:surf[i].ndiv
            surf[i].downwash[ib] = -(surf[i].kinem.u + vels[1])*sin(surf[i].kinem.alpha) - surf[i].uind[ib]*sin(surf[i].kinem.alpha) + (surf[i].kinem.hdot - vels[2])*cos(surf[i].kinem.alpha) - surf[i].wind[ib]*cos(surf[i].kinem.alpha) - surf[i].kinem.alphadot*(surf[i].x[ib] - surf[i].pvt*surf[i].c) + surf[i].cam_slope[ib]*(surf[i].uind[ib]*cos(surf[i].kinem.alpha) + (surf[i].kinem.u + vels[1])*cos(surf[i].kinem.alpha) + (surf[i].kinem.hdot - vels[2])*sin(surf[i].kinem.alpha) - surf[i].wind[ib]*sin(surf[i].kinem.alpha))
        end
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

function update_a0anda1(surf::Vector{TwoDSurf})
    for i =1:length(surf)
        surf[i].a0[1] = trapz(surf[i].downwash,surf[i].theta)
        surf[i].aterm[1] = trapz(surf[i].downwash.*cos(surf[i].theta),surf[i].theta)
        surf[i].a0[1] = -surf[i].a0[1]/(surf[i].uref*pi)
        surf[i].aterm[1] = 2.*surf[i].aterm[1]/(surf[i].uref*pi)
    end
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

function update_a2toan(surf::Vector{TwoDSurf})
    for i = 1:length(surf)
        for ia = 2:surf[i].naterm
            surf[i].aterm[ia] = trapz(surf[i].downwash.*cos(ia*surf[i].theta),surf[i].theta)
            surf[i].aterm[ia] = 2.*surf[i].aterm[ia]/(surf[i].uref*pi)
        end
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

function update_externalvel(curfield::TwoDFlowFieldMultSurf, t)
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

function update_externalvel(curfield::ThreeDFieldStrip, t)
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
	elseif (typeof(surf.kindef.alpha) == LinearDef)
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
    elseif (typeof(surf.kindef.u) == LinearDef)
        surf.kinem.u = surf.kindef.u(t)*surf.uref
        surf.kinem.udot = ForwardDiff.derivative(surf.kindef.u,t)*surf.uref*surf.uref/surf.c
    end
    # ---------------------------------------------------------------------------------------------
    return surf
end
# END kinem function
# ---------------------------------------------------------------------------------------------
# Function updating the dimensional kinematic parameters
function update_kinem(surf::Vector{TwoDSurf}, t)

    for i = 1:length(surf)
        # Pitch kinematics
        if (typeof(surf[i].kindef.alpha) == EldUpDef)
            surf[i].kinem.alpha = surf[i].kindef.alpha(t)
            surf[i].kinem.alphadot = ForwardDiff.derivative(surf[i].kindef.alpha,t)*surf[i].uref/surf[i].c
        elseif (typeof(surf[i].kindef.alpha) == EldUptstartDef)
            surf[i].kinem.alpha = surf[i].kindef.alpha(t)
            surf[i].kinem.alphadot = ForwardDiff.derivative(surf[i].kindef.alpha,t)*surf[i].uref/surf[i].c
        elseif (typeof(surf[i].kindef.alpha) == EldRampReturnDef)
            surf[i].kinem.alpha = surf[i].kindef.alpha(t)
            surf[i].kinem.alphadot = ForwardDiff.derivative(surf[i].kindef.alpha,t)*surf[i].uref/surf[i].c
        elseif (typeof(surf[i].kindef.alpha) == ConstDef)
            surf[i].kinem.alpha = surf[i].kindef.alpha(t)
            surf[i].kinem.alphadot = 0.
        elseif (typeof(surf[i].kindef.alpha) == SinDef)
            surf[i].kinem.alpha = surf[i].kindef.alpha(t)
            surf[i].kinem.alphadot = ForwardDiff.derivative(surf[i].kindef.alpha,t)*surf[i].uref/surf[i].c
        elseif (typeof(surf[i].kindef.alpha) == CosDef)
            surf[i].kinem.alpha = surf[i].kindef.alpha(t)
            surf[i].kinem.alphadot = ForwardDiff.derivative(surf[i].kindef.alpha,t)*surf[i].uref/surf[i].c
        end
        # ---------------------------------------------------------------------------------------------
        
        # Plunge kinematics
        if (typeof(surf[i].kindef.h) == EldUpDef)
            surf[i].kinem.h = surf[i].kindef.h(t)*surf[i].c
            surf[i].kinem.hdot = ForwardDiff.derivative(surf[i].kindef.h,t)*surf[i].uref
        elseif (typeof(surf[i].kindef.h) == EldUptstartDef)
            surf[i].kinem.h = surf[i].kindef.h(t)*surf[i].c
            surf[i].kinem.hdot = ForwardDiff.derivative(surf[i].kindef.h,t)*surf[i].uref
        elseif (typeof(surf[i].kindef.h) == EldUpIntDef)
            surf[i].kinem.h = surf[i].kindef.h(t)*surf[i].c
            surf[i].kinem.hdot = ForwardDiff.derivative(surf[i].kindef.h,t)*surf[i].uref
        elseif (typeof(surf[i].kindef.h) == EldUpInttstartDef)
            surf[i].kinem.h = surf[i].kindef.h(t)*surf[i].c
            surf[i].kinem.hdot = ForwardDiff.derivative(surf[i].kindef.h,t)*surf[i].uref
        elseif (typeof(surf[i].kindef.h) == EldRampReturnDef)
            surf[i].kinem.h = surf[i].kindef.h(t)*surf[i].c
            surf[i].kinem.hdot = ForwardDiff.derivative(surf[i].kindef.h,t)*surf[i].uref
        elseif (typeof(surf[i].kindef.h) == ConstDef)
            surf[i].kinem.h = surf[i].kindef.h(t)*surf[i].c
            surf[i].kinem.hdot = 0.
        elseif (typeof(surf[i].kindef.h) == SinDef)
            surf[i].kinem.h = surf[i].kindef.h(t)*surf[i].c
            surf[i].kinem.hdot = ForwardDiff.derivative(surf[i].kindef.h,t)*surf[i].uref
        elseif (typeof(surf[i].kindef.h) == CosDef)
            surf[i].kinem.h = surf[i].kindef.h(t)*surf[i].c
            surf[i].kinem.hdot = ForwardDiff.derivative(surf[i].kindef.h,t)*surf[i].uref
        end
        # ---------------------------------------------------------------------------------------------
        
        # Forward velocity
        if (typeof(surf[i].kindef.u) == EldUpDef)
            surf[i].kinem.u = surf[i].kindef.u(t)*surf[i].uref
            surf[i].kinem.udot = ForwardDiff.derivative(surf[i].kindef.u,t)*surf[i].uref*surf[i].uref/surf[i].c
        elseif (typeof(surf[i].kindef.u) == EldRampReturnDef)
            surf[i].kinem.u, surf[i].kinem.udot = surf[i].kindef.u(t)
            surf[i].kinem.u = surf[i].kinem.u*surf[i].uref
            surf[i].kinem.udot = surf[i].kinem.udot*surf[i].uref*surf[i].uref/surf[i].c
        elseif (typeof(surf[i].kindef.u) == ConstDef)
            surf[i].kinem.u = surf[i].kindef.u(t)*surf[i].uref
        surf[i].kinem.udot = 0.
        end
    end
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
function update_kinem(surf::ThreeDSurfVR, t)

    # Pitch kinematics
    if (typeof(surf.kindef.alpha) == EldUpDef)
        surf.kinem.alpha = surf.kindef.alpha(t)
        surf.kinem.alphadot = ForwardDiff.derivative(surf.kindef.alpha,t)*surf.uref/surf.cref
    elseif (typeof(surf.kindef.alpha) == EldUptstartDef)
        surf.kinem.alpha = surf.kindef.alpha(t)
        surf.kinem.alphadot = ForwardDiff.derivative(surf.kindef.alpha,t)*surf.uref/surf.cref
    elseif (typeof(surf.kindef.alpha) == EldRampReturnDef)
        surf.kinem.alpha = surf.kindef.alpha(t)
        surf.kinem.alphadot = ForwardDiff.derivative(surf.kindef.alpha,t)*surf.uref/surf.cref
    elseif (typeof(surf.kindef.alpha) == ConstDef)
        surf.kinem.alpha = surf.kindef.alpha(t)
        surf.kinem.alphadot = 0.
    elseif (typeof(surf.kindef.alpha) == SinDef)
        surf.kinem.alpha = surf.kindef.alpha(t)
        surf.kinem.alphadot = ForwardDiff.derivative(surf.kindef.alpha,t)*surf.uref/surf.cref
    elseif (typeof(surf.kindef.alpha) == CosDef)
        surf.kinem.alpha = surf.kindef.alpha(t)
        surf.kinem.alphadot = ForwardDiff.derivative(surf.kindef.alpha,t)*surf.uref/surf.cref
    end
    # ---------------------------------------------------------------------------------------------

    # Plunge kinematics
    if (typeof(surf.kindef.h) == EldUpDef)
        surf.kinem.h = surf.kindef.h(t)*surf.cref
        surf.kinem.hdot = ForwardDiff.derivative(surf.kindef.h,t)*surf.uref
    elseif (typeof(surf.kindef.h) == EldUptstartDef)
        surf.kinem.h = surf.kindef.h(t)*surf.cref
        surf.kinem.hdot = ForwardDiff.derivative(surf.kindef.h,t)*surf.uref
    elseif (typeof(surf.kindef.h) == EldUpIntDef)
            surf.kinem.h = surf.kindef.h(t)*surf.cref
            surf.kinem.hdot = ForwardDiff.derivative(surf.kindef.h,t)*surf.uref
    elseif (typeof(surf.kindef.h) == EldUpInttstartDef)
        surf.kinem.h = surf.kindef.h(t)*surf.cref
        surf.kinem.hdot = ForwardDiff.derivative(surf.kindef.h,t)*surf.uref
    elseif (typeof(surf.kindef.h) == EldRampReturnDef)
        surf.kinem.h = surf.kindef.h(t)*surf.cref
        surf.kinem.hdot = ForwardDiff.derivative(surf.kindef.h,t)*surf.uref
    elseif (typeof(surf.kindef.h) == ConstDef)
        surf.kinem.h = surf.kindef.h(t)*surf.cref
        surf.kinem.hdot = 0.
    elseif (typeof(surf.kindef.h) == SinDef)
      surf.kinem.h = surf.kindef.h(t)*surf.cref
      surf.kinem.hdot = ForwardDiff.derivative(surf.kindef.h,t)*surf.uref
    elseif (typeof(surf.kindef.h) == CosDef)
      surf.kinem.h = surf.kindef.h(t)*surf.cref
      surf.kinem.hdot = ForwardDiff.derivative(surf.kindef.h,t)*surf.uref
    end
    # ---------------------------------------------------------------------------------------------

    # Forward velocity
    if (typeof(surf.kindef.u) == EldUpDef)
        surf.kinem.u = surf.kindef.u(t)*surf.uref
        surf.kinem.udot = ForwardDiff.derivative(surf.kindef.u,t)*surf.uref*surf.uref/surf.cref
    elseif (typeof(surf.kindef.u) == EldRampReturnDef)
        surf.kinem.u, surf.kinem.udot = surf.kindef.u(t)
        surf.kinem.u = surf.kinem.u*surf.uref
        surf.kinem.udot = surf.kinem.udot*surf.uref*surf.uref/surf.cref
    elseif (typeof(surf.kindef.u) == ConstDef)
        surf.kinem.u = surf.kindef.u(t)*surf.uref
        surf.kinem.udot = 0.
    end
    # ---------------------------------------------------------------------------------------------
    return surf
end

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

function update_bv(surf::Vector{TwoDSurf})
    for i = 1:length(surf)
        gamma = zeros(surf[i].ndiv)
        for ib = 1:surf[i].ndiv
            gamma[ib] = (surf[i].a0[1]*(1 + cos(surf[i].theta[ib])))
            for ia = 1:surf[i].naterm
                gamma[ib] = gamma[ib] + surf[i].aterm[ia]*sin(ia*surf[i].theta[ib])*sin(surf[i].theta[ib])
            end
            gamma[ib] = gamma[ib]*surf[i].uref*surf[i].c
        end
                
        for ib = 2:surf[i].ndiv
            surf[i].bv[ib-1].s = (gamma[ib]+gamma[ib-1])*(surf[i].theta[2]-surf[i].theta[1])/2.
            surf[i].bv[ib-1].x = (surf[i].bnd_x[ib] + surf[i].bnd_x[ib-1])/2.
            surf[i].bv[ib-1].z = (surf[i].bnd_z[ib] + surf[i].bnd_z[ib-1])/2.
        end
    end
    return surf
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

function wakeroll(surf::Vector{TwoDSurf}, curfield::TwoDFlowFieldMultSurf, dt, shed_ind)
    #Clean induced velocities
    for i = 1:length(curfield.tev)
        for j = 1:curfield.nsurf
            curfield.tev[i][j].vx = 0
            curfield.tev[i][j].vz = 0
        end
    end
    
    for i = 1:length(curfield.lev)
        for j = 1:curfield.nsurf
            if j in shed_ind[i]
                curfield.lev[i][j].vx = 0
                curfield.lev[i][j].vz = 0
            end
        end
    end
    
    #Velocities induced by free vortices on each other
    mutual_ind([curfield.tev; curfield.lev], curfield.nsurf, shed_ind)

    #Add the influence of velocities induced by bound vortices
    for i = 1:length(curfield.tev)
        for k = 1:curfield.nsurf
            for j = 1:curfield.nsurf
                utemp, wtemp = ind_vel(surf[j].bv, curfield.tev[i][k].x, curfield.tev[i][k].z)
                curfield.tev[i][k].vx += utemp[1]
                curfield.tev[i][k].vz += wtemp[1]
            end
        end
    end
    
    for i = 1:length(curfield.lev)
        for k = 1:curfield.nsurf
            if k in shed_ind[i]
                for j = 1:curfield.nsurf
                    utemp, wtemp = ind_vel(surf[j].bv, curfield.lev[i][k].x, curfield.lev[i][k].z)
                    curfield.lev[i][k].vx += utemp[1]
                    curfield.lev[i][k].vz += wtemp[1]
                end
            end
        end
    end
        
    #Convect free vortices with their induced velocities
    for i = 1:length(curfield.tev)
        for j = 1:curfield.nsurf 
            curfield.tev[i][j].x += dt*curfield.tev[i][j].vx
            curfield.tev[i][j].z += dt*curfield.tev[i][j].vz
        end
    end
    for i = 1:length(curfield.lev)
        for j = 1:curfield.nsurf 
            if j in shed_ind[i]
                curfield.lev[i][j].x += dt*curfield.lev[i][j].vx
                curfield.lev[i][j].z += dt*curfield.lev[i][j].vz
            end
        end
    end
    return curfield
end

function wakeroll(surf::Vector{TwoDSurf}, curfield::ThreeDFieldStrip, dt, shed_ind)
    #Clean induced velocities
    for i = 1:length(curfield.tev)
        for j = 1:curfield.nspan
            curfield.tev[i][j].vx = 0
            curfield.tev[i][j].vz = 0
        end
    end
    
    for i = 1:length(curfield.lev)
        for j = 1:curfield.nspan
            if j in shed_ind[i]
                curfield.lev[i][j].vx = 0
                curfield.lev[i][j].vz = 0
            end
        end
    end
    
    #Velocities induced by free vortices on each other in each strip
    mutual_ind_llt([curfield.tev; curfield.lev], curfield.nspan, shed_ind)
        
    #Add the influence of velocities induced by bound vortices on that stip
    for i = 1:curfield.nspan
        for iv = 1:length(curfield.tev)
            utemp, wtemp = ind_vel(surf[i].bv, curfield.tev[iv][i].x, curfield.tev[iv][i].z)
            curfield.tev[iv][i].vx += utemp[1]
            curfield.tev[iv][i].vz += wtemp[1]
        end
    end

    for i = 1:curfield.nspan
        for iv = 1:length(curfield.lev)
            if i in shed_ind[iv]
                utemp, wtemp = ind_vel(surf[i].bv, curfield.lev[iv][i].x, curfield.lev[iv][i].z)
                curfield.lev[iv][i].vx += utemp[1]
                curfield.lev[iv][i].vz += wtemp[1]
            end
        end
    end
    
    #Convect free vortices with their induced velocities
    for i = 1:curfield.nspan
        for iv = 1:length(curfield.tev) 
            curfield.tev[iv][i].x += dt*curfield.tev[iv][i].vx
            curfield.tev[iv][i].z += dt*curfield.tev[iv][i].vz
        end
    end
    for i = 1:curfield.nspan
        for iv = 1:length(curfield.lev) 
            if i in shed_ind[iv]
                curfield.lev[iv][i].x += dt*curfield.lev[iv][i].vx
                curfield.lev[iv][i].z += dt*curfield.lev[iv][i].vz
            end
        end
    end
    return curfield
end

function parwakeroll(surf::Vector{TwoDSurf}, curfield::ThreeDFieldStrip, dt, shed_ind)
    #Clean induced velocities
    for i = 1:length(curfield.tev)
        for j = 1:curfield.nspan
            curfield.tev[i][j].vx = 0
            curfield.tev[i][j].vz = 0
        end
    end
    
    for i = 1:length(curfield.lev)
        for j = 1:curfield.nspan
            if j in shed_ind[i]
                curfield.lev[i][j].vx = 0
                curfield.lev[i][j].vz = 0
            end
        end
    end
    
    #Velocities induced by free vortices on each other in each strip
    parmutual_ind_llt([curfield.tev; curfield.lev], curfield.nspan, shed_ind)
        
    #Add the influence of velocities induced by bound vortices on that stip
    @sync @parallel for i = 1:curfield.nspan
        for iv = 1:length(curfield.tev)
            utemp, wtemp = ind_vel(surf[i].bv, curfield.tev[iv][i].x, curfield.tev[iv][i].z)
            curfield.tev[iv][i].vx += utemp[1]
            curfield.tev[iv][i].vz += wtemp[1]
        end
    end

    for i = 1:curfield.nspan
        for iv = 1:length(curfield.lev)
            if i in shed_ind[iv]
                utemp, wtemp = ind_vel(surf[i].bv, curfield.lev[iv][i].x, curfield.lev[iv][i].z)
                curfield.lev[iv][i].vx += utemp[1]
                curfield.lev[iv][i].vz += wtemp[1]
            end
        end
    end
    
    #Convect free vortices with their induced velocities
    for i = 1:curfield.nspan
        for iv = 1:length(curfield.tev) 
            curfield.tev[iv][i].x += dt*curfield.tev[iv][i].vx
            curfield.tev[iv][i].z += dt*curfield.tev[iv][i].vz
        end
    end
    for i = 1:curfield.nspan
        for iv = 1:length(curfield.lev) 
            if i in shed_ind[iv]
                curfield.lev[iv][i].x += dt*curfield.lev[iv][i].vx
                curfield.lev[iv][i].z += dt*curfield.lev[iv][i].vz
            end
        end
    end
    return curfield
end



function wakeroll(surf::Vector{TwoDSurf}, curfield::TwoDFlowFieldMultSurf, dt)
    #Without LEVs
    #Clean induced velocities
    for i = 1:length(curfield.tev)
        for j = 1:curfield.nsurf
            curfield.tev[i][j].vx = 0
            curfield.tev[i][j].vz = 0
        end
    end

    #Velocities induced by free vortices on each other
    mutual_ind(curfield.tev, curfield.nsurf)

    #Add the influence of velocities induced by bound vortices
    for i = 1:length(curfield.tev)
        for j = 1:curfield.nsurf
            utemp, wtemp = ind_vel(surf[j].bv, curfield.tev[i][j].x, curfield.tev[i][j].z)
            curfield.tev[i][j].vx += utemp[1]
            curfield.tev[i][j].vz += wtemp[1]
        end
    end
    
    #Convect free vortices with their induced velocities
    for i = 1:length(curfield.tev)
        for j = 1:curfield.nsurf 
            curfield.tev[i][j].x += dt*curfield.tev[i][j].vx
            curfield.tev[i][j].z += dt*curfield.tev[i][j].vz
        end
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

function place_tev(surf::Vector{TwoDSurf}, field::TwoDFlowFieldMultSurf,dt)
    ntev = length(field.tev)
    tempv = TwoDVort[]
    for i = 1:field.nsurf
        if ntev == 0
            xloc = surf[i].bnd_x[surf[i].ndiv] + 0.5*surf[i].kinem.u*dt
            zloc = surf[i].bnd_z[surf[i].ndiv]
        else
            xloc = surf[i].bnd_x[surf[i].ndiv]+(1./3.)*(field.tev[ntev][i].x - surf[i].bnd_x[surf[i].ndiv])
            zloc = surf[i].bnd_z[surf[i].ndiv]+(1./3.)*(field.tev[ntev][i].z - surf[i].bnd_z[surf[i].ndiv])
        end
        push!(tempv,TwoDVort(xloc,zloc,0.,0.02*surf[i].c,0.,0.))
    end
    push!(field.tev, tempv)
    return field
end

function place_tev(surf::Vector{TwoDSurf}, field::ThreeDFieldStrip,dt)
    ntev = length(field.tev)
    tempv = TwoDVort[]
    for i = 1:field.nspan
        if ntev == 0
            xloc = surf[i].bnd_x[surf[i].ndiv] + 0.5*surf[i].kinem.u*dt
            zloc = surf[i].bnd_z[surf[i].ndiv]
        else
            xloc = surf[i].bnd_x[surf[i].ndiv]+(1./3.)*(field.tev[ntev][i].x - surf[i].bnd_x[surf[i].ndiv])
            zloc = surf[i].bnd_z[surf[i].ndiv]+(1./3.)*(field.tev[ntev][i].z - surf[i].bnd_z[surf[i].ndiv])
        end
        push!(tempv,TwoDVort(xloc,zloc,0.,0.02*surf[i].c,0.,0.))
    end
    push!(field.tev, tempv)
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

function place_lev(surf::Vector{TwoDSurf}, field::TwoDFlowFieldMultSurf,dt, shed_ind::Vector{Vector{Int}})
    
    nlev = length(field.lev)
    tempv = TwoDVort[]
    dummyV = TwoDVort(0.,0.,0.,0.,0.,0.)
    #These dummy vortices are placed in the LEV array when that surface doesnt shed levs. These wont be used in any calculations but they still occupy memory. 
    
    for i = 1:field.nsurf
        if i in shed_ind[nlev+1]
            le_vel_x = surf[i].kinem.u - surf[i].kinem.alphadot*sin(surf[i].kinem.alpha)*surf[i].pvt*surf[i].c + surf[i].uind[1]
            le_vel_z = -surf[i].kinem.alphadot*cos(surf[i].kinem.alpha)*surf[i].pvt*surf[i].c- surf[i].kinem.hdot + surf[i].wind[1]
            
            if (surf[i].levflag[1] == 0)
                xloc = surf[i].bnd_x[1] + 0.5*le_vel_x*dt
                zloc = surf[i].bnd_z[1] + 0.5*le_vel_z*dt
            else
                xloc = surf[i].bnd_x[1]+(1./3.)*(field.lev[nlev][i].x - surf[i].bnd_x[1])
                zloc = surf[i].bnd_z[1]+(1./3.)*(field.lev[nlev][i].z - surf[i].bnd_z[1])
            end
            push!(tempv,TwoDVort(xloc,zloc,0.,0.02*surf[i].c,0.,0.))
        else
            push!(tempv,dummyV)
        end
    end

    push!(field.lev, tempv)
    return field
end

function place_lev(surf::Vector{TwoDSurf}, field::ThreeDFieldStrip, dt :: Float64)
    
    nlev = length(field.lev)
    tempv = TwoDVort[]
    dummyV = TwoDVort(0.,0.,0.,0.,0.,0.)
    #These dummy vortices are placed in the LEV array when that surface doesnt shed levs. These wont be used in any calculations but they still occupy memory. 
    
    for i = 1:field.nspan
        if (surf[i].levflag[1] == 0)
            le_vel_x = surf[i].kinem.u - surf[i].kinem.alphadot*sin(surf[i].kinem.alpha)*surf[i].pvt*surf[i].c + surf[i].uind[1]
            le_vel_z = -surf[i].kinem.alphadot*cos(surf[i].kinem.alpha)*surf[i].pvt*surf[i].c- surf[i].kinem.hdot + surf[i].wind[1]
            
            xloc = surf[i].bnd_x[1] + 0.5*le_vel_x*dt
            zloc = surf[i].bnd_z[1] + 0.5*le_vel_z*dt
        else
            xloc = surf[i].bnd_x[1]+(1./3.)*(field.lev[nlev][i].x - surf[i].bnd_x[1])
            zloc = surf[i].bnd_z[1]+(1./3.)*(field.lev[nlev][i].z - surf[i].bnd_z[1])
        end
        push!(tempv,TwoDVort(xloc,zloc,0.,0.02*surf[i].c,0.,0.))
    end
    
    push!(field.lev, tempv)
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
#Calculates first approximation of surface speed for NACA 00XX airfoil at fsep location (defined as a fraction of c)
function surfspeed(tau::Float64, fsep::Float64, surf::TwoDSurf)

	#Coefficients for NACA airfoil
	b1 = 1.48450
	b2 = -0.63
	b4 = -1.758
	b6 = 1.4215
	b8 = -0.5075
	#derivative of thickness function
	dT = tau.*(b1./(2.*sqrt(fsep)).+b2.+2.*b4.*fsep.+3.*b6.*fsep.^2.+4.*b8.*fsep.^3)

	#Delete this after singularities @ LE and TE are removed
	if(fsep <0.02)
	fsep == 0.02
	elseif(fsep > 0.98)
	fsep = 0.98
	end
	
	#First-order speed due to thickness
	q1t = tau./pi.*(1./tau.*dT.*log(fsep./(1.-fsep))+b1./sqrt(fsep).*log((1.+sqrt(fsep))./sqrt(fsep)).-2.*b4.-3./2.*b6.-4./3.*b8.-(3.*b6.+2.*b8).*fsep.-4.*b8.*fsep.^2)
		
		#Riegel's rule
		ro = (1.1019*tau^2*surf.c)/2
		cos_eta = sqrt(fsep./(fsep.+ro./2))
		q1t=q1t.*cos_eta
		
	#First-order speed due to camber
	sumAn = 0.
	fsep_theta = acos(1-2*fsep)
	for ia = 1:surf.naterm
        sumAn = sumAn + surf.aterm[ia]*sin(ia*fsep_theta)*sin(fsep_theta/2)
    end
	q1c = sqrt(2.)./sqrt(1-cos(fsep_theta)+ro)*(surf.a0[1]*cos(fsep_theta/2)+sumAn)
	
	#Total speed
	u_sh = surf.kinem.u*(1+q1t+q1c)
	
	u_sh, q1t, q1c
end

# Places a vortex at separation point location
function place_spv(surf::TwoDSurf,field::TwoDFlowField,dt, fsep::Float64, tau::Float64)
	nlev = length(field.lev)
	
	#Find idx corresponding to fsep location
		#x is defined from 0 to c
		#fsep is defined as a fraction of c
	tmp = abs(surf.x.-fsep.*surf.c) #subtract and find the smallest value
	idx = 1
	for a = 1 : length(surf.x)-1
		if(tmp[idx+1]<tmp[idx])
			idx+=1
		end
	end
		#This is more neat solution, but it crashes at some conditions for some reason.
			#E.g. for alpha = 18.0 deg
			#NACA0012 sepdef = SeparationParams(16.2,1.52,3.21,"Sheng")
		#val = minimum(tmp)
		#idx = find(tmp -> tmp==val, tmp) #find idx of the smallest value of tmp
		#idx = idx[1] #in case there are multiple minimum values	

#VELOCITY, both formulas should give similar results for steady state
	#le_vel_x = surfspeed(tau, fsep, surf)*cos(surf.kindef.alpha.amp) #contribution from hdot and alphadot needs to be added
	#le_vel_z = surfspeed(tau, fsep, surf)*sin(surf.kindef.alpha.amp) #contribution from hdot and alphadot needs to be added
	le_vel_x = surf.kinem.u - surf.kinem.alphadot*sin(surf.kinem.alpha)*abs(fsep-surf.pvt)*surf.c + surf.uind[idx]
    le_vel_z = -surf.kinem.alphadot*cos(surf.kinem.alpha)*abs(fsep-surf.pvt)*surf.c- surf.kinem.hdot + surf.wind[idx]

#POSITION
    if (surf.levflag[1] == 0)
        xloc = surf.bnd_x[idx] + 0.5*le_vel_x*dt
        zloc = surf.bnd_z[idx] + 0.5*le_vel_z*dt
   else
       xloc = surf.bnd_x[idx]+(1./3.)*(field.lev[nlev].x - surf.bnd_x[idx]) 
       zloc = surf.bnd_z[idx]+(1./3.)*(field.lev[nlev].z - surf.bnd_z[idx])
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

function mutual_ind(vorts::Vector{Vector{TwoDVort}}, nsurf :: Int)
    for i = 1:length(vorts)
        for k = 1:nsurf
            for j = i+1:length(vorts)
                for l = 1:nsurf
                    
                    dx = vorts[i][k].x - vorts[j][l].x
                    dz = vorts[i][k].z - vorts[j][l].z
                    #source- tar
                    dsq = dx*dx + dz*dz
                    
                    magitr = 1./(2*pi*sqrt(vorts[j][l].vc*vorts[j][l].vc*vorts[j][l].vc*vorts[j][l].vc + dsq*dsq))
                    magjtr = 1./(2*pi*sqrt(vorts[i][k].vc*vorts[i][k].vc*vorts[i][k].vc*vorts[i][k].vc + dsq*dsq))
                    
                    vorts[j][l].vx -= dz * vorts[i][k].s * magjtr
                    vorts[j][l].vz += dx * vorts[i][k].s * magjtr
                    
                    vorts[i][k].vx += dz * vorts[j][l].s * magitr
                    vorts[i][k].vz -= dx * vorts[j][l].s * magitr
                end
            end
        end
    end
    return vorts
end

function mutual_ind(vorts::Vector{Vector{TwoDVort}}, nsurf :: Int, shed_ind::Vector{Vector{Int}})
    nlev = length(shed_ind)
    ntev = length(vorts) - nlev
    
    for i = 1:length(vorts)
        if i <= ntev
            for k = 1:nsurf
                #If it's a TEV, no need to check shed_ind
                for j = i+1:length(vorts)
                    if j <= ntev
                        for l = 1:nsurf
                            dx = vorts[i][k].x - vorts[j][l].x
                            dz = vorts[i][k].z - vorts[j][l].z
                            #source- tar
                            dsq = dx*dx + dz*dz
                            magitr = 1./(2*pi*sqrt(vorts[j][l].vc*vorts[j][l].vc*vorts[j][l].vc*vorts[j][l].vc + dsq*dsq))
                            magjtr = 1./(2*pi*sqrt(vorts[i][k].vc*vorts[i][k].vc*vorts[i][k].vc*vorts[i][k].vc + dsq*dsq))
                            
                            vorts[j][l].vx -= dz * vorts[i][k].s * magjtr
                            vorts[j][l].vz += dx * vorts[i][k].s * magjtr
                            
                            vorts[i][k].vx += dz * vorts[j][l].s * magitr
                            vorts[i][k].vz -= dx * vorts[j][l].s * magitr
                        end
                    else
                        for l in shed_ind[j - ntev]
                            dx = vorts[i][k].x - vorts[j][l].x
                            dz = vorts[i][k].z - vorts[j][l].z
                            #source- tar
                            dsq = dx*dx + dz*dz
                            magitr = 1./(2*pi*sqrt(vorts[j][l].vc*vorts[j][l].vc*vorts[j][l].vc*vorts[j][l].vc + dsq*dsq))
                            magjtr = 1./(2*pi*sqrt(vorts[i][k].vc*vorts[i][k].vc*vorts[i][k].vc*vorts[i][k].vc + dsq*dsq))
                            
                            vorts[j][l].vx -= dz * vorts[i][k].s * magjtr
                            vorts[j][l].vz += dx * vorts[i][k].s * magjtr
                            
                            vorts[i][k].vx += dz * vorts[j][l].s * magitr
                            vorts[i][k].vz -= dx * vorts[j][l].s * magitr
                        end
                    end
                end
            end
        else
            for k in shed_ind[i - ntev]
                for j = i+1:length(vorts)
                    for l in shed_ind[j - ntev]
                        dx = vorts[i][k].x - vorts[j][l].x
                        dz = vorts[i][k].z - vorts[j][l].z
                        #source- tar
                        dsq = dx*dx + dz*dz
                        magitr = 1./(2*pi*sqrt(vorts[j][l].vc*vorts[j][l].vc*vorts[j][l].vc*vorts[j][l].vc + dsq*dsq))
                        magjtr = 1./(2*pi*sqrt(vorts[i][k].vc*vorts[i][k].vc*vorts[i][k].vc*vorts[i][k].vc + dsq*dsq))
                        
                        vorts[j][l].vx -= dz * vorts[i][k].s * magjtr
                        vorts[j][l].vz += dx * vorts[i][k].s * magjtr
                        
                        vorts[i][k].vx += dz * vorts[j][l].s * magitr
                        vorts[i][k].vz -= dx * vorts[j][l].s * magitr
                    end
                end
            end
        end
    end
    return vorts
end

function mutual_ind_llt(vorts::Vector{Vector{TwoDVort}}, nspan :: Int, shed_ind::Vector{Vector{Int}})

    #Only stripwise influence

    nlev = length(shed_ind)
    ntev = length(vorts) - nlev

    for k = 1:nspan
        for i = 1:length(vorts)
            if i <= ntev
                for j = i+1:length(vorts)
                    if j <= ntev
                        dx = vorts[i][k].x - vorts[j][k].x
                        dz = vorts[i][k].z - vorts[j][k].z
                        #source- tar
                        dsq = dx*dx + dz*dz
                        magitr = 1./(2*pi*sqrt(vorts[j][k].vc*vorts[j][k].vc*vorts[j][k].vc*vorts[j][k].vc + dsq*dsq))
                        magjtr = 1./(2*pi*sqrt(vorts[i][k].vc*vorts[i][k].vc*vorts[i][k].vc*vorts[i][k].vc + dsq*dsq))
                        
                        vorts[j][k].vx -= dz * vorts[i][k].s * magjtr
                        vorts[j][k].vz += dx * vorts[i][k].s * magjtr
                        
                        vorts[i][k].vx += dz * vorts[j][k].s * magitr
                        vorts[i][k].vz -= dx * vorts[j][k].s * magitr
                    else
                        if k in shed_ind[j - ntev]
                            dx = vorts[i][k].x - vorts[j][k].x
                            dz = vorts[i][k].z - vorts[j][k].z
                            #source- tar
                            dsq = dx*dx + dz*dz
                            magitr = 1./(2*pi*sqrt(vorts[j][k].vc*vorts[j][k].vc*vorts[j][k].vc*vorts[j][k].vc + dsq*dsq))
                            magjtr = 1./(2*pi*sqrt(vorts[i][k].vc*vorts[i][k].vc*vorts[i][k].vc*vorts[i][k].vc + dsq*dsq))
                            
                            vorts[j][k].vx -= dz * vorts[i][k].s * magjtr
                            vorts[j][k].vz += dx * vorts[i][k].s * magjtr
                            
                            vorts[i][k].vx += dz * vorts[j][k].s * magitr
                            vorts[i][k].vz -= dx * vorts[j][k].s * magitr
                        end
                    end
                end
            else
                if k in shed_ind[i - ntev]
                    for j = i+1:length(vorts)
                        if k in shed_ind[j - ntev]
                            dx = vorts[i][k].x - vorts[j][k].x
                            dz = vorts[i][k].z - vorts[j][k].z
                            #source- tar
                            dsq = dx*dx + dz*dz
                            magitr = 1./(2*pi*sqrt(vorts[j][k].vc*vorts[j][k].vc*vorts[j][k].vc*vorts[j][k].vc + dsq*dsq))
                            magjtr = 1./(2*pi*sqrt(vorts[i][k].vc*vorts[i][k].vc*vorts[i][k].vc*vorts[i][k].vc + dsq*dsq))
                            
                            vorts[j][k].vx -= dz * vorts[i][k].s * magjtr
                            vorts[j][k].vz += dx * vorts[i][k].s * magjtr
                            
                            vorts[i][k].vx += dz * vorts[j][k].s * magitr
                            vorts[i][k].vz -= dx * vorts[j][k].s * magitr
                        end
                    end
                end
            end
        end
    end

    return vorts
end


function parmutual_ind_llt(vorts::Vector{Vector{TwoDVort}}, nspan :: Int, shed_ind::Vector{Vector{Int}})

    #Only stripwise influence

    nlev = length(shed_ind)
    ntev = length(vorts) - nlev

    @sync @everywhere for k = 1:nspan
        for i = 1:length(vorts)
            if i <= ntev
                for j = i+1:length(vorts)
                    if j <= ntev
                        dx = vorts[i][k].x - vorts[j][k].x
                        dz = vorts[i][k].z - vorts[j][k].z
                        #source- tar
                        dsq = dx*dx + dz*dz
                        magitr = 1./(2*pi*sqrt(vorts[j][k].vc*vorts[j][k].vc*vorts[j][k].vc*vorts[j][k].vc + dsq*dsq))
                        magjtr = 1./(2*pi*sqrt(vorts[i][k].vc*vorts[i][k].vc*vorts[i][k].vc*vorts[i][k].vc + dsq*dsq))
                        
                        vorts[j][k].vx -= dz * vorts[i][k].s * magjtr
                        vorts[j][k].vz += dx * vorts[i][k].s * magjtr
                        
                        vorts[i][k].vx += dz * vorts[j][k].s * magitr
                        vorts[i][k].vz -= dx * vorts[j][k].s * magitr
                    else
                        if k in shed_ind[j - ntev]
                            dx = vorts[i][k].x - vorts[j][k].x
                            dz = vorts[i][k].z - vorts[j][k].z
                            #source- tar
                            dsq = dx*dx + dz*dz
                            magitr = 1./(2*pi*sqrt(vorts[j][k].vc*vorts[j][k].vc*vorts[j][l].vc*vorts[j][l].vc + dsq*dsq))
                            magjtr = 1./(2*pi*sqrt(vorts[i][k].vc*vorts[i][k].vc*vorts[i][k].vc*vorts[i][k].vc + dsq*dsq))
                            
                            vorts[j][k].vx -= dz * vorts[i][k].s * magjtr
                            vorts[j][k].vz += dx * vorts[i][k].s * magjtr
                            
                            vorts[i][k].vx += dz * vorts[j][k].s * magitr
                            vorts[i][k].vz -= dx * vorts[j][k].s * magitr
                        end
                    end
                end
            else
                if k in shed_ind[i - ntev]
                    for j = i+1:length(vorts)
                        if k in shed_ind[j - ntev]
                            dx = vorts[i][k].x - vorts[j][k].x
                            dz = vorts[i][k].z - vorts[j][k].z
                            #source- tar
                            dsq = dx*dx + dz*dz
                            magitr = 1./(2*pi*sqrt(vorts[j][k].vc*vorts[j][k].vc*vorts[j][k].vc*vorts[j][k].vc + dsq*dsq))
                            magjtr = 1./(2*pi*sqrt(vorts[i][k].vc*vorts[i][k].vc*vorts[i][k].vc*vorts[i][k].vc + dsq*dsq))
                            
                            vorts[j][k].vx -= dz * vorts[i][k].s * magjtr
                            vorts[j][k].vz += dx * vorts[i][k].s * magjtr
                            
                            vorts[i][k].vx += dz * vorts[j][k].s * magitr
                            vorts[i][k].vz -= dx * vorts[j][k].s * magitr
                        end
                    end
                end
            end
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

function update_boundpos(surf::Vector{TwoDSurf}, dt::Float64)
    for n = 1:length(surf)
        for i = 1:surf[n].ndiv
            surf[n].bnd_x[i] = surf[n].bnd_x[i] + dt*((surf[n].pvt*surf[n].c - surf[n].x[i])*sin(surf[n].kinem.alpha)*surf[n].kinem.alphadot - surf[n].kinem.u + surf[n].cam[i]*cos(surf[n].kinem.alpha)*surf[n].kinem.alphadot)
            surf[n].bnd_z[i] = surf[n].bnd_z[i] + dt*(surf[n].kinem.hdot + (surf[n].pvt*surf[n].c - surf[n].x[i])*cos(surf[n].kinem.alpha)*surf[n].kinem.alphadot - surf[n].cam[i]*sin(surf[n].kinem.alpha)*surf[n].kinem.alphadot)
        end
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
# ---------------------------------------------------------------------------------------------

#Routines for 3D vortex ring method
# Vortex Line Calculations
function vortxl(x :: Float64, y :: Float64, z :: Float64, x1 :: Float64, y1 :: Float64, z1 :: Float64, x2 :: Float64, y2 :: Float64, z2 :: Float64, Gamma :: Float64)
    #VORTXL implements the vortex line induced velocity calculation
    Eps = 1e-10

    r1v = [x-x1; y-y1; z-z1]
    r2v = [x-x2; y-y2; z-z2]

    r1r2 = cross(r1v, r2v)
    modr1r2 = norm(r1r2)^2

    r1d = norm(r1v)
    r2d = norm(r2v)

    if r1d < Eps || r2d < Eps || modr1r2 < Eps
        uind = zeros(3,1)
    else
        r0v = [x2-x1; y2-y1; z2-z1]
        r0r1 = dot(r0v, r1v)
        r0r2 = dot(r0v, r2v)

        K = Gamma/(4*pi*modr1r2)*(r0r1/r1d - r0r2/r2d)
        
        uind = K.*r1r2
    end
    return uind
end


# Vortice Ring Calculations
function voring(x :: Float64, y :: Float64, z :: Float64, i :: Int, j :: Int, Gamma :: Float64, nspan :: Int, vr :: Vector{ThreeDSurfVRingGrid})
    #VORING implements the calculation of induced velocities of a vortex ring
    nij = (i - 1)*(nspan + 1) + j
    nip1j = i*(nspan + 1) + j
    nijp1 = (i - 1)*(nspan + 1) + j + 1
    nip1jp1 = i*(nspan + 1) + j + 1
    
    uind1 = vortxl(x,y,z,vr[nij].xv, vr[nij].yv, vr[nij].zv, vr[nijp1].xv, vr[nijp1].yv, vr[nijp1].zv, Gamma)
    uind2 = vortxl(x, y, z, vr[nijp1].xv, vr[nijp1].yv, vr[nijp1].zv, vr[nip1j].xv, vr[nip1j].yv, vr[nip1j].zv, Gamma)
    uind3 = vortxl(x, y, z, vr[nip1j].xv, vr[nip1j].yv, vr[nip1j].zv, vr[nip1jp1].xv, vr[nip1jp1].yv, vr[nip1jp1].zv, Gamma)
    uind4 = vortxl(x, y, z, vr[nip1jp1].xv, vr[nip1jp1].yv, vr[nip1jp1].zv,  vr[nij].xv, vr[nij].yv, vr[nij].zv, Gamma)

    uind = uind1 + uind2 + uind3 + uind4;
    
    uind_wake = uind2 + uind4
    return uind, uind_wake
end

function voring(x :: Float64, y :: Float64, z :: Float64, i :: Int, j :: Int, Gamma :: Float64, nspan :: Int, vr :: Vector{ThreeDWakeVRingGrid})
    #VORING implements the calculation of induced velocities of a vortex ring
    nij = (i - 1)*(nspan + 1) + j
    nip1j = i*(nspan + 1) + j
    nijp1 = (i - 1)*(nspan + 1) + j + 1
    nip1jp1 = i*(nspan + 1) + j + 1
    
    uind1 = vortxl(x,y,z,vr[nij].xv, vr[nij].yv, vr[nij].zv, vr[nijp1].xv, vr[nijp1].yv, vr[nijp1].zv, Gamma)
    uind2 = vortxl(x, y, z, vr[nijp1].xv, vr[nijp1].yv, vr[nijp1].zv, vr[nip1j].xv, vr[nip1j].yv, vr[nip1j].zv, Gamma)
    uind3 = vortxl(x, y, z, vr[nip1j].xv, vr[nip1j].yv, vr[nip1j].zv, vr[nip1jp1].xv, vr[nip1jp1].yv, vr[nip1jp1].zv, Gamma)
    uind4 = vortxl(x, y, z, vr[nip1jp1].xv, vr[nip1jp1].yv, vr[nip1jp1].zv,  vr[nij].xv, vr[nij].yv, vr[nij].zv, Gamma)

    uind = uind1 + uind2 + uind3 + uind4;
    
    uind_wake = uind2 + uind4
    return uind, uind_wake
end

function voring_I(x :: Float64, y :: Float64, z :: Float64, i :: Int, j :: Int, Gamma :: Float64, nspan :: Int, vr :: Vector{ThreeDSurfVRingGrid})
    #VORING implements the calculation of induced velocities of a vortex ring
    nij = (i - 1)*(nspan + 1) + j
    nip1j = i*(nspan + 1) + j
    nijp1 = (i - 1)*(nspan + 1) + j + 1
    nip1jp1 = i*(nspan + 1) + j + 1
    
    uind1 = vortxl(x, y, z,vr[nij].xv_I, vr[nij].yv_I, vr[nij].zv_I, vr[nijp1].xv_I, vr[nijp1].yv_I, vr[nijp1].zv_I, Gamma)
    uind2 = vortxl(x, y, z, vr[nijp1].xv_I, vr[nijp1].yv_I, vr[nijp1].zv_I, vr[nip1j].xv_I, vr[nip1j].yv_I, vr[nip1j].zv_I, Gamma)
    uind3 = vortxl(x, y, z, vr[nip1j].xv_I, vr[nip1j].yv_I, vr[nip1j].zv_I, vr[nip1jp1].xv_I, vr[nip1jp1].yv_I, vr[nip1jp1].zv_I, Gamma)
    uind4 = vortxl(x, y, z, vr[nip1jp1].xv_I, vr[nip1jp1].yv_I, vr[nip1jp1].zv_I,  vr[nij].xv_I, vr[nij].yv_I, vr[nij].zv_I, Gamma)
    
    uind = uind1 + uind2 + uind3 + uind4;
    
    uind_wake = uind2 + uind4
    return uind, uind_wake
end
function voring_I(x :: Float64, y :: Float64, z :: Float64, i :: Int, j :: Int, Gamma :: Float64, nspan :: Int, vr :: Vector{ThreeDWakeVRingGrid})
    #VORING implements the calculation of induced velocities of a vortex ring
    nij = (i - 1)*(nspan + 1) + j
    nip1j = i*(nspan + 1) + j
    nijp1 = (i - 1)*(nspan + 1) + j + 1
    nip1jp1 = i*(nspan + 1) + j + 1
    
    uind1 = vortxl(x, y, z,vr[nij].xv_I, vr[nij].yv_I, vr[nij].zv_I, vr[nijp1].xv_I, vr[nijp1].yv_I, vr[nijp1].zv_I, Gamma)
    uind2 = vortxl(x, y, z, vr[nijp1].xv_I, vr[nijp1].yv_I, vr[nijp1].zv_I, vr[nip1j].xv_I, vr[nip1j].yv_I, vr[nip1j].zv_I, Gamma)
    uind3 = vortxl(x, y, z, vr[nip1j].xv_I, vr[nip1j].yv_I, vr[nip1j].zv_I, vr[nip1jp1].xv_I, vr[nip1jp1].yv_I, vr[nip1jp1].zv_I, Gamma)
    uind4 = vortxl(x, y, z, vr[nip1jp1].xv_I, vr[nip1jp1].yv_I, vr[nip1jp1].zv_I,  vr[nij].xv_I, vr[nij].yv_I, vr[nij].zv_I, Gamma)
    
    uind = uind1 + uind2 + uind3 + uind4;
    
    uind_wake = uind2 + uind4
    return uind, uind_wake
end

# ---------------------------------------------------------------------------------------------
# function calc_a03d(surf::ThreeDSurfSimple)
    
#     lhs = zeros(surf.nspan, surf.nspan)
#     rhs = zeros(surf.nspan)
    
#     for i = 1:surf.nspan
#         for n = 1:surf.nspan
#             nn = 2*n - 1
#             lhs[i,n] = sin(nn*surf.psi[i])*(sin(surf.psi[i]) + (nn*pi/(2*surf.AR)))
#         end
#         rhs[i] = pi*sin(surf.psi[i])*surf.bc[i]/(2*surf.AR)
#     end
    
#     surf.bcoeff[:] = lhs \ rhs
    
#     for i = 1:surf.nspan
#         surf.a03d[i] = 0
#         for n = 1:surf.nspan
#             nn = 2*n - 1
#             surf.a03d[i] = surf.a03d[i] - real(nn)*surf.bcoeff[n]*sin(nn*surf.psi[i])/sin(surf.psi[i])
#         end
#     end
# end
function calc_a03d(surf::ThreeDSurfSimple, shedv :: Array{Int} = Int[], levstr :: Array{Float64} = Float64[])
    
    lhs = zeros(surf.nspan, surf.nspan)
    rhs = zeros(surf.nspan)
    
    for i = 1:surf.nspan
        for n = 1:surf.nspan
            nn = 2*n - 1
            lhs[i,n] = sin(nn*surf.psi[i])*(sin(surf.psi[i]) + (nn*pi/(2*surf.AR)))
        end
        if i in shedv
            rhs[i] = pi*sin(surf.psi[i])*(surf.bc[i] + levstr[i]/pi)/(2*surf.AR)
        else
        rhs[i] = pi*sin(surf.psi[i])*surf.bc[i]/(2*surf.AR)
        end
    end
    
    surf.bcoeff[:] = lhs \ rhs
    
    for i = 1:surf.nspan
        surf.a03d[i] = 0
        for n = 1:surf.nspan
            nn = 2*n - 1
            surf.a03d[i] = surf.a03d[i] - real(nn)*surf.bcoeff[n]*sin(nn*surf.psi[i])/sin(surf.psi[i])
        end
    end
end

# function calc_a03dspl(surf::ThreeDSurfSimple)

#     nterm = 70
#     bcspl = Spline1D(surf.psi, surf.bc)
    
#     lhs = zeros(nterm, nterm)
#     rhs = zeros(nterm)
    
#     a03d = zeros(nterm)
#     psi = zeros(nterm)
#     bcoeff = zeros(nterm)
    
#     for i = 1:nterm
#         psi[i] = real(i)*(pi/2)/nterm
#         bc = Dierckx.evaluate(bcspl, psi[i])
#         for n = 1:nterm
#             nn = 2*n - 1
#             lhs[i,n] = sin(nn*psi[i])*(sin(psi[i]) + (nn*pi/(2*surf.AR)))
#         end
#         rhs[i] = pi*sin(psi[i])*bc/(2*surf.AR)
#     end

#     bcoeff[:] = lhs \ rhs
   
#     for i = 1:nterm
#         a03d[i] = 0
#         for n = 1:nterm
#             nn = 2*n - 1
#             a03d[i] = a03d[i] - real(nn)*bcoeff[n]*sin(nn*psi[i])/sin(psi[i])
#         end
#     end

#     a03dspl = Spline1D(psi, a03d)
#     for i = 1:surf.nspan
#         surf.a03d[i] = evaluate(a03dspl, surf.psi[i])
#     end
    
# end

function calc_a03dspl(surf::ThreeDSurfSimple, shedv :: Array{Int} = Int[], levstr :: Array{Float64} = Float64[])

    nterm = 70
    bcspl = Spline1D(surf.psi, surf.bc)
    
    lhs = zeros(nterm, nterm)
    rhs = zeros(nterm)
    
    a03d = zeros(nterm)
    psi = zeros(nterm)
    bcoeff = zeros(nterm)
    
    for i = 1:nterm
        psi[i] = real(i)*(pi/2)/nterm
        bc = Dierckx.evaluate(bcspl, psi[i])
        for n = 1:nterm
            nn = 2*n - 1
            lhs[i,n] = sin(nn*psi[i])*(sin(psi[i]) + (nn*pi/(2*surf.AR)))
        end
        rhs[i] = pi*sin(psi[i])*bc/(2*surf.AR)
    end

    bcoeff[:] = lhs \ rhs
   
    for i = 1:nterm
        a03d[i] = 0
        for n = 1:nterm
            nn = 2*n - 1
            a03d[i] = a03d[i] - real(nn)*bcoeff[n]*sin(nn*psi[i])/sin(psi[i])
        end
    end

    a03dspl = Spline1D(psi, a03d)
    for i = 1:surf.nspan
        surf.a03d[i] = evaluate(a03dspl, surf.psi[i])
    end
    
end


# ---------------------------------------------------------------------------------------------
function bendfirstmode(fn :: Float64, fd :: Float64, td :: Float64, c :: Float64, s :: Float64,  y :: Float64)
    
    wn = fn*2*pi
    td = td/c
    wd = 2*pi*fd
    # a = sqrt(EI/m)
    
    num1 = 0.597^2*pi^2/s^2
    a = wn/num1
    
    wna = sqrt(wn/a)
    
    
    D = td/((sin(wna*s) - sinh(wn*s))/(cosh(wn*s) + cos(wn*s))*(sinh(wna*s) - sin(wna*s)) + (cosh(wna*s) - cos(wna*s)))
    
    h_z =  D*((sin(wna*s) - sinh(wn*s))/(cosh(wn*s) + cos(wn*s))*(sinh(wna*y) - sin(wna*y)) + (cosh(wna*y) - cos(wna*y)))
    return h_z
end

# ---------------------------------------------------------------------------------------------

function calc_q1(x1 :: Vector{Float64}, x2 :: Vector{Float64}, x :: Vector{Float64})# !In the form q1(s,e,r)
    a1 = x1[1] - x[1]
    a2 = x1[2] - x[2]
    a3 = x1[3] - x[3]
    b1 = x2[1] - x[1]
    b2 = x2[2] - x[2]
    b3 = x2[3] - x[3]
    dt = dot([a1; a2; a3], [b1; b2; b3])
    cr = cross([a1; a2; a3], [b1; b2; b3])
    t1 = dot(cr,cr)
    moda = sqrt(a1^2 + a2^2 + a3^2)
    modb = sqrt(b1^2 + b2^2 + b3^2)
    t2 = moda + modb
    t3 = 1. - dt/(moda*modb)
    vx = cr[1]*t2*t3/t1
    vy = cr[2]*t2*t3/t1
    vz = cr[3]*t2*t3/t1
    return [vx; vy; vz] 
end 


function calc_q2(x1 :: Vector{Float64}, x2 :: Vector{Float64}, x :: Vector{Float64}) #In the form q2(s,t,r,out)

    a1 = x1[1] - x[1]
    a2 = x1[2] - x[2]
    a3 = x1[3] - x[3]
    moda = sqrt(a1^2 + a2^2 + a3^2)
    cr = cross([a1; a2; a3], x2)
    dotacrt = dot(cr,cr)
    dotat = dot([a1; a2; a3], x2)

    vx = cr[1]*(1.-(dotat/moda))/dotacrt
    vy = cr[2]*(1.-(dotat/moda))/dotacrt
    vz = cr[3]*(1.-(dotat/moda))/dotacrt
    return [vx; vy; vz] 
end 


# ---------------------------------------------------------------------------------------------

function update_IC(surf :: ThreeDSurfWeiss, hs_prev :: Vector{Float64}, he_prev :: Vector{Float64})
    for i = 1:surf.nlat
        if i == 1
            h_s = surf.s2d[i].kinem.h - 0.5*(surf.s2d[i+1].kinem.h - surf.s2d[i].kinem.h)
        else
            h_s = surf.s2d[i].kinem.h - 0.5*(surf.s2d[i].kinem.h - surf.s2d[i-1].kinem.h)
        end
        if i == nlat
            h_e = surf.s2d[i].kinem.h + 0.5*(surf.s2d[i].kinem.h - surf.s2d[i-1].kinem.h)
        else
            h_e = surf.s2d[i].kinem.h + 0.5*(surf.s2d[i+1].kinem.h - surf.s2d[i].kinem.h)
        end
        
        hdiff_s = hs_prev[i] - h_s
        hdiff_e = he_prev[i] - h_e
        
        surf.s1[i,3] -= hdiff_s
        surf.e1[i,3] -= hdiff_e
    end

    for ilat = 1:nlat
        surf.m[ilat,3] = 0.5*(surf.s1[ilat,3] + surf.e1[ilat,3])
        surf.s[ilat,3] = surf.s1[ilat,3]
        surf.e[ilat,3] = surf.e1[ilat,3]
        surf.cx[ilat,3] = surf.m[ilat,3]
        surf.c0[ilat,3] = surf.cx[ilat,3]
        surf.s0[ilat,3] = surf.s[ilat,3] 
        surf.e0[ilat,3] = surf.e[ilat,3]
        tan1 = surf.e[ilat,1] - surf.s[ilat,1]
        tan2 = surf.e[ilat,2] - surf.s[ilat,2]  
        tan3 = surf.e[ilat,3] - surf.s[ilat,3]
        tanmod=sqrt(tan1*tan1+tan2*tan2+tan3*tan3)
        tan1 = tan1/tanmod
        tan2 = tan2/tanmod
        tan3 = tan3/tanmod
        surf.norm[ilat,:] = cross([tan1; tan2;tan3],[1.; 0.; 0.])
        surf.dih[ilat] = atan((surf.s[ilat,3] - surf.e[ilat,3])/(surf.e[ilat,2] - surf.s[ilat,2]))
        surf.ds[ilat] = sqrt((surf.e[ilat,2] - surf.s[ilat,2])*(surf.e[ilat,2] - surf.s[ilat,2]) + (surf.e[ilat,3] - surf.s[ilat,3])*(surf.e[ilat,3] - surf.s[ilat,3]))
    end
    
    a = zeros(3)
    b = zeros(3)
    
    for i=1:surf.nlat
        for j=1:surf.nlat 
            q1terms = calc_q2(surf.s0[j,:],[-1.; 0.; 0.], surf.c0[i,:])
            q2terms = calc_q2(surf.e0[j,:],[-1.; 0.; 0.], surf.c0[i,:])
            q1t = dot(q1terms, surf.norm[i,:])
            q2t = dot(q2terms, surf.norm[i,:])
            surf.ICt[i,j] = (1/(4*pi))*(2*q2t - 2*q1t) 
            
            a[1] = surf.s[j,1] - surf.cx[i,1]
            a[2] = surf.s[j,2] - surf.cx[i,2]
            a[3] = surf.s[j,3] - surf.cx[i,3]
            b[1] = surf.e[j,1] - surf.cx[i,1]
            b[2] = surf.e[j,2] - surf.cx[i,2]
            b[3] = surf.e[j,3] - surf.cx[i,3]
            
            cr = cross(a,b)
            check = dot(cr, cr)
            
            if (check < errtiny) then
                    q1terms[:] = 0.
            else
                q1terms = calc_q1(surf.s[j,:],surf.e[j,:],surf.cx[i,:])
            end
            
            q2terms = calc_q2(surf.s[j,:], [-1.; 0.; 0.], surf.cx[i,:])
            q3terms = calc_q2(surf.e[j,:], [-1.; 0.; 0.], surf.cx[i,:])
            
            q1t =  dot(q1terms, surf.norm[i,:])
            q2t =  dot(q2terms, surf.norm[i,:])
            q3t = dot(q3terms, surf.norm[i,:])
            surf.IC[i,j] = (1/(4*pi))*(q1t - q2t + q3t)
        end
    end
end


	"""	
	```julia-repl
	Calculate separation point position f basing on static CN(alpha) data.
	To study this function step by step, refer to the notebook "1.1. Steady state separation point model". 
	
	#Arguments
	```
	- 'alpha::Vector{Float64}': vector of angles of attack
	- 'CN::Vector{Float64}': vector of normal force coefficients CN corresponding to alpha
	```
 	
	#Keyword arguments
	```
	- 'lin::Float64': max. alpha at which CN is still linear, in degrees (default lin = 7)
	```
	"""
function fFromCN(alpha::Vector{Float64}, CN::Vector{Float64}; lin::Float64 = 7.)

	#trim domain to linear part
    idx_linear = 1
    while(alpha[idx_linear]<lin)
        idx_linear+=1
    end
    
    #find b,a for y = bx + a
    a,b = linreg(alpha[1:idx_linear], CN[1:idx_linear])
    alpha0 = -a/b

	#Calculate max. slope from static data (amax) to guarantee that f will not exceed 1
    amax = 0
    for idx = 1 : idx_linear
        if(CN[idx]/(alpha[idx]-alpha0) > amax)
        amax = CN[idx]/(alpha[idx]-alpha0)
        end
    end

	#Calculate f 	
	f = zeros(length(alpha))

	for idx = 1 : length(alpha)
		f[idx] = (2*sqrt(CN[idx]/(amax*(alpha[idx]-alpha0)))-1)^2
	end 

	f
end



	"""	
	```julia-repl
	Model separation point position f. To study this function step by step, refer to the notebook "1.1 Steady state separation point model". 
	
	#Arguments
	```
	- 'alpha::Vector{Float64}': vector of angles of attack
	- 'f::Vector{Float64}': vector of separation point f positions corresponding to alpha
	- 'model::String': Sheng or Original models available.
	```
	"""
function findStaticCoeff(alpha::Vector{Float64}, f::Vector{Float64}, model::String)

	if(model == "Sheng")
		fsep =  0.6
		temp = abs(f.-fsep)
		idxf = findfirst(temp,minimum(temp))

		alpha1 = alpha[idxf]
		
		#STEP 3: OPTIMIZE 
		pcws1(x, s1) = 1.0.-0.4.*exp((abs(x).-alpha1)./s1)
		pcws2(x, s2) = 0.02.+0.58.*exp((alpha1.-abs(x))./s2)
		s1_0 = [3.5]
		s2_0 = [3.5]

		alpha_pcws1 = alpha[1:idxf]
		alpha_pcws2 = alpha[idxf:end]

		f_pcws1 = f[1:idxf]
		f_pcws2 = f[idxf:end]

		fit1 = curve_fit(pcws1, alpha_pcws1, f_pcws1, s1_0)
		fit2 = curve_fit(pcws2, alpha_pcws2, f_pcws2, s2_0)
		s1 = fit1.param[1]
		s2 = fit2.param[1]
		f=vcat(pcws1(alpha_pcws1,fit1.param), pcws2(alpha_pcws2,fit1.param)[2:end])
		
	elseif(model =="Original")
		fsep =  0.7
		temp = abs(f.-fsep)
		idxf = findfirst(temp,minimum(temp))

		alpha1 = alpha[idxf]

		#STEP 3: OPTIMIZE 
		pcws3(x, s1) = 1.0.-0.3.*exp((abs(x).-alpha1)./s1)
		pcws4(x, s2) = 0.04.+0.66.*exp((alpha1.-abs(x))./s2)	
		s1_0 = [3.5]
		s2_0 = [3.5]

		
		alpha_pcws3 = alpha[1:idxf]
		alpha_pcws4 = alpha[idxf:end]

		f_pcws3 = f[1:idxf]
		f_pcws4 = f[idxf:end]

		fit3 = curve_fit(pcws3, alpha_pcws3, f_pcws3, s1_0)
		fit4 = curve_fit(pcws4, alpha_pcws4, f_pcws4, s2_0)
		s1 = fit3.param[1]
		s2 = fit4.param[1]
		f=vcat(pcws3(alpha_pcws3,fit3.param), pcws4(alpha_pcws4,fit4.param)[2:end])
	else
		println("Wrong f model specified. Choose Sheng or Original.")
	end
		
	alpha1, s1, s2, f
end	

"""	
	```julia-repl
	Calculate f(a) characteristic for a set of static constants.
	
	Output: f(a)
	
	#Arguments
	
	```
	
	- 'alpha::Vector{Float64}': vector of angles of attack in deg
	- 'alpha1::Float64': alpha in deg for which f = 0.6 (Sheng) or f = 0.7 (Original), can be approximated as stall angle
	- 's1::Float64': static constant in deg for piecewise model, attached regime
	- 's2::Float64': static constant in deg for piecewise model, separated regime	
	- 'model::String': Two models are available: "Sheng" and "Original"
	
		
	```
	
	#Keyword argument
	
	```
	
	- 'par': Default is 0.1. Modify to larger values if no result is given by function.

	"""
function fFromConst(alpha::Vector{Float64}, alpha1::Float64, s1::Float64, s2::Float64, model::String; par=0.1)
	
	if(model == "Sheng")
		c1 = 1.0
		c2 = 0.4
		c3 = 0.02
		c4 = 0.58
	elseif(model=="Original")
		c1 = 1.0
		c2 = 0.3
		c3 = 0.04
		c4 = 0.66
	else
		println("Wrong f model specified. Choose Sheng or Original.")
	end
		
		pcws_1(x, s1, alpha1) = c1.-c2.*exp((abs(x).-alpha1)./s1)
		pcws_2(x, s2, alpha1) = c3.+c4.*exp((alpha1.-abs(x))./s2)	
	
		idxf = 1
		while(abs(alpha[idxf]-alpha1)>par)
			idxf+=1
		end
		
		f = vcat(pcws_1(alpha[1:idxf],s1,alpha1), pcws_2(alpha[idxf+1:end],s2,alpha1))
end

	"""	
	```julia-repl
	Obtain static constants for moment coefficient calculation, basing on CM(a), CN(a) and f(a). 
	
	Output: k0, k1, k2, m
	
	#Arguments
	
	```
	
	- 'CM::Vector{Float64}': vector of moment coefficients CM for a set of angles of attack
	- 'CN::Vector{Float64}': vector of normal force coefficients CN for the same set of angles of attack
	- 'f::Vector{Float64}': vector of separation point position f as a fraction of c, corresponding to the above
	
	```
	
	#Keyword argument
	
	```
	
	- 'par': Upper limit of f that is used in calculations. Default is 0.9, may be modified for a better fit.
	
	```
 	
	For more information refer to "2.1 Static constants for CM" notebook.
	```
	"""
function cmstatic(CM::Vector{Float64}, CN::Vector{Float64}, f::Vector{Float64}; par = 0.9)

fun(f, p) = p[1] + p[2]*(1.-f)+p[3]*sin(pi.*f.^p[4]) #1 - k0, 2 - k1, 3 - k2, 4 - m
init = [0., -0.135, 0.04, 2.] #initial values from [1] Z. Liu, J. C. S. Lai, J. Young, and F.-B. Tian, 
                            #“Discrete Vortex Method with Flow Separation Corrections for Flapping-Foil Power Generators,” 
                            #AIAA Journal, vol. 55, no. 2, pp. 410–418, Feb. 2017.
idx = 1
while(f[idx]>par)
    idx+=1
end
							
fit = curve_fit(fun, f[idx:end], CM[idx:end]./CN[idx:end], init)

fit.param
end

# -------------------------------------------------------------------------------------------