function calc_a03d(surf::ThreeDSurfSimple)

    lhs = zeros(surf.nspan, surf.nspan)
    rhs = zeros(surf.nspan)

    for i = 1:surf.nspan
        for n = 1:surf.nspan
            nn = 2*n - 1
            lhs[i,n] = sin(nn*surf.psi[i])*(sin(surf.psi[i]) + (nn*pi/(2*surf.AR)))
        end
        rhs[i] = pi*sin(surf.psi[i])*surf.bc[i]/(2*surf.AR)
    end

    surf.bcoeff[:] = lhs \ rhs

    for i = 1:surf.nspan
        surf.a03d[i] = 0
   
        for n = 1:surf.nspan
            nn = 2*n - 1
            surf.a03d[i] -= real(nn)*surf.bcoeff[n]*sin(nn*surf.psi[i])/sin(surf.psi[i])
        end
    end
    return surf
end

function calc_a0a13d(surf::ThreeDSurfSimple)
        
    lhs = zeros(surf.nspan, surf.nspan)
    rhs = zeros(surf.nspan)
    
    for i = 1:surf.nspan
        integ0 = simpleTrapz(surf.s2d[i].cam_slope, surf.s2d[i].theta)
        integ1 = simpleTrapz(surf.s2d[i].cam_slope.*cos.(surf.s2d[i].theta), surf.s2d[i].theta)
        for n = 1:surf.nspan
            nn = 2*n - 1
            lhs[i,n] = sin(nn*surf.psi[i])*(sin(surf.psi[i]) + nn*pi/(2*surf.AR)*
                                            (cos(surf.s2d[i].kinem.alpha) + sin(surf.s2d[i].kinem.alpha)*(integ0 - integ1)/pi)) 
        end
        rhs[i] = pi*sin(surf.psi[i])*surf.bc[i]/(2*surf.AR)
    end

    surf.bcoeff[:] = lhs \ rhs

    for i = 1:surf.nspan
        integ0 = simpleTrapz(surf.s2d[i].cam_slope, surf.s2d[i].theta)
        integ1 = simpleTrapz(surf.s2d[i].cam_slope.*cos.(surf.s2d[i].theta), surf.s2d[i].theta)
        surf.a03d[i] = 0
        surf.aterm3d[1,i] = 0

        for n = 1:surf.nspan
            nn = 2*n - 1
            surf.a03d[i] -= real(nn)*surf.bcoeff[n]*sin(nn*surf.psi[i])/sin(surf.psi[i])*(cos(surf.s2d[i].kinem.alpha) + sin(surf.s2d[i].kinem.alpha)*integ0/pi)
            surf.aterm3d[1,i] += 2*real(nn)*surf.bcoeff[n]*sin(nn*surf.psi[i])/(sin(surf.psi[i])*pi)*sin(surf.s2d[i].kinem.alpha)*integ1
        end
    end
    return surf
end

function calc_a0a13d_wlev(surf::ThreeDSurfSimple)
        
    lhs = zeros(surf.nspan, surf.nspan)
    rhs = zeros(surf.nspan)

    levcount = 0
    for i = 1:surf.nspan
        integ0 = simpleTrapz(surf.s2d[i].cam_slope, surf.s2d[i].theta)
        integ1 = simpleTrapz(surf.s2d[i].cam_slope.*cos.(surf.s2d[i].theta), surf.s2d[i].theta)
        for n = 1:surf.nspan
            nn = 2*n - 1
            lhs[i,n] = sin(nn*surf.psi[i])*(sin(surf.psi[i]) + nn*pi/(2*surf.AR)*
                                            (cos(surf.s2d[i].kinem.alpha) + sin(surf.s2d[i].kinem.alpha)*(integ0 - integ1)/pi)) 
        end
        rhs[i] = (surf.levstr[i] + pi*surf.bc[i])*sin(surf.psi[i])/(2*surf.AR)
    end

    surf.bcoeff[:] = lhs \ rhs

    for i = 1:surf.nspan
        integ0 = simpleTrapz(surf.s2d[i].cam_slope, surf.s2d[i].theta)
        integ1 = simpleTrapz(surf.s2d[i].cam_slope.*cos.(surf.s2d[i].theta), surf.s2d[i].theta)
        surf.a03d[i] = 0
        surf.aterm3d[1,i] = 0

        for n = 1:surf.nspan
            nn = 2*n - 1
            surf.a03d[i] -= real(nn)*surf.bcoeff[n]*sin(nn*surf.psi[i])/sin(surf.psi[i])*(cos(surf.s2d[i].kinem.alpha) + sin(surf.s2d[i].kinem.alpha)*integ0/pi)
            surf.aterm3d[1,i] += 2*real(nn)*surf.bcoeff[n]*sin(nn*surf.psi[i])/(sin(surf.psi[i])*pi)*sin(surf.s2d[i].kinem.alpha)*integ1
        end
    end
    return surf
end

function calc_a0a13d_wlev2(surf::ThreeDSurfSimple, levstr)
        
    lhs = zeros(surf.nspan, surf.nspan)
    rhs = zeros(surf.nspan)

    for i = 1:surf.nspan
        integ0 = simpleTrapz(surf.s2d[i].cam_slope, surf.s2d[i].theta)
        integ1 = simpleTrapz(surf.s2d[i].cam_slope.*cos.(surf.s2d[i].theta), surf.s2d[i].theta)
        for n = 1:surf.nspan
            nn = 2*n - 1
            lhs[i,n] = sin(nn*surf.psi[i])*(sin(surf.psi[i]) + nn*pi/(2*surf.AR)*
                                            (cos(surf.s2d[i].kinem.alpha) + sin(surf.s2d[i].kinem.alpha)*(integ0 - integ1)/pi)) 
        end
        
        rhs[i] = (levstr[i] + pi*surf.bc[i])*sin(surf.psi[i])/(2*surf. AR)
    end
    
    surf.bcoeff[:] = lhs \ rhs

    for i = 1:surf.nspan
        integ0 = simpleTrapz(surf.s2d[i].cam_slope, surf.s2d[i].theta)
        integ1 = simpleTrapz(surf.s2d[i].cam_slope.*cos.(surf.s2d[i].theta), surf.s2d[i].theta)
        surf.a03d[i] = 0
        surf.aterm3d[1,i] = 0

        for n = 1:surf.nspan
            nn = 2*n - 1
            surf.a03d[i] -= real(nn)*surf.bcoeff[n]*sin(nn*surf.psi[i])/sin(surf.psi[i])*(cos(surf.s2d[i].kinem.alpha) + sin(surf.s2d[i].kinem.alpha)*integ0/pi)
            surf.aterm3d[1,i] += 2*real(nn)*surf.bcoeff[n]*sin(nn*surf.psi[i])/(sin(surf.psi[i])*pi)*sin(surf.s2d[i].kinem.alpha)*integ1
        end
    end
    return surf
end


function calc_a2toan3d(surf::ThreeDSurfSimple)
    for ia = 2:surf.naterm
        for i = 1:surf.nspan
            surf.aterm3d[ia,i] = 0
            integ = simpleTrapz(surf.s2d[i].cam_slope.*cos.(ia*surf.s2d[i].theta), surf.s2d[i].theta)
            
            for n = 1:surf.nspan
                nn = 2*n - 1
                surf.aterm3d[ia,i] += 2*real(nn)*surf.bcoeff[n]*sin(nn*surf.psi[i])/(sin(surf.psi[i])*pi)*sin(surf.s2d[i].kinem.alpha)*integ
            end
        end
    end
end

# Function for calculating the wake rollup
function wakeroll(surf3d::ThreeDSurfSimple, field3d::ThreeDFieldSimple, dt)
    
    for ispan = 1:surf3d.nspan	 				  
	surf = surf3d.s2d[ispan]					    
	curfield = field3d.f2d[ispan]    

	wi = 0
        for n = 1:surf3d.nspan
            nn = 2*n -2
            wi += real(nn)*surf3d.bcoeff[n]*sin(nn*surf3d.psi[ispan])/sin(surf3d.psi[ispan])
        end


        nlev = length(curfield.lev)
        ntev = length(curfield.tev)
        nextv = length(curfield.extv)
        
        #Clean induced velocities
        for i = 1:ntev
            curfield.tev[i].vx = 0
            curfield.tev[i].vz = -wi
        end
        
        for i = 1:nlev
            curfield.lev[i].vx = 0
            curfield.lev[i].vz = -wi
        end
        
        for i = 1:nextv
            curfield.extv[i].vx = 0
            curfield.extv[i].vz = -wi
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

            #Add the influence of freestream velocities
        for i = 1:ntev
            curfield.tev[i].vx += curfield.u[1]
            curfield.tev[i].vz += curfield.w[1]
        end
        for i = 1:nlev
        curfield.lev[i].vx += curfield.u[1]
            curfield.lev[i].vz += curfield.w[1]
        end
        for i = 1:nextv
            curfield.extv[i].vx += curfield.u[1]
            curfield.extv[i].vz += curfield.w[1]
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
    end
    
    return field3d
end

function update_kinemCantilever(surf::ThreeDSurfSimple, dt::Float64, kinem::KinemCantilever, Phi::Array{Float64}, m_eta_s_etadd::Array{Float64}, m_eta_s_etad::Array{Float64}, m_eta_s_eta::Array{Float64}, cl::Array{Float64}, cm::Array{Float64})

    Qvecloc = zeros(4, surf.nspan)
    
    rho = 1.1
    q = 1. *pi/180
    
    pdyn = 0.5*rho*surf.uref^2
    
    #Update previous terms
    kinem.eta_pr = kinem.eta
    kinem.etad_pr3 = kinem.etad_pr2
    kinem.etad_pr2 = kinem.etad_pr
    kinem.etad_pr = kinem.etad
    kinem.etadd_pr3 = kinem.etadd_pr2
    kinem.etadd_pr2 = kinem.etadd_pr

    for i = 1:surf.nspan
        if i == 1
            yl = -0.5*surf.AR*surf.cref
            yu = 0.5*(surf.yle[i] + surf.yle[i+1])
        elseif i == surf.nspan
            yl = 0.5*(surf.yle[i] + surf.yle[i-1])
            yu = 0. 
        else
            yl = 0.5*(surf.yle[i] + surf.yle[i-1])
            yu = 0.5*(surf.yle[i] + surf.yle[i+1])
        end
        dy = yu - yl
        Qvecloc[:,surf.nspan-i+1] = Phi[:,:,surf.nspan-i+1]'*pdyn*[-cl[i]*dy*surf.s2d[i].c; cm[i]*dy*surf.s2d[i].c^2; 0];
    end
    
    Qvec = sum(Qvecloc,dims=2)
    
    kinem.etadd_pr = m_eta_s_etadd\(-m_eta_s_etad*kinem.etad - m_eta_s_eta*kinem.eta + Qvec)

    kinem.etad = kinem.etad_pr + (dt/12.)*(23*kinem.etadd_pr - 16*kinem.etadd_pr2 + 5*kinem.etadd_pr3)
    
    kinem.eta = kinem.eta_pr + (dt/12)*(23*kinem.etad_pr - 16*kinem.etad_pr2 + 5*kinem.etad_pr3)
    
    for i = 1:surf.nspan
        vecphys = Phi[:,:,i]*kinem.eta
        vecphys_dot = Phi[:,:,i]*kinem.etad
        j = surf.nspan - i + 1
        surf.s2d[j].kinem.h = -vecphys[1]
        surf.s2d[j].kinem.alpha = vecphys[2] + q
        surf.s2d[j].kinem.hdot = -vecphys_dot[1]
        surf.s2d[j].kinem.alphadot = vecphys_dot[2]
    end
        
end
