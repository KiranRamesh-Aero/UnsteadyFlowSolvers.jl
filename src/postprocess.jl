# ---------------------------------------------------------------------------------------------
# Function for viewing the 2D vorticity field around the airfoil
function view_vorts(surf::TwoDSurf, field::TwoDFlowField)
    scatter(map(q->q.x, field.tev),map(q->q.z,field.tev))
    scatter(map(q->q.x, field.lev),map(q->q.z,field.lev))
    plot(map(q->q.x, surf.bv),map(q->q.z,surf.bv))
end

function view_vorts(surf::TwoDSurf_2DOF, field::TwoDFlowField)
    scatter(map(q->q.x, field.tev),map(q->q.z,field.tev))
    scatter(map(q->q.x, field.lev),map(q->q.z,field.lev))
    plot(map(q->q.x, surf.bv),map(q->q.z,surf.bv))
end

function view_vorts(surf::TwoDSurfwFlap, field::TwoDFlowField)
    scatter(map(q->q.x, field.tev),map(q->q.z,field.tev))
    scatter(map(q->q.x, field.lev),map(q->q.z,field.lev))
    plot(map(q->q.x, surf.bv),map(q->q.z,surf.bv))
end

# ---------------------------------------------------------------------------------------------
# Function that calculates the aerodynamic forces
function calc_forces(surf::TwoDSurf)

    # First term in eqn (2.30) Ramesh et al. in coefficient form
    cnc = 2*pi*(surf.kinem.u*cos(surf.kinem.alpha)/surf.uref + surf.kinem.hdot*sin(surf.kinem.alpha)/surf.uref)*(surf.a0[1] + surf.aterm[1]/2.)

    # Second term in eqn (2.30) Ramesh et al. in coefficient form
    cnnc = 2*pi*(3*surf.c*surf.a0dot[1]/(4*surf.uref) + surf.c*surf.adot[1]/(4*surf.uref) + surf.c*surf.adot[2]/(8*surf.uref))

    # Suction force given in eqn (2.31) Ramesh et al.
    cs = 2*pi*surf.a0[1]*surf.a0[1]

    #The components of normal force and moment from induced velocities are calulcated in dimensional units and nondimensionalized later
    nonl=0
    nonl_m=0
    for ib = 1:surf.ndiv-1
        nonl = nonl + (surf.uind[ib]*cos(surf.kinem.alpha) - surf.wind[ib]*sin(surf.kinem.alpha))*surf.bv[ib].s
        nonl_m = nonl_m + (surf.uind[ib]*cos(surf.kinem.alpha) - surf.wind[ib]*sin(surf.kinem.alpha))*surf.x[ib]*surf.bv[ib].s
    end
    nonl = nonl*2./(surf.uref*surf.uref*surf.c)
    nonl_m = nonl_m*2./(surf.uref*surf.uref*surf.c*surf.c)

    # Normal force coefficient
    cn = cnc + cnnc + nonl

    # Lift and drag coefficients 
    cl = cn*cos(surf.kinem.alpha) + cs*sin(surf.kinem.alpha)
    cd = cn*sin(surf.kinem.alpha)-cs*cos(surf.kinem.alpha)

    #Pitching moment is clockwise or nose up positive
    cm = cn*surf.pvt - 2*pi*((surf.kinem.u*cos(surf.kinem.alpha)/surf.uref + surf.kinem.hdot*sin(surf.kinem.alpha)/surf.uref)*(surf.a0[1]/4. + surf.aterm[1]/4. - surf.aterm[2]/8.) + (surf.c/surf.uref)*(7.*surf.a0dot[1]/16. + 3.*surf.adot[1]/16. + surf.adot[2]/16. - surf.adot[3]/64.)) - nonl_m
    return cl, cd, cm
end

function calc_forces(surf::TwoDSurf_2DOF)

    # First term in eqn (2.30) Ramesh et al. in coefficient form
    cnc = 2*pi*(surf.kinem.u*cos(surf.kinem.alpha)/surf.uref + surf.kinem.hdot*sin(surf.kinem.alpha)/surf.uref)*(surf.a0[1] + surf.aterm[1]/2.)

    # Second term in eqn (2.30) Ramesh et al. in coefficient form
    cnnc = 2*pi*(3*surf.c*surf.a0dot[1]/(4*surf.uref) + surf.c*surf.adot[1]/(4*surf.uref) + surf.c*surf.adot[2]/(8*surf.uref))

    # Suction force given in eqn (2.31) Ramesh et al.
    cs = 2*pi*surf.a0[1]*surf.a0[1]

    #The components of normal force and moment from induced velocities are calulcated in dimensional units and nondimensionalized later
    nonl=0
    nonl_m=0
    for ib = 1:surf.ndiv-1
        nonl = nonl + (surf.uind[ib]*cos(surf.kinem.alpha) - surf.wind[ib]*sin(surf.kinem.alpha))*surf.bv[ib].s
        nonl_m = nonl_m + (surf.uind[ib]*cos(surf.kinem.alpha) - surf.wind[ib]*sin(surf.kinem.alpha))*surf.x[ib]*surf.bv[ib].s
    end
    nonl = nonl*2./(surf.uref*surf.uref*surf.c)
    nonl_m = nonl_m*2./(surf.uref*surf.uref*surf.c*surf.c)

    # Normal force coefficient
    cn = cnc + cnnc + nonl

    # Lift and drag coefficients 
    cl = cn*cos(surf.kinem.alpha) + cs*sin(surf.kinem.alpha)
    cd = cn*sin(surf.kinem.alpha)-cs*cos(surf.kinem.alpha)

    #Pitching moment is clockwise or nose up positive
    cm = cn*surf.pvt - 2*pi*((surf.kinem.u*cos(surf.kinem.alpha)/surf.uref + surf.kinem.hdot*sin(surf.kinem.alpha)/surf.uref)*(surf.a0[1]/4. + surf.aterm[1]/4. - surf.aterm[2]/8.) + (surf.c/surf.uref)*(7.*surf.a0dot[1]/16. + 3.*surf.adot[1]/16. + surf.adot[2]/16. - surf.adot[3]/64.)) - nonl_m
    return cl, cd, cm
end


function calc_forces(surf::TwoDSurfwFlap, dt)

    # First term in eqn (2.30) Ramesh et al. in coefficient form
    # No longer required
    # cnc = 2*pi*(surf.kinem.u*cos(surf.kinem.alpha)/surf.uref + surf.kinem.hdot*sin(surf.kinem.alpha)/surf.uref)*(surf.a0[1] + surf.aterm[1]/2.)

    # Second term in eqn (2.30) Ramesh et al. in coefficient form
    # This term can be kept as is
    cnnc = 2*pi*(3*surf.c*surf.a0dot[1]/(4*surf.uref) + surf.c*surf.adot[1]/(4*surf.uref) + surf.c*surf.adot[2]/(8*surf.uref))

    # Suction force given in eqn (2.31) Ramesh et al.
    cs = 2*pi*surf.a0[1]*surf.a0[1]

    # The components of normal force and moment from induced velocities are calulcated in dimensional units and nondimensionalized later
    nonl=0
    nonl_m=0
    nonl_cnc = 0
    nonl_m1 = 0
    nonl_m2 = 0
    for ib = 1:surf.ndiv-1
        # 1st three terms in delta p, Ucos(),hdotsin(), alphadot*eta integrated
        q1 = sqrt(1+surf.cam_slope[ib]*surf.cam_slope[ib])*(surf.kinem.u*cos(surf.kinem.alpha)+surf.kinem.hdot*sin(surf.kinem.alpha) - surf.kinem.alphadot*surf.cam[ib])
        q2 = sqrt(1+surf.cam_slope[ib+1]*surf.cam_slope[ib+1])*(surf.kinem.u*cos(surf.kinem.alpha)+surf.kinem.hdot*sin(surf.kinem.alpha) - surf.kinem.alphadot*surf.cam[ib+1]) 
        nonl_cnc = nonl_cnc + 0.5*(q1 + q2)*surf.bv[ib].s  

        # Same term as before with additional squared spatial derivative term
        q1 = sqrt(1+surf.cam_slope[ib]*surf.cam_slope[ib])*(surf.uind[ib]*cos(surf.kinem.alpha) - surf.wind[ib]*sin(surf.kinem.alpha))
        q2 = sqrt(1+surf.cam_slope[ib+1]*surf.cam_slope[ib+1])*(surf.uind[ib+1]*cos(surf.kinem.alpha) - surf.wind[ib+1]*sin(surf.kinem.alpha))
        nonl = nonl + 0.5*(q1 + q2)*surf.bv[ib].s

        # Same term as before with additional squared spatial derivative term and alphadot*eta term
        q1 = sqrt(1+surf.cam_slope[ib]*surf.cam_slope[ib])*(surf.kinem.u*cos(surf.kinem.alpha)+surf.kinem.hdot*sin(surf.kinem.alpha) - surf.kinem.alphadot*surf.cam[ib])*surf.x[ib]
       q2 = sqrt(1+surf.cam_slope[ib+1]*surf.cam_slope[ib+1])*(surf.kinem.u*cos(surf.kinem.alpha)+surf.kinem.hdot*sin(surf.kinem.alpha) - surf.kinem.alphadot*surf.cam[ib+1])*surf.x[ib+1] 
        nonl_m1 = nonl_m1 + 0.5*(q1 + q2)*surf.bv[ib].s

        # Same term as before with additional squared spatial derivative term
        q1 = sqrt(1+surf.cam_slope[ib]*surf.cam_slope[ib])*(surf.uind[ib]*cos(surf.kinem.alpha) - surf.wind[ib]*sin(surf.kinem.alpha))*surf.x[ib]
        q2 = sqrt(1+surf.cam_slope[ib+1]*surf.cam_slope[ib+1])*(surf.uind[ib+1]*cos(surf.kinem.alpha) - surf.wind[ib+1]*sin(surf.kinem.alpha))*surf.x[ib+1]
        nonl_m = nonl_m + 0.5*(q1 + q2)*surf.bv[ib].s
    end
    
    # -------------------------------------------------------------
    xbdiv = indmin(abs(surf.x[:]-surf.x_b[1]))
    m_be1a = 0
    m_be1b = 0
    m_be2a = 0
    m_be2b = 0
    me_be3 = 0
    # These are the expressions multipled by the square root term in eqn in Hinge_Moment notebook. 
    for ib = xbdiv:surf.ndiv-1
        q1 = sqrt(1+surf.cam_slope[ib]*surf.cam_slope[ib])*(surf.kinem.u*cos(surf.kinem.alpha)+surf.kinem.hdot*sin(surf.kinem.alpha) - surf.kinem.alphadot*surf.cam[ib])
        q2 = sqrt(1+surf.cam_slope[ib+1]*surf.cam_slope[ib+1])*(surf.kinem.u*cos(surf.kinem.alpha)+surf.kinem.hdot*sin(surf.kinem.alpha) - surf.kinem.alphadot*surf.cam[ib+1])
        m_be1a = m_be1a + 0.5*(q1 + q2)*surf.bv[ib].s 

        q1 = sqrt(1+surf.cam_slope[ib]*surf.cam_slope[ib])*(surf.uind[ib]*cos(surf.kinem.alpha) - surf.wind[ib]*sin(surf.kinem.alpha))
        q2 = sqrt(1+surf.cam_slope[ib+1]*surf.cam_slope[ib+1])*(surf.uind[ib+1]*cos(surf.kinem.alpha) - surf.wind[ib+1]*sin(surf.kinem.alpha))
        m_be1b = m_be1b + 0.5*(q1 + q2)*surf.bv[ib].s

        q1 = sqrt(1+surf.cam_slope[ib]*surf.cam_slope[ib])*(surf.kinem.u*cos(surf.kinem.alpha)+surf.kinem.hdot*sin(surf.kinem.alpha) - surf.kinem.alphadot*surf.cam[ib])*surf.x[ib]
        q2 = sqrt(1+surf.cam_slope[ib+1]*surf.cam_slope[ib+1])*(surf.kinem.u*cos(surf.kinem.alpha)+surf.kinem.hdot*sin(surf.kinem.alpha) - surf.kinem.alphadot*surf.cam[ib+1])*surf.x[ib+1]
        m_be2a = m_be2a + 0.5*(q1 + q2)*surf.bv[ib].s

        q1 = sqrt(1+surf.cam_slope[ib]*surf.cam_slope[ib])*(surf.uind[ib]*cos(surf.kinem.alpha) - surf.wind[ib]*sin(surf.kinem.alpha))*surf.x[ib]
        q2 = sqrt(1+surf.cam_slope[ib+1]*surf.cam_slope[ib+1])*(surf.uind[ib+1]*cos(surf.kinem.alpha) - surf.wind[ib+1]*sin(surf.kinem.alpha))*surf.x[ib+1]
        m_be2b = m_be2b + 0.5*(q1 + q2)*surf.bv[ib].s
    end 

    # This block of code relates to the time derivate, double integration term.
    intg = zeros(surf.ndiv-xbdiv+1) 
    
    for ib = xbdiv:surf.ndiv
        sumbv = 0
        sumbv_prev = 0    
    	for i_bv = 1:ib-1
    	    sumbv = sumbv + surf.bv[i_bv].s
	    sumbv_prev = sumbv_prev + surf.bv_prev[i_bv].s
	end	 
        intg[ib-xbdiv+1] = (sumbv-sumbv_prev)/dt
    end
        
    m_be3a = trapz(intg[:],surf.x[xbdiv:surf.ndiv])
    m_be3b = trapz(intg[:].*surf.x[xbdiv:surf.ndiv],surf.x[xbdiv:surf.ndiv])

    m_be = surf.x_b[1]*(m_be1a + m_be1b) - (m_be2a + m_be2b) + surf.x_b[1]*m_be3a - m_be3b 

    cm_be = m_be*2./(surf.uref*surf.uref*surf.c*surf.c)
    # ---------------------------------------------------------------------------------

    nonl = nonl*2./(surf.uref*surf.uref*surf.c)
    nonl_cnc = nonl_cnc*2./(surf.uref*surf.uref*surf.c)
    nonl_m = nonl_m*2./(surf.uref*surf.uref*surf.c*surf.c)
    nonl_m1 = nonl_m1*2./(surf.uref*surf.uref*surf.c*surf.c)

    # Normal force coefficient
    cn = nonl_cnc + cnnc + nonl

    # Lift and drag coefficients 
    cl = cn*cos(surf.kinem.alpha) + cs*sin(surf.kinem.alpha)
    cd = cn*sin(surf.kinem.alpha)-cs*cos(surf.kinem.alpha)

    #Pitching moment is clockwise or nose up positive
    cm = cn*surf.pvt - 2*pi*((surf.c/surf.uref)*(7.*surf.a0dot[1]/16. + 3.*surf.adot[1]/16. + surf.adot[2]/16. - surf.adot[3]/64.)) - nonl_m1 - nonl_m
    return cl, cd, cm, cm_be
end
# ---------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------
# Transfer moment coefficient to different chosen locations
function transfer_cm(xreq::Float64, cm::Vector{Float64}, cl::Vector{Float64}, cd::Vector{Float64},alpha::Vector{Float64},  x::Float64,c)
    for i=1:length(cm)
        cm[i] = cm[i]+cl[i]*(xreq-x)*c*cos(alpha[i])+cd[i]*(xreq-x)*c*sin(alpha[i])
    end
    return cm
end
# ---------------------------------------------------------------------------------------------
    
