# ---------------------------------------------------------------------------------------------
# Function for viewing the 2D vorticity field around the airfoil
function view_vorts(surf::TwoDSurf, field::TwoDFlowField)
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

function calc_forces(surf::TwoDSurfwFlap)

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
        nonl_cnc = nonl_cnc + sqrt(1+surf.cam_slope[ib]*surf.cam_slope[ib])*(surf.kinem.u*cos(surf.kinem.alpha)+surf.kinem.hdot*sin(surf.kinem.alpha) - surf.kinem.alphadot*surf.cam[ib])*surf.bv[ib].s  
        # Same term as before with additional squared spatial derivative term
        nonl = nonl + sqrt(1+surf.cam_slope[ib]*surf.cam_slope[ib])*(surf.uind[ib]*cos(surf.kinem.alpha) - surf.wind[ib]*sin(surf.kinem.alpha))*surf.bv[ib].s
        # Same term as before with additional squared spatial derivative term and alphadot*eta term
        nonl_m1 = nonl_m1 + surf.x[ib]*sqrt(1+surf.cam_slope[ib]*surf.cam_slope[ib])*(surf.kinem.u*cos(surf.kinem.alpha)+surf.kinem.hdot*sin(surf.kinem.alpha) - surf.kinem.alphadot*surf.cam[ib])*surf.x[ib]surf.bv[ib].s
        # Same term as before with additional squared spatial derivative term
        nonl_m = nonl_m + surf.x[in]*sqrt(1+surf.cam_slope[ib]*surf.cam_slope[ib])*(surf.uind[ib]*cos(surf.kinem.alpha) - surf.wind[ib]*sin(surf.kinem.alpha))*surf.x[ib]*surf.bv[ib].s
    end
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
    return cl, cd, cm
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
    
