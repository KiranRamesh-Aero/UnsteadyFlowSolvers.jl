# ---------------------------------------------------------------------------------------------
# Function for viewing the 2D vorticity field around the airfoil
function view_vorts(surf::TwoDSurf, field::TwoDFlowField)
    scatter(map(q->q.x, field.tev),map(q->q.z,field.tev),s=20,c=map(q->q.s,field.tev),edgecolors="none")
    sc = scatter(map(q->q.x, field.lev),map(q->q.z,field.lev),s=20,c=map(q->q.s,field.lev),edgecolors="none")
    plot(map(q->q.x, surf.bv),map(q->q.z,surf.bv),color = "black",linewidth=2.0)
end

function view_vorts(surf::Vector{TwoDSurf}, field::TwoDFlowFieldMultSurf)
    scatter([map(q->q.x, field.tev[i]) for i = 1:length(field.tev)] , [map(q->q.z,field.tev[i]) for i = 1:length(field.tev)], s=20, c=[map(q->q.s,field.tev[i]) for i = 1:length(field.tev)], edgecolors="none")
    scatter([map(q->q.x, field.lev[i]) for i = 1:length(field.lev)] , [map(q->q.z,field.lev[i]) for i = 1:length(field.lev)], s=20, c=[map(q->q.s,field.lev[i]) for i = 1:length(field.lev)], edgecolors="none")
    for i = 1:field.nsurf
        plot(map(q->q.x, surf[i].bv), map(q->q.z,surf[i].bv), color = "black", linewidth=2.0)
    end
    PyPlot.show()
end

function view_vorts(surf::TwoDSurf_2DOF, field::TwoDFlowField)
    scatter(map(q->q.x, field.tev),map(q->q.z,field.tev),s=20,c=map(q->q.z,field.tev),cmap=ColorMap("jet"))
    sc = scatter(map(q->q.x, field.lev),map(q->q.z,field.lev),s=20,c=map(q->q.z,field.lev),cmap=ColorMap("jet"))
    plot(map(q->q.x, surf.bv),map(q->q.z,surf.bv),color = "black",linewidth=1.5)
    colorbar(sc)
end

function view_vorts(surf::TwoDFreeSurf, field::TwoDFlowField)
    scatter(map(q->q.x, field.tev),map(q->q.z,field.tev),s=20,c=map(q->q.s,field.tev),cmap=ColorMap("jet"),edgecolors="none")
    scatter(map(q->q.x, field.lev),map(q->q.z,field.lev),s=20,c=map(q->q.s,field.lev),cmap=ColorMap("jet"),edgecolors="none")
    scatter(map(q->q.x, field.extv),map(q->q.z,field.extv),s=40,c=map(q->q.s,field.extv),cmap=ColorMap("jet"),edgecolors="none")
    plot(map(q->q.x, surf.bv),map(q->q.z,surf.bv),color = "black",linewidth=2.0)

end

function view_vorts(surf::TwoDSurfwFlap, field::TwoDFlowField)
     scatter(map(q->q.x, field.tev),map(q->q.z,field.tev),s=20,c=map(q->q.s,field.tev),cmap=ColorMap("jet"),edgecolors="none")
    sc = scatter(map(q->q.x, field.lev),map(q->q.z,field.lev),s=20,c=map(q->q.s,field.lev),cmap=ColorMap("jet"),edgecolors="none")
    plot(map(q->q.x, surf.bv),map(q->q.z,surf.bv),color = "black",linewidth=2.0)
colorbar(sc)
end

function view_vorts(surf::TwoDSurfLV, field::TwoDFlowField)
     scatter(map(q->q.x, field.tev),map(q->q.z,field.tev),s=20,c=map(q->q.s,field.tev),cmap=ColorMap("jet"),edgecolors="none")
    sc = scatter(map(q->q.x, field.lev),map(q->q.z,field.lev),s=20,c=map(q->q.s,field.lev),cmap=ColorMap("jet"),edgecolors="none")
    plot(map(q->q.xv, surf.lv),map(q->q.zv,surf.lv),color = "black",linewidth=2.0)
colorbar(sc)
end


function view_vorts(field::Vector{TwoDVort})
    scatter(map(q->q.x, field),map(q->q.z,field),s=20,c=map(q->q.s,field),cmap=ColorMap("jet"),edgecolors="none")
    colorbar(sc)
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

function calc_forces(surf::Vector{TwoDSurf})

    nsurf = length(surf)
    cl = zeros(nsurf)
    cd = zeros(nsurf)
    cm = zeros(nsurf)
    
    for i = 1:nsurf
        # First term in eqn (2.30) Ramesh et al. in coefficient form
        cnc = 2*pi*(surf[i].kinem.u*cos(surf[i].kinem.alpha)/surf[i].uref + surf[i].kinem.hdot*sin(surf[i].kinem.alpha)/surf[i].uref)*(surf[i].a0[1] + surf[i].aterm[1]/2.)
        
        # Second term in eqn (2.30) Ramesh et al. in coefficient form
        cnnc = 2*pi*(3*surf[i].c*surf[i].a0dot[1]/(4*surf[i].uref) + surf[i].c*surf[i].adot[1]/(4*surf[i].uref) + surf[i].c*surf[i].adot[2]/(8*surf[i].uref))
        
        # Suction force given in eqn (2.31) Ramesh et al.
        cs = 2*pi*surf[i].a0[1]*surf[i].a0[1]
        
        #The components of normal force and moment from induced velocities are calulcated in dimensional units and nondimensionalized later
        nonl=0
        nonl_m=0
        for ib = 1:surf[i].ndiv-1
            nonl = nonl + (surf[i].uind[ib]*cos(surf[i].kinem.alpha) - surf[i].wind[ib]*sin(surf[i].kinem.alpha))*surf[i].bv[ib].s
            nonl_m = nonl_m + (surf[i].uind[ib]*cos(surf[i].kinem.alpha) - surf[i].wind[ib]*sin(surf[i].kinem.alpha))*surf[i].x[ib]*surf[i].bv[ib].s
    end
        nonl = nonl*2./(surf[i].uref*surf[i].uref*surf[i].c)
        nonl_m = nonl_m*2./(surf[i].uref*surf[i].uref*surf[i].c*surf[i].c)
        
        # Normal force coefficient
        cn = cnc + cnnc + nonl
        
        # Lift and drag coefficients
        cl[i] = cn*cos(surf[i].kinem.alpha) + cs*sin(surf[i].kinem.alpha)
        cd[i] = cn*sin(surf[i].kinem.alpha) - cs*cos(surf[i].kinem.alpha)
        
        #Pitching moment is clockwise or nose up positive
        cm[i] = cn*surf[i].pvt - 2*pi*((surf[i].kinem.u*cos(surf[i].kinem.alpha)/surf[i].uref + surf[i].kinem.hdot*sin(surf[i].kinem.alpha)/surf[i].uref)*(surf[i].a0[1]/4. + surf[i].aterm[1]/4. - surf[i].aterm[2]/8.) + (surf[i].c/surf[i].uref)*(7.*surf[i].a0dot[1]/16. + 3.*surf[i].adot[1]/16. + surf[i].adot[2]/16. - surf[i].adot[3]/64.)) - nonl_m
    end
    return cl, cd, cm
end

function calc_forces_more(surf::TwoDSurf)

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
    cm_p = 2*pi*((surf.kinem.u*cos(surf.kinem.alpha)/surf.uref + surf.kinem.hdot*sin(surf.kinem.alpha)/surf.uref)*(surf.a0[1]/4. + surf.aterm[1]/4. - surf.aterm[2]/8.) + (surf.c/surf.uref)*(7.*surf.a0dot[1]/16. + 3.*surf.adot[1]/16. + surf.adot[2]/16. - surf.adot[3]/64.))
    cm = cn*surf.pvt -  cm_p - nonl_m
    
    return cl, cd, cm, surf.a0[1] + 0.5*surf.aterm[1], cn, cs, cnc, cnnc, nonl, cn*surf.pvt, cm_p, nonl_m
end

function calc_forces_E(surf::TwoDSurf, lev :: Float64, dt :: Float64)

    # First term in eqn (2.30) Ramesh et al. in coefficient form
    cnc = 2*pi*(surf.kinem.u*cos(surf.kinem.alpha)/surf.uref + surf.kinem.hdot*sin(surf.kinem.alpha)/surf.uref)*(surf.a0[1] + surf.aterm[1]/2.)

    # Second term in eqn (2.30) Ramesh et al. in coefficient form
    cnnc = 2*pi*(3*surf.c*surf.a0dot[1]/(4*surf.uref) + surf.c*surf.adot[1]/(4*surf.uref) + surf.c*surf.adot[2]/(8*surf.uref)) + (2*lev/(dt*surf.uref*surf.uref)) 

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
    cm = cn*surf.pvt - 2*pi*((surf.kinem.u*cos(surf.kinem.alpha)/surf.uref + surf.kinem.hdot*sin(surf.kinem.alpha)/surf.uref)*(surf.a0[1]/4. + surf.aterm[1]/4. - surf.aterm[2]/8.) + (surf.c/surf.uref)*(7.*surf.a0dot[1]/16. + 3.*surf.adot[1]/16. + surf.adot[2]/16. - surf.adot[3]/64.)) - nonl_m - lev*(2*surf.pvt - 1)/(dt*surf.uref*surf.uref)
    return cl, cd, cm
end

function calc_forces_E_more(surf::TwoDSurf, lev :: Float64, dt :: Float64)

    # First term in eqn (2.30) Ramesh et al. in coefficient form
    cnc = 2*pi*(surf.kinem.u*cos(surf.kinem.alpha)/surf.uref + surf.kinem.hdot*sin(surf.kinem.alpha)/surf.uref)*(surf.a0[1] + surf.aterm[1]/2.)

    # Second term in eqn (2.30) Ramesh et al. in coefficient form
    cnnc = 2*pi*(3*surf.c*surf.a0dot[1]/(4*surf.uref) + surf.c*surf.adot[1]/(4*surf.uref) + surf.c*surf.adot[2]/(8*surf.uref)) + (2*lev/(dt*surf.uref*surf.uref)) 

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
    cm = cn*surf.pvt - 2*pi*((surf.kinem.u*cos(surf.kinem.alpha)/surf.uref + surf.kinem.hdot*sin(surf.kinem.alpha)/surf.uref)*(surf.a0[1]/4. + surf.aterm[1]/4. - surf.aterm[2]/8.) + (surf.c/surf.uref)*(7.*surf.a0dot[1]/16. + 3.*surf.adot[1]/16. + surf.adot[2]/16. - surf.adot[3]/64.)) - nonl_m - lev*(2*surf.pvt - 1)/(dt*surf.uref*surf.uref)
    return cl, cd, cm, surf.a0[1] + 0.5*surf.aterm[1], cn, cs
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


function calc_forces_E(surf::TwoDSurf_2DOF, lev :: Float64, dt :: Float64)

    # First term in eqn (2.30) Ramesh et al. in coefficient form
    cnc = 2*pi*(surf.kinem.u*cos(surf.kinem.alpha)/surf.uref + surf.kinem.hdot*sin(surf.kinem.alpha)/surf.uref)*(surf.a0[1] + surf.aterm[1]/2.)

    # Second term in eqn (2.30) Ramesh et al. in coefficient form
    cnnc = 2*pi*(3*surf.c*surf.a0dot[1]/(4*surf.uref) + surf.c*surf.adot[1]/(4*surf.uref) + surf.c*surf.adot[2]/(8*surf.uref)) + (2*lev/(dt*surf.uref*surf.uref)) 

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
    cm = cn*surf.pvt - 2*pi*((surf.kinem.u*cos(surf.kinem.alpha)/surf.uref + surf.kinem.hdot*sin(surf.kinem.alpha)/surf.uref)*(surf.a0[1]/4. + surf.aterm[1]/4. - surf.aterm[2]/8.) + (surf.c/surf.uref)*(7.*surf.a0dot[1]/16. + 3.*surf.adot[1]/16. + surf.adot[2]/16. - surf.adot[3]/64.)) - nonl_m - lev*(2*surf.pvt - 1)/(dt*surf.uref*surf.uref)
    return cl, cd, cm
end


function calc_forces(surf::TwoDFreeSurf)

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


# Transfer moment coefficient to different chosen locations
function transfer_cm(xreq::Float64, cm::Vector{Float64}, cl::Vector{Float64}, cd::Vector{Float64},alpha::Vector{Float64},  x::Float64,c)
    for i=1:length(cm)
        cm[i] = cm[i]+cl[i]*(xreq-x)*c*cos(alpha[i])+cd[i]*(xreq-x)*c*sin(alpha[i])
    end
    return cm
end
# ---------------------------------------------------------------------------------------------

function anim_flow(outfolder, freq)

    #Outfolder contains the movie files
    n = length(readdir(outfolder))

    figure(1)

    for i = 1:n
        data = readdlm(string(outfolder,"/field.$(freq*i)"))
        rtev = 2:1+Int(data[1,1])
        rlev = 2+Int(data[1,1]):1+Int(data[1,1]+data[1,2])
        rextv = 2+Int(data[1,1]+data[1,2]):1+Int(data[1,1]+data[1,2]+data[1,3])
        rbv = 2+Int(data[1,1]+data[1,2]+data[1,3]):1+Int(data[1,1]+data[1,2]+data[1,3]+data[1,4])
        im = scatter(data[rtev,2],data[rtev,3],s=20,c=data[rtev,1],cmap=PyPlot.ColorMap("jet"),edgecolors="none");
        scatter(data[rlev,2],data[rlev,3],s=20,c=data[rlev,1],cmap=PyPlot.ColorMap("jet"),edgecolors="none");
        scatter(data[rextv,2],data[rextv,3],s=20,c=data[rextv,1],cmap=PyPlot.ColorMap("jet"),edgecolors="none");
        plot(data[rbv,2],data[rbv,3],color = "black",linewidth=2.0)
        PyPlot.savefig(outfolder * "/" * string(i) * ".png")
        PyPlot.clf()
    end

    run(`ffmpeg -r 25 -i $outfolder/%d.png $outfolder/anim.mpg`)
    run(`open $outfolder/anim.mpg`)
end
# ---------------------------------------------------------------------------------------------



function write_stamp(surf :: TwoDSurf, curfield :: TwoDFlowField, t :: Float64, kelv_enf :: Float64, g:: JLD.JldGroup)

    cl, cd, cm, gamma, cn, cs, cnc, cnnc, nonl, cm_n, cm_pvt, nonl_m = calc_forces_more(surf)
    
    tevmat = zeros(length(curfield.tev), 3)
    for i = 1:length(curfield.tev)
        tevmat[i,:] = [curfield.tev[i].s curfield.tev[i].x curfield.tev[i].z]
    end
    g["tev"] = tevmat
    levmat = zeros(length(curfield.tev), 3)
    for i = 1:length(curfield.lev)
        levmat[i,:] = [curfield.lev[i].s curfield.lev[i].x curfield.lev[i].z]
    end
    g["lev"] = levmat
    bvmat = zeros(length(surf.bv), 3)
    for i = 1:length(surf.bv)
        bvmat[i,:] = [surf.bv[i].s surf.bv[i].x surf.bv[i].z]
    end
    g["bv"] = bvmat
    extvmat = zeros(length(curfield.extv), 3)
    for i = 1:length(curfield.extv)
        extvmat[i,:] = [curfield.extv[i].s curfield.extv[i].x curfield.extv[i].z]
    end
    g["extv"] = extvmat
    g["t"] = t
    g["alpha"] = surf.kinem.alpha
    g["h"] = surf.kinem.h
    g["u"] = surf.kinem.u
    g["alphadot"] = surf.kinem.alphadot
    g["hdot"] = surf.kinem.hdot
    g["udot"] = surf.kinem.udot
    g["a0"] = surf.a0[1]
    g["a0dot"] = surf.a0dot[1]
    g["aterm"] = surf.aterm
    g["adot"] = surf.adot
    g["levflag"] = surf.levflag[1]
    g["kelv_enf"] = kelv_enf
    g["cl"] = cl
    g["cd"] = cd
    g["cm"] = cm
    g["gamma"] = gamma
    g["cn"] = cn
    g["cs"] = cs
    g["cnc"] = cnc
    g["cnnc"] = cnnc
    g["nonl_n"] = nonl
    g["cm_n"] = cm_n
    g["cm_pvt"] = cm_pvt
    g["nonl_m"] = nonl_m
end
# ---------------------------------------------------------------------------------------------
#Postprocessing for the 3D vortex ring method
function get_gridprop(surf:: ThreeDSurfVR, prop::String)
    #Obtains the required property in a matrix form
    mat = zeros(surf.nchord, surf.nspan)
    if prop == "xc"
        for i = 1:surf.nchord
            for j = 1:surf.nspan
                np = (i - 1)*surf.nspan + j
                mat[i,j] = surf.vr_p[np].xc
            end
        end
    elseif prop == "yc"
        for i = 1:surf.nchord
            for j = 1:surf.nspan
                np = (i - 1)*surf.nspan + j
                mat[i,j] = surf.vr_p[np].yc
            end
        end
    elseif prop == "zc"
        for i = 1:surf.nchord
            for j = 1:surf.nspan
                np = (i - 1)*surf.nspan + j
                mat[i,j] = surf.vr_p[np].zc
            end
        end
    elseif prop == "s"
        for i = 1:surf.nchord
            for j = 1:surf.nspan
                np = (i - 1)*surf.nspan + j
                mat[i,j] = surf.vr_p[np].s
            end
        end
    elseif prop == "xc_I"
        for i = 1:surf.nchord
            for j = 1:surf.nspan
                np = (i - 1)*surf.nspan + j
                mat[i,j] = surf.vr_p[np].xc_I
            end
        end
    elseif prop == "yc_I"
        for i = 1:surf.nchord
            for j = 1:surf.nspan
                np = (i - 1)*surf.nspan + j
                mat[i,j] = surf.vr_p[np].yc_I
            end
        end
    elseif prop == "zc_I"
        for i = 1:surf.nchord
            for j = 1:surf.nspan
                np = (i - 1)*surf.nspan + j
                mat[i,j] = surf.vr_p[np].zc
            end
        end
    
    else
        error("Invalid property specification")
    end
    return mat
end

function get_gridprop(field:: ThreeDFlowFieldVR, vort::String, prop::String)
    #Obtains the required property in a matrix form
 
    
    if vort == "tev"
        ntev = Int(length(field.tev_s)/field.nspan)
        mat = zeros(ntev, field.nspan)
        if prop == "s"
            for i = 1:ntev
                for j = 1:field.nspan
                    np = (i - 1)*field.nspan + j
                    mat[i,j] = field.tev_s[np]
                end
            end
        else
            error("Invalid property specification")
        end
    else
        error("Invalid property specification")
    end
    return mat
end

        
                     
                
function plotgrid(surf :: ThreeDSurfVR, field :: ThreeDFlowFieldVR)
    #Plot the grid
    #Plot chordwise lines
    x = map(q -> q.xv_I, surf.vr_g[:])
    y = map(q -> q.yv_I, surf.vr_g[:])
    z = map(q -> q.zv_I, surf.vr_g[:])
    x = reshape(x, surf.nspan+1, surf.nchord+1)'
    y = reshape(y, surf.nspan+1, surf.nchord+1)'
    z = reshape(z, surf.nspan+1, surf.nchord+1)'

    PyPlot.plot_wireframe(x,y,z)

    ntev = Int(length(curfield.tev)/(surf.nspan + 1)) - 1
    x = map(q -> q.xv_I, curfield.tev)
    y = map(q -> q.yv_I, curfield.tev)
    z = map(q -> q.zv_I, curfield.tev)
    x = reshape(x, surf.nspan+1, ntev+1)'
    y = reshape(y, surf.nspan+1, ntev+1)'
    z = reshape(z, surf.nspan+1, ntev+1)'

    PyPlot.plot_wireframe(x,y,z)

    
    PyPlot.show()
    
end

# ---------------------------------------------------------------------------------------------
