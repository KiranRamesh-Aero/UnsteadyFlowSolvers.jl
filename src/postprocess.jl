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

function view_vorts(surf::ThreeDSurfSimple, field::ThreeDFieldStrip)
    fig = figure()
    ax = PyPlot.Axes3D(fig)
    for i = 1:surf.nspan
        plot(map(q->q.x, surf.s2d[i].bv), ones(surf.ndiv-1)*surf.yle[i], map(q->q.z,surf.s2d[i].bv), color = "black",linewidth=2.0)
        
        x = map(p->p.x, map(q->q[:][i], field.tev))
        z = map(p->p.z, map(q->q[:][i], field.tev))
        y = ones(length(field.tev))*surf.yle[i]
        s = map(p->p.s, map(q->q[:][i], field.tev))
        scatter3D(x, y, z, s=20, c=s, cmap=ColorMap("jet"), edgecolors="none")
        
        if length(field.lev) > 0
            x = map(p->p.x, map(q->q[:][i], field.lev))
            z = map(p->p.z, map(q->q[:][i], field.lev))
            y = ones(length(field.lev))*surf.yle[i]
            s = map(p->p.s, map(q->q[:][i], field.lev))
            scatter3D(x, y, z, s=20, c=s, cmap=ColorMap("jet"), edgecolors="none")
        end
    end
    PyPlot.show()
end    

function view_vorts(surf::ThreeDSurfWeiss, field::ThreeDFieldStrip)
    fig = figure()
    ax = PyPlot.Axes3D(fig)
    for i = 1:surf.nlat
        plot(map(q->q.x, surf.s2d[i].bv), ones(surf.ndiv-1)*surf.yle[i], map(q->q.z,surf.s2d[i].bv), color = "black",linewidth=2.0)
        
        x = map(p->p.x, map(q->q[:][i], field.tev))
        z = map(p->p.z, map(q->q[:][i], field.tev))
        y = ones(length(field.tev))*surf.yle[i]
        s = map(p->p.s, map(q->q[:][i], field.tev))
        scatter3D(x, y, z, s=20, c=s, cmap=ColorMap("jet"), edgecolors="none")
        
        if length(field.lev) > 0
            x = map(p->p.x, map(q->q[:][i], field.lev))
            z = map(p->p.z, map(q->q[:][i], field.lev))
            y = ones(length(field.lev))*surf.yle[i]
            s = map(p->p.s, map(q->q[:][i], field.lev))
            scatter3D(x, y, z, s=20, c=s, cmap=ColorMap("jet"), edgecolors="none")
        end
    end
    PyPlot.show()
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
  #  return cl, cd, cm
	return cl, cd, cm, surf.a0[1] + 0.5*surf.aterm[1], cn, cs, cnc, cnnc, nonl, cn*surf.pvt, nonl_m
   
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

function calc_forces(surf::TwoDSurf, fsep :: Float64, sepdef::SeparationParams)


    # First term in eqn (2.30) Ramesh et al. in coefficient form
    cnc = 2*pi*(surf.kinem.u*cos(surf.kinem.alpha)/surf.uref + surf.kinem.hdot*sin(surf.kinem.alpha)/surf.uref)*(surf.a0[1] + surf.aterm[1]/2.)
 
    # Second term in eqn (2.30) Ramesh et al. in coefficient form
    cnnc = 2*pi*(3*surf.c*surf.a0dot[1]/(4*surf.uref) + surf.c*surf.adot[1]/(4*surf.uref) + surf.c*surf.adot[2]/(8*surf.uref))

    # Suction force given in eqn (2.31) Ramesh et al. and accounting for separation
	
	#Cs model 
	if(sepdef.model == "Sheng")
		fcrit = 0.6
	else	
		fcrit = 0.7
	end
	
	if(sepdef.cs_model == "Sheng_continuous")
		eta = 1.0
		E_0 = 0.2
	
		cs = eta*2*pi*surf.a0[1]*surf.a0[1]*(sqrt(fsep)-E_0)	
	
	elseif(sepdef.cs_model == "Sheng_piecewise")
	
		eta = 1.0
		E_0 = 0.2
	
		if(fsep<fcrit)
			cs = eta*2*pi*surf.a0[1]*surf.a0[1]*((fsep)^(3/2)-E_0)
		else
			cs = eta*2*pi*surf.a0[1]*surf.a0[1]*(sqrt(fsep)-E_0)	
		end
	
	elseif(sepdef.cs_model == "piecewise")
	
		eta = 1.0
		E_0 = 0.0
		
		if(fsep<fcrit)
			cs = eta*2*pi*surf.a0[1]*surf.a0[1]*((fsep)^(3/2)-E_0)
		else
			cs = eta*2*pi*surf.a0[1]*surf.a0[1]*(sqrt(fsep)-E_0)	
		end
	
	
	elseif(sepdef.cs_model == "continuous")
	
		eta = 1.0
		E_0 = 0.0
	
		cs = eta*2*pi*surf.a0[1]*surf.a0[1]*(sqrt(fsep)-E_0)	
	else
	println("Wrong Cs model specified. Choose Sheng_continuous, Sheng_piecewise, piecewise or continuous")
	end

    #The components of normal force and moment from induced velocities are calulcated in dimensional units and nondimensionalized later
    nonl=0
    nonl_m=0
    for ib = 1:surf.ndiv-1
        nonl = nonl + (surf.uind[ib]*cos(surf.kinem.alpha) - surf.wind[ib]*sin(surf.kinem.alpha))*surf.bv[ib].s
        nonl_m = nonl_m + (surf.uind[ib]*cos(surf.kinem.alpha) - surf.wind[ib]*sin(surf.kinem.alpha))*surf.x[ib]*surf.bv[ib].s
	end

	nonl = nonl*2./(surf.uref*surf.uref*surf.c)
    nonl_m = nonl_m*2./(surf.uref*surf.uref*surf.c*surf.c)

    #Account for separation
    cnsep = (cnc+nonl)*((1 + sqrt(fsep))/2.)^2
    
    # Normal force coefficient
    cn = cnsep + cnnc
    
    # Lift and drag coefficients
    cl = cn*cos(surf.kinem.alpha) + cs*sin(surf.kinem.alpha)
    cd = cn*sin(surf.kinem.alpha)-cs*cos(surf.kinem.alpha)
 
 cm = cn*surf.pvt - 2*pi*((surf.kinem.u*cos(surf.kinem.alpha)/surf.uref + surf.kinem.hdot*sin(surf.kinem.alpha)/surf.uref)*(surf.a0[1]/4. + surf.aterm[1]/4. - surf.aterm[2]/8.)
 + (surf.c/surf.uref)*(7.*surf.a0dot[1]/16. + 3.*surf.adot[1]/16. + surf.adot[2]/16. - surf.adot[3]/64.)) - nonl_m
    
	
	if(sepdef.cm_model == "Dymore")
				#Pitching moment is clockwise or nose up positive
		cm1 = (sepdef.k0+surf.pvt + surf.c*sepdef.k1*(1 - fsep) + sepdef.k2*surf.c*sin(pi*(fsep)^sepdef.m))*cnsep + surf.pvt*cnnc 
	  #  cm2 = -((1 + sqrt(fsep))^2)*2*pi*(surf.kinem.u*cos(surf.kinem.alpha)/surf.uref + surf.kinem.hdot*sin(surf.kinem.alpha)/surf.uref)*(surf.a0[1]/4. + surf.aterm[1]/4. - surf.aterm[2]/8.)
		cm2 = -2*pi*(surf.kinem.u*cos(surf.kinem.alpha)/surf.uref + surf.kinem.hdot*sin(surf.kinem.alpha)/surf.uref)*(surf.a0[1]/4. + surf.aterm[1]/4. - surf.aterm[2]/8.)
		#cm3 = -2*pi*surf.c/surf.uref*(7.*surf.a0dot[1]/16. + 3.*surf.adot[1]/16. + surf.adot[2]/16. - surf.adot[3]/64.)
		cm3 = -2*pi*surf.c/surf.uref*(7.*surf.a0dot[1]/16. + 3.*surf.adot[1]/16. + surf.adot[2]/16. - surf.adot[3]/64.)
		#cm4 = -((1 + sqrt(fsep))^2)*nonl_m
		cm4 = -nonl_m
	elseif (sepdef.cm_model == "Liu")
		#Pitching moment is clockwise or nose up positive
		cm1 = (sepdef.k0+surf.pvt + surf.c*sepdef.k1*(1 - fsep) + sepdef.k2*surf.c*sin(pi*(fsep)^sepdef.m))*cnsep + surf.pvt*cnnc  
		cm2 = -((1 + sqrt(fsep))^2)*2*pi*(surf.kinem.u*cos(surf.kinem.alpha)/surf.uref + surf.kinem.hdot*sin(surf.kinem.alpha)/surf.uref)*(surf.a0[1]/4. + surf.aterm[1]/4. - surf.aterm[2]/8.)
		cm3 = -2*pi*surf.c/surf.uref*(7.*surf.a0dot[1]/16. + 3.*surf.adot[1]/16. + surf.adot[2]/16. - surf.adot[3]/64.)
		cm4 = -((1 + sqrt(fsep))^2)*nonl_m
	else
		println("Wrong model specified.")
	end
	                                      
    cm = cm1 + cm2 + cm3 + cm4 

  #  return cl, cd, cm
   # return cl, cd, cm, surf.a0[1] + 0.5*surf.aterm[1], cn, cs, cnc, cnnc, nonl, cn*surf.pvt, cm_p, nonl_m
   return cl, cd, cm, surf.a0[1] + 0.5*surf.aterm[1], cn, cs, cnc, cnnc, nonl, cn*surf.pvt, nonl_m, cm1, cm2, cm3, cm4

end

function calc_forces(surf::TwoDSurf_2DOF, fsep :: Float64)

    #Constants for moment calculation
    k1 = -0.135
    k2 = 0.04
    
    # First term in eqn (2.30) Ramesh et al. in coefficient form
    cnc = 2*pi*(surf.kinem.u*cos(surf.kinem.alpha)/surf.uref + surf.kinem.hdot*sin(surf.kinem.alpha)/surf.uref)*(surf.a0[1] + surf.aterm[1]/2.)

    # Second term in eqn (2.30) Ramesh et al. in coefficient form
    cnnc = 2*pi*(3*surf.c*surf.a0dot[1]/(4*surf.uref) + surf.c*surf.adot[1]/(4*surf.uref) + surf.c*surf.adot[2]/(8*surf.uref))

    # Suction force given in eqn (2.31) Ramesh et al. and accounting for separation
    cs = sqrt(fsep)*2*pi*surf.a0[1]*surf.a0[1]
    
    
    #The components of normal force and moment from induced velocities are calulcated in dimensional units and nondimensionalized later
    nonl=0
    nonl_m=0
    for ib = 1:surf.ndiv-1
        nonl = nonl + (surf.uind[ib]*cos(surf.kinem.alpha) - surf.wind[ib]*sin(surf.kinem.alpha))*surf.bv[ib].s
        nonl_m = nonl_m + (surf.uind[ib]*cos(surf.kinem.alpha) - surf.wind[ib]*sin(surf.kinem.alpha))*surf.x[ib]*surf.bv[ib].s
    end
    nonl = nonl*2./(surf.uref*surf.uref*surf.c)
    nonl_m = nonl_m*2./(surf.uref*surf.uref*surf.c*surf.c)

    #Account for separation
    cnsep = (cnc+nonl)*((1 + sqrt(fsep))/2.)^2
    
    # Normal force coefficient
    cn = cnsep + cnnc
    
    # Lift and drag coefficients
    cl = cn*cos(surf.kinem.alpha) + cs*sin(surf.kinem.alpha)
    cd = cn*sin(surf.kinem.alpha)-cs*cos(surf.kinem.alpha)

    #Pitching moment is clockwise or nose up positive
    cm1 = (surf.pvt + k1*(1 - fsep) + k2*sin(pi*(fsep)^2))*cnsep + surf.pvt*cnnc
    cm2 = -((1 + sqrt(fsep))^2)*2*pi*(surf.kinem.u*cos(surf.kinem.alpha)/surf.uref + surf.kinem.hdot*sin(surf.kinem.alpha)/surf.uref)*(surf.a0[1]/4. + surf.aterm[1]/4. - surf.aterm[2]/8.)
    cm3 = -2*pi*surf.c/surf.uref*(7.*surf.a0dot[1]/16. + 3.*surf.adot[1]/16. + surf.adot[2]/16. - surf.adot[3]/64.)
    cm4 = -((1 + sqrt(fsep))^2)*nonl_m 

                                      
    cm = cm1 + cm2 + cm3 + cm4 
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

function calc_forces_more(surf::Vector{TwoDSurf})

    nsurf = length(surf)
    cl = zeros(nsurf)
    cd = zeros(nsurf)
    cm = zeros(nsurf)
    cnc = zeros(nsurf)
    cnnc = zeros(nsurf)
    cn = zeros(nsurf)
    cs = zeros(nsurf)
    nonl = zeros(nsurf)
    cm_n = zeros(nsurf)
    cm_p = zeros(nsurf)
    nonl_m = zeros(nsurf)
    bc = zeros(nsurf)
    
    for i = 1:nsurf
        # First term in eqn (2.30) Ramesh et al. in coefficient form
        cnc[i] = 2*pi*(surf[i].kinem.u*cos(surf[i].kinem.alpha)/surf[i].uref + surf[i].kinem.hdot*sin(surf[i].kinem.alpha)/surf[i].uref)*(surf[i].a0[1] + surf[i].aterm[1]/2.)
        
        # Second term in eqn (2.30) Ramesh et al. in coefficient form
        cnnc[i] = 2*pi*(3*surf[i].c*surf[i].a0dot[1]/(4*surf[i].uref) + surf[i].c*surf[i].adot[1]/(4*surf[i].uref) + surf[i].c*surf[i].adot[2]/(8*surf[i].uref))
        
        # Suction force given in eqn (2.31) Ramesh et al.
        cs[i] = 2*pi*surf[i].a0[1]*surf[i].a0[1]
        
        #Bound circulation
        bc[i] = surf[i].uref*surf[i].c*pi*(surf[i].a0[1] + 0.5*surf[i].aterm[1])
        
        #The components of normal force and moment from induced velocities are calulcated in dimensional units and nondimensionalized later
        nonl[i]=0
        nonl_m[i]=0
        for ib = 1:surf[i].ndiv-1
            nonl[i] = nonl[i] + (surf[i].uind[ib]*cos(surf[i].kinem.alpha) - surf[i].wind[ib]*sin(surf[i].kinem.alpha))*surf[i].bv[ib].s
            nonl_m[i] = nonl_m[i] + (surf[i].uind[ib]*cos(surf[i].kinem.alpha) - surf[i].wind[ib]*sin(surf[i].kinem.alpha))*surf[i].x[ib]*surf[i].bv[ib].s
        end
        nonl[i] = nonl[i]*2./(surf[i].uref*surf[i].uref*surf[i].c)
        nonl_m[i] = nonl_m[i]*2./(surf[i].uref*surf[i].uref*surf[i].c*surf[i].c)
        
        # Normal force coefficient
        cn[i] = cnc[i] + cnnc[i] + nonl[i]
        
        # Lift and drag coefficients
        cl[i] = cn[i]*cos(surf[i].kinem.alpha) + cs[i]*sin(surf[i].kinem.alpha)
        cd[i] = cn[i]*sin(surf[i].kinem.alpha) - cs[i]*cos(surf[i].kinem.alpha)
        
        #Pitching moment is clockwise or nose up positive
        cm_n[i] = cn[i]*surf[i].pvt
        cm_p[i] = 2*pi*((surf[i].kinem.u*cos(surf[i].kinem.alpha)/surf[i].uref + surf[i].kinem.hdot*sin(surf[i].kinem.alpha)/surf[i].uref)*(surf[i].a0[1]/4. + surf[i].aterm[1]/4. - surf[i].aterm[2]/8.) + (surf[i].c/surf[i].uref)*(7.*surf[i].a0dot[1]/16. + 3.*surf[i].adot[1]/16. + surf[i].adot[2]/16. - surf[i].adot[3]/64.))
        cm[i] = cn[i]*surf[i].pvt -  cm_p[i] - nonl_m[i]
    end
    
    return cl, cd, cm, bc, cn, cs, cnc, cnnc, nonl, cm_n, cm_p, nonl_m
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

function calc_forces_E(surf::Vector{TwoDSurf}, lev :: Vector{Float64}, dt :: Float64, shedv :: Vector{Int})
    nsurf = length(surf)
    cl = zeros(nsurf)
    cd = zeros(nsurf)
    cm = zeros(nsurf)

    for i = 1:nsurf
        # First term in eqn (2.30) Ramesh et al. in coefficient form
        cnc = 2*pi*(surf[i].kinem.u*cos(surf[i].kinem.alpha)/surf[i].uref + surf[i].kinem.hdot*sin(surf[i].kinem.alpha)/surf[i].uref)*(surf[i].a0[1] + surf[i].aterm[1]/2.)
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
        
        if i in shedv
            # Second term in eqn (2.30) Ramesh et al. in coefficient form
            cnnc = 2*pi*(3*surf[i].c*surf[i].a0dot[1]/(4*surf[i].uref) + surf[i].c*surf[i].adot[1]/(4*surf[i].uref) + surf[i].c*surf[i].adot[2]/(8*surf[i].uref)) + (2*lev[i]/(dt*surf[i].uref*surf[i].uref)) 
        else
            cnnc = 2*pi*(3*surf[i].c*surf[i].a0dot[1]/(4*surf[i].uref) + surf[i].c*surf[i].adot[1]/(4*surf[i].uref) + surf[i].c*surf[i].adot[2]/(8*surf[i].uref))
        end
        # Normal force coefficient
        cn = cnc + cnnc + nonl
        
        # Lift and drag coefficients
        cl[i] = cn*cos(surf[i].kinem.alpha) + cs*sin(surf[i].kinem.alpha)
        cd[i] = cn*sin(surf[i].kinem.alpha)-cs*cos(surf[i].kinem.alpha)
        
        #Pitching moment is clockwise or nose up positive
        cm[i] = cn*surf[i].pvt - 2*pi*((surf[i].kinem.u*cos(surf[i].kinem.alpha)/surf[i].uref + surf[i].kinem.hdot*sin(surf[i].kinem.alpha)/surf[i].uref)*(surf[i].a0[1]/4. + surf[i].aterm[1]/4. - surf[i].aterm[2]/8.) + (surf[i].c/surf[i].uref)*(7.*surf[i].a0dot[1]/16. + 3.*surf[i].adot[1]/16. + surf[i].adot[2]/16. - surf[i].adot[3]/64.)) - nonl_m - lev[i]*(2*surf[i].pvt - 1)/(dt*surf[i].uref*surf[i].uref)
    end
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
function write_stamp(surf :: ThreeDSurfSimple, curfield :: ThreeDFieldStrip, t :: Float64, kelv_enf :: Vector{Float64},  g:: JLD.JldGroup)
    
    cl, cd, cm, gamma, cn, cs, cnc, cnnc, nonl, cm_n, cm_pvt, nonl_m = calc_forces_more(surf.s2d)
    
    tevmat = zeros(length(curfield.tev), surf.nspan, 3)
    for i = 1:length(curfield.tev)
        for j = 1:surf.nspan
            tevmat[i,j,:] = [curfield.tev[i][j].s curfield.tev[i][j].x curfield.tev[i][j].z]
        end
    end
    g["tev"] = tevmat
    levmat = zeros(length(curfield.lev), surf.nspan, 3)
    for i = 1:length(curfield.lev)
        for j = 1:surf.nspan
            levmat[i,j,:] = [curfield.lev[i][j].s curfield.lev[i][j].x curfield.lev[i][j].z]
        end
    end
    g["lev"] = levmat
    bvmat = zeros(length(surf.s2d[1].bv), surf.nspan, 3)
    for i = 1:length(surf.s2d[1].bv)
        for j = 1:surf.nspan
            bvmat[i,j,:] = [surf.s2d[j].bv[i].s surf.s2d[j].bv[i].x surf.s2d[j].bv[i].z]
        end
    end
    g["bv"] = bvmat
    extvmat = zeros(length(curfield.extv), 3)
    for i = 1:length(curfield.extv)
        extvmat[i,:] = [curfield.extv[i].s curfield.extv[i].x curfield.extv[i].z]
    end
    g["extv"] = extvmat

    #a03d has been addded to a0 so that the calc_forces routine can be used as is
    a0 = map(q->q.a0[1], surf.s2d) - surf.a03d 
    a0dot = map(q->q.a0dot[1], surf.s2d) - surf.a03ddot
    aterm = hcat(map(q->q.aterm[:], surf.s2d)...)
    levflag = map(q->q.levflag[1], surf.s2d)

    alpha = map(q->q.kinem.alpha, surf.s2d)
    h = map(q->q.kinem.h, surf.s2d)
    u = map(q->q.kinem.u, surf.s2d)
    alphadot = map(q->q.kinem.alphadot, surf.s2d)
    hdot = map(q->q.kinem.hdot, surf.s2d)
    udot = map(q->q.kinem.udot, surf.s2d)
    
    g["t"] = t
    
    g["alpha"] = alpha
    g["h"] = h
    g["u"] = u
    g["alphadot"] = alphadot
    g["hdot"] = hdot
    g["udot"] = udot
    
    g["a0"] = a0
    g["a0dot"] = a0dot
    g["a03d"] = surf.a03d
    g["a03ddot"] = surf.a03ddot
    g["aterm"] = aterm
    g["levflag"] = levflag
    g["kelv_enf"] = kelv_enf
    
    g["cl"] = cl
    g["cd"] = cd
    g["cm"] = cm
    g["bcirc"] = gamma
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

function forces_harmonic(mat :: Array{Float64,2}, ncyc :: Int, k::Float64)
    
    T = pi/k
    nsteps = length(mat[:,1])
    range = Int(round((ncyc-1)*nsteps/ncyc))+1:nsteps 
    tbyT = (mat[range,1]-mat[range[1],1])/T

    ncol = size(mat,2)
    matnew = Array(Float64,length(tbyT),ncol)
    
    matnew[:,1] = tbyT
    for i = 2:ncol
        matnew[:,i] = mat[range,i]
    end
    
    return matnew
end

# ---------------------------------------------------------------------------------------------

#function rms(x_1::Vector{Float64}, y_1::Vector{Float64}, x_2::Vector{Float64}, y_2::Vector{Float64}; par::Float64=100)
function rms(x_1::Vector{Float64}, y_1::Vector{Float64},x_2::Vector{Float64},y_2::Vector{Float64}; par=100)
	knots_1 = (x_1,) #syntax required by interpolate function
	knots_2 = (x_2,)
	itp_1 = interpolate(knots_1, y_1, Gridded(Linear())) 
	itp_2 = interpolate(knots_2, y_2, Gridded(Linear()))

	x_max_rms = min(x_1[end],x_2[end])
	x_min_rms = max(x_1[1],x_2[1])
	x_rms = linspace(x_min_rms, x_max_rms, par)
	diff = itp_2[x_rms].-itp_1[x_rms]

	s = 0.0
	for a in diff
		s +=a*a
	end

	return sqrt(s/length(diff))
end