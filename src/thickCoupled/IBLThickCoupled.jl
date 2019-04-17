using PyPlot
function IBLThickCoupled(surf::TwoDSurfThick, curfield::TwoDFlowField, ncell::Int64, interactiveMode::Bool, nsteps::Int64 = 300, dtstar::Float64 = 0.015, startflag = 0, writeflag = 0, writeInterval = 1000., delvort = delNone(); maxwrite = 50, nround=6, wakerollup=1)
    # If a restart directory is provided, read in the simulation data
    if startflag == 0
        mat = zeros(0, 12)
        t = 0.
        tv = 0.
        del = zeros(ncell-1)
        E   = zeros(ncell-1)
        thick_orig = zeros(length(surf.thick))
        thick_orig_slope = zeros(length(surf.thick_slope))
        thick_orig[1:end] = surf.thick[1:end]
        thick_orig_slope[1:end] = surf.thick_slope[1:end]
        qu = zeros(length(surf.thick))

    elseif startflag == 1
        dirvec = readdir()
        dirresults = map(x->(v = tryparse(Float64,x); typeof(v) == Nothing ? 0.0 : v),dirvec)
        latestTime = maximum(dirresults)
        mat = DelimitedFiles.readdlm("resultsSummary")
        t = mat[end,1]
    else
        throw("invalid start flag, should be 0 or 1")
    end
    mat = mat'

    dt = dtstar*surf.c/surf.uref

    # if writeflag is on, determine the timesteps to write at
    if writeflag == 1
        writeArray = Int64[]
        tTot = nsteps*dt
        for i = 1:maxwrite
            tcur = writeInterval*real(i)
            if t > tTot
                break
            else
                push!(writeArray, Int(round(tcur/dt)))
            end
        end
    end

    # initial momentum and energy shape factors

    vcore = 0.02*surf.c

    int_wax = zeros(surf.ndiv)
    int_c = zeros(surf.ndiv)
    int_t = zeros(surf.ndiv)
    #thick_org = surf.thick

    del, E, xfvm, qu, ql, qu0, ql0 = initViscous(ncell)
    xfvm = xfvm/pi
    println("Determining Intitial step size ", dt)
    quf = zeros(length(qu))
    quf0 = zeros(length(qu))
    dt = 0.0005  #initStepSize(surf, curfield, t, dt, 0, writeArray, vcore, int_wax, int_c, int_t, del, E, mat, startflag, writeflag, writeInterval, delvort)

    #figure()
    #interactivePlot(del, E, zeros(length(E)), zeros(length(E)), xfvm, true)
    #interactivePlot(surf, true)
    # time loop


    for istep = 1:nsteps
        t = t + dt
        #@printf(" Main time loop %1.3f\n", t);
        mat, surf, curfield, int_wax, int_c, int_t = lautat(surf, curfield, t, dt, istep, writeArray, vcore, int_wax, int_c, int_t, mat, startflag, writeflag, writeInterval, delvort)
        qu, ql = calc_edgeVel(surf, [curfield.u[1], curfield.w[1]])
        quf, thetafine = mappingAerofoilToFVGrid(qu, surf, xfvm)
        xfvm = thetafine
        if quf0 == zeros(ncell-1)
            println("Inside the validation check at", t )
            quf0 = quf
        end

        w0u, Uu,  Utu, Uxu = inviscidInterface(del, E, quf, quf0, dt, surf, xfvm)
        #qux = diff1(surf.x, qu)
        #Utu = zeros(length(Utu))
        #qinter = Spline1D(surf.x,qux)
        #qfine = evaluate(qinter, xfvm)

        #Uxu = qfine

        if Utu == zeros(length(Utu))
            println("zero temporal derivative ", t )
        end

        w, dt, j1 ,j2 = FVMIBL(w0u, Uu, Utu, Uxu, xfvm);
        del[:] = w[:,1]
        E[:] = (w[:,2]./w[:,1]) .- 1.0
        #dt = 0.0005
        #dt = dtv
        # the plots of del and E

        #display(plot(x/pi,[del, E], xticks = 0:0.1:1, layout=(2,1), legend = false))

        # convergence of the Eigen-value

            #p1 = plot(x[2:end-1]/pi,j1)
            #p2 = plot(x[2:end-1]/pi,j2)
            #p3 = plot(x/pi,del)
            #p4 = plot(x/pi,E)
            #display(plot(p1, p2, p3, p4, xticks = 0:0.1:1, layout=(2,2), legend = false))

        #display(plot(sep, xticks = 0:10:200, legend = false))
            #sleep(0.05

        #    for i =1:length(j2)-1
        #       if j1[i] > 0.00015
        #            println("singularity (separation) detected at t=$t, x=$(xfvm[i]), j1=$(j1[i])")
        #        end
        #    end
        viscousInviscid3!(surf, quf, xfvm, del, thick_orig, 100000.0, true)


        #viscousInviscid!(surf, quf, del, thick_orig, xfvm, 10000.0, true)
        #interactivePlot(del, E, xfvm, interactiveMode)
        interactivePlot(del, E, j1, j2, xfvm, true)
        #interactivePlot(surf, true)
        #w0u[1:end] = w[1:end];
        #interactivePlot(qu, Uxu, surf.x, xfvm, true)

        #for i=1:length(j1)
            #if (j1[i]>0.0002)
                #println("Seperation detected at ", xfvm[i]/pi, " time ", t)
                #return
                #j1max = findallmax(j1)
                #println("J1 max ", j1[j1max], "J1 max position ",xfvm[findallmax(j1)]/pi)
            #end
        #end

        if(istep == 260)
            println("Break")
        end

        @printf("viscous Time :%1.10f , viscous Time step size %1.10f, step number %1.1f \n", t, dt, istep);
        #map(x -> @sprintf("Seperation occure at :%1.10f  \n",x), xtrunc[seperation]./pi);

            #tv = tv + dtv

            #sols = w
        quf0[1:end] = quf[1:end];

    end

    mat = mat'

    f = open("resultsSummary", "w")
    Serialization.serialize(f, ["#time \t", "alpha (deg) \t", "h/c \t", "u/uref \t", "A0 \t", "Cl \t", "Cd \t", "Cm \n"])
    DelimitedFiles.writedlm(f, mat)
    close(f)

mat, surf, curfield, del, E, quf, qu,thick_orig, thick_orig_slope

end

function lautat4IBL(surf::TwoDSurfThick, curfield::TwoDFlowField, t::Float64, dt::Float64, istep::Int64, writeArray::Array{Int64}, vcore::Float64, int_wax::Array{Float64}, int_c::Array{Float64}, int_t::Array{Float64}, wtu::Array{Float64}, wtl::Array{Float64})

    #Update kinematic parameters
    update_kinem(surf, t)

    #Update flow field parameters if any
    update_externalvel(curfield, t)

    #Update bound vortex positions
    update_boundpos(surf, dt)

    #Update incduced velocities on airfoil
    update_indbound(surf, curfield)

    #Set up the matrix problem
    surf, xloc_tev, zloc_tev = update_thickLHS(surf, curfield, dt, vcore)

    #Construct RHS vector
    update_thickRHS(surf, curfield)

    #Add transpiration velocity to RHS
    for i = 2:surf.ndiv-1
        surf.RHS[i-1] += 0.5*(sqrt(1 + (cam_slope[i] + thick_slope[i])^2)*wtu[i]
                              + sqrt(1 + (cam_slope[i] - thick_slope[i])^2)*wtl[i])

        surf.RHS[surf.ndiv+i-3] += 0.5*(sqrt(1 + (cam_slope[i] + thick_slope[i])^2)*wtu[i]
                              - sqrt(1 + (cam_slope[i] - thick_slope[i])^2)*wtl[i])
    end

    #Now solve the matrix problem
    #soln = surf.LHS[[1:surf.ndiv*2-3;2*surf.ndiv-1], 1:surf.naterm*2+2] \ surf.RHS[[1:surf.ndiv*2-3; 2*surf.ndiv-1]]
    soln = surf.LHS[1:surf.ndiv*2-2, 1:surf.naterm*2+1] \ surf.RHS[1:surf.ndiv*2-2]

    #Assign the solution
    for i = 1:surf.naterm
        surf.aterm[i] = soln[i]
        surf.bterm[i] = soln[i+surf.naterm]
    end
    tevstr = soln[2*surf.naterm+1]*surf.uref*surf.c
    push!(curfield.tev, TwoDVort(xloc_tev, zloc_tev, tevstr, vcore, 0., 0.))

    #Calculate adot
    update_atermdot(surf, dt)

    #Set previous values of aterm to be used for derivatives in next time step
    surf.a0prev[1] = surf.a0[1]
    for ia = 1:3
        surf.aprev[ia] = surf.aterm[ia]
    end

    #Update induced velocities to include effect of last shed vortex
    update_indbound(surf, curfield)

    #Calculate bound vortex strengths
    update_bv_src(surf)

    #Wake rollup
    wakeroll(surf, curfield, dt)

    #Force calculation
    cnc, cnnc, cn, cs, cl, cd, int_wax, int_c, int_t = calc_forces(surf, int_wax, int_c, int_t, dt)

    # qu, ql = calc_edgeVel(surf, [curfield.u[1]; curfield.w[1]])

    # vle = qu[1]

    # if vle > 0.
    #     qspl = Spline1D(surf.x, ql)
    #     stag = try
    #         roots(qspl, maxn=1)[1]
    #     catch
    #         0.
    #     end
    # else
    #     qspl = Spline1D(surf.x, qu)
    #     stag = try
    #         roots(qspl, maxn=1)[1]
    #     catch
    #         0.
    #     end
    # end

    # mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, vle,
    #                 cl, cd, cnc, cnnc, cn, cs, stag])


        return  mat, surf, curfield, int_wax, int_c, int_t
end


function inviscidInterface(del::Array{Float64,1}, E::Array{Float64,1}, q::Array{Float64,1}, qu0::Array{Float64,1}, dt::Float64, surf::TwoDSurfThick, xfvm::Array{Float64,1})


    #del, E, F ,B = init(n-1)

    #m = length(q) -1
    U0 = q[1:end]
    #x =x[1:n]
    #U0 = U0[1:n]

    U00 = qu0[1:end]

    #What is this?
    U0[U0.< 0.0] .= 1e-8
    U00[U00.< 0.0] .= 1e-8
    #println("finding the length",length(U00))
    #println("finding the m ", m)

    UxInter = Spline1D(xfvm, U0)
    Ux = derivative(UxInter, xfvm)  #diff1(xfvm, U0)
    Ut = temporalDeriv(U0, U00, dt)


    w1 = zeros(length(del))
    w2 = zeros(length(del))
    w1[1:end] = del[1:end]
    w2[1:end] = del[1:end].*(E[1:end] .+ 1.0)
    w0 = hcat(w1,w2)

    return w0, U0, Ut, Ux
end

function spatialDeriv(q::Array{Float64,1})

    #n = length(qu)
    #dx = 1.0/n

    #dudx = (qu[2:end]-qu[1:end-1])./dx

    #return dudx

    n = length(q)
    dx = 1/n

    #dqdx = ([q[2:end-1];q[end-1]]-[q[1];q[1:end-2]])

    #dqdx[1] = 2*dqdx[2]- dqdx[3]
    #dqdx[end] = 2*dqdx[end-1]- dqdx[end-2]

    #dqdx = ([q[2:end]; q[1]] - [q[end];q[1:end-1]])
    #dqdx[1] = 2*dqdx[2]- dqdx[3]
    #dqdx[end] = 2*dqdx[end-1]- dqdx[end-2]
    #dqdx = dqdx./(2*dx)

    dqdx = ([q[2:end];q[1]] - [q[1:end-1];q[end]])/(dx);
    dqdx[end] = 2*dqdx[end-1]- dqdx[end-2]

    return dqdx

end


function spatialDeriv(qu::Array{Float64,1}, surf::TwoDSurfThick)

theth = surf.theta

#dqdx = ([q[2:end-1];q[end-1]]-[q[1];q[1:end-2]])

dqdthe1 = ([qu[2:end-1];qu[end-1]] - [qu[1];qu[1:end-2]])./([theth[2:end-1];theth[end-1]] - [theth[1];theth[1:end-2]])

ydx = ([yinterPi1[2:end-1];yinterPi1[end]]-[yinterPi1[1];yinterPi1[1:end-2]])./(2*n)
dthedx1 = 2/(surf.c*sin.(theth))

dthedx1 = dthedx1'
dqdx1 = dqdthe1.*dthedx1[1:end-1]

return dqdx1

end

function temporalDeriv(qu::Array{Float64,1}, qu0::Array{Float64}, dt::Float64)

 return  (qu - qu0)./dt

end

function initViscous(ncell::Int64)

    m = ncell - 1
    del, E, F ,B = init(m)
    x = createUniformFVM(m)
    quf =  zeros(m)
    quf0 =  zeros(m)

    return del, E, x, quf, quf0

end

function initStepSize(surf::TwoDSurfThick, curfield::TwoDFlowField, t::Float64, dt::Float64, istep::Int64, writeArray::Array{Int64}, vcore::Float64, int_wax::Array{Float64}, int_c::Array{Float64}, int_t::Array{Float64}, del::Array{Float64,1}, E::Array{Float64,1}, xfvm::Array{Float64,1} ,mat, startflag = 0, writeflag = 0, writeInterval = 1000., delvort = delNone(); maxwrite = 50, nround=6, wakerollup=1)

    mat, surf, curfield, int_wax, int_c, int_t = lautat(surf, curfield, t, dt, istep, writeArray, vcore, int_wax, int_c, int_t, mat, startflag, writeflag, writeInterval, delvort)
    qu, ql = calc_edgeVel(surf, [curfield.u[1], curfield.w[1]])
    w0u, Uu,  Utu, Uxu = inviscidInterface(del, E, qu, qu, dt, surf, zeros(length(E)))
    dt = initDt(w0u, Uu)

    return dt


end


function interactivePlot(del::Array{Float64,1}, E::Array{Float64,1}, x::Array{Float64,1}, disp::Bool)

    if(disp)

        PyPlot.clf()
        subplot(211)
        #axis([0, 1, (minimum(del)), (maximum(del))])
        plot(x/pi,del)
        grid()
        subplot(212)
        #axis([0, 1, (minimum(E)), (minimum(E))])
        plot(x/pi,E)
        show()
        pause(0.005)
        grid()
    end

end

function interactivePlot(del::Array{Float64,1}, E::Array{Float64,1}, J1::Array{Float64,1}, J2::Array{Float64,1}, x::Array{Float64,1}, disp::Bool)

    if(disp)

        #p1 = plot(x/pi, del, xlabel = "δ",xlims = (0,1.0),xticks = 0:0.1:1.0)
        #p2 = plot(x/pi, E, xlabel = "E",xlims = (0,1.0),xticks = 0:0.1:1.0)
        PyPlot.clf()
        subplot(221)
        #ax =gca()
        #ax[: set_xlim]([0,2])
        #p1 = plot(x/pi, del)
        xticks([0.61])
        plot(x/pi, del)


        subplot(222)
        #axis([0, 1, (minimum(E)), (minimum(E))])
        xticks([0.61])
        plot(x/pi, E)

        subplot(223)
        #axis([0, 1, (minimum(E)), (minimum(E))])
        xticks([0.61])
        plot(x[1:end]/pi, J1)

        subplot(224)
        #axis([0, 1, (minimum(E)), (minimum(E))])
        xticks([0.61])
        plot(x[1:end]/pi,J2)

        show()
        pause(0.005)

    end

end


function interactivePlot(qu::Array{Float64,1}, qut::Array{Float64,1}, x::Array{Float64,1}, xfvm::Array{Float64}, disp::Bool)

    if(disp)


        #PyPlot.clf()

       subplot(211)
       #axis([0, 1.2, (minimum(qu)-0.1), (maximum(qu)+0.1)])
       scatter(x,qu)

       subplot(212)
      # axis([0, 1, (minimum(E)-0.1), (minimum(E)+0.1)])
       scatter(xfvm, qut)
       show()
       pause(0.01)

    end


end

function interactivePlot(surf::TwoDSurfThick, disp::Bool)

    if(disp)
        #PyPlot.clf()

       subplot(211)
       #axis([0, 1, (minimum(qu)-0.1), (maximum(qu)+0.1)])
       scatter(surf.x, surf.thick)

       subplot(212)
      # axis([0, 1, (minimum(E)-0.1), (minimum(E)+0.1)])
       scatter(surf.x, surf.thick_slope)
       show()
       pause(0.001)

    end


end

#function interactivePlot(surf::TwoDSurfThick, disp::Bool)

    #if(disp)


        #PyPlot.clf()

       #subplot(211)
       #axis([0, 1, (minimum(qu)-0.1), (maximum(qu)+0.1)])
       #PyPlot.scatter(surf.x, surf.thick)

      # subplot(212)
      # axis([0, 1, (minimum(E)-0.1), (minimum(E)+0.1)])
      # plot(x[1:end-1],E)
      # show()
      # pause(0.01)

    #end


#end




function mappingAerofoilToFVGrid(qu::Array{Float64,1}, theta::Array{Float64,1}, xfvm::Array{Float64,1})


    n = length(xfvm)
    dx = (xfvm[end]-xfvm[1])/(n)
    quFine = smoothEdgeVelocity(qu, theta, n, 20)

    #dqfinedx = (qfine[2:end] - qfine[1:end-1])./(dx)
    #dqfinedx[end] = 2*dqfinedx[end-1]- dqfinedx[end-2]

    #dqfinedx = ([qfine[2:end-1];qfine[end-1]]-[qfine[1];qfine[1:end-2]])

    #dqfinedx[1] = 2*dqfinedx[2]- dqfinedx[3]
    #dqfinedx[end] = 2*dqfinedx[end-1]- dqfinedx[end-2]

    #dqfinedx = dqfinedx./(2*dx)

    return quFine #, dqfinedx
end


function reverseMappingAerofoilToAeroGrid(Ue::Array{Float64,1}, surf::TwoDSurfThick, del::Array{Float64,1}, xfvm::Array{Float64,1})


    #n = length(xfvm)
    #dx = pi/(n)

    qinter = Spline1D(xfvm,Ue)
    delinter = Spline1D(xfvm[1:end-1],del)
    qfine = evaluate(qinter, surf.x)
    delfine = evaluate(delinter, surf.x)

    #dqfinedx = (qfine[2:end] - qfine[1:end-1])./(dx)
    #dqfinedx[end] = 2*dqfinedx[end-1]- dqfinedx[end-2]

    #dqfinedx = ([qfine[2:end-1];qfine[end-1]]-[qfine[1];qfine[1:end-2]])

    #dqfinedx[1] = 2*dqfinedx[2]- dqfinedx[3]
    #dqfinedx[end] = 2*dqfinedx[end-1]- dqfinedx[end-2]

    #dqfinedx = dqfinedx./(2*dx)

    return qfine, delfine #, dqfinedx
end

function createUniformFVM(ncell::Int64)

    #x = collect(0:ncell)*pi/(ncell)
    xfvm = zeros(ncell)
     x = range(0,stop=pi,length=ncell)

     for i=1:length(x)
       xfvm[i] = x[i]
       end
    return xfvm
end

# direct viscous-inviscid coupling

function vII(surf::TwoDSurfThick, U::Array{Float64,1}, del::Array{Float64}, thick_orig::Array{Float64,1}, Re::Float64)

    thickness = zeros(length(thick_orig))
    thickness = thick_orig + ((U.* del)./(sqrt(Re)))

    return thickness

end


function viscousInviscid!(surf::TwoDSurfThick, U::Array{Float64,1}, del::Array{Float64}, thick_orig::Array{Float64,1},  xfvm::Array{Float64,1}, Re::Float64, mode::Bool)

    if(mode)
    thickFVM_orig = UnsteadyFlowSolvers.mappingAerofoilToFVGrid(thick_orig, surf, xfvm)
    thickFV = UnsteadyFlowSolvers.vII(surf, U, [del;del[end]], thickFVM_orig, Re)
    thickFV_slope = UnsteadyFlowSolvers.spatialDeriv(thickFV)

    thickAero, thickAero_slope = UnsteadyFlowSolvers.reverseMappingAerofoilToAeroGrid(U, surf, del, thickFV, thickFV_slope, xfvm)

        #println("interaction mode is on")
        #qua, dela = reverseMappingAerofoilToAeroGrid(U, surf, del, xfvm)
        #thicknew = vII(surf, qua, dela, thick_orig, Re)
        #thicknew = vII(surf, qua, dela, thick_orig, Re)
        surf.thick[1:end] = thickAero[1:end]

        #thickInter = Spline1D(surf.x,surf.thick)
        surf.thick_slope[1:end] = thickAero_slope[1:end]
    else
        surf.thick[:] = thick_orig[:]
    end

end

function reverseMappingAerofoilToAeroGrid(surf::TwoDSurfThick, del::Array{Float64,1}, xfvm::Array{Float64,1})


    #n = length(xfvm)
    #dx = pi/(n)

    #qinter = Spline1D(xfvm,Ue)
    delinter = Spline1D(xfvm/pi, del)
    #thickinter =  Spline1D(xfvm,thick)
    #dtdxinter =  derivative(thickinter, surf.x)

    #qfine = evaluate(qinter, surf.x)
    delfine = evaluate(delinter, surf.x)
    #thickfine =  evaluate(thickinter, surf.x)
    #dtdxfine =  evaluate(dtdxinter, surf.x)



    return delfine  #, dqfinedx
end

function reverseMappingAerofoilToAeroGrid(Ue::Array{Float64,1}, surf::TwoDSurfThick, del::Array{Float64,1}, thick::Array{Float64,1}, thick_slope::Array{Float64,1}, xfvm::Array{Float64,1})


    #n = length(xfvm)
    #dx = pi/(n)

    qinter = Spline1D(xfvm,Ue)
    delinter = Spline1D(xfvm,del)
    thickinter =  Spline1D(xfvm,thick)
    dtdxinter =  Spline1D(xfvm[1:end-1], thick_slope)

    #qfine = evaluate(qinter, surf.x)
    #delfine = evaluate(delinter, surf.x)
    thickfine =  evaluate(thickinter, surf.x)
    thick_slopeine =  evaluate(dtdxinter, surf.x)



    return thickfine, thick_slopeine   #, dqfinedx
end



function viscousInviscid2!(surf::TwoDSurfThick, qu::Array{Float64,1}, del::Array{Float64}, thick_orig::Array{Float64,1},  xfvm::Array{Float64,1}, Re::Float64, mode::Bool)

    if(mode)

        newThickness = zeros(length(thick_orig))
        newThickness_slope = zeros(length(thick_orig))

        delaero = reverseMappingAerofoilToAeroGrid(surf, del, xfvm)
        newThickness = thick_orig + (qu .* delaero)./(sqrt(Re))
        #newThickness_slope = diff1(surf.theta, newThickness)
        newThicknessInter = Spline1D(surf.x, newThickness)
        newThickness_slope =  derivative(newThicknessInter, surf.x)

        surf.thick[1:end] = newThickness[1:end]
        surf.thick_slope[1:end] = newThickness_slope[1:end]
    end

end


function diff1(x::Array{Float64,1}, f::Array{Float64,1})

    fp = zeros(length(x))
    fpp = zeros(length(x))
    dx = x[2:end] - x[1:end-1]
    df = f[2:end] - f[1:end-1]

    fpp[1:end-1] = atan.(df,dx)

    dx1 = x[2:end-1] - x[1:end-2]
    dx2 = x[3:end] - x[2:end-1]

    ang = (dx2.*fpp[1:end-2] + dx1.*fpp[2:end-1])./(dx1 .+ dx2)
    fp[2:end-1] = tan.(ang)

    fp[1] = 2.0*tan.(fpp[1])- fp[2]
    fp[end] = 2.0*tan.(fpp[end-1])- fp[end-1]

    return fp


end


function smoothEdgeVelocity(qu::Array{Float64,1}, theta::Array{Float64,1}, ncell::Int64, xCoarse::Int64)


    thetacoarse = collect(range(0, stop = pi, length = xCoarse))
    thetafine = collect(range(0, stop = pi, length = ncell))
    quCoarseInter = Spline1D(theta, qu)
    quCoarse = evaluate(quCoarseInter, thetacoarse)
    quFineInt = Spline1D(thetacoarse, quCoarse)
    quFine = evaluate(quFineInt, thetafine)

    return quFine, thetafine

end

function smoothEdgeGradient(qu::Array{Float64,1}, xfvm::Array{Float64,1}, xCoarse::Int64)


    thetacoarse = collect(range(0, stop = pi, length = xCoarse))
    #thetafine = collect(range(0, stop = pi, length = ncell))
    quCoarseInter = Spline1D(xfvm, qu)
    quCoarse = evaluate(quCoarseInter, thetacoarse)
    quFineInt = Spline1D(thetacoarse, quCoarse)
    quFine = evaluate(quFineInt, xfvm)

    return quFine

end

function viscousInviscid3!(surf::TwoDSurfThick, quf::Array{Float64,1}, xfvm::Array{Float64}, del::Array{Float64,1}, thick_orig::Array{Float64}, Re::Float64, mode::Bool)



    if (mode)
        newThickness = zeros(length(surf.x));
        newThickness_slope = zeros(length(surf.x));

        qufInter = Spline1D(xfvm, quf)
        qufAero = evaluate(qufInter, surf.theta)

        delInter = Spline1D(xfvm, del)
        delAero = evaluate(delInter, surf.theta)

        thetacoarse = collect(range(0, stop = pi, length = 20))
        thetafine = collect(range(0, stop = pi, length = length(surf.theta)))

        quCoarseInter = Spline1D(surf.theta, delAero)
        quCoarse = evaluate(quCoarseInter, thetacoarse)
        quFineInt = Spline1D(thetacoarse, quCoarse)
        delAerosmooth = evaluate(quFineInt, surf.theta)


    #delaero = UnsteadyFlowSolvers.reverseMappingAerofoilToAeroGrid(quf, surf, del, thick_orig, xfvm)
        newThickness = thick_orig + (qufAero .* delAero)/ (sqrt(Re))
        newThickness[newThickness.<0.0] .= 0.0
        #newThickness_slopeInter = Spline1D(surf.x, newThickness)

        newthickInter = Spline1D(surf.x, newThickness)

        #xcoarse = gridGen(70)

        #newThickCoarse = evaluate(newthickInter, xcoarse)
        #newThickfineInter = Spline1D(xcoarse, newThickCoarse)
        #newThickfine = evaluate(newThickfineInter, surf.x)
        #newthickInter = Spline1D(surf.x, newThickfine)
        newThickness_slope = diff1(surf.x, newThickness)#derivative(newthickInter, surf.x) #


        surf.thick[:] = newThickness[:]
        surf.thick_slope[:] = newThickness_slope[:]
    end

end


function gridGen(ncell::Int64)

    xcoarse = zeros(ncell)
    dtheta = pi/(length(xcoarse)-1)

    for i=1:length(xcoarse)-1
    xcoarse[i+1] = 0.5*(1-cos(dtheta*i))
    end

    return xcoarse
end

function findallmax(arr)
    max_positions = Vector{Int}()
    min_val = typemin(eltype(arr))
    for i in eachindex(arr)
        if arr[i] > min_val
            min_val = arr[i]
            empty!(max_positions)
            push!(max_positions, i)
        elseif arr[i] == min_val
            push!(max_positions, i)
        end
    end
    max_positions
end

#iterIBLConstruct! = function (x, res, surf::TwoDSurf, curfield::TwoDFlowField, quf::Vector{Float64}, quf0::Vector{Float64}, dt::Float64)

struct iterIBLsolve
    surf :: TwoDSurfThick
    curfield :: TwoDFlowField
    quf :: Array{Float64}
    quf0 :: Array{Float64}
    dt :: Float64
    xfvm :: Array{Float64}
    E :: Array{Float64}
end

function (iter::iterIBLsolve)(x::Array{Float64})
    
    res = zeros(iter.surf.ndiv*2-2+2*iter.surf.ndiv)

    #Assign iterands and calculate residuals
    iter.surf.aterm[:] = x[1:iter.surf.naterm]
    iter.surf.bterm[:] = x[iter.surf.naterm+1:2*iter.surf.naterm]
    iter.curfield.tev[end].s = x[2*iter.surf.naterm+1]
    delu = x[2*iter.surf.naterm+2:2*iter.surf.naterm+1+iter.surf.ndiv]
    dell = x[2*iter.surf.naterm+iter.surf.ndiv+2:2*iter.surf.naterm+1+2*iter.surf.ndiv]

    wtu = zeros(iter.surf.ndiv)
    wtl = zeros(iter.surf.ndiv)

    qu, ql = calc_edgeVel(iter.surf, [iter.curfield.u[1], iter.curfield.w[1]])

    wtu[2:end] = diff(qu.*delu)./diff(iter.surf.x)
    wtu[1] = 2*wtu[2] - wtu[3]
    wtl[:] = wtu[:]

    RHStransp = zeros(iter.surf.ndiv*2-2)

    #Add transpiration velocity to RHS
    for i = 2:iter.surf.ndiv-1
        RHStransp[i-1] += 0.5*(sqrt(1 + (iter.surf.cam_slope[i] + iter.surf.thick_slope[i])^2)*wtu[i]
                              + sqrt(1 + (iter.surf.cam_slope[i] - iter.surf.thick_slope[i])^2)*wtl[i])

        RHStransp[iter.surf.ndiv+i-3] += 0.5*(sqrt(1 + (iter.surf.cam_slope[i] + iter.surf.thick_slope[i])^2)*wtu[i]
                                        - sqrt(1 + (iter.surf.cam_slope[i] - iter.surf.thick_slope[i])^2)*wtl[i])
    end

    res[1:iter.surf.ndiv*2-2] = iter.surf.LHS[1:iter.surf.ndiv*2-2, 1:iter.surf.naterm*2+1]*x[1:2*iter.surf.naterm+1] - iter.surf.RHS[1:iter.surf.ndiv*2-2] - RHStransp[:]

    #Update induced velocities to include effect of last shed vortex
    update_indbound(iter.surf, iter.curfield)

    iter.quf[:], thetafine = mappingAerofoilToFVGrid(qu, iter.surf.theta, iter.xfvm)

    iter.xfvm[:] = thetafine

    w0u, Uu,  Utu, Uxu = inviscidInterface(delu, iter.E, iter.quf, iter.quf0, iter.dt, iter.surf, iter.xfvm)

    w, j1 ,j2 = FVMIBL(w0u, Uu, Utu, Uxu, iter.xfvm, iter.dt);
    delfvm = w[:,1]
    Efvm = (w[:,2]./w[:,1]) .- 1.0

    #Update transpiration velocity
    delInter = Spline1D(xfvm, delfvm)
    delu[:] = evaluate(delInter, iter.surf.theta)
    dell[:] = delu[:]

    EInter = Spline1D(xfvm, Efvm)
    iter.E[:] = evaluate(EInter, iter.surf.theta)
    dell[:] = delu[:]


    res[2*iter.surf.ndiv-1:end] = x[2*iter.surf.ndiv-1:end] .- [delu; dell]

    return res
end

#iterIBLsolve!(x, res) = iterIBLConstruct!(x, res, surf, curfield, quf, quf0, dt)


