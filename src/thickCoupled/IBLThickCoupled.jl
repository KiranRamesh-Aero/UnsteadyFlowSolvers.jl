using PyPlot
function IBLThickCoupled(surf::TwoDSurfThick, curfield::TwoDFlowField, ncell::Int64, nsteps::Int64 = 300, dtstar::Float64 = 0.015, startflag = 0, writeflag = 0, writeInterval = 1000., delvort = delNone(); maxwrite = 50, nround=6, wakerollup=1)
    # If a restart directory is provided, read in the simulation data
    if startflag == 0
        mat = zeros(0, 11)
        t = 0.
        tv = 0.
        del = zeros(ncell-1)
        E   = zeros(ncell-1)
        thick_orig = zeros(length(surf.thick))
        thick_orig[1:end] = surf.thick[1:end]

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

    del, E, xfvm, qu, ql, qu0, ql0 = initViscous(ncell)
    xfvm = xfvm/pi
    println("Determining Intitial step size ", dt)
    quf = zeros(length(qu))
    quf0 = zeros(length(qu))
    dt = 0.001  #initStepSize(surf, curfield, t, dt, 0, writeArray, vcore, int_wax, int_c, int_t, del, E, mat, startflag, writeflag, writeInterval, delvort)

    figure()
    #interactivePlot(del, E, zeros(length(E)-1), zeros(length(E)-1), xfvm, true)
    #interactivePlot(qu, surf.x, true)
    # time loop
    for istep = 1:nsteps
        t = t + dt
        #@printf(" Main time loop %1.3f\n", t);
        mat, surf, curfield, int_wax, int_c, int_t = lautat(surf, curfield, t, dt, istep, writeArray, vcore, int_wax, int_c, int_t, mat, startflag, writeflag, writeInterval, delvort)
        qu, ql = calc_edgeVel(surf, [curfield.u[1], curfield.w[1]])
        quf = mappingAerofoilToFVGrid(qu, surf, xfvm)

        if quf0 == zeros(ncell-1)
            println("Inside the validation check at", t )
            quf0 = quf
        end

        w0u, Uu,  Utu, Uxu = inviscidInterface(del, E, quf, quf0, dt, surf)

        if Utu == zeros(length(Utu))
            println("zero temporal derivative ", t )
        end

        w, dt, j1 ,j2 = FVMIBL(w0u, Uu, Utu, Uxu);
        del[:] = w[:,1]
        E[:] = (w[:,2]./w[:,1]) .- 1.0

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

            for i =1:length(j2)-1
               if j1[i] > 0.00015
                    println("singularity (separation) detected at t=$t, x=$(xfvm[i]), j1=$(j1[i])")
                end
            end

        viscousInviscid2!(surf, qu, del, thick_orig,  xfvm, 10000.0, true)

        #viscousInviscid!(surf, quf, del, thick_orig, xfvm, 10000.0, true)

        #interactivePlot(del, E, j1, j2 ,xfvm, true)
        #w0u[1:end] = w[1:end];
        interactivePlot(qu, Uxu, surf.x, xfvm, true)



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

mat, surf, curfield, del

end

function lautat(surf::TwoDSurfThick, curfield::TwoDFlowField, t::Float64, dt::Float64, istep::Int64, writeArray::Array{Int64}, vcore::Float64, int_wax::Array{Float64}, int_c::Array{Float64}, int_t::Array{Float64}, mat, startflag = 0, writeflag = 0, writeInterval = 1000., delvort = delNone(); maxwrite = 50, nround=6, wakerollup=1)

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

    #Now solve the matrix problem
    #soln = surf.LHS[[1:surf.ndiv*2-3;2*surf.ndiv-1], 1:surf.naterm*2+2] \ surf.RHS[[1:surf.ndiv*2-3; 2*surf.ndiv-1]]
    soln = surf.LHS[1:surf.ndiv*2-3, 1:surf.naterm*2+2] \ surf.RHS[1:surf.ndiv*2-3]

    #Assign the solution
    surf.a0[1] = soln[1]
    for i = 1:surf.naterm
        surf.aterm[i] = soln[i+1]
        surf.bterm[i] = soln[i+surf.naterm+1]
    end
    tevstr = soln[2*surf.naterm+2]*surf.uref*surf.c
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

    # write flow details if required
    if writeflag == 1
        if istep in writeArray
            dirname = "$(round(t,sigdigits=nround))"
            writeStamp(dirname, t, surf, curfield)
        end
    end

    mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, surf.a0[1],
                    cl, cd, cnc, cnnc, cn, cs])


        return  mat, surf, curfield, int_wax, int_c, int_t
end


function inviscidInterface(del::Array{Float64,1}, E::Array{Float64,1}, q::Array{Float64,1}, qu0::Array{Float64,1}, dt::Float64, surf::TwoDSurfThick)


    #del, E, F ,B = init(n-1)

    #m = length(q) -1
    U0 = q[2:end]
    #x =x[1:n]
    #U0 = U0[1:n]

    U00 = qu0[2:end]

    U0[U0.< 0.0] .= 1e-8
    U00[U00.< 0.0] .= 1e-8
    #println("finding the length",length(U00))
    #println("finding the m ", m)

    Ux =  spatialDeriv([q[1];U0])
    Ut =   zeros(length(U0)) #temporalDeriv(U0, U00, dt)

    w1 = zeros(length(del))
    w2 = zeros(length(del))
    w1[1:end] = del[1:end]
    w2[1:end] = del[1:end].*(E[1:end].+1.0)
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

    dqdx = ([q[2:end-1];q[end-1]]-[q[1];q[1:end-2]])

    dqdx[1] = 2*dqdx[2]- dqdx[3]
    dqdx[end] = 2*dqdx[end-1]- dqdx[end-2]

    dqdx = dqdx./(2*dx)
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
    qu =  zeros(m)
    ql =  zeros(m)
    qu0 = zeros(m)
    ql0 = zeros(m)

    return del, E, x, qu, ql, qu0, ql0

end

function initStepSize(surf::TwoDSurfThick, curfield::TwoDFlowField, t::Float64, dt::Float64, istep::Int64, writeArray::Array{Int64}, vcore::Float64, int_wax::Array{Float64}, int_c::Array{Float64}, int_t::Array{Float64}, del::Array{Float64,1}, E::Array{Float64,1}, mat, startflag = 0, writeflag = 0, writeInterval = 1000., delvort = delNone(); maxwrite = 50, nround=6, wakerollup=1)

    mat, surf, curfield, int_wax, int_c, int_t = lautat(surf, curfield, t, dt, istep, writeArray, vcore, int_wax, int_c, int_t, mat, startflag, writeflag, writeInterval, delvort)
    qu, ql = calc_edgeVel(surf, [curfield.u[1], curfield.w[1]])
    w0u, Uu,  Utu, Uxu = inviscidInterface(del, E, qu, qu, dt, surf)
    dt = initDt(w0u, Uu)

    return dt


end


function interactivePlot(del::Array{Float64,1}, E::Array{Float64,1}, x::Array{Float64,1}, disp::Bool)

    if(disp)

        PyPlot.clf()
        subplot(211)
        #axis([0, 1, (minimum(del)), (maximum(del))])
        plot(x[1:end-1],del)

        subplot(212)
        #axis([0, 1, (minimum(E)), (minimum(E))])
        plot(x[1:end-1],E)
        show()
        pause(0.005)

    end

end

function interactivePlot(del::Array{Float64,1}, E::Array{Float64,1}, J1::Array{Float64,1}, J2::Array{Float64,1}, x::Array{Float64,1}, disp::Bool)

    if(disp)

        PyPlot.clf()
        subplot(221)
        #axis([0, 1, (minimum(del)), (maximum(del))])
        plot(x[1:end-1],del)

        subplot(222)
        #axis([0, 1, (minimum(E)), (minimum(E))])
        plot(x[1:end-1],E)

        subplot(223)
        #axis([0, 1, (minimum(E)), (minimum(E))])
        plot(x[2:end-1],J1)

        subplot(224)
        #axis([0, 1, (minimum(E)), (minimum(E))])
        plot(x[2:end-1],J2)

        show()
        pause(0.005)

    end

end




function interactivePlot(qu::Array{Float64,1}, qut::Array{Float64,1}, x::Array{Float64,1}, xfvm::Array{Float64}, disp::Bool)

    if(disp)


        #PyPlot.clf()

       subplot(211)
       #axis([0, 1.2, (minimum(qu)-0.1), (maximum(qu)+0.1)])
       plot(x,qu)

       subplot(212)
      # axis([0, 1, (minimum(E)-0.1), (minimum(E)+0.1)])
       plot(xfvm[1:end-1], qut)
       show()
       pause(0.01)

    end


end

function interactivePlot(surf::TwoDSurfThick, disp::Bool)

    if(disp)
        #PyPlot.clf()

       #subplot(211)
       #axis([0, 1, (minimum(qu)-0.1), (maximum(qu)+0.1)])
       PyPlot.scatter(surf.x, surf.thick)

      # subplot(212)
      # axis([0, 1, (minimum(E)-0.1), (minimum(E)+0.1)])
      # plot(x[1:end-1],E)
       show()
       pause(0.01)

    end


end

function mappingAerofoilToFVGrid(qu::Array{Float64,1}, surf::TwoDSurfThick, xfvm::Array{Float64,1})


    n = length(xfvm)
    dx = pi/(n)
    qinter = Spline1D(surf.x,qu)
    qfine = evaluate(qinter, xfvm)

    #dqfinedx = (qfine[2:end] - qfine[1:end-1])./(dx)
    #dqfinedx[end] = 2*dqfinedx[end-1]- dqfinedx[end-2]

    #dqfinedx = ([qfine[2:end-1];qfine[end-1]]-[qfine[1];qfine[1:end-2]])

    #dqfinedx[1] = 2*dqfinedx[2]- dqfinedx[3]
    #dqfinedx[end] = 2*dqfinedx[end-1]- dqfinedx[end-2]

    #dqfinedx = dqfinedx./(2*dx)

    return qfine #, dqfinedx
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

    x = collect(0:ncell)*pi/(ncell)

    return x
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
    delinter = Spline1D(xfvm[1:end-1],del)
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
    delinter = Spline1D(xfvm[1:end-1],del)
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
        newThicknessInter = Spline1D(surf.x, newThickness)
        newThickness_slope =  derivative(newThicknessInter, surf.x)

        surf.thick[1:end] = newThickness[1:end]
        surf.thick_slope[1:end] = newThickness_slope[1:end]
    end

end
