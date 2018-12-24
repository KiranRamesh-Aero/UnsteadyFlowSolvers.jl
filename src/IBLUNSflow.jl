function ldvmLinIterative(surf::TwoDSurf, curfield::TwoDFlowField, dt::Float64 ,t::Float64,vcore::Float64, T1::Array{Float64,1}, T2::Array{Float64,1}, T3::Array{Float64,1}, delvort = delNone())

    # If a restart directory is provided, read in the simulation data

        #Update external flowfield
        update_externalvel(curfield, t)

        #Update kinematic parameters
        update_kinem(surf, t)

        #Update bound vortex positions
        update_boundpos(surf, dt)

        #Update induced velocities on airfoil
        update_indbound(surf, curfield)

        #Calculate downwash
        update_downwash(surf, [curfield.u[1],curfield.w[1]])

        #The integrals I1 and J1 are based on this downwash which doesnt
        #include vortives shed at current step (equal to term T1)
        #I1 has units if circulation, J1 is dimensionless
        T1[:] = surf.downwash[:]
        I1 = surf.c*simpleTrapz(T1.*(cos.(surf.theta) .- 1. ), surf.theta)
        J1 = -simpleTrapz(T1,surf.theta)/(surf.uref*pi)

        # T2 depends on recenetly shed TEV
        ntev = length(curfield.tev)

        if ntev == 0
            xloc_tev = surf.bnd_x[surf.ndiv] + 0.5*surf.kinem.u*dt
            zloc_tev = surf.bnd_z[surf.ndiv]
        else
            xloc_tev = surf.bnd_x[surf.ndiv]+(1. /3.)*(curfield.tev[ntev].x - surf.bnd_x[surf.ndiv])
            zloc_tev = surf.bnd_z[surf.ndiv]+(1. /3.)*(curfield.tev[ntev].z - surf.bnd_z[surf.ndiv])
        end

        for ib = 1:surf.ndiv
            xdist = surf.bnd_x[ib] - xloc_tev
            zdist = surf.bnd_z[ib] - zloc_tev
            distsq = xdist*xdist + zdist*zdist
            T2[ib] = (surf.cam_slope[ib]*zdist + xdist)/(2*pi*sqrt(distsq^2 + vcore^4))
        end

        #sig_prev = sum(map(q->q.s, curfield.tev)) + sum(map(q->q.s, curfield.lev))
        sig_prev = -surf.uref*surf.c*pi*(surf.a0prev[1] + surf.aprev[1]/2. )

        I2 = simpleTrapz(T2.*(cos.(surf.theta) .- 1. ), surf.theta)
        J2 = -simpleTrapz(T2, surf.theta)/(pi*surf.uref)

        tevstr = -(I1 + sig_prev)/(1 + I2)

        #Calc first 3 fourier coefficients and derivatives
        surf.a0[1] = J1 + J2*tevstr
        for ia = 1:3
            surf.aterm[ia] = 2. *(simpleTrapz(T1.*cos.(ia*surf.theta), surf.theta) + tevstr*simpleTrapz(T2.*cos.(ia*surf.theta), surf.theta))/(pi*surf.uref)
        end

        #Calculate adot
        surf.a0dot[1] = (surf.a0[1] - surf.a0prev[1])/dt
        for ia = 1:3
            surf.adot[ia] = (surf.aterm[ia]-surf.aprev[ia])/dt
        end

        #Check if LEV shedding is true
        if abs(surf.a0[1]) > surf.lespcrit[1]
            if (surf.a0[1] >= 0.)
                lesp_cond = surf.lespcrit[1]
            else
                lesp_cond = -surf.lespcrit[1]
            end

            # T3 depends on recenetly shed LEV
            nlev = length(curfield.lev)
            if surf.levflag[1] == 0
                le_vel_x = surf.kinem.u - surf.kinem.alphadot*sin(surf.kinem.alpha)*surf.pvt*surf.c + surf.uind[1]
                le_vel_z = -surf.kinem.alphadot*cos(surf.kinem.alpha)*surf.pvt*surf.c- surf.kinem.hdot + surf.wind[1]
                xloc_lev = surf.bnd_x[1] + 0.5*le_vel_x*dt
                zloc_lev = surf.bnd_z[1] + 0.5*le_vel_z*dt
            else
                xloc_lev = surf.bnd_x[1] + (1. /3.)*(curfield.lev[nlev].x - surf.bnd_x[1])
                zloc_lev = surf.bnd_z[1]+(1. /3.)*(curfield.lev[nlev].z - surf.bnd_z[1])
            end

            for ib = 1:surf.ndiv
                xdist = surf.bnd_x[ib] - xloc_lev
                zdist = surf.bnd_z[ib] - zloc_lev
                distsq = xdist*xdist + zdist*zdist
                T3[ib] = (surf.cam_slope[ib]*zdist + xdist)/(2*pi*sqrt(distsq^2 + vcore^4))
            end
            I3 = simpleTrapz(T3.*(cos.(surf.theta) .- 1. ), surf.theta)
            J3 = -simpleTrapz(T3, surf.theta)/(pi*surf.uref)

            det = J3*(I2 + 1. ) - J2*(I3 + 1. )

            tevstr = (-J3*(I1 + sig_prev) + (I3 + 1)*(J1 - lesp_cond))/det
            levstr = (J2*(I1 + sig_prev) - (I2 + 1)*(J1 - lesp_cond))/det

            #Recalculate required fourier terms
            surf.a0[1] = J1 + J2*tevstr + J3*levstr
            for ia = 1:3
                surf.aterm[ia] = 2. *(simpleTrapz(T1.*cos.(ia*surf.theta), surf.theta) +
                                      tevstr*simpleTrapz(T2.*cos.(ia*surf.theta), surf.theta) +
                                      levstr*simpleTrapz(T3.*cos.(ia*surf.theta), surf.theta))/(pi*surf.uref)
            end

            push!(curfield.tev, TwoDVort(xloc_tev, zloc_tev, tevstr, vcore, 0., 0.))
            push!(curfield.lev, TwoDVort(xloc_lev, zloc_lev, levstr, vcore, 0., 0.))

            for ia = 4:surf.naterm
                surf.aterm[ia] = 2. *(simpleTrapz(T1.*cos.(ia*surf.theta), surf.theta) +
                                      tevstr*simpleTrapz(T2.*cos.(ia*surf.theta), surf.theta) +
                                      levstr*simpleTrapz(T3.*cos.(ia*surf.theta), surf.theta))/(pi*surf.uref)
            end

            surf.levflag[1] = 1
        else
            push!(curfield.tev, TwoDVort(xloc_tev, zloc_tev, tevstr, vcore, 0., 0.))

            for ia = 4:surf.naterm
                surf.aterm[ia] = 2. *(simpleTrapz(T1.*cos.(ia*surf.theta), surf.theta) +
                                      tevstr*simpleTrapz(T2.*cos.(ia*surf.theta), surf.theta))/(pi*surf.uref)
            end

            surf.levflag[1] = 0
        end

#Set previous values of aterm to be used for derivatives in next time step
surf.a0prev[1] = surf.a0[1]
for ia = 1:3
    surf.aprev[ia] = surf.aterm[ia]
end

#Calculate bound vortex strengths
update_bv(surf)

# Delete or merge vortices if required
controlVortCount(delvort, surf.bnd_x[Int(round(surf.ndiv/2))], surf.bnd_z[Int(round(surf.ndiv/2))], curfield)

# free wake rollup
wakeroll(surf, curfield, dt)

# Calculate force and moment coefficients
cl, cd, cm = calc_forces(surf, [curfield.u[1], curfield.w[1]])

surf, curfield, cl, cd, cm

end

function matchEdgeVelocity()


end

function invisicViscousCoupledSolver(surf::TwoDSurf, curfield::TwoDFlowField, ncell::Int64, nsteps::Int64 =5000, dtstar::Float64 = 0.015, startflag = 0, writeflag = 0, writeInterval = 1000., delvort = delNone(); maxwrite = 50, nround=6)

    if startflag == 0
        mat = zeros(0, 8)
        t = 0.
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

    #dt = dtstar*surf.c/surf.uref
     dt=0.0002
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

    vcore = 0.02*surf.c

    T1 = zeros(surf.ndiv)
    T2 = zeros(surf.ndiv)
    T3 = zeros(surf.ndiv)

  # initiliase the viscous parameters
    invis, fluxSplitPar, soln, opCond = initiliaseData(ncell)


    for istep = 1:nsteps
        #Update current time
        t = t + dt
        surf, curfield, cl, cd, cm = UNSflow.ldvmLinIterative(surf, curfield, dt ,t, vcore, T1, T2,T3)
        update(surf,invis,curfield,opCond)
        soln, dt_v, dt_vv = solveBL(invis,fluxSplitPar,soln,opCond,dt,t,ncell)
        UNSflow.identifyFlowSeperation(surf,fluxSplitPar,soln,invis.ue_us,opCond.x,opCond.Re,ncell)


        if writeflag == 1
            if istep in writeArray
                dirname = "$(round(t,sigdigits=nround))"
                writeStamp(dirname, t, surf, curfield)
            end
        end

    # for writing in resultsSummary
        mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, surf.a0[1], cl, cd, cm])

    end

    println("The main time loop completed")
    mat = mat'
    f = open("resultsSummary", "w")
    Serialization.serialize(f, ["#time \t", "alpha (deg) \t", "h/c \t", "u/uref \t", "A0 \t", "Cl \t", "Cd \t", "Cm \n"])
    DelimitedFiles.writedlm(f, mat)
    close(f)

    mat, surf, curfield, soln, fluxSplitPar, invis
end

function solveBL(invis::InvisicidTransport, fluxSplitPar::FluxSplittingParameters, soln::Solutions, opCond::OperationalConditions, dt::Float64 ,t::Float64,ncell::Int64)

    flux  = zeros(2,2,ncell+2)
    fluxt  = zeros(2,2,ncell+2)
    rhst= zeros(2,ncell+2)
    rhs= zeros(2,ncell+2)
    uet = zeros(ncell+2)
    uex = zeros(ncell+2)
    ue= zeros(ncell+2)
    ue_t0 = zeros(ncell+2)
    xbl= zeros(ncell+2)

    ue[:] = invis.ue_us[:]
    ue_t0[:] = invis.ue_us_t0[:]
    xbl[:]= opCond.x[:]
    dx= opCond.dx

    if ue_t0==zeros(ncell+2)
        uet = zeros(ncell+2)
    else
        UNSflow.derivatives(ue,ue_t0, dt, ncell)
    end
    uex = UNSflow.derivatives(ue, xbl, ncell)  # correct the time derivative


    soln.sol[1,1] = 2*soln.sol[1,2] - soln.sol[1,3]
    soln.sol[1,ncell+2] = 2*soln.sol[1,ncell+1] - soln.sol[1,ncell]
    soln.sol[2,1] = 2*soln.sol[2,2] - soln.sol[2,3]
    soln.sol[2,ncell+2] = 2*soln.sol[2,ncell+1] - soln.sol[2,ncell]

    #soln.sol[:,1] = soln.sol[:,2]
    #soln.sol[1,1] =0.025
    #soln.sol[2,1] =0.0354
    #soln.sol[:,ncell+2] = soln.sol[:,ncell+1]

    UNSflow.correlateFunction(fluxSplitPar,soln.sol,ncell)
    UNSflow.calcEigenValues(soln,fluxSplitPar,ue,ncell);
    dt_v=UNSflow.calcDt(opCond, soln, ncell)
    dt_vv=dt_v
    dt_v=dt
    UNSflow.calcFluxes(fluxSplitPar, flux, ue, soln,ncell)
    rhst =UNSflow.rhs(soln,fluxSplitPar,uet,uex,ue,ncell)
    UNSflow.fluxSplittingSchema(flux,soln,rhst,dt_v,dx,ncell)

    #soln.solt[1,1] =0.025
    #soln.solt[2,1] =0.0354
    soln.solt[1,1] = 2*soln.solt[1,2] - soln.solt[1,3]
    soln.solt[2,1] = 2*soln.solt[2,2] - soln.solt[2,3]
    soln.solt[1,ncell+2] = 2*soln.solt[1,ncell+1] - soln.solt[1,ncell]
    soln.solt[2,ncell+2] = 2*soln.solt[2,ncell+1] - soln.solt[2,ncell]

    UNSflow.correlateFunction(fluxSplitPar,soln.solt,ncell);
    UNSflow.calcEigenValues(soln,fluxSplitPar,ue,ncell);
    UNSflow.calcFluxes(fluxSplitPar,fluxt,ue,soln,ncell)
    rhs =UNSflow.rhs(soln,fluxSplitPar,uet,uex,ue,ncell)
    UNSflow.fluxSplittingSchema(flux,fluxt,soln,rhs,dt_v,dx,ncell)
    UNSflow.correlateFunction(fluxSplitPar,soln.sol,ncell);


    invis.ue_us_t0[:]=ue[:]
    dt=dt_v

    return soln, dt_v, dt_vv
end


function initiliaseData(ncell::Int64)

    invisicidPara = InvisicidTransport(ncell)
    fluxSplittingPara = FluxSplittingParameters(ncell)
    soln =Solutions(ncell, fluxSplittingPara)
    opCond= OperationalConditions(0.3, 10000.0, ncell)

    return invisicidPara, fluxSplittingPara, soln, opCond

end


function identifyFlowSeperation(surf::TwoDSurf,fluxSplit::FluxSplittingParameters,soln::Solutions,ue::Array{Float64,1},x::Array{Float64,1},re::Float64, ncell::Int64)


    for i = 2:ncell+1
        #soln.cf[i] = fluxSplit.B[i]*ue[i]*surf.c/(fluxSplit.del[i]*surf.uref*re)
        soln.Csep[i] = abs(fluxSplit.del[i+1] - fluxSplit.del[i])/abs(fluxSplit.del[i] - fluxSplit.del[i-1])
        if abs(soln.Csep[i]) > 10.
            println("singularity (separation) detected at x=$(x[i])")
        end
    end

end

function identifyFlowSeperationCirc(fluxSplit::FluxSplittingParameters,soln::Solutions,ue::Array{Float64,1},x::Array{Float64,1})


    for i = 2:size(soln.sol,2)-1
        #soln.cf[i] = fluxSplit.B[i]*ue[i]*surf.c/(fluxSplit.del[i]*surf.uref*re)
        soln.Csep[i] = abs(fluxSplit.del[i+1] - fluxSplit.del[i])/abs(fluxSplit.del[i] - fluxSplit.del[i-1])
        if abs(soln.Csep[i]) > 10.
            println("singularity (separation) detected at x=$(x[i])")
        end
    end

end



function viscousDt(invisicidPara::InvisicidTransport, fluxSplittingPara::FluxSplittingParameters, soln::Solutions, opCond::OperationalConditions)

    correlateFunction(fluxSplittingPara,soln);
    calcEigenValues(soln,fluxSplittingPara,invisicidPara.ue_us);
    dt=calcDt(opCond, soln)

  return dt
 end


function update(surf::TwoDSurf, invisidPara::InvisicidTransport, curfield::TwoDFlowField, opCond::OperationalConditions)

    q_u,q_l = UNSflow.calc_edgeVel(surf, [curfield.u[1], curfield.w[1]])

    ncell= length(q_u)
    x = zeros(ncell+2)
    dx= 1.0/(convert(Float64,ncell+2))

    for i=2:ncell+1
        x[i] = (convert(Float64,i)-1.5)*dx
        invisidPara.ue_us[i] = q_u[i-1]
        invisidPara.ue_ls[i] = q_l[i-1]
    end

    x[1] =2*x[2]-x[3]
    x[ncell+2] = 2*x[ncell+1]-x[ncell]

    #x[1] =0.
    #x[ncell+2] = 1.
    #invisidPara.ue_us[1]= 2.0*invisidPara.ue_us[2]-invisidPara.ue_us[3]
    #invisidPara.ue_us[ncell+2]= 2.0*invisidPara.ue_us[ncell+1]-invisidPara.ue_us[ncell+1]
    #invisidPara.ue_us[1]= 0.09512534545251
    #invisidPara.ue_us[2]= 0.009512534545251
    #invisidPara.ue_us[ncell+2]= 0.9947793650842671

    invisidPara.ue_ls[1]= 2.0*invisidPara.ue_ls[2]-invisidPara.ue_ls[3]
    invisidPara.ue_ls[ncell+2]= 2.0*invisidPara.ue_ls[ncell+1]-invisidPara.ue_ls[ncell+1]

    opCond.dx=dx
    opCond.x[:]= x[:]

 end


 function updateCirc(ue::Array{Float64,1}, invisidPara::InvisicidTransport)

  #q_u,q_l = UNSflow.calc_edgeVel(surf, [curfield.u[1], curfield.w[1]])

   for i=1:length(ue)
   invisidPara.ue_us[i] = ue[i]
   invisidPara.ue_ls[i] = ue[i]
   end

   #boundaryCorrection(invisidPara.ue_us)
   #boundaryCorrection(invisidPara.ue_ls)

 end

function boundaryCorrection(u::Array{Float64,1})

n_bl = length(u)
#println(n_bl)
u[1]= 2*u[2]-u[3]
u[n_bl] = 2*u[n_bl-1]-u[n_bl-2]
#println(u[n_bl])
end

function boundaryCorrection(u::Array{Float64,2})

    for i =1:ndims(u)
        n_bl = length(u[i,:])
        #println(n_bl)
        u[i,1]= 2*u[i,2]-u[i,3]
        u[i,n_bl] = 2*u[i,n_bl-1]-u[i,n_bl-2]
        #println(u[n_bl])
    end
end



function rhs(sol::Solutions,fluxSp::FluxSplittingParameters,uet::Array{Float64,1},uex::Array{Float64,1},ue::Array{Float64,1},ncell::Int64)


 rhs =zeros(2,ncell+2)
 for i = 2:ncell+1
     rhs[1,i] = fluxSp.B[i]/(2*fluxSp.del[i]) - (fluxSp.del[i]*uet[i])/ue[i] - (fluxSp.E[i]+1.0)*fluxSp.del[i]*uex[i]
     rhs[2,i] = fluxSp.S[i]/fluxSp.del[i] - (2.0*fluxSp.E[i]*fluxSp.del[i]*uet[i])/ue[i] - (2.0*fluxSp.F[i]*fluxSp.del[i]*uex[i])
 end
# The boundary values do not need
boundaryCorrection(rhs[1,:])
boundaryCorrection(rhs[2,:])

return rhs
end

function fluxSplittingSchema(flux::Array{Float64,3},sol::Solutions,rhs::Array{Float64,2},dt::Float64,dx::Float64,ncell::Int64)



    for i = 2:ncell
        sol.solt[1,i] = sol.sol[1,i] - (dt/dx)*(flux[1,1,i] - flux[1,1,i-1] + flux[2,1,i+1]
        -flux[2,1,i]) + dt*rhs[1,i]
        sol.solt[2,i] = sol.sol[2,i] - (dt/dx)*(flux[1,2,i] - flux[1,2,i-1] + flux[2,2,i+1]
        -flux[2,2,i]) + dt*rhs[2,i]
     end
    i = ncell+1
        sol.solt[1,i] = sol.sol[1,i] - (dt/dx)*(flux[1,1,i] - flux[1,1,i-1]) + dt*rhs[1,i]
        sol.solt[2,i] = sol.sol[2,i] - (dt/dx)*(flux[1,2,i] - flux[1,2,i-1]) + dt*rhs[2,i]

end

function fluxSplittingSchema(flux::Array{Float64,3},fluxt::Array{Float64,3},sol::Solutions,rhs::Array{Float64,2},dt::Float64,dx::Float64, ncell::Int64)

    i = 2
    sol.sol[1,i] = 0.5*(sol.sol[1,i] + sol.solt[1,i]) - (0.5*dt/dx)*(flux[1,1,i] - flux[1,1,i-1]
    - flux[2,1,i] + 2*flux[2,1,i+1] - flux[2,1,i+2] + fluxt[1,1,i] - fluxt[1,1,i-1]
    + fluxt[2,1,i+1] - fluxt[2,1,i]) + 0.5*dt*rhs[1,i]
    sol.sol[2,i] = 0.5*(sol.sol[2,i] + sol.solt[2,i]) - (0.5*dt/dx)*(flux[1,2,i] - flux[1,2,i-1]
    - flux[2,2,i] + 2*flux[2,2,i+1]-flux[2,2,i+2]+fluxt[1,2,i] - fluxt[1,2,i-1]
    + fluxt[2,2,i+1] - fluxt[2,2,i])+0.5*dt*rhs[2,i]

    for i = 3:ncell
        sol.sol[1,i] = 0.5*(sol.sol[1,i] + sol.solt[1,i]) - (0.5*dt/dx)*(flux[1,1,i] - 2*flux[1,1,i-1] + flux[1,1,i-2]
        - flux[2,1,i] + 2*flux[2,1,i+1] - flux[2,1,i+2] + fluxt[1,1,i] - fluxt[1,1,i-1]
        + fluxt[2,1,i+1] - fluxt[2,1,i]) + 0.5*dt*rhs[1,i]
        sol.sol[2,i] = 0.5*(sol.sol[2,i] + sol.solt[2,i]) - (0.5*dt/dx)*(flux[1,2,i] - 2*flux[1,2,i-1] + flux[1,2,i-2]
        - flux[2,2,i] + 2*flux[2,2,i+1] - flux[2,2,i+2] + fluxt[1,2,i] - fluxt[1,2,i-1]
        + fluxt[2,2,i+1] - fluxt[2,2,i]) + 0.5*dt*rhs[2,i]
    end
    i = ncell+1
    sol.sol[1,i] = 0.5*(sol.sol[1,i] + sol.solt[1,i]) - (0.5*dt/dx)*(flux[1,1,i] - 2*flux[1,1,i-1] + flux[1,1,i-2]
    +fluxt[1,1,i] - fluxt[1,1,i-1]) + 0.5*dt*rhs[1,i]
    sol.sol[2,i] = 0.5*(sol.sol[2,i] + sol.solt[2,i]) - (0.5*dt/dx)*(flux[1,2,i] - 2*flux[1,2,i-1] + flux[1,2,i-2]
    +fluxt[1,2,i] - fluxt[1,2,i-1]) + 0.5*dt*rhs[2,i]

end


function derivatives(dA::Array{Float64,1}, dB::Array{Float64,1},ncell::Int64)

#n_ele = length(dA)
    der =zeros(ncell+2)
    if ((ndims(dA)==1) && (ndims(dB)==1))
        if (length(dA)==length(dB))
            for i=2:ncell+1
                der[i] = (dA[i]-dA[i-1])/(dB[i]-dB[i-1])
            end
        else
            error("Input arrays should have identical dimensions")
        end
    else
        error("Function only acccepts one-dimensional arrays")
    end

   der[1]=2*der[2]-der[3]
   der[ncell+2] = 2*der[ncell+1]-der[ncell]
   return der
end

function derivatives(dA::Array{Float64,1},dA0::Array{Float64,1}, dt::Float64, ncell::Int64)

    #n_ele = length(dA)
    der =zeros(ncell+2)
        if ((ndims(dA)==1) && (ndims(dA0)==1))
            if (length(dA)==length(dA0))
                for i=2:ncell+1
                    der[i] = (dA[i]-dA0[i])/(dt)
                end
            else
                error("Input arrays should have identical dimensions")
            end
        else
            error("Function only acccepts one-dimensional arrays")
        end

    der[1]=2*der[2]-der[3]
    der[ncell+2] = 2*der[ncell+1]-der[ncell]
    return der
end

function correlateFunction(fluxSplitPara::FluxSplittingParameters, sol::Array{Float64,2},ncell::Int64)


    for i=1:ncell+2
        fluxSplitPara.del[i] = sol[1,i]
        fluxSplitPara.E[i] = sol[2,i]/(fluxSplitPara.del[i]) - 1.0
        fluxSplitPara.F[i] = 4.8274*fluxSplitPara.E[i]^4 - 5.9816*fluxSplitPara.E[i]^3 + 4.0274*fluxSplitPara.E[i]^2 + 0.23247*fluxSplitPara.E[i] + 0.15174

        if fluxSplitPara.E[i] < -0.0616
            fluxSplitPara.B[i] = -225.86*fluxSplitPara.E[i]^3 - 3016.6*fluxSplitPara.E[i]^2 - 208.68*fluxSplitPara.E[i] - 17.915
        elseif fluxSplitPara.E[i] > -0.0395
            fluxSplitPara.B[i] = 131.9*fluxSplitPara.E[i]^3 - 167.32*fluxSplitPara.E[i]^2 + 76.642*fluxSplitPara.E[i] - 11.068
        else
            fluxSplitPara.B[i] = 0.5*(-225.86*fluxSplitPara.E[i]^3 - 3016.6*fluxSplitPara.E[i]^2 - 208.68*fluxSplitPara.E[i] - 17.915
            + 131.9*fluxSplitPara.E[i]^3 - 167.32*fluxSplitPara.E[i]^2 + 76.642*fluxSplitPara.E[i] - 11.068)
        end
        if fluxSplitPara.E[i] < -0.0582
            fluxSplitPara.S[i] = 451.55*fluxSplitPara.E[i]^3 + 2010*fluxSplitPara.E[i]^2 + 138.96*fluxSplitPara.E[i] + 11.296
        elseif fluxSplitPara.E[i] > -0.042
            fluxSplitPara.S[i] = -96.739*fluxSplitPara.E[i]^3 + 117.74*fluxSplitPara.E[i]^2 - 46.432*fluxSplitPara.E[i] + 6.8074
        else
            fluxSplitPara.S[i] = 0.5*(451.55*fluxSplitPara.E[i]^3 + 2010*fluxSplitPara.E[i]^2 + 138.96*fluxSplitPara.E[i] + 11.296
            - 96.739*fluxSplitPara.E[i]^3 + 117.74*fluxSplitPara.E[i]^2 - 46.432*fluxSplitPara.E[i] + 6.8074)
        end

        fluxSplitPara.dfde[i] = 4*4.8274*fluxSplitPara.E[i]^3 - 3*5.9816*fluxSplitPara.E[i]^2 + 2*4.0274*fluxSplitPara.E[i] + 0.23247
    end

end




#The description according to Numerical Recipies in Fortrant 90 :
# Given arrays x and y of length N containing a tabulated function, i.e., yi = f(xi), with x1 <
# x2 < ... < xN , and given values yp1 and ypn for the first derivative of the interpolating
# function at points 1 and N, respectively, this routine returns an array y2 of length N
# that contains the second derivatives of the interpolating function at the tabulated points
# xi. If yp1 and/or ypn are equal to 1 × 1030 or larger, the routine is signaled to set the
# corresponding boundary condition for a natural spline, with zero second derivative on that
# boundary

function spline(x::Vector{Float64}, y::Vector{Float64}, n::Int64, ypn::Float64, y2::Vector{Float64})

nmax::Int64 =3000
p,qn,sig,un :: Float64
u :: Vector{Float64,nmax}
tol =99e30

if (yp1>tol)
y[1]=0
u[1]=0
else
  y2[1]=-0.5
  u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1)
end

for i=2:n-1
sig= (x[i]-x[i-1])/(x[i+1]-x[i-1])
p=sig*y2[i-1]+2
y2[i]=(sig-1.)/p
u[i]=(6.0*((y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]))/(x[i+1]-x[i-1])-sig*u[i-1])/p
end

if(ypn>tol)
  qn=0.
  un=0.
else
  qn=0.5
  un = (3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]))
end
y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0)

for j in n-1:-1:1
 y2[k]=y2[k]*y2[k+1]+u[k]
end

return y2

end

#The following description according to Numerical Recipies in Fortrant 90 :
#Given the arrays xa and ya, which tabulate a function (with the xai ’s in increasing or
#decreasing order), and given the array y2a, which is the output from spline above, and
#given a value of x, this routine returns a cubic-spline interpolated value. The arrays xa, ya
# and y2a are all of the same size.


function splint(xa::Vector{Float64},ya::Vector{Float64},y2a::Vector{Float64},n::Int64,x::Float64,y::Float64)

khi::Int64=n
klo::Int64=1
k::Int64
a,b,h::Float64

klo=max(min(locate(xa,x),n-1),1)

khi=klo+1
h= xa[khi]-xa[klo]

if(h===0.0)
  error("h must be distinct")
end
  a=(xa[khi]-x)/h
  b=(x-xa[klo])/h
  y=a*ya[klo]+b*ya[khi]+((a^3-a)*y2a[klo]+(b^3-b)*y2a[khi])*(h^2)/6.0_sp

  return y

end

# The description according to the Numerical Recipies in Fortrant 90
# Given an array xx(1:N), and given a value x, returns a value j such that x is between
# xx(j) and xx(j + 1). xx must be monotonic, either increasing or decreasing. j = 0 or
# j = N is returned to indicate that x is out of range.


function locate(xa::Vector{Float64}, x::Float64)

locate::Int64
ascend::Bool
n,jl,jm,ju::Int64

n=length(xa)[1]

ascend(xa[n]>=xa[1])

jl=0
ju=u+1

while ((ju-jl)<=1)

jm=(jl+ju)/2   #replace the integer devision with dev(jl+ju,2)

if(ascend===(x>=xa[jm]))
jl=jm
else
  ju=jm
end
end

if (x === xx[1]) #then Then set the output, being careful with the endpoints.
locate=1
elseif (x === xx[n]) then
locate=n-1
else
locate=jl
end

return locate
end
