# function transpCoupled(surf::TwoDSurfThickBL, curfield::TwoDFlowField, ncell::Int64, nsteps::Int64 = 300, dtstar::Float64 = 0.015, startflag = 0, writeflag = 0, writeInterval = 1000., delvort = delNone(); maxwrite = 50, nround=6)

#     # If a restart directory is provided, read in the simulation data
#     if startflag == 0
#         mat = zeros(0, 12)
#         t = 0.

#     elseif startflag == 1
#         dirvec = readdir()
#         dirresults = map(x->(v = tryparse(Float64,x); typeof(v) == Nothing ? 0.0 : v),dirvec)
#         latestTime = maximum(dirresults)
#         mat = DelimitedFiles.readdlm("resultsSummary")
#         t = mat[end,1]
#     else
#         throw("invalid start flag, should be 0 or 1")
#     end
#     mat = mat'

#     dt = dtstar*surf.c/surf.uref

#     # if writeflag is on, determine the timesteps to write at
#     if writeflag == 1
#         writeArray = Int64[]
#         tTot = nsteps*dt
#         for i = 1:maxwrite
#             tcur = writeInterval*real(i)
#             if t > tTot
#                 break
#             else
#                 push!(writeArray, Int(round(tcur/dt)))
#             end
#         end
#     end

#     vcore = 0.02*surf.c

#     int_wax = zeros(surf.ndiv)
#     int_c = zeros(surf.ndiv)
#     int_t = zeros(surf.ndiv)

#     for istep = 1:nsteps

#         t = t + dt

#         #Update kinematic parameters
#         update_kinem(surf, t)

#         #Update flow field parameters if any
#         update_externalvel(curfield, t)

#         #Update bound vortex positions
#         update_boundpos(surf, dt)

#         #Update incduced velocities on airfoil
#         update_indbound(surf, curfield)

#         #Set up the matrix problem
#         surf, xloc_tev, zloc_tev = update_thickLHS(surf, curfield, dt, vcore)

#         #Construct RHS vector
#         update_thickRHS(surf, curfield)

#         #Solve inviscid problem
#         invsoln = surf.LHS[1:surf.ndiv*2-2, 1:surf.naterm*2+1] \ surf.RHS[1:surf.ndiv*2-2]



#         #Assign the solution
#         for i = 1:surf.naterm
#             surf.aterm[i] = invsoln[i]
#             surf.bterm[i] = invsoln[i+surf.naterm]
#         end
#         tevstr = invsoln[2*surf.naterm+1]*surf.uref*surf.c
#         push!(curfield.tev, TwoDVort(xloc_tev, zloc_tev, tevstr, vcore, 0., 0.))

#         avisc = zeros(surf.naterm)
#         bvisc = zeros(surf.naterm)

#         #Update induced velocities to include effect of last shed vortex
#         update_indbound(surf, curfield)

#         res = 1.

#         qux = zeros(surf.ndiv)
#         qlx = zeros(surf.ndiv)
#         qut = zeros(surf.ndiv)
#         qlt = zeros(surf.ndiv)

#         iter_delu = zeros(surf.ndiv)
#         iter_dell = zeros(surf.ndiv)
#         iter_Eu = zeros(surf.ndiv)
#         iter_El = zeros(surf.ndiv)
#         wtu = zeros(surf.ndiv)
#         wtl = zeros(surf.ndiv)

#         #Iterate for viscous solution and interaction

#         iter = 0
#         while iter < 1

#             iter += 1

#             surf.aterm[:] = invsoln[1:surf.naterm] .+ avisc[:]
#             surf.bterm[:] = invsoln[surf.naterm+1:2*surf.naterm] .+ bvisc[:]

#             surf.qu[:], surf.ql[:]= calc_edgeVel(surf, [curfield.u[1], curfield.w[1]])

#             qux[2:end] = diff(surf.qu)./diff(surf.theta)
#             qux[1] = 2*qux[2] - qux[3]
#             qlx[2:end] = diff(surf.ql)./diff(surf.theta)
#             qlx[1] = 2*qlx[2] - qlx[3]

#             if iter == 3
#                 error("now")
#             end

#             if istep == 1
#                 surf.quprev[:] = surf.qu[:]
#                 surf.qlprev[:] = surf.ql[:]
#             end

#             qut[:] = (surf.qu[:] .- surf.quprev[:])./dt
#             qlt[:] = (surf.ql[:] .- surf.qlprev[:])./dt

#             w0 = [surf.delu surf.delu.*(surf.Eu .+ 1)]
#             w, j1 ,j2 = FVMIBLorig(w0, surf.qu, qut, qux, surf.theta, dt)

#             iter_delu = w[:,1]
#             iter_Eu = (w[:,2]./w[:,1]) .- 1.0
#             iter_dell[:] = surf.delu[:]
#             iter_El[:] = surf.Eu[:]

#             plot(surf.x, surf.delu)

#             wtu[2:end] = (1/100)*diff(surf.qu.*iter_delu)./diff(surf.x)
#             wtu[1] = wtu[2]
#             wtl[:] = -wtu[:]

#             figure()
#             plot(surf.x, wtu)
#             error("here")


#             #wtu[wtu .> 0.1] .= 0.1
#             #wtu[wtu .< -0.1] .= -0.1

#             RHStransp = zeros(surf.ndiv*2-2)

#             #Add transpiration velocity to RHS
#             for i = 2:surf.ndiv-1
#                 RHStransp[i-1] = 0.5*(sqrt(1 + (surf.cam_slope[i] + surf.thick_slope[i])^2)*wtu[i]
#                                        + sqrt(1 + (surf.cam_slope[i] - surf.thick_slope[i])^2)*wtl[i])

#                 RHStransp[surf.ndiv+i-3] = 0.5*(sqrt(1 + (surf.cam_slope[i] + surf.thick_slope[i])^2)*wtu[i]
#                                                       - sqrt(1 + (surf.cam_slope[i] - surf.thick_slope[i])^2)*wtl[i])
#             end

#             plot(surf.x[2:end-1], RHStransp[surf.ndiv-1:2*surf.ndiv-4])
#             println(RHStransp)

#             #Solve viscous problem
#             viscsoln = surf.LHS[1:surf.ndiv*2-2, 1:surf.naterm*2] \ RHStransp[1:surf.ndiv*2-2]

#             res = sqrt(sum((viscsoln[:] .- [avisc;bvisc]).^2))

#             println(res)

#             avisc[:] = viscsoln[1:surf.naterm]
#             bvisc[:] = viscsoln[surf.naterm+1:end]

#             #            println(bvisc)

#         end

#         #Assign bl
#         surf.delu = iter_delu[:]
#         surf.dell = iter_dell[:]
#         surf.Eu[:] = iter_Eu[:]
#         surf.El[:] = iter_El[:]

#         #Calculate adot
#         update_atermdot(surf, dt)

#         #Set previous values of aterm to be used for derivatives in next time step
#         surf.a0prev[1] = surf.a0[1]
#         for ia = 1:3
#             surf.aprev[ia] = surf.aterm[ia]
#         end

#         surf.quprev[:] = surf.qu[:]
#         surf.qlprev[:] = surf.ql[:]

#         #Calculate bound vortex strengths
#         update_bv_src(surf)

#         #Add effect of transpiration to sources and vortices

#         #Wake rollup
#         wakeroll(surf, curfield, dt)

#         #Force calculation
#         cnc, cnnc, cn, cs, cl, cd, int_wax, int_c, int_t = calc_forces(surf, int_wax, int_c, int_t, dt)

#         vle = surf.qu[1]

#         if vle > 0.
#             qspl = Spline1D(surf.x, surf.ql)
#             stag = try
#                 roots(qspl, maxn=1)[1]
# OA            catch
#                 0.
#             end
#         else
#             qspl = Spline1D(surf.x, surf.qu)
#             stag = try
#                 roots(qspl, maxn=1)[1]
#             catch
#                 0.
#             end
#         end

#         mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, vle,
#         cl, cd, cnc, cnnc, cn, cs, stag])

#         println("here")
#     end

#     mat = mat'



#     f = open("resultsSummary", "w")
#     Serialization.serialize(f, ["#time \t", "alpha (deg) \t", "h/c \t", "u/uref \t", "A0 \t", "Cl \t", "Cd \t", "Cm \n"])
#     DelimitedFiles.writedlm(f, mat)
#     close(f)

#     return mat, surf, curfield

# end




# function transpTogether(surf::TwoDSurfThickBL, curfield::TwoDFlowField, ncell::Int64, nsteps::Int64 = 300, dtstar::Float64 = 0.015, startflag = 0, writeflag = 0, writeInterval = 1000., delvort = delNone(); maxwrite = 50, nround=6)

#     # If a restart directory is provided, read in the simulation data
#     if startflag == 0
#         mat = zeros(0, 12)
#         t = 0.

#     elseif startflag == 1
#         dirvec = readdir()
#         dirresults = map(x->(v = tryparse(Float64,x); typeof(v) == Nothing ? 0.0 : v),dirvec)
#         latestTime = maximum(dirresults)
#         mat = DelimitedFiles.readdlm("resultsSummary")
#         t = mat[end,1]
#     else
#         throw("invalid start flag, should be 0 or 1")
#     end
#     mat = mat'

#     dt = dtstar*surf.c/surf.uref

#     Re = 1000


#     # if writeflag is on, determine the timesteps to write at
#     if writeflag == 1
#         writeArray = Int64[]
#         tTot = nsteps*dt
#         for i = 1:maxwrite
#             tcur = writeInterval*real(i)
#             if t > tTot
#                 break
#             else
#                 push!(writeArray, Int(round(tcur/dt)))
#             end
#         end
#     end

#     vcore = 0.02*surf.c

#     int_wax = zeros(surf.ndiv)
#     int_c = zeros(surf.ndiv)
#     int_t = zeros(surf.ndiv)

#     quprevf = zeros(surf.nfvm)

#     for istep = 1:nsteps

#         t = t + dt

#         #Update kinematic parameters
#         update_kinem(surf, t)

#         #Update flow field parameters if any
#         update_externalvel(curfield, t)

#         #Update bound vortex positions
#         update_boundpos(surf, dt)

#         #Update incduced velocities on airfoil
#         update_indbound(surf, curfield)

#         #Set up the matrix problem
#         surf, xloc_tev, zloc_tev = update_thickLHS(surf, curfield, dt, vcore)

#         #Construct RHS vector
#         update_thickRHS(surf, curfield)

#         res = 1.

#         RHStransp = zeros(surf.ndiv*2-2)

#         iter_delu = zeros(surf.ndiv)
#         iter_Eu = zeros(surf.ndiv)
#         wtu = zeros(surf.ndiv)
#         wtl = zeros(surf.ndiv)

#         quf = zeros(surf.nfvm)

#         #Iterate for viscous solution and interaction

#         iter = 0

#         while (iter < 2 || res > 5e-3)

#             iter += 1

#             if iter > 1
#                 pop!(curfield.tev)
#             end

#             #plot(surf.x[2:end-1], surf.RHS[surf.ndiv-1:2*surf.ndiv-4] + RHStransp[surf.ndiv-1:2*surf.ndiv-4])


#             soln = surf.LHS[1:surf.ndiv*2-2, 1:surf.naterm*2+1] \ (surf.RHS[1:surf.ndiv*2-2] + RHStransp[:])

#             res = sqrt(sum((soln[1:end-1] .- [surf.aterm; surf.bterm]).^2))

#             #Assign the solution
#             for i = 1:surf.naterm
#                 surf.aterm[i] = soln[i]
#                 surf.bterm[i] = soln[i+surf.naterm]
#             end
#             tevstr = soln[2*surf.naterm+1]*surf.uref*surf.c

#             push!(curfield.tev, TwoDVort(xloc_tev, zloc_tev, tevstr, vcore, 0., 0.))
#             #Update induced velocities to include effect of last shed vortex
#             update_indbound(surf, curfield)

#             #Dont update qu during iteration - just the inviscid value.
#             #This basically means no strong coupling.
#             if iter == 1
#                 surf.qu[:], surf.ql[:] = calc_edgeVel(surf, [curfield.u[1], curfield.w[1]])
#             end

#             #Transform problem to cOordinate along surface
#             srange = collect(range(0, stop=surf.su[end], length=surf.nfvm))
#             qInter = Spline1D(surf.su, surf.qu)
#             quf = evaluate(qInter, srange)
#             #qxs = derivative(qInter, srange)
#             qutf = zeros(surf.nfvm)
#             quxf = zeros(surf.nfvm)

#             quxf[2:end] = diff(quf)./diff(srange)
#             quxf[1] = 2*quxf[2] - quxf[3]

#             smoothEdges!(quxf, 10)

#             if istep == 1
#                 quprevf[:] = quf[:]
#             end

#             qutf[:] .= (quf[:] .- quprevf[:])./dt

#             delInter = Spline1D(surf.su, surf.delu)
#             EInter = Spline1D(surf.su, surf.Eu)
#             delf = evaluate(delInter, srange)
#             Ef = evaluate(EInter, srange)

#             w0 = [delf delf.*(Ef .+ 1)]

#             _, w0, tt = FVMIBLorig(w0, quf, qutf, quxf, srange, t-dt, t)



#             delf[:] = w0[:,1]

#             delInter = Spline1D(srange, delf)
#             Ef[:] = (w0[:,2]./w0[:,1]) .- 1.0


#             EInter = Spline1D(srange, Ef)

#             iter_delu[:] = evaluate(delInter, surf.su)
#             iter_Eu[:] = evaluate(EInter, surf.su)

#             #println(iter_delu)

#             #plot(surf.x, iter_delu)

#             smoothEdges!(iter_delu, 5)

#             wtu[2:end] = (1/sqrt(Re))*diff(surf.qu.*iter_delu)./diff(surf.x)
#             wtu[1] = wtu[2]

#             smoothEdges!(wtu, 10)

#             negind = 1000
#             for i = Int(floor(surf.ndiv/2)):surf.ndiv
#                 if wtu[i] < 0.
#                     negind = i
#                     break
#                 end
#             end
#             nsm = 15
#             for i = negind-nsm+1:surf.ndiv
#                 wtu[i] = wtu[negind-nsm] + (surf.x[i] - surf.x[negind-nsm])*(0 - wtu[negind-nsm])/(surf.x[end] - surf.x[negind-nsm])
#             end

#             # wtu[wtu .< 0] .= 0.

#             # maxind = argmax(wtu[Int(floor(surf.ndiv/2)):surf.ndiv]) + Int(floor(surf.ndiv/2)) - 1

#             # wtuSpl = Spline1D([surf.x[1:maxind]; surf.c], [wtu[1:maxind]; 0.])
#             # for i = maxind+1:surf.ndiv
#             #     wtu[i] = evaluate(wtuSpl, surf.x[i])
#             # end

#             # wtl[:] = -wtu[:]

#             if istep ==70 && iter == 1
#                 figure(1)
#                 plot(surf.x, iter_delu)
#                 figure(2)
#                 plot(srange, quxf)
#                 figure(3)
#                 plot(surf.x, wtu)
#                 error("here")
#             end

#             RHStransp[:] .= 0.

#             #Add transpiration velocity to RHS
#             for i = 2:surf.ndiv-1
#                 RHStransp[i-1] = 0.5*(wtu[i] + wtl[i])
#                 RHStransp[surf.ndiv+i-3] = 0.5*(wtu[i] - wtl[i])
#             end

#             println(istep, "   ", res)
#         end

#         quprevf[:] = quf[:]

#         #Assign bl
#         surf.delu[:] = iter_delu[:]
#         #surf.delu[end] = surf.delu[end-1]
#         #surf.Eu[end] = surf.Eu[end-1]

#         surf.Eu[:] = iter_Eu[:]

#         #Calculate adot
#         update_atermdot(surf, dt)

#         #Set previous values of aterm to be used for derivatives in next time step
#         surf.a0prev[1] = surf.a0[1]
#         for ia = 1:3
#             surf.aprev[ia] = surf.aterm[ia]
#         end

#         #Calculate bound vortex strengths
#         update_bv_src(surf)

#         #Add effect of transpiration to sources and vortices

#         #Wake rollup
#         wakeroll(surf, curfield, dt)

#         #Force calculation
#         cnc, cnnc, cn, cs, cl, cd, int_wax, int_c, int_t = calc_forces(surf, int_wax, int_c, int_t, dt)

#         vle = surf.qu[1]

#         if vle > 0.
#             qspl = Spline1D(surf.x, surf.ql)
#             stag = try
#                 roots(qspl, maxn=1)[1]
# OA            catch
#                 0.
#             end
#         else
#             qspl = Spline1D(surf.x, surf.qu)
#             stag = try
#                 roots(qspl, maxn=1)[1]
#             catch
#                 0.
#             end
#         end

#         mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, vle,
#         cl, cd, cnc, cnnc, cn, cs, stag])

#         println("here")
#     end

#     mat = mat'



#     f = open("resultsSummary", "w")
#     Serialization.serialize(f, ["#time \t", "alpha (deg) \t", "h/c \t", "u/uref \t", "A0 \t", "Cl \t", "Cd \t", "Cm \n"])
#     DelimitedFiles.writedlm(f, mat)
#     close(f)

#     return mat, surf, curfield

# end


# function transpTogetherWake(surf::TwoDSurfThickBL, curfield::TwoDFlowField, ncell::Int64, nsteps::Int64 = 300, dtstar::Float64 = 0.015, startflag = 0, writeflag = 0, writeInterval = 1000., delvort = delNone(); maxwrite = 50, nround=6)

#     # If a restart directory is provided, read in the simulation data
#     if startflag == 0
#         mat = zeros(0, 12)
#         t = 0.

#     elseif startflag == 1
#         dirvec = readdir()
#         dirresults = map(x->(v = tryparse(Float64,x); typeof(v) == Nothing ? 0.0 : v),dirvec)
#         latestTime = maximum(dirresults)
#         mat = DelimitedFiles.readdlm("resultsSummary")
#         t = mat[end,1]
#     else
#         throw("invalid start flag, should be 0 or 1")
#     end
#     mat = mat'

#     dt = dtstar*surf.c/surf.uref

#     Re = 1000


#     # if writeflag is on, determine the timesteps to write at
#     if writeflag == 1
#         writeArray = Int64[]
#         tTot = nsteps*dt
#         for i = 1:maxwrite
#             tcur = writeInterval*real(i)
#             if t > tTot
#                 break
#             else
#                 push!(writeArray, Int(round(tcur/dt)))
#             end
#         end
#     end

#     vcore = 0.02*surf.c

#     int_wax = zeros(surf.ndiv)
#     int_c = zeros(surf.ndiv)
#     int_t = zeros(surf.ndiv)

#     quprevf = zeros(2*surf.nfvm)

#     #Precalculate BL surface - extends to a length c beyond TE
#     ste = surf.su[end]
#     srange = collect(range(0, stop=ste, length=surf.nfvm))
#     for i = 1:surf.nfvm
#         scur = ste + surf.c*i/surf.nfvm
#         push!(srange, scur)
#     end

#     delF = zeros(2*surf.nfvm)
#     EF = zeros(2*surf.nfvm)

#     delInter = Spline1D(surf.su, surf.delu)
#     EInter = Spline1D(surf.su, surf.Eu)
#     delF[1:surf.nfvm] = evaluate(delInter, srange[1:surf.nfvm])
#     delF[surf.nfvm+1:2*surf.nfvm] .= delF[surf.nfvm]
#     EF[1:surf.nfvm] = evaluate(EInter, srange[1:surf.nfvm])
#     EF[surf.nfvm+1:2*surf.nfvm] .= EF[surf.nfvm]


#     for istep = 1:nsteps

#         t = t + dt

#         #Update kinematic parameters
#         update_kinem(surf, t)

#         #Update flow field parameters if any
#         update_externalvel(curfield, t)

#         #Update bound vortex positions
#         update_boundpos(surf, dt)

#         #Update incduced velocities on airfoil
#         update_indbound(surf, curfield)

#         #Set up the matrix problem
#         surf, xloc_tev, zloc_tev = update_thickLHS(surf, curfield, dt, vcore)

#         #Construct RHS vector
#         update_thickRHS(surf, curfield)

#         res = 1.

#         RHStransp = zeros(surf.ndiv*2-2)

#         iter_delu = zeros(surf.ndiv)
#         iter_Eu = zeros(surf.ndiv)
#         wtu = zeros(surf.ndiv)
#         wtl = zeros(surf.ndiv)

#         quf = zeros(2*surf.nfvm)

#         #Iterate for viscous solution and interaction

#         delF_iter = delF
#         EF_iter = EF

#         iter = 0

#         while (iter < 2 || res > 5e-3)

#             iter += 1

#             if iter > 1
#                 pop!(curfield.tev)
#             end

#             #plot(surf.x[2:end-1], surf.RHS[surf.ndiv-1:2*surf.ndiv-4] + RHStransp[surf.ndiv-1:2*surf.ndiv-4])


#             soln = surf.LHS[1:surf.ndiv*2-2, 1:surf.naterm*2+1] \ (surf.RHS[1:surf.ndiv*2-2] + RHStransp[:])

#             res = sqrt(sum((soln[1:end-1] .- [surf.aterm; surf.bterm]).^2))

#             #Assign the solution
#             for i = 1:surf.naterm
#                 surf.aterm[i] = soln[i]
#                 surf.bterm[i] = soln[i+surf.naterm]
#             end
#             tevstr = soln[2*surf.naterm+1]*surf.uref*surf.c

#             push!(curfield.tev, TwoDVort(xloc_tev, zloc_tev, tevstr, vcore, 0., 0.))
#             #Update induced velocities to include effect of last shed vortex
#             update_indbound(surf, curfield)

#             #Dont update qu during iteration - just the inviscid value.
#             #This basically means no strong coupling.
#             if iter == 1
#                 surf.qu[:], surf.ql[:] = calc_edgeVel(surf, [curfield.u[1], curfield.w[1]])
#             end

#             surf.qu[surf.qu .< 0.] .= 1e-6

#             #Transform problem to cOordinate along surface

#             qInter = Spline1D(surf.su, surf.qu)
#             quf[1:surf.nfvm] = evaluate(qInter, srange[1:surf.nfvm])

#             uinf = surf.kinem.u
#             ute = surf.qu[end]




#             if uinf  <= ute
#                 lt = (log(ute/uinf - 1.) - 1/log(0.02))/(srange[end] - srange[surf.nfvm])
#                 for i = 1:surf.nfvm
#                     quf[i+surf.nfvm] = uinf + (ute - uinf)*exp(-lt*(srange[i+surf.nfvm] - srange[surf.nfvm]))
#                 end
#             else
#                 lt = (log(1. - ute/uinf) - log(0.02))/(srange[end] - srange[surf.nfvm])
#                 for i = 1:surf.nfvm
#                     quf[i+surf.nfvm] = uinf - (uinf - ute)*exp(-lt*(srange[i+surf.nfvm] - srange[surf.nfvm]))
#                 end
#             end


#             #qxs = derivative(qInter, srange)

#             qutf = zeros(2*surf.nfvm)
#             quxf = zeros(2*surf.nfvm)

#             quxf[2:end] = diff(quf)./diff(srange)
#             quxf[1] = 2*quxf[2] - quxf[3]

#             smoothEdges!(quxf, 10)

#             if istep == 1
#                 quprevf[:] = quf[:]
#             end

#             qutf[:] .= (quf[:] .- quprevf[:])./dt

#             w0 = [delF delF.*(EF .+ 1)]

#             _, w0, tt = FVMIBLorig(w0, quf, qutf, quxf, srange, t-dt, t)


#             delF_iter[:] = w0[:,1]

#             #Smooth out the delta
#             nsm = 15
#             delSm = Spline1D([srange[1:surf.nfvm-nsm]; srange[surf.nfvm+nsm:end]], [delF_iter[1:surf.nfvm-nsm]; delF_iter[surf.nfvm+nsm:end]])
#             delF_iter[surf.nfvm-nsm+1:surf.nfvm+nsm-1] = evaluate(delSm, srange[surf.nfvm-nsm+1:surf.nfvm+nsm-1])


#             delInter = Spline1D(srange, delF_iter)
#             EF_iter[:] = (w0[:,2]./w0[:,1]) .- 1.0

#             EInter = Spline1D(srange, EF_iter)

#             surf.delu[:] = evaluate(delInter, surf.su)
#             surf.Eu[:] = evaluate(EInter, surf.su)

#             #println(iter_delu)

#             #plot(surf.x, iter_delu)

#             #smoothEdges!(iter_delu, 5)

#             wtu[2:end] = (1/sqrt(Re))*diff(surf.qu.*surf.delu)./diff(surf.x)
#             wtu[1] = wtu[2]

#             smoothEdges!(wtu, 10)

#             # negind = 1000
#             # for i = Int(floor(surf.ndiv/2)):surf.ndiv
#             #     if wtu[i] < 0.
#             #         negind = i
#             #         break
#             #     end
#             # end
#             # nsm = 15
#             # for i = negind-nsm+1:surf.ndiv
#             #     wtu[i] = wtu[negind-nsm] + (surf.x[i] - surf.x[negind-nsm])*(0 - wtu[negind-nsm])/(surf.x[end] - surf.x[negind-nsm])
#             # end

#             # wtu[wtu .< 0] .= 0.

#             # maxind = argmax(wtu[Int(floor(surf.ndiv/2)):surf.ndiv]) + Int(floor(surf.ndiv/2)) - 1

#             # wtuSpl = Spline1D([surf.x[1:maxind]; surf.c], [wtu[1:maxind]; 0.])
#             # for i = maxind+1:surf.ndiv
#             #     wtu[i] = evaluate(wtuSpl, surf.x[i])
#             # end

#             # wtl[:] = -wtu[:]

#             if istep ==2&& iter == 1
#                 figure(1)
#                 plot(srange, delF_iter)
#                 figure(2)
#                 plot(srange, quxf)
#                 figure(3)
#                 plot(surf.x, wtu)
#                 error("here")
#             end

#             RHStransp[:] .= 0.

#             #Add transpiration velocity to RHS
#             for i = 2:surf.ndiv-1
#                 RHStransp[i-1] = 0.5*(wtu[i] + wtl[i])
#                 RHStransp[surf.ndiv+i-3] = 0.5*(wtu[i] - wtl[i])
#             end

#             println(istep, "   ", res)
#         end

#         quprevf[:] = quf[:]

#         #Assign bl
#         delF[:] = delF_iter[:]
#         #surf.delu[end] = surf.delu[end-1]
#         #surf.Eu[end] = surf.Eu[end-1]

#         EF[:] = EF_iter[:]

#         #Calculate adot
#         update_atermdot(surf, dt)

#         #Set previous values of aterm to be used for derivatives in next time step
#         surf.a0prev[1] = surf.a0[1]
#         for ia = 1:3
#             surf.aprev[ia] = surf.aterm[ia]
#         end

#         #Calculate bound vortex strengths
#         update_bv_src(surf)

#         #Add effect of transpiration to sources and vortices

#         #Wake rollup
#         wakeroll(surf, curfield, dt)

#         #Force calculation
#         cnc, cnnc, cn, cs, cl, cd, int_wax, int_c, int_t = calc_forces(surf, int_wax, int_c, int_t, dt)

#         vle = surf.qu[1]

#         if vle > 0.
#             qspl = Spline1D(surf.x, surf.ql)
#             stag = try
#                 roots(qspl, maxn=1)[1]
# OA            catch
#                 0.
#             end
#         else
#             qspl = Spline1D(surf.x, surf.qu)
#             stag = try
#                 roots(qspl, maxn=1)[1]
#             catch
#                 0.
#             end
#         end

#         mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, vle,
#         cl, cd, cnc, cnnc, cn, cs, stag])

#         println("here")
#     end

#     mat = mat'



#     f = open("resultsSummary", "w")
#     Serialization.serialize(f, ["#time \t", "alpha (deg) \t", "h/c \t", "u/uref \t", "A0 \t", "Cl \t", "Cd \t", "Cm \n"])
#     DelimitedFiles.writedlm(f, mat)
#     close(f)

#     return mat, surf, curfield

# end



function transpSimul(surf::TwoDSurfThickBL, curfield::TwoDFlowField, ncell::Int64, nsteps::Int64 = 300, dtstar::Float64 = 0.015, startflag = 0, writeflag = 0, writeInterval = 1000., delvort = delNone(); maxwrite = 50, nround=6)

    # If a restart directory is provided, read in the simulation data
    if startflag == 0
        mat = zeros(0, 12)
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

    dt = dtstar*surf.c/surf.uref

    Re = 1000


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

    int_wax = zeros(surf.ndiv)
    int_c = zeros(surf.ndiv)
    int_t = zeros(surf.ndiv)

    quprev = zeros(surf.ndiv)


    for istep = 1:nsteps

        t = t + dt

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


        #Solve for all unknowns together through a newton iteration
        #start with an inviscid solution and BL from prev time step
        xinv = surf.LHS[1:surf.ndiv*2-2, 1:surf.naterm*2+1] \ surf.RHS[1:surf.ndiv*2-2]

        #Starting solution for del and E
        #srange = collect(range(0, stop=surf.su[end], length=surf.nfvm))
        #delInter = Spline1D(surf.su, surf.delu)
        #EInter = Spline1D(surf.su, surf.Eu)
        #delstart = evaluate(delInter, srange)
        #Estart = evaluate(EInter, srange)

        xinit = zeros(2*surf.naterm+1+2*surf.ndiv)

        xinit[1:2*surf.naterm+1] = xinv[:]
        xinit[2*surf.naterm+2:2*surf.naterm+surf.ndiv+1] = surf.delu[:]
        xinit[2*surf.naterm+surf.ndiv+2:end] = surf.Eu[:]


        #cache1=DiffCache(surf.theta)

        # IBLsimul!(Fvec, xvec) = transResidual!(Fvec, xvec, surf.naterm, surf.uref, surf.theta, xloc_tev, zloc_tev, vcore, surf.bnd_x_u, surf.bnd_z_u, surf.bnd_x_l, surf.bnd_z_l, surf.uind_u, surf.wind_u, surf.uind_l, surf.wind_l, [curfield.u[1]; curfield.w[1]], surf.kinem.alpha, surf.kinem.alphadot, surf.kinem.u, surf.kinem.hdot, surf.cam, surf.thick, surf.cam_slope, surf.thick_slope, surf.pvt, surf.c, surf.su, surf.delu, surf.Eu, surf.LHS, surf.RHS, quprev, dt, t, Re)

        resfn!(F, xvec) = transResidual!(F, xvec, surf.naterm, surf.uref, surf.theta, xloc_tev, zloc_tev, vcore, surf.bnd_x_u, surf.bnd_z_u, surf.bnd_x_l, surf.bnd_z_l, surf.uind_u, surf.wind_u, surf.uind_l, surf.wind_l, [curfield.u[1]; curfield.w[1]], surf.kinem.alpha, surf.kinem.alphadot, surf.kinem.u, surf.kinem.hdot, surf.cam, surf.thick, surf.cam_slope, surf.thick_slope, surf.pvt, surf.c, surf.su, surf.delu, surf.Eu, surf.LHS, surf.RHS, quprev, dt, t, Re)

        #resjac = xvec -> ForwardDiff.gradient(resfn, xvec)

        #Custom Newton iteration to solve system of equations



        #resjac(xinit)
        #println(resfn(xinit))

        soln = nlsolve(resfn!, xinit, iterations=5)
                       #, method=:newton, ftol=1e-3, xtol=1e-3)

        println("soln")


        error("here")

        #Iterate for viscous solution and interaction


            # if istep ==70 && iter == 1
            #     figure(1)
            #     plot(surf.x, iter_delu)
            #     figure(2)
            #     plot(srange, quxf)
            #     figure(3)
            #     plot(surf.x, wtu)
            #     error("here")
            # end


        #quprev[:] = qu[:]

        #Assign bl
        surf.delu[:] = iter_delu[:]
        #surf.delu[end] = surf.delu[end-1]
        #surf.Eu[end] = surf.Eu[end-1]

        surf.Eu[:] = iter_Eu[:]

        #Calculate adot
        update_atermdot(surf, dt)

        #Set previous values of aterm to be used for derivatives in next time step
        surf.a0prev[1] = surf.a0[1]
        for ia = 1:3
            surf.aprev[ia] = surf.aterm[ia]
        end

        #Calculate bound vortex strengths
        update_bv_src(surf)

        #Add effect of transpiration to sources and vortices

        #Wake rollup
        wakeroll(surf, curfield, dt)

        #Force calculation
        cnc, cnnc, cn, cs, cl, cd, int_wax, int_c, int_t = calc_forces(surf, int_wax, int_c, int_t, dt)

        vle = surf.qu[1]

        if vle > 0.
            qspl = Spline1D(surf.x, surf.ql)
            stag = try
                roots(qspl, maxn=1)[1]
OA            catch
                0.
            end
        else
            qspl = Spline1D(surf.x, surf.qu)
            stag = try
                roots(qspl, maxn=1)[1]
            catch
                0.
            end
        end

        mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, vle,
        cl, cd, cnc, cnnc, cn, cs, stag])

        println("here")
    end

    mat = mat'



    f = open("resultsSummary", "w")
    Serialization.serialize(f, ["#time \t", "alpha (deg) \t", "h/c \t", "u/uref \t", "A0 \t", "Cl \t", "Cd \t", "Cm \n"])
    DelimitedFiles.writedlm(f, mat)
    close(f)

    return mat, surf, curfield

end
