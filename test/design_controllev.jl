#workspace()
#include("../src/UNSflow.jl")
#using UNSflow

#Base motion is an Eldredge pitch-up of SD7003 at Re=20k with amp=30 deg, K=0.2, smoothing =0.8
#Plunge is superimposed so that LEV formation is prevented
#Plunge is an integrated Eldredge pitch-up (Derivative of plunge is an Eldrege pitch-up)

function lesp_design_max(h_amp::Float64, alphadef::MotionDef)
    println(h_amp)
    hdef = EldUpInttstartDef(h_amp,alphadef.K*h_amp/alphadef.amp,alphadef.a, 5.)
  udef = ConstDef(1.)
  full_kinem = KinemDef(alphadef, hdef, udef)
  lespcrit = [20;]
  pvt = 0.25

  surf = TwoDSurf(1., 1., "sd7003_fine.dat", pvt, 70, 35, "Prescribed", full_kinem,lespcrit)

  curfield = TwoDFlowField()

  nsteps =round(Int,5.5/0.015)+1

  mat, surf, curfield = ldvm(surf, curfield, nsteps)

  return abs(maximum(mat[:,5]) - 0.25)
end

function design_solve(alphadef::MotionDef)
    iter_h = zeros(10)
    ld = zeros(10)
    iter_h[1] = 0.
    iter_h[2] = 0.1
    ld[1] = lesp_design_max(iter_h[1],alphadef)
    iter_max = 11
    iter = 1
    eps = 1e-08
    
    
    while (ld[iter] > eps)
        if (iter > iter_max)
            error("Iteration has failed")
        end
        iter = iter + 1
        ld[iter] = lesp_design_max(iter_h[iter],alphadef)
        println(ld[iter],ld[iter-1],iter_h[iter],iter_h[iter-1]) 
        dld = (ld[iter] - ld[iter-1])/(iter_h[iter] - iter_h[iter-1])
        iter_h[iter+1] = iter_h[iter] - ld[iter]/dld
    end
    
    return iter_h[iter]
end



alphadef = EldUptstartDef(30*pi/180,0.2,0.8, 5.)
h_amp = design_solve(alphadef)

hdef = EldUpInttstartDef(h_amp,alphadef.K*h_amp/alphadef.amp,alphadef.a, 5.)

dt = 0.015
nsteps =round(Int,7.5/dt)+1

mat = zeros(nsteps,4)
for i = 1:nsteps
    t = (i-1)*dt
    mat[i,1] = t
    mat[i,2] = alphadef(t)*180/pi
    mat[i,3] = hdef(t)
    mat[i,4] = 1.
end

# mat[:,1] = mat[:,1] + 4
# mat = mat[2:end,:]
# mat[:,1] = mat[:,1]-0.01
# din = zeros(267,4)
# din[:,1] = 0:0.015:4
# din[:,2] = 0.
# din[:,3] = 0.
# din[:,4] = 1.

# dful = [din;mat]
writedlm("motion_design2.dat",mat)
