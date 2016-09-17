#workspace()
include("../src/UNSflow.jl")
using UNSflow

#Base motion is an Eldredge pitch-up of SD7003 at Re=20k with amp=30 deg, K=0.2, smoothing =0.8
#Plunge is superimposed so that LEV formation is prevented
#Plunge is an integrated Eldredge pitch-up (Derivative of plunge is an Eldrege pitch-up)




alphadef = EldUpDef(30*pi/180,0.2,0.8)
h_amp = design_solve(alphadef)

hdef = EldUpIntDef(h_amp,alphadef.K*h_amp/alphadef.amp,alphadef.a)

dt = 0.015
nsteps =round(Int,3.5/dt)+1

mat = zeros(nsteps,4)
for i = 1:nsteps
    t = (i-1)*dt
    mat[i,1] = t
    mat[i,2] = alphadef(t)*180/pi
    mat[i,3] = hdef(t)
    mat[i,4] = 1.
end
