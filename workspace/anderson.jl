workspace()
include("../src/UNSflow.jl")
using UNSflow

#Reproduce the case from Anderson et al.


alpha_mean = 4*pi/180.
alpha_amp = 22.5*pi/180.
k = 3.93
h_amp = 0.05
phi = 90*pi/180

w = 2*k
T = (2*pi/w)
t_tot = 2*T
dt = 0.015*0.2/k

hdef = CosDef(0.,h_amp,w,0.)
alphadef = CosDef(alpha_mean, alpha_amp, w, phi)

udef = ConstDef(1.)
full_kinem = KinemDef(alphadef, hdef, udef)

pvt = 0.25 #leading edge

lespcrit = [0.21;]

surf = TwoDSurf(1., 1., "sd7003_fine.dat", pvt, 70, 35, "Prescribed", full_kinem, lespcrit)

curfield = TwoDFlowField()

nsteps =round(Int,t_tot/dt)+1

ldvm(surf, curfield, nsteps, dt)

data = readdlm("results.dat")
PyPlot.figure
plot(data[:,1],data[:,6])
PyPlot.axis([0, 2, -50, 50])
