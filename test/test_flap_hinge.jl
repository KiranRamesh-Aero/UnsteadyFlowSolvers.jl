workspace()
include("../src/UNSflow.jl")
using UNSflow

#Kinematics
alpha_amp = 0.
alpha_mean = 0.
alpha_zl = 0. #Flat plate
h_amp = 0. 
k = 3.93
beta_amp = 5*pi/180.

#for UNSflow

w = 2*k 
T = (2*pi/w)
ncyc = 4
t_tot = ncyc*T 

dt = 0.015*0.2/(k*beta_amp) 
nsteps =round(Int,t_tot/dt)+1

alphadef = ConstDef(alpha_amp)
hdef = ConstDef(h_amp)
udef = ConstDef(1.)
ndef = CosDef(0., beta_amp, w, 0.)

full_kinem = KinemDefwFlap(alphadef, hdef, udef, ndef)

pvt = 0.2 #Doesnt matter for flap only

lespcrit = [21;] #high value to turn off LEV shedding

x_b = [0.75;]

surf = TwoDSurfwFlap(1., 1., "FlatPlate", pvt, 70, 35, "Prescribed", full_kinem,x_b, lespcrit)

curfield = TwoDFlowField()


mat, surf, curfield = ldvm(surf, curfield, nsteps, dt)

#for Theodorsen
theo_in = TheoDefwFlap(alpha_amp, h_amp, alpha_mean, alpha_zl, k, 0., pvt, beta_amp, x_b[1], 0.)

(t_theo, cl, cm_alpha, cm_beta) = theodorsen(theo_in)


#Plots to compare

#UNSflow
range = round(Int,(ncyc-1)*nsteps/ncyc)+1:nsteps
tbyT = (mat[range,1]-mat[range[1]])/T

plot(tbyT, mat[range,9])
plot(t_theo, real(cm_beta))
