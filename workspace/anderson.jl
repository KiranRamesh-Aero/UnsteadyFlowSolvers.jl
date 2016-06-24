# workspace()
# include("../src/UNSflow.jl")
# using UNSflow

#Reproduce the case from Anderson et al. (Pg 57 case 4)
c = 1. 
u_ref = 1.
h_amp = 0.75
theta_amp = 30*pi/180.
phi = 75.*pi/180
alpha_mean = 0

st = 0.36
A = 2*h_amp
f = st*u_ref/A
w = 2*pi*f
k = w/2.

pvt = 1./3.    

T = (2*pi/w)
ncyc = 4
t_tot = ncyc*T
dt = 0.015*0.2*4/(k*h_amp)

hdef = SinDef(0.,h_amp,w,0.)
alphadef = SinDef(alpha_mean, theta_amp, w, phi)
udef = ConstDef(1.)
full_kinem = KinemDef(alphadef, hdef, udef)

lespcrit = [0.2;]

surf = TwoDSurf(c, u_ref, "FlatPlate", pvt, 70, 35, "Prescribed", full_kinem, lespcrit)

curfield = TwoDFlowField()

nsteps =round(Int,t_tot/dt)+1

ldvm(surf, curfield, nsteps, dt)

data = readdlm("results.dat")
writedlm("and_test.dat",data)

data = readdlm("and_test.dat")

range = round(Int,(ncyc-1)*nsteps/ncyc)+1:nsteps

tbyT = (data[range,1]-data[range[1]])/T

hdot = zeros(length(range))
alphadot = zeros(length(range))
cp = zeros(length(range))
ct = zeros(length(range))
eta = zeros(length(range))
    
for i = 1:length(range)
    hdot[i] = ForwardDiff.derivative(hdef,data[range[i]])*u_ref
    alphadot[i] = ForwardDiff.derivative(alphadef,data[range[i]])*u_ref/c
end

cl = data[range,6]
cd = data[range,7]
cm = data[range,8]

ct = (1./T)*trapz(-cd,tbyT)
cp = (1./T)*trapz(cl.*hdot+cm.*alphadot,tbyT)

eta = ct/cp
                              
