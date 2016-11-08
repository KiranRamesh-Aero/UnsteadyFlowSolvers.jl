#workspace()
#include("../src/UNSflow.jl")
#using UNSflow

# Modification of Laura's code to find the best fit of critical LESP value


c = 1. 
u_ref = 1.


h_amp = 0.75
theta_amp = 5*pi/180.
phi = 90.*pi/180
alpha_mean = 0.
A = 2*h_amp
pvt = 1./3. 

###
st = 0.5
###

f = st*u_ref/A
w = 2*pi*f
k = w/2.
T = (2*pi/w)
ncyc = 3.
t_tot = ncyc*T
dt = 0.015*0.2*4/(k*h_amp)

nsteps =round(Int,t_tot/dt)+1
range = round(Int,(ncyc-1)*nsteps/ncyc)+1:nsteps

hdef = SinDef(0.,h_amp,w,0.)
alphadef = SinDef(alpha_mean, theta_amp, w, phi)
udef = ConstDef(1.)
full_kinem = KinemDef(alphadef, hdef, udef)

out = zeros(20,4)
index = 0

for lespc = 0.1:0.05:0.5
    lespcrit = [lespc;]
    index +=1
    
    surf = TwoDSurf(c, u_ref, "FlatPlate", pvt, 70, 35, "Prescribed", full_kinem, lespcrit)
    
    curfield = TwoDFlowField()
    
    ldvm(surf, curfield, nsteps, dt)
    
    data = readdlm("results.dat")
    
    tbyT = (data[range,1]-data[range[1]])/T
    

    LESP = data[range,5]
    cl = data[range,6]
    cd = data[range,7]
    cm = data[range,8]
    
    ct = (1./T)*trapz(-cd,tbyT)

    data = readdlm("6aExp_ct.csv",',')
    ig = Spline1D(data[:,1],data[:,2])
    ct_exp = evaluate(ig,st)
    
    out[index,1] = lespcrit[1]
    out[index,2] = ct
    out[index,3] = ct_exp
    out[index,4] = abs(ct-ct_exp)/ct_exp*100
end

writedlm("lesp_6a_0.5.dat",out)


