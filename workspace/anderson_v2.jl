workspace()
include("../src/UNSflow.jl")
using UNSflow

# Reproduces results Fig 6. in Anderson et al.

c = 1. 
u_ref = 1.
h_amp = 0.75
theta_amp = 5*pi/180.
alpha_max = 5*pi/180.
phi = 90.*pi/180
alpha_mean = 0.
A = 2*h_amp
pvt = 1./3. 

cp = zeros(10)
ct = zeros(10)
eta = zeros(10)
st_le = zeros(10)

for j = 1:10

st = (j-0.99)*0.6/10.
f = st*u_ref/A
w = 2*pi*f
k = w/2.
T = (2*pi/w)
ncyc = 1.
t_tot = ncyc*T
dt = 0.015*0.2*4/(k*h_amp)

nsteps =round(Int,t_tot/dt)+1
range = round(Int,(ncyc-1)*nsteps/ncyc)+1:nsteps

hdot = zeros(length(range))
alphadot = zeros(length(range))

hdef = SinDef(0.,h_amp,w,0.)
alphadef = SinDef(alpha_mean, theta_amp, w, phi)
udef = ConstDef(1.)
full_kinem = KinemDef(alphadef, hdef, udef)

lespcrit = [2;]

surf = TwoDSurf(c, u_ref, "FlatPlate", pvt, 70, 35, "Prescribed", full_kinem, lespcrit)

curfield = TwoDFlowField()

ldvm(surf, curfield, nsteps, dt)

data = readdlm("results.dat")
writedlm("and_test.dat",data)

data = readdlm("and_test.dat")

tbyT = (data[range,1]-data[range[1]])/T

    
for i = 1:length(range)
    hdot[i] = ForwardDiff.derivative(hdef,data[range[i]])*u_ref
    alphadot[i] = ForwardDiff.derivative(alphadef,data[range[i]])*u_ref/c
end

LESP = data[range,5]
cl = data[range,6]
cd = data[range,7]
cm = data[range,8]

ct[j] = (1./T)*trapz(-cd,tbyT)
cp[j] = (1./T)*(trapz(-cl.*hdot,tbyT)+trapz(cm.*alphadot,tbyT))

eta[j] = ct[j]/cp[j]

st_le[j] = st

global ct 
global cp 
global eta 
global st_le
global LESP
global tbyT

end

data_6anum = readdlm("Fig6aNum.csv")
data_6aexp = readdlm("Fig6aExp.csv")
data_6bnum = readdlm("Fig6bNum.csv")
data_6bexp = readdlm("Fig6bExp.csv")
data_6cnum = readdlm("Fig6cNum.csv")
data_6cexp = readdlm("Fig6cExp.csv")

f = figure("pyplot_subplot_mixed",figsize=(15,5))

ax = subplot(221)
plot(tbyT,LESP,"k-")
xlabel("t*")
ylabel("LESP")

ax = subplot(222)
scatter(st_le, ct,color="r",label="LDVM",marker="^")
scatter(data_6anum[:,1],data_6anum[:,2],color="b",label="Numerical",marker="o")
scatter(data_6aexp[:,1],data_6aexp[:,2],color="g",label="Experimental",marker="x")
xlabel("st_le")
ylabel("Ct")
legend(loc="upper left",fancybox="true")

ax = subplot(223)
scatter(st_le, cp,color="r",marker="^")
scatter(data_6bnum[:,1],data_6bnum[:,2],color="b",marker="o")
scatter(data_6bexp[:,1],data_6bexp[:,2],color="g",marker="x")
xlabel("st_le")
ylabel("Cp")

ax = subplot(224)
scatter(st_le, eta,color="r",marker="^")
scatter(data_6cnum[:,1],data_6cnum[:,2],color="b",marker="o")
scatter(data_6cexp[:,1],data_6cexp[:,2],color="g",marker="x")
xlabel("st_le")
ylabel("eta")
