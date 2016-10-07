#workspace()
#include("../src/UNSflow.jl")
#using UNSflow

#This repeoduces results from Case 3A in Ramesh et al. (2014)
alpha_mean = 4*pi/180.
alpha_amp = 19.9*pi/180.
k = 0.393
h_amp = 0.5
phi = 69.8*pi/180

w = 2*k #Frequency

T = (2*pi/w) #Period
ncyc = 4
t_tot = ncyc*T  #Desired total time
dt = 0.015*0.2/(k*h_amp) #time step calculation 

nsteps =round(Int,t_tot/dt)+1

# alphadef = CosDef(alpha_mean, alpha_amp, w, phi)
# hdef = CosDef(0.,h_amp,w,0.)
# udef = ConstDef(1.)
# full_kinem = KinemDef(alphadef, hdef, udef)

# pvt = 0.25 #leading edge

# lespcrit = [0.21;]

# surf = TwoDSurf(1., 1., "sd7003_fine.dat", pvt, 70, 35, "Prescribed", full_kinem, lespcrit)

# curfield = TwoDFlowField()

# surf, curfield = ldvm(surf, curfield, nsteps, dt)

# data = readdlm("results.dat")
# writedlm("uns_jfm3c.dat",data)

# # #Present plots
data = readdlm("uns_jfm3c.dat")
exp = readdlm("exp_mg_3e.csv")
comp = readdlm("cl_mg_3e.csv")

range = round(Int,(ncyc-1)*nsteps/ncyc)+1:nsteps

tbyT = (data[range,1]-data[range[1]])/T

# #Tranfer cm from LE to QC
# cm = transfer_cm(0.25,data[:,8],data[:,6],data[:,7],data[:,2],0.,1.)

f = figure("pyplot_subplot_mixed", figsize=(15,5))

ax = subplot(121)
plot(tbyT, data[range,5], "k-")
xlabel("t*")
ylabel("LESP")
xlim(0.0,1.0)
ylim(-0.3,0.3)
xticks(0:0.2:1)
yticks(-0.3:0.1:0.3)

ax2 = ax[:twinx]() # Create another axis on top of the current axis
plot(tbyT, data[range,2]*180/pi, "r--")
font2 =Dict("color"=>"red")
ylabel("h/c",fontdict=font2)
ylim(-60,60)
yticks(-60:20:60, color="r")

ax = subplot(122)
plot(tbyT, data[range,6], "k-")
plot(comp[:,1],comp[:,2],"b--")
plot(exp[:,1],exp[:,2],"g-.")
xlabel("t*")
ylabel("Cl", color="k")
xlim(0.0,1.0)
ylim(-0.25,1.25)
xticks(0:0.2:1)
yticks(-0.25:0.25:1.25)

ax2 = ax[:twinx]() # Create another axis on top of the current axis
plot(tbyT, data[range,2]*180/pi, "r--")
font2 = Dict("color"=>"red")
ylabel("h/c",fontdict=font2)
ylim(-60,60)
yticks(-60:20:60, color="r")

# ax = subplot(223)
# plot(data[:,1], data[:,7], "k-")
# plot(comp[:,1]+0.1,comp[:,6],"b--")
# plot(exp[:,1],exp[:,4]*2,"g-.")
# xlabel("t*")
# ylabel("Cd", color="k")
# xlim(0.0,7.0)
# ylim(-0.5,1.5)
# xticks(0:1:7)
# yticks(-0.5:0.5:1.5)

# ax2 = ax[:twinx]() # Create another axis on top of the current axis
# plot(data[:,1], data[:,2]*180/pi, "r--")
# font2 = Dict("color"=>"red")
# ylabel("alpha",fontdict=font2)
# ylim(-10,30)
# yticks(-10:10:30, color="r")

# ax = subplot(224)
# plot(data[:,1], cm, "k-")
# plot(comp[:,1]+0.1,-comp[:,7],"b--")
# plot(exp[:,1],-exp[:,5]*2,"g-.")
# xlabel("t*")
# ylabel("Cm", color="k")
# xlim(0.0,7.0)
# ylim(-0.8,0.8)
# xticks(0:1:7)
# yticks(-0.8:0.4:0.8)

# ax2 = ax[:twinx]() # Create another axis on top of the current axis
# plot(data[:,1], data[:,2]*180/pi, "r--")
# font2 = Dict("color"=>"red")
# ylabel("alpha",fontdict=font2)
# ylim(-10,30)
# yticks(-10:10:30, color="r")
