#workspace()
#include("../src/UNSflow.jl")
#using UNSflow

#This repeoduces results from Case 3A in Ramesh et al. (2014)
alpha_mean = 0*pi/180.
alpha_amp = 76.33*pi/180.
k = 0.377
h_amp = 1.0
phi = 90*pi/180

w = 2*k #Frequency

T = (2*pi/w) #Period
ncyc = 3
t_tot = ncyc*T  #Desired total time
dt = 0.015*0.2/(k*alpha_amp) #time step calculation 

nsteps =round(Int,t_tot/dt)+1

# alphadef = CosDef(alpha_mean,alpha_amp,w,phi)
# hdef = CosDef(0.,h_amp,w,0.)
# udef = ConstDef(1.)
# full_kinem = KinemDef(alphadef, hdef, udef)

# pvt = 1./3. #leading edge

# lespcrit = [0.19;]

# surf = TwoDSurf(1., 1., "FlatPlate", pvt, 70, 35, "Prescribed", full_kinem, lespcrit)

# curfield = TwoDFlowField()



# surf, curfield = ldvm(surf, curfield, nsteps, dt)

#data = readdlm("results.dat")
#writedlm("uns_jfm4.dat",data)

# #Present plots
data = readdlm("uns_jfm4.dat")
comp_cl = readdlm("cl_kd.csv")
comp_cd = readdlm("cd_kd.csv")
comp_cm = readdlm("cm_kd.csv")

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
ylim(-0.36,0.36)
xticks(0:0.2:1)
yticks(-0.36:0.08:0.36)

ax2 = ax[:twinx]() # Create another axis on top of the current axis
plot(tbyT, data[range,2]*180/pi, "r--")
font2 =Dict("color"=>"red")
ylabel("alpha",fontdict=font2)
ylim(-90,90)
yticks(-90:30:90, color="r")

ax = subplot(122)
plot(tbyT, data[range,6], "k-")
plot(comp_cl[:,1],comp_cl[:,2],"b--")
xlabel("t*")
ylabel("Cl", color="k")
xlim(0.0,1.0)
ylim(-4.5,4.5)
xticks(0:0.2:1)
yticks(-4.5:1.5:4.5)

ax2 = ax[:twinx]() # Create another axis on top of the current axis
plot(tbyT, data[range,2]*180/pi, "r--")
font2 = Dict("color"=>"red")
ylabel("alpha",fontdict=font2)
ylim(-90,90)
yticks(-90:30:90, color="r")

ax = subplot(223)
plot(data[:,1], data[:,7], "k-")
plot(comp_cd[:,1],comp_cd[:,2],"b--")
xlabel("t*")
ylabel("Cd", color="k")
xlim(0.0,7.0)
ylim(-1.5,7.5)
xticks(0:1:7)
yticks(-1.5:1.5:7.5)

ax2 = ax[:twinx]() # Create another axis on top of the current axis
plot(data[:,1], data[:,2]*180/pi, "r--")
font2 = Dict("color"=>"red")
ylabel("alpha",fontdict=font2)
ylim(-90,90)
yticks(-90:30:90, color="r")

ax = subplot(224)
plot(data[:,1], cm, "k-")
plot(comp_cm[:,1],comp_cm[:,2],"b--")
xlabel("t*")
ylabel("Cm", color="k")
xlim(0.0,1.0)
ylim(-1.5,1.5)
xticks(0:0.2:1)
yticks(-1.5:0.5:1.5)

ax2 = ax[:twinx]() # Create another axis on top of the current axis
plot(data[:,1], data[:,2]*180/pi, "r--")
font2 = Dict("color"=>"red")
ylabel("alpha",fontdict=font2)
ylim(-90,90)
yticks(-90:30:90, color="r")
