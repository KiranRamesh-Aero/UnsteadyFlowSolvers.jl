#workspace()
#include("../src/UNSflow.jl")
#using UNSflow
#
#Uncomment the lines above when running for the first time

#This repeoduces results from Case 3A in Ramesh et al. (2014)
alpha_mean = 0*pi/180.
alpha_amp = 00*pi/180.
k = 3.93
h_amp = 0.05
phi = 0*pi/180

w = 2*k #Frequency
T = (2*pi/w) #Period
n_cyc = 2
t_tot = n_cyc*T  #Desired total time
dt = 0.015*0.2*4/k #time step calculation 


#CosDef arguments are mean, amplitude, w, phase (all in radians)
alphadef = CosDef(alpha_mean, alpha_amp, w, phi) 
hdef = CosDef(0., h_amp, w, 0.)
udef = ConstDef(1.)
full_kinem = KinemDef(alphadef, hdef, udef)

pvt = 0.2

lesp_crit = [50;]

surf = TwoDSurf(1., 1., "FlatPlate", pvt, 70, 35, "Prescribed", full_kinem, lesp_crit)

curfield = TwoDFlowField()

nsteps =round(Int,t_tot/dt)+1

theodorsen(surf, nsteps, dt)

ldvm(surf, curfield, nsteps, dt)

# #Present plots
# data = readdlm("theo.dat")

# range = round(Int,((n_cyc-1)*nsteps/n_cyc):nsteps

# tbyT = (data[range,1]-data[range[1]])/period

# alfa_f = data[range,2]*180/pi

# h_v = data[range,3]

# cl = data[range,4]

# cd = data[range,5]

# cm = data[range,6]

# # f = figure("pyplot_subplot_mixed", figsize=(12,8))

# # ax = subplot(221)
# # plot(data[:,1], data[:,5], "k-")
# # xlabel("t*")
# # font1 = ["color"=>"black"]
# # ylabel("LESP", font1)
# # xlim(0.0,7.0)
# # ylim(-0.2,0.8)
# # xticks(0:1:7)
# # yticks(-0.2:0.2:0.8)

# # ax2 = ax[:twinx]() # Create another axis on top of the current axis
# # plot(data[:,1], data[:,2]*180/pi, "r--")
# # font2 = ["color"=>"red"]
# # ylabel("alpha",fontdict=font2)
# # ylim(-10,40)
# # yticks(-10:10:40, color="r")

# # ax = subplot(222)
# # plot(data[:,1], data[:,6], "k-")
# # xlabel("t*")
# # ylabel("Cl", color="k")
# # xlim(0.0,7.0)
# # ylim(-1,4)
# # xticks(0:1:7)
# # yticks(-2:1:4)

# # ax2 = ax[:twinx]() # Create another axis on top of the current axis
# # plot(data[:,1], data[:,2]*180/pi, "r--")
# # font2 = ["color"=>"red"]
# # ylabel("alpha",fontdict=font2)
# # ylim(-10,40)
# # yticks(-10:10:40, color="r")

# # ax = subplot(223)
# # plot(data[:,1], data[:,7], "k-")
# # xlabel("t*")
# # ylabel("Cd", color="k")
# # xlim(0.0,7.0)
# # ylim(-0.25,1)
# # xticks(0:1:7)
# # yticks(-0.5:0.25:1)

# # ax2 = ax[:twinx]() # Create another axis on top of the current axis
# # plot(data[:,1], data[:,2]*180/pi, "r--")
# # font2 = ["color"=>"red"]
# # ylabel("alpha",fontdict=font2)
# # ylim(-10,40)
# # yticks(-10:10:40, color="r")

# # ax = subplot(224)
# # plot(data[:,1], -data[:,8], "k-")
# # xlabel("t*")
# # ylabel("-Cm", color="k")
# # xlim(0.0,7.0)
# # ylim(-0.5,2)
# # xticks(0:1:7)
# # yticks(-0.5:0.5:2)

# # ax2 = ax[:twinx]() # Create another axis on top of the current axis
# # plot(data[:,1], data[:,2]*180/pi, "r--")
# # font2 = ["color"=>"red"]
# # ylabel("alpha",fontdict=font2)
# # ylim(-10,40)
# # yticks(-10:10:40, color="r")
