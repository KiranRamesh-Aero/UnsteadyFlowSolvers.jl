workspace()
include("../src/UNSflow.jl")
using UNSflow
#
#Uncomment the lines above when running for the first time

#This repeoduces results from Case 3A in Ramesh et al. (2014)
alpha_mean = 4*pi/180.
alpha_amp = 22.5*pi/180.
k = 3.93
h_amp = 0.05
phi = 90*pi/180

w = 2*k #Frequency
T = (2*pi/w) #Period
t_tot = 2*T  #Desired total time
dt = 0.015*0.2*4/k #time step calculation 


#CosDef arguments are mean, amplitude, k (red. freq.), phase (all in radians)
alphadef = CosDef(4*pi/180, 22.5*pi/180, 3.93, 90*pi/180) 
hdef = CosDef(0., 0.05, 3.93, 0.)
udef = ConstDef(1.)
full_kinem = KinemDef(alphadef, hdef, udef)

pvt = 0.25

lespcrit = [0.21;]

pvt = 0.0 #leading edge

surf = TwoDSurf(1., 1., "sd7003_fine.dat", pvt, 70, 35, "Prescribed", full_kinem, lespcrit)

curfield = TwoDFlowField()

nsteps =round(Int,t_tot/dt)+1

ldvm(surf, curfield, nsteps, dt)

#Present plots
data = readdlm("results.dat")

range = round(Int,2*nsteps/3):nsteps

tbyT = (data[range,1]-data[range[1]])/period

alfa_f = data[range,2]*180/pi

h_v = data[range,3]

lesp = data[range,5]

cl = data[range,5]

cd = data[range,6]

cm = data[range,7]

# f = figure("pyplot_subplot_mixed", figsize=(12,8))

# ax = subplot(221)
# plot(data[:,1], data[:,5], "k-")
# xlabel("t*")
# font1 = ["color"=>"black"]
# ylabel("LESP", font1)
# xlim(0.0,7.0)
# ylim(-0.2,0.8)
# xticks(0:1:7)
# yticks(-0.2:0.2:0.8)

# ax2 = ax[:twinx]() # Create another axis on top of the current axis
# plot(data[:,1], data[:,2]*180/pi, "r--")
# font2 = ["color"=>"red"]
# ylabel("alpha",fontdict=font2)
# ylim(-10,40)
# yticks(-10:10:40, color="r")

# ax = subplot(222)
# plot(data[:,1], data[:,6], "k-")
# xlabel("t*")
# ylabel("Cl", color="k")
# xlim(0.0,7.0)
# ylim(-1,4)
# xticks(0:1:7)
# yticks(-2:1:4)

# ax2 = ax[:twinx]() # Create another axis on top of the current axis
# plot(data[:,1], data[:,2]*180/pi, "r--")
# font2 = ["color"=>"red"]
# ylabel("alpha",fontdict=font2)
# ylim(-10,40)
# yticks(-10:10:40, color="r")

# ax = subplot(223)
# plot(data[:,1], data[:,7], "k-")
# xlabel("t*")
# ylabel("Cd", color="k")
# xlim(0.0,7.0)
# ylim(-0.25,1)
# xticks(0:1:7)
# yticks(-0.5:0.25:1)

# ax2 = ax[:twinx]() # Create another axis on top of the current axis
# plot(data[:,1], data[:,2]*180/pi, "r--")
# font2 = ["color"=>"red"]
# ylabel("alpha",fontdict=font2)
# ylim(-10,40)
# yticks(-10:10:40, color="r")

# ax = subplot(224)
# plot(data[:,1], -data[:,8], "k-")
# xlabel("t*")
# ylabel("-Cm", color="k")
# xlim(0.0,7.0)
# ylim(-0.5,2)
# xticks(0:1:7)
# yticks(-0.5:0.5:2)

# ax2 = ax[:twinx]() # Create another axis on top of the current axis
# plot(data[:,1], data[:,2]*180/pi, "r--")
# font2 = ["color"=>"red"]
# ylabel("alpha",fontdict=font2)
# ylim(-10,40)
# yticks(-10:10:40, color="r")
