#workspace()
#include("../src/UNSflow.jl")
#using UNSflow
#
#Uncomment the lines above when running for the first time

#This repeoduces results from Case 3A in Ramesh et al. (2014)
alpha_init = 10*pi/180
alphadot_init = 0.
h_init = 0.
hdot_init = 0.
u = 0.4667
udot = 0

kinem = KinemPar2DOF(alpha_init, h_init, alphadot_init, hdot_init, u, udot, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.)

x_alpha = 0.05
r_alpha = 0.5
kappa = 0.05
w_alpha = 1.
w_h = 1.
w_alphadot = 0.
w_hdot = 1.
cubic_h_1 = 1.
cubic_h_3 = 0.
cubic_alpha_1 = 1.
cubic_alpha_3 = 0.

strpar = TwoDOFPar(x_alpha, r_alpha, kappa, w_alpha, w_h, w_alphadot, w_hdot, cubic_h_1, cubic_h_3, cubic_alpha_1, cubic_alpha_3)
                   
dt = 0.015

pvt = 0.25

lespcrit = [21;]

pvt = 0.35 #leading edge

surf = TwoDSurf_2DOF(1., 1., "FlatPlate", pvt, 70, 35, strpar, kinem, lespcrit)

curfield = TwoDFlowField()

nsteps = 5000

ldvm(surf, curfield, nsteps, dt)

#Present plots
data = readdlm("results.dat")


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
