#workspace()
#include("../src/UNSflow.jl")
#using UNSflow

alphadef = EldUpDef(90*pi/180,0.2,11)
hdef = ConstDef(0.)
udef = ConstDef(1.)
full_kinem = KinemDef(alphadef, hdef, udef)

pvt = 1.0 #leading edge

lespcrit = [0.11;]

surf = TwoDSurf(1., 1., "sd7003_fine.dat", pvt, 70, 35, "Prescribed", full_kinem, lespcrit)

curfield = TwoDFlowField()

nsteps =round(Int,5/0.015)+1

ldvm(surf, curfield, nsteps)

#Present plots
data = readdlm("results.dat")
comp = readdlm("comp_eld_90.dat")

#Tranfer cm from LE to QC
cm = transfer_cm(0.25,data[:,8],data[:,6],data[:,7],data[:,2],1.,1.)

f = figure("pyplot_subplot_mixed", figsize=(15,8))

ax = subplot(221)
plot(data[:,1], data[:,5], "k-")
xlabel("t*")
ylabel("LESP")
xlim(0.0,5.0)
ylim(-0.3,0.3)
xticks(0:1:5)
yticks(-0.3:0.1:0.3)

ax2 = ax[:twinx]() # Create another axis on top of the current axis
plot(data[:,1], data[:,2]*180/pi, "r--")
font2 =Dict("color"=>"red")
ylabel("alpha",fontdict=font2)
ylim(0,90)
yticks(0:15:90, color="r")

ax = subplot(222)
plot(data[:,1], data[:,6], "k-")
plot(comp[:,1],comp[:,12]*2,"b--")
xlabel("t*")
ylabel("Cl", color="k")
xlim(0.0,5.0)
ylim(-1,5)
xticks(0:1:5)
yticks(-1:1:5)

ax2 = ax[:twinx]() # Create another axis on top of the current axis
plot(data[:,1], data[:,2]*180/pi, "r--")
font2 = Dict("color"=>"red")
ylabel("alpha",fontdict=font2)
ylim(0,90)
yticks(0:15:90, color="r")

ax = subplot(223)
plot(data[:,1], data[:,7], "k-")
plot(comp[:,1],comp[:,13]*2,"b--")
xlabel("t*")
ylabel("Cd", color="k")
xlim(0.0,8.0)
ylim(-1,5)
xticks(0:1:5)
yticks(-1:1:5)

ax2 = ax[:twinx]() # Create another axis on top of the current axis
plot(data[:,1], data[:,2]*180/pi, "r--")
font2 = Dict("color"=>"red")
ylabel("alpha",fontdict=font2)
ylim(0,90)
yticks(0:15:90, color="r")

ax = subplot(224)
plot(data[:,1], cm, "k-")
plot(comp[:,1],comp[:,14]*2,"b--")
xlabel("t*")
ylabel("Cm", color="k")
xlim(0.0,5.0)
ylim(-0.25,1.25)
xticks(0:1:5)
yticks(-0.25:0.25:1.25)

ax2 = ax[:twinx]() # Create another axis on top of the current axis
plot(data[:,1], data[:,2]*180/pi, "r--")
font2 = Dict("color"=>"red")
ylabel("alpha",fontdict=font2)
ylim(0,90)
yticks(0:15:90, color="r")
