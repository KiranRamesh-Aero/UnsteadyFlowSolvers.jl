#workspace()
#include("../src/UNSflow.jl")
#using UNSflow

alphadef = EldRampReturnDef(25*pi/180,0.11,11)
hdef = ConstDef(0.)
udef = ConstDef(1.)
full_kinem = KinemDef(alphadef, hdef, udef)

pvt = 0.0 #leading edge

lespcrit = [0.18;]

surf = TwoDSurf(1., 1., "sd7003_fine.dat", pvt, 70, 35, "Prescribed", full_kinem, lespcrit)

curfield = TwoDFlowField()

nsteps =round(Int,7/0.015)+1

ldvm(surf, curfield, nsteps)

#Present plots
data = readdlm("results.dat")
exp = readdlm("exp_25pr.dat")
comp = readdlm("comp_25pr.dat")

#Tranfer cm from LE to QC
cm = transfer_cm(0.25,data[:,8],data[:,6],data[:,7],data[:,2],0.,1.)

f = figure("pyplot_subplot_mixed", figsize=(15,8))

ax = subplot(221)
plot(data[:,1], data[:,5], "k-")
xlabel("t*")
ylabel("LESP")
xlim(0.0,7.0)
ylim(-0.1,0.3)
xticks(0:1:7)
yticks(-0.1:0.1:0.3)

ax2 = ax[:twinx]() # Create another axis on top of the current axis
plot(data[:,1], data[:,2]*180/pi, "r--")
font2 =Dict("color"=>"red")
ylabel("alpha",fontdict=font2)
ylim(-10,30)
yticks(-10:10:30, color="r")

ax = subplot(222)
plot(data[:,1], data[:,6], "k-")
plot(comp[:,1]+0.1,comp[:,5],"b--")
plot(exp[:,1],exp[:,3]*2,"g-.")
xlabel("t*")
ylabel("Cl", color="k")
xlim(0.0,7.0)
ylim(-1,3)
xticks(0:1:7)
yticks(-1:1:3)

ax2 = ax[:twinx]() # Create another axis on top of the current axis
plot(data[:,1], data[:,2]*180/pi, "r--")
font2 = Dict("color"=>"red")
ylabel("alpha",fontdict=font2)
ylim(-10,30)
yticks(-10:10:30, color="r")

ax = subplot(223)
plot(data[:,1], data[:,7], "k-")
plot(comp[:,1]+0.1,comp[:,6],"b--")
plot(exp[:,1],exp[:,4]*2,"g-.")
xlabel("t*")
ylabel("Cd", color="k")
xlim(0.0,7.0)
ylim(-0.5,1.5)
xticks(0:1:7)
yticks(-0.5:0.5:1.5)

ax2 = ax[:twinx]() # Create another axis on top of the current axis
plot(data[:,1], data[:,2]*180/pi, "r--")
font2 = Dict("color"=>"red")
ylabel("alpha",fontdict=font2)
ylim(-10,30)
yticks(-10:10:30, color="r")

ax = subplot(224)
plot(data[:,1], data[:,8], "k-")
plot(comp[:,1]+0.1,-comp[:,7],"b--")
plot(exp[:,1],-exp[:,5]*2,"g-.")
xlabel("t*")
ylabel("Cm", color="k")
xlim(0.0,7.0)
ylim(-0.8,0.8)
xticks(0:1:7)
yticks(-0.8:0.4:0.8)

ax2 = ax[:twinx]() # Create another axis on top of the current axis
plot(data[:,1], data[:,2]*180/pi, "r--")
font2 = Dict("color"=>"red")
ylabel("alpha",fontdict=font2)
ylim(-10,30)
yticks(-10:10:30, color="r")
