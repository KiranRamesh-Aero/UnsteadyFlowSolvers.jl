workspace()
include("../src/UNSflow.jl")
using UNSflow


alphadef = EldRampReturnDef(25,0.11,11)
hdef = ConstDef(0.)
udef = ConstDef(1.)
full_kinem = KinemDef(alphadef, hdef, udef)

pvt = 0.0 #leading edge

surf = TwoDSurf(1., 1., "sd7003_fine.dat", pvt, 70, 35, "Prescribed", full_kinem)

curfield = TwoDFlowField()

nsteps =round(Int,7/0.015)+1

lautat(surf, curfield, nsteps)

#Present plots
data = readdlm("results.dat")

f = figure("pyplot_subplot_mixed", figsize=(12,8))

ax = subplot(221)
plot(data[:,1], data[:,5], "k-")
xlabel("t*")
font1 = ["color"=>"black"]
ylabel("LESP", font1)
xlim(0.0,7.0)
ylim(-0.2,0.8)
xticks(0:1:7)
yticks(-0.2:0.2:0.8)

ax2 = ax[:twinx]() # Create another axis on top of the current axis
plot(data[:,1], data[:,2]*180/pi, "r--")
font2 = ["color"=>"red"]
ylabel("alpha",fontdict=font2)
ylim(-10,40)
yticks(-10:10:40, color="r")

ax = subplot(222)
plot(data[:,1], data[:,6], "k-")
xlabel("t*")
ylabel("Cl", color="k")
xlim(0.0,7.0)
ylim(-1,4)
xticks(0:1:7)
yticks(-2:1:4)

ax2 = ax[:twinx]() # Create another axis on top of the current axis
plot(data[:,1], data[:,2]*180/pi, "r--")
font2 = ["color"=>"red"]
ylabel("alpha",fontdict=font2)
ylim(-10,40)
yticks(-10:10:40, color="r")

ax = subplot(223)
plot(data[:,1], data[:,7], "k-")
xlabel("t*")
ylabel("Cd", color="k")
xlim(0.0,7.0)
ylim(-0.25,1)
xticks(0:1:7)
yticks(-0.5:0.25:1)

ax2 = ax[:twinx]() # Create another axis on top of the current axis
plot(data[:,1], data[:,2]*180/pi, "r--")
font2 = ["color"=>"red"]
ylabel("alpha",fontdict=font2)
ylim(-10,40)
yticks(-10:10:40, color="r")

ax = subplot(224)
plot(data[:,1], -data[:,8], "k-")
xlabel("t*")
ylabel("-Cm", color="k")
xlim(0.0,7.0)
ylim(-0.5,2)
xticks(0:1:7)
yticks(-0.5:0.5:2)

ax2 = ax[:twinx]() # Create another axis on top of the current axis
plot(data[:,1], data[:,2]*180/pi, "r--")
font2 = ["color"=>"red"]
ylabel("alpha",fontdict=font2)
ylim(-10,40)
yticks(-10:10:40, color="r")
