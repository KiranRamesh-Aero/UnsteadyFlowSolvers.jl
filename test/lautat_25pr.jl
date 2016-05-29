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

f = figure("pyplot_subplot_mixed", figsize=(5,5))

ax = PyPlot.subplot(221)
plot(data[:,1], data[:,5], "k-")
PyPlot.xlabel("t*")
PyPlot.ylabel("LESP", color="k")
PyPlot.xlim(0.0,7.0)
PyPlot.ylim(-0.2,0.8)
PyPlot.xticks(0:1:7)
PyPlot.yticks(-0.2:0.2:0.8)

ax2 = PyPlot.twinx(ax)
PyPlot.axes(ax2)
plot(data[:,1], data[:,2]*180/pi, "r--")
PyPlot.ylabel("alpha (deg)", color='r')
PyPlot.ylim(-10,40)
PyPlot.yticks(-10:10:40, color="r")

ax = PyPlot.subplot(222)
plot(data[:,1], data[:,6], "k-")
PyPlot.xlabel("t*")
PyPlot.ylabel("Cl", color="k")
PyPlot.xlim(0.0,7.0)
PyPlot.ylim(-2,4)
PyPlot.xticks(0:1:7)
PyPlot.yticks(-2:1:4)



