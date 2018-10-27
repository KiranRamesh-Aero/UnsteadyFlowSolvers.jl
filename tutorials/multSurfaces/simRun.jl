push!(LOAD_PATH,"../../src/")
import UNSflow


alphadef1 = UNSflow.ConstDef(45. *pi/180)
alphadef2 = UNSflow.ConstDef(5. *pi/180)

hdef = UNSflow.ConstDef(0.)

udef = UNSflow.ConstDef(1.)

full_kinem1 = UNSflow.KinemDef(alphadef1, hdef, udef)
full_kinem2 = UNSflow.KinemDef(alphadef2, hdef, udef)

pvt = 0.0

geometry = "FlatPlate"

lespcrit = [0.11;]

surf1 = UNSflow.TwoDSurf(geometry, pvt, full_kinem1, lespcrit)
surf2 = UNSflow.TwoDSurf(geometry, pvt, full_kinem2, lespcrit, initpos = [15.5 ;-0.5 ])

surf = [surf1;surf2]

curfield = UNSflow.TwoDFlowField()

dtstar = UNSflow.find_tstep(alphadef)

t_tot = 5.

nsteps = Int(round(t_tot/dtstar))+1

startflag = 0

writeflag = 1

writeInterval = t_tot/18.

#delvort = delSpalart(500, 12, 1e-5)
delvort = UNSflow.delNone()

mat, surf, curfield = UNSflow.ldvm(surf, curfield, nsteps, dtstar,startflag, writeflag, writeInterval, delvort)

UNSflow.makeVortPlots2D()

UNSflow.makeForcePlots2D()

UNSflow.cleanWrite()
