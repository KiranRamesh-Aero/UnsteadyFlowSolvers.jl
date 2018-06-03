push!(LOAD_PATH,"../../src/")
using UNSflow

AR = 4.

#alphadef = EldUpDef(90.*pi/180, 0.4, 0.8)
alphadef = ConstDef(45.*pi/180)

hdef = ConstDef(0.)

udef = ConstDef(1.)

kinem = KinemDef3D(alphadef, hdef, udef)

pvt = 0.0

geometry = "FlatPlate"

lespcrit = [0.11;]

surf = ThreeDSurfSimple(AR, kinem, geometry, pvt, lespcrit)

field = ThreeDFieldSimple()

dtstar = find_tstep(alphadef)

t_tot = 10.0

nsteps = Int(round(t_tot/dtstar))+1

startflag = 0

writeflag = 1

writeInterval = t_tot/10.

#delvort = delSpalart(500, 12, 1e-5)
delvort = delNone()

#maxwrite = 100
#nround=6

mat, surf, curfield = QSLLTldvmLin(surf, field, nsteps, dtstar,startflag, writeflag, writeInterval, delvort)


makeForcePlots3Dstrip()

makeVortPlots3Dstrip()

makeTevstrPlots3Dstrip()

cleanWrite()
