
alphadef = UNSflow.ConstDef(4.0 * pi/180)

hdef = UNSflow.SinDef(0., 0.05, 3.93, 0.)

udef = UNSflow.ConstDef(1.)

full_kinem = UNSflow.KinemDef(alphadef, hdef, udef)

pvt = 0.25

geometry = "sd7003.dat"

surf = UNSflow.TwoDSurf(geometry, pvt, full_kinem)

@test typeof(surf)==TwoDSurf

curfield = UNSflow.TwoDFlowField()

@test typeof(curfield) == TwoDFlowField


dtstar = UNSflow.find_tstep(hdef)

t_tot = 5. *pi/hdef.k

nsteps =Int(round(t_tot/dtstar))+1

startflag = 0

writeflag = 1

writeInterval = t_tot/20.

#delvort = delSpalart(500, 12, 1e-5)
delvort = UNSflow.delNone()

mat, surf, curfield = UNSflow.lautat(surf, curfield, nsteps, dtstar,startflag, writeflag, writeInterval, delvort, wakerollup=1)

q_u,q_l=UNSflow.calc_edgeVel(surf,[curfield.u[1],curfield.w[1]])

@test sum(q_u)>0
@test sum(q_l)>0


@test sum(q_u)!=NaN
@test sum(q_l)!=NaN
