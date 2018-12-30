alphadef = UNSflow.ConstDef(5. *pi/180)
hdef = UNSflow.ConstDef(0.)
udef = UNSflow.ConstDef(1.)
full_kinem = UNSflow.KinemDef(alphadef, hdef, udef)
pvt = 0.25
geometry = "FlatPlate"
lespc = [10.15;]

#The value of leading-edge radius (rho) must be specified according to geometry
surf = UNSflow.TwoDSurf(geometry, pvt, full_kinem, lespc, rho=0.016)

@test typeof(surf)==TwoDSurf

curfield = UNSflow.TwoDFlowField()

@test typeof(curfield) == TwoDFlowField

dtstar = UNSflow.find_tstep(alphadef)

@test dtstar== 0.0150

t_tot = 3.

nsteps =Int(round(t_tot/dtstar))+1

@test nsteps==201

startflag = 0

writeflag = 0

writeInterval = t_tot/10.

#delvort = delSpalart(500, 12, 1e-5)
delvort = UNSflow.delNone()

mat, surf, curfield = UNSflow.ldvmLin(surf, curfield, nsteps, dtstar,startflag, writeflag, writeInterval, delvort)

q_u,q_l=UNSflow.calc_edgeVel(surf,[curfield.u[1],curfield.w[1]])

#@test sum(q_u)>0
#@test sum(q_l)>0 

@test sum(q_u)==0
@test sum(q_l)==0


@test sum(q_u)!=NaN
@test sum(q_l)!=NaN
