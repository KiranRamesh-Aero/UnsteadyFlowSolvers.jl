alphadef = ConstDef(5. *pi/180)
hdef = ConstDef(0.)
udef = ConstDef(1.)
full_kinem = KinemDef(alphadef, hdef, udef)
pvt = 0.25
geometry = "FlatPlate"
lespc = [10.15;]

surf = TwoDSurf(geometry, pvt, full_kinem, lespc, rho=0.016)

@test typeof(surf)==TwoDSurf

curfield = TwoDFlowField()

@test typeof(curfield) == TwoDFlowField

dtstar = find_tstep(alphadef)

t_tot = 7.

nsteps = Int(round(t_tot/dtstar))+1

startflag = 0

writeflag = 0

writeInterval = t_tot/10.

delvort = delNone()

mat, surf, curfield = ldvmLin(surf, curfield, nsteps, dtstar,startflag, writeflag, writeInterval, delvort)

cl_expected = 2*pi*alphadef.amp
cl_out = mat[end,6]

@test abs(cl_expected - cl_out) < 0.5*cl_expected
