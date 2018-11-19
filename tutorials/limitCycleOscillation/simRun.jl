push!(LOAD_PATH,"../../src/")
import UNSflow

alpha_init = 10. *pi/180
alphadot_init = 0.
h_init = 0.
hdot_init = 0.
u = 0.467
udot = 0
kinem = UNSflow.KinemPar2DOF(alpha_init, h_init, alphadot_init, hdot_init, u)

x_alpha = 0.05
r_alpha = 0.5
kappa = 0.05
w_alpha = 1.
w_h = 1.
w_alphadot = 0.
w_hdot = 0.
cubic_h_1 = 1.
cubic_h_3 = 0.
cubic_alpha_1 = 1.
cubic_alpha_3 = 0.
strpar = UNSflow.TwoDOFPar(x_alpha, r_alpha, kappa, w_alpha, w_h, w_alphadot, w_hdot, cubic_h_1, cubic_h_3, cubic_alpha_1, cubic_alpha_3)

#Dummy kinematic definitions to initialise surface
alphadef = UNSflow.ConstDef(alpha_init)
hdef = UNSflow.ConstDef(h_init)
udef = UNSflow.ConstDef(1.) #This is relative to uref
startkinem = UNSflow.KinemDef(alphadef, hdef, udef)

lespcrit = [0.11;]

pvt = 0.35

c = 1.

surf = UNSflow.TwoDSurf("FlatPlate", pvt, startkinem, lespcrit, c=c, uref=u)

curfield = UNSflow.TwoDFlowField()

dtstar = 0.015

nsteps = 50000

t_tot = nsteps * dtstar / u

startflag = 0

writeflag = 0
writeInterval = dtstar * nsteps/20.

writeInterval = t_tot/20.

delvort = UNSflow.delSpalart(500, 12, 1e-5)

mat, surf, curfield = UNSflow.ldvm2DOF(surf, curfield, strpar, kinem, nsteps, dtstar,startflag, writeflag, writeInterval, delvort)

UNSflow.makeForcePlots2D()

#cleanWrite()
