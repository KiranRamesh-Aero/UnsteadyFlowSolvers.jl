push!(LOAD_PATH,"../../src/")
using UNSflow

alpha_init = 10*pi/180
alphadot_init = 0.
h_init = 0.
hdot_init = 0.
u = 0.467
udot = 0
kinem = KinemPar2DOF(alpha_init, h_init, alphadot_init, hdot_init, u)

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
strpar = TwoDOFPar(x_alpha, r_alpha, kappa, w_alpha, w_h, w_alphadot, w_hdot, cubic_h_1, cubic_h_3, cubic_alpha_1, cubic_alpha_3)

lespcrit = [0.11;]

pvt = 0.35

c = 1.

surf = TwoDSurf2DOF(c, u, "FlatPlate", pvt, strpar, kinem, lespcrit)

curfield = TwoDFlowField()

dtstar = 0.015

nsteps = 50000

startflag = 0

writeflag = 0
writeInterval = dtstar * nsteps/20.

delvort = delSpalart(500, 12, 1e-5)

mat, surf, curfield = ldvm(surf, curfield, nsteps, dtstar,startflag, writeflag, writeInterval, delvort)

makeForcePlots()

#cleanWrite()
