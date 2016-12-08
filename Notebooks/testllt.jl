include("../src/UNSflow.jl")
using UNSflow

cref = 1.
bref = 10.
sref = 10.

patch1 = patch(0., -5., 0., 0., "FlatPlate", 1., 0., 0.11, 10)
patch2 = patch(0., 0., 0., 0., "FlatPlate", 1., 0., 0.11, 10)
patch3 = patch(0., 5., 0., 0., "FlatPlate", 1., 0., 0.11, 5)
patchdata = [patch1; patch2; patch3]

alphadef = ConstDef(15.*pi/180)
hdef = ConstDef(0.)
udef = ConstDef(1.)
kin = KinemDef3D(alphadef, hdef, udef)

surf = ThreeDSurf(cref, bref, sref, patchdata, kin, 1., 70, 35, 21)
field = ThreeDFlowField()
dtstar = 0.015
nsteps =round(Int,5./dtstar) + 1
@time mat3d, surf3d, field3d = QSLLT_ldvm(surf, field, nsteps, dtstar)






