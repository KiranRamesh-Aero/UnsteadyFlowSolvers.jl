workspace()
include("UNSflow.jl")
using UNSflow

outfile = open("results.dat", "w")

alphadef = EldUpDef(45,0.2,0.8)
hdef = ConstDef(0.)
udef = ConstDef(1.)
full_kinem = KinemDef(alphadef, hdef, udef)
surf = TwoDSurf(1., 1., "sd7003_fine.dat", 0.35, 70, 35, "Prescribed", full_kinem)
#surf = TwoDSurf(1., 1., "FlatPlate", 0.35, 70, 35, "Prescribed", full_kinem)

dtstar = 0.015
dt = dtstar*surf.c/surf.uref
nsteps = 500
t = 0.

curfield = TwoDFlowField()

#Intialise flowfield
for istep = 1:nsteps
        #Udpate current time
    t = t + dt

    #Update kinematic parameters
    update_kinem(surf, t)

    #Update bound vortex positions
    update_boundpos(surf, dt)

    #Add a TEV with dummy strength
    place_tev(surf,curfield,dt)

    kelv = KelvinCondition(surf,curfield)
    #Solve for TEV strength to satisfy Kelvin condition
    curfield.tev[length(curfield.tev)].s = secant_method(kelv, 0., -0.01)

    #Update adot
    update_a2a3adot(surf,dt)

    #Check for LEV and shed if yes

    #Update rest of Fourier terms
    update_a2toan(surf)

    #Calculate bound vortex strengths
    update_bv(surf)

    wakeroll(surf, curfield, dt)

    write(outfile, join((t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, surf.a0[1])," "), "\n")

end

close(outfile)

