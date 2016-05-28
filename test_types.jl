workspace()
include("UNSflow.jl")
using UNSflow

alphadef = EldUpDef(45,0.2,0.8)
hdef = ConstDef(0.)
udef = ConstDef(1.)
full_kinem = KinemDef(alphadef, hdef, udef)
surf = TwoDSurf(1., 1., "sd7003_fine.dat", 0.35, 70, 35, "Prescribed", full_kinem)


dtstar = 0.015
dt = dtstar*surf.c/surf.uref
nsteps = 100
t = 0.

curfield = TwoDFlowField()

#Intialise flowfield

for istep = 1:100
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

    #Calculate velocities induced by all free vortices by each other
    for i = 1:length(curfield.tev)
        curfield.tev[i].vx = 0
        curfield.tev[i].vz = 0
    end
    for i = 1:length(curfield.lev)
        curfield.lev[i].vx = 0
        curfield.lev[i].vz = 0
    end

    mutual_ind([curfield.tev; curfield.lev])

    #Add the influence of velocities induced by bound vortices
    utemp = zeros(length(curfield.tev)+length(curfield.lev))
    wtemp = zeros(length(curfield.tev)+length(curfield.lev))
    utemp, wtemp = ind_vel(surf.bv, [map(q -> q.x, curfield.tev); map(q -> q.x, curfield.lev)], [map(q -> q.z, curfield.tev); map(q -> q.z, curfield.lev)])

    for i = 1:length(curfield.tev)
        curfield.tev[i].vx += utemp[i]
        curfield.tev[i].vz += wtemp[i]
    end
    for i = length(curfield.tev)+1:length(utemp)
        curfield.lev[i-length(curfield.tev)].vx += utemp[i]
        curfield.lev[i-length(curfield.tev)].vz += wtemp[i]
    end

    #Convect free vortices with their induced velocities
    for i = 1:length(curfield.tev)
        curfield.tev[i].x += dt*curfield.tev[i].vx
        curfield.tev[i].z += dt*curfield.tev[i].vz
    end
    for i = 1:length(curfield.lev)
        curfield.lev[i].x += dt*curfield.lev[i].vx
        curfield.lev[i].z += dt*curfield.lev[i].vz
    end

end

