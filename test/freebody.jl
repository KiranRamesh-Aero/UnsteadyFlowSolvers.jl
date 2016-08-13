workspace()
include("../src/UNSflow.jl")
using UNSflow
#
#Uncomment the lines above when running for the first time

#This is a case study of bringing down a drone using a concentrated
#vortex. The drone may be hovering or flying forward/back. The angle
#and velocity of the shot may be provided. The trajectory of the drone
#is observed

#Air
rho = 1.23
accl_g = 9.8

#Drone properties
c = 0.2 #m
uref = 5. #m/s
lespcrit = [0.1;] #Relatively sharp leading edge
m = 10. #kg
kappa = pi*rho*c*c/(4*m)

#Initial conditions 
alpha_init = 0*pi/180
alphadot_init = 0.
h_init = 10.
hdot_init = 0.
u = 0.
udot = 0.

I_g = (m/5)*(c^2+(0.2*c)^2) #Assuming its an ellipse with a = c and b = 0.2*c 
x_g = 0.5*c
pvt = x_g
r_g = 2*sqrt(I_g/(m*c*c))
strpar = TwoDFreePar(r_g, x_g, kappa)

#Shot properties
sh_vel = 100
sh_angle = 33*pi/180.
sh_vcore = 2*c
sh_str = 10.
sh_xpos = -15.
sh_zpos = 0.

#xvel = 100*cos(angle*pi/180)
#zvel = 100*sin(angle*pi/180)

#Sim properties
dt = 0.015
kinem = KinemPar2DFree(alpha_init, h_init, alphadot_init, hdot_init, u, udot, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.)
surf = TwoDFreeSurf(c, uref, "FlatPlate", pvt, 70, 35, strpar, kinem, lespcrit)
curfield = TwoDFlowField()
push!(curfield.extv, TwoDVort(sh_xpos, sh_zpos, sh_str, sh_vcore, sh_vel*cos(sh_angle), sh_vel*sin(sh_angle)))   
nsteps = 50

cf = pi*surf.c*accl_g/(2*surf.strpar.kappa*surf.uref*surf.uref) #The drone has a lifting force equalling its weight

mat, surf, frames = ldvm(surf, curfield, nsteps, dt, cf)

figure(1)
scat1 = scatter(map(q->q.x, frames[1].tev),map(q->q.z,frames[1].tev),s=20,c=map(q->q.s,frames[1].tev),cmap=PyPlot.ColorMap("jet"),edgecolors="none")
scat2 = scatter(map(q->q.x, frames[1].lev),map(q->q.z,frames[1].lev),s=20,c=map(q->q.s,frames[1].lev),cmap=PyPlot.ColorMap("jet"),edgecolors="none")
scat3 = scatter(map(q->q.x, frames[1].extv),map(q->q.z,frames[1].extv),s=40,c=map(q->q.s,frames[1].extv),cmap=PyPlot.ColorMap("jet"),edgecolors="none")
p1 = plot(map(q->q.x, frames[1].bv),map(q->q.z,frames[1].bv),color = "black",linewidth=2.0)

