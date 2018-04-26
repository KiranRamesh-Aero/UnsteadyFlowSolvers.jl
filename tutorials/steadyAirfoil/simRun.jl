using UNSflow
#= define motion kinematics for the airfoil - Choose from several
#predefined kinematic types for pitch, heave and surge

# the documentation for MotionDef will show all available types =#

# sinusoidal pitch with 10 deg amplitude, reduced frequency = 1.0 
#alphadef = SinDef(0., 10*pi/180, 1.0, 0.0)
# constant pitch of 5 deg
alphadef = ConstDef(5.*pi/180)

# no plunge
hdef = ConstDef(0.)

# constant freestream
udef = ConstDef(1.)

# combine kinematics to form the full definition
full_kinem = KinemDef(alphadef, hdef, udef)

# define airfoil properties

# pitch axis location from 0-1
pvt = 0.25

# can define flat plate or provide airfoil file in working directory
# (in XFOIL format, starting directly with coordinates)
geometry = "FlatPlate"

# combined terms to form the full definition of 2D surface
surf = TwoDSurf(geometry, pvt, full_kinem)

curfield = TwoDFlowField()

dtstar = find_tstep(alphadef)

t_tot = 10.

nsteps =Int(round(t_tot/dtstar))+1

startflag = 0

writeflag = 1

writeInterval = t_tot/10.

#delvort = delSpalart(500, 12, 1e-5)
delvort = delNone()

mat, surf, curfield = lautatRoll(surf, curfield, nsteps, dtstar,startflag, writeflag, writeInterval, delvort)

