This tutorial simulates using potential flow theory an airfoil in
constant freestream at a constant angle of attack.

The simulation can be run by entering in the julia prompt from the
current directory,

```
include("simRun.jl")
```

The simulation parameters are input in simRun.jl as follows. First,
the UNSflow library is loaded

```
push!(LOAD_PATH,"../../src/")
using UNSflow
```

Motion kinematics for the airfoil are to be defined in term of pitch
(rotation), plunge/heave (vertical translation) and surge (horizontal
translation). The documentation for type MotionDef will show all available
kinematic definitons.

```
?MotionDef
```

Here, we define a constant pitch angle of 5 deg, no plunge, and a
constant freestream velocity. 

```
alphadef = ConstDef(5.*pi/180)
hdef = ConstDef(0.)
udef = ConstDef(1.)
full_kinem = KinemDef(alphadef, hdef, udef)
```

The surface is defined in terms of it's kinematics and geometry and a
2D flowfield is initialised for the simulation. The airfoil properties
- geometry and pitch axis location are to be provided. A flat plat can
be defined, or a custom geometry from a coordinate file in the working
directory. The coordinates must be in XFOIL format (trailing edge ->
lower surface -> leading edge -> upper surface -> trailing edge) with
no headers (starting directly from coordinates).

```
pvt = 0.25
geometry = "FlatPlate"
surf = TwoDSurf(geometry, pvt, full_kinem)
curfield = TwoDFlowField()
```

The simulation time step for typical kinematics can be calculated
using the find_tstep(kin::MotionDef) function. The total run time is
defined using the number of time steps. 

```
dtstar = find_tstep(alphadef)
t_tot = 10.
nsteps =Int(round(t_tot/dtstar))+1
```

writeflag=1 is used to write the simulation details (vortex strengths
and locations) at a frequency defined using writeInterval. startflag=1
can be used to restart and continue a simulation which is not used here.

```
startflag = 0
writeflag = 1
writeInterval = t_tot/10.
delvort = delNone()

```

```
?lautat
```

to see all solver options.

Documentation for all custom types seen in simRun.jl can be accessed
the same way.

```
include("simRun.jl")
```

to define the simulation paramters and run.

```
include("postRun.jl")
```

to create vorticity maps from output timestamp directories. writeflag
must be set to 1 in simRun.jl for this. PyPlot libraries will need to
be configured correctly on your system. 

```
cleanWrite()
```

to clear all timestamp directories
