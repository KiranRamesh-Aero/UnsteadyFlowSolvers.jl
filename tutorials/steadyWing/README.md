This tutorial simulates using potential flow theory a finite wing in
constant freestream at a constant angle of attack.

Refer to the *steadyAirfoil* for general guideline in setting up the kinematics, 
surface and solver definitons. 

The simulation can be run by entering in the julia prompt from the
current directory,

```
include("simRun.jl")
```

Here, we define a wing of aspect ratio 6, constant pitch angle of 5
deg and a constant freestream velocity. 

```
push!(LOAD_PATH,"../../src/")
using UNSflow

AR = 6. 
alphadef = UNSflow.ConstDef(10. *pi/180)
hdef = UNSflow.ConstDef(0.)
udef = UNSflow.ConstDef(1.)
full_kinem = UNSflow.KinemDef(alphadef, hdef, udef)
```

Sectional properties are defined and 3D surface and flowfield are initialised. 

```
pvt = 0.25
geometry = "FlatPlate"
surf = UNSflow.ThreeDSurfSimple(AR, geometry, pvt, full_kinem)
curfield = UNSflow.ThreeDFieldSimple()
```

The simulation parameters are defined. 

```
dtstar = 0.015
t_tot = 8. 
nsteps =Int(round(t_tot/dtstar))+1
```

`writeInterval` is set so that 8 time instants are written 
for creating the vorticity plots. 

```
startflag = 0
writeflag = 1
writeInterval = t_tot/8.
delvort = delNone()
```

Finally, the input parameters are provided to the 3D unsteady solver
function `QSLLT_lautat`. This solver combines the lautat 2D solver
(see "steadyAirfoil" and "plungingAirfoil" tutorials) with a lifting
line theory to model spanwise variation of circulation. The 2D theory
is described in Ramesh, K. et al., "An unsteady airfoil theory applied
to pitching motions validated against experiment and computation",
Theor. Comput. Fluid Dyn. (2013) 27:
843. [Weblink](https://doi.org/10.1007/s00162-012-0292-8)

`QSLLT_lautat` returns a matrix containing time variation of wing
force coefficients + various sectional quantities of interest, and the
final surface and flowfield data structures at the end of the
simulation. The time variation matrix is also written as a file to the
current directory, *resultsSummary*.

```
mat, surf, curfield = UNSflow.QSLLT_lautat(surf, curfield, nsteps, dtstar,startflag, writeflag, writeInterval, delvort)
```

`makeForcePlots3D()` can be used to create plots of the time variation
of force coefficients during the simulation; they are written to a
directory *forcePlots*. 

`makeVortPlots3D()` can be used to create vorticity maps from output
timestamp directories (written to *vortPlots*). Sectional vorticity
maps at the wing root, mid and tip are plotted. This can be used only
if `writeflag=1` was used and data directories exist in the current
directory. Since the images for various times are created with the
same axes, they can be easily compiled into a movie of the simulation.

`makeInfoPlots3D()` writes a slide containing the vorticity maps,
kinematics, forces and sectional variations at the times specified
through `writeInterval`. These plots again can be easily compiled into a
movie of the simulation.

`cleanWrite()` clear all
timestamp directories from the current directory.

```
makeForcePlots3D()
makeVortPlots3D()
makeInfoPlots3D()
cleanWrite()
```

