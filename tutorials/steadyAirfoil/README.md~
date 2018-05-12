This tutorial simulates an airfoil in constant freestream
undergoing unsteady motion using the LAUTAT (Large-Angle Unsteady
Thin-Airfoil Theory) solver.
		
The theory is described in Ramesh, K. et al., "An unsteady airfoil
theory applied to pitching motions validated against experiment and
computation", Theor. Comput. Fluid Dyn. (2013) 27:
843. [Weblink](https://doi.org/10.1007/s00162-012-0292-8)

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
