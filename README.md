UNSflow is a set of low-order unsteady flow solvers managed by
Dr. Kiran Ramesh at the University of Glasgow. The solvers are written
in Julia, are based on the discrete-vortex method, cater to various
applications, and aim to simulate the aerodynamics (and aeroelastic
behavior) of airfoils/wings undergoing arbitrary unsteady motion.

[Julia](http://julialang.org) is a new high-level, high-performance dynamic programming
language for technical computing. Julia’s novel features are a
sophisticated compiler, distributed parallel execution, numerical
accuracy, and an extensive mathematical function library. Julia is
fully open-source under the MIT license, and integrates mature,
best-of-breed open source C and Fortran libraries for various
computing algorithms. 

### Installing Julia and UNSflow
[Download](http://julialang.org/downloads/) Julia vesion 0.5

Run the application, and at the prompt type, `Pkg.clone("https://github.com/KiranUofG/UNSflow")`

This will install UNSflow and the packages it depends on.

### Getting started
The simplest way to explore the codes and their capabilites is to work
out the examples provided. The examples are provided through [Jupyter](http://jupyter.org/) browser-based "notebooks". 
You can launch notebooks within Julia (using the IJulia package) by running: 
```julia
using IJulia 
notebook()
```

### Validated codes in UNSflow (published in journals):

##### Large-Angle Unsteady Thin-Airfoil Theory (LAUTAT)

A time-stepping method based on principles from traditional
thin-airfoil theory. No small-anlge approximations are made and
nonplanar wakes are modelled. An unsteady potential flow solution for
an airfoil undergoing any arbitrary motion can be generated with this code.

##### LESP-Modulated Discrete-Vortex Method (LDVM)

The Leading Edge Suction Parameter is a theoretical construct
associated with Leading Edge Vortex (LEV) formation in unsteady
flows. The LAUTAT code is augmented with this concept which is used to
predict and modulate LEV shedding. The LDVM code can be used to predict
the aerodynamics resulting from high-frequency maneuvers
(e.g. flapping wing flight) of an airfoil where LEVs play a
significant role.

##### LDVM-2DOF

The LDVM is coupled with a structural model to simulate the
aeroelastic behaviour of an airfoil that is elastically supported by
translational and torsional springs in heave and pitch
respectively. Problems that can be studied include flutter, limit cycles
resulting from LEV shedding.


### Experimental codes (development in progress):








