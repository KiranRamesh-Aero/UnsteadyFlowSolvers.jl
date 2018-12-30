# UNSflow

[![Build Status](https://travis-ci.com/desanga/UNSflow.jl.svg?branch=master)](https://travis-ci.com/desanga/UNSflow.jl)
[![Coverage Status](https://coveralls.io/repos/github/desanga/UNSflow.jl/badge.svg?branch=master)](https://coveralls.io/github/desanga/UNSflow.jl?branch=master)
[![codecov](https://codecov.io/gh/desanga/UNSflow.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/desanga/UNSflow.jl)

UNSflow is a library of low-order solvers for unsteady aerodynamics
and aeroelasticity managed by Dr. Kiran Ramesh at the University of Glasgow. The solvers.
are written in Julia, are based on the discrete-vortex method and cater to various
applications.

This project is currently supported by EPSRC grant EP/R008035/1.

[Julia](http://julialang.org) is a new high-level, high-performance dynamic programming
language for technical computing. Julia’s novel features are a
sophisticated compiler, distributed parallel execution, numerical
accuracy, and an extensive mathematical function library. Julia is
fully open-source under the MIT license, and integrates mature,
best-of-breed open source C and Fortran libraries for various
computing algorithms.

### Installing Julia and UNSflow
[Download](http://julialang.org/downloads/) the latest stable version of Julia.

Run the application, and at the prompt type, `Pkg.clone("https://github.com/KiranUofG/UNSflow")`

This will install UNSflow and the packages it depends on.

### Getting started
The tutorials provide some sample simulations that can be performed with the code. To run
your own case, copy a tutorial directory and modify the input file "simRun.jl" as necessary.  
