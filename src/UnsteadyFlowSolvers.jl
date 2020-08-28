#__precompile__(true)

module UnsteadyFlowSolvers

import Dierckx: Spline1D, derivative, evaluate

import ForwardDiff

import DelimitedFiles

import Serialization

import NLsolve: nlsolve, not_in_place

import Statistics: mean

import LinearAlgebra

import Dates

import Plots: @layout

import PyPlot: plot, scatter, figure, xlabel, ylabel, xlim, ylim,
    xticks, yticks, subplot, subplot2grid, legend, axis, savefig,
    close, tight_layout, clf, title, show, pause

import LaTeXStrings: @L_str

import DataStructures

export
    # kinematics types and funtions
    MotionDef,
    KinemPar,
    KinemDef,
    EldUpDef,
    EldUptstartDef,
    ConstDef,
    EldRampReturnDef,
    EldRampReturntstartDef,
    EldUpIntDef,
    EldUpInttstartDef,
    SinDef,
    CosDef,

    # 2D low-order solver types
    TwoDSurf,
    TwoDOFPar,
    KinemPar2DOF,
    TwoDVort,
    TwoDFlowField,
    KelvinCondition,
    KelvinKutta,

    # vortex count control utility
    delVortDef,
    delNone,
    delSpalart,

    # utility functions
    simpleTrapz,
    camber_calc,
    find_tstep,
    simpleInterp,
    cleanWrite,

    #2D low-order solver methods
    lautat,
    ldvm,
    ldvmLin,
    ldvm2DOF,
    LVE,

    # 2D plot output functions
    makeForcePlots2D,
    makeVortPlots2D,
    subPlot,

    # XFOIL Wrapper
    xfoilWrapper,

    # User Interface
    runUI

### source files

# kinematic types
include("kinem.jl")

# utility functions
include("utils.jl")

# vortex count control utility
include("delVort.jl")

# low-order 2D solvers
include("lowOrder2D/typedefs.jl")            # type definitions
include("lowOrder2D/calcs.jl")               # calculation functions
include("lowOrder2D/solvers.jl")             # solver methods
include("lowOrder2D/LVE.jl")                 # Lumped Vortex Element method
include("lowOrder2D/postprocess.jl")

# 2D plotting functions
include("plots/plots2D.jl")

# XFOIL Wrapper
include("xfoil/XfoilWrapper.jl")

# User Interface
include("UI/UI_main.jl")
include("UI/UI_utils.jl")
include("UI/s1Cmds.jl")
include("UI/s2Cmds.jl")

end
