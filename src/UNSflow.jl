#__precompile__(true)

module UNSflow

import Dierckx: Spline1D, derivative, evaluate

import ForwardDiff

import DelimitedFiles

import Serialization

import NLsolve: nlsolve, not_in_place

import Statistics: mean

import PyPlot: plot, scatter, figure, xlabel, ylabel, xlim, ylim,
    xticks, yticks, subplot, subplot2grid, legend, axis, savefig,
    close, tight_layout


import Plots: @layout


import LaTeXStrings: @L_str

#For use in development and debugging
import Revise

export
    # kinematics types and funtions
    MotionDef,
    KinemPar,
    KinemDef,
    EldUpDef,
    EldUptstartDef,
    ConstDef,
    EldRampReturnDef,
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

    # 3D low-order solver types
    ThreeDSurfSimple,
    KinemDef3D,
    ThreeDFieldSimple,

    # utility functions
    simpleTrapz,
    camber_calc,
    find_tstep,
    simpleInterp,
    cleanWrite,

    # 2D low-order solver function
    update_boundpos,
    update_kinem,
    update_indbound,
    update_downwash,
    update_a0anda1,
    place_tev,
    place_lev,
    update_a2toan,
    mutual_ind,
    update_a2a3adot,
    update_bv,
    ind_vel,
    wakeroll,
    update_adot,
    update_externalvel,
    controlVortCount,
    update_kinem,
    update_kinem2DOF,

    #2D low-order solver methods
    lautat,
    #lautatRoll,
    ldvm,
    ldvmLin,

    #3D low-order solver methods
    QSLLT_lautat,
    QSLLT_ldvm,

    # Postprocessing functions
    calc_forces,
    writeStamp,

    # 2D plotting functions
    viewVort2D,
    viewVortConnect2D,

    # 2D plot output functions
    makeForcePlots2D,
    makeVortPlots2D,

    # 3D plot output functions
    makeForcePlots3D,
    makeVortPlots3D,
    makeInfoPlots3D,

    invisicidTransport,
    fluxSplittingParameters,
    solutions,
    operationalConditions

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
include("lowOrder2D/postprocess.jl")         # postprocessing functions

# low-order 3D solvers
include("lowOrder3D/typedefs.jl")            # type definitions
include("lowOrder3D/calcs.jl")               # calculation functions
include("lowOrder3D/solvers.jl")             # solver methods
include("lowOrder3D/postprocess.jl")         # postprocessing functions

# 2D plotting functions
include("plots/plots2D.jl")
include("plots/plots3D.jl")


# IBL solver functions
include("IBLTypedefs.jl")
include("IBLUNSflow.jl")
include("fluxSplitting/StegerWarming.jl")
#include("testscript.jl")
#include("testingfile.jl")
end
