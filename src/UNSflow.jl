__precompile__(true)

module UNSflow

using Dierckx: Spline1D, derivative, evaluate

using ForwardDiff

using NLsolve: nlsolve, not_in_place

using PyCall
pygui(:tk)

using PyPlot: plot, scatter, figure, xlabel, ylabel, xlim, ylim,
xticks, yticks, subplot, legend, axis, savefig, close

using LaTeXStrings: @L_str

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
    TwoDSurfThick,
    TwoDOFPar,
    KinemPar2DOF,
    TwoDSurf2DOF,
    TwoDVort,
    TwoDFlowField,
    KelvinCondition,
    KelvinCondition2DOF,
    KelvinKutta,
    KelvinKutta2DOF,

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

    #2D low-order solver methods
    lautat,
    lautatRoll,
    ldvm,
    ldvmLin,

    # Postprocessing functions
    calc_forces,
    writeStamp,
    calc_q_cp,

    # 2D plotting functions
    viewVort2D,
    viewVortConnect2D,

    # 2D plot output functions
    makeForcePlots,
    makeVortPlots2D


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

# 2D plotting functions
include("plots/plots2D.jl")

end
