module UNSflow

using Dierckx
export Spline1D, derivative, evaluate

using ForwardDiff
export derivative

using PyPlot
export plot, scatter, figure, xlabel, ylabel, xlim, ylim, xticks, yticks, subplot, axes, legend, markers, axis

using NLsolve
export nlsolve, not_in_place



export camber_calc, update_boundpos, update_kinem, update_indbound,
update_downwash, update_a0anda1, place_tev, place_lev, update_a2toan,
mutual_ind, trapz, update_a2a3adot, update_bv, ind_vel, view_vorts,
wakeroll, lautat, lautat_wakeroll, ldvm, calc_forces, design_solve,
lesp_design_max, transfer_cm, theodorsen, anim_flow,
drone_trajectory_problem, simple_LLT, find_tstep, ldvm_E, update_adot,
ldvm_more, ldvm_E_more, calc_forces_more, calc_forces_E_more, interp,
QScorrect_lautat, QScorrect_ldvm, QSLLT_lautat, QSLLT_ldvm,
lautat_wakeroll_more, LLT_ldvm, update_externalvel

export KinemPar, KinemParwFlap, KinemPar2DOF, KinemPar2DFree,
MotionDef, KinemDef, KinemDefwFlap, EldUpDef, EldUptstartDef,
ConstDef, TwoDOFPar, TwoDFreePar, TwoDSurf, TwoDSurfLV, TwoDSurfwFlap,
TwoDSurf_2DOF, TwoDFreeSurf, TwoDVort, TwoDFlowField, KelvinCondition,
KelvinKutta, EldRampReturnDef, EldUpIntDef, EldUpInttstartDef, SinDef,
CosDef, TheoDef, TheoDefwFlap, TwoDFlowData, DelVortDef, patch,
ThreeDSurfSimple, KinemDef3D, ThreeDVort, TwoDLVPanel, ThreeDFieldSimple,
KelvinConditionLLTldvm, KelvinKuttaLLTldvm

include("types.jl")
include("calculations.jl")
include("solvers.jl")
include("postprocess.jl")

end
