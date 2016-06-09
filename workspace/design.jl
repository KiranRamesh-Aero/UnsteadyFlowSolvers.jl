workspace()
include("../src/UNSflow.jl")
using UNSflow

immutable EldUpIntDef <: MotionDef
    amp :: Float64
    K :: Float64
    a :: Float64
end

function call(eld::EldUpIntDef, t)

    dt = 0.015
    sm = pi*pi*eld.K/(2*(30*pi/180)*(1 - eld.a))
    t1 = 1.
    t2 = t1 + ((30*pi/180)/(2*eld.K))

    tmpt = 0
    prev_h = 0
    amp = 0
    while tmpt <= t
      hdot = ((eld.K/sm)*log(cosh(sm*(t - t1))/cosh(sm*(t - t2))))+(30*pi/360)
      hdot = hdot*eld.h_amp/(30*pi/180)
      amp = prev_h + hdot*dt
      prev_h = amp
      tmpt = tmpt + dt
    end
    amp
end

function lesp_design_max!(h_amp,lesp_diff)
  println("ha",h_amp)
  alphadef = EldUpDef(30,0.2,0.8)
  hdef = EldUpIntDef(h_amp[1],0.2,0.8)
  udef = ConstDef(1.)
  full_kinem = KinemDef(alphadef, hdef, udef)
  lespcrit = [20;]
  pvt = 0.25

  surf = TwoDSurf(1., 1., "sd7003_fine.dat", pvt, 70, 35, "Prescribed", full_kinem,lespcrit)

  curfield = TwoDFlowField()

  nsteps =round(Int,3/0.015)+1

  ldvm(surf, curfield, nsteps)

  data =  readdlm("results.dat")

  lesp_diff = maximum(data[:,5]) - 0.21
  println("ld",lesp_diff)
end
