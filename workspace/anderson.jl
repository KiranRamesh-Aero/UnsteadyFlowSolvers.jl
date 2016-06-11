workspace()
include("../src/UNSflow.jl")
using UNSflow

#Reproduce the case from Anderson et al.

immutable SinDef <: MotionDef
  mean :: Float64
  amp :: Float64
  w :: Float64
  phi :: Float64
end

immutable CosDef <: MotionDef
  mean :: Float64
  amp :: Float64
  w :: Float64
  phi :: Float64
end

function call(kin::SinDef, t)
  (mean*pi/180) + (kin.amp*pi/180)*sin(kin.w*t + kin.phi*pi/180)
end

function call(kin::CosDef, t)
  (mean*pi/180) + (kin.amp*pi/180)*cos(kin.w*t + kin.phi*pi/180)
end

alpha_mean = 4
alpha_amp = 22.5
k = 3.93
h_amp = 0.05
phi = 90

w = 2*k
T = (2*pi/w)
t_tot = 3*T
dt = 0.015

hdef = CosDef(0.,h_amp,w,0.)
alphadef = SinDef(alpha_mean, alpha_amp, w, phi)

udef = ConstDef(1.)
full_kinem = KinemDef(alphadef, hdef, udef)

pvt = 1./3. #leading edge

lespcrit = [0.21;]

surf = TwoDSurf(1., 1., "sd7003_fine.dat", pvt, 70, 35, "Prescribed", full_kinem, lespcrit)

curfield = TwoDFlowField()

nsteps =round(Int,t_tot/dt)+1

ldvm(surf, curfield, nsteps)
