workspace()
include("../src/UNSflow.jl")
using UNSflow


dt = 0.015;
nsteps = 1000

a1 = EldUpDef(30.*pi/180,0.2,0.8)
a2 = EldUpIntDef(0.1,0.2*0.1/(30*pi/180),0.8)
a3 = SinDef(0,20.*pi/180,1.,0.)
a4 = CosDef(0,20.*pi/180,1.,0.)
a5 = EldRampReturnDef(45*pi/180,0.2,11)

t = zeros(nsteps)
al1 = zeros(nsteps)
al2 = zeros(nsteps)
al3 = zeros(nsteps)
al4 = zeros(nsteps)
al5 = zeros(nsteps)
dal1 = zeros(nsteps)
dal2 = zeros(nsteps)
dal3 = zeros(nsteps)
dal4 = zeros(nsteps)
dal5 = zeros(nsteps)

for i = 1:nsteps
  t[i] = (i-1)*dt
  al1[i] = a1(t[i])
  dal1[i] = ForwardDiff.derivative(a1,t[i])
  al2[i] = a2(t[i])
  dal2[i] = ForwardDiff.derivative(a2,t[i])
  al3[i] = a3(t[i])
  dal3[i] = ForwardDiff.derivative(a3,t[i])
  al4[i] = a4(t[i])
  dal4[i] = ForwardDiff.derivative(a4,t[i])
  al5[i] = a5(t[i])
  dal5[i] = ForwardDiff.derivative(a5,t[i])
end
