abstract type MotionDef end

mutable struct KinemPar
    alpha :: Float64
    h :: Float64
    alphadot :: Float64
    hdot :: Float64
    u :: Float64
    udot :: Float64
end

mutable struct KinemDef
    alpha :: MotionDef
    h :: MotionDef
    u :: MotionDef
end

struct EldUpDef <: MotionDef
    amp :: Float64
    K :: Float64
    a :: Float64
end

function (eld::EldUpDef)(t)
    sm = pi*pi*eld.K/(2*(eld.amp)*(1 - eld.a))
    t1 = 1.
    t2 = t1 + ((eld.amp)/(2*eld.K))
    ((eld.K/sm)*log(cosh(sm*(t - t1))/cosh(sm*(t - t2))))+(eld.amp/2)
end


struct EldUptstartDef <: MotionDef
    amp :: Float64
    K :: Float64
    a :: Float64
    tstart :: Float64
end

function (eld::EldUptstartDef)(t)
    sm = pi*pi*eld.K/(2*(eld.amp)*(1 - eld.a))
    t1 = eld.tstart
    t2 = t1 + ((eld.amp)/(2*eld.K))
    ((eld.K/sm)*log(cosh(sm*(t - t1))/cosh(sm*(t - t2))))+(eld.amp/2)
end


struct EldRampReturnDef <: MotionDef
    amp :: Float64
    K :: Float64
    a :: Float64
end

function (eld::EldRampReturnDef)(tt)
    fr = eld.K/(pi*abs(eld.amp));
    t1 = 1.
    t2 = t1 + (1. /(2*pi*fr));
    t3 = t2 + ((1/(4*fr)) - (1/(2*pi*fr)));
    t4 = t3 + (1. /(2*pi*fr));
    t5 = t4+1.;

    nstep = round(Int,t5/0.015) + 1
    g = zeros(nstep)
    t = zeros(nstep)

    for i = 1:nstep
        t[i] = (i-1.)*0.015
        g[i] = log((cosh(eld.a*(t[i] - t1))*cosh(eld.a*(t[i] - t4)))/(cosh(eld.a*(t[i] - t2))*cosh(eld.a*(t[i] - t3))))
    end
    maxg = maximum(g);

    gg = log((cosh(eld.a*(tt - t1))*cosh(eld.a*(tt - t4)))/(cosh(eld.a*(tt - t2))*cosh(eld.a*(tt - t3))))

    return eld.amp*gg/(maxg);

end

struct EldRampReturntstartDef <: MotionDef
    amp :: Float64
    K :: Float64
    a :: Float64
    tstart :: Float64
end

function (eld::EldRampReturntstartDef)(tt)
    fr = eld.K/(pi*abs(eld.amp));
    t1 = eld.tstart
    t2 = t1 + (1. /(2*pi*fr));
    t3 = t2 + ((1/(4*fr)) - (1/(2*pi*fr)));
    t4 = t3 + (1. /(2*pi*fr));
    t5 = t4+t1;

    nstep = round(Int,t5/0.015) + 1
    g = zeros(nstep)
    t = zeros(nstep)

    for i = 1:nstep
        t[i] = (i-1.)*0.015
        g[i] = log((cosh(eld.a*(t[i] - t1))*cosh(eld.a*(t[i] - t4)))/(cosh(eld.a*(t[i] - t2))*cosh(eld.a*(t[i] - t3))))
    end
    maxg = maximum(g);

    gg = log((cosh(eld.a*(tt - t1))*cosh(eld.a*(tt - t4)))/(cosh(eld.a*(tt - t2))*cosh(eld.a*(tt - t3))))

    return eld.amp*gg/(maxg);

end

struct ConstDef <: MotionDef
    amp :: Float64
end

function (cons::ConstDef)(t)
    cons.amp
end


struct LinearDef <: MotionDef
    tstart :: Float64
    vstart :: Float64
    vend :: Float64
    len :: Float64
end

function (lin::LinearDef)(t)
    if t < lin.tstart
        lin.vstart
    elseif t > lin.tstart + lin.len
        lin.vend
    else
        lin.vstart + (lin.vend - lin.vstart)/lin.len*(t - lin.tstart)
    end
end

struct BendingDef <: MotionDef
    spl :: Spline1D
    scale :: Float64
    k :: Float64
    phi :: Float64
end


struct SinDef <: MotionDef
  mean :: Float64
  amp :: Float64
  k :: Float64
  phi :: Float64
end

function (kin::SinDef)(t)
  (kin.mean) + (kin.amp)*sin(2*kin.k*t + kin.phi)
end


struct CosDef <: MotionDef
  mean :: Float64
  amp :: Float64
  k :: Float64
  phi :: Float64
end

function (kin::CosDef)(t)
  (kin.mean) + (kin.amp)*cos(2*kin.k*t + kin.phi)
end

struct StepGustDef <: MotionDef
    amp :: Float64
    tstart :: Float64
    tgust :: Float64
end
 function (sgust::StepGustDef)(t)
    if t >= sgust.tstart && t <= sgust.tstart+sgust.tgust
        amp = (sgust.amp)
    else
        amp = 0.
    end
    amp
end

struct EldUpIntDef <: MotionDef
    amp :: Float64
    K :: Float64
    a :: Float64
end

function (eld::EldUpIntDef)(t)

    dt = 0.015
    nsteps = t/dt + 1
    nsteps = round(Int,nsteps)
    dt = t/(nsteps-1)
    sm = pi*pi*eld.K/(2*eld.amp*(1 - eld.a))
    t1 = 1.
    t2 = t1 + ((eld.amp)/(2*eld.K))


    prev_h = 0
    amp = 0
    for i = 1:nsteps
      tmpt = (i-1)*dt
      if (eld.amp == 0.)
      	 hdot = 0.
      else
         hdot = ((eld.K/sm)*log(cosh(sm*(tmpt - t1))/cosh(sm*(tmpt - t2))))+(eld.amp/2.)
      end
      amp = prev_h + hdot*dt
      prev_h = amp
    end
    if (nsteps == 1)
      amp = 0.
    end
    amp
end


struct EldUpInttstartDef <: MotionDef
    amp :: Float64
    K :: Float64
    a :: Float64
    tstart :: Float64
end

function (eld::EldUpInttstartDef)(t)

    dt = 0.015
    nsteps = t/dt + 1
    nsteps = round(Int,nsteps)
    dt = t/(nsteps-1)
    sm = pi*pi*eld.K/(2*eld.amp*(1 - eld.a))
    t1 = eld.tstart
    t2 = t1 + ((eld.amp)/(2*eld.K))


    prev_h = 0
    amp = 0
    for i = 1:nsteps
      tmpt = (i-1)*dt
	  if (eld.amp == 0.)
      	 hdot = 0.
      else
         hdot = ((eld.K/sm)*log(cosh(sm*(tmpt - t1))/cosh(sm*(tmpt - t2))))+(eld.amp/2.)
      end
      amp = prev_h + hdot*dt
      prev_h = amp
    end
    if (nsteps == 1)
      amp = 0.
    end
    amp
end

# Takes a vector of time-ordered data from an input file
struct FileDef <: MotionDef
    tvec :: Vector{Float64}
    var :: Vector{Float64}
end

function (file::FileDef)(t)
    # Find index of closest t to tvec
    idx = argmin(abs.(file.tvec.- t))
    # return var at index
    file.var[idx]
end
