type KinemDef3D
    alpha :: MotionDef
    h :: MotionDef
    u :: MotionDef

    function KinemDef3D(alpha :: MotionDef, h::MotionDef, u::MotionDef)
        new(alpha, h, u )
    end
end

immutable ThreeDFieldSimple
    f2d :: Vector{TwoDFlowField}
    function ThreeDFieldSimple()
        f2d = TwoDFlowField[]
        new(f2d)
    end
end

immutable ThreeDSurfSimple
    cref :: Float64
    AR :: Float64
    uref :: Float64
    pvt :: Float64
    lespcrit :: Vector{Float64}
    coord_file :: String
    ndiv :: Int8
    nspan :: Int8
    naterm :: Int8
    kindef :: KinemDef3D
    psi :: Vector{Float64}
    yle :: Vector{Float64}
    s2d :: Vector{TwoDSurf}
    a03d :: Vector{Float64}
    bc :: Vector{Float64}
    nshed :: Vector{Float64}
    bcoeff :: Vector{Float64}
    levstr :: Vector{Float64}
    fc :: Array{Float64}
    aterm3d :: Array{Float64}

    function ThreeDSurfSimple(AR, kindef, coord_file, pvt, lespcrit = [10.;]; nspan = 10, cref = 1., uref=1., ndiv=70, naterm=35)

        bref = AR*cref

        psi = zeros(nspan)
        yle = zeros(nspan)

        s2d = TwoDSurf[]

        for i = 1:nspan
            psi[i] = real(i)*(pi/2)/nspan
            yle[i] = -bref*cos(psi[i])/2.
        end

        #This code should be made more general to allow more motion types and combinations
        if typeof(kindef.h) == BendingDef
            for i = 1:nspan
                h_amp = evaluate(kindef.h.spl, yle[i])*kindef.h.scale
                h2d = CosDef(0., h_amp, kindef.h.k, kindef.h.phi)
                kinem2d = KinemDef(kindef.alpha, h2d, kindef.u)
                push!(s2d, TwoDSurf(coord_file, pvt,  kinem2d, lespcrit, c=cref, uref=uref, ndiv=ndiv, naterm=naterm))
            end
        else
            for i = 1:nspan
                kinem2d = KinemDef(kindef.alpha, kindef.h, kindef.u)
                lespc = lespcrit[1]
                push!(s2d, TwoDSurf(coord_file, pvt,  kinem2d, [lespc;], c=cref, uref=uref, ndiv=ndiv, naterm=naterm))
            end
        end

        a03d = zeros(nspan)
        aterm3d = zeros(naterm, nspan)

        bc = zeros(nspan)
        nshed = [0.;]
        bcoeff = zeros(nspan)
        levstr = zeros(nspan)
        fc = zeros(nspan,3)

        new(cref, AR, uref, pvt, lespcrit, coord_file,  ndiv, nspan, naterm, kindef,
        psi, yle, s2d, a03d, bc, nshed, bcoeff, levstr, fc, aterm3d)

    end
end

immutable KelvinConditionLLT
    surf :: ThreeDSurfSimple
    field :: ThreeDFieldSimple
end

function (kelv::KelvinConditionLLT)(tev_iter::Array{Float64})
    val = zeros(kelv.surf.nspan)

    #Assume symmetry condition for now
    for i = 1:kelv.surf.nspan
        kelv.field.f2d[i].tev[end].s = tev_iter[i]

        #Update incduced velocities on airfoil
        update_indbound(kelv.surf.s2d[i], kelv.field.f2d[i])

        #Calculate downwash
        update_downwash(kelv.surf.s2d[i], [kelv.field.f2d[i].u[1], kelv.field.f2d[i].w[1]])

        #Calculate first two fourier coefficients
        update_a0anda1(kelv.surf.s2d[i])

        kelv.surf.bc[i] = kelv.surf.s2d[i].a0[1] + 0.5*kelv.surf.s2d[i].aterm[1]
        end

    calc_a0a13d(kelv.surf)

    for i = 1:kelv.surf.nspan
        val[i] = kelv.surf.s2d[i].uref*kelv.surf.s2d[i].c*pi*(kelv.surf.bc[i]
        + kelv.surf.a03d[i]) + 0.5*kelv.surf.aterm3d[1,i]

        nlev = length(kelv.field.f2d[i].lev)
        ntev = length(kelv.field.f2d[i].tev)

        for iv = 1:ntev
            val[i] = val[i] + kelv.field.f2d[i].tev[iv].s
        end
        for iv = 1:nlev
            val[i] = val[i] + kelv.field.f2d[i].lev[iv].s
        end
    end

    return val
end

immutable KelvinKuttaLLT
    surf :: ThreeDSurfSimple
    field :: ThreeDFieldSimple
    nshed :: Int
end

function (kelv::KelvinKuttaLLT)(tev_iter::Array{Float64})
    val = zeros(kelv.surf.nspan + kelv.nshed)

    #Assume symmetry condition for now
    for i = 1:kelv.surf.nspan
        kelv.field.f2d[i].tev[end].s = tev_iter[i]
    end

    cntr = kelv.surf.nspan + 1
    for i = 1:kelv.surf.nspan
        if kelv.surf.s2d[i].levflag == 1
            kelv.field.f2d[i].lev[end].s = tev_iter[cntr]
            cntr += 1
        end
    end
    for i = 1:kelv.surf.nspan
        #Update incduced velocities on airfoil
        update_indbound(kelv.surf.s2d[i], kelv.field.f2d[i])

        #Calculate downwash
        update_downwash(kelv.surf.s2d[i], [kelv.field.f2d[i].u[1], kelv.field.f2d[i].w[1]])

        #Calculate first two fourier coefficients
        update_a0anda1(kelv.surf.s2d[i])

        kelv.surf.bc[i] = kelv.surf.s2d[i].a0[1] + 0.5*kelv.surf.s2d[i].aterm[1]
    end

    calc_a0a13d(kelv.surf)



    for i = 1:kelv.surf.nspan
        val[i] = kelv.surf.s2d[i].uref*kelv.surf.s2d[i].c*pi*(kelv.surf.bc[i]
        + kelv.surf.a03d[i]) + 0.5*kelv.surf.aterm3d[1,i]

        nlev = length(kelv.field.f2d[i].lev)
        ntev = length(kelv.field.f2d[i].tev)

        for iv = 1:ntev
            val[i] = val[i] + kelv.field.f2d[i].tev[iv].s
        end
        for iv = 1:nlev
            val[i] = val[i] + kelv.field.f2d[i].lev[iv].s
        end
    end

    cntr = kelv.surf.nspan + 1
    for i = 1:kelv.surf.nspan
        if kelv.surf.s2d[i].levflag == 1
            if kelv.surf.s2d[i].a0[1] > 0
                lesp_cond = kelv.surf.s2d[i].lespcrit[1]
            else
                lesp_cond = -kelv.surf.s2d[i].lespcrit[1]
            end
            val[cntr] = kelv.surf.s2d[i].a0[1] + kelv.surf.a03d[1] - lesp_cond
            cntr += 1
        end
    end

    return val
end

#===============================================================================
    ThreeDVector

    Initial code: HJAB 2018
------------------------------------------------------------------------------=#
type ThreeDVector
    x :: Float64
    y :: Float64
    z :: Float64
end

import Base.convert
""" Convert type ThreeDVector to Array{Float64}"""
function convert(::Type(Array{Float64, 1}), a::ThreeDVector)
    return [a.x, a.y, a.z]
end

""" Convert type Array{Float64} to ThreeDVector"""
function convert(::Type{ThreeDVector}, a::Array{Float64, 1})
    @assert(size(a)[1] == 3)
    b = ThreeDVector(a[1], a[2], a[3])
    return b
end

""" Convert type Array{Int64} to ThreeDVector"""
function convert(::Type{ThreeDVector}, a::Array{Float64, 1})
    @assert(size(a)[1] == 3)
    b = ThreeDVector(Float64(a[1]), Float64(a[2]), Float64(a[3]))
    return b
end

import Base.+
function +(a::ThreeDVector, b::ThreeDVector)
    c = ThreeDVector([
        a.x + b.x,
        a.y + b.y,
        a.z + b.z ])
    return c
end

import Base.-
function -(a::ThreeDVector, b::ThreeDVector)
    c = ThreeDVector([a.x - b.x, a.y - b.y, a.z - b.z])
    return c
end

import Base.*
function *(a::ThreeDVector, b::Float64)
    c = ThreeDVector([
        a.x * b,
        a.y * b,
        a.z * b ])
    return c
end

function *(a::Float64, b::ThreeDVector)
    return b * a
end

import Base./
function /(a::ThreeDVector, b::Float64)
    c = ThreeDVector([
        a.x / b,
        a.y / b,
        a.z / b ])
    return c
end

"""Cross product of two ThreeDVectors"""
function cross(a::ThreeDVector, b::ThreeDVector)
    c = ThreeDVector(
        [a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x ])
    return c
end

"""Dot product of two ThreeDVectors"""
function dot(a::ThreeDVector, b::ThreeDVector)
    return a.x * b.x + a.y * b.y + a.z * b.z
end

function abs(a::ThreeDVector)
    return sqrt(a.x^2 + a.y^2 + a.z^2)
end

""" Set vector a to zero"""
function zero!(a::ThreeDVector)
    a.x = 0.0
    a.y = 0.0
    a.z = 0.0
    return Void
end

"""Return vector of length 1 with same direction"""
function unit(a::ThreeDVector)
    b = abs(a)
    return a / b
end

"""Set a vector normalised to length 1 with same direction"""
function unit!(a::ThreeDVector)
    b = abs(a)
    a.x /= b
    a.y /= b
    a.z /= b
    return Void
end

import Base.getindex
function getindex(a::ThreeDVector, i::Int64)
    if i == 1 return a.x
    elseif i == 2 return a.y
    elseif i == 3 return a.z
    else throw(Core.BoundsError)
    end
end

function size(a::ThreeDVector)
    return [3]
end

#= END ThreeDVector ----------------------------------------------------------=#

#===============================================================================
    ThreeDVortexParticle

    Initial code: HJAB 2018
------------------------------------------------------------------------------=#
type ThreeDVortexParticle
    coord :: ThreeDVector
    vorticity :: ThreeDVector
    size :: Float64

    velocity :: ThreeDVector
    vorticity_time_derivative :: ThreeDVector
end

"""
According to input vortex particle field particles, and kernal functions g and
f, compute the change for a timestep dt and return this as a new
particles_updated array of particles.
"""
function one_step(
    particles::Array{ThreeDVortexParticle},
    dt::Float64,
    g_function::Function,
    f_function::Function
    )

    delta_x = Array{ThreeDVector, 1}(size(particles))
    delta_vort = Array{ThreeDVector, 1}(size(particles))
    particles_updated = deepcopy(particles)

    function get_dx(particle_idx::Int64, particles::Array{ThreeDVortexParticle})
        location = particles[particle_idx].coord
        v = map(x->induced_velocity(x, location, g_function),
                        particles[vcat(1:particle_idx-1, particle_idx+1:end)])
        dx = sum(v) * dt
        return dx
    end
    function get_dvort(particle_idx::Int64, particles::Array{ThreeDVortexParticle})
        location = particles[particle_idx].coord
        vo = map(x->rate_of_change_of_vorticity(x, particles[particle_idx], g_function, f_function),
                    particles[vcat(1:particle_idx-1, particle_idx+1:end)])
        dvo = sum(vo) * dt
        return dvo
    end

    delta_x = map(i->get_dx(i, particles), 1:size(particles)[1])
    delta_vort = map(i->get_dvort(i, particles), 1:size(particles)[1])

    for i = 1 : size(particles)[1]
        particles_updated[i].coord += delta_x[i]
        particles_updated[i].vorticity += delta_vort[i]
    end
    return particles_updated
end

"""
Compute the velocity induced by ThreeDParticle particle at ThreeDCoordinate
coordinate, given a g function.
"""
function induced_velocity(
    particle::ThreeDVortexParticle,
    coordinate::ThreeDVector,
    g_function::Function
     )
     # Robertson 2010 Eq. 3
     rad = particle.coord - coordinate
     a = g_function(abs(rad)/particle.size) / (4. * pi)
     den = abs(rad)^3
     c = cross(rad, particle.vorticity)
     vel :: ThreeDVector = c * a / den
     return vel
end

"""
The vorticity of a particle over time changes as vortices stretch
 and what not. This RETURNS (doesn't change the value of) the rate of change
 of particle j with respect to time as domega_x / dt, domega_y / dt,
 domega_z / dt
"""
function rate_of_change_of_vorticity(
    particle_j::ThreeDVortexParticle,
    particle_k::ThreeDVortexParticle,
    g_function::Function,
    f_function::Function
    )
    # Robertson 2010 Eq. 4
    sigma_k = particle_k.size
    om_j = particle_j.vorticity
    om_k = particle_k.vorticity
    r = particle_j.coord - particle_k.coord
    g = g_function(abs(r) / particle_k.size)
    f = f_function(abs(r) / particle_k.size)
    om_cross = cross(particle_j.vorticity, particle_k.vorticity)
    r_cross_k = cross(r, particle_k.vorticity)
    om_dot_r = dot(particle_j.vorticity, r)

    t1 = 1. / (4 * pi * sigma_k ^ 3)
    t21n = - g * om_cross
    t21d = (abs(r) / sigma_k) ^ 3
    t21 = t21n / t21d
    t221 = 1 / (abs(r)^2)
    t222 = 3 * g / (abs(r) / sigma_k) ^ 3 - f
    t223 = om_dot_r * r_cross_k
    t22 = t221 * t222 * t223
    t2 = t21 + t22
    t = t1 * t2
    return t
end

"""Take a set of vortex particles and the value of g for each of them,
 and computes the velocity of each. Returns a vector of velocities. """
function particle_velocities(
    particles::Array{ThreeDVortexParticle},
    particle_g_function::Function)

    @assert(particles.size() == particles_g_values.size())  # Same size arrays

    vel = zeros(particles.size())
    i :: Int64
    v0 :: ThreeDVector
    for i = 1 : particles.size()
        x0 = particles[i].coord
        zero!(v0)
        vel[i] = mapreduce(x->induced_velocity(x, x0, particle_g_function),
            +, v0, particles)
    end
    return vel
end

"""Take a set of vortex particles and the value of g and f for each of them,
 and computes the rate of change of vorticity for each of them. """
function particle_rate_of_change_of_vorticity(
    particles::Array{ThreeDVortexParticle},
    particle_g_function::Function,
    particle_f_function::Function)

    @assert(particles.size() == particles_g_values.size())  # Same size arrays

    dvort = zeros(particles.size())
    i :: Int64
    v0 :: ThreeDVector
    for i = 1 : particles.size()
        x0 = particles[i].coord
        zero!(v0)
        dvort[i] = mapreduce(
            x->rate_of_change_of_vorticity(dvort[i], x[1],
                                        particle_g_values, particle_f_values),
            +, v0, particles)
    end
    return vel
end
#= END ThreeDVortexParticle --------------------------------------------------=#
