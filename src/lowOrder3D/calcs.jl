#===============================================================================
calcs.jl

Calculations for lowOrder3D methods.
===============================================================================#

#===============================================================================
    ThreeDSurfSimple methods
------------------------------------------------------------------------------=#
function calc_a0a13d(surf::ThreeDSurfSimple)

    lhs = zeros(surf.nspan, surf.nspan)
    rhs = zeros(surf.nspan)


    for i = 1:surf.nspan
        integ = simpleTrapz(surf.s2d[i].cam_slope, surf.s2d[i].theta)
        -simpleTrapz(surf.s2d[i].cam_slope.*cos.(surf.s2d[i].theta), surf.s2d[i].theta)

        for n = 1:surf.nspan
            nn = 2*n - 1
            lhs[i,n] = sin(nn*surf.psi[i])*(sin(surf.psi[i]) + (nn*pi/(2*surf.AR))
            *(cos(surf.s2d[i].kinem.alpha) + sin(surf.s2d[i].kinem.alpha)*integ/pi))

        end
        rhs[i] = pi*sin(surf.psi[i])*surf.bc[i]/(2*surf.AR)
    end

    surf.bcoeff[:] = lhs \ rhs

    for i = 1:surf.nspan
        surf.a03d[i] = 0
        surf.aterm3d[1,i] = 0
        integ0 = simpleTrapz(surf.s2d[i].cam_slope, surf.s2d[i].theta)
        integ1 = simpleTrapz(surf.s2d[i].cam_slope.*cos.(surf.s2d[i].theta), surf.s2d[i].theta)

        for n = 1:surf.nspan
            nn = 2*n - 1
            surf.a03d[i] -= real(nn)*surf.bcoeff[n]*sin(nn*surf.psi[i])/sin(surf.psi[i])*(cos(surf.s2d[i].kinem.alpha) + sin(surf.s2d[i].kinem.alpha)*integ0/pi)

            surf.aterm3d[1,i] += 2*real(nn)*surf.bcoeff[n]*sin(nn*surf.psi[i])/(sin(surf.psi[i])*pi)*sin(surf.s2d[i].kinem.alpha)*integ1


        end
    end
end

function calc_a2toan3d(surf::ThreeDSurfSimple)
    for ia = 2:surf.naterm
        for i = 1:surf.nspan
            surf.aterm3d[ia,i] = 0
            integ = simpleTrapz(surf.s2d[i].cam_slope.*cos.(ia*surf.s2d[i].theta), surf.s2d[i].theta)

            for n = 1:surf.nspan
                nn = 2*n - 1
                surf.aterm3d[ia,i] += 2*real(nn)*surf.bcoeff[n]*sin(nn*surf.psi[i])/(sin(surf.psi[i])*pi)*sin(surf.s2d[i].kinem.alpha)*integ
            end
        end
    end
end
 # END ThreeDSurfSimple methods ===============================================#

#===============================================================================
    ThreeDVortexParticle methods
------------------------------------------------------------------------------=#
"""The induced velocity at a point 'coordinate' due to a vortex particle.
reduction_factor_fn is the function g that avoids velocity singularities."""
function ind_vel(
    particle::ThreeDVortexParticle,
    coordinate::ThreeDVector,
    reduction_factor_fn::Function
     )
     if particle.coord == coordinate
         return ThreeDVector(0., 0., 0.)
     end
     # Robertson 2010 Eq. 3 with additional - sign.
     rad = particle.coord - coordinate
     a = reduction_factor_fn(abs(rad)/particle.size) / (4. * pi)
     den = abs(rad)^3
     c = cross(rad, particle.vorticity)
     vel :: ThreeDVector = c * a / den
     return vel
end

function ind_vel(
    particles :: Vector{ThreeDVortexParticle},
    coordinate :: ThreeDVector,
    reduction_factor_fn :: Function
     )
     z = ThreeDVector(0., 0., 0.)
     v = mapreduce(x->ind_vel(x, coordinate, reduction_factor_fn),
        +, z, particles)
     return v
 end

 function ind_vel(
     particles :: Vector{ThreeDVortexParticle},
     reduction_factor_fn :: Function
      )
      v = map(x->ind_vel(particles, x.coord, reduction_factor_fn), particles)
      return v
  end

""" The time rate of change of vorticity of a particle
"dvort_induced_on_particle" due to another particle "inducing particle" """
function ind_dvortdt(
    dvort_induced_on_particle :: ThreeDVortexParticle,
    inducing_particle :: ThreeDVortexParticle,
    reduction_factor_fn :: Function,
    vorticity_fraction_fn :: Function
    )
    if dvort_induced_on_particle.coord == inducing_particle.coord
        return ThreeDVector(0., 0., 0.)
    end

    # Robertson 2010 Eq. 4
    sigma_k = inducing_particle.size
    om_j = dvort_induced_on_particle.vorticity
    om_k = inducing_particle.vorticity
    r = dvort_induced_on_particle.coord - inducing_particle.coord
    g = reduction_factor_fn(abs(r) / inducing_particle.size)
    f = vorticity_fraction_fn(abs(r) / inducing_particle.size)
    om_cross = cross(dvort_induced_on_particle.vorticity,
                            inducing_particle.vorticity)
    r_cross_k = cross(r, inducing_particle.vorticity)
    om_dot_r = dot(dvort_induced_on_particle.vorticity, r)

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

function ind_dvortdt(
    dvort_induced_on_particle :: ThreeDVortexParticle,
    inducing_particle :: Vector{ThreeDVortexParticle},
    reduction_factor_fn :: Function,
    vorticity_fraction_fn :: Function
    )
    z = ThreeDVector(0., 0., 0.)
    v = mapreduce(x->ind_dvortdt(dvort_induced_on_particle, x,
        reduction_factor_fn, vorticity_fraction_fn), +, z, particles)
    return v
end

function ind_dvortdt(
    particles :: Vector{ThreeDVortexParticle},
    reduction_factor_fn :: Function,
    vorticity_fraction_fn :: Function
    )
    v = map(x->ind_dvortdt(x, particles,
        reduction_factor_fn, vorticity_fraction_fn), particles)
    return v
end

"""
Compute the velocity and rate of change of vorticity of a set of vortex
particles on itself.
"""
function mutual_ind(
    particles::Vector{ThreeDVortexParticle},
    reduction_factor_fn :: Function,
    vorticity_fraction_fn :: Function
    )

    vel = ind_vel(particles, reduction_factor_fn)
    dvortdt = ind_dvortdt(particles, reduction_factor_fn, vorticity_fraction_fn)
    for (particle, v, dvort) in zip(particles, vel, dvortdt)
        particle.velocity = v
        particle.vorticity_time_derivative = dvort
    end
    return particles
end

"""
Integrate the particle velocities and rate of change of vorticities over time
dt to compute new coordinate and vorticity. Assumes derivatives are already
computed using mutual_ind.
"""
function integration_step!(
    particles :: Vector{ThreeDVortexParticle},
    delta_t :: Float64
    )
    for particle in particles
        particle.coord += particle.velocity * dt
        particle.vorticity += particle.vorticity_time_derivative * dt
    end
    return particles
end

"""
Use the forward Euler method to calculate the location of vortex particles
at time dt in the particles future. Updates particles.
"""
function euler_forward_step!(
    particles :: Vector{ThreeDVortexParticle},
    delta_t :: Float64,
    reduction_factor_fn :: Function,
    vorticity_fraction_fn :: Function
    )

    mutual_ind(particles, reduction_factor_fn, vorticity_fraction_fn)
    integration_step!(particles, dt)
    return particles
end

"""
Use the explicit midpoint method to calculate the location of vortex particles
at time dt in the particles future. Updates particles.
"""
function explicit_midpoint_step!(
    particles :: Vector{ThreeDVortexParticle},
    delta_t :: Float64,
    reduction_factor_fn :: Function,
    vorticity_fraction_fn :: Function
    )

    mutual_ind(particles, reduction_factor_fn, vorticity_fraction_fn)
    particles_cpy = deepcopy(particles)
    integration_step!(particles_cpy, dt / 2.)
    mutual_ind(particles_cpy, reduction_factor_fn, vorticity_fraction_fn)
    for (particle, particlecpy) in zip(particles, particles_cpy)
        particle.velocity = particlecpy.velocity
        particle.vorticity_time_derivative =
            particlecpy.vorticity_time_derivative
    end
    integration_step!(particles, dt)
    return particles
end
 # END ThreeDVortexParticle methods ===========================================#


#===============================================================================
     ThreeDVortexParticle kernels functions

     These methods return the reduction_factor_fn and vorticity_fraction_fn
     used by ThreeDVortexParticle in velocity etc. calculations.

     These are taken from  Robertson, Joo and Reich: "Vortex Particle
     Aerodynamic Modelling of Perching Manoeuvres with Micro Air Vehicles",
     51st AIAA/ASME/ASCE/AHS/ASC Structures, Structural Dynamics, and
     Materials Conference
------------------------------------------------------------------------------=#
"""
Returns (reduction_factor_fn :: Function, vorticity_fraction_fn :: Function)
for a singular vortex particle. Used by vortex particle functions such
as mutual_ind, ind_vel and ind_dvortdt.

These functions describle vortex particle interaction.
See Robertson, Joo and Reich: "Vortex Particle
Aerodynamic Modelling of Perching Manoeuvres with Micro Air Vehicles",
51st AIAA/ASME/ASCE/AHS/ASC Structures, Structural Dynamics, and
Materials Conference for more detail on the mathematical detail implemented.
"""
function threed_singularity_kernels()
    function vortex_f_singularity(rho::Float64)
        return 0.0
    end
    function vortex_g_singularity(rho::Float64)
        return 1.0
    end
    return vortex_g_singularity, vortex_f_singularity
end

"""
Returns (reduction_factor_fn :: Function, vorticity_fraction_fn :: Function)
for a planetary vortex particle. Used by vortex particle functions such
as mutual_ind, ind_vel and ind_dvortdt.

These functions describle vortex particle interaction.
See Robertson, Joo and Reich: "Vortex Particle
Aerodynamic Modelling of Perching Manoeuvres with Micro Air Vehicles",
51st AIAA/ASME/ASCE/AHS/ASC Structures, Structural Dynamics, and
Materials Conference for more detail on the mathematical detail implemented.
"""
function threed_planetary_kernels()
    function vortex_f_planetary(rho::Float64)
        return rho <  1 ? 3 : 0
    end
    function vortex_g_planetary(rho::Float64)
        return (min(1, rho)) ^ 3
    end
    return vortex_g_planetary, vortex_f_planetary
end

"""
Returns (reduction_factor_fn :: Function, vorticity_fraction_fn :: Function)
for a exponential vortex particle. Used by vortex particle functions such
as mutual_ind, ind_vel and ind_dvortdt.

These functions describle vortex particle interaction.
See Robertson, Joo and Reich: "Vortex Particle
Aerodynamic Modelling of Perching Manoeuvres with Micro Air Vehicles",
51st AIAA/ASME/ASCE/AHS/ASC Structures, Structural Dynamics, and
Materials Conference for more detail on the mathematical detail implemented.
"""
function threed_exponential_kernels()
    function vortex_f_exponential(rho::Float64)
        return 3.0 * exp(- rho ^ 3)
    end
    function vortex_g_exponential(rho::Float64)
        return 1 - vortex_f_exponential(rho) / 3.0
    end
    return vortex_g_exponential, vortex_f_exponential
end

"""
Returns (reduction_factor_fn :: Function, vorticity_fraction_fn :: Function)
for a winkelmans and Leonard vortex particle. Used by vortex particle functions
such as mutual_ind, ind_vel and ind_dvortdt.

These functions describle vortex particle interaction.
See Robertson, Joo and Reich: "Vortex Particle
Aerodynamic Modelling of Perching Manoeuvres with Micro Air Vehicles",
51st AIAA/ASME/ASCE/AHS/ASC Structures, Structural Dynamics, and
Materials Conference for more detail on the mathematical detail implemented.
"""
function threed_winckelmans_kernels()
    function vortex_f_winckelmans(rho::Float64)
        a = 15. / 2.
        b = rho ^ 2 + 1
        return a / (b ^ (7./2.))
    end
    function vortex_g_winckelmans(rho::Float64)
        a = rho^2 + 2.5
        b = rho^3
        c = rho^2 + 1
        d = a * b
        e = c ^(5./2.)
        return d / e
    end
    return  vortex_g_winckelmans, vortex_f_winckelmans
end

"""
Returns (reduction_factor_fn :: Function, vorticity_fraction_fn :: Function)
for a hyperbolic tangent vortex particle. Used by vortex particle functions such
as mutual_ind, ind_vel and ind_dvortdt.

These functions describle vortex particle interaction.
See Robertson, Joo and Reich: "Vortex Particle
Aerodynamic Modelling of Perching Manoeuvres with Micro Air Vehicles",
51st AIAA/ASME/ASCE/AHS/ASC Structures, Structural Dynamics, and
Materials Conference for more detail on the mathematical detail implemented.
"""
function threed_tanh_kernels()
    function vortex_f_tanh(rho::Float64)
        return 3 * sech(rho^3)^2
    end
    function vortex_g_tanh(rho::Float64)
        return tanh(rho^3)
    end
    return  vortex_g_tanh, vortex_f_tanh
end

"""
Returns (reduction_factor_fn :: Function, vorticity_fraction_fn :: Function)
for a Gaussian vortex particle. Used by vortex particle functions such
as mutual_ind, ind_vel and ind_dvortdt.

These functions describle vortex particle interaction.
See Robertson, Joo and Reich: "Vortex Particle
Aerodynamic Modelling of Perching Manoeuvres with Micro Air Vehicles",
51st AIAA/ASME/ASCE/AHS/ASC Structures, Structural Dynamics, and
Materials Conference for more detail on the mathematical detail implemented.
"""
function threed_gaussian_kernels()
    function vortex_f_gaussian(rho::Float64)
        return sqrt(2 / pi) * exp((- rho ^2) / 2)
    end
    function vortex_g_gaussian(rho::Float64)
        return erf(rho / sqrt(2.)) - rho * vortex_f_gaussian(rho)
    end
    return  vortex_g_gaussian, vortex_f_gaussian
end

"""
Returns (reduction_factor_fn :: Function, vorticity_fraction_fn :: Function)
for a super Gaussian vortex particle. Used by vortex particle functions such
as mutual_ind, ind_vel and ind_dvortdt.

These functions describle vortex particle interaction.
See Robertson, Joo and Reich: "Vortex Particle
Aerodynamic Modelling of Perching Manoeuvres with Micro Air Vehicles",
51st AIAA/ASME/ASCE/AHS/ASC Structures, Structural Dynamics, and
Materials Conference for more detail on the mathematical detail implemented.
"""
function threed_super_gaussian_kernels()
    function vortex_f_super_gaussian(rho::Float64)
        return sqrt(2 / pi) * (2.5 - rho^2 / 2) * exp((- rho ^2) / 2)
    end
    function vortex_g_super_gaussian(rho::Float64)
        return erf(rho / sqrt(2.)) - ((2 - rho^2) / (5 - rho^2)) *
            vortex_f_super_gaussian(rho)
    end
    return  vortex_g_super_gaussian, vortex_f_super_gaussian
end
  # END ThreeDVortexParticle kernels functions ================================#
