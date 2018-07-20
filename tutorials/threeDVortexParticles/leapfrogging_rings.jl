#===============================================================================
leapfrogging_rings.jl

A demo of leapfrogging vortex rings using vortex particles.

HJA Bird 2018
h.bird.1@research.gla.ac.uk
===============================================================================#

#=--------------------------- Dependencies -----------------------------------=#
using Dierckx
using SpecialFunctions
import Base.size
include("../../src/kinem.jl")
include("../../src/utils.jl")
include("../../src/delVort.jl")
include("../../src/lowOrder2D/typedefs.jl")
include("../../src/lowOrder2D/calcs.jl")
include("../../src/lowOrder3D/typedefs.jl")
include("../../src/lowOrder3D/calcs.jl")
include("ThreeDVortexParticleFlowFeatures.jl")
using WriteVTK  # We'll use this package to output to VTK for visualisation.

#-------------------------- User (that's you) parameters ---------------------=#
# ODE integration parameters
const num_steps = 2
const dt = 0.1
# Options: euler_forward_step!, explicit_midpoint_step!
ode_method = euler_forward_step!
# Data saving parameters
basepath = "./output/leapfrogging_rings_"   # Where to write our output files
const save_every = 10                       # Save every 10 steps.
# Initial conditions
const ring_strength = [1., 1.]
const ring_particles = [3, 3]
const ring_radii = [1., 1.]
const ring_locations = [0., 1.]
# options: singular, planetary, exponential, winckelmans, tanh, gaussian,
# and super_gaussian
reduction_factor_fn, vorticity_fraction_fn = threed_winckelmans_kernels()

#=---------------------- Automated problem setup -----------------------------=#
num_rings = size(ring_particles)[1]
particles = Vector{ThreeDVortexParticle}()
for i = 1 : num_rings
    centre = ThreeDVector(ring_locations[i], 0., 0.)
    normal = ThreeDVector(1., 0., 0.)
    radius = ring_radii[i]
    strength = ring_strength[i]
    n_ring_particles = ring_particles[i]
    ring = vortex_particle_ring(
        centre, normal, radius, strength, n_ring_particles
        )
    particles = [particles; ring]
end
num_particles = size(particles)[1]

#=-------------------------- ODE time integration ----------------------------=#
for i = 1 : num_steps
    # Save the current state to vtk if required
    if (i - 1) % save_every == 0
        points = zeros(3, num_particles)
        point_vorticity = zeros(3, num_particles)
        cells = Array{MeshCell, 1}(num_particles)
        for j = 1 : num_particles
            points[:,j] = [particles[j].coord.x,
                particles[j].coord.y, particles[j].coord.z]
            point_vorticity[:, j] = [particles[j].vorticity.x,
                particles[j].vorticity.y, particles[j].vorticity.z]
            cells[j] = MeshCell(VTKCellTypes.VTK_VERTEX, [j])
        end
        vtkfile = vtk_grid(string(basepath, i), points, cells)
        vtk_point_data(vtkfile, point_vorticity, "vorticity")
        outfiles = vtk_save(vtkfile)
    end

    # Calculate the next iteration
    particles = ode_method(particles, dt,
        reduction_factor_fn, vorticity_fraction_fn)
end

#=------------------- Now repeat until it doesn't blow up --------------------=#
