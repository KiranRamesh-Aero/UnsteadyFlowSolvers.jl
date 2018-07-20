#===============================================================================
orbiting_particles.jl

2 vortex particles orbiting each other - as simple as it gets.

HJA Bird 2018
h.bird.1@research.gla.ac.uk
===============================================================================#


#=--------------------------- Dependencies -----------------------------------=#
# Plotting in Julia on my PC is broken, so I have to load dependencies like
# this. If yours works, try import UNSflow.
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
using WriteVTK  # We'll use this package to output to VTK for visualisation.

#-------------------------- User (that's you) parameters ---------------------=#
# ODE integration parameters
const num_steps = 2000
const dt = 0.1
# Options: euler_forward_step!, explicit_midpoint_step!
ode_method = explicit_midpoint_step!
# Data saving parameters
basepath = "./output/orbiting_particles_"   # Where to write our output files
const save_every = 10                       # Save every 10 steps.
# Initial conditions
particle1 = ThreeDVortexParticle(
    ThreeDVector(-1., 0., 0.),  # Coordinates
    ThreeDVector(0., 0., 1.),   # Vorticity
    0.1,                        # Radius
    # These don't matter, but have to be filled in anyway.
    ThreeDVector(-1., 0., 0.),  # Velocity
    ThreeDVector(0., 0., 1.)    # Rate of change of vorticity
)

particle2 = deepcopy(particle1)
particle2.coord = [1., 0., 0.]      # We can also assign to our second particle.
particle2.vorticity = [0., 0., 0.5]
particle2.size = 0.08

# options: singular, planetary, exponential, winckelmans, tanh, gaussian,
# and super_gaussian
reduction_factor_fn, vorticity_fraction_fn = threed_super_gaussian_kernels()

#=-------------------------- ODE time integration ----------------------------=#
particles = [particle1, particle2]
num_particles = size(particles)[1]
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
