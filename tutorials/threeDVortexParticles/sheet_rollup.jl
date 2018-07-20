#===============================================================================
sheet_rollup.jl

A demo of a sheet of vortex particles rollling up.

HJA Bird 2018
h.bird.1@research.gla.ac.uk
===============================================================================#

#=--------------------------- Dependencies -----------------------------------=#
# Plotting in Julia on my PC is broken, so I have to load dependencies like
# this. If yours works, try import UNSflow instead.
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

#---------------------------- User parameters --------------------------------=#
# ODE integration parameters
const num_steps = 300
const dt = 0.05
# Options: euler_forward_step!, explicit_midpoint_step!
ode_method = euler_forward_step!
# Data saving parameters
basepath = "./output/sheet_rollup_"   # Where to write our output files
const save_every = 10                       # Save every 10 steps.

# We can make a sheet of vortex particles defined by functions. Lets make
# a flat sheet for now, because the function that does the work here
# currently only works for nearly flat stuff.
# Define a function that turns x,y to TheeDVector
function f1(x, y)
    return ThreeDVector(x, y, 0)
end
# We also want it to be bounded
bounds = [-1, 2, -0.5, 0.5] # minx, maxx, miny, maxy
# And it needs to have a known number of particles
np_x = 40
np_y = 15
# and finally we define a continious vorticy density function function.
# In theory, out of plane vorticity is not ok unless we have thickness, but
# we'll conveniently gloss over that. We define vorticity here in global coord.
function f2(x, y)
    return ThreeDVector(0.05 * y, 0.1 * (x^2)^0.25, 0.0)
end

# options: singular, planetary, exponential, winckelmans, tanh, gaussian,
# and super_gaussian
reduction_factor_fn, vorticity_fraction_fn = threed_winckelmans_kernels()

#=---------------------- Automated problem setup -----------------------------=#
particles = vortex_particle_sheet(
    f1,
    f2,
    bounds,
    np_x,
    np_y
)
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
