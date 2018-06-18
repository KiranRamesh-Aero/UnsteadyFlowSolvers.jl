# Code to test if Fourier series can capture the shape of thickness slope

ndiv = 140
naterm = 120

theta = [linspace(-pi,pi,ndiv);]
dtheta = 2*pi/(ndiv-1)

x = 0.5*(1-cos.(theta))

coord_file = "NACA0012"

m = parse(Int, coord_file[5])/100.
p = parse(Int, coord_file[6])/10.
th = parse(Int, coord_file[7:8])/100.

b1 = 0.2969; b2 = -0.1260; b3 = -0.3516; b4 = 0.2843; b5 = -0.1015
#
thick = zeros(ndiv)
thick_slope = zeros(ndiv)
#
bterm = zeros(naterm)
#b0 = zeros(2)
#
 for i = 2:ndiv
     thick[i] = 5*th*(b1*sqrt(x[i]) + b2*x[i] + b3*x[i]^2 + b4*x[i]^3 + b5*x[i]^4)
     thick_slope[i] = 5*th*(b1/(2*sqrt(x[i])) + b2 + 2*b3*x[i] + 3*b4*x[i]^2 + 4*b5*x[i]^3)
end
 thick[1] = 5*th*(b1*sqrt(x[1]) + b2*x[1] + b3*x[1]^2 + b4*x[1]^3 + b5*x[1]^4)
 thick_slope[1] = 2*thick_slope[2] - thick_slope[3]

# #Cylinder
# for i = 1:ndiv
#     thick[i] = 0.5*sin(theta[i])
# end
# thickspl = Spline1D(x, thick)
# thick_slope = derivative(thickspl, x)

#Jowkowski Airfoil
# tau = 0.12
# tau1 = 4*tau/(3*sqrt(3.))
#
# for i = 1:ndiv
#     thick[i] = 0.5*tau1*sqrt(1. - 4*(x[i] - 0.5)^2)*(1 - 2*(x[i] - 0.5))
# end
# thickspl = Spline1D(x, thick)
# thick_slope = derivative(thickspl, x)
# for i = 2:ndiv
#     thick_slope[i] = tau1*sqrt(1-x[i])*((1-2*x[i])/sqrt(x[i]) - 2*sqrt(x[i]))
# end
# thick_slope[1] = 2*thick_slope[2] - thick_slope[3]


LHS = zeros(ndiv-2, ndiv-2)
RHS = zeros(ndiv-2)

#for i = 2:ndiv-1 #TE vorticity is zero, nothing to solve for
#    for j = 2:ndiv-1
#        LHS[i,j] = plot
