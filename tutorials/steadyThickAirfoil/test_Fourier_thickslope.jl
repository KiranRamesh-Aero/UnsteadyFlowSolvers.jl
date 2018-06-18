# Code to test if Fourier series can capture the shape of thickness slope

ndiv = 140
naterm = 120

theta = [linspace(0,pi,ndiv);]

x = 0.5*(1-cos.(theta))

 coord_file = "NACA0012"

 m = parse(Int, coord_file[5])/100.
 p = parse(Int, coord_file[6])/10.
 th = parse(Int, coord_file[7:8])/100.
#
 b1 = 0.2969; b2 = -0.1260; b3 = -0.3516; b4 = 0.2843; b5 = -0.1015
#
thick = zeros(ndiv)
thick_slope = zeros(ndiv)
#
bterm = zeros(naterm)
b0 = zeros(2)
#
 for i = 2:ndiv
     thick[i] = 5*th*(b1*sqrt(x[i]) + b2*x[i] + b3*x[i]^2 + b4*x[i]^3 + b5*x[i]^4)
     thick_slope[i] = 5*th*(b1/(2*sqrt(x[i])) + b2 + 2*b3*x[i] + 3*b4*x[i]^2 + 4*b5*x[i]^3)
end
 thick[1] = 5*th*(b1*sqrt(x[1]) + b2*x[1] + b3*x[1]^2 + b4*x[1]^3 + b5*x[1]^4)
 thick_slope[1] = 2*thick_slope[2] - thick_slope[3]
 rho = 1.1019*th*th


# #Cylinder
# for i = 1:ndiv
#     thick[i] = 0.5*sin(theta[i])
# end
# thickspl = Spline1D(x, thick)
# thick_slope = derivative(thickspl, x)

#Jowkowski Airfoil
# tau = 0.12
# tau1 = 4*tau/(3*sqrt(3.))
# rho = 2*tau1*tau1
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

LHS = zeros(ndiv-2, naterm+2)
RHS = zeros(ndiv-2)

for i = 2:ndiv-1
    for n = 1:naterm
        LHS[i-1,n] = sin(n*theta[i])
    end
    LHS[i-1,naterm+1] = (1+cos(theta[i]))/sin(theta[i])
    LHS[i-1,naterm+2] = (1-cos(theta[i]))/sin(theta[i])
    RHS[i-1] = thick_slope[i]
end

soln = LHS \ RHS
bterm = soln[1:end-2]
b0 = soln[end-1:end]
u1t = zeros(ndiv)
w1t = zeros(ndiv)
for i = 2:ndiv-1
    w1t[i] = b0[1]*(1+cos(theta[i]))/sin(theta[i])
    w1t[i] += b0[2]*(1-cos(theta[i]))/sin(theta[i])
    u1t[i] = b0[1] - b0[2]
    for n = 1:naterm
        u1t[i] -= bterm[n]*cos(n*theta[i])
        w1t[i] += bterm[n]*sin(n*theta[i])
    end
end
w1t[1] = 2*w1t[2] - w1t[3]
w1t[ndiv] = 2*w1t[ndiv-1] - w1t[ndiv-2]
u1t[1] = 2*u1t[2] - u1t[3]
u1t[ndiv] = 2*u1t[ndiv-1] - u1t[ndiv-2]

t1 = zeros(ndiv)
for i = 1:ndiv
    t1[i] = thick[i]*u1t[i]
end
t1spl = Spline1D(x, t1,k=1)
dt1 = derivative(t1spl,x)

LHS = zeros(ndiv-2, naterm+2)
RHS = zeros(ndiv-2)

for i = 2:ndiv-1
    for n = 1:naterm
        LHS[i-1,n] = sin(n*theta[i])
    end
    LHS[i-1,naterm+1] = (1+cos(theta[i]))/sin(theta[i])
    LHS[i-1,naterm+2] = (1-cos(theta[i]))/sin(theta[i])
    RHS[i-1] = dt1[i]
end

soln = LHS \ RHS
bterm2 = soln[1:end-2]
b02 = soln[end-1:end]
u2t = zeros(ndiv)
w2t = zeros(ndiv)
for i = 2:ndiv-1
    w2t[i] = b02[1]*(1+cos(theta[i]))/sin(theta[i])
    w2t[i] += b02[2]*(1-cos(theta[i]))/sin(theta[i])
    u2t[i] = b02[1] - b02[2]
    for n = 1:naterm
        u2t[i] -= bterm2[n]*cos(n*theta[i])
        w2t[i] += bterm2[n]*sin(n*theta[i])
    end
end
w2t[1] = 0.
w2t[ndiv] = 0.
u2t[1] = 2*u2t[2] - u2t[3]
u2t[ndiv] = 2*u2t[ndiv-1] - u2t[ndiv-2]


thsl_spl = Spline1D(x, thick_slope)
thick_dslope = derivative(thsl_spl, x)

q1 = 1+u1t
q2 = 1 + u1t + u2t + thick.*thick_dslope + 0.5*thick_slope.^2
q1c = sqrt(x./(x + 0.5*rho)).*q1
q2c = sqrt(x./(x + 0.5*rho)).*(q2 + rho./(4*x))
