#workspace()
#include("../src/UNSflow.jl")
#using UNSflow


alphadef = EldRampReturnDef(25,0.11,11)
hdef = ConstDef(0.)
udef = ConstDef(1.)
full_kinem = KinemDef(alphadef, hdef, udef)

pvt = 0.0 #leading edge

surf = TwoDSurf(1., 1., "sd7003_fine.dat", pvt, 70, 35, "Prescribed", full_kinem)

curfield = TwoDFlowField()

#Form panels

#Strategy - Should have greater density of panels near leading edge -


#Create splines of x and y coordinates as a function of s

airfoil = surf.coord_file
in_air = readdlm(airfoil);
x = in_air[:,1];
y = in_air[:,2];
n = length(x);

s = zeros(n)

for i = 2:n
    s[i] = s[i-1] + sqrt((x[i] - x[i-1])^2 + (y[i] - y[i-1])^2)
end

x_spl = Spline1D(s,x);
y_spl = Spline1D(s,y);


type TwoDPanel
    x1 :: Float64
    z1 :: Float64
    x2 :: Float64
    z2 :: Float64
    s :: Float64
end

sp = TwoDPanel[]

n_p = 200;

for i = 1:n_p
    s_p = (i-1)*s[n]/n_p
    s_n = i*s[n]/n_p
    x1 = evaluate(x_spl,s_p)
    x2 = evaluate(x_spl,s_n)
    z1 = evaluate(y_spl,s_p)
    z2 = evaluate(y_spl,s_n)
    push!(sp,TwoDPanel(x1,z1,x2,z2,0.))
end

for i = 2:n_p
    plot([x[i-1] x[i]],[y[i-1] y[i]])
end
