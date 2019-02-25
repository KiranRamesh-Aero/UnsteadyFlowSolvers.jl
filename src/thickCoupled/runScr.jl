using Printf
using Plots

n =200
del, E, F ,B = init(n)


x = collect(0:n+1)*pi/(n+1)
U0 = 2.0*sin.(x)
x =x[2:n+1]
U0 = U0[2:n+1]

Ut = zeros(n)
Ux= 2.0*cos.(x)

w1 = del
w2 = del.*(E.+1.0)
w0 = hcat(w1,w2)

function runSim(w0,U0,Ut,Ux,tEnd)

tStart = 0;
t1 =tStart;
deltaT = 0.005;
#tEnd = 1.5;
w = zeros(n,2)
sols = zeros(n,2)
dt =0.001

while t1 < tEnd

w, dt, j1 ,j2 = FVMIBL(w0,U0,Ut,Ux);
ww = w
del = ww[:,1]
E = (ww[:,2]./ww[:,1]) .- 1.0
ncell = length(del)
Csep = zeros(ncell-2)

Csep = abs.(diff(del[2:end]))./abs.(diff(del[1:end-1]));
#xtrunc = x[2:end-1]
sep = abs.(diff(Csep))
for i =1:length(j2)-1
   if j1[i] > 0.0002
        println("singularity (separation) detected at t=$t1, x=$(x[i]./pi), j1=$(j1[i])")
    end
end


# the plots of del and E

#display(plot(x/pi,[del, E], xticks = 0:0.1:1, layout=(2,1), legend = false))

# convergence of the Eigen-value
p1 = plot(x[1:end-1]/pi,j1)
p2 = plot(x[1:end-1]/pi,j2)

p3 = plot(x/pi,del)
p4 = plot(x/pi,E)
display(plot(p1, p2, p3, p4, xticks = 0:0.1:1, layout=(2,2), legend = false))

#display(plot(sep, xticks = 0:10:200, legend = false))
sleep(0.05)
w0 = w;

@printf("Time :%1.10f , Time step size %1.10f \n", t1, dt);
#map(x -> @sprintf("Seperation occure at :%1.10f  \n",x), xtrunc[seperation]./pi);

t1 = t1+dt

sols = w
end

return w
end
