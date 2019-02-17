using Printf
using Plots

n =100
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

while t1 < tEnd

w, dt = FVMIBL(w0,U0,Ut,Ux);
ww = w
del = ww[:,1]
E = (ww[:,2]./ww[:,1]) .- 1.0
ncell = length(del)
Csep = zeros(ncell-1)
#xtrunc = x[2:end-1]
#for i = 2:ncell-1
    #Csep[i] = (del[i] - del[i-1])/(del[i+1] - del[i])
    #if Csep[i] > 5.
    #    println("singularity (separation) detected at t=$t1, x=$(x[i]./pi), Csep=$(Csep[i])")
    #end
#end
display(plot(x/pi,[del, E], xticks = 0:0.1:1, layout=(2,1), legend = false))
sleep(0.05)
w0 = w;

@printf("Time :%1.10f , Time step size %1.10f \n", t1, dt);
#map(x -> @sprintf("Seperation occure at :%1.10f  \n",x), xtrunc[seperation]./pi);

t1 = t1+dt
end

return del, E
end
