#=!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!UNSflow 
!Created by Kiran Ramesh in 2016
!University of Glasgow
!Code to solve unsteady flow problems
!Papers on LESP and LAUTAT theories at :: http://www.mae.ncsu.edu/apa/kiran_ramesh/Publications.html

!Copyright 2016 Kiran Ramesh 

!     This program is free software: you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation, either version 3 of
!     the License, or (at your option) any later version.

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public
!     License along with this program.  If not, see
!     <http://www.gnu.org/licenses/>.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! =#



function camber_calc(x::Vector,airfoil;c=1.)
#Determine camberslope on airfoil from airfoil input file

n_div = length(x);
in_air = readdlm(airfoil);
xcoord = in_air[:,1];
ycoord = in_air[:,2];
ncoord = length(xcoord);
xcoord_sum=Float64[0 for i in 1:ncoord];
for icoord=1:ncoord-1
    xcoord_sum[icoord+1] = xcoord_sum[icoord]+abs(xcoord[icoord+1]-xcoord[icoord]);
end
y_spl = Spline1D(xcoord_sum,ycoord);
y_ans=Float64[0. for i in 1:2*n_div];
#cam=Float64[0. for i in 1:n_div];
y_ans[1:n_div] = evaluate(y_spl,x[1:n_div]/c);
y_ans[n_div+1:2*n_div] = evaluate(y_spl,(x[n_div]/c)+(x[(n_div+1:2*n_div)-n_div]/c));
cam[1:n_div] =[(y_ans[i_div]+y_ans[(2*n_div)+1-i_div])*c/2 for i_div=n_div:-1:1];
cam[1] = 0;
cam_spl = Spline1D(x,cam);
cam_slope[1:n_div]=derivative(cam_spl,x);
return cam, cam_slope
end    


function eld_fn(t::AbstractVector;K=0.2,amp=45,t1=1.,tf=1.,a=11.)
fr = K/(pi*abs(amp)*pi/180);
t2=t1+(1./(2*pi*fr));
t3 = t2+((1/(4*fr))-(1/(2*pi*fr)));
t4 = t3+(1./(2*pi*fr));
t5 = t4+tf;

g = log((cosh(a*(t-t1)).*cosh(a*(t-t4)))./(cosh(a*(t-t2)).*cosh(a*(t-t3))));
maxg = maximum(g);
return res = amp*g/maxg;
end		

function const_fn(t::Vector;amp=45)
res = amp;
end



#Assume input is given here
using Dierckx
c = 1;
u_ref = 1;
pvt = 0.25;
cm_pvt = 0.25;
foil_name = "sd7003_fine.dat";
re_ref = 10000;
lesp_crit = 0.18;
#motion_filename = "motion_25pr.dat";
force_filename = "force_test_25pr.dat";
flow_filename = "flow_test_25pr.dat";
np_flow = 1;

#Program variables
eps=10e-6     #Tolerance or iteration
iter_max=100  #Max. iterations
v_core=0.02   #Non-dimensional core radius of point vortices 
n_div=70      #No. of divisions along chord on airfoil
n_aterm=45    #Number of fourier terms used to compute vorticity at a location on chord  
del_dist=10

#Generate motion
t = [0:0.015:10;];
alpha=eld_fn(t);
h = const_fn(t,amp=0);
u = const_fn(t,amp=1);

n_step = length(t);

#Convert motion inputs to dimensional quantities
t=t*c/u_ref;
alpha=alpha*pi/180;
vh=h*c;
u=u*u_ref;
v_core=v_core*c;
del_dist=del_dist*c;

#Defining chordwise divisions
dtheta=pi/(n_div-1);
theta=(collect(1:n_div)-1)*dtheta;
x=0.5c*(1-cos(theta));
cam=Float64[0. for i in 1:n_div]

cam_slope=Float64[0. for i in 1:n_div]
#cam[1:n_div]],[cam_slope[1:n_div] = 
camber_calc(x,foil_name);






