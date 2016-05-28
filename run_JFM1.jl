workspace()
using Dierckx
include("UNSflow.jl")
using UNSflow
#using PyPlot

c = 1;
u_ref = 1;
pvt = 0.0;
cm_pvt = 0.25;
foil_name = "sd7003_fine.dat";
re_ref = 10000;
lesp_crit = 0.18;
#motion_filename = "motion_25pr.dat";
force_filename = "force_test_25pr.dat";
flow_filename = "flow_test_25pr.dat";
np_flow = 1;

#Program variables
eps = 10e-6     #Tolerance or iteration
iter_max = 100  #Max. iterations
v_core = 0.02   #Non-dimensional core radius of point vortices
v_core = v_core*c;
ndiv = 70      #No. of divisions along chord on airfoil
naterm = 45    #Number of fourier terms used to compute vorticity at a location on chord
del_dist = 10
del_dist = del_dist*c;

#Generate motion
t = [0:0.014123742:10;];
t = t*c/u_ref;
nstep = length(t);

alpha = Array(Float64,nstep)
h = Array(Float64,nstep)
alphadot = Array(Float64,nstep)
hdot = Array(Float64,nstep)
u = Array(Float64,nstep)

x = Array(Float64,ndiv)
theta = Array(Float64,ndiv)
cam = Array(Float64,ndiv)
cam_slope = Array(Float64,ndiv)
bnd_x = Array(Float64,ndiv)
bnd_z = Array(Float64,ndiv)
bv_s = Array(Float64,ndiv)
bv_x = Array(Float64,ndiv)
bv_z = Array(Float64,ndiv)

cl = Array(Float64,nstep)
cd = Array(Float64,nstep)
cm = Array(Float64,nstep)
cl[1] = 0
cd[1] = 0
cm[1] = 0


alpha = eld_fn(t);
[h[i] = 0. for i = 1:nstep]
[u[i] = 1. for i = 1:nstep]

for i = 1:nstep
    h[i] = h[i]*c;
    u[i]=u[i]*u_ref;
    alpha[i] = alpha[i]*pi/180
end


alphadot[2:nstep] = diff(alpha)./diff(t)
hdot[2:nstep] = diff(h)./diff(t)
alphadot[1] = 2*alphadot[2] - alphadot[3]
hdot[1] = 2*hdot[2] - hdot[3]

#Defining chordwise divisions
dtheta = pi/(ndiv - 1);
for i = 1:ndiv
    theta[i] = (real(i) - 1)*dtheta;
    x[i] = 0.5c*(1 - cos(theta[i]));
end

cam,cam_slope = camber_calc(x,foil_name);

#Initial conditions - Figure out how to include all this
n_lev = 0
n_tev=0

aterm = Array(Float64,naterm)
aterm_prev = Array(Float64,naterm)

#levflag=0
#dist_wind=0
kelv_enf=0
# uind_lev_prev(:)=0
# uind_tev_prev(:)=0
# wind_lev_prev(:)=0
# wind_tev_prev(:)=0
# uind_lev_prprev(:)=0
# uind_tev_prprev(:)=0
# wind_lev_prprev(:)=0
# wind_tev_prprev(:)=0

#Place initial bound vortices
bnd_x, bnd_z = start_bound(alpha[1],h[1],pvt,ndiv,c,x,cam)

 tev_x = Float64[];
    tev_z = Float64[];
    lev_x = Float64[];
    lev_z = Float64[];
    tev_s = Float64[];
    lev_s = Float64[];
    lev_vc = Float64[];
    tev_vc = Float64[];
    bv_vc = Float64[];
    uind_tev = Float64[];
    wind_tev = Float64[];
    uind_lev = Float64[];
    wind_lev = Float64[];

uind = Array(Float64,ndiv)
    wind = Array(Float64,ndiv)
    dwash = Array(Float64,ndiv)

adot = Array(Float64,3)
a0 = 0
a0_prev = 0
    aterm_prev[1:3] = 0
    gam = 0
  kelv_enf = 0
   tev_iter = 0
#Iterate over time steps
for ist = 2:nstep

    #Update bound vortex positions
    bnd_x, bnd_z = update_bound(bnd_x,bnd_z,alphadot[ist],alpha[ist],h[ist],hdot[ist],pvt,ndiv,c,x,cam,u[ist],t[ist]-t[ist-1])

    #Shed a TEV at every time step

    tev_x, tev_z, tev_vc = add_tev!(tev_x, tev_z, tev_vc, bnd_x, bnd_z, ndiv, u[ist], t[ist]-t[ist-1], c)

    push!(tev_s,0.01) #Initial value for iteration

    #push!(tev_s,fzero(tev_iter -> kelv(),0.))

    #For downwash

#    dwash = downwash(uind,wind,lev_s,lev_x,lev_z,lev_vc,tev_s,tev_x,tev_z,tev_vc,bnd_x,bnd_z,u[ist],alpha[ist],alphadot[ist],hdot[ist],pvt,c,cam_slope,ndiv,x)

 #   gam, a0 = bcirc(a0,aterm,dwash,theta,ndiv,dtheta,u_ref,c)


#    kelv = kelv1d(a0,aterm,dwash,tev_iter, tev_s, lev_s,kelv_enf,theta,dtheta,ndiv,u_ref,c)

    tev_s[length(tev_s)] = fzero(tev_iter -> kelv1d(tev_iter,kelv_enf,theta,dtheta,ndiv,u_ref,c,lev_s,lev_x,lev_z,lev_vc,tev_s,tev_x,tev_z,tev_vc,bnd_x,bnd_z,u[ist],alpha[ist],alphadot[ist],hdot[ist],pvt,cam_slope,x),[10, -10])

    uind, wind, dwash = downwash(lev_s,lev_x,lev_z,lev_vc,tev_s,tev_x,tev_z,tev_vc,bnd_x,bnd_z,u[ist],alpha[ist],alphadot[ist],hdot[ist],pvt,c,cam_slope,ndiv,x)

    a0, aterm[1], gam = bcirc(theta,ndiv,dtheta,u_ref,c,dwash)

    a0dot, adot = calc_adot(a0,a0_prev,aterm,aterm_prev,dwash,theta,dtheta,t[ist]-t[ist-1],u_ref,ndiv)

    aterm = calc_aterm(aterm,naterm,theta,dtheta,ndiv,u_ref,lev_s,lev_x,lev_z,lev_vc,tev_s,tev_x,tev_z,tev_vc,bnd_x,bnd_z,u[ist],alpha[ist],alphadot[ist],hdot[ist],pvt,c,cam_slope,x,dwash)

    a0_prev = a0
    aterm_prev[1:3]=aterm[1:3]

    bv_s, bv_x, bv_z, bv_vc = calc_bndvorts(a0,aterm,naterm,theta,dtheta,ndiv,bnd_x,bnd_z,u_ref,c)

    #!Calculate bound_vortex strengths

    # !Wake rollup
    uind_tev = Array(Float64,length(tev_s))
    uind_lev = Array(Float64,length(lev_s))
    wind_tev = Array(Float64,length(tev_s))
    wind_lev = Array(Float64,length(lev_s))

    uind_tev[1:length(tev_s)] = 0
    wind_tev[1:length(tev_s)] = 0
    uind_lev[1:length(lev_s)] = 0
    wind_lev[1:length(lev_s)] = 0

    uind_tev, wind_tev = ind_vel(uind_tev,wind_tev,[lev_s; tev_s; bv_s[2:ndiv]],[lev_x; tev_x; bv_x[2:ndiv]],[lev_z; tev_z; bv_z[2:ndiv]],[lev_vc; tev_vc; bv_vc[2:ndiv]],tev_x,tev_z)

    uind_lev, wind_lev = ind_vel(uind_lev,wind_lev,[lev_s; tev_s; bv_s[2:ndiv]],[lev_x; tev_x; bv_x[2:ndiv]],[lev_z; tev_z; bv_z[2:ndiv]],[lev_vc; tev_vc; bv_vc[2:ndiv]],lev_x,lev_z)

   tev_x, tev_z = convect_vort(tev_x,tev_z,uind_tev,wind_tev,t[ist]-t[ist-1])
   lev_x, lev_z = convect_vort(lev_x,lev_z,uind_lev,wind_lev,t[ist]-t[ist-1])

    cl[ist], cd[ist], cm[ist] = force_calc(u[ist],alpha[ist],hdot[ist],u_ref,c,a0,aterm,a0dot, adot,uind,wind,ndiv,bv_s,x,cm_pvt)


end

