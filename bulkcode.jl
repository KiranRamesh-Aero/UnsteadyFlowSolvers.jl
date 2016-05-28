function cam_calc(x,airfoil;c=1)
#Determine camberslope on airfoil from airfoil input file
using Dierckx

n_div = length(x);
in_air = readdlm(airfoil);
xcoord = in_air[:,1];
ycoord = in_air[:,2];
ncoord = length(xcoord);
xcoord_sum[1] = 0	  			
xcoord_sum[2:ncoord] = xcoord_sum[1:ncoord-1]+abs(xcoord[2:ncoord]-xcoord[1:ncoord-1]);
y_spl = Spline1D(xcoord_sum,ycoord)
y_ans[1:n_div] = evaluate(y_spl,x[1:n_div]/c);
y_ans[n_div+1:2*n_div] = evaluate(y_spl,(x[n_div]/c)+(x[(n_div+1:2*n_div)-n_div]/c);
cam_calc[1:n_div] = (y_ans[n_div:-1:1]+y_ans[(2*n_div)+1-(n_div:-1:1)])*c/2;
cam_calc[1] = 0.0
end    