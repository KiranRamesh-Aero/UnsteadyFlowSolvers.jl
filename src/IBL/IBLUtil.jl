#The description according to Numerical Recipies in Fortrant 90 :
# Given arrays x and y of length N containing a tabulated function, i.e., yi = f(xi), with x1 <
# x2 < ... < xN , and given values yp1 and ypn for the first derivative of the interpolating
# function at points 1 and N, respectively, this routine returns an array y2 of length N
# that contains the second derivatives of the interpolating function at the tabulated points
# xi. If yp1 and/or ypn are equal to 1 × 1030 or larger, the routine is signaled to set the
# corresponding boundary condition for a natural spline, with zero second derivative on that
# boundary

function spline(x::Vector{Float64}, y::Vector{Float64}, n::Int64, ypn::Float64, y2::Vector{Float64})

nmax::Int64 =3000
p,qn,sig,un :: Float64
u :: Vector{Float64,nmax}
tol =99e30

if (yp1>tol)
y[1]=0
u[1]=0
else
  y2[1]=-0.5
  u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1)
end

for i=2:n-1
sig= (x[i]-x[i-1])/(x[i+1]-x[i-1])
p=sig*y2[i-1]+2
y2[i]=(sig-1.)/p
u[i]=(6.0*((y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]))/(x[i+1]-x[i-1])-sig*u[i-1])/p
end

if(ypn>tol)
  qn=0.
  un=0.
else
  qn=0.5
  un = (3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]))
end
y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0)

for j in n-1:-1:1
 y2[k]=y2[k]*y2[k+1]+u[k]
end

return y2

end

#The following description according to Numerical Recipies in Fortrant 90 :
#Given the arrays xa and ya, which tabulate a function (with the xai ’s in increasing or
#decreasing order), and given the array y2a, which is the output from spline above, and
#given a value of x, this routine returns a cubic-spline interpolated value. The arrays xa, ya
# and y2a are all of the same size.


function splint(xa::Vector{Float64},ya::Vector{Float64},y2a::Vector{Float64},n::Int64,x::Float64,y::Float64)

khi::Int64=n
klo::Int64=1
k::Int64
a,b,h::Float64

klo=max(min(locate(xa,x),n-1),1)

khi=klo+1
h= xa[khi]-xa[klo]

if(h===0.0)
  error("h must be distinct")
end
  a=(xa[khi]-x)/h
  b=(x-xa[klo])/h
  y=a*ya[klo]+b*ya[khi]+((a^3-a)*y2a[klo]+(b^3-b)*y2a[khi])*(h^2)/6.0_sp

  return y

end

# The description according to the Numerical Recipies in Fortrant 90
# Given an array xx(1:N), and given a value x, returns a value j such that x is between
# xx(j) and xx(j + 1). xx must be monotonic, either increasing or decreasing. j = 0 or
# j = N is returned to indicate that x is out of range.


function locate(xa::Vector{Float64}, x::Float64)

locate::Int64
ascend::Bool
n,jl,jm,ju::Int64

n=length(xa)[1]

ascend(xa[n]>=xa[1])

jl=0
ju=u+1

while ((ju-jl)<=1)

jm=(jl+ju)/2   #replace the integer devision with dev(jl+ju,2)

if(ascend===(x>=xa[jm]))
jl=jm
else
  ju=jm
end
end

if (x === xx[1]) #then Then set the output, being careful with the endpoints.
locate=1
elseif (x === xx[n]) then
locate=n-1
else
locate=jl
end

return locate
end
