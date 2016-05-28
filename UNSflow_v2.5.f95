!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!UNSflow
!Created by Kiran Ramesh in 2013
!Applied Aerodynamics Group, NC State University
!Code to solve 2D unsteady flow problems : Using unsteady thin-airfoil theory augmented with LEV model
!Papers on LESP and LAUTAT theories at :: http://www.mae.ncsu.edu/apa/kiran_ramesh/Publications.html
!Major change - LESP_crit is calibrated in accordance with vel. magnitude at leading edge
!This version is optimized for speed - some base quantities are not calculated/printed
!These can be obtained using UNSflow_v2.0
!LEVs and TEVs can be deleted from the solution when they are a set distance from airfoil

!Copyright 2013 Kiran Ramesh

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program UNSflow

implicit none

double precision :: pi, eps
double precision :: c, u_ref, pvt, v_core, cm_pvt, del_dist
integer, parameter :: expect_vort=3000
integer :: n_div, n_step, n_aterm, n_coord, n_lev, n_tev, n_calib
integer :: i_div, i_step, i_coord, io1, iter_max, iter, i_aterm, i_calib
integer :: i_lev, i_tev, j_lev, j_tev, levflag
character(len=50) :: foil_name, motion_filename, flow_filename, force_filename, input_filename, cp_filename
character(len=50) :: lesp_calib
double precision, dimension(100000) :: t, alpha, h, u, alphadot, hdot
double precision, dimension(300) :: x, theta, cam, cam_slope, gamma
double precision :: dtheta, xreq
double precision, dimension(1000) :: recalib, lespcalib, lesp_splined
double precision, dimension(0:50) :: aterm
double precision, dimension(0:3) :: aterm_prev, adot
double precision, dimension(300,3) :: bound, bound_int
double precision, dimension(0:100) :: tev_iter, lev_iter, kelv, kutta
double precision, dimension(expect_vort,3) :: tev, lev
double precision, dimension(3,3,expect_vort,expect_vort) :: xdist, zdist
double precision :: lesp, lesp_crit, le_vel_x, le_vel_y, re_le, re_ref, lesp_cond, vmag
double precision, dimension(expect_vort) :: uind_tev, wind_tev, uind_lev, wind_lev
double precision, dimension(expect_vort) :: uind_tev_prev, wind_tev_prev, uind_lev_prev, wind_lev_prev
double precision, dimension(expect_vort) :: uind_tev_prprev, wind_tev_prprev, uind_lev_prprev, wind_lev_prprev
double precision, dimension(2000) :: xcoord, ycoord, xcoord_sum, ysplined, ycoord_ans
double precision, dimension(300) :: uind, wind, downwash
double precision :: bound_circ, cnc, cnnc, cn, cs, cl, cd, cm, nonl, nonl_m
double precision :: dkelv, dkelv_tev, dkelv_lev, dkutta_lev, dkutta_tev
double precision :: kelv_lev, kelv_tev, kelv_enf, kutta_lev, kutta_tev
double precision :: bound_int_zdist, bound_int_xdist
double precision :: dist,  dist_wind
integer :: n_pts_flow, n_pts_cp
double precision :: t1, t2, dt
real :: version
!double precision :: lera=0.008

pi=acos(-1.)

!Program variables
eps=10e-6     !Tolerance or iteration
iter_max=100  !Max. iterations
v_core=0.02   !Non-dimensional core radius of point vortices
n_div=70     !No. of divisions along chord on airfoil
n_aterm=45    !Number of fourier terms used to compute vorticity at a location on chord
del_dist=10

version=2.5
write(*,1001)version
1001 format( &
& /'================================================================'&
& /'UNSflow Version',F5.2&
& /'Copyright (C) 2013 Kiran Ramesh'&
& //'This software comes with ABSOLUTELY NO WARRANTY,' &
& /'subject to the GNU General Public License.'&
& /'================================================================')

write(*,*)'Enter name of input file for analysis (Default:input_UNSflow.dat)'
read(*,103)input_filename
103 format(A30)
if(input_filename .ne. "") then
   open(unit=3, file=input_filename,status='unknown')
else
   open(unit=3, file='input_UNSflow.dat',status='unknown')
end if

call cpu_time(t1)

!Read input file
read(3,*)c
read(3,*)u_ref
read(3,*)pvt
read(3,*)cm_pvt
read(3,*)foil_name
read(3,*)re_ref
read(3,*)lesp_crit
read(3,*)motion_filename
read(3,*)force_filename
read(3,*)flow_filename, n_pts_flow
write(*,*)motion_filename
!Open output files
if (flow_filename .ne. 'nil') then
   open(unit=13, file=flow_filename, status='replace')
   write(13,*)'NaN ','NaN ','NaN '
end if
open(unit=17, file=force_filename, status='replace')


!Reading motion file
open(unit=7, file=motion_filename, status='unknown')
i_step=0
do
   read(7,*,iostat=io1)t(i_step+1),alpha(i_step+1),h(i_step+1),u(i_step+1)
   if (io1<0) then
      n_step=i_step
      exit
   end if
   i_step=i_step+1
end do

!Convert motion inputs to dimensional quantities
t(1:n_step)=t(1:n_step)*c/u_ref
alpha(1:n_step)=alpha(1:n_step)*pi/180
h(1:n_step)=h(1:n_step)*c
u(1:n_step)=u(1:n_step)*u_ref
v_core=v_core*c
del_dist=del_dist*c

!Defining chordwise divisions
dtheta=pi/(n_div-1)
do i_div=1,n_div
   theta(i_div)=(i_div-1)*dtheta
   x(i_div)=(c/2.)*(1-cos(theta(i_div)))
end do

!Determine camberslope on airfoil from airfoil input file
cam_slope(1:n_div)=0
cam(1:n_div)=0
call calc_camberslope

!Calculating derivatives (Backward differentiation)
do i_step=2,n_step
   alphadot(i_step)=(alpha(i_step)-alpha(i_step-1))/(t(i_step)-t(i_step-1))
   hdot(i_step)=(h(i_step)-h(i_step-1))/(t(i_step)-t(i_step-1))
end do
alphadot(1)=alphadot(2)
hdot(1)=hdot(2)

!Initial conditions
n_lev=0
n_tev=0
aterm(0:n_aterm)=0
aterm_prev(0:n_aterm)=0
levflag=0
dist_wind=0
kelv_enf=0
uind_lev_prev(:)=0
uind_tev_prev(:)=0
wind_lev_prev(:)=0
wind_tev_prev(:)=0
uind_lev_prprev(:)=0
uind_tev_prprev(:)=0
wind_lev_prprev(:)=0
wind_tev_prprev(:)=0

!Iteration
do i_step=2,n_step
   write(*,*)i_step, n_tev, n_lev

   !Calculate bound vortex positions at this time step
      dist_wind=dist_wind+(u(i_step-1)*(t(i_step)-t(i_step-1)))
      do i_div=1,n_div
         bound(i_div,2)=-((c-pvt*c)+((pvt*c-x(i_div))*cos(alpha(i_step)))+dist_wind)&
              &+(cam(i_div)*sin(alpha(i_step)))
         bound(i_div,3)=h(i_step)+((pvt*c-x(i_div))*sin(alpha(i_step)))+(cam(i_div)&
              &*cos(alpha(i_step)))
      end do

   !TEV shed at every time step
   n_tev=n_tev+1
   tev_iter(0)=0
   tev_iter(1)=-0.01
   if (n_tev==1) then
      tev(n_tev,2)=bound(n_div,2)+(0.5*u(i_step)*(t(i_step)-t(i_step-1)))
      tev(n_tev,3)=bound(n_div,3)
   else
      tev(n_tev,2)=bound(n_div,2)+((1./3.)*(tev(n_tev-1,2)-bound(n_div,2)))
      tev(n_tev,3)=bound(n_div,3)+((1./3.)*(tev(n_tev-1,3)-bound(n_div,3)))
   end if

   !Precalculate xdist and zdist arrays
   do i_div=1,n_div
      do i_tev=1,n_tev
         xdist(1,2,i_div,i_tev)=tev(i_tev,2)-bound(i_div,2)
         zdist(1,2,i_div,i_tev)=tev(i_tev,3)-bound(i_div,3)
      end do
   end do
   do i_div=1,n_div
      do i_lev=1,n_lev
         xdist(1,3,i_div,i_lev)=lev(i_lev,2)-bound(i_div,2)
         zdist(1,3,i_div,i_lev)=lev(i_lev,3)-bound(i_div,3)
      end do
   end do
   do i_tev=1,n_tev
      do j_tev=1,n_tev
         xdist(2,2,i_tev,j_tev)=tev(j_tev,2)-tev(i_tev,2)
         zdist(2,2,i_tev,j_tev)=tev(j_tev,3)-tev(i_tev,3)
      end do
   end do
   do i_tev=1,n_tev
      do i_lev=1,n_lev
         xdist(2,3,i_tev,i_lev)=lev(i_lev,2)-tev(i_tev,2)
         zdist(2,3,i_tev,i_lev)=lev(i_lev,3)-tev(i_tev,3)
      end do
   end do
   do i_lev=1,n_lev
      do j_lev=1,n_lev
         xdist(3,3,i_lev,j_lev)=lev(j_lev,2)-lev(i_lev,2)
         zdist(3,3,i_lev,j_lev)=lev(j_lev,3)-lev(i_lev,3)
      end do
   end do

   !Iterating to find A0 value assuming no LEV is formed
   iter=0
   do
      iter=iter+1
      if (iter>iter_max) then
         write(*,*)i_step,'1D iteration failed'
         stop
      end if
      tev(n_tev,1)=tev_iter(iter)
      call calc_downwash_boundcirc
      kelv(iter)=kelv_enf
      do i_lev=1,n_lev
         kelv(iter)=kelv(iter)+lev(i_lev,1)
      end do
      do i_tev=1,n_tev
         kelv(iter)=kelv(iter)+tev(i_tev,1)
      end do
      kelv(iter)=kelv(iter)+bound_circ
      if (abs(kelv(iter))<eps) then
         exit
      end if
      dkelv=(kelv(iter)-kelv(iter-1))/(tev_iter(iter)-tev_iter(iter-1))
      tev_iter(iter+1)=tev_iter(iter)-(kelv(iter)/dkelv)
   end do
   aterm(2:3)=0
   do i_aterm=2,3
      do i_div=2,n_div
         aterm(i_aterm)=aterm(i_aterm)+((((downwash(i_div)*cos(i_aterm*theta(i_div)))&
              &+(downwash(i_div-1)*cos(i_aterm*theta(i_div-1))))/2)*dtheta)
      end do
      aterm(i_aterm)=(2./(u_ref*pi))*aterm(i_aterm)
   end do
   adot(0)=(aterm(0)-aterm_prev(0))/(t(i_step)-t(i_step-1))
   adot(1)=(aterm(1)-aterm_prev(1))/(t(i_step)-t(i_step-1))
   adot(2)=(aterm(2)-aterm_prev(2))/(t(i_step)-t(i_step-1))
   adot(3)=(aterm(3)-aterm_prev(3))/(t(i_step)-t(i_step-1))

   le_vel_x=(u(i_step))-(alphadot(i_step)*sin(alpha(i_step))*pvt*c)+uind(1)
   le_vel_y=-(alphadot(i_step)*cos(alpha(i_step))*pvt*c)-(hdot(i_step))+wind(1)
   vmag=sqrt(le_vel_x*le_vel_x+le_vel_y*le_vel_y)
   re_le=re_ref*vmag/u_ref
   lesp=aterm(0)





   !2D iteration if LESP_crit is exceeded
   if (abs(lesp)>lesp_crit) then
      if (lesp>0) then
         lesp_cond=lesp_crit
      else
         lesp_cond=-lesp_crit
      end if
      n_lev=n_lev+1
      tev_iter(0)=0
      tev_iter(1)=-0.01
      lev_iter(0)=0
      lev_iter(1)=0.01
      if (levflag==0) then
         lev(n_lev,2)=bound(1,2)+(0.5*le_vel_x*(t(i_step)-t(i_step-1)))
         lev(n_lev,3)=bound(1,3)+(0.5*le_vel_y*(t(i_step)-t(i_step-1)))
      else
         lev(n_lev,2)=bound(1,2)+((1./3.)*(lev(n_lev-1,2)-bound(1,2)))
         lev(n_lev,3)=bound(1,3)+((1./3.)*(lev(n_lev-1,3)-bound(1,3)))
      end if
      levflag=1

      !Updating the xdist and zdist arrays
      i_lev=n_lev
      do i_div=1,n_div
         xdist(1,3,i_div,i_lev)=lev(i_lev,2)-bound(i_div,2)
         zdist(1,3,i_div,i_lev)=lev(i_lev,3)-bound(i_div,3)
      end do
      i_lev=n_lev
      do i_tev=1,n_tev
         xdist(2,3,i_tev,i_lev)=lev(i_lev,2)-tev(i_tev,2)
         zdist(2,3,i_tev,i_lev)=lev(i_lev,3)-tev(i_tev,3)
      end do
      j_lev=n_lev
      do i_lev=1,n_lev
         xdist(3,3,i_lev,j_lev)=lev(j_lev,2)-lev(i_lev,2)
         zdist(3,3,i_lev,j_lev)=lev(j_lev,3)-lev(i_lev,3)
      end do
      i_lev=n_lev
      do j_lev=1,n_lev
         xdist(3,3,i_lev,j_lev)=lev(j_lev,2)-lev(i_lev,2)
         zdist(3,3,i_lev,j_lev)=lev(j_lev,3)-lev(i_lev,3)
      end do

      iter=0
      do
         iter=iter+1
         if (iter>iter_max) then
            write(*,*)i_step,'2D NR iteration failed'
            stop
         end if
         !Advancing with tev strength
         lev(n_lev,1)=lev_iter(iter-1)
         tev(n_tev,1)=tev_iter(iter)
         call calc_downwash_boundcirc
         kelv_tev=kelv_enf
         do i_lev=1,n_lev
            kelv_tev=kelv_tev+lev(i_lev,1)
         end do
         do i_tev=1,n_tev
            kelv_tev=kelv_tev+tev(i_tev,1)
         end do
         kelv_tev=kelv_tev+bound_circ
         kutta_tev=aterm(0)-lesp_cond
         dkelv_tev=(kelv_tev-kelv(iter-1))/(tev_iter(iter)-tev_iter(iter-1))
         dkutta_tev=(kutta_tev-kutta(iter-1))/(tev_iter(iter)-tev_iter(iter-1))
         !Advancing with lev strength
         lev(n_lev,1)=lev_iter(iter)
         tev(n_tev,1)=tev_iter(iter-1)
         call calc_downwash_boundcirc
         kelv_lev=kelv_enf
         do i_lev=1,n_lev
            kelv_lev=kelv_lev+lev(i_lev,1)
         end do
         do i_tev=1,n_tev
            kelv_lev=kelv_lev+tev(i_tev,1)
         end do
         kelv_lev=kelv_lev+bound_circ
         kutta_lev=aterm(0)-lesp_cond
         dkelv_lev=(kelv_lev-kelv(iter-1))/(lev_iter(iter)-lev_iter(iter-1))
         dkutta_lev=(kutta_lev-kutta(iter-1))/(lev_iter(iter)-lev_iter(iter-1))
         !Advancing with both
         lev(n_lev,1)=lev_iter(iter)
         tev(n_tev,1)=tev_iter(iter)
         call calc_downwash_boundcirc
         kelv(iter)=kelv_enf
         do i_lev=1,n_lev
            kelv(iter)=kelv(iter)+lev(i_lev,1)
         end do
         do i_tev=1,n_tev
            kelv(iter)=kelv(iter)+tev(i_tev,1)
         end do
         kelv(iter)=kelv(iter)+bound_circ
         kutta(iter)=aterm(0)-lesp_cond
         if (abs(kelv(iter))<eps .and. abs(kutta(iter))<eps) then
            exit
         end if
         tev_iter(iter+1)=tev_iter(iter)-((1/(dkelv_tev*dkutta_lev-dkelv_lev*dkutta_tev))*&
              &((dkutta_lev*kelv(iter))-(dkelv_lev*kutta(iter))))
         lev_iter(iter+1)=lev_iter(iter)-((1/(dkelv_tev*dkutta_lev-dkelv_lev*dkutta_tev))*&
              &((-dkutta_tev*kelv(iter))+(dkelv_tev*kutta(iter))))
      end do
   else
      levflag=0
   end if

   !To remove any massive starting vortices
   !if (i_step==2) then
   !   tev(1,1)=0
   !end if
   !Calculate fourier terms and bound vorticity

   aterm(2:n_aterm)=0
   do i_aterm=2,n_aterm
      do i_div=2,n_div
         aterm(i_aterm)=aterm(i_aterm)+((((downwash(i_div)*cos(i_aterm*theta(i_div)))&
              &+(downwash(i_div-1)*cos(i_aterm*theta(i_div-1))))/2)*dtheta)
      end do
      aterm(i_aterm)=(2./(u_ref*pi))*aterm(i_aterm)
   end do

   !Set previous values of aterm to be used for derivatives in next time step
   aterm_prev(0:3)=aterm(0:3)

   !Calculate bound_vortex strengths
   do i_div=1,n_div
      gamma(i_div)=(aterm(0)*(1+cos(theta(i_div))))
      do i_aterm=1,n_aterm
         gamma(i_div)=gamma(i_div)+(aterm(i_aterm)*sin(i_aterm*theta(i_div))*sin(theta(i_div)))
      end do
      gamma(i_div)=gamma(i_div)*u_ref*c
   end do

   do i_div=2,n_div
      bound_int(i_div,1)=((gamma(i_div)+gamma(i_div-1))/2)*dtheta
      bound_int(i_div,2)=(bound(i_div,2)+bound(i_div-1,2))/2
      bound_int(i_div,3)=(bound(i_div,3)+bound(i_div-1,3))/2
   end do

   !Wake rollup
   uind_tev(1:n_tev)=0
   wind_tev(1:n_tev)=0
   do i_tev=1,n_tev
      do j_tev=1,n_tev
         if (i_tev .ne. j_tev) then
            dist=xdist(2,2,i_tev,j_tev)**2+zdist(2,2,i_tev,j_tev)**2
            uind_tev(i_tev)=uind_tev(i_tev)+&
                 &((tev(j_tev,1)*(-zdist(2,2,i_tev,j_tev)))/(2*pi*sqrt(v_core**4+dist**2)))
            wind_tev(i_tev)=wind_tev(i_tev)+&
                 &((-tev(j_tev,1)*(-xdist(2,2,i_tev,j_tev)))/(2*pi*sqrt(v_core**4+dist**2)))
         end if
      end do
      do i_lev=1,n_lev
         dist=xdist(2,3,i_tev,i_lev)**2+zdist(2,3,i_tev,i_lev)**2
         uind_tev(i_tev)=uind_tev(i_tev)+&
              &((lev(i_lev,1)*(-zdist(2,3,i_tev,i_lev)))/(2*pi*sqrt(v_core**4+dist**2)))
         wind_tev(i_tev)=wind_tev(i_tev)+&
              &((-lev(i_lev,1)*(-xdist(2,3,i_tev,i_lev)))/(2*pi*sqrt(v_core**4+dist**2)))
      end do
      do i_div=2,n_div
         bound_int_xdist=tev(i_tev,2)-bound_int(i_div,2)
         bound_int_zdist=tev(i_tev,3)-bound_int(i_div,3)
         dist=bound_int_xdist**2+bound_int_zdist**2
         uind_tev(i_tev)=uind_tev(i_tev)+((bound_int(i_div,1)*bound_int_zdist)/(2*pi*sqrt(v_core**4+dist**2)))
         wind_tev(i_tev)=wind_tev(i_tev)+((-bound_int(i_div,1)*bound_int_xdist)/(2*pi*sqrt(v_core**4+dist**2)))
      end do
   end do
   uind_lev(1:n_lev)=0
   wind_lev(1:n_lev)=0
   do i_lev=1,n_lev
      do j_lev=1,n_lev
         if (i_lev .ne. j_lev) then
            dist=xdist(3,3,i_lev,j_lev)**2+zdist(3,3,i_lev,j_lev)**2
            uind_lev(i_lev)=uind_lev(i_lev)+((lev(j_lev,1)*(-zdist(3,3,i_lev,j_lev)))/(2*pi*sqrt(v_core**4+dist**2)))
            wind_lev(i_lev)=wind_lev(i_lev)+((-lev(j_lev,1)*(-xdist(3,3,i_lev,j_lev)))/(2*pi*sqrt(v_core**4+dist**2)))
         end if
      end do
      do i_tev=1,n_tev
         dist=xdist(2,3,i_tev,i_lev)**2+zdist(2,3,i_tev,i_lev)**2
         uind_lev(i_lev)=uind_lev(i_lev)+((tev(i_tev,1)*zdist(2,3,i_tev,i_lev))/(2*pi*sqrt(v_core**4+dist**2)))
         wind_lev(i_lev)=wind_lev(i_lev)+((-tev(i_tev,1)*xdist(2,3,i_tev,i_lev))/(2*pi*sqrt(v_core**4+dist**2)))
      end do
      do i_div=2,n_div
         bound_int_xdist=lev(i_lev,2)-bound_int(i_div,2)
         bound_int_zdist=lev(i_lev,3)-bound_int(i_div,3)
         dist=bound_int_xdist**2+bound_int_zdist**2
         uind_lev(i_lev)=uind_lev(i_lev)+((bound_int(i_div,1)*bound_int_zdist)/(2*pi*sqrt(v_core**4+dist**2)))
         wind_lev(i_lev)=wind_lev(i_lev)+((-bound_int(i_div,1)*bound_int_xdist)/(2*pi*sqrt(v_core**4+dist**2)))
      end do
   end do

   !Adam-bashforth method
   dt=t(i_step)-t(i_step-1)
   ! do i_tev=1,n_tev
   !    tev(i_tev,2)=tev(i_tev,2)+((dt/12)*((23*uind_tev(i_tev))-&
   !         &(16*uind_tev_prev(i_tev))+(5*uind_tev_prprev(i_tev))))
   !    tev(i_tev,3)=tev(i_tev,3)+((dt/12)*((23*wind_tev(i_tev))-&
   !         &(16*wind_tev_prev(i_tev))+(5*wind_tev_prprev(i_tev))))
   ! end do
   ! do i_lev=1,n_lev
   !    lev(i_lev,2)=lev(i_lev,2)+((dt/12)*((23*uind_lev(i_lev))-&
   !         &(16*uind_lev_prev(i_lev))+(5*uind_lev_prprev(i_lev))))
   !    lev(i_lev,3)=lev(i_lev,3)+((dt/12)*((23*wind_lev(i_lev))-&
   !         &(16*wind_lev_prev(i_lev))+(5*wind_lev_prprev(i_lev))))
   ! end do

   do i_tev=1,n_tev
      tev(i_tev,2)=tev(i_tev,2)+(dt*uind_tev(i_tev))
      tev(i_tev,3)=tev(i_tev,3)+(dt*wind_tev(i_tev))
   end do
   do i_lev=1,n_lev
      lev(i_lev,2)=lev(i_lev,2)+(dt*uind_lev(i_lev))
      lev(i_lev,3)=lev(i_lev,3)+(dt*wind_lev(i_lev))
   end do





   uind_lev_prprev(1:n_lev)=uind_lev_prev(1:n_lev)
   uind_lev_prev(1:n_lev)=uind_lev(1:n_lev)
   wind_lev_prprev(1:n_lev)=wind_lev_prev(1:n_lev)
   wind_lev_prev(1:n_lev)=wind_lev(1:n_lev)
   uind_tev_prprev(1:n_tev)=uind_tev_prev(1:n_tev)
   uind_tev_prev(1:n_tev)=uind_tev(1:n_tev)
   wind_tev_prprev(1:n_tev)=wind_tev_prev(1:n_tev)
   wind_tev_prev(1:n_tev)=wind_tev(1:n_tev)

   !Remove TEVs and LEVs thats have crossed a certain distance and update kelvin condition
   if (tev(1,2)-bound(n_div,2)>del_dist) then
      do i_tev=1,n_tev-1
         tev(i_tev,:)=tev(i_tev+1,:)
      end do
      n_tev=n_tev-1
      kelv_enf=kelv_enf+tev(1,1)
   end if
   if (n_lev>0 .and. lev(1,2)-bound(n_div,2)>del_dist) then
      do i_lev=1,n_lev-1
         lev(i_lev,:)=lev(i_lev+1,:)
      end do
      n_lev=n_lev-1
      kelv_enf=kelv_enf+lev(1,1)
   end if



   !Load coefficient calculation (nondimensional units)
   cnc=(2*pi*((u(i_step)*cos(alpha(i_step))/u_ref)+(hdot(i_step)*sin(alpha(i_step))/u_ref))*&
        &(aterm(0)+(aterm(1)/2)))
   cnnc=(2*pi*((3*c*adot(0)/(4*u_ref))+(c*adot(1)/(4*u_ref))+&
        &(c*adot(2)/(8*u_ref))))
   cs=2*pi*aterm(0)*aterm(0)
   !The components of normal force and moment from induced velocities are
   !calulcated in dimensional units and nondimensionalized later
   nonl=0
   nonl_m=0
   do i_div=2,n_div
      nonl=nonl+(((uind(i_div)*cos(alpha(i_step)))-(wind(i_div)&
           &*sin(alpha(i_step))))*bound_int(i_div,1))
      nonl_m=nonl_m+(((uind(i_div)*cos(alpha(i_step)))-(wind(i_div)&
           &*sin(alpha(i_step))))*(x(i_div))*bound_int(i_div,1))
   end do
   nonl=nonl*(2/(u_ref*u_ref*c))
   nonl_m=nonl_m*(2/(u_ref*u_ref*c*c))

   cn=cnc+cnnc+nonl
   cl=cn*cos(alpha(i_step))+cs*sin(alpha(i_step))
   cd=cn*sin(alpha(i_step))-cs*cos(alpha(i_step))
   !Pitching moment is clockwise or nose up positive
   !cm=-((2*pi*((u(i_step)*cos(alpha(i_step))/u_ref)+(hdot(i_step)*sin(alpha(i_step))/u_ref))&
   !     &*(((aterm(0)/4)+(aterm(1)/4)-(aterm(2)/8))&
   !     &-(cm_pvt*(aterm(0)+(aterm(1)/2)))))+&
   !     &((pi*c/u_ref)*(((7*adot(0)/8)+(3*adot(1)/8)+(adot(2)/8)-(adot(3)/32))&
   !     &-(cm_pvt*((3*adot(0)/2)+(adot(1)/2)+(adot(2)/4)))))+(nonl_m))
   cm=cn*cm_pvt-(2*pi*(((u(i_step)*cos(alpha(i_step))/u_ref)+(hdot(i_step)*sin(alpha(i_step))/u_ref))&
        &*((aterm(0)/4)+(aterm(1)/4)-(aterm(2)/8))+(c/u_ref)*((7*adot(0)/16)+(3*adot(1)/16)&
        &+(adot(2)/16)-(adot(3)/64))))-nonl_m


   !Nondimensional quantities and coefficients are output
   write(17,999)t(i_step),alpha(i_step)*180/pi,h(i_step),u(i_step)&
        &,bound_circ,aterm(0),cn,cs,cl,cd,cm
999 format (E16.8,E16.8,E16.8,E16.8,E16.8,E16.8,E16.8,&
         &E16.8,E16.8,E16.8,E16.8,E16.8,E16.8,E16.8)

   if (flow_filename .ne. 'nil') then
      if (mod(i_step,n_pts_flow)==0) then
         do i_lev=1,n_lev
            write(13,*)lev(i_lev,1),lev(i_lev,2),lev(i_lev,3)
         end do
         do i_tev=1,n_tev
            write(13,*)tev(i_tev,1),tev(i_tev,2),tev(i_tev,3)
         end do
         do i_div=2,n_div
            write(13,*)bound_int(i_div,1),bound_int(i_div,2),bound_int(i_div,3)
         end do
         write(13,*)'NaN ','NaN ','NaN '
      end if
   end if

end do

call cpu_time(t2)
write(*,*)t2-t1,n_lev, n_tev

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine splint(xa,ya,y2a,n,x,y)
integer n
double precision x,y,xa(n),y2a(n),ya(n)
integer k,khi,klo
double precision a,b,h
klo=1
khi=n

1 if (khi-klo.gt.1) then
  k=(khi+klo)/2
  if (xa(k).gt.x) then
    khi=k
  else
    klo=k
  end if
goto 1
end if
h=xa(khi)-xa(klo)
if (h.eq.0) write(*,*) 'bad xa input'
a=(xa(khi)-x)/h
b=(x-xa(klo))/h
y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6
return
end subroutine splint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine spline(x,y,n,yp1,ypn,y2) !From numerical recipes
integer n,nmax
double precision :: yp1,ypn,x(n),y(n),y2(n)
parameter(nmax=3000)
integer i,k
double precision p,qn,sig,un,u(nmax)
if (yp1.gt..99e30) then
  y2(1)=0
  u(1)=0
else
  y2(1)=-0.5
  u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1  )
end if
do i=2,n-1
  sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
  p=sig*y2(i-1)+2
  y2(i)=(sig-1.)/p
  u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))&
  &/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
end do
if (ypn.gt..99e30) then
  qn=0
  un=0
else
  qn=0.5
  un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
end if
y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
do k=n-1,1,-1
  y2(k)=y2(k)*y2(k+1)+u(k)
end do
return
end subroutine spline
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_camberslope
  !Constructing camber slope from airfoil file
  if (foil_name=='flat_plate') then
     cam_slope(:)=0
  else
     open(unit=5, file=foil_name, status='unknown')
     i_coord=0
     do
        read(5,*,iostat=io1)xcoord(i_coord+1),ycoord(i_coord+1)
        if (io1<0) then
           n_coord=i_coord
           exit
        end if
        i_coord=i_coord+1
     end do
     xcoord_sum(1)=0
     do i_coord=2,n_coord
        xcoord_sum(i_coord)=xcoord_sum(i_coord-1)+abs(xcoord(i_coord)-xcoord(i_coord-1))
     end do
     call spline(xcoord_sum,ycoord,n_coord,dble(0),dble(0),ysplined)
     do i_div=1,n_div
        xreq=x(i_div)/c
        !xreq goes from 0 to 1
        call splint(xcoord_sum,ycoord,ysplined,n_coord,xreq,ycoord_ans(i_div))
     end do
     do i_div=n_div+1,2*n_div
        xreq=(x(n_div)/c)+(x(i_div-n_div)/c)
        !xreq goes from 1 to 2
        call splint (xcoord_sum,ycoord,ysplined,n_coord,xreq,ycoord_ans(i_div))
     end do
     do i_div=n_div,1,-1
        cam(n_div+1-i_div)=(ycoord_ans(i_div)+ycoord_ans((2*n_div)+1-i_div))/2
     end do
     !Scale the camber according to chord length
     cam(1:n_div)=cam(1:n_div)*c
     cam(1) = 0.0
     !Determine camber slope - both x and cam are in dimensional units
     cam_slope(1)=(cam(2)-cam(1))/(x(2)-x(1))
     do i_div=2,n_div
        cam_slope(i_div)=(cam(i_div)-cam(i_div-1))/(x(i_div)-x(i_div-1))
     end do
  end if
end subroutine calc_camberslope
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_downwash_boundcirc
  uind(1:n_div)=0
  wind(1:n_div)=0
  do i_div=1,n_div
     do i_lev=1,n_lev
        dist=xdist(1,3,i_div,i_lev)**2+zdist(1,3,i_div,i_lev)**2
        uind(i_div)=uind(i_div)+(lev(i_lev,1)*&
             &(-zdist(1,3,i_div,i_lev))/(2*pi*sqrt(v_core**4+dist**2)))
        wind(i_div)=wind(i_div)+(-lev(i_lev,1)*&
             &(-xdist(1,3,i_div,i_lev))/(2*pi*sqrt(v_core**4+dist**2)))
     end do
     do i_tev=1,n_tev
        dist=xdist(1,2,i_div,i_tev)**2+zdist(1,2,i_div,i_tev)**2
        uind(i_div)=uind(i_div)+(tev(i_tev,1)*&
             &(-zdist(1,2,i_div,i_tev))/(2*pi*sqrt(v_core**4+dist**2)))
        wind(i_div)=wind(i_div)+(-tev(i_tev,1)*&
             &(-xdist(1,2,i_div,i_tev))/(2*pi*sqrt(v_core**4+dist**2)))
     end do
     downwash(i_div)=(-u(i_step)*sin(alpha(i_step)))+(-uind(i_div)*sin(alpha(i_step)))&
          &+(hdot(i_step)*cos(alpha(i_step)))+(-wind(i_div)*cos(alpha(i_step)))&
          &+(-alphadot(i_step)*(x(i_div)-pvt*c))+(cam_slope(i_div)*((uind(i_div)&
          &*cos(alpha(i_step)))+(u(i_step)*cos(alpha(i_step)))+(hdot(i_step)*sin(alpha(i_step)))&
          &+(-wind(i_div)*sin(alpha(i_step)))))
  end do
  aterm(0:1)=0
  do i_div=2,n_div
     aterm(0)=aterm(0)+(((downwash(i_div)&
          +downwash(i_div-1))/2)*dtheta)
     aterm(1)=aterm(1)+(((downwash(i_div)&
          *cos(theta(i_div))+downwash(i_div-1)*cos(theta(i_div-1)))/2)*dtheta)
  end do
  aterm(0)=(-1./(u_ref*pi))*aterm(0)
  aterm(1)=(2./(u_ref*pi))*aterm(1)
  bound_circ=u_ref*c*pi*(aterm(0)+(aterm(1)/2.))
end subroutine calc_downwash_boundcirc

end program UNSflow
