program signal
implicit none

integer,parameter :: nf = 14, maxnmds = 20000

real(8),parameter :: time1 = 477, pi= 3.1415926535d0, dts = 5d-5
real(8),parameter :: dt = 2d0*pi*dts
real(8),parameter :: e1 = 0.05, omega = 260

integer,parameter :: tau = 143 

character(len=23) :: fname1
character(len=24) :: fname2
character(len=53) :: fname3,fname4

complex(8) :: integral,integrand,pul
complex(8),dimension(0:maxnmds,0:nf) :: pol, f_pol_oo, f_pol_ox,f_pol_xo
complex(8),dimension(0:maxnmds,0:nf) :: f_pol_xx,pulset,prod

integer :: i,it,j,k,switch,x,y,voidi,nmds,time

real(8),dimension(0:nf) :: f_re1,f_im1,f_re2,f_im2,delay
real(8) :: gaussian,intensity,voidr,treal,timag!,time

pol = dcmplx(0d0,0d0)
f_pol_oo = dcmplx(0d0,0d0)
f_pol_ox = dcmplx(0d0,0d0)
f_pol_xo = dcmplx(0d0,0d0)
f_pol_xx = dcmplx(0d0,0d0)

!goto 100

do i = 0, nf
   x = 700 + i
   y = 800 + i

   write(fname1,'(a9,i2.2,a12)') '../e3-map',i,'/polariz.out'
   open(x,file=fname1)

   write(fname2,'(a9,i2.2,a13)') '../e3-map',i,'/intensity.in'
   open(y,file=fname2)

   read(y,*)
   read(y,*) voidr, nmds, voidr, voidr, voidr
   
   read(x,*)
   read(x,*)
   do it = 1, nmds
      read(x,*) voidi, voidr,treal!, timag
      read(x,*) timag

      f_pol_xx(it,i) = dcmplx(treal,timag)
   end do

   close(y)
   close(x)
end do

do i = 0, nf
   x = 700 + i
   y = 800 + i
   
   write(fname1,'(a9,i2.2,a12)') '../e2-map',i,'/polariz.out'
   open(x,file=fname1)

   write(fname2,'(a9,i2.2,a13)') '../e2-map',i,'/intensity.in'
   open(y,file=fname2)
   
   read(y,*)
   read(y,*) voidr, nmds, voidr, voidr, voidr
   
   read(x,*)
   read(x,*)
   do it = 1, nmds
      read(x,*) voidi, voidr, treal!, timag
      read(x,*) timag

      f_pol_xo(it,i) = dcmplx(treal,timag)
   end do

   close(y)
   close(x)
end do

!100 continue

do i = 0, nf
   x = 700 + i
   y = 800 + i

   write(fname1,'(a9,i2.2,a12)') '../e1-map',i,'/polariz.out'
   open(x,file=fname1)
   
   write(fname2,'(a9,i2.2,a13)') '../e1-map',i,'/intensity.in'
   open(y,file=fname2)

   read(y,*)
   read(y,*) voidr, nmds, voidr, voidr, voidr
   
   read(x,*)
   read(x,*)
   do it = 1, nmds
      read(x,*) voidi, voidr, treal!, timag
      read(x,*) timag

      f_pol_ox(it,i) = dcmplx(treal,timag)
   end do

   close(y)
   close(x)
end do

do i = 0, nf
   x = 700 + i
   y = 800 + i

   write(fname1,'(a9,i2.2,a12)') '../e0-map',i,'/polariz.out'
   open(x,file=fname1)

   write(fname2,'(a9,i2.2,a13)') '../e0-map',i,'/intensity.in'
   open(y,file=fname2)

   read(y,*)
   read(y,*) voidr, nmds, voidr, voidr, voidr
   
   read(x,*)
   read(x,*)
   do it = 1, nmds
      read(x,*) delay(i), voidr, treal!, timag
      read(x,*) timag

      f_pol_oo(it,i) = dcmplx(treal,timag)
   end do

   close(y)
   close(x)
end do

open(102,file='intensity.out')

do i = 0, nmds
   do j = 0, nf
      pol(i,j) =  f_pol_oo(i,j) - f_pol_ox(i,j) + f_pol_xx(i,j) - f_pol_xo(i,j)
   end do
end do

pulset = 0d0

open(103,file='pol000.out')
open(104,file='pol001.out')
open(105,file='pol110.out')
open(106,file='pol111.out')
open(107,file='poltot.out')
do i = 0, nmds
   write(103,'(i6,62f22.16)') i, (dble(f_pol_xx(i,j)),j=0,nf),(aimag(f_pol_xx(i,j)),j=0,nf)
   write(104,'(i6,62f22.16)') i, (dble(f_pol_xo(i,j)),j=0,nf),(aimag(f_pol_xo(i,j)),j=0,nf)
   write(105,'(i6,62f22.16)') i, (dble(f_pol_ox(i,j)),j=0,nf),(aimag(f_pol_ox(i,j)),j=0,nf)
   write(106,'(i6,62f22.16)') i, (dble(f_pol_oo(i,j)),j=0,nf),(aimag(f_pol_oo(i,j)),j=0,nf)
   write(107,'(i6,62f22.16)') i, (dble(pol(i,j)),j=0,nf), (aimag(pol(i,j)),j=0,nf)
end do
close(107)
close(106)
close(105)
close(104)
close(103)

open(108,file='pulse.out')
do i = 0, nf
   integral = dcmplx(0,0)
   
!   time = (0.15 + 0.04*i)!0.99596334853845002735338323789331d0)
   time = delay(i)/dt
!   time = 477 + 127*i
   do it = 1, nmds
      call int_red(pi,tau,it,dt,time,e1,omega,pol(it,i),integrand,switch)
      call pulse(pi,tau,it,dt,time,e1,omega,pul)
      
      prod(it,i) = pul*dconjg(pol(it,i))
      pulset(it,i) = pul

      integral = integral + integrand

   end do

   if (switch == 0) then
      intensity = 2*omega*imag(integral)*dt
   else if (switch == 1) then
      intensity = imag(integral)*dt
   else
      intensity = 2d0*omega*imag(integral)*dt
   end if
  
   write(102,*) (time-time1)*dt/0.00149891, intensity, delay(i)
end do

do it = 0, nmds
   write(108,'(i6,62f20.16)') it, (dble(pulset(it,j)),j=0,nf), (aimag(pulset(it,j)),j=0,nf)
   write(200,'(i6,62f20.16)') it, (dble(prod(it,j)),j=0,nf), (aimag(prod(it,j)),j=0,nf)
end do

close(108)
close(102)

contains

subroutine int_full(pi,tau,it,dt,time,e1,omega,pol,integrand,switch)
implicit none

complex(8),intent(in) :: pol
complex(8),intent(out) :: integrand

integer,intent(in) :: it
integer,intent(out) :: switch

real(8),intent(in) :: pi,tau,dt,time,e1,omega

real(8) :: gaussian,cofac,term

gaussian = dsqrt(4d0*dlog(2d0)/(pi*tau**2))
gaussian = gaussian*dexp(-4d0*dlog(2d0)*((it-0.5d0)*dt - time)**2/tau**2)

cofac = -8d0*dlog(2d0)*((it-0.5d0)*dt - time)/tau**2

term = omega*((it-0.5d0)*dt-time)

integrand = gaussian*e1*dcmplx(dcos(term),-dsin(term))
integrand = integrand*dcmplx(cofac,-omega)
integrand = integrand*pol

switch = 1

end subroutine int_full

subroutine int_red(pi,tau,i,dt,time,e1,omega,pol,integrand,switch)
implicit none

complex(8),intent(in) :: pol
complex(8),intent(out) :: integrand

integer,intent(in) :: i,tau,time
integer,intent(out) :: switch

real(8),intent(in) :: pi,dt,e1,omega!,time

real(8) :: gaussian,it

it = i!*0.994959272d0

gaussian = dsqrt(4d0*dlog(2d0)/(pi*(dt*tau)**2))
gaussian = gaussian*dexp(-4d0*dlog(2d0)*((it-0.5- time))**2/tau**2)

integrand = gaussian*e1* ( dcmplx( dcos(omega*(it-0.5-time)*dt) , -dsin(omega*(it-0.5-time)*dt) ) &
                        )! - dcmplx( dcos(omega*(it-0.5-time)*dt) , -dsin(omega*(it-0.5-time)*dt) ) )

integrand = integrand*dconjg(pol)

switch = 2

end subroutine int_red

subroutine pulse(pi,tau,it,dt,time,e1,omega,pul)
implicit none

integer,intent(in) :: it,tau,time

real(8),intent(in) :: pi,dt,e1,omega!,time

complex(8),intent(out) :: pul

real(8) :: gaussian

gaussian = dsqrt(4d0*dlog(2d0)/(pi*(dt*tau)**2))
gaussian = gaussian*dexp(-4d0*dlog(2d0)*((it-0.5d0- time)**2)/tau**2)
pul = gaussian*e1*dcmplx ( dcos(omega*dt*(it-0.5d0-time)) , -dsin(omega*dt*(it-0.5d0-time)) )
!pul = pul + gaussian*e1*dcmplx ( dcos(omega*dt*(it+0.5d0-time)) , -dsin(omega*dt*(it+0.5d0-time)) )
!pul = pul/2d0
end subroutine pulse

subroutine int_reformatted(pi,tau,it,dt,time,e1,omega,pol,integrand,switch)
implicit none

complex(8),intent(in) :: pol
complex(8),intent(out) :: integrand

integer,intent(in) :: it,time
integer,intent(out) :: switch

real(8),intent(in) :: pi,tau,dt,e1,omega

real(8) :: gaussian

gaussian = dsqrt(4d0*dlog(2d0)/(pi*tau**2))
gaussian = gaussian*dexp(-4d0*dlog(2d0)*((it-0.5d0)*dt - time)**2/tau**2)

integrand = gaussian*e1*dcos(omega*((it-0.5d0)*dt-time))*dconjg(pol)

switch = 0

end subroutine int_reformatted

end program signal
