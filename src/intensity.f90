program signal
implicit none

integer,parameter :: maxnmds = 20000

real(8),parameter :: time1 = 0.30, pi= 3.1415926535d0, dts = 5d-5
real(8),parameter :: dt = 2d0*pi*dts

character(len=1) :: dir
character(len=20) :: fname1
character(len=21) :: fname2
character(len=5) :: fname3,fname4

complex(8) :: integral,integrand,pul
complex(8),dimension(:,:),allocatable :: pol, f_pol_000, f_pol_010, f_pol_001
complex(8),dimension(:,:),allocatable :: f_pol_011, f_pol_100, f_pol_110, f_pol_101
complex(8),dimension(:,:),allocatable :: f_pol_111,pulset,prod

integer :: i,it,j,k,switch,x,y,voidi,nmds,nf

real(8),dimension(:),allocatable :: f_re1,f_im1,f_re2,f_im2,delay, time, tau
real(8) :: gaussian,intensity,voidr,treal,timag,omega,e1

print *, 'number of points to resolve'
read(*,*) nf

allocate(pol(0:maxnmds,0:nf),f_pol_000(0:maxnmds,0:nf),f_pol_010(0:maxnmds,0:nf))
allocate(f_pol_001(0:maxnmds,0:nf),f_pol_011(0:maxnmds,0:nf),f_pol_100(0:maxnmds,0:nf))
allocate(f_pol_110(0:maxnmds,0:nf),f_pol_101(0:maxnmds,0:nf),f_pol_111(0:maxnmds,0:nf))
allocate(pulset(0:maxnmds,0:nf),prod(0:maxnmds,0:nf))

allocate(f_re1(0:nf),f_im1(0:nf),f_re2(0:nf),f_im2(0:nf))
allocate(delay(0:nf),time(0:nf),tau(0:nf))

pol = dcmplx(0d0,0d0)
f_pol_000 = dcmplx(0d0,0d0)
f_pol_010 = dcmplx(0d0,0d0)
f_pol_001 = dcmplx(0d0,0d0)
f_pol_011 = dcmplx(0d0,0d0)
f_pol_100 = dcmplx(0d0,0d0)
f_pol_110 = dcmplx(0d0,0d0)
f_pol_101 = dcmplx(0d0,0d0)
f_pol_111 = dcmplx(0d0,0d0)

!goto 100
do j = 0, 7
   do i = 0, nf
      x = 700 + i
      y = 800 + i
      
      write(dir,'(i1.1)') j
      
      write(fname1,'(a1,i1.1,a4,i2.2,a12)') 'e',j,'-sec',i,'/polariz.out'
      open(x,file=fname1)

      write(fname2,'(a1,i1.1,a4,i2.2,a13)') 'e',j,'-sec',i,'/intensity.in'
      open(y,file=fname2)

      read(y,*)
      read(y,*) e1, nmds, tau(i), omega, time(i)
      
      read(x,*)
      read(x,*)
      do it = 1, nmds
         read(x,*) delay(i), voidr,treal!, timag
         read(x,*) timag
         
         select case (j)
         case(0)   
            f_pol_000(it,i) = dcmplx(treal,timag)
         case(1)
            f_pol_010(it,i) = dcmplx(treal,timag)
         case(2)
            f_pol_001(it,i) = dcmplx(treal,timag)
         case(3)
            f_pol_011(it,i) = dcmplx(treal,timag)
         case(4)
            f_pol_100(it,i) = dcmplx(treal,timag)
         case(5)
            f_pol_110(it,i) = dcmplx(treal,timag)
         case(6)
            f_pol_101(it,i) = dcmplx(treal,timag)
         case(7)
            f_pol_111(it,i) = dcmplx(treal,timag)
         case default
            print *, 'reading error'
         end select
      end do

      close(y)
      close(x)
   end do
end do

open(102,file='intensity.log')

do i = 0, nmds
   do j = 0, nf
      !ALL
      pol(i,j) = -f_pol_000(i,j)
      pol(i,j) = pol(i,j) + f_pol_100(i,j) + f_pol_010(i,j) + f_pol_001(i,j)
      pol(i,j) = pol(i,j) - (f_pol_110(i,j) + f_pol_101(i,j) + f_pol_011(i,j))
      pol(i,j) = pol(i,j) + f_pol_111(i,j)
      !RWA
      !pol(i,j) = -f_pol_000(i,j) + f_pol_111(i,j) - f_pol_110(i,j) - f_pol_101(i,j)
   end do
end do

pulset = 0d0

do j = 0, 7
   select case (j)
   case(0)
      open(103,file='000.out')
      do it = 0, nmds
         write(103,'(i6,62f22.16)') it, (dble(f_pol_000(it,i)),i=0,nf),(aimag(f_pol_000(it,i)),i=0,nf)
      end do
      close(103)
   case(1)
      open(103,file='010.out')
      do it = 0, nmds
         write(103,'(i6,62f22.16)') it, (dble(f_pol_010(it,i)),i=0,nf),(aimag(f_pol_010(it,i)),i=0,nf)
      end do
      close(103)
   case(2)
      open(103,file='001.out')
      do it = 0, nmds
         write(103,'(i6,62f22.16)') it, (dble(f_pol_001(it,i)),i=0,nf),(aimag(f_pol_001(it,i)),i=0,nf)
      end do
      close(103)
   case(3)
      open(103,file='011.out')
      do it = 0, nmds
         write(103,'(i6,62f22.16)') it, (dble(f_pol_011(it,i)),i=0,nf),(aimag(f_pol_011(it,i)),i=0,nf)
      end do
      close(103)
   case(4)
      open(103,file='100.out')
      do it = 0, nmds
         write(103,'(i6,62f22.16)') it, (dble(f_pol_100(it,i)),i=0,nf),(aimag(f_pol_100(it,i)),i=0,nf)
      end do
      close(103)
   case(5)
      open(103,file='110.out')
      do it = 0, nmds
         write(103,'(i6,62f22.16)') it, (dble(f_pol_110(it,i)),i=0,nf),(aimag(f_pol_110(it,i)),i=0,nf)
      end do
      close(103)
   case(6)
      open(103,file='101.out')
      do it = 0, nmds
         write(103,'(i6,62f22.16)') it, (dble(f_pol_101(it,i)),i=0,nf),(aimag(f_pol_101(it,i)),i=0,nf)
      end do
      close(103)
   case(7)
      open(103,file='111.out')
      do it = 0, nmds
         write(103,'(i6,62f22.16)') it, (dble(f_pol_111(it,i)),i=0,nf),(aimag(f_pol_111(it,i)),i=0,nf)
      end do
      close(103)
   case default
      print *, 'error'
   end select
end do

pul = cmplx(0d0,0d0)
open(108,file='pulse.out')
do i = 0, nf
   integral = dcmplx(0,0)
   
   do it = 1, nmds
      call pulse(pi,tau(i),it,dt,time(i),e1,omega,pul)
      call int_red(pul,pol(it,i),integrand,switch)
!      call int_reformatted(pi,tau(i),it,dt,time(i),e1,omega,pol(it,i),integrand,switch)

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
  
   write(102,*) (time(i)-time1)/0.00149891, intensity, delay(i)
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

subroutine int_red(field,pol,integrand,switch)
implicit none

complex(8),intent(in) :: pol,field
complex(8),intent(out) :: integrand

integer,intent(out) :: switch

integrand = field*dconjg(pol)

switch = 2

end subroutine int_red

subroutine pulse(pi,tau,it,dt,time,e1,omega,pul)
implicit none

integer,intent(in) :: it

real(8),intent(in) :: pi,dt,e1,omega,time,tau

complex(8),intent(out) :: pul

real(8) :: gaussian,ln2t4,arg

arg = omega*((it - 0.5d0)*dt - time)
ln2t4 = 4d0*dlog(2d0)

gaussian = dsqrt(ln2t4/(pi*(tau**2)))*dexp(-ln2t4*((arg/omega)**2)/tau**2)

pul = gaussian*e1*dcmplx(dcos(arg),-dsin(arg))
end subroutine pulse

subroutine int_reformatted(pi,tau,it,dt,time,e1,omega,pol,integrand,switch)
implicit none

complex(8),intent(in) :: pol
complex(8),intent(out) :: integrand

integer,intent(in) :: it
integer,intent(out) :: switch

real(8),intent(in) :: pi,tau,dt,e1,omega,time

real(8) :: gaussian,ln2t4,arg

arg = omega*((it - 0.5d0)*dt - time)
ln2t4 = 4d0*log(2d0)

gaussian = sqrt(ln2t4/(pi*tau**2))*exp(-ln2t4*((arg/omega)**2)/tau**2)

integrand = gaussian*e1*cos(arg)*dconjg(pol)

switch = 0

end subroutine int_reformatted

end program signal
