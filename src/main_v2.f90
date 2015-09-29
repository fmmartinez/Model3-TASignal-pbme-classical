program popmodel3
use ifport, only:srand,rand
implicit none

integer,parameter :: nmap = 3

real(8),parameter :: pi=3.1415926535d0, oc=37.7d0, kc=sqrt(10d0)*oc

complex(8) :: et

real(8) :: coeff,a1,a2,fact
real(8) :: pc,qc,av1,av2,fc,fc1,fc2,hmtrace,etotal,eclas,equan
real(8),dimension(1:3) :: rm,pm
real(8),dimension(1:3,1:3) :: hm
real(8),dimension(:),allocatable :: pol,x,p,fx,polt

integer :: i,j,ng,nb,nd,basispc,stp,cont,p_i,p_j,p_k,omc,nfile,step1
integer :: np,nosc,nmcs,nmds,seed_dimension,bath,init,mcs,it,is,ib,ie,je
integer :: overflow
integer,dimension(:),allocatable :: g

logical :: overflowcheck

real(8) :: delta,ome_max,dt,kondo,lumda_d,eg,eb,ed,mu,e0,e1,beta,check,vomega
real(8) :: dt2,uj,qbeta,lambdacheck,gaussian,rtemp,step2,tau1,omega1,tau2,omega2
real(8) :: time3,arg
real(8),dimension(:),allocatable :: tau,time,omega,ome,c2,kosc

call iniconc()

call srand(seed_dimension)

allocate(ome(1:nosc),c2(1:nosc),kosc(1:nosc))
allocate(tau(1:np),omega(1:np),time(1:np),g(1:np))

call iniconq_d()

do i = 1, np-1
   tau(i)   = tau1
   omega(i) = omega1
end do
tau(np) = tau2
omega(np) = omega2

time(1) = 0.3d0
time(2) = 0.3d0
time(3) = (time3 + nfile*step2)

nmds = nmds + nfile*step1

allocate(pol(1:nmds+1),polt(1:nmds+1))
allocate(x(1:nosc),p(1:nosc),fx(1:nosc))

pol  = cmplx(0d0,0d0)
polt  = cmplx(0d0,0d0)

dt  = 2d0*pi*dt
dt2 = 0.5d0*dt

g(1) = p_i
g(2) = p_j
g(3) = p_k

omc = 0
overflow = 0
overflowcheck = .false.
MC: do mcs = 1, nmcs
!sampling bath oscillators
   do is=1,nosc
      uj = 0.5d0*beta*dsqrt(kosc(is))
      
      qbeta = beta/(uj/tanh(uj))
      
      p(is) = cmplx(gauss_noise2()/dsqrt(qbeta),0d0)
      x(is) = cmplx(gauss_noise2()/dsqrt(qbeta*kosc(is)),0d0)
      
      if(bath == 1) x(is)=x(is)+c2(is)/kosc(is)
   end do
  
!sampling oscillator coupled 
   uj = 0.5d0*beta*dsqrt(oc**2)
   qbeta = beta/(uj/tanh(uj))
   pc = cmplx(gauss_noise2()/dsqrt(qbeta),0d0)
   qc = cmplx(gauss_noise2()/dsqrt(qbeta*oc**2),0d0)
!sampling mapping variables
   if (init == 3) then
      do i = 1, nmap
         rm(i) = cmplx(gauss_noise2()/sqrt(2d0),0d0)
         pm(i) = cmplx(gauss_noise2()/sqrt(2d0),0d0)
      end do
   else
      rm = 0d0
      pm = 0d0
   end if
 
   coeff = rm(1)**2 + pm(1)**2 - 0.5d0

   call get_traceless_force_bath(kosc,x,c2,rm,pm,fx)
   call get_traceless_force_coupledosc(oc,qc,kc,rm,pm,fc1,fc2)
   
   ib = 1

   call get_facts_pol(mu,coeff,rm,pm,fact)
 
   pol(ib) = fact
   
   a1 = cmplx(0d0,0d0)
   a2 = cmplx(0d0,0d0)
   do is=1, nosc 
      a1 = a1 + 2.d0*c2(is)**2/ome(is)**2
      a2 = a2 + 2.d0*c2(is)*x(is)
   end do

   av1 = 2.d0*kc**2/oc**2
   av2 = 2.d0*kc*qc
   
   MD: do it = 1, nmds
      gaussian=sqrt(4d0*log(2d0)/(pi*tau(1)**2))*exp(-4d0*log(2d0)*((it-0.5d0)*dt-time(1))**2/(tau(1)**2))
      arg = omega(1)*((it-0.5d0)*dt-time(1))
      et = g(1)*gaussian*e0*cmplx(cos(arg),sin(arg))
      
      gaussian=sqrt(4d0*log(2d0)/(pi*tau(2)**2))*exp(-4d0*log(2d0)*((it-0.5d0)*dt-time(2))**2/(tau(2)**2))
      arg = omega(2)*((it-0.5d0)*dt-time(2))
      et = et + g(2)*gaussian*e0*cmplx(cos(arg),-sin(arg))
      
      gaussian=sqrt(4d0*log(2d0)/(pi*tau(3)**2))*exp(-4d0*log(2d0)*((it-0.5d0)*dt-time(3))**2/(tau(3)**2))
      arg = omega(3)*((it-0.5d0)*dt-time(3))
      et = et + g(3)*gaussian*e1*cmplx(cos(arg),-sin(arg))
   
      do is = 1, nosc
         p(is) = p(is) + dt2*fx(is)
      end do
      pc = pc + dt2*(fc1+fc2)
      
      call get_hm(delta,mu,et,a1,a2,av1,av2,pc,oc,qc,hm)
      call make_hm_traceless(hm,hmtrace)
      
      call evolve_pm(nmap,dt2,hm,rm,pm)

      do is = 1, nosc
         x(is) = x(is) + dt*p(is)
      end do
      qc = qc + dt*pc

      a2=cmplx(0d0,0d0)
      do is = 1, nosc
          a2 = a2 + 2.d0*c2(is)*x(is)
      end do
      av2 = 2.d0*kc*qc 

      call get_hm(delta,mu,et,a1,a2,av1,av2,pc,oc,qc,hm)
      call make_hm_traceless(hm,hmtrace)

      call evolve_rm(nmap,dt,hm,pm,rm)

      call evolve_pm(nmap,dt2,hm,rm,pm)

      call get_traceless_force_bath(kosc,x,c2,rm,pm,fx)
      call get_traceless_force_coupledosc(oc,qc,kc,rm,pm,fc1,fc2) 

      do is = 1, nosc
         p(is) = p(is) + dt2*fx(is)
      end do
      pc = pc + dt2*(fc1+fc2)

      ib = it + 1
      
      call get_facts_pol(mu,coeff,rm,pm,fact)
      
      polt(ib)  = fact

      if ((polt(ib) /= polt(ib)).or.(polt(ib)-1 == polt(ib))) then
         print *, it, qc
         print *, hm
         print *, rm
         print *, pm
         print *, 'overflow occurred at', mcs
         overflow = overflow + 1
         overflowcheck = .true.
         exit
      end if

!      if (mcs == 3) then
!         etotal = 0d0
!         do is = 1, nosc
!            etotal = 0.5d0*(p(is)**2 + kosc(is)*x(is)**2)
!         end do
!         eclas = etotal
!         etotal = etotal + hmtrace/3d0
!         do ie = 1, 3
!            do je = 1, 3
!               etotal = etotal + 0.5d0*hm(ie,je)*(pm(ie)*pm(je) + rm(ie)*rm(je))
!            end do
!         end do
!         equan = etotal - hmtrace/3d0 - eclas
!         write(69,'(i5,4f20.12)') it,eclas,hmtrace/3d0,equan,etotal
!         write(70,'(i5,9f20.12)') it, hm
!         write(71,'(i5,6f20.12)') it, rm, pm
!         if (it == nmds) stop
!      end if
   end do MD

   if (overflowcheck == .false.) then
      pol = pol + polt
      omc = omc + 1
   else
      overflowcheck = .false.
   end if
end do MC

open(333,file='polariz.out')
do ib = 1, nmds+1
   pol(ib) = pol(ib)/dble(omc)
   write(333,*) time(3), ib-1, dble(pol(ib))!, aimag(pol(ib))
end do
close(333)

print *, 'MC cycles', nmcs, 'with', overflow, 'overflows'

deallocate(ome,c2,kosc)
deallocate(pol,polt)
deallocate(x,p)

contains

subroutine evolve_rm(nmap,dt,hm,pm,rm)
implicit none

integer :: i, j
integer,intent(in) :: nmap

real(8),intent(in) :: dt
real(8),dimension(:),intent(in) :: pm
real(8),dimension(:),intent(inout) :: rm
real(8),dimension(:,:),intent(in) :: hm

do i = 1, nmap
   do j = 1, nmap
      rm(i) = rm(i) + dt*hm(i,j)*pm(j)
   end do
end do

end subroutine evolve_rm

subroutine evolve_pm(nmap,dt2,hm,rm,pm)
implicit none

integer :: i, j
integer,intent(in) :: nmap

real(8),intent(in) :: dt2
real(8),dimension(:),intent(in) :: rm
real(8),dimension(:),intent(inout) :: pm
real(8),dimension(:,:),intent(in) :: hm

do i = 1, nmap
   do j = 1, nmap
      pm(i) = pm(i) - dt2*hm(i,j)*rm(j)
   end do
end do

end subroutine evolve_pm

subroutine get_facts_pol(mu,coeff,rm,pm,fact)
implicit none

real(8),intent(in) :: coeff
real(8),intent(out) :: fact
real(8),dimension(:),intent(in) :: rm,pm

real(8), intent(in) :: mu

fact = mu*coeff*2d0*(rm(1)*rm(2) + pm(1)*pm(2))

end subroutine get_facts_pol

function kronecker_delta(i,j) result (d)
implicit none

integer,intent(in) :: i, j

real(8) :: d

if (i == j) then
   d = 1d0
else
   d = 0d0
end if

end function kronecker_delta

subroutine get_force_bath(kosc,x,c2,rm,pm,f)
implicit none

integer :: j,n

real(8),dimension(:),intent(in) :: kosc,x,c2,rm,pm
real(8),dimension(:),intent(out) :: f

n = size(x)

f = 0d0
do j = 1, n
   f(j) = -kosc(j)*x(j) - c2(j)*(rm(3)**2 + pm(3)**2 - 1d0)
end do

end subroutine get_force_bath

subroutine get_traceless_force_bath(kosc,x,c2,rm,pm,f)
implicit none

integer :: j,n

real(8),dimension(:),intent(in) :: x,rm,pm
real(8),dimension(:),intent(out) :: f

real(8) :: trace
real(8),dimension(:),intent(in) :: kosc,c2

n = size(x)

f = 0d0
do j = 1, n
   trace = -2d0*c2(j)/3d0
   f(j) = -kosc(j)*x(j) - trace*(0.5d0*(rm(1)**2 + pm(1)**2 + rm(2)**2 + pm(2)**2 - 2d0*rm(3)**2 - 2d0*pm(3)**2) - 1d0)
end do

end subroutine get_traceless_force_bath

subroutine get_force_coupledosc(oc,qc,kc,rm,pm,f1,f2)

real(8),intent(in) :: oc,qc,kc
real(8),intent(in),dimension(:) :: rm,pm
real(8),intent(out) :: f1,f2

!f1 = (-oc**2*qc)*(rm(1)**2+rm(2)**2+rm(3)**2+pm(1)**2+pm(2)**2+pm(3)**2-3d0) 
!f2 = (rm(1)**2 + pm(1)**2 - 1d0)
!f2 = f2 + (2d0*kc)*(rm(2)**2 + pm(2)**2 - 1d0)
!f2 = f2 + (kc)*(rm(3)**2 + pm(3)**2 - 1d0)

f1 = -oc**2*qc
f2 = kc*(2d0*(rm(2)**2+pm(2)**2-1d0) + (rm(3)**2+pm(3)**2-1d0))

end subroutine get_force_coupledosc

subroutine get_traceless_force_coupledosc(oc,qc,kc,rm,pm,f1,f2)

real(8),intent(in) :: oc,kc

real(8),intent(in) :: qc
real(8),intent(in),dimension(:) :: rm,pm
real(8),intent(out) :: f1,f2

f1 = -oc**2*qc + kc
f2 = kc*((rm(2)**2+pm(2)**2-1d0) + (rm(1)**2+pm(1)**2-1d0))

end subroutine get_traceless_force_coupledosc

!subroutine update_hm(a1,a2,av1,av2,pc,oc,qc,hm)
!implicit none
!
!real(8),parameter :: eg=0, eb=240, ed=240
!
!real(8) :: ev
!real(8),intent(in) :: a1,a2,av1,av2,pc,oc,qc
!real(8),dimension(:,:),intent(inout) :: hm
!
!ev = 0.5d0*(pc**2 + (oc*qc)**2)
!
!!only diagonal part is updated
!!1 x 1
!hm(1,1) = 0d0
!hm(1,1) = eg + ev
!!2 x 2
!hm(2,2) = 0d0
!hm(2,2) = eb + ev + av1 - av2
!!3 x 3
!hm(3,3) = 0d0
!hm(3,3) = ed + ev + a1 + a2 + 0.25d0*av1 - 0.5d0*av2
!
!end subroutine update_hm

subroutine get_hm(delta,mu,et,a1,a2,av1,av2,pc,oc,qc,hm)
implicit none

real(8),parameter :: eg=0, eb=240, ed=240

complex(8) :: et

real(8) :: ev
real(8),intent(in) :: a1,a2,pc,qc,av1,av2
real(8),dimension(:,:),intent(out) :: hm

real(8),intent(in) :: delta,mu,oc

ev = 0.5d0*(pc**2 + (oc*qc)**2)

hm = cmplx(0d0,0d0)
!1 x 1
hm(1,1) = eg + ev
!1 x 2
hm(1,2) = -mu*et
!1 x 3
!this part is zero
!2 x 1
hm(2,1) = hm(1,2)
!2 x 2
hm(2,2) = eb + ev + av1 - av2
!2 x 3
hm(2,3) = delta
!3 x 1
!this part is zero
!3 x 2
hm(3,2) = hm(2,3)
!3 x 3
hm(3,3) = ed + ev + a1 + a2 + 0.25d0*av1 - 0.5d0*av2

end subroutine get_hm

subroutine make_hm_traceless(hm,trace)
implicit none

real(8),intent(inout) :: trace
real(8),dimension(:,:),intent(inout) :: hm

trace = hm(1,1) + hm(2,2) + hm(3,3)
hm(1,1) = hm(1,1) - trace/3d0
hm(2,2) = hm(2,2) - trace/3d0
hm(3,3) = hm(3,3) - trace/3d0
end subroutine make_hm_traceless

subroutine iniconc()
implicit none

integer :: i

open (666,file='map.in')
read(666,*)
read(666,*) np,delta,kondo,nosc,ome_max
read(666,*)
read(666,*) nmcs,nmds,seed_dimension,dt,lumda_d
read(666,*)
read(666,*) eg,eb,ed,mu,e0,e1,beta,vomega
read(666,*)
read(666,*) tau1,omega1,tau2,omega2,time3,step1,step2
read(666,*)
read(666,*) bath,init,nfile
read(666,*)
read(666,*) ng,nb,nd
read(666,*)
read(666,*) p_i, p_j, p_k
close(666)

!call random_seed(size=seed_dimension)
!allocate (seed(seed_dimension))
!do i=1,seed_dimension
!  seed(i) = 3*2**i-1
!enddo
!call random_seed(put=seed)

end subroutine iniconc

subroutine iniconq_d()
implicit none

integer :: i

check=0
do i=1,nosc
!  ome(i)=tan(pi*i/(4*nosc))
   ome(i)=tan(i*datan(ome_max)/nosc)
!  ome(i)=i**2*ome_max/nosc**2
!  rho(i)=sqrt(ome(i)*ome_max)/(1+ome(i)**2)
!  c2(i)=0.5d0*ome(i)*dsqrt(lumda_d/nosc)
   c2(i)=ome(i)*dsqrt(datan(ome_max)*lumda_d/(pi*nosc))
!  c2(i)=ome(i)*sqrt(lumda_d*rho(i)*2/(pi*nosc))
   check=check+c2(i)**2/ome(i)**2
   kosc(i)=ome(i)**2 
end do
write(6,*) check

end subroutine iniconq_d

function gauss_noise2() result(g)
implicit none

real(8),parameter :: pi2=2.0*3.141592654

real(8) :: g,z1,z2

!call random_number(z1)
!call random_number(z2)
z1 = rand()
z2 = rand()
g = sqrt(-2.d0*log(z1))*cos(pi2*z2)

end function gauss_noise2

end program popmodel3
