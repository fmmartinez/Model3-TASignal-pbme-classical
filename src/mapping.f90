module m_map
implicit none

private

public iniconq_d,get_preh,sampling_class,sampling_mapng,get_coeff,get_fact,get_a
public get_traceless_force_bath, get_traceless_force_coupledosc
public get_pulsefield
public get_hm,make_hm_traceless
public update_p,update_x,update_pm,update_rm,update_a2
public get_total_energy

real(8),parameter :: pi=3.1415926535d0

contains

subroutine get_hm(delta,mu,et,a1,a2,av1,av2,pc,oc,qc,hm)
implicit none

real(8),parameter :: eg=0, eb=240, ed=240

complex(8) :: ev
complex(8),intent(in) :: et,a1,a2,qc,av1,av2
complex(8),dimension(:,:),intent(out) :: hm

real(8),intent(in) :: delta,mu,pc,oc

ev = 0.5d0*(pc**2 + (oc*qc)**2)

hm = 0d0
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

subroutine make_hm_traceless(hm,tracen)
implicit none

complex(8),intent(inout) :: tracen
complex(8),dimension(:,:),intent(inout) :: hm

tracen = (hm(1,1) + hm(2,2) + hm(3,3))/3d0
hm(1,1) = hm(1,1) - tracen
hm(2,2) = hm(2,2) - tracen
hm(3,3) = hm(3,3) - tracen
end subroutine make_hm_traceless

subroutine get_traceless_force_coupledosc(oc,qc,kc,rm,pm,f1,f2)
implicit none

complex(8),intent(in) :: qc
complex(8),intent(in),dimension(:) :: rm,pm
complex(8),intent(out) :: f1,f2

real(8),intent(in) :: oc,kc

f1 = -oc**2*qc + kc
f2 = kc*((rm(2)**2+pm(2)**2-1d0) + (rm(1)**2+pm(1)**2-1d0))

end subroutine get_traceless_force_coupledosc

subroutine get_traceless_force_bath(kosc,x,c2,rm,pm,f)
implicit none

complex(8),dimension(:),intent(out) :: f
complex(8),dimension(:),intent(in) :: x,rm,pm

integer :: j,n

real(8) :: trace
real(8),dimension(:),intent(in) :: kosc,c2

n = size(kosc)

f = 0d0
do j = 1, n
   trace = -2d0*c2(j)/3d0
   f(j) = -kosc(j)*x(j) - trace*(0.5d0*(rm(1)**2 + pm(1)**2 + rm(2)**2 + pm(2)**2 - 2d0*rm(3)**2 - 2d0*pm(3)**2) - 1d0)
end do

end subroutine get_traceless_force_bath


subroutine get_total_energy(nosc,nmap,kosc,p,x,hm,trace,rm,pm,h,hcl,hma)
implicit none

complex(8),intent(out) :: h,hcl,hma
complex(8),intent(in) :: trace
complex(8),intent(in),dimension(:) :: x,p,rm,pm
complex(8),intent(in),dimension(:,:) :: hm

integer :: i,j
integer,intent(in) :: nmap,nosc

real(8),intent(in),dimension(:) :: kosc

h = cmplx(0d0,0d0)

!classical
do i = 1, nosc
   h = h + 0.5d0*(p(i)**2 + kosc(i)*x(i)**2)
end do

hcl = h

!mapping
do i = 1, nmap
   do j = 1, nmap
      h = h + 0.5d0*hm(i,j)*(rm(i)*rm(j) + pm(i)*pm(j))
   end do
end do

h = h + trace

hma = h - hcl

end subroutine get_total_energy

subroutine get_coeff(ng,beta,omega,rm,pm,coeff)
implicit none

complex(8),intent(out) :: coeff
complex(8),dimension(:),intent(in) :: rm,pm

integer :: i
integer,intent(in) :: ng

real(8),intent(in) :: beta,omega
real(8) :: z
real(8),dimension(:),allocatable :: exp_be,prob

allocate(exp_be(1:ng),prob(1:ng))

exp_be = 0d0
z = 0d0
do i = 1, ng
   !exp_be(i) = exp(-beta*omega*(i - 0.5d0))
   !exp_be(i) = exp(-1.44d0*(i - 0.5d0))
   exp_be(i) = exp(-2.29d-1*(i - 0.5d0))
   z = z + exp_be(i)
end do

prob = exp_be/z

coeff = cmplx(0d0,0d0)
!only diagonal because in the ngxng matrix by construction their lambdas are
! (1,0,0...), (0,1,0...), (0,0,1...), etc
do i = 1, ng
   coeff = coeff + (rm(i)**2 + pm(i)**2 - 0.5d0)*prob(i)
end do

deallocate(exp_be)
deallocate(prob)
end subroutine get_coeff

subroutine get_fact(ng,nb,coeff,llgb,llbg,mu,rm,pm,fact)
implicit none

complex(8),intent(in) :: coeff
complex(8),intent(out) :: fact
complex(8),dimension(:),intent(in) :: rm,pm

integer :: a,b
integer,intent(in) :: ng,nb

real(8),intent(in) :: mu
real(8),dimension(:,:),intent(in) :: llgb,llbg

fact = cmplx(0d0,0d0)

do a = 1,ng
   do b = ng+1,ng+nb
      fact = fact + llgb(a,b)*(rm(a)*rm(b) + pm(a)*pm(b))
   end do
end do

do a = ng+1,ng+nb
   do b = 1,ng
      fact = fact + llbg(a,b)*(rm(a)*rm(b) + pm(a)*pm(b))
   end do
end do

fact = fact*coeff*mu
end subroutine get_fact

subroutine get_lambda_eigenvectors(ng,nb,nd,eg,eb,ed,delta,omega,&
                                    sgg,sgb,sgd,sbg,sbb,sbd,sdg,sdb,sdd,lambda,hsf)
use m_vib
implicit none

integer,parameter :: ip = 10000

character(len=2) :: c_nb,c_nt
character(len=9) :: fmt1,fmt2

integer :: i, j, nt
integer,intent(in) :: nd,ng,nb

real(8) :: omgsq,cg,cb,cd,uint,lint,alpha
real(8),intent(in) :: eg,eb,ed,delta,omega
real(8),dimension(1:3,1:3) :: hs
real(8),dimension(:),allocatable :: eigval
real(8),dimension(:,:),intent(out) :: lambda,sgg,sgb,sgd,sbg,sbb,sbd,sdg,sdb,sdd,hsf
real(8),dimension(:,:),allocatable :: kg,kb,kd,vg,vb,vd,hxp,hxp_temp,sxp

!lapack variables
integer,parameter :: lwork = 30000
real(8),dimension(1:lwork) :: work
integer :: info

!generating some constants
nt = ng + nb + nd

allocate(kg(1:ng,1:ng),kb(1:nb,1:nb),kd(1:nd,1:nd))
allocate(vg(1:ng,1:ng),vb(1:nb,1:nb),vd(1:nd,1:nd))
allocate(hxp(1:nt,1:nt),hxp_temp(1:nt,1:nt),sxp(1:nt,1:nt))
allocate(eigval(1:nt))

omgsq = omega**2

cg = 0d0
cb = 2d0*sqrt(10d0)/(omega)
cd = cb/2d0

uint = cb + 5.5d0
lint = cg - 5.5d0

alpha = sqrt(omega)

!original subsystem hamiltonian
hs = 0d0
hs(1,1) = eg
hs(2,2) = eb
hs(2,3) = delta
hs(3,2) = delta
hs(3,3) = ed

!initialization
kg = 0d0
vg = 0d0
sgg = 0d0

sgb = 0d0

sgd = 0d0
!
sbg = 0d0

kb = 0d0
vb = 0d0
sbb = 0d0

sbd = 0d0
!
sdg = 0d0

sdb = 0d0

kd = 0d0
vd = 0d0
sdd = 0d0

!construction of integrals for the extended hamiltonian
! 1,1
do i = 1, ng
   do j = 1, ng
      kg(i,j) = integrate_t_phid2p(ip,lint,uint,i,cg,j,cg,alpha)
      vg(i,j) = integrate_t_pqcsqp(ip,lint,uint,i,cg,j,cg,alpha,cg)
   end do
   sgg(i,i) = 1d0
end do
! 1,2
do i = 1, ng
   do j = 1, nb
      sgb(i,j) = integrate_t_phiphi(ip,lint,uint,i,cg,j,cb,alpha)
   end do
end do
! 1,3
! stays at 0, will not be used
! 2,1
do i = 1, nb
   do j = 1, ng
      sbg(i,j) = integrate_t_phiphi(ip,lint,uint,i,cb,j,cg,alpha)
   end do
end do
! 2,2
do i = 1, nb
   do j = 1, nb
      kb(i,j) = integrate_t_phid2p(ip,lint,uint,i,cb,j,cb,alpha)
      vb(i,j) = integrate_t_pqcsqp(ip,lint,uint,i,cb,j,cb,alpha,cb)
   end do
   sbb(i,i) = 1d0
end do
! 2,3
do i = 1, nb
   do j = 1, nd
      sbd(i,j) = integrate_t_phiphi(ip,lint,uint,i,cb,j,cd,alpha)
   end do
end do
! 3,1
! stays at 0, will not be used
! 3,2
do i = 1, nd
   do j = 1, nb
      sdb(i,j) = integrate_t_phiphi(ip,lint,uint,i,cd,j,cb,alpha)
   end do
end do
! 3,3
do i = 1, nd
   do j = 1, nd
      kd(i,j) = integrate_t_phid2p(ip,lint,uint,i,cd,j,cd,alpha)
      vd(i,j) = integrate_t_pqcsqp(ip,lint,uint,i,cd,j,cd,alpha,cd)
   end do
   sdd(i,i) = 1d0
end do

!expansion of subsystem hamiltonian with basis functions
hxp(1:ng,1:ng)       = hs(1,1)*sgg(1:ng,1:ng) + 0.5d0*(kg(1:ng,1:ng) + omgsq*vg(1:ng,1:ng))
hxp(1:ng,ng+1:ng+nb) = hs(1,2)*sgb(1:ng,1:nb)
hxp(1:ng,ng+nb+1:nt) = hs(1,3)*sgd(1:ng,1:nd)

hxp(ng+1:ng+nb,1:ng)       = hs(2,1)*sbg(1:nb,1:ng)
hxp(ng+1:ng+nb,ng+1:ng+nb) = hs(2,2)*sbb(1:nb,1:nb) + 0.5d0*(kb(1:nb,1:nb) + omgsq*vb(1:nb,1:nb))
hxp(ng+1:ng+nb,ng+nb+1:nt) = hs(2,3)*sbd(1:nb,1:nd)

hxp(ng+nb+1:nt,1:ng)       = hs(3,1)*sdg(1:nd,1:ng)
hxp(ng+nb+1:nt,ng+1:ng+nb) = hs(3,2)*sdb(1:nd,1:nb)
hxp(ng+nb+1:nt,ng+nb+1:nt) = hs(3,3)*sdd(1:nd,1:nd) + 0.5d0*(kd(1:nd,1:nd) + omgsq*vd(1:nd,1:nd))

!overlap?
sxp(1:nt,1:nt) = 0d0
do i = 1, nt
   sxp(i,i) = 1d0
end do

if (nb > 9) then
   write(c_nb,'(i2)') nb
else
   write(c_nb,'(i1)') nb
end if

if (nt > 9) then
   write(c_nt,'(i2)') nt
else
   write(c_nt,'(i1)') nt
end if

fmt1 = '('//trim(c_nb)//'f10.5)'
fmt2 = '('//trim(c_nt)//'f10.5)'

write(*,*) 'hxp'
write(*,fmt2) hxp

write(*,*) 'sbg'
write(*,fmt1) sbg

!Diagonalization
hxp_temp = hxp
!call dsyev('v','u',nt,hxp_temp,nt,eigval,work,lwork,info)
call dsygv(1,'v','u',nt,hxp_temp,nt,sxp,nt,eigval,work,lwork,info)
if (info /=0) then
   print *, info, 'info from Diagonalization Hc=ec'
   stop
end if

write(*,*) 'diagonalization result'
write(*,*) 'eigenvectors'
do i = 1, nt
   write(*,*) i, 'E=',eigval(i)
   write(*,fmt2) hxp_temp(1:nt,i)
end do

lambda(1:nt,1:nt) = hxp_temp(1:nt,1:nt)

hsf = 0d0

hsf = matmul(matmul(transpose(lambda),hxp),lambda)

write(*,*) 'hsf'
write(*,fmt2) hsf

deallocate(hxp,hxp_temp,eigval)

end subroutine get_lambda_eigenvectors

subroutine update_a2(c2,x,a2)
implicit none

complex(8),intent(out) :: a2
complex(8),dimension(:),intent(in) :: x

integer :: i,n

real(8),dimension(:),intent(in) :: c2

n = size(c2)

a2 = cmplx(0d0,0d0)
do i = 1, n
   a2 = a2 + 2d0*c2(i)*x(i)
end do
end subroutine update_a2

subroutine update_hmp(ng,nb,nmap,ed,a1,c2,x,hc,hm)
implicit none

complex(8) :: a2
complex(8),intent(in) :: a1
complex(8),dimension(:),intent(in) :: x
complex(8),dimension(:,:),intent(inout) :: hm

integer :: a,b,i,n
integer,intent(in) :: ng,nb,nmap

real(8),intent(in) :: ed
real(8),dimension(:),intent(in) :: c2
real(8),dimension(:,:),intent(in) :: hc

n = size(x)

a2 = cmplx(0d0,0d0)
do i = 1, n
   a2 = a2 + 2d0*c2(i)*x(i)
end do

!only part that is updated
!3 x 3
hm(ng+nb+1:nmap,ng+nb+1:nmap) = 0d0
do a = ng+nb+1, nmap
   do b = ng+nb+1, nmap
      if (a == b) hm(a,b) = hc(a,b) + (ed + a1 + a2)
   end do
end do

end subroutine update_hmp

subroutine update_x(dt,p,x)
implicit none

complex(8),dimension(:),intent(in) :: p
complex(8),dimension(:),intent(inout) :: x

integer :: i,n

real(8) :: dt

n = size(x)

do i = 1, n
   x(i) = x(i) + dt*p(i)
end do

end subroutine update_x

subroutine update_pm(dt,hm,rm,pm)
implicit none

complex(8),dimension(:),intent(in) :: rm
complex(8),dimension(:),intent(inout) :: pm
complex(8),dimension(:,:),intent(in) :: hm

integer :: i,j,n

real(8),intent(in) :: dt

n = size(pm)

do i = 1, n
   do j = 1, n
      pm(i) = pm(i) - dt*hm(i,j)*rm(j)
   end do
end do

end subroutine update_pm

subroutine update_rm(dt,hm,pm,rm)
implicit none

complex(8),dimension(:),intent(in) :: pm
complex(8),dimension(:),intent(inout) :: rm
complex(8),dimension(:,:),intent(in) :: hm

integer :: i,j,n

real(8),intent(in) :: dt 

n = size(rm)

do i = 1, n
   do j = 1, n
      rm(i) = rm(i) + dt*hm(i,j)*pm(j)
   end do
end do
end subroutine update_rm

subroutine update_p(dt,f,p)
implicit none

complex(8),dimension(:),intent(in) :: f
complex(8),dimension(:),intent(inout) :: p

integer :: i,n

real(8),intent(in) :: dt

n = size(p)

do i = 1, n
   p(i) = p(i) + dt*f(i)
end do

end subroutine update_p

subroutine get_hm2(nmap,ng,nb,mu,et,a1,a2,hs,hm)
implicit none

integer :: i
integer,intent(in) :: nmap,ng,nb

complex(8),intent(in) :: et,a1,a2
complex(8),dimension(:,:),intent(out) :: hm

real(8),intent(in) :: mu
real(8),dimension(:,:),intent(in) :: hs

hm = hs
hm(1:ng,ng+1:ng+nb) = hs(1:ng,ng+1:ng+nb)*(-mu*et)
hm(ng+1:ng+nb,1:ng) = hs(ng+1:ng+nb,1:ng)*(-mu*et)
do i = ng+nb+1, nmap
      hm(i,i) = hm(i,i) + (a1+a2)
end do
end subroutine get_hm2

subroutine get_pulsefield(np,tau,it,dt,time,g,e0,e1,omega,et)
implicit none

complex(8) :: etc
complex(8),intent(out) :: et

integer :: i
integer,intent(in) :: np,it
integer,dimension(:),intent(in) :: g

real(8),intent(in) :: dt,e0,e1
real(8),dimension(:),intent(in) :: tau,time,omega

et = cmplx(0d0,0d0)

do i = 1,np
   call get_fieldcomponent(i,tau(i),it,dt,time(i),g(i),E0,E1,omega(i),etc)
   et = et + etc
end do

end subroutine get_pulsefield

subroutine get_fieldcomponent(i,tau,it,dt,time,g,E0,E1,omega,etc)
implicit none

complex(8) :: temp
complex(8),intent(out) :: etc

integer :: it
integer,intent(in) :: i,g

real(8) :: gaussian,arg
real(8),intent(in) :: dt,E0,E1,tau,time,omega

gaussian = get_pulseenvelope(tau,it,dt,time)

arg = omega*((it - 0.5d0)*dt - time)
select case(i)
   case(1)
      temp = E0*cmplx(cos(arg), sin(arg))
   case(2)
      temp = E0*cmplx(cos(arg),-sin(arg))
   case(3)
      temp = E1*cmplx(cos(arg),-sin(arg))
end select

etc = g*gaussian*temp

end subroutine get_fieldcomponent

function get_pulseenvelope(tau,it,dt,time) result(f)
implicit none

integer,intent(in) :: it

real(8) :: f,ln2t4
real(8),intent(in) :: tau,dt,time

ln2t4 = 4d0*dlog(2d0)

f = dsqrt(ln2t4/(pi*tau**2))*dexp(-ln2t4*((it - 0.5d0)*dt - time)**2/tau**2)

end function get_pulseenvelope

subroutine get_force(nmap,ng,nb,lld,kosc,x,c2,rm,pm,f)
implicit none

complex(8),dimension(:),intent(in) :: x,rm,pm
complex(8),dimension(:),intent(out) :: f

integer :: a,b,i,j,n
integer,intent(in) :: nmap,ng,nb

real(8),dimension(:),intent(in) :: kosc,c2
real(8),dimension(:,:),intent(in) :: lld

n = size(x)

f = cmplx(0d0,0d0)
do j = 1, n
   f(j) = -kosc(j)*x(j)
   !exclude lambdas from 1 to ng, because I know those won't contribute
   do a = ng+1, nmap
      do b = ng+1, nmap
         if (a == b) then
            f(j) = f(j) - c2(j)*lld(a,b)*(rm(a)*rm(b) + pm(a)*pm(b) - 1d0)
         else
            f(j) = f(j) - c2(j)*lld(a,b)*(rm(a)*rm(b) + pm(a)*pm(b))
         end if
      end do
   end do
end do

end subroutine get_force

subroutine get_force_traceless(nmap,ng,nb,lld,kosc,x,c2,rm,pm,f,fcla,ftra,fqua)
implicit none

complex(8),dimension(:),allocatable :: c
complex(8),dimension(:),intent(in) :: rm,pm,x
complex(8),dimension(:),intent(out) :: f,fcla,ftra,fqua

integer :: a,b,i,j,n
integer,intent(in) :: nmap,ng,nb

real(8),dimension(:),intent(in) :: kosc,c2
real(8),dimension(:,:),intent(in) :: lld
!real(8),dimension(:,:),allocatable :: mdh

!allocate(dh(1:nmap,1:nmap))
allocate(c(1:nmap))

n = size(c2)

!getting product for faster calculation
c = cmplx(0d0,0d0)
do a = 1, nmap
   c(a) = 0.5d0*(rm(a)**2 + pm(a)**2)
end do

f = cmplx(0d0,0d0)
fcla = cmplx(0d0,0d0)
ftra = cmplx(0d0,0d0)
fqua = cmplx(0d0,0d0)

do j = 1, n
   fcla(j) = -kosc(j)*x(j)
   
!   mdh = (lld)*c2(j)*2d0
   
   ftra(j) = (nmap-ng-nb)*(-2d0*c2(j))/nmap
!   do a = 1, nmap
!      trace = trace + dh(a,a)
!   end do

!   tn = trace/nmap
   !for force trace is substracted, in hamiltonian the trace is added (F = -Div V)
!   f(j) = f(j) - tn
!   do a = 1, nmap
!      f(j) = f(j) - (dh(a,a) - tn)*c(a)
!   end do

!   fmap(j) = f(j) - fclas(j)
   do a = 1, ng+nb
      fqua(j) = fqua(j) + (-ftra(j))*c(a)
   end do
   do a = ng+nb+1,nmap
      fqua(j) = fqua(j) + (-2d0*c2(j)-ftra(j))*c(a)
   end do

   f(j) = fcla(j) + ftra(j) + fqua(j)
end do

end subroutine get_force_traceless

subroutine get_force_total(ng,nb,nmap,kosc,x,c2,rm,pm,f)
implicit none

complex(8),dimension(:),intent(in) :: x,rm,pm
complex(8),dimension(:),intent(out) :: f

integer :: i,n
integer,intent(in) :: ng,nb,nmap

real(8),dimension(:),intent(in) :: kosc,c2

n = size(f)

do i = 1, n
   f(i) = get_force_oneosc(kosc(i),x(i)) + get_force_ssosc(ng,nb,nmap,c2(i),rm,pm)
end do

end subroutine get_force_total

function get_force_ssosc(ng,nb,nmap,c2,rm,pm) result(f)
implicit none

complex(8),dimension(:),intent(in) :: rm,pm

integer :: a
integer,intent(in) :: ng,nb,nmap

real(8),intent(in) :: c2

complex(8) :: f

f = 0d0
do a = ng+nb+1, nmap
   f = f - c2*(rm(a)**2 + pm(a)**2 - 1d0)
end do
end function get_force_ssosc

function get_force_oneosc(k,x) result(f)

complex(8),intent(in) :: x

real(8),intent(in) :: k

complex(8) :: f

f = -k*x

end function get_force_oneosc

subroutine get_a(c2,ome,x,a1,a2)
implicit none

complex(8),dimension(:),intent(in) :: x

complex(8), intent(out) :: a1,a2

real(8),dimension(:),intent(in) :: c2,ome

integer :: i,n

n = size(c2)

a1 = cmplx(0d0,0d0)
a2 = cmplx(0d0,0d0)
do i = 1, n
   a1 = a1 + 2d0*c2(i)**2/ome(i)**2
   a2 = a2 + 2d0*c2(i)*x(i)
end do

end subroutine

subroutine sampling_mapng(init,rm,pm)
implicit none

real(8),parameter :: sqrt2 = 1.414213562d0

complex(8),dimension(:),intent(out) :: rm, pm

integer,intent(in) :: init

integer :: i, n

real(8) :: gauss

n = size(rm)

select case(init)
   case (:2)
      rm = cmplx(0d0,0d0)
      pm = cmplx(0d0,0d0)
   case (3)
      do i = 1, n
         call gauss_noise2(gauss)
         rm(i) = cmplx(gauss/sqrt2,0d0)
         
         call gauss_noise2(gauss)
         pm(i) = cmplx(gauss/sqrt2,0d0)
      end do
   case (4:)
      rm = cmplx(0d0,0d0)
      pm = cmplx(0d0,0d0)
end select

end subroutine sampling_mapng

subroutine sampling_class(bath,beta,kosc,c2,x,p)
implicit none

complex(8),dimension(:),intent(out) :: x,p

integer,intent(in) :: bath

real(8),intent(in) :: beta

real(8),dimension(:),intent(in) :: kosc,c2

integer :: i, nosc

real(8) :: uj,qbeta,gauss

nosc = size(kosc)

do i = 1, nosc
   uj = 0.5d0*beta*sqrt(kosc(i))
   qbeta = beta/(uj/tanh(uj))
   
   call gauss_noise2(gauss)
   p(i) = cmplx(gauss/sqrt(qbeta),0d0)

   call gauss_noise2(gauss)
   x(i) = cmplx(gauss/sqrt(qbeta*kosc(i)),0d0)

   if (bath == 1) x(i) = x(i) + c2(i)/kosc(i)
end do

end subroutine sampling_class

subroutine get_preh(ng,nb,nd,eg,eb,ed,delta,omega,hs)
use m_vib
implicit none

integer,parameter :: ip = 10000

character(len=2) :: c_nt
character(len=9) :: fmt1

integer :: i,j,nt
integer,intent(in) :: ng,nb,nd

real(8) :: cg,cb,cd,uint,lint,alpha
real(8),intent(in) :: eg,eb,ed,delta,omega
real(8),dimension(:,:),intent(out) :: hs

nt = ng + nb + nd

cg = 0d0
cb = 2d0*sqrt(10d0)/omega
cd = cb/2d0

uint = cb + 6d0
lint = cb - 6d0

alpha = sqrt(omega)

hs = 0d0

!fill g|g
do i = 1, ng
   hs(i,i) = eg + (i - 0.5d0)*omega
end do
!fill g|b
do i = 1, ng
   do j = ng+1, ng+nb
      hs(i,j) = integrate_t_phiphi(ip,lint,uint,i,cg,j-ng,cb,alpha)
   end do
end do
!fill g|d
!not necessary
!fill b|g
do i = ng+1, ng+nb
   do j = 1, ng
      hs(i,j) = integrate_t_phiphi(ip,lint,uint,i-ng,cb,j,cg,alpha)
   end do
end do
!fill b|b
do i = ng+1, ng+nb
   hs(i,i) = eb + (i -ng - 0.5d0)*omega
end do
!fill b|d
do i = ng+1, ng+nb
   do j = ng+nb+1, nt
      hs(i,j) = delta*integrate_t_phiphi(ip,lint,uint,i-ng,cb,j-ng-nb,cd,alpha)
   end do
end do
!fill d|g
!not necessary
!fill d|b
do i = ng+nb+1, nt
   do j = ng+1, ng+nb
      hs(i,j) = delta*integrate_t_phiphi(ip,lint,uint,i-ng-nb,cd,j-ng,cb,alpha)
   end do
end do
!fill d|d
do i = ng+nb+1, nt
   hs(i,i) = ed + (i -ng -nb - 0.5d0)*omega
end do

if (nt > 9) then
   write(c_nt,'(i2)') nt
else
   write(c_nt,'(i1)') nt
end if

fmt1 = '('//trim(c_nt)//'f10.5)'

!print fmt1, hs
end subroutine get_preh

subroutine iniconq_d(nosc,lumda_d,ome_max,ome,c2,kosc)
implicit none

integer :: i
integer,intent(in) :: nosc

real(8),intent(in) :: ome_max,lumda_d
real(8),dimension(:),intent(inout) :: ome,c2,kosc

do i = 1, nosc
   ome(i)  = tan(i*atan(ome_max)/nosc)
   c2(i)   = ome(i)*sqrt(atan(ome_max)*lumda_d/(pi*nosc))
   kosc(i) = ome(i)**2
end do

end subroutine iniconq_d


subroutine gauss_noise2(gauss)
implicit none

real(8),parameter :: twopi = 2d0*3.141592654
real(8) :: z1,z2
real(8),intent(out) :: gauss

call random_number(z1)
call random_number(z2)
gauss = sqrt(-2d0*log(z1))*cos(twopi*z2)

end subroutine gauss_noise2

end module m_map
