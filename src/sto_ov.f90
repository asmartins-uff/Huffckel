! -----------------------------------------------------------------------------
!                             MODULE OVERLAP                              
!
! This module contains all subroutines needed to calculate the overlap between
! The basis orbitals (S, P and D) and to build the Huckel Hamiltonian.
!
! LAST MODIFICATION: 20/12/2012
! -----------------------------------------------------------------------------

module sto_ov
implicit none

! Common Variables to the module:

public :: calcovlp, overlap, arrsel, twocenter

contains

! -----------------------------------------------------------------------------
! SUBROUTINE ARRSEL: Load the arrays that control the cases inside the OVERLAP
! Subroutine.
! -----------------------------------------------------------------------------

subroutine arrsel(nnsel,llsel)
implicit none

!----------------------- Variables and Arrays ---------------------------------
integer, parameter :: lmax = 2, nmax = 5
integer            :: llsel(0:lmax,0:lmax), nnsel(nmax,nmax)
! -----------------------------------------------------------------------------

! The NNSEL(n1,n2) array specific the principal quantum numbers involve.

nnsel(1,1) = 1
nnsel(1,2) = 2
nnsel(1,3) = 3
nnsel(1,4) = 4
nnsel(1,5) = 5
nnsel(2,1) = 2
nnsel(2,2) = 6
nnsel(2,3) = 7
nnsel(2,4) = 8
nnsel(2,5) = 9
nnsel(3,1) = 3
nnsel(3,2) = 7
nnsel(3,3) = 10
nnsel(3,4) = 11
nnsel(3,5) = 12
nnsel(4,1) = 4
nnsel(4,2) = 8
nnsel(4,3) = 11
nnsel(4,4) = 13
nnsel(4,5) = 14
nnsel(5,1) = 5
nnsel(5,2) = 9
nnsel(5,3) = 12
nnsel(5,4) = 14
nnsel(5,5) = 15

! Cases for the overlap between the orbitals:
! llsel = 1 -> <s1|s2>
! llsel = 2 -> <s1|p2>
! llsel = 3 -> <s1|d2>
! llsel = 4 -> <p1|p2>
! llsel = 5 -> <p1|d2>
! llsel = 6 -> <d1|d2>

llsel(0,0) = 0
llsel(0,1) = 1
llsel(0,2) = 2
llsel(1,0) = 3
llsel(1,1) = 4
llsel(1,2) = 5
llsel(2,0) = 6
llsel(2,1) = 7
llsel(2,2) = 8

end subroutine arrsel

! -----------------------------------------------------------------------------
! Compute the overlap matrix among the orbitals of the atoms. In the symmetric
! Case, that is, for systems where the atoms are arranjed in a ordered array,
! The overlap matrix includes only the atoms that belong to the basis of the
! Unit cell.
! -----------------------------------------------------------------------------

subroutine calcovlp(i,j,k,n,nat,numtp,nzt,nval,lval,itype,zeta,cf,ui,uj,&
                    llsel,nnsel,ovout,ipbc)

!----------------------- Variables and Arrays ---------------------------------
! (I) Extern and parameters:
integer, parameter  :: nmax = 5, lmax = 2
integer, intent(in) :: i, j, k, n, nat, numtp, ipbc
integer, intent(in) :: llsel(0:lmax,0:lmax), nnsel(nmax,nmax)
integer, intent(in) :: nzt(numtp,3), nval(numtp,3), lval(numtp,3), itype(nat)
real*8,  intent(in) :: zeta(numtp,3,2), cf(numtp,3,2), ui(3), uj(3)
! (II) Intern:
real*8, intent(out) :: ovout(25)
real*8              :: rho, ovlp(0:lmax), ovin(3), zti, ztj, fc
real*8              :: cl, cm, cn, cl2, cm2, cn2, xij, yij, zij
integer             :: ni, nj, li, lj, iz, jz, llopt, nnopt, kj, lp
! (III) Commons:
real*8              :: sdx, sdy, sdz
! -----------------------------------------------------------------------------

common /sdx / sdx
common /sdy / sdy
common /sdz / sdz

! Initializing OVOUT:

ovout(:) = 0.0d0

! Convertion factor for the distances:

fc = 0.529177d0

! Assigning the ni and li quantum numbers of (I,K) orbital:

ni  = nval(itype(i),k)
li  = lval(itype(i),k)

! Assigning the nj and lj quantum numbers of (J,N) orbital:

nj  = nval(itype(j),n)
lj  = lval(itype(j),n)

! Distance |rj - ri| with PBC (in Bohrs):

  xij = (uj(1) - ui(1)) - dnint((uj(1) - ui(1))/sdx)*sdx*ipbc
  yij = (uj(2) - ui(2)) - dnint((uj(2) - ui(2))/sdy)*sdy*ipbc
  zij = (uj(3) - ui(3)) - dnint((uj(3) - ui(3))/sdz)*sdz*ipbc
  rho = dsqrt(xij**2 + yij**2 + zij**2)/fc

! Loop over the zetas:

ovin(:) = 0.0d0

do iz = 1, nzt(itype(i),k)
  do jz = 1, nzt(itype(j),n)

  zti = zeta(itype(i),k,iz)
  ztj = zeta(itype(j),n,jz)

! Overlap <ni,li|nj,lj>:

  call overlap(ni,nj,li,lj,nnsel,llsel,zti,ztj,rho,ovlp)

  ovin(1) = ovin(1) + cf(itype(i),k,iz)*cf(itype(j),n,jz)*ovlp(0)
  ovin(2) = ovin(2) + cf(itype(i),k,iz)*cf(itype(j),n,jz)*ovlp(1)
  ovin(3) = ovin(3) + cf(itype(i),k,iz)*cf(itype(j),n,jz)*ovlp(2)

  enddo
enddo

! Applying the Two-Center Slater and Koster Rules (SIGNALS TO BE VERYFIED):
! Direction cosines:
!
!   cl  : (xj-xi)/|rj-ri|
!   cm  : (yj-yi)/|rj-ri|
!   cn  : (zj-zi)/|rj-ri|
!   cl2,cm2,cn2: squares of the direction cosines cl,cm,cn

cl = xij/(rho*fc)
cm = yij/(rho*fc)
cn = zij/(rho*fc)
cl2 = cl*cl
cm2 = cm*cm
cn2 = cn*cn
   
call twocenter(li,lj,ovin,ovout,cl,cl2,cm,cm2,cn,cn2)

! End of the subroutine

end subroutine calcovlp

! -----------------------------------------------------------------------------
! Compute the overlap integral between two Slater-Type orbitals (STO), denoted
! by O_q(kp,t). The STO are normalized, with quantum numbers (na, la, m) and
! (nb, lb, m) and screening constants (ka, kb). We have: kp = ka/kb, t = kb*rho 
! And q = (na,nb,la,lb,m) the quantum numbers. The overlaps were tabuled by 
! Michael Barnett and obey the following convention:
! n1, n2: the principal quantum numbers (n1 >= n2)
! l1, l2: azimuthal quantum numbers
!      m: magnetic quantum number
!    rho: distance between the a and b sites (in BOHRS).
!  alpha: (alpha)
! alphan: alpha**n (n=2,3,4,5,...)
! -----------------------------------------------------------------------------

subroutine overlap(n1,n2,l1,l2,nnsel,llsel,zeta1,zeta2,rho,ovlp)
implicit none

!----------------------- Variables and Arrays ---------------------------------
! (I) Extern:
real(8)             :: zeta1, zeta2, rho
integer             :: n1, n2, l1, l2
integer, intent(in) :: llsel(0:2,0:2), nnsel(5,5)
! (II) Intern:
real(8) :: ovlp(0:2), pi, k, t, f, skt, k2, k3, k4, k5, k6, k7, k8, k9, t2, t3, t4, t5, t6, t7, t8, t9, t10, sqrtk
real(8) :: alpha, alpha2, alpha3, alpha4, alpha5, alpha6, exp_minus_t
integer :: na, nb, la, lb, llopt, nnopt
! -----------------------------------------------------------------------------

! Initializing the overlaps and other parameters:

ovlp(:) = 0
la = l1
lb = l2
na = n1
nb = n2
skt = dfloat(n1)/dfloat(n2)

! The overlap tables impose n1 >= n2. When n1 < n2, we have to interchange the
! Orbitals in order to get the overlap: in this case, we have to redefine K
! And T:

k = zeta1/zeta2
t = zeta2*rho
f = 1.0

! Redefining K and T in order to impose n1 > n2 or, when n1 = n2, to correct
! The overlap signal when l1 < l2:

if(skt.lt.1) then
t = k*t
k = 1./k
f = (-1.)**(l1 + l2)
la = l2
lb = l1
endif

if((skt.eq.1).and.(l1.lt.l2)) then
t = k*t
k = 1./k
f = (-1.)**(l1 + l2)
endif

! Selectioning llopt and nnopt:

llopt = llsel(la,lb)
nnopt = nnsel(n1,n2)

k2 = k*k
k3 = k2*k
k4 = k3*k
k5 = k4*k
k6 = k5*k 
k7 = k6*k 
k8 = k7*k 
k9 = k8*k 

t2 = t*t
t3 = t2*t
t4 = t3*t
t5 = t4*t
t6 = t5*t
t7 = t6*t
t8 = t7*t
t9 = t8*t
t10 = t9*t

alpha = (alpha)
alpha2 = alpha*alpha
alpha3 = alpha2*alpha
alpha4 = alpha3*alpha
alpha5 = alpha4*alpha
alpha6 = alpha5*alpha

sqrtk = sqrtk
exp_minus_t = exp_minus_t

select case(nnopt)

! ------------------------------ n1 = n2 = 1 -----------------------------------
  case(1) 
! ------------------------------------------------------------------------------

select case(llopt)

!*********************************
case(0)   ! l1 = l2 = 0 => <s1|s2>
!*********************************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((3.+ 3.*t + t2)*dexp_minus_t)/3.
else
ovlp(0) = ovlp(0) + (8*k*((-4*k + k*(alpha)*t)*exp_minus_t + &
                    (4*k + (alpha)*t)*dexp(-(k*t)))*sqrtk)/(alpha3*t)
endif

end select

! ------------------------------ n1 = 2 and n2 = 1 -----------------------------
  case(2) 
! ------------------------------------------------------------------------------

select case(llopt)

!*********************************
case(0)   ! l1 = l2 = 0 => <s1|s2>
!*********************************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((9 + 9*t + 4*t2 + t3)*exp_minus_t)/(6*sqrt(3.))
else
ovlp(0) = ovlp(0) + (8*k2*((-4 - 20*k2 + (alpha)*(1 + 3*k2)*t)*exp_minus_t + (4 + 20*k2 &
                  + 8*k*(alpha)*t + alpha2*t2)*exp(-(k*t)))*sqrtk)/(sqrt(3.)*alpha4*t)
endif

!*******************
case(3)   ! <p1|s2>
!********************

if(k.eq.1) then
ovlp(0) = ovlp(0) + (t*(3 + 3*t + t2)*exp_minus_t)/6
else
ovlp(0) = ovlp(0) + (8*k2*((-24*k - 24*k*t + 4*k*(alpha)*t2)*exp_minus_t + (24*k + 24*k2*t &
                  + 8*k*(alpha)*t2 + alpha2*t3)*exp(-(k*t)))*sqrtk)/(alpha4*t2)
endif

end select

! ------------------------------ n1 = 3 and n2 = 1 -----------------------------
  case(3) 
! ------------------------------------------------------------------------------

select case(llopt)

!*********************************
case(0)   ! l1 = l2 = 0 => <s1|s2>
!*********************************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((30 + 30*t + 15*t2 + 5*t3 + t4)*exp_minus_t*sqrt(2./5.))/30
else
ovlp(0) = ovlp(0) + (8*k3*((-72*k - 120*k3 + 12*k*(alpha)*(1 + k2)*t)*exp_minus_t &
                  + (72*k + 120*k3 + 12*(alpha)*(1 + 5*k2)*t + 12*k*alpha2*t2 &
                  + alpha3*t3)*exp(-(k*t)))*sqrt(2./5.)*sqrtk)/(3*alpha5*t)
endif

!*********************
case(3)   ! <p1|s2>
!*********************

if(k.eq.1) then
ovlp(0) = ovlp(0) + (t*(30 + 30*t + 13*t2 + 3*t3)*exp_minus_t*sqrt(2./15.))/30
else
ovlp(0) = ovlp(0) + (8*k3*((-24 - 168*k2 - 24*(1 + 7*k2)*t + 4*(alpha)*(1 + 5*k2)*t2)*exp_minus_t &
                  + (24 + 168*k2 + 24*k*(1 + 7*k2)*t + 8*(alpha)*(1 + 8*k2)*t2 + 12*k*(-1 &
                  + k2)**2*t3 + alpha3*t4)*exp(-(k*t)))*sqrt(2./15.)*sqrtk)/(alpha5*t2)
endif

!*********************
case(6)   ! <d1|s2>
!**********************

if(k.eq.1) then
ovlp(0) = ovlp(0) + (t2*(3 + 3*t + t2)*exp_minus_t*sqrt(2.))/30
else
ovlp(0) = ovlp(0) + (8*k3*((-576*k - 576*k*t + 24*(-3 + k)*k*(3 + k)*t2 + 24*k*(alpha)*t3)*exp_minus_t &
                  + (576*k + 576*k2*t + 24*k*(-3 + 11*k2)*t2 + 72*k2*(alpha)*t3 + 12*k*(-1 &
                  + k2)**2*t4 + alpha3*t5)*exp(-(k*t)))*sqrt(2.)*sqrtk)/(3*alpha5*t3)
endif

end select

! ------------------------------ n1 = 4 and n2 = 1 -----------------------------
  case(4) 
! ------------------------------------------------------------------------------

select case(llopt)

!*********************************
case(0)   ! l1 = l2 = 0 => <s1|s2>
!*********************************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((225 + 225*t + 120*t2 + 45*t3 + 12*t4 + 2*t5)*exp_minus_t)/(90*sqrt(35.))
else
ovlp(0) = ovlp(0) + (8*k4*((-72 - 1008*k2 - 840*k4 + 12*(alpha)*(1 + 10*k2 + 5*k4)*t)*exp_minus_t + &
                    (72 + 1008*k2 + 840*k4 + 96*k*(alpha)*(3 + 5*k2)*t + 24*alpha2*(1 + 5*k2)*t2 &
                  + 16*k*alpha3*t3 + alpha4*t4)*exp(-(k*t)))*sqrtk)/(3*sqrt(35.)*alpha6*t)
endif

!*********************
case(3)   ! <p1|s2>
!*********************

if(k.eq.1) then
ovlp(0) = ovlp(0) + (t*(75 + 75*t + 36*t2 + 11*t3 + 2*t4)*exp_minus_t)/(30*sqrt(105.))
else
ovlp(0) = ovlp(0) + (8*k4*((-576*k - 1344*k3 - 192*k*(3 + 7*k2)*t + 24*k*(alpha)*(3 + 5*k2)*t2)*exp_minus_t &
                  + (576*k + 1344*k3 + 192*k2*(3 + 7*k2)*t + 24*k*(alpha)*(9 + 23*k2)*t2 + &
                  4*alpha2*(5 + 31*k2)*t3 + 16*k*alpha3*t4 + &
                  alpha4*t5)*exp(-(k*t)))*sqrtk)/(sqrt(105.)*alpha6*t2)
endif


!*********************
case(6)   ! <d1|s2>
!**********************

if(k.eq.1) then
ovlp(0) = ovlp(0) + (t2*(21 + 21*t + 9*t2 + 2*t3)*exp_minus_t)/(90*sqrt(7.))
else
ovlp(0) = ovlp(0) + (8*k4*((-576 - 5184*k2 - 576*(1 + 9*k2)*t + 24*(-9 - 78*k2 + 7*k4)*t2 + &
                 24*(alpha)*(1 + 7*k2)*t3)*exp_minus_t + (576 + 5184*k2 + 576*k*(1 + 9*k2)*t + &
                 24*(-3 - 18*k2 + 101*k4)*t2 + 24*k*(alpha)*(3 + 29*k2)*t3 + 12*alpha2*(1 & 
               + 11*k2)*t4 + 16*k*alpha3*t5 + alpha4*t6)*exp(-(k*t)))*sqrtk)/(3*sqrt(7.)*(-1 &
               + k2)**6*t3)
endif

end select

! ------------------------------ n1 = 5 and n2 = 1 -----------------------------
  case(5) 
! ------------------------------------------------------------------------------

select case(llopt)

!*********************************
case(0)   ! l1 = l2 = 0 => <s1|s2>
!*********************************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((945 + 945*t + 525*t2 + 210*t3 + 63*t4 + 14*t5 + 2*t6)*exp_minus_t*sqrt(2./7.))/1890.
else
ovlp(0) = ovlp(0) + (8*k5*((-2880*k - 13440*k3 - 6720*k5 + 120*k*(alpha)*(3 + k2)*(1 + 3*k2)*t)*exp_minus_t & 
                  + (2880*k + 13440*k3 + 6720*k5 + 120*(alpha)*(3 + 42*k2 + 35*k4)*t + 240*k*(-1 &
                  + k2)**2*(3 + 5*k2)*t2 + 40*alpha3*(1 + 5*k2)*t3 + 20*k*alpha4*t4  &
                  + alpha5*t5)*exp(-(k*t)))*sqrt(2./7.)*sqrtk)/(45*alpha7*t)
endif

!*********************
case(3)   ! <p1|s2>
!*********************

if(k.eq.1) then
ovlp(0) = ovlp(0) + (t*(1575 + 1575*t + 798*t2 + 273*t3 + 66*t4 + 10*t5)*exp_minus_t*sqrt(2./21.))/3150. 
else
ovlp(0) = ovlp(0) + (8*k5*((-576 - 10368*k2 - 12096*k4 - 576*(1 + 18*k2 + 21*k4)*t + 24*(-1 &
                  + k2)*(3 + 42*k2 + 35*k4)*t2)*exp_minus_t + (576 + 10368*k2 + 12096*k4 &
                  + 576*k*(1 + 18*k2 + 21*k4)*t + 24*(alpha)*(9 + 174*k2 + 217*k4)*t2 + 48*k*(-1 &
                  + k2)**2*(13 + 27*k2)*t3 + 12*alpha3*(3 + 17*k2)*t4 + 20*k*alpha4*t5 &
                  + alpha5*t6)*exp(-(k*t)))*sqrt(2./21.)*sqrtk)/(15*alpha7*t2)
endif

!*********************
case(6)   ! <d1|s2>
!**********************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((210*t3 + 210*t4 + 99*t5 + 29*t6 + 5*t7)*exp_minus_t*sqrt(2./35.))/(945*t)
else
ovlp(0) = ovlp(0) + (8*k5*((-17280*k - 51840*k3 - 17280*k*(1 + 3*k2)*t + 192*k*(-33 - 94*k2 + 7*k4)*t2 &
                  + 192*k*(alpha)*(3 + 7*k2)*t3)*exp_minus_t + (17280*k + 51840*k3 + 17280*k2*(1 + 3*k2)*t &
                  + 768*k*(-3 + k2 + 32*k4)*t2 + 384*k2*(alpha)*(6 + 19*k2)*t3 + 48*k*(-1 &
                  + k2)**2*(9 + 31*k2)*t4 + 4*alpha3*(7 + 53*k2)*t5 + 20*k*alpha4*t6 &
                  + alpha5*t7)*exp(-(k*t)))*sqrt(2./35.)*sqrtk)/(9*alpha7*t3)
endif


end select

! ------------------------------ n1 = n2 = 2 -----------------------------------
  case(6) 
! ------------------------------------------------------------------------------

select case(llopt)

!*********************************
case(0)   ! l1 = l2 = 0 => <s1|s2>
!*********************************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((45 + 45*t + 20*t2 + 5*t3 + t4)*exp_minus_t)/45
else
ovlp(0) = ovlp(0) + (8*k2*((20 + 152*k2 + 20*k4 - 8*(alpha)*(1 + 5*k2)*t + alpha2*(1 &
                  + 3*k2)*t2)*exp_minus_t + (-20 - 152*k2 - 20*k4 - 8*k*(alpha)*(5 + k2)*t &
                  - alpha2*(3 + k2)*t2)*exp(-(k*t)))*sqrtk)/(3*alpha5*t)
endif

!*********************************
case(1,3)   ! <s1|p2> or <p1|s2>
!*********************************

if(k.eq.1) then
ovlp(0) = ovlp(0) +(t*(15 + 15*t + 7*t2 + 2*t3)*exp_minus_t)/(30*sqrt(3.))
else
ovlp(0) = ovlp(0) + (-8*k2*((-168*k - 24*k3 - 24*k*(7 + k2)*t + 4*k*(alpha)*(11 + k2)*t2 - 4*k*(-1 &
                  + k2)**2*t3)*exp_minus_t + (168*k + 24*k3 + 24*k2*(7 + k2)*t + 8*k*(alpha)*(5 &
                  + k2)*t2 + alpha2*(3 + k2)*t3)*exp(-(k*t)))*sqrtk)/(sqrt(3.)*alpha5*t2)
endif

!*********************************
case(4)   ! l1 = l2 = 1 => <p1|p2>
!*********************************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((15*t + 15*t2 + 3*t3 - 2*t4 - t5)*exp_minus_t)/(15*t)
ovlp(1) = ovlp(1) - ((-15*t - 15*t2 - 6*t3 - t4)*exp_minus_t)/(15*t)
else
ovlp(0) = ovlp(0) + (-32*k2*((96*k + 96*k*t + 48*k*t2 - k*(alpha)*(11 + k2)*t3 + k*(-1 &
                  + k2)**2*t4)*exp_minus_t + (-96*k - 96*k2*t - 48*k3*t2 - (alpha)*(1 + 11*k2)*t3 &
                  - k*alpha2*t4)*exp(-(k*t)))*sqrtk)/(alpha5*t3)
ovlp(1) = ovlp(1) + (32*k2*((48*k + 48*k*t - 12*k*(alpha)*t2 + k*alpha2*t3)*exp_minus_t &
                  + (-48*k - 48*k2*t - 12*k*(alpha)*t2 &
                  - alpha2*t3)*exp(-(k*t)))*sqrtk)/(alpha5*t3)
endif

end select

! ------------------------------ n1 = 3 and n2 = 2 -----------------------------
  case(7) 
! ------------------------------------------------------------------------------

select case(llopt)

!*********************************
case(0)   ! l1 = l2 = 0 => <s1|s2>
!*********************************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((225 + 225*t + 105*t2 + 30*t3 + 6*t4 + t5)*exp_minus_t*sqrt(2./15.))/90
else
ovlp(0) = ovlp(0) + (-8*k3*((-504*k - 1296*k3 - 120*k5 + 48*k*(alpha)*(3 + 5*k2)*t &
                  - 12*k*alpha2*(1 + k2)*t2)*exp_minus_t + (504*k + 1296*k3 + 120*k5 &
                  + 12*(alpha)*(5 + 38*k2 + 5*k4)*t + 12*k*alpha2*(5 + k2)*t2 &
                  + alpha3*(3 + k2)*t3)*exp(-(k*t)))*sqrt(2./15.)*sqrtk)/(3*alpha6*t)
endif

!*******************
case(1)   ! <s1|p2>
!*******************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((-15*t2 - 15*t3 - 9*t4 - 4*t5 - t6)*exp_minus_t*sqrt(2./5.))/(90*t)
else
ovlp(0) = ovlp(0) + (32*k3*((-144*k - 336*k3 - 48*k*(3 + 7*k2)*t + 12*k*(alpha)*(3 + 5*k2)*t2 &
                  - 3*k*alpha2*(1 + k2)*t3)*exp_minus_t + (144*k + 336*k3 + 48*k2*(3 + 7*k2)*t &
                  + 36*k*(alpha)*(1 + 3*k2)*t2 + 2*alpha2*(1 + 8*k2)*t3 &
                  + k*alpha3*t4)*exp(-(k*t)))*sqrt(2./5.)*sqrtk)/(3*alpha6*t2)
endif

!*********************
case(3)   ! <p1|s2>
!*********************

if(k.eq.1) then
ovlp(0) = ovlp(0) + (t*(45 + 45*t + 20*t2 + 5*t3 + t4)*exp_minus_t*sqrt(2./5.))/90
else
ovlp(0) = ovlp(0) + (-8*k3*((-168 - 1584*k2 - 168*k4 - 24*(7 + 66*k2 + 7*k4)*t &
                  + 4*(alpha)*(11 + 80*k2 + 5*k4)*t2 - 4*alpha2*(1 + 5*k2)*t3)*exp_minus_t &
                  + (168 + 1584*k2 + 168*k4 + 24*k*(7 + 66*k2 + 7*k4)*t + 8*(alpha)*(5 + 59*k2 &
                  + 8*k4)*t2 + 12*k*alpha2*(5 + k2)*t3 + alpha3*(3 &
                  + k2)*t4)*exp(-(k*t)))*sqrt(2./5.)*sqrtk)/(3*alpha6*t2)
endif

!*********************************
case(4)   ! l1 = l2 = 1 => <p1|p2>
!*********************************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((75*t + 75*t2 + 24*t3 - t4 - 3*t5 - t6)*exp_minus_t*sqrt(2./15.))/(30*t)
ovlp(1) = ovlp(1) - ((-75*t - 75*t2 - 33*t3 - 8*t4 - t5)*exp_minus_t*sqrt(2./15.))/(30*t)
else
ovlp(0) = ovlp(0) + (32*k3*((-96 - 864*k2 - 96*(1 + 9*k2)*t - 48*(1 + 9*k2)*t2 &
                  + (alpha)*(11 + 80*k2 + 5*k4)*t3 - alpha2*(1 + 5*k2)*t4)*exp_minus_t &
                  + (96 + 864*k2 + 96*k*(1 + 9*k2)*t + 48*k2*(1 + 9*k2)*t2 + 6*k*(-1 &
                  + k2)*(5 + 19*k2)*t3 + 2*alpha2*(1 + 8*k2)*t4 + k*(-1 &
                  + k2)**3*t5)*exp(-(k*t)))*sqrt(2./15.)*sqrtk)/(alpha6*t3)

ovlp(1) = ovlp(1) + (-32*k3*((-48 - 432*k2 - 48*(1 + 9*k2)*t + 12*(alpha)*(1 + 7*k2)*t2 &
                  - alpha2*(1 + 5*k2)*t3)*exp_minus_t + (48 + 432*k2 + 48*k*(1 + 9*k2)*t &
                  + 12*(alpha)*(1 + 11*k2)*t2 + 18*k*alpha2*t3 + (-1 &
                  + k2)**3*t4)*exp(-(k*t)))*sqrt(2./15.)*sqrtk)/(alpha6*t3)

endif

!********************
case(6)   ! <d1|s2>
!********************

if(k.eq.1) then
ovlp(0) = ovlp(0) + (t2*(6 + 6*t + 3*t2 + t3)*exp_minus_t*sqrt(2./3.))/90
else
ovlp(0) = ovlp(0) + (-8*k3*((-5184*k - 576*k3 - 576*k*(9 + k2)*t + 24*k*(-87 + 6*k2 &
                  + k4)*t2 + 24*k*(alpha)*(15 + k2)*t3 - 24*k*alpha2*t4)*exp_minus_t &
                  + (5184*k + 576*k3 + 576*k2*(9 + k2)*t + 24*k*(-21 + 90*k2 + 11*k4)*t2 &
                  + 72*k2*(alpha)*(7 + k2)*t3 + 12*k*alpha2*(5 + k2)*t4 &
                  + alpha3*(3 + k2)*t5)*exp(-(k*t)))*sqrt(2./3.)*sqrtk)/(3*alpha6*t3)

endif

!*********************
case(7)  ! <d1|p2> 
!*********************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((30*t2 + 30*t3 + 9*t4 - t5 - t6)*exp_minus_t*sqrt(2.))/(90*t)
ovlp(1) = ovlp(1) + ((15*t2 + 15*t3 + 6*t4 + t5)*exp_minus_t*sqrt(2./3.))/(30*t)
else
ovlp(0) = ovlp(0) + (32*k3*((-4320*k - 4320*k*t + 96*k*(-21 + k2)*t2 + 96*k*(-6 + k2)*t3 &
                  + 6*k*(alpha)*(15 + k2)*t4 - 6*k*alpha2*t5)*exp_minus_t + (4320*k  &
                  + 4320*k2*t + 48*k*(-3 + 43*k2)*t2 + 48*k2*(-3 + 13*k2)*t3 + 18*k*(-1 &
                  + k2)*(1 + 7*k2)*t4 + 2*alpha2*(1 + 8*k2)*t5 + k*(-1 &
                  + k2)**3*t6)*exp(-(k*t)))*sqrt(2.)*sqrtk)/(3*alpha6*t4)

ovlp(1) = ovlp(1) + (-32*k3*((-1440*k - 1440*k*t + 96*k*(-6 + k2)*t2 + 96*k*(alpha)*t3 &
                  - 6*k*alpha2*t4)*exp_minus_t + (1440*k + 1440*k2*t + 48*k*(-3 + 13*k2)*t2 &
                  + 144*k2*(alpha)*t3 + 18*k*alpha2*t4 + (-1 &
                  + k2)**3*t5)*exp(-(k*t)))*sqrt(2./3.)*sqrtk)/(alpha6*t4)

endif

end select

! ------------------------------ n1 = 4 and n2 = 2 -----------------------------
  case(8) 
! ------------------------------------------------------------------------------

select case(llopt)

!*********************************
case(0)   ! l1 = l2 = 0 => <s1|s2>
!*********************************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((4725 + 4725*t + 2310*t2 + 735*t3 + 168*t4 + 28*t5 + 4*t6)*exp_minus_t)/&
                    (630*sqrt(105.))
else
ovlp(0) = ovlp(0) + (-8*k4*((-504 - 9432*k2 - 12264*k4 - 840*k6 + 48*(alpha)*(3 + 42*k2 &
                  + 35*k4)*t - 12*alpha2*(1 + 10*k2 + 5*k4)*t2)*exp_minus_t + (504 + 9432*k2 &
                  + 12264*k4 + 840*k6 + 96*k*(alpha)*(21 + 54*k2 + 5*k4)*t + 24*alpha2*(5 &
                  + 38*k2 + 5*k4)*t2 + 16*k*alpha3*(5 + k2)*t3 + alpha4*(3 + &
                    k2)*t4)*exp(-(k*t)))*sqrtk)/(3*sqrt(105.)*alpha7*t)
endif

!*******************
case(1)   ! <s1|p2>
!*******************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((-21*t4 - 21*t5 - 9*t6 - 2*t7)*exp_minus_t)/(315*t*sqrt(35.))
else
ovlp(0) = ovlp(0) + (32*k4*((-144 - 2592*k2 - 3024*k4 - 144*(1 + 18*k2 + 21*k4)*t + 12*(-1 &
                  + k2)*(3 + 42*k2 + 35*k4)*t2 - 3*alpha2*(1 + 10*k2 + 5*k4)*t3)*exp_minus_t &
                  + (144 + 2592*k2 + 3024*k4 + 144*k*(1 + 18*k2 + 21*k4)*t + 12*(alpha)*(3 + 66*k2 &
                  + 91*k4)*t2 + 12*k*alpha2*(7 + 17*k2)*t3 + 3*alpha3*(1 + 7*k2)*t4 &
                  + k*alpha4*t5)*exp(-(k*t)))*sqrtk)/(3*Sqrt(35.)*alpha7*t2)
endif

!*********************
case(3)   ! <p1|s2>
!*********************

if(k.eq.1) then
ovlp(0) = ovlp(0) + (t*(1050 + 1050*t + 483*t2 + 133*t3 + 25*t4 + 4*t5)*exp_minus_t)/(630*sqrt(35.))
else
ovlp(0) = ovlp(0) + (-8*k4*((-5184*k - 16512*k3 - 1344*k5 - 192*k*(27 + 86*k2 + 7*k4)*t &
                  + 120*k*(alpha)*(9 + 22*k2 + k4)*t2 - 24*k*alpha2*(3 + 5*k2)*t3)*exp_minus_t &
                  + (5184*k + 16512*k3 + 1344*k5 + 192*k2*(27 + 86*k2 + 7*k4)*t + 24*k*(alpha)*(63 &
                  + 234*k2 + 23*k4)*t2 + 4*alpha2*(25 + 232*k2 + 31*k4)*t3 + 16*k*(-1 &
                  + k2)**3*(5 + k2)*t4 + alpha4*(3 + k2)*t5)*exp(-(k*t)))*sqrtk)/(3*sqrt(35.)*(-1 &
                  + k2)**7*t2)
endif

!*********************************
case(4)   ! l1 = l2 = 1 => <p1|p2>
!*********************************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((1575*t + 1575*t2 + 630*t3 + 105*t4 - 15*t5 - 15*t6 - 4*t7)*exp_minus_t)/&
                    (210*t*sqrt(105.))
ovlp(1) = ovlp(1) - ((-1575*t - 1575*t2 - 735*t3 - 210*t4 - 39*t5 - 4*t6)*exp_minus_t)/(210*t*sqrt(105.))
else
ovlp(0) = ovlp(0) + (32*k4*((-2880*k - 8640*k3 - 2880*k*(1 + 3*k2)*t - 1440*k*(1 + 3*k2)*t2 &
                  +  30*k*(alpha)*(9 + 22*k2 + k4)*t3 - 6*k*alpha2*(3 + 5*k2)*t4)*exp_minus_t &
                  + (2880*k + 8640*k3 + 2880*k2*(1 + 3*k2)*t + 1440*k3*(1 + 3*k2)*t2 + 30*(-1 &
                  + k2)*(1 + 22*k2 + 41*k4)*t3 + 6*k*alpha2*(13 + 35*k2)*t4 + 3*(-1 &
                  + k2)**3*(1 + 7*k2)*t5 + k*alpha4*t6)*exp(-(k*t)))*sqrtk)/(sqrt(105.)*(-1 &
                  + k2)**7*t3)

ovlp(1) = ovlp(1) + (-32*k4*((-1440*k - 4320*k3 - 1440*k*(1 + 3*k2)*t + 96*k*(alpha)*(3 + 7*k2)*t2 &
                  - 6*k*alpha2*(3 + 5*k2)*t3)*exp_minus_t + (1440*k + 4320*k3 + 1440*k2*(1 + 3*k2)*t &
                  + 48*k*(alpha)*(9 + 31*k2)*t2 + 6*alpha2*(5 + 43*k2)*t3 + 24*k*(-1 &
                  + k2)**3*t4 + alpha4*t5)*exp(-(k*t)))*sqrtk)/(sqrt(105.)*alpha7*t3)

endif

!********************
case(6)   ! <d1|s2>
!********************

if(k.eq.1) then
ovlp(0) = ovlp(0) + (t2*(168 + 168*t + 75*t2 + 19*t3 + 4*t4)*exp_minus_t)/(630*sqrt(21.))
else
ovlp(0) = ovlp(0) + (-8*k4*((-5184 - 58752*k2 - 5184*k4 - 1728*(3 + 34*k2 + 3*k4)*t &
                  + 24*(-87 - 939*k2 + 59*k4 + 7*k6)*t2 + 24*(alpha)*(15 + 138*k2 &
                  + 7*k4)*t3 - 24*alpha2*(1 + 7*k2)*t4)*exp_minus_t + (5184 + 58752*k2 &
                  + 5184*k4 + 1728*k*(3 + 34*k2 + 3*k4)*t + 24*(-21 - 177*k2 + 1057*k4 + 101*k6)*t2 &
                  + 24*k*(alpha)*(21 + 270*k2 + 29*k4)*t3 + 12*alpha2*(5 + 80*k2 &
                  + 11*k4)*t4 + 16*k*alpha3*(5 + k2)*t5 + alpha4*(3 &
                  + k2)*t6)*exp(-(k*t)))*sqrtk)/(3*Sqrt(21.)*alpha7*t3)

endif

!*********************
case(7)  ! <d1|p2> 
!*********************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((630*t2 + 630*t3 + 237*t4 + 27*t5 - 9*t6 - 4*t7)*exp_minus_t)/(630*t*sqrt(7.))
ovlp(1) = ovlp(1) + ((315*t2 + 315*t3 + 138*t4 + 33*t5 + 4*t6)*exp_minus_t)/(210*t*sqrt(21.))
else
ovlp(0) = ovlp(0) + (32*k4*((-4320 - 47520*k2 - 4320*(1 + 11*k2)*t + 288*(-7 - 76*k2 + 3*k4)*t2 &
                  + 288*(-2 - 21*k2 + 3*k4)*t3 + 6*(alpha)*(15 + 138*k2 + 7*k4)*t4 - 6*(-1 &
                  + k2)**2*(1 + 7*k2)*t5)*exp_minus_t + (4320 + 47520*k2 + 4320*k*(1 + 11*k2)*t + 144*(-1 &
                  + 2*k2 + 159*k4)*t2 + 144*k*(-1 - 8*k2 + 49*k4)*t3 + 6*(alpha)*(3 + 66*k2 &
                  + 251*k4)*t4 + 6*k*alpha2*(11 + 37*k2)*t5 + 3*alpha3*(1 + 7*k2)*t6 &
                  + k*alpha4*t7)*exp(-(k*t)))*sqrtk)/(3*sqrt(7.)*alpha7*t4)
 
ovlp(1) = ovlp(1) + (-32*k4*((-1440 - 15840*k2 - 1440*(1 + 11*k2)*t + 288*(-2 - 21*k2 + 3*k4)*t2 &
                  + 96*(alpha)*(1 + 9*k2)*t3 - 6*alpha2*(1 + 7*k2)*t4)*exp_minus_t + (1440 &
                  + 15840*k2 + 1440*k*(1 + 11*k2)*t + 144*(-1 - 8*k2 + 49*k4)*t2 + 48*k*(alpha)*(3 &
                  + 37*k2)*t3 + 18*alpha2*(1 + 15*k2)*t4 + 24*k*alpha3*t5 + (-1 &
                  + k2)**4*t6)*exp(-(k*t)))*sqrtk)/(sqrt(21.)*alpha7*t4)

endif

end select

! ------------------------------ n1 = 5 and n2 = 2 -----------------------------
  case(9) 
! ------------------------------------------------------------------------------

select case(llopt)

!*********************************
case(0)   ! l1 = l2 = 0 => <s1|s2>
!*********************************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((6615 + 6615*t + 3360*t2 + 1155*t3 + 294*t4 + 56*t5 + 8*t6 &
                  + t7)*exp_minus_t*sqrt(2./21.))/3780.
else
ovlp(0) = ovlp(0) + (-8*k5*((-25920*k - 162240*k3 - 127680*k5 - 6720*k7 + 1920*k*(alpha)*(3 + 14*k2 &
                  + 7*k4)*t - 120*k*alpha2*(3 + k2)*(1 + 3*k2)*t2)*exp_minus_t + (25920*k + 162240*k3 &
                  + 127680*k5 + 6720*k7 + 120*(alpha)*(21 + 393*k2 + 511*k4 + 35*k6)*t + 240*k*(-1 &
                  + k2)**2*(21 + 54*k2 + 5*k4)*t2 + 40*alpha3*(5 + 38*k2 + 5*k4)*t3 + 20*k*(-1 &
                  + k2)**4*(5 + k2)*t4 + alpha5*(3 + k2)*t5)*exp(-(k*t)))*sqrt(2./21.)*sqrtk)/(45*(-1 &
                  + k2)**8*t)

endif

!*******************
case(1)   ! <s1|p2>
!*******************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((315*t2 + 315*t3 + 105*t4 - 12*t6 - 5*t7 - t8)*exp_minus_t*sqrt(2./7.))/(3780*t)
else
ovlp(0) = ovlp(0) + (32*k5*((-7200*k - 43200*k3 - 30240*k5 - 1440*k*(5 + 30*k2 + 21*k4)*t + 480*k*(-1 &
                  + k2)*(3 + 14*k2 + 7*k4)*t2 - 30*k*alpha2*(3 + k2)*(1 + 3*k2)*t3)*exp_minus_t &
                  + (7200*k + 43200*k3 + 30240*k5 + 1440*k2*(5 + 30*k2 + 21*k4)*t + 240*k*(alpha)*(9 &
                  + 62*k2 + 49*k4)*t2 + 120*alpha2*(1 + 18*k2 + 21*k4)*t3 + 30*k*alpha3*(5 &
                  + 11*k2)*t4 + 2*alpha4*(2 + 13*k2)*t5 + k*(-1 &
                  + k2)**5*t6)*exp(-(k*t)))*sqrt(2./7.)*sqrtk)/(45*alpha8*t2)
endif

!*********************
case(3)   ! <p1|s2>
!*********************

if(k.eq.1) then
ovlp(0) = ovlp(0) + (t*(7875 + 7875*t + 3759*t2 + 1134*t3 + 240*t4 + 37*t5 + 5*t6)*exp_minus_t*sqrt(2./7.))/18900.
else
ovlp(0) = ovlp(0) + (-8*k5*((-5184 - 116928*k2 - 188352*k4 - 12096*k6 - 576*(9 + 203*k2 + 327*k4 &
                  + 21*k6)*t + 120*(alpha)*(9 + 165*k2 + 203*k4 + 7*k6)*t2 - 24*alpha2*(3 &
                  + 42*k2 + 35*k4)*t3)*exp_minus_t + (5184 + 116928*k2 + 188352*k4 + 12096*k6 + 576*k*(9 &
                  + 203*k2 + 327*k4 + 21*k6)*t + 24*(alpha)*(63 + 1611*k2 + 2909*k4 + 217*k6)*t2 &
                  + 48*k*alpha2*(1 + 3*k2)*(91 + 9*k2)*t3 + 12*alpha3*(15 + 128*k2 &
                  + 17*k4)*t4 + 20*k*alpha4*(5 + k2)*t5 + alpha5*(3 &
                  + k2)*t6)*exp(-(k*t)))*sqrt(2./7.)*sqrtk)/(45*alpha8*t2)

endif

!*********************************
case(4)   ! l1 = l2 = 1 => <p1|p2>
!*********************************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((11025*t + 11025*t2 + 5040*t3 + 1365*t4 + 180*t5 - 30*t6 - 22*t7 &
                  - 5*t8)*exp_minus_t*sqrt(2./21.))/(6300*t)
ovlp(1) = ovlp(1) - ((-11025*t - 11025*t2 - 5355*t3 - 1680*t4 - 372*t5 - 57*t6 &
                  - 5*t7)*exp_minus_t*sqrt(2./21.))/(6300*t)
else
ovlp(0) = ovlp(0) + (32*k5*((-2880 - 63360*k2 - 95040*k4 - 2880*(1 + 22*k2 + 33*k4)*t - 1440*(1 + 22*k2 &
                  + 33*k4)*t2 + 30*(alpha)*(9 + 165*k2 + 203*k4 + 7*k6)*t3 - 6*alpha2*(3 &
                  + 42*k2 + 35*k4)*t4)*exp_minus_t + (2880 + 63360*k2 + 95040*k4 + 2880*k*(1 + 22*k2 &
                  + 33*k4)*t + 1440*k2*(1 + 22*k2 + 33*k4)*t2 + 240*k*(alpha)*(7 + 54*k2 &
                  + 59*k4)*t3 + 12*alpha2*(9 + 166*k2 + 225*k4)*t4 + 48*k*alpha3*(3 &
                  + 7*k2)*t5 + 2*alpha4*(2 + 13*k2)*t6 + k*(-1 &
                  + k2)**5*t7)*exp(-(k*t)))*sqrt(2./21.)*sqrtk)/(15*alpha8*t3)

ovlp(1) = ovlp(1) + (-32*k5*((-1440 - 31680*k2 - 47520*k4 - 1440*(1 + 22*k2 + 33*k4)*t + 288*(alpha)*(1 &
                  + 18*k2 + 21*k4)*t2 - 6*alpha2*(3 + 42*k2 + 35*k4)*t3)*exp_minus_t + (1440 &
                  + 31680*k2 + 47520*k4 + 1440*k*(1 + 22*k2 + 33*k4)*t + 144*(alpha)*(3 + 74*k2 &
                  + 123*k4)*t2 + 96*k*alpha2*(13 + 37*k2)*t3 + 6*alpha3*(9 + 71*k2)*t4 &
                  + 30*k*alpha4*t5 + alpha5*t6)*exp(-(k*t)))*sqrt(2./21.)*sqrtk)/(15*alpha8*t3)

endif

!********************
case(6)   ! <d1|s2>
!********************

if(k.eq.1) then
ovlp(0) = ovlp(0) + (t2*(1365 + 1365*t + 624*t2 + 169*t3 + 31*t4 + 5*t5)*exp_minus_t*sqrt(2./105.))/3780.
else
ovlp(0) = ovlp(0) + (-8*k5*((-190080*k - 725760*k3 - 51840*k5 - 17280*k*(11 + 42*k2 + 3*k4)*t + 192*k*(-387 &
                  - 1379*k2 + 79*k4 + 7*k6)*t2 + 192*k*(alpha)*(57 + 176*k2 + 7*k4)*t3 - 192*k*(-1 &
                  + k2)**2*(3 + 7*k2)*t4)*exp_minus_t + (190080*k + 725760*k3 + 51840*k5 + 17280*k2*(11 + 42*k2 &
                  + 3*k4)*t + 768*k*(-27 - 4*k2 + 419*k4 + 32*k6)*t2 + 384*k2*(alpha)*(54 + 227*k2 &
                  + 19*k4)*t3 + 48*k*alpha2*(63 + 306*k2 + 31*k4)*t4 + 4*alpha3*(35 + 392*k2 &
                  + 53*k4)*t5 + 20*k*alpha4*(5 + k2)*t6 + alpha5*(3 &
                  + k2)*t7)*exp(-(k*t)))*sqrt(2./105.)*sqrtk)/(9*alpha8*t3)

endif

!*********************
case(7)  ! <d1|p2>
!*********************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((4410*t2 + 4410*t3 + 1875*t4 + 405*t5 + 21*t6 - 16*t7 - 5*t8)*exp_minus_t*sqrt(2./35.))/(3780*t)
ovlp(1) = ovlp(1) + ((2205*t2 + 2205*t3 + 1020*t4 + 285*t5 + 51*t6 + 5*t7)*exp_minus_t*sqrt(2./105.))/(1260*t)
else
ovlp(0) = ovlp(0) + (32*k5*((-155520*k - 570240*k3 - 51840*k*(3 + 11*k2)*t + 2880*k*(-25 - 90*k2 + 3*k4)*t2 &
                  + 2880*k*(-7 - 24*k2 + 3*k4)*t3 + 48*k*(alpha)*(57 + 176*k2 + 7*k4)*t4 - 48*k*(-1 &
                  + k2)**2*(3 + 7*k2)*t5)*exp_minus_t + (155520*k + 570240*k3 + 51840*k2*(3 + 11*k2)*t + 5760*k*(-1 &
                  + 9*k2 + 48*k4)*t2 + 5760*k2*(-1 + 15*k4)*t3 + 96*k*(alpha)*(9 + 92*k2 + 199*k4)*t4 &
                  + 12*alpha2*(7 + 138*k2 + 255*k4)*t5 + 12*k*alpha3*(11 + 29*k2)*t6 + 2*(-1 &
                + k2)**4*(2 + 13*k2)*t7 + k*alpha5*t8)*exp(-(k*t)))*sqrt(2./35.)*sqrtk)/(9*alpha8*t4)

ovlp(1) = ovlp(1) + (-32*k5*((-51840*k - 190080*k3 - 17280*k*(3 + 11*k2)*t + 2880*k*(-7 - 24*k2 + 3*k4)*t2 &
                  + 2880*k*(alpha)*(1 + 3*k2)*t3 - 48*k*alpha2*(3 + 7*k2)*t4)*exp_minus_t + (51840*k &
                  + 190080*k3 + 17280*k2*(3 + 11*k2)*t + 5760*k*(-1 + 15*k4)*t2 + 5760*k2*(alpha)*(1 &
                  + 4*k2)*t3 + 96*k*alpha2*(9 + 41*k2)*t4 + 6*alpha3*(7 + 73*k2)*t5 + 30*k*(-1 &
                  + k2)**4*t6 + alpha5*t7)*exp(-(k*t)))*sqrt(2./105.)*sqrtk)/(3*alpha8*t4)

endif

end select

! ------------------------------ n1 = 4 and n2 = 3 -----------------------------
  case(11) 
! ------------------------------------------------------------------------------

select case(llopt)

!*********************************
case(0)   ! l1 = l2 = 0 => <s1|s2>
!*********************************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((11025 + 11025*t + 5250*t2 + 1575*t3 + 336*t4 + 56*t5 + 8*t6 &
                  + t7)*exp_minus_t*sqrt(2./7.))/6300 
else
ovlp(0) = ovlp(0) + (32*k4*((-336 - 8112*k2 - 15408*k4 - 3024*k6 + 6*(alpha)*(21 + 393*k2 &
                  + 511*k4 + 35*k6)*t - 6*alpha2*(3 + 42*k2 + 35*k4)*t2 + alpha3*(1 &
                  + 10*k2 + 5*k4)*t3)*exp_minus_t + (336 + 8112*k2 + 15408*k4 + 3024*k6 + 192*k*(-1 &
                  + k2)*(7 + 26*k2 + 7*k4)*t + 12*alpha2*(5 + 54*k2 + 21*k4)*t2 + 8*k*(-1 &
                  + k2)**3*(5 + 3*k2)*t3 + alpha4*(1 &
                  + k2)*t4)*exp(-(k*t)))*sqrt(2./7.)*sqrtk)/(15*alpha8*t)
endif

!*******************
case(1)   ! <s1|p2>
!*******************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((-1575*t2 - 1575*t3 - 861*t4 - 336*t5 - 96*t6 - 19*t7 - 3*t8)*&
                    exp_minus_t*sqrt(2./21.))/(6300*t)
else
ovlp(0) = ovlp(0) + (-32*k4*((-1296 - 29232*k2 - 47088*k4 - 3024*k6 - 144*(9 + 203*k2 + 327*k4 &
                  + 21*k6)*t + 12*(alpha)*(33 + 609*k2 + 763*k4 + 35*k6)*t2 - 18*alpha2*(3 &
                  + 42*k2 + 35*k4)*t3 + 3*alpha3*(1 + 10*k2 + 5*k4)*t4)*exp_minus_t + (1296 &
                  + 29232*k2 + 47088*k4 + 3024*k6 + 144*k*(9 + 203*k2 + 327*k4 + 21*k6)*t + 12*(-1 &
                  + k2)*(21 + 609*k2 + 1199*k4 + 91*k6)*t2 + 12*k*alpha2*(49 + 174*k2 &
                  + 17*k4)*t3 + 3*alpha3*(5 + 52*k2 + 7*k4)*t4 + k*alpha4*(5 &
                  + k2)*t5)*exp(-(k*t)))*sqrt(2./21.)*sqrtk)/(15*alpha8*t2)
endif

!*******************
case(2)   ! <s1|d2>
!*******************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((-105*t3 - 105*t4 - 36*t5 - t6 + 3*t7 + t8)*exp_minus_t*sqrt(2./35.))/(1260*t)
else
ovlp(0) = ovlp(0) + (32*k4*((-1440 - 31680*k2 - 47520*k4 - 1440*(1 + 22*k2 + 33*k4)*t + 48*(-13 &
                  - 271*k2 - 339*k4 + 63*k6)*t2 + 144*(alpha)*(1 + 18*k2 + 21*k4)*t3 - 6*(-1 &
                  + k2)**2*(3 + 42*k2 + 35*k4)*t4 + alpha3*(1 + 10*k2 + 5*k4)*t5)*exp_minus_t &
                  + (1440 + 31680*k2 + 47520*k4 + 1440*k*(1 + 22*k2 + 33*k4)*t + 96*(-1 - 22*k2 + 87*k4 &
                  + 216*k6)*t2 + 96*k*(alpha)*(1 + 28*k2 + 51*k4)*t3 + 6*alpha2*(1 + 46*k2 &
                  + 113*k4)*t4 + 2*k*alpha3*(5 + 27*k2)*t5 + 2*k2*(-1 &
                  + k2)**4*t6)*exp(-(k*t)))*sqrt(2./35.)*sqrtk)/(3*alpha8*t3)
endif

!*********************
case(3)   ! <p1|s2>
!*********************

if(k.eq.1) then
ovlp(0) = ovlp(0) + (t*(1575 + 1575*t + 735*t2 + 210*t3 + 42*t4 + 7*t5 + t6)*exp_minus_t*sqrt(2./21.))/2100
else
ovlp(0) = ovlp(0) + (32*k4*((-4320*k - 18240*k3 - 4320*k5 - 480*k*(9 + 38*k2 + 9*k4)*t + 96*k*(-1 &
                  + k2)*(12 + 41*k2 + 7*k4)*t2 - 4*k*alpha2*(33 + 82*k2 + 5*k4)*t3 &
                  + 2*k*alpha3*(3 + 5*k2)*t4)*exp_minus_t + (4320*k + 18240*k3 + 4320*k5 + 480*k2*(9 &
                  + 38*k2 + 9*k4)*t + 48*k*(alpha)*(21 + 108*k2 + 31*k4)*t2 + 2*alpha2*(25 &
                  + 326*k2 + 129*k4)*t3 + 8*k*alpha3*(5 + 3*k2)*t4 + alpha4*(1 &
                  + k2)*t5)*exp(-(k*t)))*sqrt(2./21.)*sqrtk)/(5*alpha8*t2)
endif

!*********************************
case(4)   ! l1 = l2 = 1 => <p1|p2>
!*********************************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((11025*t + 11025*t2 + 4410*t3 + 735*t4 - 60*t5 - 60*t6 - 16*t7 & 
                  - 3*t8)*exp_minus_t*sqrt(2./7.))/(6300*t)
ovlp(1) = ovlp(1) - ((-11025*t - 11025*t2 - 5145*t3 - 1470*t4 - 282*t5 - 37*t6 - 3*t7)*&
                    exp_minus_t*sqrt(2./7.))/(6300*t)
else
ovlp(0) = ovlp(0) + (-32*k4*((-31680*k - 120960*k3 - 8640*k5 - 2880*k*(11 + 42*k2 + 3*k4)*t - 1440*k*(11 &
                  + 42*k2 + 3*k4)*t2 + 240*k*(alpha)*(15 + 50*k2 + 7*k4)*t3 - 12*k*alpha2*(33 &
                  + 82*k2 + 5*k4)*t4 + 6*k*alpha3*(3 + 5*k2)*t5)*exp_minus_t + (31680*k + 120960*k3 &
                  + 8640*k5 + 2880*k2*(11 + 42*k2 + 3*k4)*t + 1440*k3*(11 + 42*k2 + 3*k4)*t2 + 30*(-1 &
                  + k2)*(7 + 203*k2 + 517*k4 + 41*k6)*t3 + 6*k*alpha2*(91 + 354*k2 &
                  + 35*k4)*t4 + 3*alpha3*(5 + 52*k2 + 7*k4)*t5 + k*alpha4*(5 &
                  + k2)*t6)*exp(-(k*t)))*sqrt(2./7.)*sqrtk)/(15*alpha8*t3)

ovlp(1) = ovlp(1) + (32*k4*((-15840*k - 60480*k3 - 4320*k5 - 1440*k*(11 + 42*k2 + 3*k4)*t + 96*k*(-1 &
                  + k2)*(42 + 131*k2 + 7*k4)*t2 - 144*k*alpha2*(3 + 7*k2)*t3 + 6*k*(-1 &
                  + k2)**3*(3 + 5*k2)*t4)*exp_minus_t + (15840*k + 60480*k3 + 4320*k5 + 1440*k2*(11 &
                  + 42*k2 + 3*k4)*t + 48*k*(alpha)*(81 + 368*k2 + 31*k4)*t2 + 6*alpha2*(35 &
                  + 402*k2 + 43*k4)*t3 + 24*k*alpha3*(7 + k2)*t4 + alpha4*(5 &
                  + k2)*t5)*exp(-(k*t)))*sqrt(2./7.)*sqrtk)/(15*alpha8*t3)

endif

!***************************************
case(5)   ! l1 = 1 and l2 = 2 => <p1|d2>
!***************************************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((-630*t2 - 630*t3 - 285*t4 - 75*t5 - 9*t6 + 2*t7 + t8)*exp_minus_t*sqrt(2./105.))/(420*t)
ovlp(1) = ovlp(1) - ((315*t2 + 315*t3 + 150*t4 + 45*t5 + 9*t6 + t7)*exp_minus_t*sqrt(2./35.))/(420*t)
else
ovlp(0) = ovlp(0) + (64*k4*((-25920*k - 95040*k3 - 8640*k*(3 + 11*k2)*t + 240*k*(-51 - 182*k2 &
                  + 9*k4)*t2 + 240*k*(-15 - 50*k2 + 9*k4)*t3 + 24*k*(alpha)*(27 + 86*k2 &
                  + 7*k4)*t4 - 2*k*alpha2*(33 + 82*k2 + 5*k4)*t5 + k*alpha3*(3 &
                  + 5*k2)*t6)*exp_minus_t + (25920*k + 95040*k3 + 8640*k2*(3 + 11*k2)*t + 240*k*(-3 &
                  + 38*k2 + 189*k4)*t2 + 240*k2*(-3 + 2*k2 + 57*k4)*t3 + 24*k*(alpha)*(3 &
                  + 44*k2 + 113*k4)*t4 + 4*alpha2*(1 + 3*k2)*(1 + 29*k2)*t5 + k*(-1 &
                  + k2)**3*(5 + 27*k2)*t6 + k2*alpha4*t7)*exp(-(k*t)))*sqrt(2./105.)*sqrtk)/((-1 &
                  + k2)**8*t4)

ovlp(1) = ovlp(1) + (-64*k4*((-8640*k - 31680*k3 - 2880*k*(3 + 11*k2)*t + 240*k*(-15 - 50*k2 + 9*k4)*t2 &
                  + 720*k*(alpha)*(1 + 3*k2)*t3 - 24*k*alpha2*(3 + 7*k2)*t4 + k*alpha3*(3 &
                  + 5*k2)*t5)*exp_minus_t + (8640*k + 31680*k3 + 2880*k2*(3 + 11*k2)*t + 240*k*(-3 + 2*k2 &
                  + 57*k4)*t2 + 240*k2*(alpha)*(3 + 13*k2)*t3 + 24*k*alpha2*(3 + 17*k2)*t4 &
                + 2*alpha3*(1 + 15*k2)*t5 + k*alpha4*t6)*exp(-(k*t)))*sqrt(2./35.)*sqrtk)/((-1 &
                  + k2)**8*t4)

endif

!********************
case(6)   ! <d1|s2>
!********************

if(k.eq.1) then
ovlp(0) = ovlp(0) + (t2*(105 + 105*t + 54*t2 + 19*t3 + 5*t4 + t5)*exp_minus_t*sqrt(2./35.))/1260
else
ovlp(0) = ovlp(0) + (32*k4*((-4320 - 60480*k2 - 15840*k4 - 1440*(3 + 42*k2 + 11*k4)*t + 96*(-19 &
                  - 253*k2 - 17*k4 + 9*k6)*t2 + 96*(alpha)*(4 + 47*k2 + 9*k4)*t3 - 6*(-1 &
                  + k2)**2*(7 + 66*k2 + 7*k4)*t4 + 2*alpha3*(1 + 7*k2)*t5)*exp_minus_t + (4320 &
                  + 60480*k2 + 15840*k4 + 1440*k*(3 + 42*k2 + 11*k4)*t + 48*(-7 - 79*k2 + 499*k4 &
                  + 147*k6)*t2 + 48*k*(alpha)*(7 + 116*k2 + 37*k4)*t3 + 30*alpha2*(1 + 22*k2 &
                  + 9*k4)*t4 + 8*k*alpha3*(5 + 3*k2)*t5 + alpha4*(1 &
                  + k2)*t6)*exp(-(k*t)))*sqrt(2./35.)*sqrtk)/(3*alpha8*t3)

endif

!*********************
case(7)  ! <d1|p2> 
!*********************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((3150*t2 + 3150*t3 + 1281*t4 + 231*t5 - 3*t6 - 10*t7 - 3*t8)*exp_minus_t*&
                    sqrt(2./105.))/(1260*t)
ovlp(1) = ovlp(1) + ((1575*t2 + 1575*t3 + 714*t4 + 189*t5 + 31*t6 + 3*t7)*exp_minus_t*& 
                    sqrt(2./35.))/(1260*t)
else
ovlp(0) = ovlp(0) + (-32*k4*((-47520 - 630720*k2 - 47520*k4 - 4320*(11 + 146*k2 + 11*k4)*t + 288*(-78 &
                  - 1021*k2 - 24*k4 + 3*k6)*t2 + 288*(-23 - 291*k2 + 31*k4 + 3*k6)*t3 + 6*(-1 &
                  + k2)*(201 + 2325*k2 + 347*k4 + 7*k6)*t4 - 18*alpha2*(7 + 66*k2 + 7*k4)*t5 &
                  + 6*alpha3*(1 + 7*k2)*t6)*exp_minus_t + (47520 + 630720*k2 + 47520*k4 + 4320*k*(11 &
                  + 146*k2 + 11*k4)*t + 144*(-9 + 17*k2 + 2073*k4 + 159*k6)*t2 + 144*k*(-9 - 93*k2 &
                  + 613*k4 + 49*k6)*t3 + 6*(alpha)*(21 + 609*k2 + 2959*k4 + 251*k6)*t4 + 6*k*(-1 &
                  + k2)**2*(77 + 366*k2 + 37*k4)*t5 + 3*alpha3*(5 + 52*k2 + 7*k4)*t6 + k*(-1 &
                  + k2)**4*(5 + k2)*t7)*exp(-(k*t)))*sqrt(2./105.)*sqrtk)/(3*alpha8*t4)
 
ovlp(1) = ovlp(1) + (32*k4*((-15840 - 210240*k2 - 15840*k4 - 1440*(11 + 146*k2 + 11*k4)*t + 288*(-23 &
                  - 291*k2 + 31*k4 + 3*k6)*t2 + 96*(alpha)*(14 + 157*k2 + 9*k4)*t3 - 6*(-1 &
                  + k2)**2*(23 + 210*k2 + 7*k4)*t4 + 6*alpha3*(1 + 7*k2)*t5)*exp_minus_t + (15840 &
                  + 210240*k2 + 15840*k4 + 1440*k*(11 + 146*k2 + 11*k4)*t + 144*(-9 - 93*k2 + 613*k4 &
                  + 49*k6)*t2 + 48*k*(alpha)*(27 + 416*k2 + 37*k4)*t3 + 18*alpha2*(7 &
                  + 138*k2 + 15*k4)*t4 + 24*k*alpha3*(7 + k2)*t5 + alpha4*(5 &
                  + k2)*t6)*exp(-(k*t)))*sqrt(2./35.)*sqrtk)/(3*alpha8*t4)

endif

!*********************
case(8)  ! <d1|d2> 
!*********************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((2205*t + 2205*t2 + 690*t3 - 45*t4 - 84*t5 - 20*t6 + t8)*exp_minus_t*sqrt(2./7.))/(1260*t)
ovlp(1) = ovlp(1) - ((-735*t - 735*t2 - 255*t3 - 10*t4 + 20*t5 + 7*t6 + t7)*exp_minus_t*sqrt(2./7.))/(420*t)
ovlp(2) = ovlp(2) + ((735*t + 735*t2 + 330*t3 + 85*t4 + 13*t5 + t6)*exp_minus_t*sqrt(2./7.))/(420*t)
else
ovlp(0) = ovlp(0) + (64*k4*((-51840 - 673920*k2 - 51840*(1 + 13*k2)*t + 1440*(-17 - 218*k2 + 11*k4)*t2 &
                  + 1440*(-5 - 62*k2 + 11*k4)*t3 - 48*(31 + 367*k2 - 127*k4 + 9*k6)*t4 + 3*(-1 &
                  + k2)*(73 + 821*k2 + 59*k4 + 7*k6)*t5 - 3*alpha2*(7 + 66*k2 + 7*k4)*t6 &
                  + alpha3*(1 + 7*k2)*t7)*exp_minus_t + (51840 + 673920*k2 + 51840*k*(1 + 13*k2)*t &
                  + 1440*(-1 + 2*k2 + 223*k4)*t2 + 1440*k*(-1 - 10*k2 + 67*k4)*t3 + 48*(1 + 7*k2 &
                  - 157*k4 + 429*k6)*t4 + 120*k*(alpha)*(1 + 4*k2 + 27*k4)*t5 + 6*alpha2*(1 &
                  + 18*k2 + 61*k4)*t6 + k*alpha3*(5 + 27*k2)*t7 + k2*(-1 &
                  + k2)**4*t8)*exp(-(k*t)))*sqrt(2./7.)*sqrtk)/(3*alpha8*t5)
 
ovlp(1) = ovlp(1) + (-64*k4*((-11520 - 149760*k2 - 11520*(1 + 13*k2)*t + 480*(-11 - 140*k2 + 11*k4)*t2 &
                  + 480*(-3 - 36*k2 + 11*k4)*t3 + 240*(alpha)*(1 + 11*k2)*t4 - alpha2*(23 &
                  + 210*k2 + 7*k4)*t5 + alpha3*(1 + 7*k2)*t6)*exp_minus_t + (11520 + 149760*k2 &
                  + 11520*k*(1 + 13*k2)*t + 480*(-1 - 4*k2 + 145*k4)*t2 + 480*k*(-1 - 12*k2 + 41*k4)*t3 &
                  + 240*k2*(alpha)*(1 + 15*k2)*t4 + 8*k*alpha2*(7 + 53*k2)*t5 + 2*(-1 &
                  + k2)**3*(1 + 15*k2)*t6 + k*alpha4*t7)*exp(-(k*t)))*sqrt(2./7.)*sqrtk)/((-1 &
                  + k2)**8*t5)

ovlp(2) = ovlp(2) + (64*k4*((-2880 - 37440*k2 - 2880*(1 + 13*k2)*t + 240*(-5 - 62*k2 + 11*k4)*t2 &
                  + 240*(alpha)*(1 + 11*k2)*t3 - 24*alpha2*(1 + 9*k2)*t4 + alpha3*(1 &
                  + 7*k2)*t5)*exp_minus_t + (2880 + 37440*k2 + 2880*k*(1 + 13*k2)*t + 240*(-1 - 10*k2 &
                  + 67*k4)*t2 + 240*k*(alpha)*(1 + 15*k2)*t3 + 24*alpha2*(1 + 19*k2)*t4 &
               + 32*k*alpha3*t5 + alpha4*t6)*exp(-(k*t)))*sqrt(2./7.)*sqrtk)/(alpha8*t5)

endif

end select

! ------------------------------ n1 = n2 = 3 -----------------------------------
  case(10) 
! ------------------------------------------------------------------------------

select case(llopt)

!*********************************
case(0)  ! l1 = l2 = 0 => <s1|s2>
!*********************************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((1575 + 1575*t + 735*t2 + 210*t3 + 42*t4 + 7*t5 + t6)*exp_minus_t)/1575
else
ovlp(0) = ovlp(0) + (64*k3*((-336*k - 1248*k3 - 336*k5 + 6*k*(alpha)*(21 + 54*k2 + 5*k4)*t & 
         - 6*k*alpha2*(3 + 5*k2)*t2 + k*alpha3*(1 + k2)*t3)*exp_minus_t  & 
         + (336*k + 1248*k3 + 336*k5 + 6*(alpha)*(5 + 54*k2 + 21*k4)*t + 6*k*alpha2*(5 + 3*k2)*t2 &
         + alpha3*(1 + k2)*t3)*exp(-(k*t)))*sqrtk)/(15*alpha7*t)
endif

!*****************************************
case(1,3)  ! l1 = 1 and l2 = 0  => <p1|s2> 
!*****************************************

if(k.eq.1) then
ovlp(0) = ovlp(0) + (t*(525 + 525*t + 252*t2 + 77*t3 + 17*t4 + 3*t5)*exp_minus_t)/(1575*sqrt(3.))
else
ovlp(0) = ovlp(0) + (64*k3*((-336 - 4128*k2 - 1296*k4 - 48*(7 + 86*k2 + 27*k4)*t + &
                    36*(alpha)*(3 + 30*k2 + 7*k4)*t2 - 2*alpha2*(8 + 59*k2 + 5*k4)*t3 + &
                    alpha3*(1 + 5*k2)*t4)*exp_minus_t + (336 + 4128*k2 + 1296*k4 + 48*k*(7 + 86*k2 + &
                    27*k4)*t + 12*(alpha)*(5 + 82*k2 + 33*k4)*t2 + 18*k*alpha2*(5 + 3*k2)*t3 + &
                    3*alpha3*(1 + k2)*t4)*exp(-(k*t)))*sqrtk)/(15*sqrt(3.)*alpha7*t2)
endif

!***************************************
case(2,6)  ! l1 = 2 and l2 = 0  => <d1|s2> 
!***************************************

if(k.eq.1) then
ovlp(0) = ovlp(0) + (t4*(3 + 3*t + t2)*exp_minus_t)/(315*sqrt(5.))
else
ovlp(0) = ovlp(0) + (64*k3*((-4320*k - 1440*k3 - 1440*k*(3 + k2)*t + 96*k*(-19 - 2*k2 + k4)*t2 &
                  + 96*k*(alpha)*(4 + k2)*t3 - 6*k*alpha2*(7 + k2)*t4 &
                  + 2*k*alpha3*t5)*exp_minus_t + (4320*k + 1440*k3 + 1440*k2*(3 + k2)*t &
                  + 48*k*(-7 + 34*k2 + 13*k4)*t2 + 48*k2*(alpha)*(7 + 3*k2)*t3 &
                  + 6*k*alpha2*(5 + 3*k2)*t4 + alpha3*(1 + k2)*t5) &
                  *exp(-(k*t)))*sqrtk)/(3*sqrt(5.)*alpha7*t3)
endif

!***************************************
case(4)  ! l1 = l2 = 1  => <p1|p2> 
!***************************************

if(k.eq.1) then
ovlp(0) = ovlp(0) + (( 1575*t + 1575*t2 + 567*t3 + 42*t4 - 34*t5 - 13*t6 - 3*t7)*exp_minus_t)/(1575*t)
ovlp(1) = ovlp(1) - ((-1575*t - 1575*t2 - 714*t3 - 189*t4 - 31*t5 - 3*t6)*exp_minus_t)/(1575*t)
else
ovlp(0) = ovlp(0) + (-64*k3*((-864 - 9792*k2 - 864*k4 - 288*(3 + 34*k2 + 3*k4)*t &
                  - 144*(3 + 34*k2 + 3*k4)*t2 + 6*(alpha)*(19 + 186*k2 + 35*k4)*t3 &
                  - 2*alpha2*(8 + 59*k2 + 5*k4)*t4 + alpha3*(1 + 5*k2)*t5)*exp_minus_t &
                  + (864 + 9792*k2 + 864*k4 + 288*k*(3 + 34*k2 + 3*k4)*t + 144*k2*(3 + 34*k2 &
                  + 3*k4)*t2 + 6*k*(alpha)*(35 + 186*k2 + 19*k4)*t3 + 2*alpha2*(5 + 59*k2 &
                  + 8*k4)*t4 + k*alpha3*(5 + k2)*t5)*exp(-(k*t)))*sqrtk)/(15*alpha7*t3)

ovlp(1) = ovlp(1) + (64*k3*((-432 - 4896*k2 - 432*k4 - 144*(3 + 34*k2 + 3*k4)*t &
                  + 12*(alpha)*(11 + 102*k2 + 7*k4)*t2 - 18*alpha2*(1 + 7*k2)*t3 &
                  + alpha3*(1 + 5*k2)*t4)*exp_minus_t + (432 + 4896*k2 + 432*k4 &
                  + 144*k*(3 + 34*k2 + 3*k4)*t + 12*(alpha)*(7 + 102*k2 + 11*k4)*t2 &
                  + 18*k*alpha2*(7 + k2)*t3 + alpha3*(5 + k2)*t4)*exp(-(k*t)))*&
                    sqrtk)/(15*alpha7*t3)
endif

!***************************************
case(5,7)  ! l1 = 2 and l2 = 1  => <d1|p2> 
!***************************************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((420*t2 + 420*t3 + 162*t4 + 22*t5 - 5*t6 - 3*t7)*exp_minus_t)/(315*t*sqrt(15.))
ovlp(1) = ovlp(1) + ((210*t2 + 210*t3 + 93*t4 + 23*t5 + 3*t6)*exp_minus_t)/(315*t*sqrt(5.))
else
ovlp(0) = ovlp(0) + (-64*k3*((-47520*k - 4320*k3 - 4320*k*(11 + k2)*t + 96*k*(-234 - 7*k2 + k4)*t2 &
                  + 96*k*(-69 + 8*k2 + k4)*t3 + 6*k*(alpha)*(201 + 38*k2 + k4)*t4 &
                  - 18*k*alpha2*(7 + k2)*t5 + 6*k*alpha3*t6)*exp_minus_t + (47520*k + 4320*k3 &
                  + 4320*k2*(11 + k2)*t + 48*k*(-27 + 464*k2 + 43*k4)*t2 + 48*k2*(-27 + 134*k2 &
                  + 13*k4)*t3 + 18*k*(alpha)*(7 + 66*k2 + 7*k4)*t4 + 2*alpha2*(5 + 59*k2 &
                  + 8*k4)*t5 + k*alpha3*(5 + k2)*t6)*exp(-(k*t)))*sqrtk)/(3*Sqrt(15.)*alpha7*t4)

ovlp(1) = ovlp(1) + (64*k3*((-15840*k - 1440*k3 - 1440*k*(11 + k2)*t + 96*k*(-69 + 8*k2 + k4)*t2 &
                  + 96*k*(alpha)*(14 + k2)*t3 - 6*k*alpha2*(23 + k2)*t4 &
                  + 6*k*alpha3*t5)*exp_minus_t + (15840*k + 1440*k3 + 1440*k2*(11 + k2)*t &
                  + 48*k*(-27 + 134*k2 + 13*k4)*t2 + 144*k2*(alpha)*(9 + k2)*t3 &
                  + 18*k*alpha2*(7 + k2)*t4 + alpha3*(5 + k2)*t5)&
                  *exp(-(k*t)))*sqrtk)/(3*Sqrt(5.)*alpha7*t4)
endif

!***************************************
case(8)  ! l1 = l2 = 2  => <d1|d2> 
!***************************************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((315*t + 315*t2 + 75*t3 - 30*t4 - 18*t5 - t6 + t7)*exp_minus_t)/(315*t)
ovlp(1) = ovlp(1) - ((-105*t - 105*t2 - 30*t3 + 5*t4 + 5*t5 + t6)*exp_minus_t)/(105*t)
ovlp(2) = ovlp(2) + ((105*t + 105*t2 + 45*t3 + 10*t4 + t5)*exp_minus_t)/(105*t)
else
ovlp(0) = ovlp(0) + (128*k3*((-51840*k - 51840*k*t + 1440*k*(-17 + k2)*t2 + 1440*k*(-5 + k2)*t3 &
                  - 48*k*(31 - 12*k2 + k4)*t4 + 3*k*(alpha)*(73 + 6*k2 + k4)*t5 &
                  - 3*k*alpha2*(7 + k2)*t6 + k*alpha3*t7)*exp_minus_t + (51840*k &
                  + 51840*k2*t + 1440*k*(-1 + 17*k2)*t2 + 1440*k2*(-1 + 5*k2)*t3 + 48*k*(1 &
                  - 12*k2 + 31*k4)*t4 + 3*(alpha)*(1 + 6*k2 + 73*k4)*t5 + 3*k*alpha2*(1 &
                  + 7*k2)*t6 + k2*alpha3*t7)*exp(-(k*t)))*sqrtk)/(3*alpha7*t5)

ovlp(1) = ovlp(1) + (-128*k3*((-11520*k - 11520*k*t + 480*k*(-11 + k2)*t2 + 480*k*(-3 + k2)*t3 &
                  + 240*k*(alpha)*t4 - k*alpha2*(23 + k2)*t5 + k*alpha3*t6)*exp_minus_t &
                  + (11520*k + 11520*k2*t + 480*k*(-1 + 11*k2)*t2 + 480*k2*(-1 + 3*k2)*t3 &
                  + 240*k3*(alpha)*t4 + alpha2*(1 + 23*k2)*t5 &
                  + k*alpha3*t6)*exp(-(k*t)))*sqrtk)/(alpha7*t5)

ovlp(2) = ovlp(2) + (128*k3*((-2880*k - 2880*k*t + 240*k*(-5 + k2)*t2 + 240*k*(alpha)*t3 &
                  - 24*k*alpha2*t4 + k*alpha3*t5)*exp_minus_t + (2880*k + 2880*k2*t &
                  + 240*k*(-1 + 5*k2)*t2 + 240*k2*(alpha)*t3 + 24*k*alpha2*t4 &
                  + alpha3*t5)*exp(-(k*t)))*sqrtk)/(alpha7*t5)
endif


end select  

! ------------------------------ n1 = 5 and n2 = 3 -----------------------------
  case(12) 
! ------------------------------------------------------------------------------

select case(llopt)

!*********************************
case(0)   ! l1 = l2 = 0 => <s1|s2>
!*********************************

if(k.eq.1) then
ovlp(0) = ovlp(0) +  ((79380 + 79380*t + 38745*t2 + 12285*t3 + 2835*t4 + 504*t5 + 72*t6 + 9*t7 &
                  + t8)*exp_minus_t)/(17010*sqrt(35.))
else

ovlp(0) = ovlp(0) + (64*k5*((-21600*k - 175200*k3 - 203040*k5 - 30240*k7 + 240*k*(alpha)*(27 + 169*k2 + 133*k4 &
                  + 7*k6)*t - 240*k*alpha2*(3 + 14*k2 + 7*k4)*t2 + 10*k*alpha3*(3 + k2)*(1 &
                  + 3*k2)*t3)*exp_minus_t + (21600*k + 175200*k3 + 203040*k5 + 30240*k7 + 240*(alpha)*(7 + 169*k2 &
                  + 321*k4 + 63*k6)*t + 480*k*alpha2*(7 + 26*k2 + 7*k4)*t2 + 20*alpha3*(5 + 54*k2 &
                  + 21*k4)*t3 + 10*k*alpha4*(5 + 3*k2)*t4 + alpha5*(1 &
                  + k2)*t5)*exp(-(k*t)))*sqrtk)/(45*sqrt(35.)*alpha9*t)

endif

!*******************
case(1)   ! <s1|p2>
!*******************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((-315*t4 - 315*t5 - 144*t6 - 39*t7 - 7*t8 - t9)*exp_minus_t)/(5670*t*sqrt(105.))

else

ovlp(0) = ovlp(0) + (-64*k5*((-79200*k - 597600*k3 - 583200*k5 - 30240*k7 - 1440*k*(55 + 415*k2 + 405*k4 &
                  + 21*k6)*t + 3360*k*(alpha)*(6 + 37*k2 + 28*k4 + k6)*t2 - 720*k*alpha2*(3 &
                  + 14*k2 + 7*k4)*t3 + 30*k*alpha3*(3 + k2)*(1 + 3*k2)*t4)*exp_minus_t + (79200*k &
                  + 597600*k3 + 583200*k5 + 30240*k7 + 1440*k2*(55 + 415*k2 + 405*k4 + 21*k6)*t + 240*k*(-1 &
                  + k2)*(81 + 727*k2 + 823*k4 + 49*k6)*t2 + 120*alpha2*(7 + 167*k2 + 285*k4 &
                  + 21*k6)*t3 + 30*k*alpha3*(35 + 114*k2 + 11*k4)*t4 + 2*alpha4*(10 + 97*k2 &
                  + 13*k4)*t5 + k*alpha5*(5 + k2)*t6)*exp(-(k*t)))*sqrtk)/(45*sqrt(105.)*alpha9*t2)

endif

!*******************
case(2)   ! <s1|d2>
!*******************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((-945*t3 - 945*t4 - 405*t5 - 90*t6 - 6*t7 + 3*t8 + t9)*exp_minus_t)/(17010*t*sqrt(7.))

else

ovlp(0) = ovlp(0) + (128*k5*((-43200*k - 316800*k3 - 285120*k5 - 2880*k*(15 + 110*k2 + 99*k4)*t + 240*k*(-75 &
                  - 515*k2 - 369*k4 + 63*k6)*t2 + 720*k*(alpha)*(5 + 30*k2 + 21*k4)*t3 - 120*k*(-1 &
                  + k2)**2*(3 + 14*k2 + 7*k4)*t4 + 5*k*alpha3*(3 + k2)*(1 + 3*k2)*t5)*exp_minus_t &
                  + (43200*k + 316800*k3 + 285120*k5 + 2880*k2*(15 + 110*k2 + 99*k4)*t + 240*k*(-15 - 55*k2 &
                  + 435*k4 + 531*k6)*t2 + 1200*k2*(alpha)*(3 + 26*k2 + 27*k4)*t3 + 120*k*alpha2*(3 &
                  + 34*k2 + 43*k4)*t4 + 8*alpha3*(1 + 33*k2 + 66*k4)*t5 + k*alpha4*(7 &
                  + 33*k2)*t6 + k2*alpha5*t7)*exp(-(k*t)))*sqrtk)/(45*sqrt(7.)*alpha9*t3)

endif

!*********************
case(3)   ! <p1|s2>
!*********************

if(k.eq.1) then
ovlp(0) = ovlp(0) + (t*(66150 + 66150*t + 31185*t2 + 9135*t3 + 1872*t4 + 297*t5 + 41*t6 &
                  + 5*t7)*exp_minus_t)/(28350*sqrt(105.))
else

ovlp(0) = ovlp(0) + (64*k5*((-4320 - 119520*k2 - 258720*k4 - 47520*k6 - 480*(9 + 249*k2 + 539*k4 + 99*k6)*t &
                  + 288*(alpha)*(4 + 93*k2 + 162*k4 + 21*k6)*t2 - 4*alpha2*(33 + 609*k2 + 763*k4 &
                  + 35*k6)*t3 + 2*alpha3*(3 + 42*k2 + 35*k4)*t4)*exp_minus_t + (4320 + 119520*k2 + 258720*k4 &
                  + 47520*k6 + 480*k*(9 + 249*k2 + 539*k4 + 99*k6)*t + 48*(alpha)*(21 + 687*k2 + 1723*k4 &
                  + 369*k6)*t2 + 32*k*alpha2*(91 + 398*k2 + 111*k4)*t3 + 2*alpha3*(45 + 542*k2 &
                  + 213*k4)*t4 + 10*k*alpha4*(5 + 3*k2)*t5 + alpha5*(1 &
                  + k2)*t6)*exp(-(k*t)))*sqrtk)/(15*sqrt(105.)*alpha9*t2)

endif

!*********************************
case(4)   ! l1 = l2 = 1 => <p1|p2>
!*********************************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((132300*t + 132300*t2 + 57645*t3 + 13545*t4 + 1305*t5 - 270*t6 - 142*t7 &
                  - 31*t8 - 5*t9)*exp_minus_t)/(28350*t*sqrt(35.))
ovlp(1) = ovlp(1) - ((-132300*t - 132300*t2 - 63315*t3 - 19215*t4 - 4104*t5 - 639*t6 - 71*t7 &
                  - 5*t8)*exp_minus_t)/(28350*t*sqrt(35.))
else

ovlp(0) = ovlp(0) + (-64*k5*((-31680 - 838080*k2 - 1615680*k4 - 95040*k6 - 2880*(11 + 291*k2 + 561*k4 + 33*k6)*t &
                  - 1440*(11 + 291*k2 + 561*k4 + 33*k6)*t2 + 720*(alpha)*(5 + 115*k2 + 195*k4 + 21*k6)*t3 &
                  - 12*alpha2*(33 + 609*k2 + 763*k4 + 35*k6)*t4 + 6*alpha3*(3 + 42*k2 &
                  + 35*k4)*t5)*exp_minus_t + (31680 + 838080*k2 + 1615680*k4 + 95040*k6 + 2880*k*(11 + 291*k2 + 561*k4 &
                  + 33*k6)*t + 1440*k2*(11 + 291*k2 + 561*k4 + 33*k6)*t2 + 240*k*(alpha)*(63 + 629*k2 &
                  + 929*k4 + 59*k6)*t3 + 36*alpha2*(21 + 513*k2 + 991*k4 + 75*k6)*t4 + 48*k*(-1 &
                  + k2)**3*(21 + 72*k2 + 7*k4)*t5 + 2*alpha4*(10 + 97*k2 + 13*k4)*t6 + k*alpha5*(5 &
                  + k2)*t7)*exp(-(k*t)))*sqrtk)/(45*sqrt(35.)*alpha9*t3)

ovlp(1) = ovlp(1) + (64*k5*((-15840 - 419040*k2 - 807840*k4 - 47520*k6 - 1440*(11 + 291*k2 + 561*k4 + 33*k6)*t &
                + 288*(alpha)*(14 + 313*k2 + 492*k4 + 21*k6)*t2 - 432*alpha2*(1 + 18*k2 + 21*k4)*t3 &
                + 6*alpha3*(3 + 42*k2 + 35*k4)*t4)*exp_minus_t + (15840 + 419040*k2 + 807840*k4 + 47520*k6 &
                + 1440*k*(11 + 291*k2 + 561*k4 + 33*k6)*t + 144*(alpha)*(27 + 829*k2 + 1821*k4 + 123*k6)*t2 &
                + 96*k*alpha2*(117 + 446*k2 + 37*k4)*t3 + 6*alpha3*(63 + 666*k2 + 71*k4)*t4 &
                + 30*k*alpha4*(7 + k2)*t5 + alpha5*(5 + k2)*t6)*exp(-(k*t)))*sqrtk)/(45*sqrt(35.)*(-1 &
                + k2)**9*t3)

endif

!***************************************
case(5)   ! l1 = 1 and l2 = 2 => <p1|d2>
!***************************************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((-13230*t2 - 13230*t3 - 6885*t4 - 2475*t5 - 594*t6 - 63*t7 + 11*t8 &
                  + 5*t9)*exp_minus_t)/(28350*t*sqrt(21.))
ovlp(1) = ovlp(1) - ((6615*t2 + 6615*t3 + 3375*t4 + 1170*t5 + 294*t6 + 51*t7 &
                  + 5*t8)*exp_minus_t)/(28350*t*sqrt(7.))
else

ovlp(0) = ovlp(0) + (128*k5*((-25920 - 673920*k2 - 1235520*k4 - 8640*(3 + 78*k2 + 143*k4)*t + 240*(-51 - 1311*k2 &
                  - 2321*k4 + 99*k6)*t2 + 240*(-15 - 375*k2 - 605*k4 + 99*k6)*t3 + 72*(alpha)*(9 + 203*k2 &
                  + 327*k4 + 21*k6)*t4 - 2*alpha2*(33 + 609*k2 + 763*k4 + 35*k6)*t5 + alpha3*(3 &
                  + 42*k2 + 35*k4)*t6)*exp_minus_t + (25920 + 673920*k2 + 1235520*k4 + 8640*k*(3 + 78*k2 + 143*k4)*t &
                  + 240*(-3 - 39*k2 + 1151*k4 + 2475*k6)*t2 + 240*k*(-3 - 75*k2 + 215*k4 + 759*k6)*t3 + 24*(-1 &
                  + k2)*(3 + 141*k2 + 1069*k4 + 1587*k6)*t4 + 8*k*alpha2*(47 + 466*k2 + 687*k4)*t5 + (-1 &
                  + k2)**3*(9 + 254*k2 + 537*k4)*t6 + k*alpha4*(7 + 33*k2)*t7 + k2*(-1 &
                  + k2)**5*t8)*exp(-(k*t)))*sqrtk)/(15*sqrt(21.)*alpha9*t4)

ovlp(1) = ovlp(1) + (-128*k5*((-8640 - 224640*k2 - 411840*k4 - 2880*(3 + 78*k2 + 143*k4)*t + 240*(-15 - 375*k2 &
                  - 605*k4 + 99*k6)*t2 + 720*(alpha)*(1 + 22*k2 + 33*k4)*t3 - 72*alpha2*(1 + 18*k2 &
                  + 21*k4)*t4 + alpha3*(3 + 42*k2 + 35*k4)*t5)*exp_minus_t + (8640 + 224640*k2 + 411840*k4 &
                  + 2880*k*(3 + 78*k2 + 143*k4)*t + 240*(-3 - 75*k2 + 215*k4 + 759*k6)*t2 + 240*k*(alpha)*(3 &
                  + 90*k2 + 187*k4)*t3 + 24*alpha2*(3 + 114*k2 + 283*k4)*t4 + 8*k*alpha3*(19 &
                  + 81*k2)*t5 + alpha4*(3 + 37*k2)*t6 + k*(-1 &
                  + k2)**5*t7)*exp(-(k*t)))*sqrtk)/(15*sqrt(7.)*alpha9*t4)

endif

!********************
case(6)   ! <d1|s2>
!********************

if(k.eq.1) then
ovlp(0) = ovlp(0) + (t2*(6615 + 6615*t + 3105*t2 + 900*t3 + 186*t4 + 33*t5 + 5*t6)*exp_minus_t)/(85050*sqrt(7.))

else

ovlp(0) = ovlp(0) + (64*k5*((-190080*k - 910080*k3 - 190080*k5 - 5760*k*(33 + 158*k2 + 33*k4)*t + 960*k*(-81 &
                  - 361*k2 - 15*k4 + 9*k6)*t2 + 2880*k*(alpha)*(5 + 20*k2 + 3*k4)*t3 - 48*k*(-1 &
                  + k2)**2*(27 + 86*k2 + 7*k4)*t4 + 16*k*alpha3*(3 + 7*k2)*t5)*exp_minus_t + (190080*k &
                  + 910080*k3 + 190080*k5 + 5760*k2*(33 + 158*k2 + 33*k4)*t + 1920*k*(-9 - 7*k2 + 195*k4 &
                  + 45*k6)*t2 + 1920*k2*(alpha)*(9 + 49*k2 + 12*k4)*t3 + 96*k*alpha2*(21 + 138*k2 &
                  + 41*k4)*t4 + 2*alpha3*(35 + 546*k2 + 219*k4)*t5 + 10*k*alpha4*(5 + 3*k2)*t6 &
                  + alpha5*(1 + k2)*t7)*exp(-(k*t)))*sqrtk)/(45*sqrt(7.)*alpha9*t3)

endif

!*********************
case(7)  ! <d1|p2> 
!*********************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((39690*t2 + 39690*t3 + 17055*t4 + 3825*t5 + 342*t6 - 51*t7 - 23*t8 &
                  - 5*t9)*exp_minus_t)/(28350*t*sqrt(21.))
ovlp(1) = ovlp(1) + ((19845*t2 + 19845*t3 + 9225*t4 + 2610*t5 + 492*t6 + 63*t7 + 5*t8)*exp_minus_t)/(28350*t*sqrt(7.))

else

ovlp(0) = ovlp(0) + (-64*k5*((-2021760*k - 9020160*k3 - 570240*k5 - 51840*k*(39 + 174*k2 + 11*k4)*t + 2880*k*(-329 &
                  - 1439*k2 - 27*k4 + 3*k6)*t2 + 2880*k*(-95 - 395*k2 + 39*k4 + 3*k6)*t3 + 48*k*(alpha)*(933 &
                  + 3661*k2 + 439*k4 + 7*k6)*t4 - 144*k*alpha2*(27 + 86*k2 + 7*k4)*t5 + 48*k*(-1 &
                  + k2)**3*(3 + 7*k2)*t6)*exp_minus_t + (2021760*k + 9020160*k3 + 570240*k5 + 51840*k2*(39 + 174*k2 &
                  + 11*k4)*t + 5760*k*(-11 + 112*k2 + 747*k4 + 48*k6)*t2 + 5760*k2*(-11 - 5*k2 + 225*k4 &
                  + 15*k6)*t3 + 96*k*(alpha)*(81 + 1057*k2 + 2863*k4 + 199*k6)*t4 + 12*alpha2*(49 &
                  + 1277*k2 + 3219*k4 + 255*k6)*t5 + 12*k*alpha3*(77 + 294*k2 + 29*k4)*t6 + 2*(-1 &
                  + k2)**4*(10 + 97*k2 + 13*k4)*t7 + k*alpha5*(5 &
                  + k2)*t8)*exp(-(k*t)))*sqrtk)/(45*sqrt(21.)*alpha9*t4)
 
ovlp(1) = ovlp(1) + (64*k5*((-673920*k - 3006720*k3 - 190080*k5 - 17280*k*(39 + 174*k2 + 11*k4)*t + 2880*k*(-95 &
                  - 395*k2 + 39*k4 + 3*k6)*t2 + 2880*k*(alpha)*(17 + 64*k2 + 3*k4)*t3 - 48*k*(-1 &
                  + k2)**2*(87 + 266*k2 + 7*k4)*t4 + 48*k*alpha3*(3 + 7*k2)*t5)*exp_minus_t + (673920*k &
                  + 3006720*k3 + 190080*k5 + 17280*k2*(39 + 174*k2 + 11*k4)*t + 5760*k*(-11 - 5*k2 + 225*k4 &
                  + 15*k6)*t2 + 5760*k2*(alpha)*(11 + 55*k2 + 4*k4)*t3 + 96*k*alpha2*(81 + 478*k2 &
                  + 41*k4)*t4 + 6*alpha3*(49 + 678*k2 + 73*k4)*t5 + 30*k*alpha4*(7 + k2)*t6 + (-1 &
                  + k2)**5*(5 + k2)*t7)*exp(-(k*t)))*sqrtk)/(45*sqrt(7.)*alpha9*t4)

endif

!*********************
case(8)  ! <d1|d2> 
!*********************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((79380*t + 79380*t2 + 29295*t3 + 2835*t4 - 1485*t5 - 666*t6 - 114*t7 + 3*t8 &
                  + 5*t9)*exp_minus_t)/(17010*t*sqrt(35.))
ovlp(1) = ovlp(1) - ((-26460*t - 26460*t2 - 10395*t3 - 1575*t4 + 240*t5 + 177*t6 + 43*t7 &
                  + 5*t8)*exp_minus_t)/(5670*t*sqrt(35.))
ovlp(2) = ovlp(2) + ((26460*t + 26460*t2 + 12285*t3 + 3465*t4 + 645*t5 + 78*t6 + 5*t7)*exp_minus_t)/(5670*t*sqrt(35.))

else

ovlp(0) = ovlp(0) + (128*k5*((-2177280*k - 9434880*k3 - 725760*k*(3 + 13*k2)*t + 5760*k*(-177 - 752*k2 + 33*k4)*t2 &
                  + 5760*k*(-51 - 206*k2 + 33*k4)*t3 - 480*k*(123 + 463*k2 - 147*k4 + 9*k6)*t4 + 24*k*(-1 &
                  + k2)*(333 + 1261*k2 + 79*k4 + 7*k6)*t5 - 24*k*alpha2*(27 + 86*k2 + 7*k4)*t6 + 8*k*(-1 &
                  + k2)**3*(3 + 7*k2)*t7)*exp_minus_t + (2177280*k + 9434880*k3 + 725760*k2*(3 + 13*k2)*t + 11520*k*(-6 &
                  + 61*k2 + 393*k4)*t2 + 23040*k2*(-3 - k2 + 60*k4)*t3 + 960*k*(3 - 13*k2 - 81*k4 &
                  + 315*k6)*t4 + 24*(alpha)*(7 + 139*k2 + 581*k4 + 2073*k6)*t5 + 24*k*alpha2*(17 &
                  + 126*k2 + 257*k4)*t6 + alpha3*(11 + 234*k2 + 555*k4)*t7 + k*alpha4*(7 &
                  + 33*k2)*t8 + k2*alpha5*t9)*exp(-(k*t)))*sqrtk)/(9*sqrt(35.)*alpha9*t5)
 
ovlp(1) = ovlp(1) + (-128*k5*((-483840*k - 2096640*k3 - 161280*k*(3 + 13*k2)*t + 1920*k*(-114 - 479*k2 + 33*k4)*t2 &
                  + 1920*k*(-30 - 115*k2 + 33*k4)*t3 + 2880*k*(alpha)*(3 + 11*k2)*t4 - 8*k*alpha2*(87 &
                  + 266*k2 + 7*k4)*t5 + 8*k*alpha3*(3 + 7*k2)*t6)*exp_minus_t + (483840*k + 2096640*k3 &
                  + 161280*k2*(3 + 13*k2)*t + 1920*k*(-12 + 59*k2 + 513*k4)*t2 + 1920*k2*(-12 - 25*k2 &
                  + 149*k4)*t3 + 1920*k3*(alpha)*(6 + 29*k2)*t4 + 8*alpha2*(7 + 266*k2 + 927*k4)*t5 &
                  + 8*k*alpha3*(17 + 83*k2)*t6 + alpha4*(3 + 37*k2)*t7 + k*(-1 &
                  + k2)**5*t8)*exp(-(k*t)))*sqrtk)/(3*sqrt(35.)*alpha9*t5)

ovlp(2) = ovlp(2) + (128*k5*((-120960*k - 524160*k3 - 40320*k*(3 + 13*k2)*t + 960*k*(-51 - 206*k2 + 33*k4)*t2 &
                  + 2880*k*(alpha)*(3 + 11*k2)*t3 - 720*k*alpha2*(1 + 3*k2)*t4 + 8*k*alpha3*(3 &
                  + 7*k2)*t5)*exp_minus_t + (120960*k + 524160*k3 + 40320*k2*(3 + 13*k2)*t + 3840*k*(-3 - k2 &
                  + 60*k4)*t2 + 1920*k2*(alpha)*(6 + 29*k2)*t3 + 480*k*alpha2*(3 + 17*k2)*t4 &
                  + 8*alpha3*(7 + 93*k2)*t5 + 40*k*alpha4*t6 + (-1 &
                  + k2)**5*t7)*exp(-(k*t)))*sqrtk)/(3*sqrt(35.)*alpha9*t5)

endif

end select

! ------------------------------ n1 = n2 = 4 -----------------------------------
  case(13) 
! ------------------------------------------------------------------------------

select case(llopt)

!*********************************
case(0)  ! l1 = l2 = 0 => <s1|s2>
!*********************************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((99225 + 99225*t + 47250*t2 + 14175*t3 + 3024*t4 + 504*t5 + 72*t6 + 9*t7 &
                  + t8)*exp_minus_t)/99225 
else
ovlp(0) = ovlp(0) + (32*k4*((3024 + 91584*k2 + 240864*k4 + 91584*k6 + 3024*k8 - 192*(alpha)*(7 &
                  + 169*k2 + 321*k4 + 63*k6)*t + 12*alpha2*(21 + 393*k2 + 511*k4 &
                  + 35*k6)*t2 - 8*alpha3*(3 + 42*k2 + 35*k4)*t3 + alpha4*(1 + 10*k2 &
                  + 5*k4)*t4)*exp_minus_t + (-3024 - 91584*k2 - 240864*k4 - 91584*k6 - 3024*k8 - 192*k*(-1 &
                  + k2)*(63 + 321*k2 + 169*k4 + 7*k6)*t - 12*alpha2*(35 + 511*k2 + 393*k4 &
                  + 21*k6)*t2 - 8*k*alpha3*(35 + 42*k2 + 3*k4)*t3 - alpha4*(5 + 10*k2 &
                  + k4)*t4)*exp(-(k*t)))*sqrtk)/(105*alpha9*t) 
endif

!*****************************************
case(1,3)  ! l1 = 1 and l2 = 0  => <p1|s2> 
!*****************************************

if(k.eq.1) then
ovlp(0) = ovlp(0) + (t*(33075 + 33075*t + 16065*t2 + 5040*t3 + 1152*t4 + 207*t5 + 31*t6 &
                  + 4*t7)*exp_minus_t)/(132300*sqrt(3.))
else
ovlp(0) = ovlp(0) + (-32*k4*((-47520*k - 258720*k3 - 119520*k5 - 4320*k7 - 480*k*(99 + 539*k2 &
                  + 249*k4 + 9*k6)*t + 96*k*(alpha)*(153 + 701*k2 + 259*k4 + 7*k6)*t2 - 32*k*(-1 &
                  + k2)**2*(69 + 242*k2 + 49*k4)*t3 + 2*k*alpha3*(87 + 218*k2 + 15*k4)*t4 &
                  - 2*k*alpha4*(3 + 5*k2)*t5)*exp_minus_t + (47520*k + 258720*k3 + 119520*k5 + 4320*k7 &
                  + 480*k2*(99 + 539*k2 + 249*k4 + 9*k6)*t + 48*k*(alpha)*(189 + 1293*k2 + 727*k4 &
                  + 31*k6)*t2 + 2*alpha2*(175 + 3059*k2 + 2397*k4 + 129*k6)*t3 + 8*k*(-1 &
                  + k2)**3*(35 + 42*k2 + 3*k4)*t4 + alpha4*(5 + 10*k2 &
                  + k4)*t5)*exp(-(k*t)))*sqrtk)/(35*sqrt(3.)*alpha9*t2)
endif

!***************************************
case(2,6)  ! l1 = 2 and l2 = 0  => <d1|s2> 
!***************************************

if(k.eq.1) then
ovlp(0) = ovlp(0) +  (t2*(-945 - 945*t - 216*t2 + 99*t3 + 75*t4 + 21*t5 + 4*t6)*exp_minus_t)/(79380*sqrt(5.))
else
ovlp(0) = ovlp(0) + (-32*k4*((-47520 - 807840*k2 - 419040*k4 - 15840*k6 - 1440*(33 + 561*k2 + 291*k4 &
                  + 11*k6)*t + 96*(-216 - 3501*k2 - 1021*k4 + 249*k6 + 9*k8)*t2 + 96*(alpha)*(51 &
                  + 747*k2 + 313*k4 + 9*k6)*t3 - 6*alpha2*(113 + 1381*k2 + 419*k4 + 7*k6)*t4 &
                  + 2*alpha3*(27 + 258*k2 + 35*k4)*t5 - 2*alpha4*(1 + 7*k2)*t6)*exp_minus_t &
                  + (47520 + 807840*k2 + 419040*k4 + 15840*k6 + 1440*k*(33 + 561*k2 + 291*k4 + 11*k6)*t &
                  + 48*(-63 - 918*k2 + 6092*k4 + 3702*k6 + 147*k8)*t2 + 48*k*(alpha)*(63 + 1311*k2 &
                  + 829*k4 + 37*k6)*t3 + 30*alpha2*(7 + 203*k2 + 165*k4 + 9*k6)*t4 + 8*k*(-1 &
                  + k2)**3*(35 + 42*k2 + 3*k4)*t5 + alpha4*(5 + 10*k2 &
                  + k4)*t6)*exp(-(k*t)))*sqrtk)/(21*sqrt(5.)*alpha9*t3)
endif

!***************************************
case(4)  ! l1 = l2 = 1  => <p1|p2> 
!***************************************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((66150*t + 66150*t2 + 27405*t3 + 5355*t4 + 90*t5 - 225*t6 - 67*t7 - 13*t8 &
                  - 2*t9)*exp_minus_t)/(66150*t)
ovlp(1) = ovlp(1) -((-66150*t - 66150*t2 - 31185*t3 - 9135*t4 - 1845*t5 - 270*t6 - 29*t7 &
                  - 2*t8)*exp_minus_t)/(66150*t) 
else
ovlp(0) = ovlp(0) + (-64*k4*((63360*k + 303360*k3 + 63360*k5 + 1920*k*(33 + 158*k2 + 33*k4)*t + 960*k*(33 &
                  + 158*k2 + 33*k4)*t2 - 40*k*(alpha)*(201 + 877*k2 + 259*k4 + 7*k6)*t3 + 8*k*(-1 &
                  + k2)**2*(141 + 488*k2 + 91*k4)*t4 - k*alpha3*(87 + 218*k2 + 15*k4)*t5 &
                  + k*alpha4*(3 + 5*k2)*t6)*exp_minus_t + (-63360*k - 303360*k3 - 63360*k5 &
                  - 1920*k2*(33 + 158*k2 + 33*k4)*t - 960*k3*(33 + 158*k2 + 33*k4)*t2 - 40*(-1 &
                  + k2)*(7 + 259*k2 + 877*k4 + 201*k6)*t3 - 8*k*alpha2*(91 + 488*k2 &
                  + 141*k4)*t4 - alpha3*(15 + 218*k2 + 87*k4)*t5 - k*alpha4*(5 &
                  + 3*k2)*t6)*exp(-(k*t)))*sqrtk)/(35*alpha9*t3)

ovlp(1) = ovlp(1) + (64*k4*((31680*k + 151680*k3 + 31680*k5 + 960*k*(33 + 158*k2 + 33*k4)*t - 240*k*(-1 &
                  + k2)*(39 + 158*k2 + 27*k4)*t2 + 40*k*alpha2*(33 + 104*k2 + 7*k4)*t3 &
                  - 32*k*alpha3*(3 + 7*k2)*t4 + k*alpha4*(3 + 5*k2)*t5)*exp_minus_t + (-31680*k &
                  - 151680*k3 - 31680*k5 - 960*k2*(33 + 158*k2 + 33*k4)*t - 240*k*(alpha)*(27 &
                  + 158*k2 + 39*k4)*t2 - 40*alpha2*(7 + 104*k2 + 33*k4)*t3 - 32*k*(-1 &
                  + k2)**3*(7 + 3*k2)*t4 - alpha4*(5 + 3*k2)*t5)*exp(-(k*t)))*sqrtk)/(35*(-1 &
                  + k2)**9*t3)

endif

!***************************************
case(5,7)  ! l1 = 2 and l2 = 1  => <d1|p2> 
!***************************************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((13230*t2 + 13230*t3 + 5805*t4 + 1395*t5 + 162*t6 - 9*t7 - 8*t8 &
                  - 2*t9)*exp_minus_t)/(13230*t*sqrt(15.))
ovlp(1) = ovlp(1) +  ((6615*t2 + 6615*t3 + 3105*t4 + 900*t5 + 177*t6 + 24*t7 &
                  + 2*t8)*exp_minus_t)/(13230*t*sqrt(5.))
else
ovlp(0) = ovlp(0) + (64*k4*((-95040 - 1503360*k2 - 336960*k4 - 8640*(11 + 174*k2 + 39*k4)*t + 240*(-189 &
                  - 2949*k2 - 479*k4 + 33*k6)*t2 + 240*(-57 - 861*k2 - 11*k4 + 33*k6)*t3 + 24*(-1 &
                  + k2)*(113 + 1601*k2 + 499*k4 + 27*k6)*t4 - 4*alpha2*(87 + 1053*k2 + 293*k4 &
                  + 7*k6)*t5 + alpha3*(27 + 258*k2 + 35*k4)*t6 - alpha4*(1 &
                  + 7*k2)*t7)*exp_minus_t + (95040 + 1503360*k2 + 336960*k4 + 8640*k*(11 + 174*k2 + 39*k4)*t &
                  + 240*(-9 + 15*k2 + 2909*k4 + 669*k6)*t2 + 240*k*(-9 - 117*k2 + 821*k4 + 201*k6)*t3 &
                  + 24*(alpha)*(7 + 259*k2 + 1581*k4 + 393*k6)*t4 + 8*k*alpha2*(77 + 496*k2 &
                  + 147*k4)*t5 + alpha3*(15 + 218*k2 + 87*k4)*t6 + k*alpha4*(5 &
                  + 3*k2)*t7)*exp(-(k*t)))*sqrtk)/(7*sqrt(15.)*alpha9*t4)

ovlp(1) = ovlp(1) + (64*k4*((31680 + 501120*k2 + 112320*k4 + 2880*(11 + 174*k2 + 39*k4)*t - 240*(-57 &
                  - 861*k2 - 11*k4 + 33*k6)*t2 - 240*(alpha)*(13 + 178*k2 + 33*k4)*t3 + 24*(-1 &
                  + k2)**2*(17 + 196*k2 + 27*k4)*t4 - 2*alpha3*(15 + 138*k2 + 7*k4)*t5 + (-1 &
                  + k2)**4*(1 + 7*k2)*t6)*exp_minus_t + (-31680 - 501120*k2 - 112320*k4 - 2880*k*(11 &
                  + 174*k2 + 39*k4)*t - 240*(-9 - 117*k2 + 821*k4 + 201*k6)*t2 - 240*k*(alpha)*(9 &
                  + 170*k2 + 45*k4)*t3 - 24*alpha2*(7 + 176*k2 + 57*k4)*t4 - 32*k*(-1 &
                  + k2)**3*(7 + 3*k2)*t5 - alpha4*(5 &
                  + 3*k2)*t6)*exp(-(k*t)))*sqrtk)/(7*sqrt(5.)*alpha9*t4)
endif

!***************************************
case(8)  ! l1 = l2 = 2  => <d1|d2> 
!***************************************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((39690*t + 39690*t2 + 13905*t3 + 675*t4 - 972*t5 - 315*t6 - 33*t7 &
                  + 3*t8 + 2*t9)*exp_minus_t)/(39690*t)
ovlp(1) = ovlp(1) - ((-13230*t - 13230*t2 - 4995*t3 - 585*t4 + 195*t5 + 96*t6 + 19*t7 &
                  + 2*t8)*exp_minus_t)/(13230*t)
ovlp(2) = ovlp(2) + ((13230*t + 13230*t2 + 6075*t3 + 1665*t4 + 294*t5 + 33*t6 &
                  + 2*t7)*exp_minus_t)/(13230*t)
else
ovlp(0) = ovlp(0) + (-64*k4*((-673920 - 10264320*k2 - 673920*k4 - 51840*(13 + 198*k2 + 13*k4)*t &
                  + 1440*(-223 - 3351*k2 - 21*k4 + 11*k6)*t2 + 1440*(-67 - 975*k2 + 135*k4 &
                  + 11*k6)*t3 - 48*(429 + 5964*k2 - 2006*k4 + 84*k6 + 9*k8)*t4 + 120*(-1 &
                  + k2)*(27 + 363*k2 + 49*k4 + 9*k6)*t5 - 6*alpha2*(61 + 725*k2 &
                  + 167*k4 + 7*k6)*t6 + alpha3*(27 + 258*k2 + 35*k4)*t7 - alpha4*(1 &
                  + 7*k2)*t8)*exp_minus_t + (673920 + 10264320*k2 + 673920*k4 + 51840*k*(13 + 198*k2 &
                  + 13*k4)*t + 1440*(-11 + 21*k2 + 3351*k4 + 223*k6)*t2 + 1440*k*(-11 - 135*k2 &
                  + 975*k4 + 67*k6)*t3 + 48*(9 + 84*k2 - 2006*k4 + 5964*k6 + 429*k8)*t4 + 120*k*(-1 &
                  + k2)*(9 + 49*k2 + 363*k4 + 27*k6)*t5 + 6*alpha2*(7 + 167*k2 + 725*k4 &
                  + 61*k6)*t6 + k*alpha3*(35 + 258*k2 + 27*k4)*t7 + k2*alpha4*(7 &
                  + k2)*t8)*exp(-(k*t)))*sqrtk)/(21*alpha9*t5)

ovlp(1) = ovlp(1) + (64*k4*((-149760 - 2280960*k2 - 149760*k4 - 11520*(13 + 198*k2 + 13*k4)*t + 480*(-145 &
                  - 2163*k2 + 57*k4 + 11*k6)*t2 + 480*(-41 - 579*k2 + 161*k4 + 11*k6)*t3 + 240*(-1 &
                  + k2)*(15 + 198*k2 + 11*k4)*t4 - 8*alpha2*(53 + 604*k2 + 63*k4)*t5 + 2*(-1 &
                  + k2)**3*(15 + 138*k2 + 7*k4)*t6 - alpha4*(1 + 7*k2)*t7)*exp_minus_t + (149760 &
                  + 2280960*k2 + 149760*k4 + 11520*k*(13 + 198*k2 + 13*k4)*t + 480*(-11 - 57*k2 + 2163*k4 &
                  + 145*k6)*t2 + 480*k*(-11 - 161*k2 + 579*k4 + 41*k6)*t3 + 240*k2*(alpha)*(11 &
                  + 198*k2 + 15*k4)*t4 + 8*k*alpha2*(63 + 604*k2 + 53*k4)*t5 + 2*alpha3*(7 &
                  + 138*k2 + 15*k4)*t6 + k*alpha4*(7 + k2)*t7)*exp(-(k*t)))*sqrtk)/(7*(-1 &
                  + k2)**9*t5)

ovlp(2) = ovlp(2) + (64*k4*((37440 + 570240*k2 + 37440*k4 + 2880*(13 + 198*k2 + 13*k4)*t - 240*(-67 &
                  - 975*k2 + 135*k4 + 11*k6)*t2 - 240*(alpha)*(15 + 198*k2 + 11*k4)*t3 + 24*(-1 &
                  + k2)**2*(19 + 212*k2 + 9*k4)*t4 - 32*alpha3*(1 + 9*k2)*t5 + alpha4*(1 &
                  + 7*k2)*t6)*exp_minus_t + (-37440 - 570240*k2 - 37440*k4 - 2880*k*(13 + 198*k2 + 13*k4)*t &
                  - 240*(-11 - 135*k2 + 975*k4 + 67*k6)*t2 - 240*k*(alpha)*(11 + 198*k2 + 15*k4)*t3 &
                  - 24*alpha2*(9 + 212*k2 + 19*k4)*t4 - 32*k*alpha3*(9 + k2)*t5 - (-1 &
                  + k2)**4*(7 + k2)*t6)*exp(-(k*t)))*sqrtk)/(7*alpha9*t5)

endif

end select  

! ------------------------------ n1 = 5 and n2 = 4 -----------------------------
  case(14) 
! ------------------------------------------------------------------------------

select case(llopt)

!*********************************
case(0)   ! l1 = l2 = 0 => <s1|s2>
!*********************************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((893025 + 893025*t + 429975*t2 + 132300*t3 + 29295*t4 + 5040*t5 + 720*t6 + 90*t7 &
                  + 10*t8 + t9)*exp_minus_t*sqrt(2./5.))/595350.
else

ovlp(0) = ovlp(0) + (-32*k5*((-237600*k - 2428800*k3 - 3921600*k5 - 1123200*k7 - 30240*k9 + 1920*k*(alpha)*(45 &
                  + 365*k2 + 423*k4 + 63*k6)*t - 480*k*alpha2*(27 + 169*k2 + 133*k4 + 7*k6)*t2 &
                  + 320*k*alpha3*(3 + 14*k2 + 7*k4)*t3 - 10*k*alpha4*(3 + k2)*(1 + 3*k2)*t4)*exp_minus_t &
                  + (237600*k + 2428800*k3 + 3921600*k5 + 1123200*k7 + 30240*k9 + 240*(alpha)*(63 + 1908*k2 &
                  + 5018*k4 + 1908*k6 + 63*k8)*t + 480*k*alpha2*(63 + 321*k2 + 169*k4 + 7*k6)*t2 + 20*(-1 &
                  + k2)**3*(35 + 511*k2 + 393*k4 + 21*k6)*t3 + 10*k*alpha4*(35 + 42*k2 + 3*k4)*t4 + (-1 &
                  + k2)**5*(5 + 10*k2 + k4)*t5)*exp(-(k*t)))*sqrt(2./5.)*sqrtk)/(315*alpha10*t)
endif

!*******************
case(1)   ! <s1|p2>
!*******************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((-66150*t2 - 66150*t3 - 34965*t4 - 12915*t5 - 3600*t6 - 765*t7 - 125*t8 - 17*t9 &
                  - 2*t10)*exp_minus_t*sqrt(2./15.))/(396900*t)
else
ovlp(0) = ovlp(0) + (64*k5*((-158400*k - 1473600*k3 - 1953600*k5 - 285120*k7 - 960*k*(165 + 1535*k2 + 2035*k4 &
                  + 297*k6)*t + 240*k*(alpha)*(195 + 1535*k2 + 1665*k4 + 189*k6)*t2 - 40*k*alpha2*(165 &
                  + 1025*k2 + 791*k4 + 35*k6)*t3 + 160*k*alpha3*(3 + 14*k2 + 7*k4)*t4 - 5*k*(-1 &
                  + k2)**4*(3 + k2)*(1 + 3*k2)*t5)*exp_minus_t + (158400*k + 1473600*k3 + 1953600*k5 + 285120*k7 &
                  + 960*k2*(165 + 1535*k2 + 2035*k4 + 297*k6)*t + 1200*k*(alpha)*(27 + 307*k2 + 481*k4 &
                  + 81*k6)*t2 + 160*alpha2*(7 + 214*k2 + 511*k4 + 108*k6)*t3 + 200*k*alpha3*(7 &
                  + 32*k2 + 9*k4)*t4 + 4*alpha4*(5 + 68*k2 + 27*k4)*t5 + k*alpha5*(5 &
                  + 3*k2)*t6)*exp(-(k*t)))*sqrt(2./15.)*sqrtk)/(105*alpha10*t2)
endif

!*******************
case(2)   ! <s1|d2>
!*******************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((-13230*t3 - 13230*t4 - 5535*t5 - 1125*t6 - 45*t7 + 36*t8 + 11*t9 &
                  + 2*t10)*exp_minus_t*sqrt(2.))/(1190700*t)
else

ovlp(0) = ovlp(0) + (-64*k5*((-561600*k - 4968000*k3 - 5797440*k5 - 285120*k7 - 8640*k*(65 + 575*k2 + 671*k4 &
                  + 33*k6)*t + 240*k*(-1005 - 8360*k2 - 7978*k4 + 1152*k6 + 63*k8)*t2 + 240*k*(alpha)*(225 &
                  + 1685*k2 + 1611*k4 + 63*k6)*t3 - 120*k*alpha2*(57 + 349*k2 + 259*k4 + 7*k6)*t4 &
                  + 160*k*alpha3*(3 + 14*k2 + 7*k4)*t5 - 5*k*alpha4*(3 + k2)*(1 + 3*k2)*t6)*exp_minus_t &
                  + (561600*k + 4968000*k3 + 5797440*k5 + 285120*k7 + 8640*k2*(65 + 575*k2 + 671*k4 + 33*k6)*t &
                  + 240*k*(-165 - 820*k2 + 6250*k4 + 10332*k6 + 531*k8)*t2 + 1200*k2*(alpha)*(33 + 353*k2 &
                  + 483*k4 + 27*k6)*t3 + 120*k*alpha2*(27 + 389*k2 + 661*k4 + 43*k6)*t4 + 8*(-1 &
                  + k2)**3*(7 + 302*k2 + 825*k4 + 66*k6)*t5 + k*alpha4*(49 + 318*k2 + 33*k4)*t6 &
                  + k2*alpha5*(7 + k2)*t7)*exp(-(k*t)))*sqrt(2.)*sqrtk)/(315*alpha10*t3)
endif

!*********************
case(3)   ! <p1|s2>
!*********************

if(k.eq.1) then
ovlp(0) = ovlp(0) + (t*(99225 + 99225*t + 47250*t2 + 14175*t3 + 3024*t4 + 504*t5 + 72*t6 + 9*t7 &
                  + t8)*exp_minus_t*sqrt(2./15.))/198450
else

ovlp(0) = ovlp(0) + (-32*k5*((-47520 - 1584000*k2 - 4478400*k4 - 1584000*k6 - 47520*k8 - 1440*(33 + 1100*k2 &
                  + 3110*k4 + 1100*k6 + 33*k8)*t + 96*(alpha)*(153 + 4398*k2 + 10408*k4 + 2898*k6 &
                  + 63*k8)*t2 - 96*alpha2*(23 + 541*k2 + 969*k4 + 147*k6)*t3 + 6*alpha3*(29 &
                  + 537*k2 + 679*k4 + 35*k6)*t4 - 2*alpha4*(3 + 42*k2 + 35*k4)*t5)*exp_minus_t + (47520 &
                  + 1584000*k2 + 4478400*k4 + 1584000*k6 + 47520*k8 + 1440*k*(33 + 1100*k2 + 3110*k4 + 1100*k6 &
                  + 33*k8)*t + 48*(alpha)*(189 + 7704*k2 + 25834*k4 + 10704*k6 + 369*k8)*t2 + 96*k*(-1 &
                  + k2)**2*(273 + 1611*k2 + 879*k4 + 37*k6)*t3 + 6*alpha3*(105 + 1701*k2 + 1323*k4 &
                  + 71*k6)*t4 + 10*k*alpha4*(35 + 42*k2 + 3*k4)*t5 + alpha5*(5 + 10*k2 &
                  + k4)*t6)*exp(-(k*t)))*sqrt(2./15.)*sqrtk)/(105*alpha10*t2)
endif

!*********************************
case(4)   ! l1 = l2 = 1 => <p1|p2>
!*********************************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((595350*t + 595350*t2 + 257985*t3 + 59535*t4 + 5805*t5 - 810*t6 - 426*t7 - 93*t8 &
                  - 15*t9 - 2*t10)*exp_minus_t*sqrt(2./5.))/(396900*t)
ovlp(1) = ovlp(1) - ((-595350*t - 595350*t2 - 284445*t3 - 85995*t4 - 18360*t5 - 2925*t6 - 357*t7 - 33*t8 &
                  - 2*t9)*exp_minus_t*sqrt(2./5.))/(396900*t) 
else

ovlp(0) = ovlp(0) + (64*k5*((-63360 - 1987200*k2 - 4867200*k4 - 823680*k6 - 5760*(11 + 345*k2 + 845*k4 &
                  + 143*k6)*t - 2880*(11 + 345*k2 + 845*k4 + 143*k6)*t2 + 120*(alpha)*(67 + 1882*k2 &
                  + 4232*k4 + 966*k6 + 21*k8)*t3 - 24*alpha2*(47 + 1099*k2 + 1941*k4 + 273*k6)*t4 &
                  + 3*alpha3*(29 + 537*k2 + 679*k4 + 35*k6)*t5 - alpha4*(3 + 42*k2 &
                  + 35*k4)*t6)*exp_minus_t + (63360 + 1987200*k2 + 4867200*k4 + 823680*k6 + 5760*k*(11 + 345*k2 &
                  + 845*k4 + 143*k6)*t + 2880*k2*(11 + 345*k2 + 845*k4 + 143*k6)*t2 + 240*k*(alpha)*(105 &
                  + 1321*k2 + 2587*k4 + 467*k6)*t3 + 48*alpha2*(21 + 657*k2 + 1743*k4 + 379*k6)*t4 &
                  + 96*k*alpha3*(14 + 67*k2 + 19*k4)*t5 + 4*alpha4*(5 + 68*k2 + 27*k4)*t6 + k*(-1 &
                  + k2)**5*(5 + 3*k2)*t7)*exp(-(k*t)))*sqrt(2./5.)*sqrtk)/(105*alpha10*t3)

ovlp(1) = ovlp(1) + (-64*k5*((-31680 - 993600*k2 - 2433600*k4 - 411840*k6 - 2880*(11 + 345*k2 + 845*k4 + 143*k6)*t &
                  + 240*(alpha)*(39 + 1059*k2 + 2189*k4 + 297*k6)*t2 - 120*alpha2*(11 + 247*k2 + 393*k4 &
                  + 21*k6)*t3 + 96*alpha3*(1 + 18*k2 + 21*k4)*t4 - alpha4*(3 + 42*k2 &
                  + 35*k4)*t5)*exp_minus_t + (31680 + 993600*k2 + 2433600*k4 + 411840*k6 + 2880*k*(11 + 345*k2 &
                  + 845*k4 + 143*k6)*t + 240*(alpha)*(27 + 1011*k2 + 2881*k4 + 561*k6)*t2 + 480*k*(-1 &
                  + k2)**2*(39 + 194*k2 + 47*k4)*t3 + 24*alpha3*(21 + 288*k2 + 91*k4)*t4 + 40*k*(-1 &
                  + k2)**4*(7 + 3*k2)*t5 + alpha5*(5 + 3*k2)*t6)*exp(-(k*t)))*sqrt(2./5.)*sqrtk)/(105*(-1 &
                  + k2)**10*t3)
endif

!***************************************
case(5)   ! l1 = 1 and l2 = 2 => <p1|d2>
!***************************************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((-79380*t2 - 79380*t3 - 37530*t4 - 11070*t5 - 2151*t6 - 225*t7 + 12*t8 + 9*t9 &
                  + 2*t10)*exp_minus_t*sqrt(2./3.))/(396900*t)
ovlp(1) = ovlp(1) - ((39690*t2 + 39690*t3 + 19305*t4 + 6075*t5 + 1371*t6 + 228*t7 + 27*t8 &
                  + 2*t9)*exp_minus_t*sqrt(2.))/(396900*t)
else

ovlp(0) = ovlp(0) + (-64*k5*((-336960 - 10238400*k2 - 23025600*k4 - 1235520*k6 - 8640*(39 + 1185*k2 + 2665*k4 &
                  + 143*k6)*t + 720*(-223 - 6700*k2 - 14570*k4 - 44*k6 + 33*k8)*t2 + 720*(-67 - 1960*k2 &
                  - 3910*k4 + 528*k6 + 33*k8)*t3 + 24*(alpha)*(393 + 10638*k2 + 21848*k4 + 2898*k6 &
                  + 63*k8)*t4 - 24*alpha2*(49 + 1133*k2 + 1947*k4 + 231*k6)*t5 + 3*alpha3*(29 &
                  + 537*k2 + 679*k4 + 35*k6)*t6 - alpha4*(3 + 42*k2 + 35*k4)*t7)*exp_minus_t + (336960 &
                  + 10238400*k2 + 23025600*k4 + 1235520*k6 + 8640*k*(39 + 1185*k2 + 2665*k4 + 143*k6)*t + 720*(-11 &
                  - 176*k2 + 5690*k4 + 15176*k6 + 825*k8)*t2 + 720*k*(-11 - 332*k2 + 950*k4 + 4516*k6 &
                  + 253*k8)*t3 + 24*(alpha)*(27 + 1572*k2 + 14602*k4 + 27012*k6 + 1587*k8)*t4 + 24*k*(-1 &
                  + k2)**2*(141 + 1787*k2 + 3443*k4 + 229*k6)*t5 + 3*alpha3*(21 + 777*k2 + 2223*k4 &
                  + 179*k6)*t6 + k*alpha4*(49 + 318*k2 + 33*k4)*t7 + k2*alpha5*(7 &
                  + k2)*t8)*exp(-(k*t)))*sqrt(2./3.)*sqrtk)/(105.*alpha10*t4)

ovlp(1) = ovlp(1) + (64*k5*((-112320 - 3412800*k2 - 7675200*k4 - 411840*k6 - 2880*(39 + 1185*k2 + 2665*k4 &
                  + 143*k6)*t + 720*(-67 - 1960*k2 - 3910*k4 + 528*k6 + 33*k8)*t2 + 240*(alpha)*(45 &
                  + 1185*k2 + 2255*k4 + 99*k6)*t3 - 72*alpha2*(19 + 423*k2 + 657*k4 + 21*k6)*t4 &
                  + 96*alpha3*(1 + 18*k2 + 21*k4)*t5 - alpha4*(3 + 42*k2 + 35*k4)*t6)*exp_minus_t &
                  + (112320 + 3412800*k2 + 7675200*k4 + 411840*k6 + 2880*k*(39 + 1185*k2 + 2665*k4 + 143*k6)*t &
                  + 720*(-11 - 332*k2 + 950*k4 + 4516*k6 + 253*k8)*t2 + 240*k*(alpha)*(33 + 1185*k2 &
                  + 3075*k4 + 187*k6)*t3 + 24*alpha2*(27 + 1269*k2 + 4021*k4 + 283*k6)*t4 + 24*k*(-1 &
                  + k2)**3*(57 + 316*k2 + 27*k4)*t5 + alpha4*(21 + 342*k2 + 37*k4)*t6 + k*(-1 &
                  + k2)**5*(7 + k2)*t7)*exp(-(k*t)))*sqrt(2.)*sqrtk)/(105.*alpha10*t4)

endif

!********************
case(6)   ! <d1|s2>
!********************

if(k.eq.1) then
ovlp(0) = ovlp(0) + (t2*(6615 + 6615*t + 3510*t2 + 1305*t3 + 369*t4 + 81*t5 + 14*t6 &
                  + 2*t7)*exp_minus_t*sqrt(2.))/1190700
else

ovlp(0) = ovlp(0) + (-32*k5*((-2471040*k - 14601600*k3 - 5961600*k5 - 190080*k7 - 17280*k*(143 + 845*k2 + 345*k4 &
                  + 11*k6)*t + 960*k*(-1089 - 6010*k2 - 1280*k4 + 306*k6 + 9*k8)*t2 + 960*k*(alpha)*(231 &
                  + 1171*k2 + 381*k4 + 9*k6)*t3 - 48*k*alpha2*(543 + 2281*k2 + 529*k4 + 7*k6)*t4 &
                  + 80*k*alpha3*(21 + 68*k2 + 7*k4)*t5 - 16*k*alpha4*(3 + 7*k2)*t6)*exp_minus_t &
                  + (2471040*k + 14601600*k3 + 5961600*k5 + 190080*k7 + 17280*k2*(143 + 845*k2 + 345*k4 + 11*k6)*t &
                  + 1920*k*(-99 - 154*k2 + 2890*k4 + 1350*k6 + 45*k8)*t2 + 1920*k2*(alpha)*(99 + 682*k2 &
                  + 327*k4 + 12*k6)*t3 + 96*k*alpha2*(189 + 1623*k2 + 947*k4 + 41*k6)*t4 + 2*(-1 &
                  + k2)**3*(245 + 5089*k2 + 4047*k4 + 219*k6)*t5 + 10*k*alpha4*(35 + 42*k2 + 3*k4)*t6 &
                  + alpha5*(5 + 10*k2 + k4)*t7)*exp(-(k*t)))*sqrt(2.)*sqrtk)/(315.*alpha10*t3)

endif

!*********************
case(7)  ! <d1|p2> 
!*********************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((132300*t2 + 132300*t3 + 58590*t4 + 14490*t5 + 1935*t6 + 45*t7 - 38*t8 - 11*t9 &
                  - 2*t10)*exp_minus_t*sqrt(2./3.))/(396900*t)
ovlp(1) = ovlp(1) + ((66150*t2 + 66150*t3 + 31185*t4 + 9135*t5 + 1845*t6 + 270*t7 + 29*t8 &
                  + 2*t9)*exp_minus_t*sqrt(2.))/(396900*t)
else

ovlp(0) = ovlp(0) + (64*k5*((-4717440*k - 25401600*k3 - 4717440*k5 - 362880*k*(13 + 70*k2 + 13*k4)*t + 960*k*(-2325 &
               - 12275*k2 - 1627*k4 + 99*k6)*t2 + 960*k*(-687 - 3455*k2 + 11*k4 + 99*k6)*t3 + 240*k*(-1 &
               + k2)*(501 + 2441*k2 + 615*k4 + 27*k6)*t4 - 16*k*alpha2*(831 + 3452*k2 + 743*k4 & 
               + 14*k6)*t5 + 40*k*alpha3*(21 + 68*k2 + 7*k4)*t6 - 8*k*alpha4*(3 &
               + 7*k2)*t7)*exp_minus_t + (4717440*k + 25401600*k3 + 4717440*k5 + 362880*k2*(13 + 70*k2 + 13*k4)*t &
               + 1920*k*(-66 + 751*k2 + 6200*k4 + 1179*k6)*t2 + 3840*k2*(-33 - 34*k2 + 895*k4 + 180*k6)*t3 &
               + 480*k*(alpha)*(27 + 439*k2 + 1477*k4 + 297*k6)*t4 + 16*alpha2*(49 + 1633*k2 + 5467*k4 &
               + 1251*k6)*t5 + 16*k*alpha3*(77 + 406*k2 + 117*k4)*t6 + 4*alpha4*(5 + 68*k2 &
               + 27*k4)*t7 + k*alpha5*(5 + 3*k2)*t8)*exp(-(k*t)))*sqrt(2./3.)*sqrtk)/(105*alpha10*t4)
 
ovlp(1) = ovlp(1) + (-64*k5*((-1572480*k - 8467200*k3 - 1572480*k5 - 120960*k*(13 + 70*k2 + 13*k4)*t + 960*k*(-687 &
                  - 3455*k2 + 11*k4 + 99*k6)*t2 + 960*k*(alpha)*(141 + 656*k2 + 99*k4)*t3 - 720*k*(-1 &
                  + k2)**2*(21 + 82*k2 + 9*k4)*t4 + 16*k*alpha3*(57 + 176*k2 + 7*k4)*t5 - 8*k*(-1 &
                  + k2)**4*(3 + 7*k2)*t6)*exp_minus_t + (1572480*k + 8467200*k3 + 1572480*k5 + 120960*k2*(13 &
                  + 70*k2 + 13*k4)*t + 3840*k*(-33 - 34*k2 + 895*k4 + 180*k6)*t2 + 1920*k2*(alpha)*(66 &
                  + 407*k2 + 87*k4)*t3 + 480*k*alpha2*(27 + 202*k2 + 51*k4)*t4 + 8*alpha3*(49 &
                  + 872*k2 + 279*k4)*t5 + 40*k*alpha4*(7 + 3*k2)*t6 + alpha5*(5 &
                  + 3*k2)*t7)*exp(-(k*t)))*sqrt(2.)*sqrtk)/(105*alpha10*t4)

endif

!*********************
case(8)  ! <d1|d2> 
!*********************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((357210*t + 357210*t2 + 137025*t3 + 17955*t4 - 4185*t5 - 2232*t6 - 414*t7 - 27*t8 &
                  + 5*t9 + 2*t10)*exp_minus_t*sqrt(2./5.))/(238140*t)
ovlp(1) = ovlp(1) - ((-119070*t - 119070*t2 - 48195*t3 - 8505*t4 + 390*t5 + 579*t6 + 153*t7 + 23*t8 &
                  + 2*t9)*exp_minus_t*sqrt(2./5.))/(79380*t) 
ovlp(2) = ovlp(2) + ((119070*t + 119070*t2 + 55755*t3 + 16065*t4 + 3135*t5 + 426*t6 + 39*t7 &
                  + 2*t8)*exp_minus_t*sqrt(2./5.))/(79380*t)

else

ovlp(0) = ovlp(0) + (-64*k5*((-32659200*k - 166924800*k3 - 9434880*k5 - 725760*k*(45 + 230*k2 + 13*k4)*t &
                  + 17280*k*(-893 - 4475*k2 - 19*k4 + 11*k6)*t2 + 17280*k*(-263 - 1255*k2 + 163*k4 + 11*k6)*t3 &
                  - 480*k*(1965 + 8740*k2 - 2758*k4 + 108*k6 + 9*k8)*t4 + 240*k*(alpha)*(579 + 2639*k2 &
                  + 321*k4 + 45*k6)*t5 - 48*k*alpha2*(288 + 1171*k2 + 214*k4 + 7*k6)*t6 + 40*k*(-1 &
                  + k2)**3*(21 + 68*k2 + 7*k4)*t7 - 8*k*alpha4*(3 + 7*k2)*t8)*exp_minus_t + (32659200*k &
                  + 166924800*k3 + 9434880*k5 + 725760*k2*(45 + 230*k2 + 13*k4)*t + 34560*k*(-26 + 295*k2 &
                  + 2288*k4 + 131*k6)*t2 + 69120*k2*(-13 - 10*k2 + 339*k4 + 20*k6)*t3 + 960*k*(33 - 148*k2 &
                  - 1280*k4 + 5112*k6 + 315*k8)*t4 + 24*(alpha)*(63 + 1578*k2 + 8248*k4 + 32838*k6 &
                  + 2073*k8)*t5 + 24*k*alpha2*(153 + 1471*k2 + 3719*k4 + 257*k6)*t6 + alpha3*(77 &
                  + 2161*k2 + 6807*k4 + 555*k6)*t7 + k*alpha4*(49 + 318*k2 + 33*k4)*t8 + k2*(-1 &
                  + k2)**5*(7 + k2)*t9)*exp(-(k*t)))*sqrt(2./5.)*sqrtk)/(63*alpha10*t5) 
 
ovlp(1) = ovlp(1) + (64*k5*((-7257600*k - 37094400*k3 - 2096640*k5 - 161280*k*(45 + 230*k2 + 13*k4)*t + 5760*k*(-578 &
                  - 2865*k2 + 72*k4 + 11*k6)*t2 + 1920*k*(-474 - 2155*k2 + 580*k4 + 33*k6)*t3 + 960*k*(-1 &
                  + k2)*(159 + 704*k2 + 33*k4)*t4 - 240*k*alpha2*(65 + 250*k2 + 21*k4)*t5 + 16*k*(-1 &
                  + k2)**3*(57 + 176*k2 + 7*k4)*t6 - 8*k*alpha4*(3 + 7*k2)*t7)*exp_minus_t + (7257600*k &
                  + 37094400*k3 + 2096640*k5 + 161280*k2*(45 + 230*k2 + 13*k4)*t + 5760*k*(-52 + 275*k2 + 2966*k4 &
                  + 171*k6)*t2 + 1920*k2*(-156 - 435*k2 + 2458*k4 + 149*k6)*t3 + 1920*k3*(alpha)*(78 &
                  + 453*k2 + 29*k4)*t4 + 24*alpha2*(21 + 987*k2 + 4283*k4 + 309*k6)*t5 + 8*k*(-1 &
                  + k2)**3*(153 + 964*k2 + 83*k4)*t6 + alpha4*(21 + 342*k2 + 37*k4)*t7 + k*(-1 &
                  + k2)**5*(7 + k2)*t8)*exp(-(k*t)))*sqrt(2./5.)*sqrtk)/(21*alpha10*t5)

ovlp(2) = ovlp(2) + (-64*k5*((-1814400*k - 9273600*k3 - 524160*k5 - 40320*k*(45 + 230*k2 + 13*k4)*t + 2880*k*(-263 &
                  - 1255*k2 + 163*k4 + 11*k6)*t2 + 960*k*(alpha)*(159 + 704*k2 + 33*k4)*t3 - 720*k*(-1 &
                  + k2)**2*(23 + 86*k2 + 3*k4)*t4 + 960*k*alpha3*(1 + 3*k2)*t5 - 8*k*alpha4*(3 &
                  + 7*k2)*t6)*exp_minus_t + (1814400*k + 9273600*k3 + 524160*k5 + 40320*k2*(45 + 230*k2 + 13*k4)*t &
                  + 11520*k*(-13 - 10*k2 + 339*k4 + 20*k6)*t2 + 1920*k2*(alpha)*(78 + 453*k2 + 29*k4)*t3 &
                  + 480*k*alpha2*(33 + 230*k2 + 17*k4)*t4 + 24*alpha3*(21 + 348*k2 + 31*k4)*t5 &
                  + 40*k*alpha4*(9 + k2)*t6 + alpha5*(7 &
                  + k2)*t7)*exp(-(k*t)))*sqrt(2./5.)*sqrtk)/(21*alpha10*t5)

endif

end select

! ------------------------------ n1 = 5 and n2 = 5 -----------------------------
  case(15) 
! ------------------------------------------------------------------------------

select case(llopt)

!*********************************
case(0)   ! l1 = l2 = 0 => <s1|s2>
!*********************************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((9823275 + 9823275*t + 4729725*t2 + 1455300*t3 + 322245*t4 + 55440*t5 + 7920*t6 + 990*t7 &
                  + 110*t8 + 11*t9 + t10)*exp_minus_t)/9823275.
else
ovlp(0) = ovlp(0) + (128*k5*((-285120*k - 3590400*k3 - 7731840*k5 - 3590400*k7 - 285120*k9 + 240*k*(alpha)*(495 &
                  + 5060*k2 + 8170*k4 + 2340*k6 + 63*k8)*t - 480*k*alpha2*(45 + 365*k2 + 423*k4 &
                  + 63*k6)*t2 + 80*k*alpha3*(27 + 169*k2 + 133*k4 + 7*k6)*t3 - 40*k*alpha4*(3 &
                  + 14*k2 + 7*k4)*t4 + k*alpha5*(3 + k2)*(1 + 3*k2)*t5)*exp_minus_t + (285120*k + 3590400*k3 &
                  + 7731840*k5 + 3590400*k7 + 285120*k9 + 240*(alpha)*(63 + 2340*k2 + 8170*k4 + 5060*k6 &
                  + 495*k8)*t + 480*k*alpha2*(63 + 423*k2 + 365*k4 + 45*k6)*t2 + 80*alpha3*(7 &
                  + 133*k2 + 169*k4 + 27*k6)*t3 + 40*k*alpha4*(7 + 14*k2 + 3*k4)*t4 + alpha5*(3 &
                  + k2)*(1 + 3*k2)*t5)*exp(-(k*t)))*sqrtk)/(945*alpha11*t)
endif

!*******************
case(1)   ! <s1|p2>
!*******************

if(k.eq.1) then
ovlp(0) = ovlp(0) +  ((3274425*t2 + 3274425*t3 + 1600830*t4 + 509355*t5 + 118800*t6 + 21780*t7 + 3300*t8 &
                  + 429*t9 + 49*t10 + 5*t**11)*exp_minus_t)/(16372125*t*sqrt(3.))
else

ovlp(0) = ovlp(0) + (128*k5*((-285120 - 11278080*k2 - 40579200*k4 - 23212800*k6 - 2059200*k8 - 2880*(99 &
                  + 3916*k2 + 14090*k4 + 8060*k6 + 715*k8)*t + 1200*(alpha)*(81 + 2820*k2 + 8750*k4 &
                  + 4180*k6 + 297*k8)*t2 - 480*alpha2*(36 + 1051*k2 + 2571*k4 + 801*k6 + 21*k8)*t3 &
                  + 600*alpha3*(3 + 71*k2 + 129*k4 + 21*k6)*t4 - 4*alpha4*(27 + 501*k2 + 637*k4 &
                  + 35*k6)*t5 + alpha5*(3 + 42*k2 + 35*k4)*t6)*exp_minus_t + (285120 + 11278080*k2 &
                  + 40579200*k4 + 23212800*k6 + 2059200*k8 + 2880*k*(99 + 3916*k2 + 14090*k4 + 8060*k6 &
                  + 715*k8)*t + 240*(alpha)*(189 + 9396*k2 + 40790*k4 + 27460*k6 + 2805*k8)*t2 &
                  + 480*k*alpha2*(273 + 2097*k2 + 1875*k4 + 235*k6)*t3 + 840*alpha3*(3 + 63*k2 &
                  + 81*k4 + 13*k6)*t4 + 200*k*alpha4*(7 + 14*k2 + 3*k4)*t5 + 5*alpha5*(3 &
                  + k2)*(1 + 3*k2)*t6)*exp(-(k*t)))*sqrtk)/(1575*sqrt(3.)*alpha11*t2)

endif

!*******************
case(2)   ! <s1|d2>
!*******************

if(k.eq.1) then
ovlp(0) = ovlp(0) + (t2*(-291060 - 291060*t - 114345*t2 - 17325*t3 + 2475*t4 + 1782*t5 + 437*t6 + 74*t7 &
                  + 10*t8)*exp_minus_t)/(19646550*sqrt(5.))
else

ovlp(0) = ovlp(0) + (128*k5*((-17297280*k - 124225920*k3 - 82857600*k5 - 7862400*k7 - 120960*k*(143 + 1027*k2 &
                  + 685*k4 + 65*k6)*t + 960*k*(-7821 - 52744*k2 - 23770*k4 + 3200*k6 + 495*k8)*t2 + 4800*k*(-1 &
                  + k2)*(363 + 2285*k2 + 1285*k4 + 99*k6)*t3 - 1200*k*alpha2*(201 + 1081*k2 + 483*k4 &
                  + 27*k6)*t4 + 16*k*alpha3*(1287 + 5554*k2 + 1531*k4 + 28*k6)*t5 - 8*k*alpha4*(129 &
                  + 422*k2 + 49*k4)*t6 + 8*k*alpha5*(3 + 7*k2)*t7)*exp_minus_t + (17297280*k + 124225920*k3 &
                  + 82857600*k5 + 7862400*k7 + 120960*k2*(143 + 1027*k2 + 685*k4 + 65*k6)*t + 3840*k*(-297 &
                  - 737*k2 + 11329*k4 + 8965*k6 + 900*k8)*t2 + 1920*k2*(alpha)*(594 + 5071*k2 + 3980*k4 &
                  + 435*k6)*t3 + 480*k*alpha2*(189 + 2061*k2 + 1975*k4 + 255*k6)*t4 + 40*alpha3*(49 &
                  + 1309*k2 + 1723*k4 + 279*k6)*t5 + 200*k*alpha4*(7 + 14*k2 + 3*k4)*t6 + 5*(-1 &
                  + k2)**5*(3 + k2)*(1 + 3*k2)*t7)*exp(-(k*t)))*sqrtk)/(945*sqrt(5.)*alpha11*t3)

endif

!*********************
case(3)   ! <p1|s2>
!*********************

if(k.eq.1) then
ovlp(0) = ovlp(0) +  ((3274425*t2 + 3274425*t3 + 1600830*t4 + 509355*t5 + 118800*t6 + 21780*t7 + 3300*t8 &
                  + 429*t9 + 49*t10 + 5*t**11)*exp_minus_t)/(16372125*t*sqrt(3.))

else

ovlp(0) = ovlp(0) + (128*k5*((-285120 - 11278080*k2 - 40579200*k4 - 23212800*k6 - 2059200*k8 - 2880*(99 &
                  + 3916*k2 + 14090*k4 + 8060*k6 + 715*k8)*t + 1200*(alpha)*(81 + 2820*k2 + 8750*k4 &
                  + 4180*k6 + 297*k8)*t2 - 480*alpha2*(36 + 1051*k2 + 2571*k4 + 801*k6 + 21*k8)*t3 &
                  + 600*alpha3*(3 + 71*k2 + 129*k4 + 21*k6)*t4 - 4*alpha4*(27 + 501*k2 + 637*k4 &
                  + 35*k6)*t5 + alpha5*(3 + 42*k2 + 35*k4)*t6)*exp_minus_t + (285120 + 11278080*k2 &
                  + 40579200*k4 + 23212800*k6 + 2059200*k8 + 2880*k*(99 + 3916*k2 + 14090*k4 + 8060*k6 &
                  + 715*k8)*t + 240*(alpha)*(189 + 9396*k2 + 40790*k4 + 27460*k6 + 2805*k8)*t2 &
                  + 480*k*alpha2*(273 + 2097*k2 + 1875*k4 + 235*k6)*t3 + 840*alpha3*(3 + 63*k2 &
                  + 81*k4 + 13*k6)*t4 + 200*k*alpha4*(7 + 14*k2 + 3*k4)*t5 + 5*alpha5*(3 &
                  + k2)*(1 + 3*k2)*t6)*exp(-(k*t)))*sqrtk)/(1575*sqrt(3.)*alpha11*t2)

endif

!*********************************
case(4)   ! l1 = l2 = 1 => <p1|p2>
!*********************************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((16372125*t + 16372125*t2 + 7203735*t3 + 1746360*t4 + 212355*t5 - 5940*t6 - 8052*t7 &
                  - 1914*t8 - 318*t9 - 43*t10 - 5*t**11)*exp_minus_t)/(16372125*t)
ovlp(1) = ovlp(1) - ((-16372125*t - 16372125*t2 - 7858620*t3 - 2401245*t4 - 522720*t5 - 86130*t6 - 11154*t7 &
                  - 1155*t8 - 93*t9 - 5*t10)*exp_minus_t)/(16372125*t)
else

ovlp(0) = ovlp(0) + (-128*k5*((-823680 - 30251520*k2 - 92678400*k4 - 30251520*k6 - 823680*k8 - 5760*(143 + 5252*k2 &
                  + 16090*k4 + 5252*k6 + 143*k8)*t - 2880*(143 + 5252*k2 + 16090*k4 + 5252*k6 + 143*k8)*t2 &
                  + 240*(alpha)*(467 + 15700*k2 + 45410*k4 + 17908*k6 + 1155*k8)*t3 - 48*alpha2*(379 &
                  + 10914*k2 + 25944*k4 + 7374*k6 + 189*k8)*t4 + 96*alpha3*(19 + 448*k2 + 807*k4 &
                  + 126*k6)*t5 - 4*alpha4*(27 + 501*k2 + 637*k4 + 35*k6)*t6 + alpha5*(3 + 42*k2 &
                  + 35*k4)*t7)*exp_minus_t + (823680 + 30251520*k2 + 92678400*k4 + 30251520*k6 + 823680*k8 + 5760*k*(143 &
                  + 5252*k2 + 16090*k4 + 5252*k6 + 143*k8)*t + 2880*k2*(143 + 5252*k2 + 16090*k4 + 5252*k6 &
                  + 143*k8)*t2 + 240*k*(alpha)*(1155 + 17908*k2 + 45410*k4 + 15700*k6 + 467*k8)*t3 + 48*(-1 &
                  + k2)**2*(189 + 7374*k2 + 25944*k4 + 10914*k6 + 379*k8)*t4 + 96*k*alpha3*(126 + 807*k2 &
                  + 448*k4 + 19*k6)*t5 + 4*alpha4*(35 + 637*k2 + 501*k4 + 27*k6)*t6 + k*alpha5*(35 &
                  + 42*k2 + 3*k4)*t7)*exp(-(k*t)))*sqrtk)/(1575*alpha11*t3)

ovlp(1) = ovlp(1) + (128*k5*((-411840 - 15125760*k2 - 46339200*k4 - 15125760*k6 - 411840*k8 - 2880*(143 + 5252*k2 &
                  + 16090*k4 + 5252*k6 + 143*k8)*t + 720*(alpha)*(187 + 6060*k2 + 16090*k4 + 4444*k6 &
                  + 99*k8)*t2 - 480*alpha2*(47 + 1287*k2 + 2717*k4 + 429*k6)*t3 + 24*alpha3*(91 &
                  + 2047*k2 + 3273*k4 + 189*k6)*t4 - 120*alpha4*(1 + 18*k2 + 21*k4)*t5 + alpha5*(3 &
                  + 42*k2 + 35*k4)*t6)*exp_minus_t + (411840 + 15125760*k2 + 46339200*k4 + 15125760*k6 + 411840*k8 &
                  + 2880*k*(143 + 5252*k2 + 16090*k4 + 5252*k6 + 143*k8)*t + 720*(alpha)*(99 + 4444*k2 &
                  + 16090*k4 + 6060*k6 + 187*k8)*t2 + 480*k*alpha2*(429 + 2717*k2 + 1287*k4 + 47*k6)*t3 &
                  + 24*alpha3*(189 + 3273*k2 + 2047*k4 + 91*k6)*t4 + 120*k*alpha4*(21 + 18*k2 &
                  + k4)*t5 + alpha5*(35 + 42*k2 + 3*k4)*t6)*exp(-(k*t)))*sqrtk)/(1575*alpha11*t3)

endif

!***************************************
case(5)   ! l1 = 1 and l2 = 2 => <p1|d2>
!***************************************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((5239080*t2 + 5239080*t3 + 2411640*t4 + 665280*t5 + 117315*t6 + 11979*t7 + 42*t8 &
                  - 255*t9 - 62*t10 - 10*t**11)*exp_minus_t)/(6548850*t*sqrt(15.))
ovlp(1) = ovlp(1) + ((2619540*t2 + 2619540*t3 + 1257795*t4 + 384615*t5 + 83655*t6 + 13662*t7 + 1713*t8 &
                  + 162*t9 + 10*t10)*exp_minus_t)/(6548850*t*sqrt(5.))

else

ovlp(0) = ovlp(0) + (-128*k5*((-70761600*k - 455414400*k3 - 165836160*k5 - 4717440*k7 - 362880*k*(195 + 1255*k2 &
                  + 457*k4 + 13*k6)*t + 2880*k*(-11713 - 73982*k2 - 22496*k4 + 638*k6 + 33*k8)*t2 + 2880*k*(-3523 &
                  - 21272*k2 - 3302*k4 + 1184*k6 + 33*k8)*t3 + 240*k*(alpha)*(8259 + 49060*k2 + 21386*k4 &
                  + 1908*k6 + 27*k8)*t4 - 240*k*alpha2*(1055 + 5535*k2 + 2229*k4 + 141*k6)*t5 + 24*k*(-1 &
                  + k2)**3*(869 + 3723*k2 + 987*k4 + 21*k6)*t6 - 8*k*alpha4*(129 + 422*k2 + 49*k4)*t7 &
                  + 8*k*alpha5*(3 + 7*k2)*t8)*exp_minus_t + (70761600*k + 455414400*k3 + 165836160*k5 + 4717440*k7 &
                  + 362880*k2*(195 + 1255*k2 + 457*k4 + 13*k6)*t + 5760*k*(-286 + 3601*k2 + 36385*k4 + 13667*k6 &
                  + 393*k8)*t2 + 11520*k2*(-143 - 247*k2 + 5015*k4 + 2035*k6 + 60*k8)*t3 + 480*k*(alpha)*(297 &
                  + 5896*k2 + 24350*k4 + 9480*k6 + 297*k8)*t4 + 48*alpha2*(147 + 6102*k2 + 26412*k4 &
                  + 11722*k6 + 417*k8)*t5 + 48*k*alpha3*(231 + 1617*k2 + 913*k4 + 39*k6)*t6 + 4*(-1 &
                  + k2)**4*(35 + 637*k2 + 501*k4 + 27*k6)*t7 + k*alpha5*(35 + 42*k2 &
                  + 3*k4)*t8)*exp(-(k*t)))*sqrtk)/(315*sqrt(15.)*alpha11*t4)

ovlp(1) = ovlp(1) + (128*k5*((-23587200*k - 151804800*k3 - 55278720*k5 - 1572480*k7 - 120960*k*(195 + 1255*k2 &
                  + 457*k4 + 13*k6)*t + 2880*k*(-3523 - 21272*k2 - 3302*k4 + 1184*k6 + 33*k8)*t2 + 2880*k*(-1 &
                  + k2)*(793 + 4495*k2 + 1399*k4 + 33*k6)*t3 - 240*k*alpha2*(1257 + 6137*k2 + 1539*k4 &
                  + 27*k6)*t4 + 240*k*alpha3*(101 + 402*k2 + 57*k4)*t5 - 24*k*alpha4*(47 + 146*k2 &
                  + 7*k4)*t6 + 8*k*alpha5*(3 + 7*k2)*t7)*exp_minus_t + (23587200*k + 151804800*k3 + 55278720*k5 &
                  + 1572480*k7 + 120960*k2*(195 + 1255*k2 + 457*k4 + 13*k6)*t + 11520*k*(-143 - 247*k2 + 5015*k4 &
                  + 2035*k6 + 60*k8)*t2 + 5760*k2*(alpha)*(286 + 2145*k2 + 900*k4 + 29*k6)*t3 &
                  + 480*k*alpha2*(297 + 2761*k2 + 1371*k4 + 51*k6)*t4 + 24*alpha3*(147 + 3279*k2 &
                  + 2081*k4 + 93*k6)*t5 + 120*k*alpha4*(21 + 18*k2 + k4)*t6 + alpha5*(35 + 42*k2 &
                  + 3*k4)*t7)*exp(-(k*t)))*sqrtk)/(315*sqrt(5.)*alpha11*t4)

endif

!********************
case(6)   ! <d1|s2>
!********************

if(k.eq.1) then
ovlp(0) = ovlp(0) + (t2*(-291060 - 291060*t - 114345*t2 - 17325*t3 + 2475*t4 + 1782*t5 + 437*t6 + 74*t7 &
                  + 10*t8)*exp_minus_t)/(19646550*sqrt(5.))
else

ovlp(0) = ovlp(0) + (128*k5*((-17297280*k - 124225920*k3 - 82857600*k5 - 7862400*k7 - 120960*k*(143 + 1027*k2 &
                  + 685*k4 + 65*k6)*t + 960*k*(-7821 - 52744*k2 - 23770*k4 + 3200*k6 + 495*k8)*t2 + 4800*k*(-1 &
                  + k2)*(363 + 2285*k2 + 1285*k4 + 99*k6)*t3 - 1200*k*alpha2*(201 + 1081*k2 + 483*k4 &
                  + 27*k6)*t4 + 16*k*alpha3*(1287 + 5554*k2 + 1531*k4 + 28*k6)*t5 - 8*k*alpha4*(129 &
                  + 422*k2 + 49*k4)*t6 + 8*k*alpha5*(3 + 7*k2)*t7)*exp_minus_t + (17297280*k + 124225920*k3 &
                  + 82857600*k5 + 7862400*k7 + 120960*k2*(143 + 1027*k2 + 685*k4 + 65*k6)*t + 3840*k*(-297 &
                  - 737*k2 + 11329*k4 + 8965*k6 + 900*k8)*t2 + 1920*k2*(alpha)*(594 + 5071*k2 + 3980*k4 &
                  + 435*k6)*t3 + 480*k*alpha2*(189 + 2061*k2 + 1975*k4 + 255*k6)*t4 + 40*alpha3*(49 &
                  + 1309*k2 + 1723*k4 + 279*k6)*t5 + 200*k*alpha4*(7 + 14*k2 + 3*k4)*t6 + 5*(-1 &
                  + k2)**5*(3 + k2)*(1 + 3*k2)*t7)*exp(-(k*t)))*sqrtk)/(945*sqrt(5.)*alpha11*t3)

endif

!*********************
case(7)  ! <d1|p2> 
!*********************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((5239080*t2 + 5239080*t3 + 2411640*t4 + 665280*t5 + 117315*t6 + 11979*t7 + 42*t8 &
                  - 255*t9 - 62*t10 - 10*t**11)*exp_minus_t)/(6548850*t*sqrt(15.))
ovlp(1) = ovlp(1) + ((2619540*t2 + 2619540*t3 + 1257795*t4 + 384615*t5 + 83655*t6 + 13662*t7 + 1713*t8 &
                  + 162*t9 + 10*t10)*exp_minus_t)/(6548850*t*sqrt(5.))

else

ovlp(0) = ovlp(0) + (-128*k5*((-70761600*k - 455414400*k3 - 165836160*k5 - 4717440*k7 - 362880*k*(195 + 1255*k2 &
                  + 457*k4 + 13*k6)*t + 2880*k*(-11713 - 73982*k2 - 22496*k4 + 638*k6 + 33*k8)*t2 + 2880*k*(-3523 &
                  - 21272*k2 - 3302*k4 + 1184*k6 + 33*k8)*t3 + 240*k*(alpha)*(8259 + 49060*k2 + 21386*k4 &
                  + 1908*k6 + 27*k8)*t4 - 240*k*alpha2*(1055 + 5535*k2 + 2229*k4 + 141*k6)*t5 + 24*k*(-1 &
                  + k2)**3*(869 + 3723*k2 + 987*k4 + 21*k6)*t6 - 8*k*alpha4*(129 + 422*k2 + 49*k4)*t7 &
                  + 8*k*alpha5*(3 + 7*k2)*t8)*exp_minus_t + (70761600*k + 455414400*k3 + 165836160*k5 + 4717440*k7 &
                  + 362880*k2*(195 + 1255*k2 + 457*k4 + 13*k6)*t + 5760*k*(-286 + 3601*k2 + 36385*k4 + 13667*k6 &
                  + 393*k8)*t2 + 11520*k2*(-143 - 247*k2 + 5015*k4 + 2035*k6 + 60*k8)*t3 + 480*k*(alpha)*(297 &
                  + 5896*k2 + 24350*k4 + 9480*k6 + 297*k8)*t4 + 48*alpha2*(147 + 6102*k2 + 26412*k4 &
                  + 11722*k6 + 417*k8)*t5 + 48*k*alpha3*(231 + 1617*k2 + 913*k4 + 39*k6)*t6 + 4*(-1 &
                  + k2)**4*(35 + 637*k2 + 501*k4 + 27*k6)*t7 + k*alpha5*(35 + 42*k2 &
                  + 3*k4)*t8)*exp(-(k*t)))*sqrtk)/(315*sqrt(15.)*alpha11*t4)

ovlp(1) = ovlp(1) + (128*k5*((-23587200*k - 151804800*k3 - 55278720*k5 - 1572480*k7 - 120960*k*(195 + 1255*k2 &
                  + 457*k4 + 13*k6)*t + 2880*k*(-3523 - 21272*k2 - 3302*k4 + 1184*k6 + 33*k8)*t2 + 2880*k*(-1 &
                  + k2)*(793 + 4495*k2 + 1399*k4 + 33*k6)*t3 - 240*k*alpha2*(1257 + 6137*k2 + 1539*k4 &
                  + 27*k6)*t4 + 240*k*alpha3*(101 + 402*k2 + 57*k4)*t5 - 24*k*alpha4*(47 + 146*k2 &
                  + 7*k4)*t6 + 8*k*alpha5*(3 + 7*k2)*t7)*exp_minus_t + (23587200*k + 151804800*k3 + 55278720*k5 &
                  + 1572480*k7 + 120960*k2*(195 + 1255*k2 + 457*k4 + 13*k6)*t + 11520*k*(-143 - 247*k2 + 5015*k4 &
                  + 2035*k6 + 60*k8)*t2 + 5760*k2*(alpha)*(286 + 2145*k2 + 900*k4 + 29*k6)*t3 &
                  + 480*k*alpha2*(297 + 2761*k2 + 1371*k4 + 51*k6)*t4 + 24*alpha3*(147 + 3279*k2 &
                  + 2081*k4 + 93*k6)*t5 + 120*k*alpha4*(21 + 18*k2 + k4)*t6 + alpha5*(35 + 42*k2 &
                  + 3*k4)*t7)*exp(-(k*t)))*sqrtk)/(315*sqrt(5.)*alpha11*t4)

endif

!*********************
case(8)  ! <d1|d2> 
!*********************

if(k.eq.1) then
ovlp(0) = ovlp(0) + ((9823275*t + 9823275*t2 + 3939705*t3 + 665280*t4 - 41580*t5 - 45045*t6 - 10125*t7 &
                  - 1116*t8 - 23*t9 + 19*t10 + 5*t**11)*exp_minus_t)/(9823275*t)
ovlp(1) = ovlp(1) - ((-3274425*t - 3274425*t2 - 1372140*t3 - 280665*t4 - 10395*t5 + 10395*t6 + 3360*t7 + 588*t8 &
                  + 69*t9 + 5*t10)*exp_minus_t)/(3274425*t)
ovlp(2) = ovlp(2) + ((3274425*t + 3274425*t2 + 1548855*t3 + 457380*t4 + 93555*t5 + 13860*t6 + 1500*t7 + 114*t8 &
                  + 5*t9)*exp_minus_t)/(3274425*t)

else

ovlp(0) = ovlp(0) +  (1024*k5*((-65318400*k - 391910400*k3 - 65318400*k5 - 65318400*k*(1 + 6*k2 + k4)*t &
                  + 120960*k*(-257 - 1513*k2 - 163*k4 + 13*k6)*t2 + 120960*k*(-77 - 433*k2 + 17*k4 &
                  + 13*k6)*t3 - 480*k*(4131 + 21724*k2 - 5342*k4 - 452*k6 + 99*k8)*t4 + 30*k*(alpha)*(10299 &
                  + 55840*k2 + 12494*k4 + 1944*k6 + 63*k8)*t5 - 30*k*alpha2*(1155 + 5795*k2 + 1857*k4 &
                  + 153*k6)*t6 + k*alpha3*(2673 + 11291*k2 + 2759*k4 + 77*k6)*t7 - k*alpha4*(129 &
                  + 422*k2 + 49*k4)*t8 + k*alpha5*(3 + 7*k2)*t9)*exp_minus_t + (65318400*k + 391910400*k3 &
                  + 65318400*k5 + 65318400*k2*(1 + 6*k2 + k4)*t + 120960*k*(-13 + 163*k2 + 1513*k4 + 257*k6)*t2 &
                  + 120960*k2*(-13 - 17*k2 + 433*k4 + 77*k6)*t3 + 480*k*(99 - 452*k2 - 5342*k4 + 21724*k6 &
                  + 4131*k8)*t4 + 30*(alpha)*(63 + 1944*k2 + 12494*k4 + 55840*k6 + 10299*k8)*t5 + 30*k*(-1 &
                  + k2)**2*(153 + 1857*k2 + 5795*k4 + 1155*k6)*t6 + alpha3*(77 + 2759*k2 + 11291*k4 &
                  + 2673*k6)*t7 + k*alpha4*(49 + 422*k2 + 129*k4)*t8 + k2*alpha5*(7 &
                  + 3*k2)*t9)*exp(-(k*t)))*sqrtk)/(945*alpha11*t5) 
 
ovlp(1) = ovlp(1) + (-1024*k5*((-14515200*k - 87091200*k3 - 14515200*k5 - 14515200*k*(1 + 6*k2 + k4)*t + 40320*k*(-167 &
                - 973*k2 - 73*k4 + 13*k6)*t2 + 40320*k*(-47 - 253*k2 + 47*k4 + 13*k6)*t3 + 20160*k*(-1 &
                + k2)*(17 + 90*k2 + 13*k4)*t4 - 210*k*alpha2*(193 + 913*k2 + 171*k4 + 3*k6)*t5 &
                + 30*k*alpha3*(103 + 406*k2 + 51*k4)*t6 - 3*k*alpha4*(47 + 146*k2 + 7*k4)*t7 &
                + k*alpha5*(3 + 7*k2)*t8)*exp_minus_t + (14515200*k + 87091200*k3 + 14515200*k5 + 14515200*k2*(1 &
                + 6*k2 + k4)*t + 40320*k*(-13 + 73*k2 + 973*k4 + 167*k6)*t2 + 40320*k2*(-13 - 47*k2 + 253*k4 &
                + 47*k6)*t3 + 20160*k3*(alpha)*(13 + 90*k2 + 17*k4)*t4 + 210*alpha2*(3 + 171*k2 &
                + 913*k4 + 193*k6)*t5 + 30*k*alpha3*(51 + 406*k2 + 103*k4)*t6 + 3*alpha4*(7 &
                + 146*k2 + 47*k4)*t7 + k*alpha5*(7 + 3*k2)*t8)*exp(-(k*t)))*sqrtk)/(315*alpha11*t5)

ovlp(2) = ovlp(2) + (1024*k5*((-3628800*k - 21772800*k3 - 3628800*k5 - 3628800*k*(1 + 6*k2 + k4)*t + 20160*k*(-77 &
                - 433*k2 + 17*k4 + 13*k6)*t2 + 20160*k*(alpha)*(17 + 90*k2 + 13*k4)*t3 - 240*k*(-1 &
                + k2)**2*(183 + 838*k2 + 99*k4)*t4 + 30*k*alpha3*(113 + 426*k2 + 21*k4)*t5 - 150*k*(-1 &
                + k2)**4*(1 + 3*k2)*t6 + k*alpha5*(3 + 7*k2)*t7)*exp_minus_t + (3628800*k + 21772800*k3 &
                + 3628800*k5 + 3628800*k2*(1 + 6*k2 + k4)*t + 20160*k*(-13 - 17*k2 + 433*k4 + 77*k6)*t2 &
                + 20160*k2*(alpha)*(13 + 90*k2 + 17*k4)*t3 + 240*k*alpha2*(99 + 838*k2 + 183*k4)*t4 &
                + 30*alpha3*(21 + 426*k2 + 113*k4)*t5 + 150*k*alpha4*(3 + k2)*t6 + alpha5*(7 &
                + 3*k2)*t7)*exp(-(k*t)))*sqrtk)/(315*alpha11*t5)

endif

end select

! ---------------- End of Select case over NNSEL -------------------------------

end select

! ------------------------------------------------------------------------------

! Correcting the signal, in the case of interchange due the n1 >= n2 imposition:

ovlp(0) = f*ovlp(0)
ovlp(1) = f*ovlp(1)
ovlp(2) = f*ovlp(2)

! End of the subroutine:

end subroutine overlap

!-----------------------------------------------------------------------
!
! Given two l-shells li,lj, centered at atoms ri and rj, this subroutine
! generates all the < ljm | lim > by applying the two-center SK formulas
! 
! Input:
!
!   ovin: sigma, pi & delta components for current <li|lj> interaction
!         (they may correspond to the overlap or the Hamiltonian)
!   cl  : (xj-xi)/|rj-ri|
!   cm  : (yj-yi)/|rj-ri|
!   cn  : (zj-zi)/|rj-ri|
!   cl2,cm2,cn2: squares of the direction cosines cl,cm,cn
!   
!
! Output:
!   ovout: overlap/Hamiltonian elements are given in ovout
!
! Care with signs:
!
! ----------------------------------------------------------------------

subroutine twocenter(li,lj,ovin,ovout,cl,cl2,cm,cm2,cn,cn2)

      implicit none

      integer , intent(in) :: li,lj
      real*8,  intent(in)  :: ovin(3) , cl , cl2 , cm , cm2 , cn ,cn2
      real*8,  intent(out) :: ovout(25)

      integer  :: itype
      integer  :: kao2(0:2,0:2)
      real*8   :: sqr3
      data kao2 / 0 , 3 , 6 , 1 , 4 , 7 , 2 , 5 , 8 /

! ----------------------------------------------------------------------

      sqr3 = sqrt(3.)
      itype = kao2(li,lj)

      select case(itype)

! --  itype = 0 => (S(I):S(J))
      case(0)
        ovout(1)= ovin(1)

! --  itype = 1 => (S(I):P(J))
      case(1)
        ovout(1) = cl*ovin(1)                           ! s,px
        ovout(2) = cm*ovin(1)                           ! s,py
        ovout(3) = cn*ovin(1)                           ! s,pz

! --  itype = 2 => (S(I):D(J))
      case(2)
        ovout(1) = 0.5d0*sqr3*(cl2-cm2)   * ovin(1)     ! s,dx2-y2
        ovout(2) = (cn2-0.5d0*(cl2+cm2))  * ovin(1)     ! s,dz2
        ovout(3) = sqr3*cl*cm             * ovin(1)     ! s,dxy
        ovout(4) = sqr3*cn*cl             * ovin(1)     ! s,dzx
        ovout(5) = sqr3*cm*cn             * ovin(1)     ! s,dyz

! --  itype = 3 => (P(I):S(J))  ! Signal changed here!!!
      case(3)
        ovout(1) =+cl*ovin(1)                           ! px,s
        ovout(2) =+cm*ovin(1)                           ! py,s
        ovout(3) =+cn*ovin(1)                           ! pz,s

! --  itype = 4 => (P(I):P(J))
      case(4)
        ovout(1) =   cl2*ovin(1) + (1.0d0-cl2)*ovin(2)  ! px,px
        ovout(2) = cm*cl*ovin(1) - cm*cl* ovin(2)       ! py,px
        ovout(3) = cn*cl*ovin(1) - cn*cl* ovin(2)       ! pz,px
        ovout(4) = ovout(2)                             ! px,py
        ovout(5) =   cm2*ovin(1) + (1.0d0-cm2)*ovin(2)  ! py,py
        ovout(6) = cn*cm*ovin(1) - cn*cm*ovin(2)        ! pz,py
        ovout(7) = ovout(3)                             ! px,pz
        ovout(8) = ovout(6)                             ! py,pz
        ovout(9) =   cn2*ovin(1) + (1.0d0-cn2)*ovin(2)  ! pz,pz

! --  itype = 5 => (P(I):D(J))
      case(5)
        ovout( 1)= 0.5d0*sqr3*cl*(cl2-cm2)*ovin(1)  + cl*(1.0d0-cl2+cm2)*ovin(2)  ! px,dx2
        ovout( 2)= 0.5d0*sqr3*cm*(cl2-cm2)*ovin(1)  - cm*(1.0d0+cl2-cm2)*ovin(2)  ! py,dx2
        ovout( 3)= 0.5d0*sqr3*cn*(cl2-cm2)*ovin(1)  - cn*(cl2-cm2)*ovin(2)        ! pz,dx2
        ovout( 4)= cl*(cn2-0.5d0*(cl2+cm2))*ovin(1) - sqr3*cl*cn2*ovin(2)         ! px,dz2
        ovout( 5)= cm*(cn2-0.5d0*(cl2+cm2))*ovin(1) - sqr3*cm*cn2* ovin(2)        ! py,dz2
        ovout( 6)= cn*(cn2-0.5d0*(cl2+cm2))*ovin(1) + sqr3*cn*(cl2+cm2)* ovin(2)  ! pz,dz2
        ovout( 7)= sqr3*cl2*cm*ovin(1) + cm*(1.0d0-2.0d0*cl2)* ovin(2)            ! pxdxy
        ovout( 8)= sqr3*cm2*cl*ovin(1) + cl*(1.0d0-2.0d0*cm2)* ovin(2)            ! pydxy
        ovout( 9)= sqr3*cl*cm*cn*ovin(1) - 2.0d0*cl*cm*cn* ovin(2)                ! pzdxy
        ovout(10)= sqr3*cl2*cn*ovin(1) + cn*(1.0d0-2.0d0*cl2)* ovin(2)            ! px,dzx
        ovout(11)= ovout(9)                                                       ! py,dzx
        ovout(12)= sqr3*cn2*cl*ovin(1) + cl*(1.0d0-2.0d0*cn2)* ovin(2)            ! pz,dzx
        ovout(13)= ovout(9)                                                       ! px,dyz
        ovout(14)= sqr3*cm2*cn*ovin(1) + cn*(1.0d0-2.0d0*cm2)* ovin(2)            ! py,dyz
        ovout(15)= sqr3*cn2*cm*ovin(1) + cm*(1.0d0-2.0d0*cn2)* ovin(2)            ! py,dyz

! --  itype = 6 => (D(I):S(J))
      case(6)
        ovout(1) =  0.5d0*sqr3*(cl2-cm2)  * ovin(1)    ! s,dx2
        ovout(2) = (cn2-0.5d0*(cl2+cm2))  * ovin(1)    ! s,dz2
        ovout(3) = sqr3*cl*cm             * ovin(1)    ! s,dxy
        ovout(4) = sqr3*cn*cl             * ovin(1)    ! s,dzx
        ovout(5) = sqr3*cm*cn             * ovin(1)    ! s,dyz

! --  itype = 7 => (D(I):P(J))   ! Signals changed here too
      case(7)
        ovout( 1) =-(-0.5d0*sqr3*cl*(cl2-cm2)*ovin(1) - cl*(1.0d0-cl2+cm2)*ovin(2))   ! dx2,px
        ovout( 2) =-(-cl*(cn2-0.5d0*(cl2+cm2))*ovin(1) + sqr3*cl*cn2*ovin(2)      )   ! dz2,px
        ovout( 3) =-(-sqr3*cl2*cm*ovin(1) - cm*(1.0d0-2.0d0*cl2)*ovin(2)          )   ! dxy,px
        ovout( 4) =-(-sqr3*cl2*cn*ovin(1) - cn*(1.0d0-2.0d0*cl2)*ovin(2)          )   ! dzx,px
        ovout( 5) =-(-sqr3*cl*cm*cn*ovin(1) + 2.0d0*cl*cm*cn*ovin(2)              )   ! dyz,px
        ovout( 6) =-(-0.5d0*sqr3*cm*(cl2-cm2)*ovin(1) + cm*(1.0d0+cl2-cm2)*ovin(2))   ! dx2,py
        ovout( 7) =-(-cm*(cn2-0.5d0*(cl2+cm2))*ovin(1) + sqr3*cm*cn2*ovin(2)      )   ! dz2,py
        ovout( 8) =-(-sqr3*cm2*cl*ovin(1) - cl*(1.0d0-2.0d0*cm2)*ovin(2)          )   ! dxy,py
        ovout( 9) = ovout(5)                                                          ! dzx,py
        ovout(10) =-(-sqr3*cm2*cn*ovin(1) - cn*(1.0d0-2.0d0*cm2)*ovin(2)          )   ! dyz,py
        ovout(11) =-(-0.5d0*sqr3*cn*(cl2-cm2)*ovin(1) + cn*(cl2-cm2)*ovin(2)      )   ! dx2,pz
        ovout(12) =-(-cn*(cn2-0.5d0*(cl2+cm2))*ovin(1) - sqr3*cn*(cl2+cm2)*ovin(2))   ! dz2,pz
        ovout(13) = ovout(5)                                                          ! dxy,pz
        ovout(14) =-(-sqr3*cn2*cl*ovin(1) - cl*(1.0d0-2.0d0*cn2)*ovin(2)          )   ! dzx,pz
        ovout(15) =-(-sqr3*cn2*cm*ovin(1) - cm*(1.0d0-2.0d0*cn2)*ovin(2)          )   ! dyz,pz

! --  itype = 8 =>  (D(I):D(J))
      case(8)
        ovout( 1)= 0.75d0*((cl2-cm2)**2)*ovin(1) + (cl2+cm2-(cl2-cm2)**2)*ovin(2) &
                  + (cn2+0.25*(cl2-cm2)**2)*ovin(3)                                      ! dx2,dx2
        ovout( 2)= 0.5d0*sqr3*(cl2-cm2)*(cn2-0.5d0*(cl2+cm2))*ovin(1) &
                  + sqr3*cn2*(cm2-cl2)*ovin(2) + (0.25*sqr3)*(1.+cn2)*(cl2-cm2)*ovin(3)  ! dz2,dx2
        ovout( 3)= 1.5d0*cl*cm*(cl2-cm2)*ovin(1) + 2.0d0*cl*cm*(cm2-cl2)*ovin(2) &
                  + 0.5d0*cl*cm*(cl2-cm2)*ovin(3)                                        ! dxy,dx2
        ovout( 4)= 1.5d0*cn*cl*(cl2-cm2)*ovin(1) + cn*cl*(1. - 2.*(cl2-cm2))*ovin(2) &
                  - cn*cl*(1.0d0-0.5d0*(cl2-cm2))*ovin(3)                                ! dzx,dx2
        ovout( 5)= 1.5d0*cm*cn*(cl2-cm2)*ovin(1) - cm*cn*(1. + 2.*(cl2-cm2))*ovin(2) &
                  + cm*cn*(1.0d0+0.5d0*(cl2-cm2))*ovin(3)                                ! dyz,dx2
        ovout( 6)= ovout(2)                                                              ! dx2,dz2 = dz2,dx2
        ovout( 7)= (cn2-0.5d0*(cl2+cm2))**2*ovin(1) + 3.0d0*cn2*(cl2+cm2)*ovin(2) &
                  +0.75d0*(cl2+cm2)**2*ovin(3)                                          ! dz2,dz2

        ovout( 8)= sqr3*cl*cm*(cn2-0.5*(cl2+cm2))*ovin(1) - sqr3*2.*cl*cm*cn2*ovin(2) &
                  +0.5d0*sqr3*cl*cm*(1.0d0+cn2)*ovin(3)                                 ! dxy,dz2
        ovout( 9)= sqr3*cl*cn*(cn2-0.5*(cl2+cm2))*ovin(1) &
                  +sqr3*cl*cn*(cl2+cm2-cn2)*ovin(2) - 0.5*sqr3*cl*cn*(cl2+cm2)*ovin(3)  ! dzx,dz2
        ovout(10) = sqr3*cm*cn*(cn2-0.5d0*(cl2+cm2))*ovin(1) &
                   +sqr3*cm*cn*(cl2+cm2-cn2)*ovin(2)- 0.5*sqr3*cm*cn*(cl2+cm2)*ovin(3)  ! dyz,dz2
        ovout(11)= ovout(3)                                                             ! dx2,dxy = dxy,dx2
        ovout(12)= ovout(8)                                                             ! dz2,dxy = dxy,dz2
        ovout(13)= 3.0d0*cl2*cm2*ovin(1) + (cl2+cm2-4.0d0*cl2*cm2)*ovin(2) &
                  + (cn2+cl2*cm2)*ovin(3)                                               ! dxy,dxy
        ovout(14)= 3.0d0*cn*cl2*cm*ovin(1) + cn*cm*(1.0d0-4.0d0*cl2)*ovin(2) &
                  + cn*cm*(cl2-1.0d0)*ovin(3)                                           ! dzx,dxy
        ovout(15)= 3.0d0*cm2*cn*cl*ovin(1) + cn*cl*(1.0d0-4.0d0*cm2)*ovin(2) &
                  + cn*cl*(cm2-1.0d0)*ovin(3)                                           ! dyz,dxy
        ovout(16)= ovout(4)                                                             ! dx2,dzx = dzx,dx2
        ovout(17)= ovout(9)                                                             ! dz2,dzx = dzx,dz2
        ovout(18)= ovout(14)                                                            ! dxy,dzx = dzx,dxy
        ovout(19)= 3.0d0*cn2*cl2*ovin(1) + (cn2+cl2-4.0d0*cn2*cl2)*ovin(2) &
                  + (cm2+cn2*cl2)*ovin(3)                                               ! dzx,dzy
        ovout(20)= 3.0d0*cm*cn2*cl*ovin(1) + cm*cl*(1.0d0-4.0d0*cn2)*ovin(2) &
                  + cm*cl*(cn2-1.0d0)*ovin(3)                                           ! dyz,dzx
        ovout(21)= ovout(5)                                                             ! dx2,dyz = dyz,dx2
        ovout(22)= ovout(10)                                                            ! dz2,dyz = dyz,dz2
        ovout(23)= ovout(15)                                                            ! dxy,dyz = dyz,dxy
        ovout(24)= ovout(20)                                                            ! dzx,dyz = dyz,dzx
        ovout(25)= 3.0d0*cm2*cn2*ovin(1) + (cm2+cn2-4.0d0*cm2*cn2)*ovin(2) &
                  + (cl2+cm2*cn2)*ovin(3)                                               !dyz,dyz

      case default
       write(*,'(a,i3,/a)') 'twocenter: parameter ITYPE had an ilegal value =', itype
       stop

      end select

      return
      end subroutine twocenter

! -----------------------------------------------------------------------------
end module sto_ov
! -----------------------------------------------------------------------------
