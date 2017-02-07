! -----------------------------------------------------------------------------
! MODULE TOTEN: performs the total energy calculation                     
!
! (1) TBTE:   Main subroutine, that drives the calculation.
!
! (2) IONPOT: Calculates the repulsive contribution to the total energy due the
!             ion-ion interaction.
!
! LAST MODIFICATION: 20/08/2012
! -----------------------------------------------------------------------------

module toten  
use diagonalize
implicit none

public:: tbte, ionpot                         

CONTAINS

!-----------------------------------------------------------------------
! SUBROUTINE TBTE: Subroutine that drives the total energy computation
!-----------------------------------------------------------------------

subroutine tbte(kcalc,ibrav,numtp,nat,tnao,norbv,itype,neao,nx,ny,nz,rcut,rx,ry,rz,pv1,pv2,pv3,oe,&
                esit,keht,cf,zeta,sip,ldm,dbov,lptb,nzt,nval,lval,lpar,eshift,nelect,nkpt,kpoint)

implicit none

! ----------------- Variables and Parameters of the Subroutine ----------------
! Externals:
integer, parameter  :: lmax = 2
integer, intent(in) :: kcalc, ibrav, numtp, nat, tnao, norbv(numtp), itype(nat)
integer, intent(in) :: neao(numtp,3), nx, ny, nz
integer, intent(in) :: sip(nat), ldm(0:lmax), dbov(numtp), lptb(numtp,0:lmax)
integer, intent(in) :: nzt(numtp,3), nval(numtp,3), lval(numtp,3), nelect, nkpt
real*8 , intent(in) :: oe(numtp,9), esit(numtp,3), cf(numtp,3,2), zeta(numtp,3,2)
real*8 , intent(in) :: keht(numtp,numtp), rcut, rx(nat), ry(nat), rz(nat)
real*8 , intent(in) :: pv1(3), pv2(3), pv3(3), lpar, eshift, kpoint(nkpt,3)
! Internals:
real*8, allocatable     :: eband(:,:), hr(:,:), sr(:,:), W(:)
complex*16, allocatable :: hamk(:,:,:), ovlk(:,:,:)
integer :: i, k, flp, imax(1), imin(1)
real*8  :: eele, erep, etot

! -----------------------------------------------------------------------------

! Fermi Level position:

if(mod(nelect,2).eq.0) then
 flp = nelect/2
  else
 flp = (nelect + 1)/2
endif

! ------------- Calculating the Electronic Energy According KCALC -------------

select case(kcalc)

case(4)  ! Real supercell calculation

! Allocating HR, SR and W:

allocate(hr(tnao,tnao), sr(tnao,tnao), W(tnao))

call fcspec(ibrav,numtp,nat,nx,ny,nz,tnao,neao,norbv,itype,rcut,rx,ry,rz,pv1,pv2,pv3,oe,&
            esit,keht,cf,zeta,ldm,sip,dbov,lptb,nzt,nval,lval,lpar,eshift,nelect,hr,sr,W)

! Electronic energy:

eele = 0.0d0
do i = 1, flp
 eele = eele + W(i)
enddo

case(5)  ! Energy spectrum in the k-speace:

! Allocating EBAND, HAMK and OVLK:

allocate(eband(nkpt,tnao),hamk(nkpt,tnao,tnao),ovlk(nkpt,tnao,tnao))

 eband(:,:) = 0.0d0
hamk(:,:,:) = cmplx(0.0d0,0.0d0)
ovlk(:,:,:) = cmplx(0.0d0,0.0d0)

call bands(kcalc,ibrav,numtp,nat,tnao,norbv,itype,rx,ry,rz,pv1,pv2,pv3,oe,esit,keht,cf,zeta,&
           sip,ldm,dbov,lptb,nzt,nval,lval,lpar,eshift,nelect,nkpt,kpoint,hamk,ovlk,eband)

! Electronic energy:

eele = 0.0d0
do k = 1, nkpt
 do i = 1, flp
  eele = eele + eband(k,i)
 enddo
enddo

eele = eele/dfloat(nkpt)

case default

write(*,'(a)') 'KCALC value not recognize'

end select

! -------------------- Repulsive Energy Contribution  -------------------------

call ionpot(kcalc,nat,numtp,itype,rcut,rx,ry,rz,erep)

etot = eele + erep

write(*,'(a,f14.4)') 'Electronic Energy = ', eele
write(*,'(a,f14.4)') 'Repulsive  Energy = ', erep
write(*,'(a,f14.4)') 'Total Energy (eV) = ', etot

end subroutine tbte   

! -----------------------------------------------------------------------------
! Subroutine IONPOT: calculates the repulsive energy due the interaction among
!                    the ions.
! -----------------------------------------------------------------------------

subroutine ionpot(kcalc,nat,numtp,itype,rcut,rx,ry,rz,erep)
implicit none 

! ----------------- Variables and Parameters of the Subroutine ----------------
! Externals:
integer, intent(in)   :: kcalc, nat, numtp, itype(nat)
real*8, intent(in)    :: rcut, rx(nat), ry(nat), rz(nat)
real*8, intent(inout) :: erep
!Internals:
integer             :: i, j, k
real*8              :: chi0, d0, alpha, rij
! -----------------------------------------------------------------------------

! Parameters for Silicon (Bernstein, PRB 56, 10488 (1997)):

   d0 = 2.35d0
! chi0 = 0.0082d0
!alpha = 1.72d0


 chi0 = 4.11d0  
alpha = 1.62d0

! Loop over all pairs:

erep = 0.0d0

do i = 1, nat-1
  do j = i + 1, nat
    rij = dsqrt((rx(j)-rx(i))**2 + (ry(j)-ry(i))**2 + (rz(j)-rz(i))**2)
      if(rij.lt.rcut) then
        erep = erep + chi0*dexp(-4.d0*alpha*(rij - d0))
      endif
  enddo
enddo

end subroutine ionpot

! -----------------------------------------------------------------------------
end module toten
! -----------------------------------------------------------------------------
