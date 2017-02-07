! -----------------------------------------------------------------------------
! MODULE DIAGONALIZE: In this module we have the subroutines that diagonalizes 
! The system hamiltonian. The subroutine were taken from the numerical recipes
!
! (1) BANDS: calculates the band structure of a material.
!
! (2) SPECTRUM: computes the energy spectrum of a supercell
!
! (3)  FCSPEC: Computes the energy spectrun of a supercell for IBRAV = 6 (free
!              cell mode)
!
! (4) KPOINTS: Generates the k-points according with the bravais lattice for 
!              Band-Structure Calculations
!
! (5) MONKPACK: Generates the k-points according the Mokhost-Pack recipe to
!               Calculate the DOS (Ref.: ...)
!
! LAST MODIFICATION: 20/08/2012
! -----------------------------------------------------------------------------

module diagonalize
use sto_ov
use postproc
implicit none

public:: bands, spectrum, fcspec, kpoints, monkpack

CONTAINS

!-----------------------------------------------------------------------
! SUBROUTINE BANDS: Calculates the band structure of the system for
!                   The minimum unit cell mode
!-----------------------------------------------------------------------

subroutine bands(kcalc,ibrav,numtp,nat,tnao,norbv,itype,rx,ry,rz,pv1,pv2,pv3,oe,esit,keht,cf,zeta,&
           sip,ldm,dbov,lptb,nzt,nval,lval,lpar,eshift,nelect,nkpt,kpoint,hamk,ovlk,eband)
implicit none

! ----------------- Variables and Parameters of the Subroutine ----------------
! Externals:
integer, parameter  :: lmax = 2, nmax = 5
integer, intent(in) :: kcalc, ibrav, numtp, nat, tnao, norbv(numtp), itype(nat)
integer, intent(in) :: sip(nat), ldm(0:lmax), dbov(numtp), lptb(numtp,0:lmax)
integer, intent(in) :: nzt(numtp,3), nval(numtp,3), lval(numtp,3), nelect, nkpt
real*8 , intent(in) :: oe(numtp,9), esit(numtp,3), cf(numtp,3,2), zeta(numtp,3,2)
real*8 , intent(in) :: keht(numtp,numtp), rx(nat), ry(nat), rz(nat)
real*8 , intent(in) :: pv1(3), pv2(3), pv3(3), lpar, kpoint(nkpt,3)
real*8 , intent(inout) :: eband(nkpt,tnao)
complex*16, intent(inout) :: hamk(nkpt,tnao,tnao), ovlk(nkpt,tnao,tnao)
! Internals:
integer, parameter  :: ncell = 4
integer             :: llsel(0:lmax,0:lmax), nnsel(nmax,nmax), mdcount, ipbc
integer             :: i, j, oi, oj, ct, ix, iy, iz, jo, jd, io, id, li, lj
integer             :: izt, izoi, izoj, ik, idv, ikcount
real*8              :: ovout(25), ui(3), uj(3), rij(3), kp(nkpt,3), eshift
real*8              :: kpb(3), kdotr, pi, kh
complex*16          :: hk(tnao,tnao), sk(tnao,tnao)
! Parametros para a subrotina de diagonalizacao (LAPACK, ZHEGVD):
integer                  :: itp, ndm, lda, ldb, lwork, lrwork, liwork, info
complex*16, allocatable  :: A(:,:), B(:,:), work(:)
real*8, allocatable      :: w(:), rwork(:)
integer, allocatable     :: iwork(:)
character*1              :: jobz, uplo
! ------------------------------------------------------------------------

! PI value and Initializing the EGAP value:

pi = 3.141592654d0

! In this mode, we assign IPBC = 0 (no periodic boundary conditions):

ipbc = 0

! -----------------------------------------------------------------------------
! Setting the parameters for the lapack diagonalization subroutine:

itp = 1
 jobz = 'V'
  uplo = 'U'
   ndm = tnao
    lda = ndm
   ldb = ndm
  lwork = 2*ndm + ndm*ndm
 lrwork = 1 + 5*ndm + 2*ndm*ndm
liwork = 3 + 5*ndm
allocate(A(ndm,ndm),B(ndm,ndm),W(ndm),work(lwork),rwork(lrwork),iwork(liwork))
! -----------------------------------------------------------------------------

! Loading the arrays that control the cases for the overlap computation:

call arrsel(nnsel,llsel)

! Loop over the K-points and the divisions:

ikcount = 0

do ik = 1, nkpt     

kpb(1) = kpoint(ik,1)
kpb(2) = kpoint(ik,2)
kpb(3) = kpoint(ik,3)

!write(*,*) ik, kpb(1), kpb(2), kpb(3)

! Initializing the arrays SK and HK:

sk(:,:) = cmplx(0.0d0,0.0d0)
hk(:,:) = cmplx(0.0d0,0.0d0)

! Loading the main diagonal of the hamiltonian:

mdcount = 0
do i = 1, nat
 do j = 1, dbov(itype(i))
  mdcount = mdcount + 1
   hk(mdcount,mdcount) = cmplx(oe(itype(i),j),0.0d0)
 enddo
enddo

! Building the hamiltonian:

mdcount = 0

do i = 1, nat                ! Loop over the atoms in the central cell:

ui(1) = rx(i)
ui(2) = ry(i)
ui(3) = rz(i)

do oi = 1, norbv(itype(i))   ! Loop over the I's orbitals:

do j = 1, nat                ! Loop over the J atoms

do oj = 1, norbv(itype(j))   ! Loop over the J's orbitals:

li = lval(itype(i),oi)
lj = lval(itype(j),oj)

! ------------------------------------------------------------------------
!                       Loop over the units cells
! ------------------------------------------------------------------------

do ix = -ncell, ncell
  do iy = -ncell, ncell
    do iz = -ncell, ncell

uj(1) = rx(j) + ix*pv1(1) + iy*pv2(1) + iz*pv3(1)
uj(2) = ry(j) + ix*pv1(2) + iy*pv2(2) + iz*pv3(2)
uj(3) = rz(j) + ix*pv1(3) + iy*pv2(3) + iz*pv3(3)

! Relative position

rij(1) = uj(1) - ui(1)
rij(2) = uj(2) - ui(2)
rij(3) = uj(3) - ui(3)

! Excluding the case j = i in the central cell:

if((j.eq.i).and.(ix.eq.0).and.(iy.eq.0).and.(iz.eq.0)) then 

mdcount = mdcount + 1

! in the case of NZT = 2 the orbital is always normalized, otherwise, its value squared
! is equal to cf**2 (??????):

 sk(mdcount,mdcount) = sk(mdcount,mdcount) + cmplx(1.0d0,0.0d0)  ! Normalized orbitals  

else

call calcovlp(i,j,oi,oj,nat,numtp,nzt,nval,lval,itype,zeta,cf,ui,uj,&
              llsel,nnsel,ovout,ipbc)

! VERIFICAR CUIDADOSAMENTE ESTA PARTE!!!!!

ct = 0
do jo = 1, ldm(lj)

  jd = sip(j) + lptb(itype(j),lj) + jo - 1

  do io = 1, ldm(li)

   id = sip(i) + lptb(itype(i),li) + io - 1

   ct = ct + 1

   kdotr = dot_product(kpb,rij)
   kh = keht(itype(i),itype(j))
   sk(id,jd) = sk(id,jd) + ovout(ct)*cmplx(cos(kdotr),sin(kdotr))   ! AQUI MORA O DEMO!!!
   hk(id,jd) = hk(id,jd) + 0.5d0*kh*(esit(itype(i),oi) + esit(itype(j),oj))*ovout(ct)*cmplx(cos(kdotr),sin(kdotr))

 enddo
enddo

endif

  enddo   ! End of IZ loop
 enddo   ! End of IY loop
enddo   ! End of IX loop
! ------------------------------------------------------------------------

enddo         ! End of the loop over the J's orbitals

enddo         ! End of the Loop over the J atoms

enddo         ! End of the loop over the I's orbitals

enddo         ! End of the loop over the I atoms

! ------------------------------------------------------------------------
!                    DIAGONALIZATION OF THE HAMILTONIAN
! ------------------------------------------------------------------------

! Apllying the shift on the energies of the hamiltonian:

hk(:,:) = hk(:,:) + eshift*sk(:,:)

A(:,:) = hk(:,:)
B(:,:) = sk(:,:)

call zhegvd(itp, jobz, uplo, ndm, A, lda, B, ldb, W, work, lwork, rwork, lrwork, IWORK, LIWORK, INFO)

! Storing the bands in EBAND, HK(k) and OVLK(k) to calculate DOS:

ikcount = ikcount + 1
eband(ikcount,:)   =  W(:)
 hamk(ikcount,:,:) =  A(:,:)
 ovlk(ikcount,:,:) = sk(:,:)

W(:) = 0.0d0

enddo         ! End of the loop over the k-points

! ------------------------------------------------------------------------------
! End of the subroutine
! ------------------------------------------------------------------------------

end subroutine bands

!------------------------------------------------------------------------------
! SUBROUTINE SPECTRUM: Calculates the eigen-states for a supercell Hamiltonian
!                      HR(TNAO,TNAO).  
!------------------------------------------------------------------------------

subroutine spectrum(numtp,nat,tnao,neao,norbv,itype,nnn,nngh,lisngh,rx,ry,rz,oe,esit,keht,cf,&
                    zeta,ldm,sip,dbov,lptb,nzt,nval,lval,lpar,eshift,sdx,sdy,sdz,nelect,hr,sr,W)
implicit none

! ----------------- Variables and Parameters of the Subroutine ----------------
! Externals:
integer, parameter     :: lmax = 2, nmax = 5
integer, intent(in)    :: numtp, nat, tnao, neao(numtp,3), norbv(numtp), itype(nat), nnn
integer, intent(in)    :: ldm(0:lmax), sip(nat), dbov(numtp), lptb(numtp,0:lmax)
integer, intent(in)    :: nngh(nat), lisngh(nat,nnn), nzt(numtp,3), nval(numtp,3), lval(numtp,3), nelect
real*8 , intent(in)    :: oe(numtp,9), esit(numtp,3), cf(numtp,3,2), zeta(numtp,3,2)
real*8 , intent(in)    :: keht(numtp,numtp), rx(nat), ry(nat), rz(nat), lpar, eshift
real*8 , intent(in)    :: sdx, sdy, sdz
real*8 , intent(inout) :: hr(tnao,tnao), sr(tnao,tnao), W(tnao)
! Internals:
integer                :: llsel(0:lmax,0:lmax), nnsel(nmax,nmax), mdcount, ipbc
integer                :: i, j, oi, oj, ct, jng, jo, jd, io, id, li, lj, k
integer                :: izt, izoi, izoj, ik, idv, itest, it, jt, iot, jot
real*8                 :: ovout(25), ui(3), uj(3), rij(3), kh, sum
! Parametros para a subrotina de diagonalizacao (LAPACK, DSYGVD):
integer                :: itp, ndm, lda, ldb, lwork, liwork, info
real*8, allocatable    :: work(:)
character*1            :: jobz, uplo
integer, allocatable   :: iwork(:)
!-----------------------------------------------------------------------

! Initializing the arrays SR (overlaps) and HR (hamiltonian):

sr(:,:) = 0.0d0
hr(:,:) = 0.0d0

! Loading the main diagonal of the hamiltonian and the overlap matrices:

mdcount = 0
do i = 1, nat
 do j = 1, dbov(itype(i))
  mdcount = mdcount + 1
   hr(mdcount,mdcount) = oe(itype(i),j)
   sr(mdcount,mdcount) = 1.0d0
 enddo
enddo

! In this mode, we assign IPBC = 1 (periodic boundary conditions allowed):

ipbc = 1

! Loading the arrays that control the cases for the overlap computation:

call arrsel(nnsel,llsel)

! Building the hamiltonian:

mdcount = 0

do i = 1, nat                ! Loop over the I atoms:

ui(1) = rx(i)
ui(2) = ry(i)
ui(3) = rz(i)

do oi = 1, norbv(itype(i))   ! Loop over the I's orbitals:

do jng = 1, nngh(i)          ! Loop over the I's neighbors J

j = lisngh(i,jng)

if(j.gt.i) then

uj(1) = rx(j)
uj(2) = ry(j)
uj(3) = rz(j)

do oj = 1, norbv(itype(j))   ! Loop over the J's orbitals:

li = lval(itype(i),oi)
lj = lval(itype(j),oj)

call calcovlp(i,j,oi,oj,nat,numtp,nzt,nval,lval,itype,zeta,cf,ui,uj,&
              llsel,nnsel,ovout,ipbc)

! ---------- VERIFICAR CUIDADOSAMENTE ESTA PARTE!!!!! -----------------------

ct = 0
do jo = 1, ldm(lj)
! jd = (j - 1)*dbov(itype(j)) + lptb(itype(j),lj) + jo - 1
  jd = sip(j) + lptb(itype(j),lj) + jo - 1
  do io = 1, ldm(li)

!  id = (i - 1)*dbov(itype(i)) + lptb(itype(i),li) + io - 1
   id = sip(i) + lptb(itype(i),li) + io - 1
   ct = ct + 1

   kh = keht(itype(i),itype(j))
   sr(id,jd) = sr(id,jd) + ovout(ct)   ! AQUI MORA O DEMO!!!
   hr(id,jd) = hr(id,jd) + 0.5d0*kh*(esit(itype(i),oi) + esit(itype(j),oj))*ovout(ct)
   sr(jd,id) = sr(id,jd)
   hr(jd,id) = hr(id,jd)

 enddo
enddo

! ------------------------------------------------------------------------

enddo         ! End of the loop over the J's orbitals

endif

enddo         ! End of the Loop over the J atoms

enddo         ! End of the loop over the I's orbitals

enddo         ! End of the loop over the I atoms

! ------------------------------------------------------------------------
!               VERIFICACAO DA HAMILTONIANA E DO OVERLAP
! ------------------------------------------------------------------------

!do itest = 1, 1000
!write(*,*) 'Entre com o sitio I e J e os respectivos indices dos orbitais'
!read(5,*) it, iot, jt, jot
!write(*,*) dsqrt(rx(it) - rx(jt)**2 + (ry(it) - ry(jt))**2 + (rz(it) - rz(jt))**2)
!it = (it - 1)*9 + iot
!jt = (jt - 1)*9 + jot
!write(*,*) it, jt 
!write(*,*) hr(it,jt), sr(it,jt)
!enddo

! ------------------------------------------------------------------------
!                    DIAGONALIZATION OF THE HAMILTONIAN
! ------------------------------------------------------------------------

! Apllying the shift on the energies of the hamiltonian:

hr(:,:) = hr(:,:) + eshift*sr(:,:)

! Setting the parameters for the lapack diagonalization subroutine:

itp = 1
jobz = 'V'
uplo = 'U'
ndm = tnao
lda = ndm
ldb = ndm
lwork = 1 + 6*ndm + 2*ndm*ndm 
liwork = 3 + 5*ndm
allocate(work(lwork),iwork(liwork))

call dsygvd(itp, jobz, uplo, ndm, hr, lda, sr, ldb, W, work, lwork, iwork, liwork, info)

! ------------------------------------------------------------------------
!                    NORMALIZATION OF THE EIGENSTATES  
! ------------------------------------------------------------------------

do i = 1, tnao ! All states will be normalized

sum = 0.0d0

  do j = 1, tnao
    do k = 1, tnao
      sum = sum + hr(j,i)*hr(k,i)*sr(k,j)
    enddo
  enddo

hr(:,i) = hr(:,i)/dsqrt(sum)

enddo

! End of the Subroutine

end subroutine spectrum

!-----------------------------------------------------------------------
! SUBROUTINE FCSPEC: Calculates the band structure of the system for
!                    A free unit cell
!-----------------------------------------------------------------------

subroutine fcspec(ibrav,numtp,nat,nx,ny,nz,tnao,neao,norbv,itype,rcut,rx,ry,rz,pv1,pv2,pv3,oe,&
                  esit,keht,cf,zeta,ldm,sip,dbov,lptb,nzt,nval,lval,lpar,eshift,nelect,hr,sr,W)
implicit none

! ----------------- Variables and Parameters of the Subroutine ----------------
! Externals:
integer, parameter     :: lmax = 2, nmax = 5
integer, intent(in)    :: numtp, nat, nx, ny, nz, tnao, norbv(numtp), neao(numtp,3), nelect
integer, intent(in)    :: itype(nat), ldm(0:lmax), sip(nat), dbov(numtp), lptb(numtp,0:lmax)
integer, intent(in)    :: nzt(numtp,3), nval(numtp,3), lval(numtp,3), ibrav
real*8 , intent(in)    :: oe(numtp,9), esit(numtp,3), cf(numtp,3,2), zeta(numtp,3,2)
real*8 , intent(in)    :: keht(numtp,numtp), rx(nat), ry(nat), rz(nat), pv1(3), pv2(3), pv3(3)
real*8 , intent(in)    :: lpar, rcut
real*8 , intent(inout) :: hr(tnao,tnao), sr(tnao,tnao), W(tnao)
! Internals:
integer, parameter     :: ncell = 1
integer                :: llsel(0:lmax,0:lmax), nnsel(nmax,nmax), mdcount, ipbc
integer                :: i, j, oi, oj, ct, ix, iy, iz, jo, jd, io, id, li, lj, k
integer                :: izt, izoi, izoj, ik, idv, flp
real*8                 :: ovout(25), ui(3), uj(3), rij, kpb(3), eshift, kh
real*8                 :: sdx, sdy, sdz, sv1(3), sv2(3), sv3(3), sum
! Parametros para a subrotina de diagonalizacao (LAPACK, DSYGVD):
integer                :: itp, ndm, lda, ldb, lwork, liwork, info
real*8, allocatable    :: work(:)
character*1            :: jobz, uplo
integer, allocatable   :: iwork(:)
! ------------------------------------------------------------------------

! In this mode, we assign IPBC = 0 (no periodic boundary conditions):

ipbc = 0

! -------------------- Defining the supercell vectors -------------------------

select case(ibrav) 

case(1,2,3,4)  ! Supercells with cubic symmetry

sv1(:) = (/nx*lpar, 0.0d0, 0.0d0/)  !pv1(:)
sv2(:) = (/0.0d0, ny*lpar, 0.0d0/)  !pv2(:)
sv3(:) = (/0.0d0, 0.0d0, nz*lpar/)  !pv3(:)

case(5) 

write(*,*) 'SUPERCELL: Case IBRAV = 5 not implemented!'

case(6)        ! Supercells with arbitrary shape

sv1(:) = nx*pv1(:)
sv2(:) = ny*pv2(:)
sv3(:) = nz*pv3(:)

case default 

write(*,*) 'SUPERCELL: Case IBRAV = 5 not implemented!'

end select

! -----------------------------------------------------------------------------
! Loading the arrays that control the cases for the overlap computation:

call arrsel(nnsel,llsel)

! Initializing the arrays SR (overlap) and HR (hamiltonian):

sr(:,:) = 0.0d0              
hr(:,:) = 0.0d0

! Loading the main diagonal of the hamiltonian:

mdcount = 0
do i = 1, nat
 do j = 1, dbov(itype(i))
  mdcount = mdcount + 1
   hr(mdcount,mdcount) = oe(itype(i),j)
 enddo
enddo

! Building the hamiltonian:

mdcount = 0

do i = 1, nat                ! Loop over the atoms in the central cell:

ui(1) = rx(i)
ui(2) = ry(i)
ui(3) = rz(i)

do oi = 1, norbv(itype(i))   ! Loop over the I's orbitals:

do j = 1, nat                ! Loop over the J atoms

do oj = 1, norbv(itype(j))   ! Loop over the J's orbitals:

li = lval(itype(i),oi)
lj = lval(itype(j),oj)

! ------------------------------------------------------------------------
!                       Loop over the neighbor's supercells
! ------------------------------------------------------------------------

do ix = -ncell, ncell
  do iy = -ncell, ncell
    do iz = -ncell, ncell

uj(1) = rx(j) + ix*sv1(1) + iy*sv2(1) + iz*sv3(1)
uj(2) = ry(j) + ix*sv1(2) + iy*sv2(2) + iz*sv3(2)
uj(3) = rz(j) + ix*sv1(3) + iy*sv2(3) + iz*sv3(3)

! Relative position

rij = dsqrt((uj(1) - ui(1))**2 + (uj(2) - ui(2))**2 + (uj(3) - ui(3))**2)

! Excluding the case j = i in the central cell:

if(rij.eq.0.0) then 

mdcount = mdcount + 1

sr(mdcount,mdcount) = sr(mdcount,mdcount) + 1.0d0  ! Normalized orbitals  

else

! -----------------------
if(rij.lt.rcut) then

call calcovlp(i,j,oi,oj,nat,numtp,nzt,nval,lval,itype,zeta,cf,ui,uj,&
              llsel,nnsel,ovout,ipbc)

! Setting the HR and SR elements:

ct = 0
do jo = 1, ldm(lj)
 jd = sip(j) + lptb(itype(j),lj) + jo - 1
  do io = 1, ldm(li)

   id = sip(i) + lptb(itype(i),li) + io - 1
   ct = ct + 1

   kh = keht(itype(i),itype(j))
   sr(id,jd) = sr(id,jd) + ovout(ct)
   hr(id,jd) = hr(id,jd) + 0.5d0*kh*(esit(itype(i),oi) + esit(itype(j),oj))*ovout(ct)

 enddo
enddo

endif
! -----------------------

endif

  enddo   ! End of IZ loop
 enddo   ! End of IY loop
enddo   ! End of IX loop
! ------------------------------------------------------------------------

enddo         ! End of the loop over the J's orbitals

enddo         ! End of the Loop over the J atoms

enddo         ! End of the loop over the I's orbitals

enddo         ! End of the loop over the I atoms

! ------------------------------------------------------------------------
!                    DIAGONALIZATION OF THE HAMILTONIAN
! ------------------------------------------------------------------------

! Apllying the shift on the energies of the hamiltonian:

hr(:,:) = hr(:,:) + eshift*sr(:,:)

! Setting the parameters for the lapack diagonalization subroutine:

itp = 1
jobz = 'V'
uplo = 'U'
ndm = tnao
lda = ndm
ldb = ndm
lwork = 1 + 6*ndm + 2*ndm*ndm 
liwork = 3 + 5*ndm
allocate(work(lwork),iwork(liwork))

call dsygvd(itp, jobz, uplo, ndm, hr, lda, sr, ldb, W, work, lwork, iwork, liwork, info)

! ------------------------------------------------------------------------
!                    NORMALIZATION OF THE EIGENSTATES  
! ------------------------------------------------------------------------

do i = 1, tnao  ! All states will be normalized

sum = 0.0d0

  do j = 1, tnao
    do k = 1, tnao
      sum = sum + hr(j,i)*hr(k,i)*sr(k,j)
    enddo
  enddo

hr(:,i) = hr(:,i)/dsqrt(sum)

enddo

! End of the subroutine

end subroutine fcspec

!-----------------------------------------------------------------------
! SUBROUTINE KPOINTS: calculates the primitive vector of the recipocral 
!                     Lattice, identifyed with IBRAV, and set the k-points
!                     For Band Structure Calculation
!-----------------------------------------------------------------------

subroutine kpoints(nkpt,ibrav,lpar,pi,kp,pv1,pv2,pv3)
implicit none

! ----------------- Variables and Parameters of the Subroutine ----------------
! Externals:
integer, intent(in) :: nkpt, ibrav       
real*8,  intent(in) :: lpar, pi, pv1(3), pv2(3), pv3(3)
real*8,  intent(out):: kp(nkpt,3)    
! Internals
real*8              :: b1(3), b2(3), b3(3), kcoord(nkpt,3), vuc(3)
integer             :: i, j, opt
!-----------------------------------------------------------------------

! Volum of the primitive cell

vuc(1) = pv1(1)*(pv2(2)*pv3(3) - pv2(3)*pv3(2)) + pv1(2)*(pv2(3)*pv3(1) - pv2(1)*pv3(3)) + &
         pv1(3)*(pv2(1)*pv3(2) - pv2(2)*pv3(1))

vuc(2) = pv2(1)*(pv3(2)*pv1(3) - pv3(3)*pv1(2)) + pv2(2)*(pv3(3)*pv1(1) - pv3(1)*pv1(3)) + &
         pv2(3)*(pv3(1)*pv1(2) - pv3(2)*pv1(1))

vuc(3) = pv3(1)*(pv1(2)*pv2(3) - pv1(3)*pv2(2)) + pv3(2)*(pv1(3)*pv2(1) - pv1(1)*pv2(3)) + &
         pv3(3)*(pv1(1)*pv2(2) - pv1(2)*pv2(1))

! Primitive vectors of the reciprocal lattice

b1(1) = 2*pi*(pv2(2)*pv3(3) - pv2(3)*pv3(2))/vuc(1)
b1(2) = 2*pi*(pv2(3)*pv3(1) - pv2(1)*pv3(3))/vuc(1)
b1(3) = 2*pi*(pv2(1)*pv3(2) - pv2(2)*pv3(1))/vuc(1)

b2(1) = 2*pi*(pv3(2)*pv1(3) - pv3(3)*pv1(2))/vuc(2)
b2(2) = 2*pi*(pv3(3)*pv1(1) - pv3(1)*pv1(3))/vuc(2)
b2(3) = 2*pi*(pv3(1)*pv1(2) - pv3(2)*pv1(1))/vuc(2)

b3(1) = 2*pi*(pv1(2)*pv2(3) - pv1(3)*pv2(2))/vuc(3)
b3(2) = 2*pi*(pv1(3)*pv2(1) - pv1(1)*pv2(3))/vuc(3)
b3(3) = 2*pi*(pv1(1)*pv2(2) - pv1(2)*pv2(1))/vuc(3)

! Selecting the coordinates of the K-Points according IBRAV

opt = ibrav

select case(opt)  

! * * * * * * * * * * * * * * * * * * FCC LATTICE * * * * * * * * * * * * * * *

case(1) 

kcoord(1,:) = (/0.00d0,0.00d0,0.00d0/)   ! Gamma 
kcoord(2,:) = (/0.50d0,0.00d0,0.50d0/)   ! X
kcoord(3,:) = (/0.50d0,0.25d0,0.75d0/)   ! W
kcoord(4,:) = (/0.50d0,0.50d0,0.50d0/)   ! L
kcoord(5,:) = (/0.00d0,0.00d0,0.00d0/)   ! Gamma
kcoord(6,:) = (/0.50d0,0.00d0,0.50d0/)   ! X

! * * * * * * * * * * * * * * * * * ZINCBLEND LATTICE * * * * * * * * * * * * *

case(2) 

kcoord(1,:) = (/0.00d0,0.00d0,0.00d0/)   ! Gamma 
kcoord(2,:) = (/0.50d0,0.00d0,0.50d0/)   ! X
kcoord(3,:) = (/0.50d0,0.25d0,0.75d0/)   ! W
kcoord(4,:) = (/0.50d0,0.50d0,0.50d0/)   ! L
kcoord(5,:) = (/0.00d0,0.00d0,0.00d0/)   ! Gamma
kcoord(6,:) = (/0.00d0,0.50d0,0.50d0/)   ! X

! * * * * * * * * * * * * * * * * * DIAMOND LATTICE * * * * * * * * * * * * * *

case(3) 

!kcoord(1,:) = (/0.00d0,0.00d0,0.00d0/)   ! Gamma 
!kcoord(2,:) = (/0.50d0,0.00d0,0.50d0/)   ! X
!kcoord(3,:) = (/0.50d0,0.25d0,0.75d0/)   ! W
!kcoord(4,:) = (/0.50d0,0.50d0,0.50d0/)   ! L
!kcoord(5,:) = (/0.00d0,0.00d0,0.00d0/)   ! Gamma
!kcoord(6,:) = (/0.50d0,0.00d0,0.50d0/)   ! X

kcoord(1,:) = (/0.00d0,0.00d0,0.00d0/)   ! Gamma 
kcoord(2,:) = (/0.50d0,0.00d0,0.50d0/)   ! X
kcoord(3,:) = (/0.50d0,0.50d0,0.50d0/)   ! L
kcoord(4,:) = (/0.00d0,0.00d0,0.00d0/)   ! Gamma
kcoord(5,:) = (/0.50d0,0.00d0,0.50d0/)   ! X
kcoord(6,:) = (/0.50d0,0.25d0,0.75d0/)   ! W

! * * * * * * * * * * * * * * * * * BCC LATTICE * * * * * * * * * * * * * * * *

case(4) 

kcoord(1,:) = (/0.00d0, 0.00d0, 0.00d0/)   ! Gamma 
kcoord(2,:) = (/0.50d0,-0.50d0, 0.50d0/)   ! H
kcoord(3,:) = (/0.00d0, 0.00d0, 0.50d0/)   ! N
kcoord(4,:) = (/0.00d0, 0.00d0, 0.00d0/)   ! Gamma
kcoord(5,:) = (/0.25d0, 0.25d0, 0.25d0/)   ! P    
kcoord(6,:) = (/0.50d0,-0.50d0, 0.50d0/)   ! H

! * * * * * * * * * * * * * * * * * HCP LATTICE * * * * * * * * * * * * * * * *

case(5) 

kcoord(1,:) = (/0.5d0, 0.0d0, 0.0d0/)         ! M     
kcoord(2,:) = (/0.5d0, 0.0d0, 0.5d0/)         ! L
kcoord(3,:) = (/0.0d0, 0.0d0, 0.5d0/)         ! A
kcoord(4,:) = (/0.0d0, 0.0d0, 0.0d0/)         ! Gamma
kcoord(5,:) = (/1.d0/3.d0, 1.d0/3.d0, 0.0d0/) ! K    
kcoord(6,:) = (/1.d0/3.d0, 1.d0/3.d0, 0.5d0/) ! H    

! * * * * * * * * * * * * * * * * * FREE LATTICE  * * * * * * * * * * * * * * *
! In this case, assign below the points according to the system's symmetry

case(6) 

!kcoord(1,:) = (/0.5d0, 0.0d0, 0.0d0/)         ! M     
!kcoord(2,:) = (/0.5d0, 0.0d0, 0.5d0/)         ! L
!kcoord(3,:) = (/0.0d0, 0.0d0, 0.5d0/)         ! A
!kcoord(4,:) = (/0.0d0, 0.0d0, 0.0d0/)         ! Gamma
!kcoord(5,:) = (/1.d0/3.d0, 1.d0/3.d0, 0.0d0/) ! K    
!kcoord(6,:) = (/1.d0/3.d0, 1.d0/3.d0, 0.5d0/) ! H    

!kcoord(1,:) = (/0.5d0, 0.0d0, 0.0d0/)         ! M     
!kcoord(2,:) = (/2.d0/3.d0, 1.d0/3.d0, 0.0d0/) ! K    
!kcoord(3,:) = (/0.0d0, 0.0d0, 0.0d0/)         ! Gamma
!kcoord(4,:) = (/0.5d0, 0.0d0, 0.0d0/)         ! M     
!kcoord(5,:) = (/2.d0/3.d0, 1.d0/3.d0, 0.0d0/) ! K    
!kcoord(6,:) = (/0.0d0, 0.0d0, 0.0d0/)         ! Gamma

!kcoord(1,:) = (/-1.d0/3.d0, 1.d0/3.d0, 0.0d0/) ! K
!kcoord(2,:) = (/ 0.0d0, 0.0d0, 0.0d0/)         ! Gamma
!kcoord(3,:) = (/ 0.0d0, 0.5d0, 0.0d0/)         ! M
!kcoord(4,:) = (/ 0.0d0, 0.0d0, 0.0d0/)         ! Gamma
!kcoord(5,:) = (/-1.d0/3.d0, 1.d0/3.d0, 0.0d0/) ! K
!kcoord(6,:) = (/ 0.0d0, 0.5d0, 0.0d0/)         ! M

!SILICENO e CIA
kcoord(1,:) = (/ 0.5d0, 0.0d0, 0.0d0/)         ! M
kcoord(2,:) = (/ 0.0d0, 0.0d0, 0.0d0/)         ! Gamma
kcoord(3,:) = (/ 1.d0/3.d0, 1.d0/3.d0, 0.0d0/) ! K
kcoord(4,:) = (/ 0.5d0, 0.0d0, 0.0d0/)         ! M
kcoord(5,:) = (/ 0.0d0, 0.0d0, 0.0d0/)         ! Gamma
kcoord(6,:) = (/ 1.d0/3.d0, 1.d0/3.d0, 0.0d0/) ! K

case default

write(*,*) 'IBRAV =', ibrav, '->  CASE NOT RECOGNIZED'
stop

end select

do i = 1, nkpt
 kp(i,1) = kcoord(i,1)*b1(1) + kcoord(i,2)*b2(1) + kcoord(i,3)*b3(1)
  kp(i,2) = kcoord(i,1)*b1(2) + kcoord(i,2)*b2(2) + kcoord(i,3)*b3(2)
 kp(i,3) = kcoord(i,1)*b1(3) + kcoord(i,2)*b2(3) + kcoord(i,3)*b3(3)
enddo

! End of the subroutine

end subroutine kpoints

!-----------------------------------------------------------------------
! SUBROUTINE MONKPACK: Calculates the k-points according the Monkhost-Pack
!                      Recipe.
!-----------------------------------------------------------------------

subroutine monkpack(tnkp,qp,qr,qs,lpar,pi,kpoint,pv1,pv2,pv3)
implicit none

! ----------------- Variables and Parameters of the Subroutine ----------------
! Externals:
integer, intent(in) :: tnkp, qp, qr, qs
real*8,  intent(in) :: lpar, pi, pv1(3), pv2(3), pv3(3)
real*8,  intent(out):: kpoint(tnkp,3)    
! Internals
real*8              :: b1(3), b2(3), b3(3), vruc, up, ur, us
integer             :: i, j, k, ic
!------------------------------------------------------------------------------

! Volum of the reciprocal unit cell

vruc = pv1(1)*(pv2(2)*pv3(3) - pv2(3)*pv3(2)) + pv1(2)*(pv2(3)*pv3(1) - pv2(1)*pv3(3)) + &
       pv1(3)*(pv2(1)*pv3(2) - pv2(2)*pv3(1))

! Reciprocal lattice Vectors:

b1(1) = 2*pi*(pv2(2)*pv3(3) - pv2(3)*pv3(2))/vruc
b1(2) = 2*pi*(pv2(3)*pv3(1) - pv2(1)*pv3(3))/vruc
b1(3) = 2*pi*(pv2(1)*pv3(2) - pv2(2)*pv3(1))/vruc

b2(1) = 2*pi*(pv3(2)*pv1(3) - pv3(3)*pv1(2))/vruc
b2(2) = 2*pi*(pv3(3)*pv1(1) - pv3(1)*pv1(3))/vruc
b2(3) = 2*pi*(pv3(1)*pv1(2) - pv3(2)*pv1(1))/vruc

b3(1) = 2*pi*(pv1(2)*pv2(3) - pv1(3)*pv2(2))/vruc
b3(2) = 2*pi*(pv1(3)*pv2(1) - pv1(1)*pv2(3))/vruc
b3(3) = 2*pi*(pv1(1)*pv2(2) - pv1(2)*pv2(1))/vruc

! Gamma Point:

kpoint(1,:) = 0.0d0

! Generating the K-points:

ic = 2

do i = 1, qp
  do j = 1, qr
    do k = 1, qs

      up = dfloat((2*i) - qp - 1)/dfloat(2*qp)
      ur = dfloat((2*j) - qr - 1)/dfloat(2*qr)
      us = dfloat((2*k) - qs - 1)/dfloat(2*qs)

      kpoint(ic,1) = up*b1(1) + ur*b2(1) + us*b3(1)
      kpoint(ic,2) = up*b1(2) + ur*b2(2) + us*b3(2)
      kpoint(ic,3) = up*b1(3) + ur*b2(3) + us*b3(3)

      ic = ic + 1

    enddo
  enddo
enddo

! End of the subroutine

end subroutine monkpack

! -----------------------------------------------------------------------------
end module diagonalize
! -----------------------------------------------------------------------------
