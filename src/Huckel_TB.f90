! -----------------------------------------------------------------------------
!                               Program HUCKEL_TB
! Date: 27/03/2014
! -----------------------------------------------------------------------------

program huckel_tb
use sto_ov
use supercell
use diagonalize
use toten
use postproc

implicit none

!----------------------- Variables and Arrays ---------------------------------
! (i) Parameters:
integer, parameter       :: nmax = 5, lmax = 2
real*8, parameter        :: pi = 3.1415926535898d0
! (ii) From the input file:
integer                  :: nat, nx, ny, nz, naux, numtp, ibrav, kcalc, imc, ids
integer                  :: qp, qr, qs, icp1, icp2
real*8                   :: lpar, cova, zeta1, zeta2, rho, rcut
integer, allocatable     :: zat(:)
integer, allocatable     :: norbv(:), nval(:,:), lval(:,:), nzt(:,:), neao(:,:)
real*8 , allocatable     :: zeta(:,:,:), cf(:,:,:), esit(:,:), keht(:,:)
character*2, allocatable :: zsimb(:)
character*40             :: cltit, ibravopt(6), kcalcopt(5)
character*5              :: ppopt(0:1)
! (iii) Internals
real*8,  allocatable     :: rx(:), ry(:), rz(:), oe(:,:), hr(:,:), sr(:,:), W(:)
real*8                   :: rpos(2,3), cl, cm, cn, cl2, cm2, cn2, rsq, ebin, version
real*8                   :: ovlp(0:lmax), ovin(3), ovout(25), eshift, vuc
real*8                   :: xij, yij, zij, rij, sdx, sdy, sdz, pv1(3), pv2(3), pv3(3)
integer, allocatable     :: itype(:), nngh(:), lisngh(:,:), dbov(:), ldm(:), sip(:)
integer, allocatable     :: lptb(:,:), nspc(:)
integer                  :: ni, nj, li, lj, llp, i, j, k, m, n, id, jd, ik, kk, today(3)
integer                  :: nnn, tnao, ErrorFlag, nelect, natopt(6,5), imin(1), imax(1)
! (iv) Internals - Reciprocal space
integer                  :: nspbz, ndiv, nkpt, idv, ikcount, flp
real*8                   :: emin, emax, egap
real*8, allocatable      :: eband(:,:)
real*8, allocatable      :: kp(:,:), kpmod(:), kpmodacc(:), kpoint(:,:), kpcoord(:)
complex*16, allocatable  :: hamk(:,:,:), ovlk(:,:,:)
! -----------------------------------------------------------------------------

! Commons:

common /sdx/ sdx
common /sdy/ sdy
common /sdz/ sdz

version = 1.0
call idate(today)

print*
print*, "--------------------- PROGRAM HUCKEL_TB ------------------------------"
write(*,'(a,f3.1)') '     VERSION: ', version
write(*,'(a,i2,a,i2,a,i4)') 'CURRENT DATE: ', today(1),'/',today(2),'/',today(3)
print*, "----------------------------------------------------------------------"
print*
print*, "Calculates the electronic structure of a crystalline material in the  "
print*, "Tight-Binding aproximation. The radial part of the atomic orbitals are" 
print*, "Described  as linear combination  of Slater-Type Functions (STO).  The"
print*, "Elements of the Hamiltonian matrix follows the Huckel prescription: the" 
print*, "Overlaps among the orbitals are explicitely computed."
print* 
print*, "Unit System:"
print*
print*, "DISTANCES: Angstrons (A)"
print*, " ENERGIES: electron volts (eV)"
print*, "   FORCES: eV/Angs"
print*
print*, "In the current version  of the program, the maximum allowed values for"
print*, "The N and L quantum numbers are NMAX = 5 and LMAX = 2 (s, p and d) and"
print*, "The maximal number of the orbitals for each specie is 3."

! Files managed by the main program:

open(4, file ='input.dat', status='unknown')
open(8, file ='atoms.dat', status='unknown')
open(7, file ='positions.xyz', status='unknown')

! ------------------------ Reading the input file -----------------------------

read(4,*) cltit                           ! Title of the Calculation
read(4,*) kcalc                           ! Kind of Calculation (1 = REAL, 2 = BANDS and 3 = DOS)
read(4,*) ibrav                           ! Bravais Group for the positions
read(4,*) numtp                           ! Number of atom types of the system

allocate(zat(numtp), zsimb(numtp), norbv(numtp))

read(4,*) (  zat(i), i = 1, numtp)        ! Atomic numbers of the species
read(4,*) (zsimb(i), i = 1, numtp)        ! Atomic symbols of the species
read(4,*) (norbv(i), i = 1, numtp)        ! Number of Valence orbitals

allocate(nval(numtp,3), lval(numtp,3),neao(numtp,3))

do i = 1, numtp
read(4,*) (nval(i,j), j = 1, norbv(i))    ! N Quantum Numbers of the valence orbitals
enddo

do i = 1, numtp
read(4,*) (lval(i,j), j = 1, norbv(i))    ! L Quantum Numbers of the valence orbitals
enddo

do i = 1, numtp
read(4,*) (neao(i,j), j = 1, norbv(i))    ! Number of Electrons of the J atomic Orbital
enddo

read(4,*) lpar, cova
read(4,*) nx, ny, nz
read(4,*) rcut

! Allocating KEHT:

allocate(keht(numtp,numtp))
read(4,*) (keht(i,i), i = 1, numtp)

! Building the cross terms of KEHT (geometric average):

do i = 1, numtp - 1
  do j = i + 1, numtp
    keht(i,j) = dsqrt(keht(i,i)*keht(j,j))    !0.50d0*(keht(i,i) + keht(j,j))
    keht(j,i) = keht(i,j)
  enddo
enddo

! Allocating and reading the informations about the STO's: 

allocate(nzt(numtp,3))    ! 1 = s, 2 = p and 3 = d
nzt(:,:) = 0

! Reading the number of STO for each basis orbital:

do i = 1, numtp
read(4,*) (nzt(i,j), j = 1, norbv(i))
enddo

! Reading the ZETA(NUMTP,3,2) and the COEFFICIENTS CF(NUMTP,3,2) for the STO's
! (3 refers to the maximal number of valence orbitals and 2 the maximal number
! Of ZETAS.):

allocate(zeta(numtp,3,2), cf(numtp,3,2))

zeta(:,:,:) = 0.0d0
  cf(:,:,:) = 0.0d0

do i = 1, numtp
  do j = 1, norbv(i)
     do k = 1, nzt(i,j)
       read(4,*) zeta(i,j,k), cf(i,j,k)
     enddo
  enddo
enddo

! Reading the on-site energies, ESIT(NUMTP,3), of the valence orbitals:

allocate(esit(numtp,3))
esit(:,:) = 0.0d0

do i = 1, numtp
read(4,*) (esit(i,j), j = 1, norbv(i))
enddo

! Reading the shift tho the hamiltonian values (to set up the fermi energy):

read(4,*) eshift

! Reanding the integers that control the post-processing part:

read(4,*) imc, ids

! Parameters needed to the post-processing:

read(4,*) qp, qr, qs        ! Monkhost-Pack Grid
read(4,*) icp1, icp2        ! Basis orbital index for COOP

! ------------------- SETTING THE INITIAL VARIABLES -----------------------------

! Defining the vectors related to the basis orbitals: LDM(0:lmax) gives the
! Number of orbitals for each L value and DBOV(NUMTP) gives the total number of
! Orbitals for each specie of the system:

allocate(ldm(0:lmax), dbov(numtp))

ldm(0) = 1
ldm(1) = 3
ldm(2) = 5

! Assigning the total number of valence orbitals according to the basis dimension:

do i = 1, numtp
 do j = 1, norbv(i)
  dbov(i) = dbov(i) + ldm(lval(i,j))
 enddo
enddo

! Pointer for each L value inside the basis: for each specie, this pointer
! Gives the initial position in the basis for the orbital of L quantum number.

allocate(lptb(numtp,0:lmax))

do i = 1, numtp

llp = norbv(i)
select case (llp)
case(1) ! Basis with only S orbital
lptb(i,0) = 1
case(2) ! Basis with S and P orbitals
lptb(i,0) = 1
lptb(i,1) = 2
case(3) ! Basis with S, P and D orbitals
lptb(i,0) = 1
lptb(i,1) = 2
lptb(i,2) = 5
end select

enddo

! Post-Proceeding Options:

ppopt(0) = '   NO'
ppopt(1) = '  YES'

! ------------------- GENERATION OF THE ATOMIC POSITIONS ---------------------
! For the positions generation we need two input parameters: KCALC and IBRAV.
! If KCALC = 1, the program work with the minimum unit cell. On the other hand,
! If KCALC = 1, the program generates a supercell according the number of units
! Cells Nx, Ny and Nz, readen from the input file. We have the following options
! For the bravais Lattices: ibrav = 1 (FCC), ibrav = 2 (ZINC), ibrav = 3 (DIAMOND),
! Ibrav = 4 (BCC), ibrav = 5 (HCP) or give the atomic positions (FREE):

ibravopt(1) = 'fcc'
ibravopt(2) = 'zinc'
ibravopt(3) = 'sili'
ibravopt(4) = 'bcc'
ibravopt(5) = 'hcp'
ibravopt(6) = 'free'

pv1(:) = 0.0d0
pv2(:) = 0.0d0
pv3(:) = 0.0d0

! Kind of calculations:

kcalcopt(1) = 'Energy Spectrum in Real Space'
kcalcopt(2) = 'Energy Bands in Reciprocal Space'
kcalcopt(3) = 'Density of States Calculation'
kcalcopt(4) = 'Total Energy Calculation - Real Space'
kcalcopt(5) = 'Total Energy Calculation - Reciprocal Space'

! Assigning the number of atoms natopt(IBRAV,KCALC):

natopt(:,:) = 0
natopt(1,1) = 4*nx*ny*nz
natopt(1,2) = 1
natopt(1,3) = 1
natopt(1,4) = 4*nx*ny*nz
natopt(1,5) = 1
natopt(2,1) = 8*nx*ny*nz
natopt(2,2) = 2
natopt(2,3) = 2
natopt(2,4) = 8*nx*ny*nz
natopt(2,5) = 2
natopt(3,1) = 8*nx*ny*nz
natopt(3,2) = 2
natopt(3,3) = 2
natopt(3,4) = 8*nx*ny*nz
natopt(3,5) = 2
natopt(4,1) = 2*nx*ny*nz
natopt(4,2) = 1
natopt(4,3) = 1
natopt(4,4) = 2*nx*ny*nz
natopt(4,5) = 1
natopt(5,1) = 2*nx*ny*nz
natopt(5,2) = 2
natopt(5,3) = 2
natopt(5,4) = 2*nx*ny*nz
natopt(5,5) = 2

if(ibrav.eq.6) then   ! Read the first lines of ATOMS.DAT, the number 
                      ! Of atoms in the unit cell and primitive vectors
 read(8,*) naux
 read(8,*) pv1(1), pv1(2), pv1(3)
 read(8,*) pv2(1), pv2(2), pv2(3)
 read(8,*) pv3(1), pv3(2), pv3(3)
 close(8)
natopt(6,1) = naux*nx*ny*nz
natopt(6,2) = naux
natopt(6,3) = naux

! Scaling the lattice Primitive Vectors:

pv1(:) = pv1(:)*lpar
pv2(:) = pv2(:)*lpar
pv3(:) = pv3(:)*lpar*cova

endif

! Number of Atoms:

nat = natopt(ibrav,kcalc)

allocate(rx(nat), ry(nat), rz(nat), itype(nat))

! Generating the positions and atributing the species:

call gerapos(nat,ibrav,kcalc,nx,ny,nz,itype,lpar,cova,rx,ry,rz,pv1,pv2,pv3)

! Printing the positions in a output file (XYZ format):

write(7,*) nat
write(7,*) 
do i = 1, nat
!write(7,'(a,3f14.5)') zsimb(itype(i)), rx(i), ry(i), rz(i)
write(7,'(i2,3f14.5)') itype(i), rx(i), ry(i), rz(i)
enddo

! Dimensions of the supercell (needed for KCALC = 1)

if(ibrav.eq.6) then
 sdx = nx*dsqrt(pv1(1)*pv1(1) + pv1(2)*pv1(2) + pv1(3)*pv1(3))  
  sdy = ny*dsqrt(pv2(1)*pv2(1) + pv2(2)*pv2(2) + pv2(3)*pv2(3))
   sdz = nz*dsqrt(pv3(1)*pv3(1) + pv3(2)*pv3(2) + pv3(3)*pv3(3))
    else
   sdx = nx*lpar
  sdy = ny*lpar
 sdz = nz*lpar*cova
endif

! Volum of the primitive cell

vuc = pv1(1)*(pv2(2)*pv3(3) - pv2(3)*pv3(2)) + pv1(2)*(pv2(3)*pv3(1) - pv2(1)*pv3(3)) + &
      pv1(3)*(pv2(1)*pv3(2) - pv2(2)*pv3(1))

vuc = dabs(vuc)   ! For non-orthogonal primitive vectors

! ------------------- CALCULATION OF THE NEIGHBOR LIST ------------------------
! Estimating  the number of neighbors for each atom: for systems with cristal
! Simmetry, all the sites are equivalent and the number of neighbohrs will be
! Proportional to the relative volume:

nnn = 0

if(kcalc.eq.1) then

nnn = dnint(1.0d0*nat) 

! Allocating the array that stores the information about the neighbors for each
! Atom and calling the subroutine NEIGHBORS:

allocate(nngh(nat), lisngh(nat,nnn))
call neighbors(ibrav,nat,nnn,nx,ny,nz,lpar,pv1,pv2,pv3,rx,ry,rz,rcut,nngh,lisngh)

endif

! -----------------------------------------------------------------------------
!              Calculating the total number of electrons NELECT:

nelect = 0

do i = 1, nat
  do k = 1, norbv(itype(i))
    nelect = nelect + neao(itype(i),k)
  enddo
enddo

! -----------------------------------------------------------------------------
!             Printing the information readed from the input file:

print*
print*, "----------------------------------------------------------------------"
print*, "          SUMMARY OF THE DATA READ FROM THE INPUT FILE                "
print*, "----------------------------------------------------------------------"
print*
print*, "(i) System Parameters:"
print*
write(*,'(a,a)')    ' System                                    =  ',cltit
write(*,'(a,a)')    ' Bravais Lattice                           =  ',ibravopt(ibrav)
write(*,'(a,a)')    ' Kind of Calculation                       =  ',kcalcopt(kcalc)
write(*,'(a,i8)')   ' Number of atoms (NAT)                     =', nat
write(*,'(a,i8)')   ' Number of atom types (NUMTP)              =', numtp
write(*,'(a,i8)')   ' Number of Neighbors  (NNN)                =', nnn  
write(*,'(a,3i4)')  ' Supercell Dimensions (Nx,Ny,Nz)           =', nx, ny, nz
write(*,'(a,f8.3)') ' Lattice Parameter (LPAR)                  =', lpar 
write(*,'(a,f8.3)') ' C/A Ratio (COVA)                          =', cova 
print*
print*, "(ii) Structural Supercell Parameters:"
print*
write(*,'(a,f8.3)') ' Supercell Dimension Along X (SDX)         =', sdx  
write(*,'(a,f8.3)') ' Supercell Dimension Along Y (SDY)         =', sdy  
write(*,'(a,f8.3)') ' Supercell Dimension Along Z (SDZ)         =', sdz  
write(*,'(a,f8.3)') ' Unit Cell Volum (VUC)                     =', vuc  
print*
print*, "(iii) Information about the atom types:"
print*
write(*,'(a)')   '(a) Atomic numbers (ZAT) and Symbols:          '
print*
do k = 1, numtp
write(*,'(a,a,i1,a,a,i2,a,a,a)') ' ZAT','(',k,')',' = ',zat(k),' (',zsimb(k),')'
enddo
print*
write(*,'(a)')   '(b) Valence orbitals for each atom type:       '
print*
do k = 1, numtp
write(*,'(a,a,i1,a,a,i2,a,a,a)') ' NORBV','(',k,')',' = ', norbv(k)
enddo
print*
write(*,'(a)')   '(c) Number of Electrons of the Atomic Orbitals:       '
print*
do k = 1, numtp
write(*,'(a,a,a,a,a,a,a,a,a)') ' | ' , 'Specie' ,   ' | ' ,  '  (n,l)', '  | ', 'NEAO', ' |  ',  'ESIT', '   |  '
print*
  do j = 1, norbv(k)
write(*,'(i7,a,i1,a,i1,a,i7,f10.3)') k, '       (',nval(k,j),',',lval(k,j),')', neao(k,j), esit(k,j)
     enddo
print*
enddo
write(*,'(a,i6)') ' Total Number of Electrons: ', nelect
print*
write(*,'(a)')      '(d) Informations about the STOs for each specie:'
do i = 1, numtp
print*
write(*,'(a,a,a)') ' (',zsimb(i),'):'
print*
write(*,'(a,a,a,a,a,a,a,a)') ' | ' , 'Specie' ,   ' | ' ,  ' (n,l)', ' | ', 'ZETA', '  | ', ' CJ '
print*
  do j = 1, norbv(i)
     do k = 1, nzt(i,j)
write(*,'(i7,a,i1,a,i1,a,2f8.3)') i, '      (',nval(i,j),',',lval(i,j),')', zeta(i,j,k), cf(i,j,k)
     enddo
  enddo
enddo
print*
print*, "(iv) Post-Processing Options :"
print*
write(*,'(a,a)') ' Computation of the Mulliken Charges          =', ppopt(imc)
write(*,'(a,a)') ' Computation of the Density of States         =', ppopt(ids)
write(*,'(a,i2,a,i2,a,i2)') ' Monkhost-Pack Grid                           = ', qp,' x ',qr,' x ',qs
print*

! -----------------------------------------------------------------------------
! Counting the total number of atoms of each specie:

allocate(nspc(numtp))
nspc(:) = 0

do i = 1, nat
j = itype(i)
nspc(j) = nspc(j) + 1
enddo

! Total number of atomic orbitals:

tnao = 0
do i = 1, numtp
tnao = tnao + nspc(i)*dbov(i)
enddo

write(*,*) 'TNAO =', tnao

! Creating the SIP pointer, that gives the number of atomic orbitals up the I
! Site.

allocate(sip(nat))

sip(1) = 0   

do i = 2, nat
sip(i) = sip(i-1) + dbov(itype(i-1))
enddo

! Allocating the orbital energies according the specie:

allocate(oe(numtp,9))

do i = 1, numtp
 do j = 1, norbv(i)
  do k = 1, ldm(lval(i,j))
   id = lptb(i,lval(i,j)) + k - 1
   oe(i,id) = esit(i,j)
  enddo
 enddo
enddo

! -----------------------------------------------------------------------------
!                             K-POINT GENERATION                               
! -----------------------------------------------------------------------------

select case(kcalc)

!**********************************************************************
case(2)   ! Generating the K-Points according with the Bravais Lattice:
!**********************************************************************

nspbz = 6
 ndiv = 20

allocate(kp(nspbz,3),kpmod(nspbz),kpmodacc(nspbz))

! Reading the special points in the BZ for band structure calculation:

call kpoints(nspbz,ibrav,lpar,pi,kp,pv1,pv2,pv3)

! Modules of the difference (kp(i+1) - kp(i)):

kpmodacc(1) = 0.0d0

do i = 1, nspbz - 1
 kpmod(i) = dsqrt((kp(i+1,1) - kp(i,1))**2 + (kp(i+1,2) - kp(i,2))**2 + (kp(i+1,3) - kp(i,3))**2)
  kpmodacc(i+1) = kpmodacc(i) + kpmod(i)
enddo

! Defining the total number of K-points:

nkpt = (nspbz - 1)*ndiv

! Allocating and generating the KPOINTS:

allocate(kpoint(nkpt,3),kpcoord(nkpt))

ikcount = 1 
do ik = 1, nspbz - 1
  do idv = 0, ndiv - 1
    kpoint(ikcount,1) = kp(ik,1) + idv*(kp(ik+1,1) - kp(ik,1))/dfloat(ndiv)
      kpoint(ikcount,2) = kp(ik,2) + idv*(kp(ik+1,2) - kp(ik,2))/dfloat(ndiv)
    kpoint(ikcount,3) = kp(ik,3) + idv*(kp(ik+1,3) - kp(ik,3))/dfloat(ndiv)

   kpcoord(ikcount) = kpmodacc(ik) + idv*kpmod(ik)/dfloat(ndiv)
   ikcount = ikcount + 1

  enddo
enddo

!**********************************************************************
case(3,5)   ! Generating a Monkhost-Pack K-Point grid:
!**********************************************************************

! Number of K-points of the Monkhost-Pack grid:

nkpt = 1 + qp*qr*qs

! Allocating KPOINT:

allocate(kpoint(nkpt,3))

call monkpack(nkpt,qp,qr,qs,lpar,pi,kpoint,pv1,pv2,pv3)

end select

! ---------------------- CALCULATION OF THE EIGENSTATES ------------------------
! The difference between KCALC = 2 and KCALC = 3 options is the k-point 
! Generation: in the former case we generate the k-points along the lines that
! Connect the high-symmetry points of the Brillouin zone and, in the later
! Case, the k-kpoints are generated according the Monkhost-Pack recipe.
! -----------------------------------------------------------------------------

select case(kcalc)

! * * * * * * * * * * * 
       case(1) 
! * * * * * * * * * * *

write(*,'(a)') 'KCALC = 1: Real Space Supercell Calculation'
print*

! Opening the output file for the spectrum:

open(10,file='spectrum.dat', status='unknown')

! Allocating the hamiltonian and the overlap matrices and the eingenvals arrays:

allocate(hr(tnao,tnao), sr(tnao,tnao), W(tnao))

if(ibrav.eq.6) then
call fcspec(ibrav,numtp,nat,nx,ny,nz,tnao,neao,norbv,itype,rcut,rx,ry,rz,pv1,pv2,pv3,oe,&
            esit,keht,cf,zeta,ldm,sip,dbov,lptb,nzt,nval,lval,lpar,eshift,nelect,hr,sr,W)
else
call spectrum(numtp,nat,tnao,neao,norbv,itype,nnn,nngh,lisngh,rx,ry,rz,oe,esit,keht,cf,&
              zeta,ldm,sip,dbov,lptb,nzt,nval,lval,lpar,eshift,sdx,sdy,sdz,nelect,hr,sr,W)
endif

! Writing the eigenvalues:

do i = 1, tnao
write(10,'(i8,f12.5)') i, W(i)
enddo

! Post-processing analysis:

if((imc.eq.1).or.(ids.eq.1)) then
call pp_real(nat,tnao,numtp,ibrav,neao,nelect,norbv,itype,dbov,W,hr,sr,imc,ids)
endif

! Deallocations:

deallocate(hr,sr,W)

! * * * * * * * * * * *
     case(2)
! * * * * * * * * * * *

write(*,'(a)') 'KCALC = 2: K-Space Calculation With Minimum Unit Cell'
print*

! Allocating EBAND, HAMK and OVLK:

allocate(eband(nkpt,tnao),hamk(nkpt,tnao,tnao),ovlk(nkpt,tnao,tnao))

 eband(:,:) = 0.0d0
hamk(:,:,:) = cmplx(0.0d0,0.0d0)
ovlk(:,:,:) = cmplx(0.0d0,0.0d0)

! Band-Structure Calculation:

call bands(kcalc,ibrav,numtp,nat,tnao,norbv,itype,rx,ry,rz,pv1,pv2,pv3,oe,esit,keht,cf,zeta,&
           sip,ldm,dbov,lptb,nzt,nval,lval,lpar,eshift,nelect,nkpt,kpoint,hamk,ovlk,eband,rcut)

! File for the Band Structure and determining the Fermi Level Position:

open(9, file='Bands.dat', status ='unknown')

if(mod(nelect,2).eq.0) then
 flp = nelect/2
  else
 flp = (nelect + 1)/2
endif

! Writing the bands:

do ik = 1, nkpt
  write(9,'(f5.2,7f10.4)') kpcoord(ik), (eband(ik,j), j = flp - 3, flp + 3)
enddo

! Printing the GAP and the K-point at the minimum:

imax = maxloc(eband(:,flp))
imin = minloc(eband(:,flp+1))
emin = maxval(eband(:,flp))
emax = minval(eband(:,flp+1))
egap = emax - emin

write(*,'(a,1f10.4)') 'Gap Energy                 =', egap
write(*,'(a,1f10.4)') 'Fermi Energy               =', eband(imax,flp)
write(*,'(a,1i6   )') 'IMAX =                     =', imax
write(*,'(a,3f10.4)') 'K at the Maximum of the VB =', kpoint(imax,1), kpoint(imax,2), kpoint(imax,3)
write(*,'(a,3f10.4)') 'K at the Minimum of the CB =', kpoint(imin,1), kpoint(imin,2), kpoint(imin,3)
print* 

! * * * * * * * * * * *
     case(3)
! * * * * * * * * * * *

write(*,'(a)') 'KCALC = 3: Calculation of the Total Density of States'
print*

! Allocating EBAND, HAMK and OVLK:

allocate(eband(nkpt,tnao),hamk(nkpt,tnao,tnao),ovlk(nkpt,tnao,tnao))

 eband(:,:) = 0.0d0
hamk(:,:,:) = cmplx(0.0d0,0.0d0)
ovlk(:,:,:) = cmplx(0.0d0,0.0d0)

! Band-Structure Calculation:

call bands(kcalc,ibrav,numtp,nat,tnao,norbv,itype,rx,ry,rz,pv1,pv2,pv3,oe,esit,keht,cf,zeta,&
           sip,ldm,dbov,lptb,nzt,nval,lval,lpar,eshift,nelect,nkpt,kpoint,hamk,ovlk,eband,rcut)

! Fermi Level Position and DOS calculation:  

if(mod(nelect,2).eq.0) then
 flp = nelect/2
  else
 flp = (nelect + 1)/2
endif

call pp_recp(numtp,nat,tnao,nelect,nkpt,itype,dbov,eband,hamk,ovlk,flp,imc,ids,icp1,icp2)

! * * * * * * * * * * *
     case(4)
! * * * * * * * * * * *

write(*,'(a)') 'KCALC = 4: Total Energy Calculation in the Real Space'
print*

call tbte(kcalc,ibrav,numtp,nat,tnao,norbv,itype,neao,nx,ny,nz,rcut,rx,ry,rz,pv1,pv2,pv3,oe,&
          esit,keht,cf,zeta,sip,ldm,dbov,lptb,nzt,nval,lval,lpar,eshift,nelect,nkpt,kpoint)

! * * * * * * * * * * *
     case(5)
! * * * * * * * * * * *

write(*,'(a)') 'KCALC = 5: Total Energy Calculation in the Reciprocal Space'
print*

call tbte(kcalc,ibrav,numtp,nat,tnao,norbv,itype,neao,nx,ny,nz,rcut,rx,ry,rz,pv1,pv2,pv3,oe,&
          esit,keht,cf,zeta,sip,ldm,dbov,lptb,nzt,nval,lval,lpar,eshift,nelect,nkpt,kpoint)

! * * * * * * * * * *
     case default
! * * * * * * * * * *

write(*,'(a)') 'KCALC value not recognized'

end select

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

! Deallocations:

deallocate(zat, zsimb, norbv, nval, lval, neao, nzt, zeta, cf, esit, oe, ldm, dbov, lptb)
deallocate(rx, ry, rz, itype)

! -----------------------------------------------------------------------------
end program huckel_tb
! -----------------------------------------------------------------------------
