! -----------------------------------------------------------------------------
! MODULE SUPERCELL: In this module we have the subroutines that generates the
! Atomic positions for diffetent bravais symmetries. The subroutines of the
! Module are
!
! (1) GERAPOS: Main subroutine that, according with the readen KCALC and
!              IBRAV parameters, build the atomic positions and species.
!  
! (2)     FCC: Supercell with the atoms positioned at the Face-centered Cubic
!              Lattice sites.
! 
! (3)    ZINC: Supercell with the atoms positioned at the Zinc-Blend structure,
!              with one type of atom in the basis (silicon and germanium)
! 
! (4)    SILI: In this case the atomic positions are calculated using subroutine
!              ZINC, but the atom types are setted to just one specie (Si, Ge, ...)
!
! (5) BCCPOS: Supercell with atoms positioned in the BCC structure
!
! (6) HCPPOS: Supercell with atoms positioned in the HCP structure
! 
! (7) FREEPOS: Generates a arbitrary supercell from the basis atoms and primitive
!              Vectors readen from the ATOMS.DAT input file.
!
! (8) NEIGHBORHS: Calculates, for each atom, their neighbor list
!
! LAST MODIFICATION: 21/08/2012
! -----------------------------------------------------------------------------

module supercell
implicit none

public:: gerapos, fccpos , zincpos, bccpos, hcppos, freepos, neighbors 

CONTAINS

! -----------------------------------------------------------------------------
! Main Subroutine for generating the atomic positions and species                
! -----------------------------------------------------------------------------

subroutine gerapos(nat,ibrav,kcalc,nx,ny,nz,itype,lpar,cova,rx,ry,rz,pv1,pv2,pv3)
implicit none

! ----------------- Variables and Parameters of the Subroutine ----------------
! External
integer, intent(in)  :: nat, ibrav, kcalc, nx, ny, nz
real*8,  intent(in)  :: lpar, cova
real*8,intent(inout) :: pv1(3), pv2(3), pv3(3)
real*8               :: rx(nat), ry(nat), rz(nat)
integer              :: itype(nat)
! Internal
real*8               :: cp
integer              :: i, nasp
! ------------------------------------------------------------------------

select case(ibrav)

! ***************************
case(1)   ! FCC LATTICE
! ***************************

! Primitive Vectors:

pv1 = (/0.0d0, lpar*0.5d0, lpar*0.5d0 /)
pv2 = (/lpar*0.5d0, 0.0d0, lpar*0.5d0 /)
pv3 = (/lpar*0.5d0, lpar*0.5d0, 0.0d0 /)

select case(kcalc) 

case(1,4)

! SUPERCELL

call fccpos(nat,nx,ny,nz,itype,lpar,rx,ry,rz)

case(2,3,5)

! 1-ATOM UNIT CELL:

! Positions of the atom:

rx(1) = 0.0d0
ry(1) = 0.0d0
rz(1) = 0.0d0

! Atomic Specie:

itype(1) = 1

end select

! *****************************
case(2)   ! ZINCBLEND LATTICE
! *****************************

! Primitive Vectors:

pv1 = (/0.0d0, lpar*0.5d0, lpar*0.5d0 /)
pv2 = (/lpar*0.5d0, 0.0d0, lpar*0.5d0 /)
pv3 = (/lpar*0.5d0, lpar*0.5d0, 0.0d0 /)

select case(kcalc) 

case(1,4)

! SUPERCELL

call zincpos(nat,nx,ny,nz,itype,lpar,rx,ry,rz)

case(2,3,5)

! 2-ATOMS UNIT CELL:

! Positions of the Basis atoms:

rx(1) = 0.0d0
ry(1) = 0.0d0
rz(1) = 0.0d0
rx(2) = lpar/4.0d0
ry(2) = lpar/4.0d0
rz(2) = lpar/4.0d0

! Atomic Species:

itype(1) = 1
itype(2) = 2

end select

! ******************************
case(3)  ! DIAMOND-TYPE LATTICE 
! ******************************

! Primitive Vectors:

pv1 = (/0.0d0, lpar*0.5d0, lpar*0.5d0 /)
pv2 = (/lpar*0.5d0, 0.0d0, lpar*0.5d0 /)
pv3 = (/lpar*0.5d0, lpar*0.5d0, 0.0d0 /)

select case(kcalc) 

case(1,4)

! SUPERCELL

call zincpos(nat,nx,ny,nz,itype,lpar,rx,ry,rz)
itype(:) = 1

case(2,3,5)

! 2-ATOMS UNIT CELL:

! Positions of the Basis atoms:

rx(1) = 0.0d0
ry(1) = 0.0d0
rz(1) = 0.0d0
rx(2) = lpar/4.0d0
ry(2) = lpar/4.0d0
rz(2) = lpar/4.0d0

! Atomic Species:

itype(1) = 1
itype(2) = 1

end select

! *********************************
case(4)  ! BCC LATTICE
! *********************************

! Primitive Vectors:

pv1 = (/-lpar*0.5d0, lpar*0.5d0, lpar*0.5d0 /)
pv2 = (/ lpar*0.5d0,-lpar*0.5d0, lpar*0.5d0 /)
pv3 = (/ lpar*0.5d0, lpar*0.5d0,-lpar*0.5d0 /)

select case(kcalc) 

case(1,4)

! SUPERCELL

call bccpos(nat,nx,ny,nz,itype,lpar,rx,ry,rz)
itype(:) = 1

case(2,3,5)

! 1-ATOM UNIT CELL:

! Positions of the Basis atoms:

rx(1) = 0.0d0
ry(1) = 0.0d0
rz(1) = 0.0d0

! Atomic Species:

itype(1) = 1

end select

! *********************************
case(5)  ! HCP Lattice
! *********************************

! Primitive Vectors:

pv1 = (/ lpar*0.5d0,-dsqrt(3.d0)*lpar/2.d0, 0.0d0 /)
pv2 = (/ lpar*0.5d0, dsqrt(3.d0)*lpar/2.d0, 0.0d0 /)
pv3 = (/ 0.0d0, 0.0d0, lpar*cova /)

select case(kcalc) 

case(1,4)

! SUPERCELL

call hcppos(nat,nx,ny,nz,itype,lpar,cova,rx,ry,rz)
itype(:) = 1

case(2,3,5)

! 2-ATOMS UNIT CELL:

! Positions of the Basis atoms:

rx(1) = 0.5d0*lpar
ry(1) = (dsqrt(3.d0)/6.d0)*lpar
rz(1) = lpar*cova/4.d0
rx(2) = 0.5d0*lpar
ry(2) =-(dsqrt(3.d0)/6.d0)*lpar
rz(2) = 3.d0*lpar*cova/4.d0

! Atomic Species:

itype(1) = 1
itype(2) = 1

end select


! *********************************
case(6)
! *********************************

select case(kcalc) 

case(1,4)

! SUPERCELL

call freepos(nat,nx,ny,nz,itype,rx,ry,rz,lpar,cova)

case(2,3,5)

! UNIT CELL (Read from ATOMS.DAT input file):

open(8, file ='atoms.dat', status='unknown')
read(8,*) nasp

! Reading the primitive vectors:

 read(8,*)  pv1(1), pv1(2), pv1(3)
 read(8,*)  pv2(1), pv2(2), pv2(3)
 read(8,*)  pv3(1), pv3(2), pv3(3)

! Reading the types and the coordinates of the basis vectors:

do i = 1, nasp
read(8,*) itype(i), rx(i), ry(i), rz(i)
enddo

close(8)

! Scaling the lattice Primitive Vectors:

pv1(:) = pv1(:)*lpar
pv2(:) = pv2(:)*lpar
pv3(:) = pv3(:)*lpar*cova

end select

! *********************************
case default
! *********************************

write(*,'(a)') 'INITIAL POSITIONS -> IBRAV value not recognized!!'
stop

end select

! End of the subroutine

end subroutine gerapos

! -----------------------------------------------------------------------------
! Generates a supercell with the atomic positions in the FCC lattice                
! -----------------------------------------------------------------------------

subroutine fccpos(nat,nx,ny,nz,itype,lpar,rx,ry,rz)
implicit none

! ----------------- Variables and Parameters of the Subroutine ----------------
! External
integer,   intent(in) :: nx, ny, nz, nat
real*8,    intent(in) :: lpar
real*8, intent(inout) :: rx(nat), ry(nat), rz(nat)
integer,intent(inout) :: itype(nat)
! Calculated in the subroutine 
real*8   :: sdx, sdy, sdz
! Internals
real*8  ::  base(4,3), itp(4)
integer :: iatom, ix, iy, iz, ib, jk
! ------------------------------------------------------------------------

sdx = nx*lpar
sdy = ny*lpar
sdz = nz*lpar

! Basis vectors: 

base(1,1) = 0.0d0
base(1,2) = 0.0d0
base(1,3) = 0.0d0
base(2,1) = 0.0d0
base(2,2) = (0.5d0)*lpar
base(2,3) = (0.5d0)*lpar
base(3,1) = (0.5d0)*lpar
base(3,2) = 0.0d0
base(3,3) = (0.5d0)*lpar
base(4,1) = (0.5d0)*lpar
base(4,2) = (0.5d0)*lpar
base(4,3) = 0.0d0

! Types of the basis atoms (change here in the case of more types in the basis):

itp(:) = 1

! Loop over the atoms:

iatom = 1

do ix = 0, nx - 1
  do iy = 0, ny - 1
    do iz = 0, nz - 1

! Loop over the basis vectors:

do ib = 1, 4
 rx(iatom) = ix*lpar + base(ib,1)
  ry(iatom) = iy*lpar + base(ib,2)
 rz(iatom) = iz*lpar + base(ib,3)
 itype(iatom) = itp(ib)
iatom = iatom + 1
enddo

! End of the loop over the basis vectors

    enddo
  enddo
enddo

! Shifting the center of the supercell to the (0,0,0) position            

do jk = 1, nat
 rx(jk) = rx(jk) - 0.5d0*sdx
  ry(jk) = ry(jk) - 0.5d0*sdy
 rz(jk) = rz(jk) - 0.5d0*sdz
enddo

! End of the subroutine

end subroutine fccpos

! -----------------------------------------------------------------------------
! Generates a supercell with the atomic positions in the zincblend structure,
! With TWO different atom types in the basis.
! -----------------------------------------------------------------------------

subroutine zincpos(nat,nx,ny,nz,itype,lpar,rx,ry,rz)
implicit none

! ----------------- Variables and Parameters of the Subroutine ----------------
! External
integer,   intent(in) :: nx, ny, nz, nat
real*8,    intent(in) :: lpar
real*8, intent(inout) :: rx(nat), ry(nat), rz(nat)
integer,intent(inout) :: itype(nat)
! Calculated in the subroutine 
real*8   :: sdx, sdy, sdz
! Internals
real*8  ::  base(4,3)
integer :: iatom, ix, iy, iz, ib, jk
! -----------------------------------------------------------------------------

sdx = nx*lpar
sdy = ny*lpar
sdz = nz*lpar

! Basis vectors: 

base(1,1) = 0.0d0
base(1,2) = 0.0d0
base(1,3) = 0.0d0
base(2,1) = 0.0d0
base(2,2) = (0.5d0)*lpar
base(2,3) = (0.5d0)*lpar
base(3,1) = (0.5d0)*lpar
base(3,2) = 0.0d0
base(3,3) = (0.5d0)*lpar
base(4,1) = (0.5d0)*lpar
base(4,2) = (0.5d0)*lpar
base(4,3) = 0.0d0

! Loop over the atoms:

iatom = 1

do ix = 0, nx - 1
  do iy = 0, ny - 1
    do iz = 0, nz - 1

! Loop over the basis vectors:

do ib = 1, 4

! Alpha sub-lattice:

 rx(iatom) = ix*lpar + base(ib,1)
  ry(iatom) = iy*lpar + base(ib,2)
 rz(iatom) = iz*lpar + base(ib,3)
itype(iatom) = 1

! Beta sub-lattice:

 rx(iatom + 1) = rx(iatom) + lpar/4. 
  ry(iatom + 1) =  ry(iatom) + lpar/4.
 rz(iatom + 1) = rz(iatom) + lpar/4.
itype(iatom + 1) = 2

iatom = iatom + 2

! End of the loop over the basis vectors

enddo

    enddo
  enddo
enddo

! Shifting the center of the supercell to the (0,0,0) position            

!do jk = 1, nat
! rx(jk) = rx(jk) - 0.5d0*sdx
!  ry(jk) = ry(jk) - 0.5d0*sdy
! rz(jk) = rz(jk) - 0.5d0*sdz
!enddo

! End of the ZINCPOS subroutine

end subroutine zincpos

! -----------------------------------------------------------------------------
! Generates a supercell with the atomic positions in the BCC lattice                
! -----------------------------------------------------------------------------

subroutine bccpos(nat,nx,ny,nz,itype,lpar,rx,ry,rz)
implicit none

! ----------------- Variables and Parameters of the Subroutine ----------------
! External
integer,   intent(in) :: nx, ny, nz, nat
real*8,    intent(in) :: lpar
real*8, intent(inout) :: rx(nat), ry(nat), rz(nat)
integer,intent(inout) :: itype(nat)
! Calculated in the subroutine 
real*8   :: sdx, sdy, sdz
! Internals
real*8  ::  base(2,3), itp(2)
integer :: iatom, ix, iy, iz, ib, jk
! ------------------------------------------------------------------------

sdx = nx*lpar
sdy = ny*lpar
sdz = nz*lpar

! Basis atoms: 

base(1,1) = 0.0d0
base(1,2) = 0.0d0
base(1,3) = 0.0d0
base(2,1) = 0.5d0*lpar
base(2,2) = 0.5d0*lpar
base(2,3) = 0.5d0*lpar

! Types of the basis atoms (change here in the case of more types in the basis):

itp(:) = 1

! Loop over the atoms:

iatom = 1

do ix = 0, nx - 1
  do iy = 0, ny - 1
    do iz = 0, nz - 1

! Loop over the basis vectors:

do ib = 1, 2
 rx(iatom) = ix*lpar + base(ib,1)
  ry(iatom) = iy*lpar + base(ib,2)
 rz(iatom) = iz*lpar + base(ib,3)
 itype(iatom) = itp(ib)
iatom = iatom + 1
enddo

! End of the loop over the basis vectors

    enddo
  enddo
enddo

! Shifting the center of the supercell to the (0,0,0) position            

do jk = 1, nat
 rx(jk) = rx(jk) - 0.5d0*sdx
  ry(jk) = ry(jk) - 0.5d0*sdy
 rz(jk) = rz(jk) - 0.5d0*sdz
enddo

! End of the subroutine

end subroutine bccpos

! -----------------------------------------------------------------------------
! Generates a supercell with the atomic positions in the HCP lattice                
! -----------------------------------------------------------------------------

subroutine hcppos(nat,nx,ny,nz,itype,lpar,cova,rx,ry,rz)
implicit none

! ----------------- Variables and Parameters of the Subroutine ----------------
! External
integer,   intent(in) :: nx, ny, nz, nat
real*8,    intent(in) :: lpar, cova
real*8, intent(inout) :: rx(nat), ry(nat), rz(nat)
integer,intent(inout) :: itype(nat)
! Calculated in the subroutine 
real*8                :: cp, sdx, sdy, sdz
! Internals
real*8                ::  base(2,3), itp(2)
integer               :: iatom, ix, iy, iz, ib, jk
! ------------------------------------------------------------------------

sdx = nx*lpar
sdy = ny*lpar
sdz = nz*lpar*cova

! Basis atoms: 

base(1,1) = 0.5d0*lpar
base(1,2) = (dsqrt(3.d0)/6.d0)*lpar
base(1,3) = lpar*cova/4.d0
base(2,1) = 0.5d0*lpar
base(2,2) =-(dsqrt(3.d0)/6.d0)*lpar
base(1,3) = 3.d0*lpar*cova/4.d0

! Types of the basis atoms (change here in the case of more types in the basis):

itp(:) = 1

! Loop over the atoms:

iatom = 1

do ix = 0, nx - 1
  do iy = 0, ny - 1
    do iz = 0, nz - 1

! Loop over the basis vectors:

do ib = 1, 2
 rx(iatom) = ix*lpar + base(ib,1)
  ry(iatom) = iy*lpar + base(ib,2)
 rz(iatom) = iz*lpar + base(ib,3)
 itype(iatom) = itp(ib)
iatom = iatom + 1
enddo

! End of the loop over the basis vectors

    enddo
  enddo
enddo

! Shifting the center of the supercell to the (0,0,0) position            

do jk = 1, nat
 rx(jk) = rx(jk) - 0.5d0*sdx
  ry(jk) = ry(jk) - 0.5d0*sdy
 rz(jk) = rz(jk) - 0.5d0*sdz
enddo

! End of the subroutine

end subroutine hcppos

! -----------------------------------------------------------------------------
! Generates a supercell with the atomic positions in a free cell mode, that is,
! In a general unit cell readen from the ATOMS.DAT input file
! -----------------------------------------------------------------------------

subroutine freepos(nat,nx,ny,nz,itype,rx,ry,rz,lpar,cova)
implicit none

! ----------------- Variables and Parameters of the Subroutine ----------------
! External
integer,   intent(in) :: nx, ny, nz, nat
real*8,    intent(in) :: lpar, cova
real*8, intent(inout) :: rx(nat), ry(nat), rz(nat)
integer,intent(inout) :: itype(nat)
! Internals
real*8, allocatable  ::  base(:,:), itp(:)
real*8               :: sdx, sdy, sdz
real*8               :: pv1(3), pv2(3), pv3(3), lp(3)
integer              :: nauc, iatom, i, ix, iy, iz, ib, jk
! -----------------------------------------------------------------------------

! UNIT CELL (Read from ATOMS.DAT input file):

open(8, file ='atoms.dat', status='unknown')
read(8,*) nauc

! Reading the primitive vectors:

 read(8,*)  pv1(1), pv1(2), pv1(3)
 read(8,*)  pv2(1), pv2(2), pv2(3)
 read(8,*)  pv3(1), pv3(2), pv3(3)

! Reading the types and the coordinates of the basis vectors:

allocate(base(nauc,3),itp(nauc))

do i = 1, nauc
read(8,*) itp(i), base(i,1), base(i,2), base(i,3)
enddo

close(8)

! Scaling the lattice Primitive Vectors:

pv1(:) = pv1(:)*lpar
pv2(:) = pv2(:)*lpar
pv3(:) = pv3(:)*lpar*cova

! Lattice parameters:

lp(1) = dsqrt(pv1(1)*pv1(1) + pv2(1)*pv2(1) + pv3(1)*pv3(1))
 lp(2) = dsqrt(pv1(2)*pv1(2) + pv2(2)*pv2(2) + pv3(2)*pv3(2))
  lp(3) = dsqrt(pv1(3)*pv1(3) + pv2(3)*pv2(3) + pv3(3)*pv3(3))
  sdx = nx*lp(1)
 sdy = ny*lp(2)
sdz = nz*lp(3)

! Loop over the atoms:

iatom = 1

do ix = 0, nx - 1
  do iy = 0, ny - 1
    do iz = 0, nz - 1

! Loop over the basis vectors:

do ib = 1, nauc

! Alpha sub-lattice:

 rx(iatom) = base(ib,1) + ix*pv1(1) + iy*pv2(1) + iz*pv3(1)
  ry(iatom) = base(ib,2) + ix*pv1(2) + iy*pv2(2) + iz*pv3(2)
 rz(iatom) = base(ib,3) + ix*pv1(3) + iy*pv2(3) + iz*pv3(3)

! Setting the type:

itype(iatom) = itp(ib)

iatom = iatom + 1

! End of the loop over the basis vectors

enddo

    enddo
  enddo
enddo

! Shifting the center of the supercell to the (0,0,0) position            

!do jk = 1, nat
! rx(jk) = rx(jk) - 0.5d0*sdx
!  ry(jk) = ry(jk) - 0.5d0*sdy
! rz(jk) = rz(jk) - 0.5d0*sdz
!enddo

! End of the subroutine

end subroutine freepos

! -----------------------------------------------------------------------------
! Generates, for all atoms of the system, their neighbor lists LISNGH
! -----------------------------------------------------------------------------

subroutine neighbors(ibrav,nat,nnn,nx,ny,nz,lpar,pv1,pv2,pv3,rx,ry,rz,rcut,nngh,lisngh)

! ----------------- Variables and Parameters of the Subroutine ----------------
! External
integer,   intent(in) :: ibrav, nat, nnn, nx, ny, nz
real*8,    intent(in) :: rcut, lpar, pv1(3), pv2(3), pv3(3), rx(nat), ry(nat), rz(nat)
integer,intent(inout) :: nngh(nat), lisngh(nat,nnn)
! Internals
integer, parameter    :: ncell = 1
real*8                :: rij, xi, yi, zi, xj, yj, zj, xij, yij, zij, sdx, sdy, sdz
integer               :: icx, icy, icz, i, j, k, icount
! -----------------------------------------------------------------------------

select case(ibrav)

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
case(1,2,3,4,5)  ! Supercells with Orthogonal Primitive Vectors
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

sdx = nx*lpar
sdy = ny*lpar
sdz = nz*lpar

! Loop over all (I,J) pairs:

do i = 1, nat

icount = 0

! Coordinates of i atom:

xi = rx(i)
yi = ry(i)
zi = rz(i)

  do j = 1, nat

   if(j.eq.i) cycle

! Coordinates of j atom:

  xj = rx(j)
  yj = ry(j)
  zj = rz(j)

! Distance rij = |ri - rj|:

  xij = (xj - xi) - dnint((xj - xi)/sdx)*sdx
  yij = (yj - yi) - dnint((yj - yi)/sdy)*sdy
  zij = (zj - zi) - dnint((zj - zi)/sdz)*sdz
  rij = dsqrt(xij**2 + yij**2 + zij**2)

! Verifying if J is neighbor of the I:

    if(rij.le.rcut) then
     icount = icount + 1
     lisngh(i,icount) = j
    endif

! End of the loops or all (I,J) pairs:

  enddo  ! End of the J Loop

nngh(i) = icount

if(icount.gt.nnn) then
print*
write(*,'(a)') 'ERROR: the estimated value of NNN will not take all neighbors of'
write(*,'(a)') '       the I atom. IMPROVED IT'
print*
 stop
endif

enddo  ! End of the I Loop

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
case(6)  ! Supercells with non-Orthogonal Primitive Vectors
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

! Loop over all (I,J) pairs:

do i = 1, nat

icount = 0

! Coordinates of i atom:

xi = rx(i)
yi = ry(i)
zi = rz(i)

! Loop over the copies of the supercell:

do icx = -ncell, ncell
  do icy = -ncell, ncell
    do icz = -ncell, ncell

! Loop over the J atoms:

  do j = 1, nat

   if(j.eq.i) cycle

! Coordinates of j atom:

  xj = rx(j) + icx*nx*pv1(1) + icy*ny*pv2(1) + icz*nz*pv3(1)
  yj = ry(j) + icx*nx*pv1(2) + icy*ny*pv2(2) + icz*nz*pv3(2)
  zj = rz(j) + icx*nx*pv1(3) + icy*ny*pv2(3) + icz*nz*pv3(3)

! Distance rij = |ri - rj|:

  xij = (xj - xi) 
  yij = (yj - yi) 
  zij = (zj - zi) 
  rij = dsqrt(xij**2 + yij**2 + zij**2)

! Verifying if J is neighbor of the I:

    if(rij.le.rcut) then
     icount = icount + 1
     lisngh(i,icount) = j
    endif

! End of the J loop:

  enddo

! Assigning the number of neighbors of the I atom:

nngh(i) = icount

if(icount.gt.nnn) then
print*
write(*,'(a)') 'ERROR: the estimated value of NNN will not take all neighbors of'
write(*,'(a)') '       the I atom. IMPROVED IT'
print*
stop
endif

! End of the loops over the cells:

    enddo
  enddo
enddo

! End of the I loop:

enddo

end select

end subroutine neighbors

! -----------------------------------------------------------------------------
! End of the SUPERCELL module
! -----------------------------------------------------------------------------

end module supercell  

! -----------------------------------------------------------------------------
