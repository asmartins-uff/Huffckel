! -----------------------------------------------------------------------------
! MODULE POSTPROC: performs the post-processing of the calculated electronic 
! States. 
!
! (1) PP_REAL: Calculates the Mullinken Charges of the System.
!
! (2) PP_RECP: Calculates the total density of states applying gaussian
!              Broadening to the electronic states.
!
! LAST MODIFICATION: 28/03/2014
! -----------------------------------------------------------------------------

module postproc 
implicit none

public:: pp_real, pp_recp                        

CONTAINS

!-----------------------------------------------------------------------
! SUBROUTINE PP_REAL: Computes the Mulliken Charges of the atoms and the
!                     Total and Partial density of states for supercells
!                     (Real Space) Reference: Andrew Leach, pp 77 - 80  
!-----------------------------------------------------------------------

subroutine pp_real(nat,tnao,numtp,ibrav,neao,nelect,norbv,itype,dbov,W,ham,ovl,imc,ids) 
implicit none

! ----------------- Variables and Parameters of the Subroutine ----------------
! Externals:
integer, parameter  :: lmax = 2
integer, intent(in) :: nat, tnao, numtp, ibrav, norbv(numtp), itype(nat), neao(numtp,3)
integer, intent(in) :: nelect, imc, ids, dbov(numtp)
real*8 , intent(in) :: W(tnao), ham(tnao,tnao), ovl(tnao,tnao) 
! Internals:
real*8, parameter   :: sig = 0.08d0
real*8              :: zk(numtp), sum, occ(tnao), sip(numtp)
real*8              :: erange, emax, emin, ef, df1, df2, rdf 
real*8              :: ene, deltae, broad, xx, pi, wao, yy
real*8, allocatable :: MC(:), MullQ(:), MCat(:), qmn(:,:), dos(:), pdos(:,:), atpdos(:,:)
integer             :: i, j, k, m, oi, ii, jj, id, li, ie, noccs, ndven, nauc, nuc
! -----------------------------------------------------------------------------

! Output files:

open(11,file='dos.dat',status='unknown')
open(13,file='MullQ.dat',status='unknown')
open(14,file='pdos.dat',status='unknown')

! -----------------------------------------------------------------------------
! Defining the number of atoms in the minimum unit cell, needed to normalize the
! Partial Density of States:

select case(ibrav)

case(1,4)   ! FCC and BCC lattices

nauc = 1

case(2,3,5) ! ZINCBLEND, SILICON and HCP lattices

nauc = 2

case(6)     ! Arbitrary Unit Cell

open(16,file='atoms.dat',status = 'unknown')
 read(16,*) nauc
close(16)

case default
write(*,'(a)') 'Subroutine POSTPROC: IBRAV value not recognized'

end select
! -----------------------------------------------------------------------------

! Defining PI:

pi = 3.141592654d0

! Defining the number of occupied states NOCCS:

if(mod(nelect,2).eq.0) then
 noccs = nelect/2             ! For EVEN value of NELECT
  else
 noccs = (nelect + 1)/2       ! For ODD value of NELECT
endif

! Defining the occupations:

occ(:) = 0.0d0
 do i = 1, noccs
  occ(i) = 2.0d0
 enddo

if(mod(nelect,2).eq.1) then
 occ(noccs) = 1.0d0         ! For odd nelect, the last band has occupation ONE
endif

! Setting ZK (ion charges) to the number of valence electrons for each specie:

zk(:) = 0.0d0
do i = 1, numtp
  do j = 1, norbv(i)
    zk(i) = zk(i) + 1.d0*neao(i,j)
  enddo
enddo

! Determining the energy range: if the difference between the fermi energy with
! The maximum value is greater that 1.5 times the difference of the fermi level
! With the minimum, EMAX is setted to the later difference

  emin = W(1) - 0.4d0
  emax = W(tnao) + 0.4d0
    ef = W(noccs)
   df1 = ef - emin
   df2 = emax - ef
   rdf = df2/df1
   if(rdf.gt.1.5d0) emax = ef + 1.5d0*df1

erange = emax - emin
ndven  = 2*dnint(erange/sig)
deltae = erange/dfloat(ndven)

! -----------------------------------------------------------------------------
!                              MULLIKEN CHARGES 
! -----------------------------------------------------------------------------

if(imc.eq.1) then

allocate(MC(tnao), MCat(nat), MullQ(nat), qmn(tnao,tnao))
MC(:)   = 0.0d0
qmn(:,:) = 0.0d0
MCat(:) = 0.0d0

! Calculating MC, the mulliken charge associated to each orbital of the system

do i = 1, tnao
  do m = 1, tnao 
   if(W(m).gt.emax) cycle
    do j = 1, tnao
      qmn(i,m) =     qmn(i,m) + ham(i,m)*ovl(i,j)*ham(j,m)
         MC(i) = MC(i) + occ(m)*ham(i,m)*ovl(i,j)*ham(j,m)
    enddo
  enddo
enddo

! Checkin the sum of the mulliken charges:

sum = 0.0d0
do k = 1, tnao
sum = sum + MC(k)
enddo

write(*,'(a,f16.2)') 'Sum of The Mulliken Charges =', sum

! Defining the Mulliken charges for each atom:

id = 0
do i = 1, nat
  do j = 1, dbov(itype(i))
    id = id + 1
    MCat(i) = MCat(i) + MC(id)
  enddo
enddo

! Taking the difference of Z(i) - MCat(i):

do j = 1, nat
 MullQ(j) = zk(itype(j)) - MCat(j)
enddo

! Writing the Mulliken Charges:

do i = 1, nat      
write(13,'(i6,2f12.4)') i, MCat(i), MullQ(i)
enddo

deallocate(MC,MullQ,MCat)

endif

! ------------------------------------------------------------------------
!                  TOTAL AND PROJECTED DENSITY OF STATES:
! ------------------------------------------------------------------------

if(ids.eq.1) then

! The ATPDOS here will be related to each specie and their corresponding
! Orbitlas:

id = 0
do k = 1, numtp
 id = id + dbov(k)
enddo

allocate(dos(ndven),pdos(tnao,ndven),atpdos(id,ndven))
     dos(:) = 0.0d0
  pdos(:,:) = 0.0d0
atpdos(:,:) = 0.0d0

do i = 1, tnao      ! Loop over all orbitals
    do ie = 0, ndven - 1

! DOS:

      ene = emin + ie*deltae
       xx = ene - W(i)
       broad =  (1.d0/(sig*dsqrt(2.d0*pi)))*dexp(-xx*xx/(2.d0*sig*sig))
      dos(ie+1) = dos(ie+1) + 2.0d0*deltae*broad

! Atom weighted PDOS:

     do m = 1, tnao      ! Loop over the bands
        wao = qmn(i,m)
         yy = ene - W(m)
         broad =  (1.d0/(sig*dsqrt(2.d0*pi)))*dexp(-yy*yy/(2.d0*sig*sig))
       pdos(i,ie+1) = pdos(i,ie+1) + 2.d0*wao*deltae*broad
     enddo

    enddo    ! End of the Loop over the energy divisions
enddo        ! End of the I loop

! Defining ATPDOS:

sip(1) = 0

do i = 2, numtp
sip(i) = sip(i-1) + dbov(i-1)
enddo

jj = 0
do i = 1, nat
 do oi = 1, dbov(itype(i))

  jj = jj + 1
  ii = sip(itype(i)) + oi
  atpdos(ii,:) = atpdos(ii,:) + pdos(jj,:)

 enddo
enddo

! Normalizing DOS and PDOS, in order to give the number of electrons in the
! Unit Cell:

 nuc = nat/nauc
     dos(:) = dos(:)/dfloat(nuc)
atpdos(:,:) = atpdos(:,:)/dfloat(nuc)

! Writing DOS and PDOS:

do i = 0, ndven - 1
  ene = emin + i*deltae
  write(11,'(2f12.4)')  ene, dos(i+1)
  write(14,'(19f12.4)') ene, (atpdos(j,i+1), j = 1, id)
enddo

deallocate(dos,pdos,qmn)

endif

! ------------------------------------------------------------------------

! End of the subroutine

end subroutine pp_real 

!-----------------------------------------------------------------------
! SUBROUTINE PP_RECP: Computes the Density of States and Mulliken Charges
!                     For K-dependent states.
!-----------------------------------------------------------------------

subroutine pp_recp(numtp,nat,tnao,nelect,nkpt,itype,dbov,eband,hamk,ovlk,flp,imc,ids,icp1,icp2) 
implicit none

! ----------------- Variables and Parameters of the Subroutine ----------------
! Externals:
integer, intent(in)     :: numtp, nat, tnao, nelect, nkpt, flp, imc, ids, icp1, icp2
integer, intent(in)     :: itype(nat), dbov(numtp) 
real*8 , intent(in)     :: eband(nkpt,tnao) 
complex*16              :: hamk(nkpt,tnao,tnao), ovlk(nkpt,tnao,tnao)
! Internals:
real*8,  parameter      :: sig = 0.08d0
integer                 :: i, j, k, io, jo, ko, ik, id, jd, kd, m, md, maxlv, sip(nat)
integer                 :: kc, mc, ok, mo, ie, ndven
real*8                  :: occ(tnao), erange, emax, emin, ene, deltae, broad
real*8                  :: e1, e2, xx, yy, pi, sum, sumacc(nat), chrg, wao, wop
real*8, allocatable     :: pdos(:,:), atpdos(:,:), dos(:), MullQ(:), qmn(:,:,:)
real*8, allocatable     :: ovp(:,:,:), coop(:,:,:)
complex*16, allocatable :: uk(:,:)
! ------------------------------------------------------------------------

! Output file

open(11,file='dos.dat',status='unknown')
open(13,file='pdos.dat',status='unknown')
open(15,file='coop.dat',status='unknown')

! Defining PI:

pi = 3.141592654d0

! Occupations:

    occ(:) = 0.0d0
occ(1:flp) = 2.d0

if(mod(nelect,2).eq.1) then
 occ(flp) = 1.0d0            ! For ODD nelect, the last band has occupation ONE
endif

! Maximum level:

maxlv = flp 

! Determining the energy range:

  emin = minval(eband(:,1)) - 0.2d0
  emax = maxval(eband(:,tnao)) + 0.2d0
erange = emax - emin
ndven  = 2*dnint(erange/sig)   
deltae = erange/dfloat(ndven)

allocate(pdos(tnao,ndven),dos(ndven))

    pdos(:,:) = 0.0d0

! Pointer SIP:

sip(1) = 0

do i = 2, nat
 sip(i) = sip(i-1) + dbov(itype(i-1))
enddo

! ------------------------------------------------------------------------
!                 NORMALIZATION OF THE EIGENSTATES
! ------------------------------------------------------------------------

allocate(uk(tnao,tnao))

do k = 1, nkpt  ! Loop over the k-points

uk(:,:) = 0.0d0

do m = 1, tnao  ! Loop over the bands

uk(:,m) = hamk(k,:,m)

sum = 0.0d0

  do i = 1, tnao
    do j = 1, tnao
     sum = sum + conjg(uk(i,m))*ovlk(k,i,j)*uk(j,m)
    enddo
  enddo

hamk(k,:,m) = hamk(k,:,m)/dsqrt(sum)

enddo
enddo

! ------------------------------------------------------------------------
!                           MULLIKEN CHARGES:
! ------------------------------------------------------------------------

if(imc.eq.1) then

allocate(MullQ(tnao),qmn(nkpt,tnao,tnao))

qmn(:,:,:) = 0.0d0     ! QMN(k,orb,band)
MullQ(:) = 0.0d0

do i = 1, tnao    ! I orbitals

do k = 1, nkpt    ! Loop over the K-Points

 do m = 1, tnao   ! Bands
  do j = 1, tnao  ! J orbitals 

! Mulliken charge for I orbital and charge matrix:

   MullQ(i) = MullQ(i) + occ(m)*conjg(hamk(k,i,m))*ovlk(k,i,j)*hamk(k,j,m)
   qmn(k,m,i) = qmn(k,m,i) + conjg(hamk(k,m,i))*ovlk(k,m,j)*hamk(k,j,i)

  enddo
 enddo

enddo
enddo

! Normalizing and printing the Mulliken Charges:

MullQ(:) = MullQ(:)/dfloat(nkpt)

       sum = 0.0d0
sumacc(:) = 0.0d0

  do i = 1, nat 

   id = sip(i) 

   do io = 1, dbov(itype(i))

   id = id + 1
   sumacc(i) = sumacc(i) + MullQ(id)
   sum = sum + MullQ(id)

   enddo
  enddo

write(*,'(a,a,a,a,a,a,a)') ' | ', ' Atom ', ' | ', ' Type ', ' | ', ' Mulliken Charge ', ' | '
print*
do i = 1, nat  
write(*,'(i7,i9,f16.4)') i, itype(i), sumacc(i)
enddo
write(*,'(a)')  '  --------------------------------------'
write(*,'(a,f10.4)') '          SUM  = ', sum
write(*,'(a)')  '  --------------------------------------'
print*

endif

! ------------------------------------------------------------------------
!                  TOTAL AND PROJECTED DENSITY OF STATES:
! ------------------------------------------------------------------------

if(ids.eq.1) then

! **********************************
! (a) K-WEIGHTED OVERLAP POPULATION:
! **********************************

allocate(ovp(tnao,tnao,tnao),coop(tnao,tnao,ndven))  ! OVP(mu,nu,band)
 ovp(:,:,:) = 0.0d0
coop(:,:,:) = 0.0d0
id = 0

do i = 1, nat               ! Atoms I    
do io = 1, dbov(itype(i))   ! I's orbitals

id = id + 1

  do k = 1, nkpt    ! Loop over the K-Points
    do m = 1, tnao  ! Loop over the bands

      do j = i + 1, nat          ! Atom J

      jd = sip(j)

      do jo = 1, dbov(itype(j))  ! J's orbitals

! Overlap Population:

       jd = jd + 1
!      ovp(id,jd,m) = ovp(id,jd,m) + 2.d0*occ(m)*conjg(hamk(k,id,m))*ovlk(k,id,jd)*hamk(k,jd,m)
       ovp(id,jd,m) = ovp(id,jd,m) + 4.d0*conjg(hamk(k,id,m))*ovlk(k,id,jd)*hamk(k,jd,m)
       ovp(jd,id,m) = ovp(id,jd,m)

      enddo
      enddo

    enddo            ! End of the bands loop
  enddo              ! End of the K-point Loop

enddo
enddo

! ********************
! (b) DOS CALCULATION:
! ********************

sum = 0.0d0

do i = 1, tnao      ! Loop over all orbitals
  do k = 1, nkpt    ! Loop over the K-Points
    do ie = 0, ndven - 1

! DOS:

      ene = emin + ie*deltae
       xx = ene - eband(k,i)
       broad =  (1.d0/(sig*dsqrt(2.d0*pi)))*dexp(-xx*xx/(2.d0*sig*sig))
      dos(ie+1) = dos(ie+1) + 2.0d0*deltae*broad

! Atom weighted PDOS:

       do m = 1, tnao      ! Loop over the bands

        wao = qmn(k,i,m)
         yy = ene - eband(k,m)
         broad =  (1.d0/(sig*dsqrt(2.d0*pi)))*dexp(-yy*yy/(2.d0*sig*sig))
       pdos(i,ie+1) = pdos(i,ie+1) + 2.d0*wao*deltae*broad

       do j = 1, tnao
        wop = ovp(i,j,m)
         yy = ene - eband(k,m)
         broad =  (1.d0/(sig*dsqrt(2.d0*pi)))*dexp(-yy*yy/(2.d0*sig*sig))
        coop(i,j,ie+1) = coop(i,j,ie+1) + wop*deltae*broad
       enddo

       enddo

    enddo    ! End of the Loop over the energy divisions
  enddo      ! End of the K loop
enddo        ! End of the I loop

! Normalizing DOS, PDOS and COOP:

     dos(:) = dos(:)/dfloat(nkpt)
  pdos(:,:) = pdos(:,:)/dfloat(nkpt)
coop(:,:,:) = coop(:,:,:)/dfloat(nkpt)

! Summing the PDOS for each atom:

allocate(atpdos(nat,ndven))

do   i = 1, nat
id = sip(i)
 do io = 1, dbov(itype(i))
  id = id + 1
  atpdos(i,:) = atpdos(i,:) + pdos(id,:)
 enddo
enddo

! Writing to the DOS.DAT file:

id = icp1
jd = icp2

do i = 0, ndven - 1
  ene = emin + i*deltae
  write(11,'(2f12.4)')  ene, dos(i+1)
! write(13,'(f12.4,18f10.4)')  ene, (pdos(j,i+1), j = 1, tnao)
  write(13,'(f12.4,18f10.4)')  ene, (atpdos(j,i+1), j = 1, nat) 
  write(15,'(f12.4,4f10.4)')  ene, coop(id,jd,i+1), coop(id,jd+1,i+1), coop(id,jd+2,i+1), coop(id,jd+3,i+1)
enddo

endif

! End of the subroutine 

end subroutine pp_recp

! -----------------------------------------------------------------------------
end module postproc
! -----------------------------------------------------------------------------
