  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                        
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !
  !--------------------------------------------------------------------------
  SUBROUTINE dynwan2bloch ( nbnd, nrr, irvec, ndegen, xq, cuf, eig, vph)
  !--------------------------------------------------------------------------
  !
  !
  !  WARNING: this SUBROUTINE is identical to hamwan2bloch.f90, except
  !           that here rdw is a real array, not a complex one. This is
  !           required to obtain proper phonon dispersion interpolation
  !           and corresponds to the reality of the interatomic force
  !           constants
  !
  ! -------------------------------------------------------------------------
  !
  !  From the Hamiltonian in Wannier representation, find the corresponding
  !  Hamiltonian in Bloch representation for a given k point
  !
  !  input  : number of bands nbnd
  !           number of WS vectors, coordinates and degeneracy
  !           Hamiltonian in Wannier representation chw(nbnd, nbnd, nrr)
  !           qpoint coordinate xq(3)
  !
  !  output : rotation matrix cuf (nbnd, nbnd)
  !           interpolated hamiltonian eigenvalues eig(nbnd)
  !
  !
  !--------------------------------------------------------------------------
  !
#include "f_defs.h"
  USE kinds,     ONLY : DP
  USE pwcom,     ONLY : at, bg, omega
  USE phcom,     ONLY : nq1, nq2, nq3
  USE cell_base, ONLY : alat
  USE ions_base, ONLY : amass, tau, nat, ityp
  USE elph2,     ONLY : epsi, zstar, ifc
  USE epwcom,    ONLY : lpolar, vg_ph
  USE constants_epw, ONLY : pi, twopi, ci, czero
  implicit none
  !
  !  input variables
  !
  integer :: nbnd, nrr, irvec (3, nrr), ndegen (nrr)
  ! number of bands (possibly of the optimal subspace)
  ! kpoint number for the interpolation
  ! record length and unit for direct write of rotation matrix
  ! number of WS points, crystal coordinates, degeneracy
  !
  ! Hamiltonian in wannier basis
  !
  real(kind=DP) :: xq (3)
  ! kpoint coordinates for the interpolation
  !
  ! output variables
  !
  real(kind=DP) :: eig (nbnd)
  ! interpolated hamiltonian eigenvalues for this kpoint
  complex(kind=DP) :: cuf(nbnd, nbnd)
  ! Rotation matrix, fine mesh
  !
  ! variables for lapack ZHPEVX
  !
  integer :: neig, info, ifail( nbnd ), iwork( 5*nbnd )
  real(kind=DP) :: w( nbnd ), rwork( 7*nbnd )
  complex(kind=DP) :: champ( nbnd*(nbnd+1)/2 ), &
    cwork( 2*nbnd ), cz( nbnd, nbnd)
  !
  ! work variables
  !
  complex(kind=DP) :: chf(nbnd, nbnd)
  ! Hamiltonian in Bloch basis, fine mesh
  integer :: ibnd, jbnd, ir, na, nb
  real(kind=DP) :: rdotk, massfac
  complex(kind=DP) :: cfac
  !
  !
  ! THL
  REAL(KIND=DP), PARAMETER :: eps = 1.0d-6
  INTEGER, PARAMETER       :: nrwsx = 200
  ! the maximun number of nearest neighbor
  INTEGER                  :: nrws             
  ! number of nearest neighbor
  REAL(KIND=DP)            :: atws(3,3), rws(0:3,nrwsx), r(3), r_ws(3), weight
  ! lattice vector for WS initialization
  ! nearest neighbor list, rws(0,*) = norm^2
  COMPLEX(KIND=DP)         :: ddyn_s(3,nbnd,nbnd), ddyn_l(3,nbnd,nbnd), ddyn(3,nbnd,nbnd)
  ! short-range part of derivative of DM
  ! long-range (non-analytical) part of derivative of DM
  REAL(KIND=DP)            :: vph(3,nbnd)
  INTEGER                  :: i, j, k, n1, n2, n3, m1, m2, m3, mu, nu
  REAL(KIND=DP), EXTERNAL  :: wsweight
  !
  !
  CALL start_clock ( 'DynW2B' )
  !
  !
  !------------------------------------------------------------------------------------------------
  ! THL: Compute chf (dynamical matrix) through ifc. The obtained phonon frequency now is 
  !      identical with that from matdyn.x
  !------------------------------------------------------------------------------------------------
  !
  !
  atws(:,1) = at(:,1)*DBLE(nq1)
  atws(:,2) = at(:,2)*DBLE(nq2)
  atws(:,3) = at(:,3)*DBLE(nq3)
  CALL wsinit(rws,nrwsx,nrws,atws)
  !
  ! Fourier interpolation 
  chf = czero
  ddyn_s = czero
  ddyn_l = czero
  !
  DO na = 1,nat
     DO nb = 1,nat
        !
        DO n1 = -2*nq1,2*nq1
           DO n2 = -2*nq2,2*nq2
              DO n3 = -2*nq3,2*nq3
                 !
                 DO k = 1,3
                    r(k) = n1*at(k,1)+n2*at(k,2)+n3*at(k,3)
                    r_ws(k) = r(k) + tau(k,na) - tau(k,nb)
                 ENDDO
                 weight = wsweight(r_ws,rws,nrws)
                 !
                 IF (weight .GT. 0.0d0) THEN
                    !
                    m1 = MOD(n1+1,nq1)
                    IF(m1.LE.0) m1=m1+nq1
                    m2 = MOD(n2+1,nq2)
                    IF(m2.LE.0) m2=m2+nq2
                    m3 = MOD(n3+1,nq3)
                    IF(m3.LE.0) m3=m3+nq3
                    !
                    DO i = 1,3
                       DO j = 1,3
                          !
                          mu = (na-1)*3+i
                          nu = (nb-1)*3+j
                          !
                          rdotk = twopi * (  xq(1)*DBLE(n1) + xq(2)*DBLE(n2) + xq(3)*DBLE(n3) )
                          cfac = EXP(-ci*rdotk)*weight
                          !
                          chf(mu,nu) = chf(mu,nu) + cfac * REAL(ifc(m1,m2,m3,mu,nu))
                          !
                          IF (vg_ph .EQ. 'matrix') ddyn_s(:,mu,nu) = ddyn_s(:,mu,nu) - ci * alat * r(:) * cfac * REAL(ifc(m1,m2,m3,mu,nu))
                          !
                       ENDDO
                    ENDDO
                    !
                 ENDIF
                 !
              ENDDO
           ENDDO
        ENDDO  
        !
     ENDDO 
  ENDDO 
  !
  !
  ! bring xq in cart. coordinates (needed for rgd_blk call)
  CALL cryst_to_cart (1, xq, bg, 1)
  !
  !  add the long-range term to D(q)
  IF (lpolar) THEN
     CALL rgd_blk (nq1,nq2,nq3,nat,chf,1,1,ddyn_l,xq, &  !xq has to be in 2pi/a
                  tau,epsi,zstar,bg,omega,+1.d0)
   !WRITE (6,'(a)') "Done rigid"
  ENDIF
  !
  !  divide by the square root of masses 
  !
  DO na = 1, nat
   DO nb = 1, nat
      massfac = 1.d0 / sqrt ( amass(ityp(na)) * amass(ityp(nb)) )
      !
      chf(3*(na-1)+1:3*na, 3*(nb-1)+1:3*nb) = &
         chf(3*(na-1)+1:3*na, 3*(nb-1)+1:3*nb) * massfac
      ! 
   ENDDO
  ENDDO
  !
  ! bring xq back to crystal coordinates
  CALL cryst_to_cart (1, xq, at, -1)
  !
  !---------------------------------------------------------------------
  !  STEP 4: diagonalize smooth Hamiltonian on k points of the fine grid
  !---------------------------------------------------------------------
  !
  ! champ: complex hamiltonian packed (upper triangular part for zhpevx)
  ! after hermitian-ization
  !
  DO jbnd = 1, nbnd
   DO ibnd = 1, jbnd
      champ (ibnd + (jbnd - 1) * jbnd/2 ) = &
      ( chf ( ibnd, jbnd) + conjg ( chf ( jbnd, ibnd) ) ) / 2.d0
   ENDDO
  ENDDO
  !
  CALL zhpevx ('V', 'A', 'U', nbnd, champ , 0.0, 0.0, &
               0, 0,-1.0, neig, w, cz, nbnd, cwork, &
               rwork, iwork, ifail, info)
  !
  ! rotation matrix and Ham eigenvalues
  ! [in Ry, mind when comparing with wannier code]
  !
  cuf = cz
  eig = w
  !
  !
  ! THL : Compute the phonon group velocity from dynamical matrix
  !
  IF (vg_ph .EQ. 'matrix') THEN
     !
     ddyn = ddyn_s + ddyn_l
     !  
     DO na = 1, nat
        DO nb = 1, nat
           massfac = 1.d0 / SQRT (amass(ityp(na))*amass(ityp(nb)))
           ddyn(:,3*(na-1)+1:3*na,3*(nb-1)+1:3*nb) = ddyn(:,3*(na-1)+1:3*na,3*(nb-1)+1:3*nb) * massfac
        ENDDO
     ENDDO
     !
     DO ibnd = 1, nbnd
        DO i = 1, 3
           vph(i,ibnd) = REAL( DOT_PRODUCT( cuf(:,ibnd) , MATMUL(ddyn(i,:,:),cuf(:,ibnd)) ) ) / (2.0d0*SQRT(ABS(eig(ibnd))))
        ENDDO
     ENDDO
     !
  ENDIF
  !
  !
  !
  CALL stop_clock ( 'DynW2B' )
  !   
  END SUBROUTINE dynwan2bloch
  !
  !
  !--------------------------------------------------------------------------
  SUBROUTINE dynwan2bloch_s ( nbnd, nrr, irvec, ndegen, xq, cuf, eig, vph, ntemp, ndope)
  !--------------------------------------------------------------------------
  !
  ! note: this copied version of dynwan2bloch takes into account the
  !       screening effect
  !
  ! ERROR: only cuf, eig are updated, ddyn_l, vph not updated
  !
  !--------------------------------------------------------------------------
  !
#include "f_defs.h"
  USE kinds,     ONLY : DP
  USE pwcom,     ONLY : at, bg, omega
  USE phcom,     ONLY : nq1, nq2, nq3
  USE cell_base, ONLY : alat
  USE ions_base, ONLY : amass, tau, nat, ityp
  USE elph2,     ONLY : epsi, zstar, ifc
  USE epwcom,    ONLY : lpolar, vg_ph
  USE constants_epw, ONLY : pi, twopi, ci, czero
  implicit none
  !
  !  input variables
  !
  integer :: nbnd, nrr, irvec (3, nrr), ndegen (nrr), ntemp, ndope
  ! number of bands (possibly of the optimal subspace)
  ! kpoint number for the interpolation
  ! record length and unit for direct write of rotation matrix
  ! number of WS points, crystal coordinates, degeneracy
  !
  ! Hamiltonian in wannier basis
  !
  real(kind=DP) :: xq (3)
  ! kpoint coordinates for the interpolation
  !
  ! output variables
  !
  real(kind=DP) :: eig (nbnd, ntemp, ndope)
  ! interpolated hamiltonian eigenvalues for this kpoint
  complex(kind=DP) :: cuf(nbnd, nbnd, ntemp, ndope)
  ! Rotation matrix, fine mesh
  !
  ! variables for lapack ZHPEVX
  !
  integer :: neig, info, ifail( nbnd ), iwork( 5*nbnd )
  real(kind=DP) :: w( nbnd ), rwork( 7*nbnd )
  complex(kind=DP) :: champ( nbnd*(nbnd+1)/2, ntemp, ndope ), &
    cwork( 2*nbnd ), cz( nbnd, nbnd)
  !
  ! work variables
  !
  complex(kind=DP) :: chf(nbnd, nbnd, ntemp, ndope)
  ! Hamiltonian in Bloch basis, fine mesh
  integer :: ibnd, jbnd, ir, na, nb, itemp, idope
  real(kind=DP) :: rdotk, massfac
  complex(kind=DP) :: cfac
  !
  !
  ! THL
  REAL(KIND=DP), PARAMETER :: eps = 1.0d-6
  INTEGER, PARAMETER       :: nrwsx = 200
  ! the maximun number of nearest neighbor
  INTEGER                  :: nrws             
  ! number of nearest neighbor
  REAL(KIND=DP)            :: atws(3,3), rws(0:3,nrwsx), r(3), r_ws(3), weight
  ! lattice vector for WS initialization
  ! nearest neighbor list, rws(0,*) = norm^2
  COMPLEX(KIND=DP)         :: ddyn_s(3,nbnd,nbnd), ddyn_l(3,nbnd,nbnd), ddyn(3,nbnd,nbnd)
  ! short-range part of derivative of DM
  ! long-range (non-analytical) part of derivative of DM
  REAL(KIND=DP)            :: vph(3,nbnd,ntemp,ndope)
  INTEGER                  :: i, j, k, n1, n2, n3, m1, m2, m3, mu, nu
  REAL(KIND=DP), EXTERNAL  :: wsweight
  !
  !
  CALL start_clock ( 'DynW2B' )
  !
  !
  !------------------------------------------------------------------------------------------------
  ! THL: Compute chf (dynamical matrix) through ifc. The obtained phonon frequency now is 
  !      identical with that from matdyn.x
  !------------------------------------------------------------------------------------------------
  !
  !
  atws(:,1) = at(:,1)*DBLE(nq1)
  atws(:,2) = at(:,2)*DBLE(nq2)
  atws(:,3) = at(:,3)*DBLE(nq3)
  CALL wsinit(rws,nrwsx,nrws,atws)
  !
  ! Fourier interpolation 
  chf = czero
  ddyn_s = czero
  ddyn_l = czero
  !
  DO na = 1,nat
     DO nb = 1,nat
        !
        DO n1 = -2*nq1,2*nq1
           DO n2 = -2*nq2,2*nq2
              DO n3 = -2*nq3,2*nq3
                 !
                 DO k = 1,3
                    r(k) = n1*at(k,1)+n2*at(k,2)+n3*at(k,3)
                    r_ws(k) = r(k) + tau(k,na) - tau(k,nb)
                 ENDDO
                 weight = wsweight(r_ws,rws,nrws)
                 !
                 IF (weight .GT. 0.0d0) THEN
                    !
                    m1 = MOD(n1+1,nq1)
                    IF(m1.LE.0) m1=m1+nq1
                    m2 = MOD(n2+1,nq2)
                    IF(m2.LE.0) m2=m2+nq2
                    m3 = MOD(n3+1,nq3)
                    IF(m3.LE.0) m3=m3+nq3
                    !
                    DO i = 1,3
                       DO j = 1,3
                          !
                          mu = (na-1)*3+i
                          nu = (nb-1)*3+j
                          !
                          rdotk = twopi * (  xq(1)*DBLE(n1) + xq(2)*DBLE(n2) + xq(3)*DBLE(n3) )
                          cfac = EXP(-ci*rdotk)*weight
                          !
                          chf(mu,nu,:,:) = chf(mu,nu,:,:) + cfac * REAL(ifc(m1,m2,m3,mu,nu))
                          !
                          IF (vg_ph .EQ. 'matrix') ddyn_s(:,mu,nu) = ddyn_s(:,mu,nu) - ci * alat * r(:) * cfac * REAL(ifc(m1,m2,m3,mu,nu))
                          !
                       ENDDO
                    ENDDO
                    !
                 ENDIF
                 !
              ENDDO
           ENDDO
        ENDDO  
        !
     ENDDO 
  ENDDO 
  !
  !
  ! bring xq in cart. coordinates (needed for rgd_blk call)
  CALL cryst_to_cart (1, xq, bg, 1)
  !
  !  add the long-range term to D(q)
  IF (lpolar) THEN
     CALL rgd_blk (nq1,nq2,nq3,nat,chf,ntemp,ndope,ddyn_l,xq, &  !xq has to be in 2pi/a
                  tau,epsi,zstar,bg,omega,+1.d0)
   !WRITE (6,'(a)') "Done rigid"
  ENDIF
  !
  !  divide by the square root of masses 
  !
  DO na = 1, nat
   DO nb = 1, nat
      massfac = 1.d0 / sqrt ( amass(ityp(na)) * amass(ityp(nb)) )
      !
      chf(3*(na-1)+1:3*na, 3*(nb-1)+1:3*nb, :,:) = &
         chf(3*(na-1)+1:3*na, 3*(nb-1)+1:3*nb, :,:) * massfac
      ! 
   ENDDO
  ENDDO
  !
  ! bring xq back to crystal coordinates
  CALL cryst_to_cart (1, xq, at, -1)
  !
  !---------------------------------------------------------------------
  !  STEP 4: diagonalize smooth Hamiltonian on k points of the fine grid
  !---------------------------------------------------------------------
  !
  ! champ: complex hamiltonian packed (upper triangular part for zhpevx)
  ! after hermitian-ization
  !
  DO jbnd = 1, nbnd
   DO ibnd = 1, jbnd
      champ (ibnd + (jbnd - 1) * jbnd/2, :, : ) = &
      ( chf ( ibnd, jbnd, :, :) + conjg ( chf ( jbnd, ibnd, :, :) ) ) / 2.d0
   ENDDO
  ENDDO
  !
  !
  do itemp = 1, ntemp
     do idope = 1, ndope
        CALL zhpevx ('V', 'A', 'U', nbnd, champ(:,itemp,idope) , 0.0, 0.0, &
                     0, 0,-1.0, neig, w, cz, nbnd, cwork, &
                     rwork, iwork, ifail, info)
        !
        ! rotation matrix and Ham eigenvalues
        ! [in Ry, mind when comparing with wannier code]
        !
        cuf(:,:,itemp,idope) = cz
        eig(:,itemp,idope) = w
     enddo
  enddo
  !
  !
  ! THL : Compute the phonon group velocity from dynamical matrix
  !
  IF (vg_ph .EQ. 'matrix') THEN
     !
     ddyn = ddyn_s + ddyn_l
     !  
     DO na = 1, nat
        DO nb = 1, nat
           massfac = 1.d0 / SQRT (amass(ityp(na))*amass(ityp(nb)))
           ddyn(:,3*(na-1)+1:3*na,3*(nb-1)+1:3*nb) = ddyn(:,3*(na-1)+1:3*na,3*(nb-1)+1:3*nb) * massfac
        ENDDO
     ENDDO
     !
     DO ibnd = 1, nbnd
        DO i = 1, 3
           do itemp = 1, ntemp
              do idope = 1, ndope
                 vph(i,ibnd,itemp,idope) = REAL( DOT_PRODUCT( cuf(:,ibnd,itemp,idope) , MATMUL(ddyn(i,:,:),cuf(:,ibnd,itemp,idope)) ) ) / (2.0d0*SQRT(ABS(eig(ibnd,itemp,idope))))
              enddo
           enddo
        ENDDO
     ENDDO
     !
  ENDIF
  !
  CALL stop_clock ( 'DynW2B' )
  !   
  END SUBROUTINE dynwan2bloch_s
  !
  !
  !--------------------------------------------------------------------------
  SUBROUTINE ifc_load (nrr, irvec)
  !--------------------------------------------------------------------------
  !
#INCLUDE "f_defs.h"
  !
  USE kinds,         ONLY : DP
  USE pwcom,         ONLY : at, bg
  USE phcom,         ONLY : nq1, nq2, nq3
  USE ions_base,     ONLY : nat
  USE elph2,         ONLY : rdw, ifc
  USE epwcom,        ONLY : ifc_read
  USE io_global,     ONLY : stdout
  USE io_files,      ONLY : tmp_dir, prefix
  USE constants_epw, ONLY : ci
#ifdef __PARA
  USE io_global,     ONLY : ionode_id
  USE mp_global,     ONLY : my_pool_id
#ENDIF
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)      :: nrr, irvec(3,nrr)
  ! irvec was calculated in wigner_seitz2
  !
  REAL(KIND=DP), PARAMETER :: eps = 1.0d-6
  INTEGER, PARAMETER       :: nrwsx = 200
  INTEGER                  :: nrws             
  REAL(KIND=DP)            :: atws(3,3), rws(0:3,nrwsx), r(3), r_ws(3), weight
  INTEGER                  :: ir, n1, n2, n3, m1, m2, m3, i, j, k, na, nb, mu, nu
  REAL(KIND=DP), EXTERNAL  :: wsweight
  REAL(KIND=DP)            :: nq1_tmp, nq2_tmp, nq3_tmp, i_tmp, j_tmp, na_tmp, nb_tmp, n1_tmp, n2_tmp, n3_tmp
  CHARACTER (LEN=256)      :: ifc_fmt
  !
  !
  ! If epwread is true, the subroutine dynbloch2wan will not be called
  ! and thus ifc is empty.
  atws(:,1) = at(:,1)*DBLE(nq1)
  atws(:,2) = at(:,2)*DBLE(nq2)
  atws(:,3) = at(:,3)*DBLE(nq3)
  CALL wsinit(rws,nrwsx,nrws,atws)
  !
  ! copy rdw to ifc when epwread is true
  IF (ALLOCATED(ifc)) DEALLOCATE (ifc)
  ALLOCATE (ifc(nq1,nq2,nq3,3*nat,3*nat))
  ifc = (0.0d0,0.0d0)
  !
  IF (.NOT. ifc_read) THEN
     !
     DO ir = 1, nrr
        !
        DO n1 = -2*nq1,2*nq1
           DO n2 = -2*nq2,2*nq2
              DO n3 = -2*nq3,2*nq3
                 !
                 DO k = 1,3
                    r(k) = n1*at(k,1)+n2*at(k,2)+n3*at(k,3)
                    r_ws(k) = r(k)
                 ENDDO
                 weight = wsweight(r_ws,rws,nrws)
                 !
                 IF (weight .GT. 0.0d0) THEN
                    !
                    m1 = MOD(n1+1,nq1)
                    IF(m1.LE.0) m1=m1+nq1
                    m2 = MOD(n2+1,nq2)
                    IF(m2.LE.0) m2=m2+nq2
                    m3 = MOD(n3+1,nq3)
                    IF(m3.LE.0) m3=m3+nq3
                    !
                    IF ( ABS(DBLE(irvec(1,ir))-DBLE(n1)).LT.eps .AND. &
                         ABS(DBLE(irvec(2,ir))-DBLE(n2)).LT.eps .AND. &
                         ABS(DBLE(irvec(3,ir))-DBLE(n3)).LT.eps ) THEN
                       !
                       ifc(m1,m2,m3,:,:) = rdw (:,:,ir)
                       !
                    ENDIF
                    !
                 ENDIF
                 !
              ENDDO
           ENDDO
        ENDDO  
        !
     ENDDO
     !
  ELSE
     !
     ! THL: input the force constat on coarse mesh
     !      format is the same as QE's (remove the caption in ifc file)
     !
     WRITE (stdout,'(/5x,a/)') 'Reading external interatomic force constant'
     !
     ifc_fmt = TRIM(prefix)//'.ifc'
     OPEN (69517,FILE=ifc_fmt,STATUS='old')
     READ (69517,*) nq1_tmp, nq2_tmp, nq3_tmp
     IF (NINT(nq1_tmp) .NE. nq1 .OR. NINT(nq2_tmp) .NE. nq2 .OR. NINT(nq3_tmp) .NE. nq3) &
        CALL errore ('ifc_load','coarse q mesh of input IFC file is wrong)',1)
     !
     DO i = 1,3
       DO j = 1,3
         DO na = 1,nat
           DO nb = 1,nat
             !
             READ (69517,*) i_tmp, j_tmp, na_tmp, nb_tmp
             IF (NINT(i_tmp) .NE. i .OR. NINT(j_tmp) .NE. j .OR. NINT(na_tmp) .NE. na .OR. NINT(nb_tmp) .NE. nb) &
                CALL errore ('ifc_load','ordering of force constant of input IFC file is wrong)',1)
             !
             mu = (NINT(na_tmp)-1)*3+NINT(i_tmp)
             nu = (NINT(nb_tmp)-1)*3+NINT(j_tmp)
             !
             DO n3 = 1, nq3
               DO n2 = 1, nq2
                 DO n1 = 1, nq1
                   !
                   READ (69517,*) n1_tmp, n2_tmp, n3_tmp, ifc(n1,n2,n3,mu,nu)
                   IF (NINT(n1_tmp) .NE. n1 .OR. NINT(n2_tmp) .NE. n2 .OR. NINT(n3_tmp) .NE. n3) &
                       CALL errore ('ifc_load','ordering of force constant of input IFC file is wrong)',999)
                   !
                 ENDDO
               ENDDO
             ENDDO
             !
           ENDDO
         ENDDO 
       ENDDO
     ENDDO
     !
     CLOSE (69517)
     !
  ENDIF
  !
  !
  ! THL: output the force constat on coarse mesh
  !      this file can be checked with that computed by q2r.x
  !
  ! format:
  ! Cart_i, Cart_j, Atom_a, Atom_b
  ! Index_1, Index_2, Index_3, IFC
  !
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     OPEN(99999,file='BTE/EPCHECK/ifc.dat')
     WRITE(99999,'(3i3)') nq1, nq2, nq3
     !
     DO i = 1,3
       DO j = 1,3
         DO na = 1,nat
           DO nb = 1,nat
             !
             WRITE(99999,'(4i3)') i, j, na, nb
             mu = (na-1)*3+i
             nu = (nb-1)*3+j
             !
             DO n3 = 1, nq3
               DO n2 = 1, nq2
                 DO n1 = 1, nq1
                   !
                   WRITE(99999,'(3x,i3,i3,i3,f16.12)') n1, n2, n3, REAL(ifc(n1,n2,n3,mu,nu))
                   !
                 ENDDO
               ENDDO
             ENDDO
             !
           ENDDO
         ENDDO 
       ENDDO
     ENDDO
     !
     CLOSE(99999)
     !
  ENDIF
  !
  !
  END SUBROUTINE ifc_load
