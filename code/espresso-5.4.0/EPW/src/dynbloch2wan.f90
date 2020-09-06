  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                        
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !
  !------------------------------------------------------------------------
  SUBROUTINE dynbloch2wan ( nmodes, nq, xk, dynq, nrr, irvec, wslen )
  !------------------------------------------------------------------------
  !
  !  From the Dynamical Matrix in Bloch representation (coarse mesh), 
  !  find the corresponding matrix in Wannier representation 
  !
  !  NOTA BENE: it seems to be very important that the matrix is kept real
  !  physically these are truely the interatomic force constants.
  !  If you use a complex matrix instead, you may get some spurious 
  !  oscillations when you interpolate the phonon dispersions.
  !
  !  Note also that the acoustic sum rule for the q=0 case has been imposed
  !  already in readmat_shuffle.f90
  !
  !---------------------------------------------------------------------
  !
#include "f_defs.h"
  USE kinds,         ONLY : DP
  USE pwcom,         ONLY : at, bg, celldm, omega
  USE ions_base,     ONLY : nat, tau
  USE phcom,         ONLY : nq1, nq2, nq3
  USE control_flags, ONLY : iverbosity
  USE elph2,         ONLY : rdw, epsi, zstar, ifc
  USE epwcom,        ONLY : lpolar
  USE constants_epw, ONLY : bohr2ang, twopi, ci, czero
#ifdef __PARA
  USE io_global,     ONLY : ionode_id
  USE mp_world,      ONLY : mpime
  USE mp,            ONLY : mp_barrier
  USE mp_global,     ONLY : inter_pool_comm, my_pool_id
#endif
  implicit none
  !
  !  input variables
  !
  integer :: nmodes, nq, nrr, irvec (3, nrr)
  ! number of branches
  ! number of qpoints
  ! number of WS points and coordinates
  complex(kind=DP) :: dynq (nmodes, nmodes, nq)
  ! dynamical matrix in bloch representation (Cartesian coordinates)
  real(kind=DP) :: xk (3, nq), wslen (nrr) 
  ! kpoint coordinates (cartesian in units of 2piba)
  ! WS vectors length (alat units)
  !
  !  output variables
  !
  ! work variables 
  !
  integer :: ik, ir
  real(kind=DP) :: rdotk, tmp
  complex(kind=DP) :: cfac
  !
  !
  ! THL: used in wsweight.f90
  REAL(KIND=DP), PARAMETER :: eps = 1.0d-6
  INTEGER, PARAMETER       :: nrwsx = 200
  ! the maximun number of nearest neighbor
  INTEGER                  ::  nrws             
  ! number of nearest neighbor
  REAL(KIND=DP)            :: atws(3,3), rws(0:3,nrwsx), r(3), r_ws(3), weight
  ! lattice vector for WS initialization
  ! nearest neighbor list, rws(0,*) = norm^2
  INTEGER                  :: i, j, k, na, nb, n1, n2, n3, m1, m2, m3, mu, nu
  REAL(KIND=DP), EXTERNAL  :: wsweight
  COMPLEX(KIND=DP)         :: ddyn_l(3,3*nat,3*nat)

  !
  !
  !  subtract the long-range term from D(q)
  IF (lpolar) THEN
     DO ik = 1, nq
        CALL rgd_blk (nq1,nq2,nq3,nat,dynq(1,1,ik),1,1,ddyn_l,xk(:,ik), &  !xk has to be in cart. coord.
                  tau,epsi,zstar,bg,omega,-1.d0)
      IF (iverbosity.eq.1) WRITE (6,'(a,i3)') "Done rigid ", ik
     ENDDO
  ENDIF
  !
  !
  !------------------------------------------------------------------------------------------------
  ! THL: EPW cannot prepare identical rdw (interatomic force constant) with that from QE due to 
  !      the subroutine wigner_seitz2 might be wrong in computing the weighting of q-points. 
  !      Insdead, the new subroutine wsinit and function wsweight in wsweight.f90 are added to 
  !      compute the weighting. The correct interatomic force constant (saved in variable ifc) 
  !      is computed for the Fourier interpolation, and rdw will no londer be used.
  !------------------------------------------------------------------------------------------------
  !
  ! impose hermiticity on dynamical matrix
  DO ik = 1,nq
    CALL trasl(dynq(:,:,ik),nat)
  END DO
  !
  ! bring xk in crystal coordinates
  CALL cryst_to_cart (nq, xk, at, -1)
  !
  ! Fourier transform
  ALLOCATE (ifc(nq1,nq2,nq3,3*nat,3*nat))
  ifc(:,:,:,:,:) = 0.d0
  !
  DO na = 1,nat
    DO nb = 1,nat
      DO i = 1,3
        DO j = 1,3
          !
          mu = (na-1)*3+i
          nu = (nb-1)*3+j
          !
          DO n1 = 1,nq1
            DO n2 = 1,nq2
              DO n3 = 1,nq3
                !
                DO ik = 1,nq
                  !
                  rdotk = twopi * ( xk(1,ik)*DBLE(n1-1) + xk(2,ik)*DBLE(n2-1) + xk(3,ik)*DBLE(n3-1) )
                  cfac = EXP(ci*rdotk)/DBLE(nq)
                  ifc(n1,n2,n3,mu,nu) = ifc(n1,n2,n3,mu,nu) + REAL(cfac*dynq(mu,nu,ik))
                  !
                END DO
                !
              END DO  
            END DO
          END DO 
          ! 
        END DO 
      END DO  
    END DO 
  END DO 
  !
  !
  CALL ifc_asr(ifc,nq1,nq2,nq3,nat)
  !
  ! copy ifc to rdw
  rdw (:,:,:) = 0.d0
  !
  atws(:,1) = at(:,1)*DBLE(nq1)
  atws(:,2) = at(:,2)*DBLE(nq2)
  atws(:,3) = at(:,3)*DBLE(nq3)
  CALL wsinit(rws,nrwsx,nrws,atws)
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
          END DO
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
              rdw (:,:,ir) = ifc(m1,m2,m3,:,:)
              !
            END IF
            !
          END IF
          !
        END DO
      END DO
    END DO  
    !
  END DO
  !
  !
  ! THL: output the force constat on coarse mesh
  !      check the consistency between ifc and rdw
  !
  ! format:
  ! Cart_i, Cart_j, Atom_a, Atom_b
  ! Index, IFC 
  !
  !   IF (my_pool_id .EQ. ionode_id) THEN
  !      !
  !      OPEN(99999,file='rdw_coarse.dat')
  !      DO i = 1,3
  !        DO j = 1,3
  !          DO na = 1,nat
  !            DO nb = 1,nat
  !              !
  !              WRITE(99999,'(i3,i3,i3,i3)') i, j, na, nb
  !              mu = (na-1)*3+i
  !              nu = (nb-1)*3+j
  !              !
  !              DO ir = 1,nrr
  !                WRITE(99999,'(2x,i4,f16.12)') ir, REAL(rdw(mu,nu,ir))
  !              END DO
  !              !
  !            END DO
  !           END DO
  !        END DO
  !      END DO
  !     !
  !   ENDIF
  !CLOSE(99999)
  !
  !
  ! bring xk back into cart coord
  !
  CALL cryst_to_cart (nq, xk, bg, 1)
  !
  !
  !  check spatial decay of dynamical matrix in Wannier basis
  !  the unit in r-space is angstrom, and I am plotting
  !  the matrix for the first mode only
  !
#ifdef __PARA
    IF (mpime.eq.ionode_id) THEN
#endif
      OPEN(unit=302,file='decay.dynmat')
      WRITE(302, '(/3x,a/)') '#Spatial decay of Dynamical matrix in Wannier basis'
      DO ir = 1, nrr
        !
        tmp =  maxval ( abs( rdw (:,:,ir)) )
        WRITE(302, *) wslen(ir) * celldm (1) * bohr2ang, tmp
        !
      ENDDO
      CLOSE(302)
#ifdef __PARA
    ENDIF
    CALL mp_barrier(inter_pool_comm)
#endif
  !
  END SUBROUTINE dynbloch2wan
  !-----------------------------------------------------
  !
  !--------------------------------------------------------------------------
  SUBROUTINE trasl(dyn, nat)
  !--------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  ! input variables
  INTEGER, INTENT(IN)             :: nat
  COMPLEX(KIND=DP), INTENT(INOUT) :: dyn(3*nat,3*nat)
  !
  ! work variables  
  INTEGER                         :: i, j, na, nb, mu, nu
  COMPLEX(KIND=DP)                :: dyn_hermi(3,3,nat,nat)
  COMPLEX(KIND=DP)                :: dyn_pseu(3,3,nat,nat)
  !
  ! THL: ensure the dynamical matrix is Hermitian
  !
  ! transform the format of dynamical matrix from EPW to QE
  DO na = 1,nat
    DO nb = 1,nat
      DO i = 1,3
        DO j = 1,3
          !
          mu = (na-1)*3+i
          nu = (nb-1)*3+j
          dyn_pseu(i,j,na,nb) = dyn(mu,nu)
          !
        END DO
      END DO
    END DO
  END DO
  !
  ! impose hermiticity on dynamical matrix
  DO i=1,3
    DO j=1,3
       DO na=1,nat
          DO nb=1,nat
            dyn_hermi(i,j,na,nb) = 0.5d0 * ( dyn_pseu(i,j,na,nb) + CONJG(dyn_pseu(j,i,nb,na)) )
          END DO
       END DO
    END DO
  END DO
  !
  ! transform the format of dynamical matrix from QE to EPW
  DO na = 1,nat
    DO nb = 1,nat
      DO i = 1,3
        DO j = 1,3
          !
          mu = (na-1)*3+i
          nu = (nb-1)*3+j
          dyn(mu,nu) = dyn_hermi(i,j,na,nb)
          !
        END DO
      END DO
    END DO
  END DO
  !
  !
  RETURN
  ! 
  END SUBROUTINE trasl
  !
  !
  !--------------------------------------------------------------------------
  SUBROUTINE ifc_asr(ifc_inout, nq1, nq2, nq3, nat)
  !--------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  ! input variables
  INTEGER, INTENT(IN)             :: nq1, nq2, nq3, nat
  COMPLEX(KIND=DP), INTENT(INOUT) :: ifc_inout(nq1,nq2,nq3,3*nat,3*nat)
  !
  ! work variables
  INTEGER                         :: n1, n2, n3, i, j, na, nb, mu, nu
  REAL(KIND=DP)                   :: sum_asr
  COMPLEX(KIND=DP)                :: ifc_tmp(nq1,nq2,nq3,3,3,nat,nat)
  !
  ! THL: impose acoustic sum rule on interatomic force constant
  !
  ! transform the format of ifc from EPW to QE
  DO n1 = 1,nq1
    DO n2 = 1,nq2
      DO n3 = 1,nq3
        DO na = 1,nat
          DO nb = 1,nat
            DO i = 1,3
              DO j = 1,3
                !
                mu = (na-1)*3+i
                nu = (nb-1)*3+j
                ifc_tmp(n1,n2,n3,i,j,na,nb) = ifc_inout(n1,n2,n3,mu,nu)
                !
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO
  END DO
  !
  ! acoustic sum rule 
  DO i = 1,3
    DO j = 1,3
      DO na = 1,nat
        sum_asr = 0.0d0
        DO nb = 1,nat
          DO n1 = 1,nq1
            DO n2 = 1,nq2
              DO n3 = 1,nq3
                !
                sum_asr = sum_asr + ifc_tmp(n1,n2,n3,i,j,na,nb)
                !
              END DO
            END DO
          END DO
        END DO
        !
        ifc_tmp(1,1,1,i,j,na,na) = ifc_tmp(1,1,1,i,j,na,na) - sum_asr
        !
      END DO
    END DO
  END DO
  !
  ! transform the format of ifc from EPW QE to EPW
  DO n1 = 1,nq1
    DO n2 = 1,nq2
      DO n3 = 1,nq3
        DO na = 1,nat
          DO nb = 1,nat
            DO i = 1,3
              DO j = 1,3
                !
                mu = (na-1)*3+i
                nu = (nb-1)*3+j
                ifc_inout(n1,n2,n3,mu,nu) = ifc_tmp(n1,n2,n3,i,j,na,nb)
                !
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO
  END DO
  !
  !
  RETURN
  ! 
  END SUBROUTINE ifc_asr
  !
