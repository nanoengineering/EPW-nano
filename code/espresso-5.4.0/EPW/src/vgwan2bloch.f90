!
! add by THL : compute electron and phonon group velocity using finite difference method
!
!--------------------------------------------------------------------------
SUBROUTINE vgwan2bloch_k (xkk, xxq, vel)
!--------------------------------------------------------------------------
  !
#INCLUDE "f_defs.h"
  USE kinds,         ONLY : DP
  USE epwcom,        ONLY : nbndsub
  USE cell_base,     ONLY : at, bg, alat
  USE elph2,         ONLY : ef_m, delta_egap, nrr_k, irvec, ndegen_k, chw
  USE constants_epw, ONLY : ryd2ev, bohr2ang, twopi
#ifdef __PARA
  USE io_global,     ONLY : ionode_id
  USE mp_global,     ONLY : my_pool_id
#ENDIF
  !
  IMPLICIT NONE
  !
  !
  !
  ! constant
  REAL(KIND=DP), PARAMETER   :: dk = 1.0d-4
  ! input variable
  REAL(KIND=DP), INTENT(IN)  :: xkk(3), xxq(3)
  ! output variable
  REAL(KIND=DP), INTENT(OUT) :: vel(3,nbndsub)
  ! work variable
  INTEGER                    :: ir, ibnd
  REAL(KIND=DP)              :: eig_1(nbndsub), eig_2(nbndsub), &
                                xkk_0(3), xxq_0(3), xkq_1(3), xkq_2(3)
  ! pseudo variable
  COMPLEX(KIND=DP)           :: cuf(nbndsub,nbndsub) ! cuf will not be used further
  !
  !
  xkk_0(:) = xkk(:)
  xxq_0(:) = xxq(:)
  CALL cryst_to_cart (1, xkk_0, bg, 1) ! bring xkk from cryst to cart
  CALL cryst_to_cart (1, xxq_0, bg, 1) ! bring xxq from cryst to cart
  !
  !
  DO ir = 1, 3
     !
     xkq_1(:) = xkk_0(:) + xxq_0(:)
     xkq_2(:) = xkk_0(:) + xxq_0(:)
     !
     xkq_1(ir) = xkq_1(ir) - (dk/2.0d0)
     xkq_2(ir) = xkq_2(ir) + (dk/2.0d0)
     !
     CALL cryst_to_cart (1, xkq_1, at, -1) ! bring xkq_1 from cart to cryst
     CALL hamwan2bloch (nbndsub, nrr_k, irvec, ndegen_k, xkq_1, cuf, eig_1, chw)
     CALL cryst_to_cart (1, xkq_2, at, -1) ! bring xkq_2 from cart to cryst
     CALL hamwan2bloch (nbndsub, nrr_k, irvec, ndegen_k, xkq_2, cuf, eig_2, chw)
     !
     ! calculate electron group velocity
     DO ibnd = 1, nbndsub
        vel(ir,ibnd) = (eig_2(ibnd)-eig_1(ibnd)) / ((twopi/alat)*dk) ! [Ry*Bohr]
     ENDDO
     !
  ENDDO
  !
END SUBROUTINE vgwan2bloch_k



!--------------------------------------------------------------------------
SUBROUTINE vgwan2bloch_q (xkk, xxq, vph, ntemp, ndope)
!--------------------------------------------------------------------------
  !
#INCLUDE "f_defs.h"
  USE kinds,         ONLY : DP
  USE phcom,         ONLY : nmodes
  use epwcom, only : screen_polar
  USE cell_base,     ONLY : at, bg, alat
  USE elph2,         ONLY : ef_m, delta_egap, nrr_q, irvec, ndegen_q, chw
  USE constants_epw, ONLY : ryd2ev, bohr2ang, twopi
#ifdef __PARA
  USE io_global,     ONLY : ionode_id
  USE mp_global,     ONLY : my_pool_id
#ENDIF
  !
  IMPLICIT NONE
  !
  !
  !
  ! constant
  integer :: ntemp, ndope
  REAL(KIND=DP), PARAMETER   :: dq = 1.0d-4
  ! input variable
  REAL(KIND=DP), INTENT(IN)  :: xkk(3), xxq(3)
  ! output variable
  REAL(KIND=DP), INTENT(OUT) :: vph(3,nmodes,ntemp,ndope)
  ! work variable
  INTEGER                    :: ir, imode, itemp, idope
  REAL(KIND=DP)              :: eig_1(nmodes,ntemp,ndope), eig_2(nmodes,ntemp,ndope), &
                                xkk_0(3), xxq_0(3), xkq_1(3), xkq_2(3)
  ! pseudo variable
  COMPLEX(KIND=DP)           :: cuf(nmodes,nmodes,ntemp,ndope), v(3,nmodes,ntemp,ndope) ! will not be used further
  !
  !
  xkk_0(:) = xkk(:)
  xxq_0(:) = xxq(:)
  CALL cryst_to_cart (1, xkk_0, bg, 1) ! bring xkk from cryst to cart
  CALL cryst_to_cart (1, xxq_0, bg, 1) ! bring xxq from cryst to cart
  !
  !
  DO ir = 1, 3
     !
     xkq_1(:) = xkk_0(:) + xxq_0(:)
     xkq_2(:) = xkk_0(:) + xxq_0(:)
     !
     xkq_1(ir) = xkq_1(ir) - (dq/2.0d0)
     xkq_2(ir) = xkq_2(ir) + (dq/2.0d0)
     !
     CALL cryst_to_cart (1, xkq_1, at, -1) ! bring xkq_1 from cart to cryst
     if (screen_polar) then
        CALL dynwan2bloch_s (nmodes, nrr_q, irvec, ndegen_q, xkq_1, cuf, eig_1, v, ntemp, ndope)
     else
        CALL dynwan2bloch (nmodes, nrr_q, irvec, ndegen_q, xkq_1, cuf(:,:,1,1), eig_1(:,1,1), v(:,:,1,1))
     endif
     CALL cryst_to_cart (1, xkq_2, at, -1) ! bring xkq_2 from cart to cryst
     if (screen_polar) then
        CALL dynwan2bloch_s (nmodes, nrr_q, irvec, ndegen_q, xkq_2, cuf, eig_2, v, ntemp, ndope)
     else
        CALL dynwan2bloch (nmodes, nrr_q, irvec, ndegen_q, xkq_2, cuf(:,:,1,1), eig_2(:,1,1), v(:,:,1,1))
     endif
     !
     ! calculate phonon group velocity
     DO imode = 1, nmodes
        do itemp = 1, ntemp
           do idope = 1, ndope
              !
              IF (eig_1(imode,itemp,idope) .GT. 0.0d0) THEN
                 eig_1(imode,itemp,idope) =  SQRT(ABS(eig_1(imode,itemp,idope)))
              ELSE
                 eig_1(imode,itemp,idope) =  -SQRT(ABS(eig_1(imode,itemp,idope)))
              ENDIF
              !
              IF (eig_2(imode,itemp,idope) .GT. 0.0d0) THEN
                 eig_2(imode,itemp,idope) =  SQRT(ABS(eig_2(imode,itemp,idope)))
              ELSE
                 eig_2(imode,itemp,idope) =  -SQRT(ABS(eig_2(imode,itemp,idope)))
              ENDIF
              !
              vph(ir,imode,itemp,idope) = (eig_2(imode,itemp,idope)-eig_1(imode,itemp,idope)) / ((twopi/alat)*dq) ! [Ry*Bohr]
              !
           enddo
        enddo
     ENDDO
     !
  ENDDO
  !
END SUBROUTINE vgwan2bloch_q
