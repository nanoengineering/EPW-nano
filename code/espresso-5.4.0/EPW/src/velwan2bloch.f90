  !
  ! add by THL for ADP
  !
  ! Calculate the electron group velocity on fine mech for each input k+q in turn
  !
  !--------------------------------------------------------------------------
  SUBROUTINE velwan2bloch (xkk, xxq, eig, chw, nrr, irvec, ndegen, vel)
  !--------------------------------------------------------------------------
  !
#INCLUDE "f_defs.h"
  USE mp_global, ONLY : my_pool_id
  USE io_global, ONLY : ionode_id
  USE kinds,     ONLY : DP
  use pwcom,     only : ef
  USE epwcom,    ONLY : nbndsub
  USE cell_base, ONLY : at, bg, alat
  !
  IMPLICIT NONE
  !
  ! constant
  REAL(KIND=DP), PARAMETER :: tPI        = 2.d0*3.1415926535897932d0, &
                              Ry2eV      = 13.6058d0, &                        ! Ry -> eV
                              Ry2Hz      = Ry2eV*8065.5444d0*2.99793d10*tPI, & ! Ry -> eV -> cm^-1 -> Hz
                              Bohr2m     = 5.2917721092d-11
  complex(kind=DP), parameter :: ci = (0.d0,1.d0), czero = (0.d0, 0.d0)
  !
  ! input variable
  REAL(KIND=DP) :: xkk(3), xxq(3), eig(nbndsub)
  INTEGER :: nrr, irvec(3, nrr), ndegen(nrr)     ! irvec and ndegen will be used in hamwan2bloch.f90
  COMPLEX(KIND=DP) :: chw(nbndsub, nbndsub, nrr) ! chw will be used in hamwan2bloch.f90
  !
  ! output variable
  REAL(KIND=DP) :: vel(3,nbndsub)
  !
  ! work variable
  REAL(KIND=DP) :: eig_0(nbndsub), eig_1(nbndsub), eig_2(nbndsub), &
                   xkq_0(3), xkq_1(3), xkq_2(3), xkq_norm, dxkq(3,3), dxkq_norm, &
                   xkk_0(3), xxq_0(3), d_xkq0
  INTEGER :: ibnd, idir
  !
  ! pseudo variable
  COMPLEX(KIND=DP) :: cuf(nbndsub, nbndsub)  ! cuf will not be used further
  !
  d_xkq0 = 5.0d-4
  !
  ! copy the unshifted xkq and eigen value to xkq_0 and eig_0
  xkk_0(:) = xkk(:)
  xxq_0(:) = xxq(:)
  !
  CALL cryst_to_cart (1, xkk_0, bg, 1) ! bring xkk from cryst to cart
  CALL cryst_to_cart (1, xxq_0, bg, 1) ! bring xxq from cryst to cart
  !
  xkq_0(:) = xkk_0(:) + xxq_0(:)
  eig_0(:) = eig(:)
  !
  ! calculate tiny back-shift xkq_1
  ! Gamma point will not be integrated in selfen_elec2.f90
  IF ( (xkk_0(1) .EQ. 0.0d0 .AND. xkk_0(2) .EQ. 0.0d0 .AND. xkk_0(3) .EQ. 0.0d0) .AND. &
       (xxq_0(1) .EQ. 0.0d0 .AND. xxq_0(2) .EQ. 0.0d0 .AND. xxq_0(3) .EQ. 0.0d0) ) THEN
     !
     xkq_1(:) = xkq_0(:)
     xkq_2(:) = xkq_0(:)
     eig_1(:) = eig_0(:)
     eig_2(:) = eig_0(:)
     DO ibnd = 1, nbndsub
        vel(:,ibnd) = 0.d0
     ENDDO
     !
  ELSEIF ( (xkk_0(1) .NE. 0.0d0 .OR. xkk_0(2) .NE. 0.0d0 .OR. xkk_0(3) .NE. 0.0d0) .AND. &
           (xxq_0(1) .EQ. 0.0d0 .AND. xxq_0(2) .EQ. 0.0d0 .AND. xxq_0(3) .EQ. 0.0d0) ) THEN
     !
     xkq_norm = SQRT(DOT_PRODUCT(xkq_0,xkq_0)) 
     dxkq = 0.d0
     dxkq(1,1) = d_xkq0
     dxkq(2,2) = d_xkq0
     dxkq(3,3) = d_xkq0
     do idir = 1, 3
        xkq_1(:) = xkq_0(:) - dxkq(:,idir)
        xkq_2(:) = xkq_0(:) + dxkq(:,idir)
        dxkq_norm = 2.d0*(tPI/alat)*SQRT(DOT_PRODUCT(dxkq(:,idir),dxkq(:,idir)))  ! in units [Bohr]
        !
        ! get the eigen value of back-shift xkq_1
        CALL cryst_to_cart (1, xkq_1, at, -1)                                   ! bring xkq_1 from cart to cryst
        CALL hamwan2bloch (nbndsub, nrr, irvec, ndegen, xkq_1, cuf, eig_1, chw)
        CALL cryst_to_cart (1, xkq_2, at, -1)                                   ! bring xkq_2 from cart to cryst
        CALL hamwan2bloch (nbndsub, nrr, irvec, ndegen, xkq_2, cuf, eig_2, chw)
        !
        ! calculate electron group velocity
        DO ibnd = 1, nbndsub
           vel(idir,ibnd) = (eig_2(ibnd)-eig_1(ibnd)) / dxkq_norm ! [Ry*Bohr]
        ENDDO
     enddo
     !
  ELSE
     ! for below, see ep-coupling/espresso-4.0.3/epw-3.2.8/src
     !
  ENDIF
  !
  CALL cryst_to_cart (1, xkq_0, at, -1)    ! bring xkq_0 from cart to cryst
  CALL cryst_to_cart (1, xkk_0, at, -1)    ! bring xkk_0 from cart to cryst
  CALL cryst_to_cart (1, xxq_0, at, -1)    ! bring xxq_0 from cart to cryst
  !
  ! output and check
  IF (my_pool_id .EQ. ionode_id) THEN
     !
  ENDIF
  !
  END SUBROUTINE velwan2bloch
