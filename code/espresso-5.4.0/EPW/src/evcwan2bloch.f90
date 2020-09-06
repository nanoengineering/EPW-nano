!---------------------------------------------------------------------------
SUBROUTINE evcwan2bloch (xk, irvec, ndegen, cufkk, nbnd, nrr_k, evc_f, ng)
  !---------------------------------------------------------------------------
  !
  ! this subroutine is used for electron self-energy calculation,
  ! in this case, the transfrom from wannier to bloch for electrons
  ! is done first
  !
#INCLUDE "f_defs.h"
  !
  USE kinds,         ONLY : DP
  USE elph2,         ONLY : epmatwp, nbnd_red, ibndmin, ibndmax
  USE constants_epw, ONLY : twopi, ci, czero, cone
  USE noncollin_module, ONLY : noncolin, npol
  USE epwcom,        ONLY : nbndsub, save_m_matw
  use elph2,         only : evwan
  USE io_files,      ONLY : tmp_dir, prefix
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)           :: irvec(3, nrr_k), ndegen(nrr_k), nbnd, nrr_k, ng
  COMPLEX(KIND=DP), INTENT(IN)  :: cufkk(nbnd, nbnd)
  REAL(KIND=DP), INTENT(IN)     :: xk(3)
  !
  COMPLEX(KIND=DP), INTENT(OUT) :: evc_f(ng, npol, nbnd_red)
  !
  INTEGER                       :: imode, ibnd, jbnd, ir, ig
  REAL(KIND=DP)                 :: rdotk
  COMPLEX(KIND=DP)              :: cfac, evctmp(nbnd, npol, ng), evc_g(nbnd)
  !
  ! save_m_matw
  CHARACTER(LEN=256)            :: filename
  !
  evctmp = czero
  !
  DO ir = 1, nrr_k
     !   
     !
     rdotk = twopi * DOT_PRODUCT(xk,DBLE(irvec(:,ir)))
     cfac = EXP(ci*rdotk) / DBLE(ndegen(ir))
     !
     DO ibnd = 1, nbnd
        evctmp(ibnd,1,1:ng) = evctmp(ibnd,1,1:ng) + cfac * evwan (ibnd,1,ir,1:ng)
        if (noncolin) &
           evctmp(ibnd,2,1:ng) = evctmp(ibnd,2,1:ng) + cfac * evwan (ibnd,2,ir,1:ng)
     ENDDO
     !
  ENDDO
  !
  !
  do ig = 1, ng
     !
     CALL zgemm ('n', 'c', 1, nbnd, nbnd, cone, evctmp(:,1,ig), &
                  1, cufkk, nbnd, czero, evc_g, 1)
     evc_f(ig,1,1:nbnd_red) = evc_g(ibndmin:ibndmax)

     if (noncolin) then
        CALL zgemm ('n', 'c', 1, nbnd, nbnd, cone, evctmp(:,2,ig), &
                     1, cufkk, nbnd, czero, evc_g, 1)
        evc_f(ig,2,1:nbnd_red) = evc_g(ibndmin:ibndmax)
     endif
     !
  enddo
  !
  !
END SUBROUTINE evcwan2bloch
