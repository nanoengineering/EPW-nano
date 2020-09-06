!---------------------------------------------------------------------------
SUBROUTINE ephwan2bloch2 (nbnd, nrr_q, irvec, ndegen, epmatw, xxq, cuf, cufkk, cufkq, epmatf, nmodes, ntemp, ndope)
!---------------------------------------------------------------------------
  !
#include "f_defs.h"
  !
  USE kinds,         ONLY : DP
  USE constants_epw, ONLY : twopi, ci, czero, cone
  implicit none
  !
  INTEGER, INTENT(IN)           :: nbnd, nrr_q, irvec(3, nrr_q), ndegen(nrr_q), nmodes, ntemp, ndope
  COMPLEX(KIND=DP), INTENT(IN)  :: epmatw(nbnd, nbnd, nrr_q, nmodes), cufkq(nbnd, nbnd), cufkk(nbnd, nbnd), &
                                   cuf(nmodes, nmodes, ntemp, ndope)
  REAL(KIND=DP), INTENT(IN)     :: xxq(3) 
  !
  COMPLEX(KIND=DP), INTENT(OUT) :: epmatf (nbnd, nbnd, nmodes, ntemp, ndope)
  !
  INTEGER                       :: ir, imode, ibnd, jbnd, itemp, idope
  REAL(KIND=DP)                 :: rdotk 
  COMPLEX(KIND=DP)              :: cfac, eptmp(nbnd, nbnd, nmodes)
  !
  !
  epmatf = czero
  !
  DO ir = 1, nrr_q
     !
     rdotk = twopi * DOT_PRODUCT(xxq,DBLE(irvec(:,ir)))
     cfac = EXP(ci*rdotk) / DBLE(ndegen(ir))
     !
     DO imode = 1, nmodes
        epmatf(:, :, imode, 1, 1) = epmatf(:, :, imode, 1, 1) + cfac * epmatw(:, :, ir, imode)
     ENDDO
     !
  ENDDO
  !
  !
  DO imode = 1, nmodes
    !
    CALL zgemm ('n', 'n', nbnd, nbnd, nbnd, cone, cufkq, &
         nbnd, epmatf (:,:,imode,1,1), nbnd, czero, eptmp(:,:,imode), nbnd)
    CALL zgemm ('n', 'c', nbnd, nbnd, nbnd, cone, eptmp(:,:,imode), &
         nbnd, cufkk, nbnd, czero, epmatf(:,:,imode,1,1), nbnd)
    !
  ENDDO
  !
  eptmp = epmatf(:,:,:,1,1)
  !
  DO ibnd = 1, nbnd
     DO jbnd = 1, nbnd
        do itemp = 1, ntemp
           do idope = 1, ndope
              !
              CALL zgemv ('t', nmodes, nmodes, cone, cuf(:,:,itemp,idope), nmodes, eptmp(ibnd,jbnd,:), &
                          1, czero, epmatf(ibnd,jbnd,1:nmodes,itemp,idope), 1)
              !
           enddo
        enddo
     ENDDO
  ENDDO
  !
  !
END SUBROUTINE ephwan2bloch2
