!---------------------------------------------------------------------------
SUBROUTINE eimpwan2bloch2 (nbnd, nrr_q, irvec, ndegen, eimpmatw, xxq, cufkk, cufkq, eimpmatf)
!---------------------------------------------------------------------------
  !
#include "f_defs.h"
  !
  use io_global, only : stdout
  USE kinds,         ONLY : DP
  USE constants_epw, ONLY : twopi, ci, czero, cone
  use epwcom, only : epcheck
  implicit none
  !
  INTEGER, INTENT(IN)           :: nbnd, nrr_q, irvec(3, nrr_q), ndegen(nrr_q)
  COMPLEX(KIND=DP), INTENT(IN)  :: eimpmatw(nbnd, nbnd, nrr_q), cufkq(nbnd, nbnd), cufkk(nbnd, nbnd)
  REAL(KIND=DP), INTENT(IN)     :: xxq(3) 
  !
  COMPLEX(KIND=DP), INTENT(OUT) :: eimpmatf (nbnd, nbnd)
  !
  INTEGER                       :: ir, ibnd, jbnd
  REAL(KIND=DP)                 :: rdotk 
  COMPLEX(KIND=DP)              :: cfac, eimptmp(nbnd, nbnd)
  !
  !
  eimpmatf = czero
  !
  DO ir = 1, nrr_q
     !
     rdotk = twopi * DOT_PRODUCT(xxq,DBLE(irvec(:,ir)))
     cfac = EXP(ci*rdotk) / DBLE(ndegen(ir))
     !
     eimpmatf(:, :) = eimpmatf(:, :) + cfac * eimpmatw(:, :, ir)
     !
  ENDDO
  !
  if (epcheck .and. (xxq(1) == 0) .and. (xxq(2) == 0) .and. (xxq(3) == 0)) then
     write(*,*)
     write(*,*) ' check in eimpwan2bloch2:'
     write(*,*) ' nbnd:', nbnd
     write(*,*) ' eimpmatw:', eimpmatw(10,10,1:5)
     write(*,*) ' eimpmatf:', eimpmatf(10,10)
  endif
  !
  CALL zgemm ('n', 'n', nbnd, nbnd, nbnd, cone, cufkq, &
       nbnd, eimpmatf (:,:), nbnd, czero, eimptmp(:,:), nbnd)
  CALL zgemm ('n', 'c', nbnd, nbnd, nbnd, cone, eimptmp(:,:), &
       nbnd, cufkk, nbnd, czero, eimpmatf(:,:), nbnd)
  !
  !
END SUBROUTINE eimpwan2bloch2
