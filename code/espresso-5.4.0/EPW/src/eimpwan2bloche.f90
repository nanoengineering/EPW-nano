!---------------------------------------------------------------------------
SUBROUTINE eimpwan2bloche (xkk, irvec, ndegen, nrr_q, cufkk, eimpmatf, nbnd, nrr_k)
  !---------------------------------------------------------------------------
  !
  ! this subroutine is used for electron self-energy calculation,
  ! in this case, the transfrom from wannier to bloch for electrons
  ! is done first
  !
#INCLUDE "f_defs.h"
  !
  USE kinds,         ONLY : DP
  USE elph2,         ONLY : eimpmatwp
  USE constants_epw, ONLY : twopi, ci, czero, cone
  USE epwcom,        ONLY : nbndsub, save_m_matw, epcheck
  USE io_files,      ONLY : tmp_dir, prefix
  use io_global, only : stdout
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)           :: nrr_q, irvec(3, nrr_k), ndegen(nrr_k), nbnd, nrr_k
  COMPLEX(KIND=DP), INTENT(IN)  :: cufkk(nbnd, nbnd)
  REAL(KIND=DP), INTENT(IN)     :: xkk(3)
  !
  COMPLEX(KIND=DP), INTENT(OUT) :: eimpmatf(nbnd, nbnd, nrr_q)
  !
  INTEGER                       :: ibnd, jbnd, irk, irq
  REAL(KIND=DP)                 :: rdotk
  COMPLEX(KIND=DP)              :: cfac, eimptmp(nbnd, nbnd, nrr_q)
  !
  ! save_m_matw
  COMPLEX(KIND=DP)              :: eimpmatwp0(nbndsub,nbndsub,nrr_q) !eimpmatwp(nbndsub,nbndsub,nrr_k,nrr_q)
  CHARACTER(LEN=256)            :: filename
  !
  eimptmp = czero
  IF (save_m_matw) THEN
     filename = TRIM(tmp_dir)//TRIM(prefix)//'.eimpmatwp_k'
     OPEN (91915,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=2*nbndsub*nbndsub*nrr_q*DP,STATUS='old')
  ENDIF
  !
  DO irk = 1, nrr_k
     !   
     IF (save_m_matw) THEN
        READ (91915,REC=irk) eimpmatwp0(1:nbndsub,1:nbndsub,1:nrr_q)
     ELSE
        eimpmatwp0(:,:,:) = eimpmatwp(:,:,irk,:)
     ENDIF
     !
     rdotk = twopi * DOT_PRODUCT(xkk,DBLE(irvec(:,irk)))
     cfac = EXP(ci*rdotk) / DBLE(ndegen(irk))
     !
     DO ibnd = 1, nbnd
        DO jbnd = 1, nbnd
           DO irq = 1, nrr_q                                   
              eimptmp(ibnd,jbnd,irq) = eimptmp(ibnd,jbnd,irq) + cfac * eimpmatwp0(ibnd,jbnd,irq)
           ENDDO
        ENDDO
     ENDDO
     !
  ENDDO
  !
  IF (save_m_matw) CLOSE (91915)
  !
  eimpmatf = eimptmp
  !
  if (epcheck) then
     write(*,*)
     write(*,*) ' check in eimpwan2bloche:'
     write(*,*) ' eimpmatwp:', eimpmatwp(10,10,1,1:5)
     write(*,*) ' eimpmatwp0:', eimpmatwp0(10,10,1:5)
     write(*,*) ' eimptmp:', eimptmp(10,10,1:5)
     write(*,*) ' eimpmatf:', eimpmatf(10,10,1:5)
  endif
  !
  !
END SUBROUTINE eimpwan2bloche
