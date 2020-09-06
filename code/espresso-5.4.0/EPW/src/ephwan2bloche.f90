!---------------------------------------------------------------------------
SUBROUTINE ephwan2bloche (nmodes, xkk, irvec, ndegen, nrr_q, cufkk, epmatf, nbnd, nrr_k)
  !---------------------------------------------------------------------------
  !
  ! this subroutine is used for electron self-energy calculation,
  ! in this case, the transfrom from wannier to bloch for electrons
  ! is done first
  !
#INCLUDE "f_defs.h"
  !
  USE kinds,         ONLY : DP
  USE elph2,         ONLY : epmatwp
  USE constants_epw, ONLY : twopi, ci, czero, cone
  USE epwcom,        ONLY : nbndsub, save_m_matw
  USE io_files,      ONLY : tmp_dir, prefix
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)           :: nmodes, nrr_q, irvec(3, nrr_k), ndegen(nrr_k), nbnd, nrr_k
  COMPLEX(KIND=DP), INTENT(IN)  :: cufkk(nbnd, nbnd)
  REAL(KIND=DP), INTENT(IN)     :: xkk(3)
  !
  COMPLEX(KIND=DP), INTENT(OUT) :: epmatf(nbnd, nbnd, nrr_q, nmodes)
  !
  INTEGER                       :: imode, ibnd, jbnd, irk, irq
  REAL(KIND=DP)                 :: rdotk
  COMPLEX(KIND=DP)              :: cfac, eptmp(nbnd, nbnd, nrr_q, nmodes)
  !
  ! save_m_matw
  COMPLEX(KIND=DP)              :: epmatwp0(nbndsub,nbndsub,nmodes,nrr_q) !epmatwp(nbndsub,nbndsub,nrr_k,nmodes,nrr_q)
  CHARACTER(LEN=256)            :: filename
  !
  eptmp = czero
  IF (save_m_matw) THEN
     filename = TRIM(tmp_dir)//TRIM(prefix)//'.epmatwp_k'
     OPEN (91915,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=2*nbndsub*nbndsub*nmodes*nrr_q*DP,STATUS='old')
  ENDIF
  !
  DO irk = 1, nrr_k
     !   
     IF (save_m_matw) THEN
        READ (91915,REC=irk) epmatwp0(1:nbndsub,1:nbndsub,1:nmodes,1:nrr_q)
     ELSE
        epmatwp0(:,:,:,:) = epmatwp(:,:,irk,:,:)
     ENDIF
     !
     rdotk = twopi * DOT_PRODUCT(xkk,DBLE(irvec(:,irk)))
     cfac = EXP(ci*rdotk) / DBLE(ndegen(irk))
     !
     DO ibnd = 1, nbnd
        DO jbnd = 1, nbnd
           DO irq = 1, nrr_q                                   
              eptmp(ibnd,jbnd,irq,:) = eptmp(ibnd,jbnd,irq,:) + cfac * epmatwp0(ibnd,jbnd,:,irq)
           ENDDO
        ENDDO
     ENDDO
     !
  ENDDO
  !
  IF (save_m_matw) CLOSE (91915)
  !
  epmatf = eptmp
  !
  !
END SUBROUTINE ephwan2bloche
