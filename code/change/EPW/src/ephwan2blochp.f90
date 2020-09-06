 !---------------------------------------------------------------------------
 SUBROUTINE ephwan2blochp (nmodes, xxq, irvec, ndegen, nrr_q, cuf, epmatf, nbnd, nrr_k)
 !---------------------------------------------------------------------------
 !
 ! even though this is for phonons, I use the same notations
 ! adopted for the electronic case (nmodes->nmodes etc)
 !
 !
#INCLUDE "f_defs.h"
  USE kinds,         ONLY : DP
  USE elph2,         ONLY : epmatwp
  USE constants_epw, ONLY : twopi, ci, czero, cone
  USE epwcom,        ONLY : nbndsub, save_m_matw
  USE io_files,      ONLY : tmp_dir, prefix
  !
  IMPLICIT NONE
  !
  !
  INTEGER, INTENT(IN)           :: nmodes, nrr_q, irvec(3,nrr_q), ndegen(nrr_q), nbnd, nrr_k
  COMPLEX(KIND=DP), INTENT(IN)  :: cuf(nmodes,nmodes)
  REAL(KIND=DP), INTENT(IN)     :: xxq(3)
  !
  COMPLEX(KIND=DP), INTENT(OUT) :: epmatf(nbnd,nbnd,nrr_k,nmodes)
  !
  INTEGER                       :: imode, ibnd, jbnd, irk, irq
  REAL(KIND=DP)                 :: rdotk
  COMPLEX(KIND=DP)              :: cfac, eptmp(nbnd,nbnd,nrr_k,nmodes)
  !
  ! save_m_matw
  COMPLEX(KIND=DP)              :: epmatwp0(nbnd,nbnd,nrr_k,nmodes) !epmatwp(nbndsub,nbndsub,nrr_k,nmodes,nrr_q)
  CHARACTER(LEN=256)            :: filename
  !
  !
  eptmp = czero
  IF (save_m_matw) THEN
     filename = TRIM(tmp_dir)//TRIM(prefix)//'.epmatwp_q'
     OPEN (91915,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=2*nbnd*nbnd*nrr_k*nmodes*DP,STATUS='old')
  ENDIF
  !
  DO irq = 1, nrr_q
     !   
     IF (save_m_matw) THEN
        READ (91915,REC=irq) epmatwp0(1:nbnd,1:nbnd,1:nrr_k,1:nmodes)
     ELSE
        epmatwp0(:,:,:,:) = epmatwp(:,:,:,:,irq)
     ENDIF
     !
     rdotk = twopi * DOT_PRODUCT(xxq,DBLE(irvec(:,irq)))
     cfac = EXP(ci*rdotk) / DBLE(ndegen(irq))
     !
     eptmp(:,:,:,:) = eptmp(:,:,:,:) + cfac * epmatwp0(:,:,:,:)
     !
  ENDDO
  !
  !
  ! un-rotate to Bloch space, fine grid
  ! epmatf(j) = sum_i eptmp(i) * uf(i,j)
  !
  CALL zgemm ('n', 'n', nbnd*nbnd*nrr_k, nmodes, nmodes, (1.0d0,0.0d0), eptmp, nbnd*nbnd*nrr_k, cuf, nmodes, &
                                                         (0.0d0,0.0d0), epmatf, nbnd*nbnd*nrr_k)
  !
  IF (save_m_matw) CLOSE (91915)
  !
END SUBROUTINE ephwan2blochp
