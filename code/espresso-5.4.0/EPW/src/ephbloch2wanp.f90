  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !--------------------------------------------------------------------------
  SUBROUTINE ephbloch2wanp ( nbnd, nmodes, xk, nq, irvec, wslen, nrk, nrr)
  !--------------------------------------------------------------------------
  !
  !  From the EP Matrix in Electron Bloch representation (coarse mesh), 
  !  find the corresponding matrix in Phonon Wannier representation 
  !
  !--------------------------------------------------------------------------------
  !
#include "f_defs.h"
  USE kinds,         ONLY : DP
  USE pwcom,         ONLY : at, bg, celldm
  USE epwcom,        ONLY : asr_eph, save_m_matw
  USE elph2,         ONLY : epmatwp, epmatwp_asr, epmatwe
  USE constants_epw, ONLY : bohr2ang, twopi, ci, czero
#ifdef __PARA
  USE io_files,      ONLY : tmp_dir, prefix
  USE io_global,     ONLY : ionode_id
  USE mp,            ONLY : mp_barrier
  USE mp_world,      ONLY : mpime
#endif
  implicit none
  !
  !  input variables - note irvec is dimensioned with nrr_k 
  !                    (which is assumed to be larger than nrr_q)
  !
  integer :: nbnd, nrk, nmodes, nq, nrr, irvec (3, nrk)
  ! number of electronic bands
  ! number of electronic WS points
  ! number of branches
  ! number of qpoints
  ! number of WS points and coordinates
  real(kind=DP) :: xk (3, nq), wslen (nrr) 
  ! kpoint coordinates (cartesian in units of 2piba)
  ! WS vectors length (alat units)
  !
  !  output variables
  !
  ! EP matrix in electron-wannier representation and phonon-wannier
  ! representation
  !
  ! work variables 
  !
  INTEGER :: ik, iqc, ir, ire, ipol, imode, irk, irq
  REAL(KIND=DP) :: rdotk, tmp, rvec1(3), rvec2(3), len1, len2
  COMPLEX(KIND=DP) :: cfac, epmatwe_q(nbnd,nbnd,nrk,nmodes), epmatwp_q(nbnd,nbnd,nrk,nmodes), epmatwp_k(nbnd,nbnd,nmodes,nrr)
  CHARACTER(LEN=256) :: matwe_ufmt, matwp_q_ufmt, matwp_q_tmp_ufmt, matwp_k_ufmt, matwp_k_tmp_ufmt
  !
  !----------------------------------------------------------
  !  Fourier transform to go into Wannier basis
  !----------------------------------------------------------
  !
  IF (save_m_matw) THEN
     !
     matwe_ufmt = TRIM(tmp_dir)//TRIM(prefix)//'.epmatwe_q'
     OPEN (10011,FILE=matwe_ufmt,FORM='unformatted',ACCESS='direct',RECL=2*nbnd*nbnd*nrk*nmodes*DP,STATUS='old')
     !
     matwp_k_tmp_ufmt = TRIM(tmp_dir)//TRIM(prefix)//'.epmatwp_k_tmp'
     OPEN (20022,FILE=matwp_k_tmp_ufmt,FORM='unformatted',ACCESS='direct',RECL=2*nbnd*nbnd*nmodes*DP,STATUS='replace')
     !
     matwp_q_tmp_ufmt = TRIM(tmp_dir)//TRIM(prefix)//'.epmatwp_q_tmp'
     OPEN (30033,FILE=matwp_q_tmp_ufmt,FORM='unformatted',ACCESS='direct',RECL=2*nbnd*nbnd*nrk*nmodes*DP,STATUS='replace')
     !
  ENDIF
  !
  !  D (R) = (1/nk) sum_k e^{-ikR} D (k)
  !
  ! bring xk in crystal coordinates
  !
  CALL cryst_to_cart (nq, xk, at, -1)
  ! 
  DO ir = 1, nrr
    !
    epmatwp_q = (0.0d0,0.0d0)
    !
    DO iqc = 1, nq
       !
       rdotk = twopi * DOT_PRODUCT(xk(:,iqc),DBLE(irvec(:,ir)))
       cfac = EXP(-ci*rdotk) / DBLE(nq)
       !
       IF (save_m_matw) THEN
          READ (10011,REC=iqc) epmatwe_q(1:nbnd,1:nbnd,1:nrk,1:nmodes)
          epmatwp_q(:,:,:,:) = epmatwp_q(:,:,:,:) + cfac * epmatwe_q(:,:,:,:)
       ELSE
          epmatwp(:,:,:,:,ir) = epmatwp(:,:,:,:,ir) + cfac * epmatwe(:,:,:,:,iqc)
       ENDIF
       !
    ENDDO
    !
    IF (save_m_matw) THEN
       !
       DO irk = 1, nrk
          WRITE (20022,REC=(irk-1)*nrr+ir) epmatwp_q(1:nbnd,1:nbnd,irk,1:nmodes)
       ENDDO
       !
       WRITE (30033,REC=ir) epmatwp_q(1:nbnd,1:nbnd,1:nrk,1:nmodes)
       !
    ENDIF
    !
    IF (asr_eph) THEN
       !
       epmatwp_asr = czero
       !
       DO imode = 1, nmodes
          !
          ipol = MOD(imode-1,3) + 1
          !
          IF (save_m_matw) THEN
             epmatwp_asr(:,:,:,ipol) = epmatwp_asr(:,:,:,ipol) + epmatwp_q(:,:,:,imode)
          ELSE
             epmatwp_asr(:,:,:,ipol) = epmatwp_asr(:,:,:,ipol) + epmatwp(:,:,:,imode,ir)
          ENDIF
          !
       ENDDO
       !
    ENDIF
    !
    !
    !
    !  check spatial decay of e-p matrix elements in wannier basis - electrons
    !  + phonons
    !
    !  we plot: R_e, R_p, max_{m,n,nu} |g(m,n,nu;R_e,R_p)|
    !
#ifdef __PARA
    IF (mpime.eq.ionode_id) THEN
#endif
      IF (ir.eq.1) open(unit=303,file='decay.epmat_wanep',status='unknown')
      DO ire = 1, nrk
        !
        rvec1 = dble(irvec(1,ire))*at(:,1) + &
                dble(irvec(2,ire))*at(:,2) + &
                dble(irvec(3,ire))*at(:,3)
        rvec2 = dble(irvec(1,ir))*at(:,1) + &
                dble(irvec(2,ir))*at(:,2) + &
                dble(irvec(3,ir))*at(:,3)
        len1 = sqrt(rvec1(1)**2.d0+rvec1(2)**2.d0+rvec1(3)**2.d0)
        len2 = sqrt(rvec2(1)**2.d0+rvec2(2)**2.d0+rvec2(3)**2.d0)
        !
        IF (save_m_matw) THEN
           tmp = MAXVAL(ABS(epmatwp_q(:,:,ire,:)))
        ELSE
           tmp = MAXVAL(ABS(epmatwp(:,:,ire,:,ir)))
        ENDIF
        ! 
        !
        ! rvec1 : electron-electron0 distance
        ! rvec2 : phonon - electron0 distance
        !
        WRITE(303, '(5f15.10)') len1 * celldm (1) * bohr2ang, &
                                len2 * celldm (1) * bohr2ang, tmp
      ENDDO
      IF (ir.eq.nrr) close(303)
#ifdef __PARA
    ENDIF
#endif
    !
    !
  ENDDO ! loop ir
  !
  !
  !
  IF (save_m_matw) THEN
     !
     ! epmatwp_k
     matwp_k_tmp_ufmt = TRIM(tmp_dir)//TRIM(prefix)//'.epmatwp_k_tmp'
     OPEN (40044,FILE=matwp_k_tmp_ufmt,FORM='unformatted',ACCESS='direct',RECL=2*nbnd*nbnd*nmodes*DP,STATUS='old')
     !
     matwp_k_ufmt = TRIM(tmp_dir)//TRIM(prefix)//'.epmatwp_k'
     OPEN (50055,FILE=matwp_k_ufmt,FORM='unformatted',ACCESS='direct',RECL=2*nbnd*nbnd*nmodes*nrr*DP,STATUS='replace')
     !
     DO irk = 1, nrk
        !
        DO irq = 1, nrr
           READ (40044,REC=irq+(irk-1)*nrr) epmatwp_k(1:nbnd,1:nbnd,1:nmodes,irq)
           !
           IF (asr_eph) THEN
              DO imode = 1, nmodes
                 ipol = MOD(imode-1,3) + 1
                 epmatwp_k(:,:,imode,irq) = epmatwp_k(:,:,imode,irq) - (3.0d0/nrr/nmodes) * epmatwp_asr(:,:,irk,ipol)
              ENDDO
           ENDIF
           !
        ENDDO
        !
        WRITE(50055,REC=irk) epmatwp_k(1:nbnd,1:nbnd,1:nmodes,1:nrr)
        !
     ENDDO
     !
     CALL SYSTEM ('rm '//TRIM(tmp_dir)//TRIM(prefix)//'.epmatwp_k_tmp')
     !
     ! epmatwp_q
     matwp_q_tmp_ufmt = TRIM(tmp_dir)//TRIM(prefix)//'.epmatwp_q_tmp'
     OPEN (40044,FILE=matwp_q_tmp_ufmt,FORM='unformatted',ACCESS='direct',RECL=2*nbnd*nbnd*nrk*nmodes*DP,STATUS='old')
     !
     matwp_q_ufmt = TRIM(tmp_dir)//TRIM(prefix)//'.epmatwp_q'
     OPEN (90099,FILE=matwp_q_ufmt,FORM='unformatted',ACCESS='direct',RECL=2*nbnd*nbnd*nrk*nmodes*DP,STATUS='replace')
     !
     DO irq = 1, nrr
        !
        READ (40044,REC=irq) epmatwp_q(1:nbnd,1:nbnd,1:nrk,1:nmodes)
        !
        IF (asr_eph) THEN
           DO imode = 1, nmodes
              ipol = MOD(imode-1,3) + 1
              epmatwp_q(:,:,:,imode) = epmatwp_q(:,:,:,imode) - (3.0d0/nrr/nmodes) * epmatwp_asr(:,:,:,ipol)
           ENDDO
        ENDIF
        !
        WRITE (90099,REC=irq) epmatwp_q(1:nbnd,1:nbnd,1:nrk,1:nmodes)
        !
     ENDDO
     !
     CALL SYSTEM ('rm '//TRIM(tmp_dir)//TRIM(prefix)//'.epmatwp_q_tmp')
     !
     CLOSE (10011)
     CLOSE (20022)
     CLOSE (40044)
     CLOSE (90099)
     !
  ELSE
     !
     IF (asr_eph) THEN
        DO ir = 1, nrr
           DO imode = 1, nmodes
              ipol = MOD(imode-1,3) + 1
              epmatwp(:,:,:,imode,ir) = epmatwp(:,:,:,imode,ir) - (3.0d0/nrr/nmodes) * epmatwp_asr(:,:,:,ipol)
           ENDDO
        ENDDO
     ENDIF
     !
     matwp_q_ufmt = TRIM(tmp_dir)//TRIM(prefix)//'.epmatwp_q'
     OPEN (90099,FILE=matwp_q_ufmt,FORM='unformatted',ACCESS='direct',RECL=2*nbnd*nbnd*nrk*nmodes*DP,STATUS='replace')
     DO ir = 1, nrr
        WRITE (90099,REC=ir) epmatwp(1:nbnd,1:nbnd,1:nrk,1:nmodes,ir)
     ENDDO
     CLOSE (90099)
     !
  ENDIF
  !
  ! bring xk back into cart coord
  !
  CALL cryst_to_cart (nq, xk, bg, 1)
  !
  END SUBROUTINE ephbloch2wanp
