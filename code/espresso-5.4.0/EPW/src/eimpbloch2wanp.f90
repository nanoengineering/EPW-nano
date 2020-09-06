  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !--------------------------------------------------------------------------
  SUBROUTINE eimpbloch2wanp ( nbnd, xk, nq, irvec, wslen, nrk, nrr)
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
  USE elph2,         ONLY : eimpmatwp, eimpmatwe
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
  integer :: nbnd, nrk, nq, nrr, irvec (3, nrk)
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
  COMPLEX(KIND=DP) :: cfac, eimpmatwe_q(nbnd,nbnd,nrk), eimpmatwp_q(nbnd,nbnd,nrk), eimpmatwp_k(nbnd,nbnd,nrr)
  CHARACTER(LEN=256) :: matwe_ufmt, matwp_q_ufmt, matwp_q_tmp_ufmt, matwp_k_ufmt, matwp_k_tmp_ufmt
  !
  !----------------------------------------------------------
  !  Fourier transform to go into Wannier basis
  !----------------------------------------------------------
  !
  IF (save_m_matw) THEN
     !
     matwe_ufmt = TRIM(tmp_dir)//TRIM(prefix)//'.eimpmatwe_q'
     OPEN (10011,FILE=matwe_ufmt,FORM='unformatted',ACCESS='direct',RECL=2*nbnd*nbnd*nrk*DP,STATUS='old')
     !
     matwp_k_tmp_ufmt = TRIM(tmp_dir)//TRIM(prefix)//'.eimpmatwp_k_tmp'
     OPEN (20022,FILE=matwp_k_tmp_ufmt,FORM='unformatted',ACCESS='direct',RECL=2*nbnd*nbnd*DP,STATUS='replace')
     !
     matwp_q_tmp_ufmt = TRIM(tmp_dir)//TRIM(prefix)//'.eimpmatwp_q_tmp'
     OPEN (30033,FILE=matwp_q_tmp_ufmt,FORM='unformatted',ACCESS='direct',RECL=2*nbnd*nbnd*nrk*DP,STATUS='replace')
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
    eimpmatwp_q = (0.0d0,0.0d0)
    !
    DO iqc = 1, nq
       !
       rdotk = twopi * DOT_PRODUCT(xk(:,iqc),DBLE(irvec(:,ir)))
       cfac = EXP(-ci*rdotk) / DBLE(nq)
       !
       IF (save_m_matw) THEN
          READ (10011,REC=iqc) eimpmatwe_q(1:nbnd,1:nbnd,1:nrk)
          eimpmatwp_q(:,:,:) = eimpmatwp_q(:,:,:) + cfac * eimpmatwe_q(:,:,:)
       ELSE
          eimpmatwp(:,:,:,ir) = eimpmatwp(:,:,:,ir) + cfac * eimpmatwe(:,:,:,iqc)
       ENDIF
       !
    ENDDO
    !
    IF (save_m_matw) THEN
       !
       DO irk = 1, nrk
          WRITE (20022,REC=(irk-1)*nrr+ir) eimpmatwp_q(1:nbnd,1:nbnd,irk)
       ENDDO
       !
       WRITE (30033,REC=ir) eimpmatwp_q(1:nbnd,1:nbnd,1:nrk)
       !
    ENDIF
    !
    !
    !  check spatial decay of e-imp matrix elements in wannier basis - electrons
    !  + phonons
    !
    !  we plot: R_e, R_p, max_{m,n,nu} |g(m,n,nu;R_e,R_p)|
    !
#ifdef __PARA
    IF (mpime.eq.ionode_id) THEN
#endif
      IF (ir.eq.1) open(unit=303,file='decay.eimpmat_wanep',status='unknown')
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
           tmp = MAXVAL(ABS(eimpmatwp_q(:,:,ire)))
        ELSE
           tmp = MAXVAL(ABS(eimpmatwp(:,:,ire,ir)))
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
     ! eimpmatwp_k
     matwp_k_tmp_ufmt = TRIM(tmp_dir)//TRIM(prefix)//'.eimpmatwp_k_tmp'
     OPEN (40044,FILE=matwp_k_tmp_ufmt,FORM='unformatted',ACCESS='direct',RECL=2*nbnd*nbnd*DP,STATUS='old')
     !
     matwp_k_ufmt = TRIM(tmp_dir)//TRIM(prefix)//'.eimpmatwp_k'
     OPEN (50055,FILE=matwp_k_ufmt,FORM='unformatted',ACCESS='direct',RECL=2*nbnd*nbnd*nrr*DP,STATUS='replace')
     !
     DO irk = 1, nrk
        !
        DO irq = 1, nrr
           READ (40044,REC=irq+(irk-1)*nrr) eimpmatwp_k(1:nbnd,1:nbnd,irq)
        ENDDO
        !
        WRITE(50055,REC=irk) eimpmatwp_k(1:nbnd,1:nbnd,1:nrr)
        !
     ENDDO
     !
     CALL SYSTEM ('rm '//TRIM(tmp_dir)//TRIM(prefix)//'.eimpmatwp_k_tmp')
     !
     ! eimpmatwp_q
     matwp_q_tmp_ufmt = TRIM(tmp_dir)//TRIM(prefix)//'.eimpmatwp_q_tmp'
     OPEN (40044,FILE=matwp_q_tmp_ufmt,FORM='unformatted',ACCESS='direct',RECL=2*nbnd*nbnd*nrk*DP,STATUS='old')
     !
     matwp_q_ufmt = TRIM(tmp_dir)//TRIM(prefix)//'.eimpmatwp_q'
     OPEN (90099,FILE=matwp_q_ufmt,FORM='unformatted',ACCESS='direct',RECL=2*nbnd*nbnd*nrk*DP,STATUS='replace')
     !
     DO irq = 1, nrr
        !
        READ (40044,REC=irq) eimpmatwp_q(1:nbnd,1:nbnd,1:nrk)
        !
        WRITE (90099,REC=irq) eimpmatwp_q(1:nbnd,1:nbnd,1:nrk)
        !
     ENDDO
     !
     CALL SYSTEM ('rm '//TRIM(tmp_dir)//TRIM(prefix)//'.eimpmatwp_q_tmp')
     !
     CLOSE (10011)
     CLOSE (20022)
     CLOSE (40044)
     CLOSE (90099)
     !
  ELSE
     !
     matwp_q_ufmt = TRIM(tmp_dir)//TRIM(prefix)//'.eimpmatwp_q'
     OPEN (90099,FILE=matwp_q_ufmt,FORM='unformatted',ACCESS='direct',RECL=2*nbnd*nbnd*nrk*DP,STATUS='replace')
     DO ir = 1, nrr
        WRITE (90099,REC=ir) eimpmatwp(1:nbnd,1:nbnd,1:nrk,ir)
     ENDDO
     CLOSE (90099)
     !
  ENDIF
  !
  ! bring xk back into cart coord
  !
  CALL cryst_to_cart (nq, xk, bg, 1)
  !
  END SUBROUTINE eimpbloch2wanp
