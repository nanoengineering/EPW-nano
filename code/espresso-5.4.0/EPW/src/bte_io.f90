!-------------------------------------------------------------------------------
SUBROUTINE kq_load ()
!-------------------------------------------------------------------------------
!
#INCLUDE "f_defs.h"
  USE kinds,     ONLY : DP
  USE epwcom,    ONLY : bte, nkf1, nkf2, nkf3, nqf1, nqf2, nqf3, phdrag
  USE bte_var
#ifdef __PARA
  USE mp,        ONLY : mp_bcast, mp_barrier
  USE io_global, ONLY : ionode_id, stdout
  USE mp_global, ONLY : my_pool_id, inter_pool_comm
#ENDIF
  !

  IMPLICIT NONE
  !
  INTEGER       :: ik, iq, isym, ir, is, it
  REAL(KIND=DP) :: symmat_pt(3,3), non_r, n_rot_tmp
  !
  !
  WRITE (stdout,'(/5x,a)') 'Load information of kq-mesh'
  !
  !
  ! k and q number
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     OPEN (99999,FILE='BTE/META/nk_ful',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='old')
     READ (99999,REC=1) nk_ful
     CLOSE (99999)
     !
     OPEN (99999,FILE='BTE/META/nq_ful',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='old')
     READ (99999,REC=1) nq_ful
     CLOSE (99999)
     !
     OPEN (99999,FILE='BTE/META/nk_irr',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='old')
     READ (99999,REC=1) nk_irr
     CLOSE (99999)
     !
     OPEN (99999,FILE='BTE/META/nq_irr',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='old')
     READ (99999,REC=1) nq_irr
     CLOSE (99999)
     !
     OPEN (9999,FILE='BTE/META/n_rot',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='old')
     READ (9999,REC=1) n_rot
     CLOSE (9999)
     !
  ENDIF
#ifdef __PARA
  CALL mp_bcast (nk_ful,ionode_id,inter_pool_comm)
  CALL mp_bcast (nk_irr,ionode_id,inter_pool_comm)
  CALL mp_bcast (nq_ful,ionode_id,inter_pool_comm)
  CALL mp_bcast (nq_irr,ionode_id,inter_pool_comm)
  CALL mp_bcast (n_rot,ionode_id,inter_pool_comm)
  CALL mp_barrier (inter_pool_comm)
#ENDIF
  !
  !
  IF (ALLOCATED (symmat_lat)) DEALLOCATE (symmat_lat)
  ALLOCATE (symmat_lat(3,3,n_rot))
  symmat_lat = 0.0d0
  !
  IF (ALLOCATED (ful2irr)) DEALLOCATE (ful2irr)
  ALLOCATE (ful2irr(nk_ful))
  ful2irr = 0
  !
  IF (ALLOCATED (irr2ful)) DEALLOCATE (irr2ful)
  ALLOCATE (irr2ful(nk_irr))
  irr2ful = 0
  !
  IF (ALLOCATED (xkf_irr)) DEALLOCATE (xkf_irr)
  ALLOCATE (xkf_irr(3,nk_irr))
  xkf_irr = 0.0d0
  !
  IF (ALLOCATED (ful2irr_q)) DEALLOCATE (ful2irr_q)
  ALLOCATE (ful2irr_q(nq_ful))
  ful2irr_q = 0
  !
  IF (ALLOCATED (xqf_ful)) DEALLOCATE (xqf_ful)
  ALLOCATE (xqf_ful(3,nq_ful))
  xqf_ful = 0.0d0
  !
  IF (phdrag) THEN
     !
     IF (ALLOCATED (xqf_irr)) DEALLOCATE (xqf_irr)
     ALLOCATE (xqf_irr(3,nq_irr))
     xqf_irr = 0.0d0
     !
     IF (ALLOCATED (irr2ful_q)) DEALLOCATE (irr2ful_q)
     ALLOCATE (irr2ful_q(nk_irr))
     irr2ful_q = 0
     !
     IF (ALLOCATED (xkf_ful)) DEALLOCATE (xkf_ful)
     ALLOCATE (xkf_ful(3,nk_ful))
     xkf_ful = 0.0d0
     !
  ENDIF
  !
  IF (bte .EQ. 2) THEN
     !
     IF (ALLOCATED (xqf_irr)) DEALLOCATE (xqf_irr)
     ALLOCATE (xqf_irr(3,nq_irr))
     xqf_irr = 0.0d0
     !
     IF (ALLOCATED (xkf_ful)) DEALLOCATE (xkf_ful)
     ALLOCATE (xkf_ful(3,nk_ful))
     xkf_ful = 0.0d0
     !
  ENDIF
  ! k and q number
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     OPEN (9999,FILE='BTE/META/symmat_lat',FORM='unformatted',ACCESS='direct',RECL=9*DP,STATUS='old')
     DO isym = 1, n_rot
        READ (9999,REC=isym) symmat_lat(1,1:3,isym), symmat_lat(2,1:3,isym), symmat_lat(3,1:3,isym)
     ENDDO
     CLOSE (9999)
     !
     OPEN (9999,FILE='BTE/META/ful2irr',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='old')
     DO ik = 1, nk_ful
        READ (9999,REC=ik) ful2irr(ik)
     ENDDO
     CLOSE (9999)
     !
     OPEN (9999,FILE='BTE/META/irr2ful',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='old')
     DO ik = 1, nk_irr
        READ (9999,REC=ik) irr2ful(ik)
     ENDDO
     CLOSE (9999)
     !
     OPEN (9999,FILE='BTE/META/xkf_irr_cryst',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='old')
     DO ik = 1, nk_irr
        READ (9999,REC=ik) xkf_irr(1:3,ik)
     ENDDO
     CLOSE (9999)
     !
     OPEN (9999,FILE='BTE/META/ful2irr_q',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='old')
     DO iq = 1, nq_ful
        READ (9999,REC=iq) ful2irr_q(iq)
     ENDDO
     CLOSE (9999)
     !
     OPEN (9999,FILE='BTE/META/xqf_ful_cryst',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='old')
     DO iq = 1, nq_ful
        READ (9999,REC=iq) xqf_ful(1:3,iq)
     ENDDO
     CLOSE (9999)
     !
     IF (phdrag) THEN
        !
        OPEN (9999,FILE='BTE/META/xqf_irr_cryst',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='old')
        DO iq = 1, nq_irr
           READ (9999,REC=iq) xqf_irr(1:3,iq)
        ENDDO
        CLOSE (9999)
        !
        OPEN (9999,FILE='BTE/META/irr2ful_q',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='old')
        DO iq = 1, nq_irr
           READ (9999,REC=iq) irr2ful_q(iq)
        ENDDO
        CLOSE (9999)
        !
        OPEN (9999,FILE='BTE/META/xkf_ful_cryst',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='old')
        DO ik = 1, nk_ful
           READ (9999,REC=ik) xkf_ful(1:3,ik)
        ENDDO
        CLOSE (9999)
        !
     ENDIF
     !
     IF (bte .EQ. 2) THEN
        !
        OPEN (9999,FILE='BTE/META/xqf_irr_cryst',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='old')
        DO iq = 1, nq_irr
           READ (9999,REC=iq) xqf_irr(1:3,iq)
        ENDDO
        CLOSE (9999)
        !
        OPEN (9999,FILE='BTE/META/xkf_ful_cryst',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='old')
        DO ik = 1, nk_ful
           READ (9999,REC=ik) xkf_ful(1:3,ik)
        ENDDO
        CLOSE (9999)
        !
     ENDIF
     !
  ENDIF
#ifdef __PARA
  CALL mp_bcast (symmat_lat,ionode_id,inter_pool_comm)
  CALL mp_bcast (ful2irr,ionode_id,inter_pool_comm)
  CALL mp_bcast (irr2ful,ionode_id,inter_pool_comm)
  CALL mp_bcast (xkf_irr,ionode_id,inter_pool_comm)
  CALL mp_bcast (ful2irr_q,ionode_id,inter_pool_comm)
  CALL mp_bcast (xqf_ful,ionode_id,inter_pool_comm)
  IF (phdrag) CALL mp_bcast (xqf_irr,ionode_id,inter_pool_comm)
  IF (phdrag) CALL mp_bcast (irr2ful_q,ionode_id,inter_pool_comm)
  IF (phdrag) CALL mp_bcast (xkf_ful,ionode_id,inter_pool_comm)
  IF (bte .EQ. 2) CALL mp_bcast (xqf_irr,ionode_id,inter_pool_comm)
  IF (bte .EQ. 2) CALL mp_bcast (xkf_ful,ionode_id,inter_pool_comm)
  CALL mp_barrier (inter_pool_comm)
#ENDIF
  !
  !
END SUBROUTINE kq_load



!-------------------------------------------------------------------------------
SUBROUTINE kq_red_load ()
!-------------------------------------------------------------------------------
!
#INCLUDE "f_defs.h"
  USE kinds,     ONLY : DP
  USE epwcom,    ONLY : nkf1, nkf2, nkf3, nqf1, nqf2, nqf3, phdrag
  USE bte_var
#ifdef __PARA
  USE mp,        ONLY : mp_bcast, mp_barrier
  USE io_global, ONLY : ionode_id, stdout
  USE mp_global, ONLY : my_pool_id, inter_pool_comm
#ENDIF
  !

  IMPLICIT NONE
  !
  INTEGER       :: ik, iq, isym, ir, is, it
  REAL(KIND=DP) :: symmat_pt(3,3), non_r, n_rot_tmp
  !
  !
  IF (ALLOCATED (ful2irr) .AND. .NOT. phdrag) DEALLOCATE (ful2irr) ! save memory
  !
  WRITE (stdout,'(5x,a)') 'Load information of reduced kq-mesh'
  !
  ! k and q number
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     OPEN (99999,FILE='BTE/META/nk_ful_red',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='old')
     READ (99999,REC=1) nk_ful_red
     CLOSE (99999)
     !
     OPEN (99999,FILE='BTE/META/nk_irr_red',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='old')
     READ (99999,REC=1) nk_irr_red
     CLOSE (99999)
     !
  ENDIF
#ifdef __PARA
  CALL mp_bcast (nk_ful_red,ionode_id,inter_pool_comm)
  CALL mp_bcast (nk_irr_red,ionode_id,inter_pool_comm)
  CALL mp_barrier (inter_pool_comm)
#ENDIF
  !
  !
  IF (ALLOCATED (rful2ful)) DEALLOCATE (rful2ful)
  ALLOCATE (rful2ful(nk_ful_red))
  rful2ful = 0
  !
  IF (ALLOCATED (ful2rful)) DEALLOCATE (ful2rful)
  ALLOCATE (ful2rful(nk_ful))
  ful2rful = 0
  !
  IF (ALLOCATED (rirr2irr)) DEALLOCATE (rirr2irr)
  ALLOCATE (rirr2irr(nk_irr_red))
  rirr2irr = 0
  !
  IF (ALLOCATED (rful2rirr)) DEALLOCATE (rful2rirr)
  ALLOCATE (rful2rirr(nk_ful_red))
  rful2rirr = 0
  !
  ! k and q number
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     OPEN (9999,FILE='BTE/META/rful2ful',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='old')
     DO ik = 1, nk_ful_red
        READ (9999,REC=ik) rful2ful(ik)
     ENDDO
     CLOSE (9999)
     !
     OPEN (9999,FILE='BTE/META/ful2rful',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='old')
     DO ik = 1, nk_ful
        READ (9999,REC=ik) ful2rful(ik)
     ENDDO
     CLOSE (9999)
     !
     OPEN (9999,FILE='BTE/META/rful2rirr',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='old')
     DO ik = 1, nk_ful_red
        READ (9999,REC=ik) rful2rirr(ik)
     ENDDO
     CLOSE (9999)
     !
     OPEN (9999,FILE='BTE/META/rirr2irr',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='old')
     DO ik = 1, nk_irr_red
        READ (9999,REC=ik) rirr2irr(ik)
     ENDDO
     CLOSE (9999)
     !
     !
  ENDIF
#ifdef __PARA
  CALL mp_bcast (rful2ful,ionode_id,inter_pool_comm)
  CALL mp_bcast (ful2rful,ionode_id,inter_pool_comm)
  CALL mp_bcast (rirr2irr,ionode_id,inter_pool_comm)
  CALL mp_bcast (rful2rirr,ionode_id,inter_pool_comm)
  CALL mp_barrier (inter_pool_comm)
#ENDIF
  !
  !
END SUBROUTINE kq_red_load


!-------------------------------------------------------------------------------
SUBROUTINE meta_save ()
!-------------------------------------------------------------------------------
!
! write the needed properties
!
#INCLUDE "f_defs.h"
  USE kinds,     ONLY : DP
  USE phcom,     ONLY : nmodes
  USE pwcom,     ONLY : ef 
  USE epwcom,    ONLY : nbndsub, bte, nkf1, nkf2, nkf3, nqf1, nqf2, nqf3, neptemp, nepdope, efermi_read, &
                        fermi_energy, phdrag, eimp_mode, alloy_pot
  USE elph2,     ONLY : ibndmin, ibndmax, nbnd_red, ef_epw, dope_ef, &
                        etf_all, vel_all, wf_ful, vph_ful, sigmai_mode_all_abs, sigmai_mode_all_emi, &
                        n_hole, n_elec, n_intr, cbnd_emin, vbnd_emax, cfsthick, vfsthick, wf_irr, &
                        vph_irr, wf_all, vph_all, sigmai_mode_all_ela_intra, sigmai_mode_all_ela_inter, &
                        sigmai_mode_all_alloy_intra, sigmai_mode_all_alloy_inter, &
                        sigmai_mode_all_intra, sigmai_mode_all_inter
  USE bte_var
#ifdef __PARA
  USE mp,            ONLY : mp_bcast, mp_barrier
  USE io_global,     ONLY : ionode_id, stdout
  USE mp_global,     ONLY : my_pool_id, inter_pool_comm
#ENDIF
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, iq, ibnd, imode, ir, itemp, idope, irec_sigmai
  !
  !
  WRITE (stdout,'(/5x,a)') 'Save files for BTE calculation'
  !
  !
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     IF (bte .GE. -100) THEN
        !
        OPEN (99999,FILE='BTE/META/parameters',FORM='unformatted',ACCESS='direct',RECL=(3*4)+(4*DP)+(4*neptemp*nepdope*DP),STATUS='replace')
        WRITE (99999,REC=1) ibndmin, ibndmax, nbnd_red, vbnd_emax, cbnd_emin, vfsthick, cfsthick, &
                            ef_epw(1:neptemp,1:nepdope), n_hole(1:neptemp,1:nepdope), n_elec(1:neptemp,1:nepdope), n_intr(1:neptemp,1:nepdope)
        CLOSE (99999)
        !
     ENDIF
     !
     IF (bte .EQ. 0 .OR. bte .EQ. 1 .OR. bte .EQ. 19 .OR. bte .EQ. 18) THEN
        !
        OPEN (99999,FILE='BTE/META/sigmai_mode',FORM='unformatted',ACCESS='direct',RECL=4*DP,STATUS='replace')
        if (eimp_mode > 0) &
           OPEN (99998,FILE='BTE/META/sigmai_eimp_mode',FORM='unformatted',ACCESS='direct',RECL=DP,STATUS='replace')
           IF (alloy_pot) &
           OPEN (77777,FILE='BTE/META/sigmai_alel_mode',FORM='unformatted',ACCESS='direct',RECL=DP,STATUS='replace')
        DO itemp = 1, neptemp
           DO idope = 1, nepdope
              DO ik = 1, nk_irr_red
                 DO ibnd = 1, nbnd_red
                    DO imode = 1, nmodes
                       WRITE (99999,REC=(itemp-1)*nepdope*nk_irr_red*nbnd_red*nmodes+(idope-1)*nk_irr_red*nbnd_red*nmodes+(ik-1)*nbnd_red*nmodes+(ibnd-1)*nmodes+imode) &
                             sigmai_mode_all_abs(itemp,idope,imode,ibnd,ik), sigmai_mode_all_emi(itemp,idope,imode,ibnd,ik), &
sigmai_mode_all_intra(itemp,idope,imode,ibnd,ik), sigmai_mode_all_inter(itemp,idope,imode,ibnd,ik)

                       irec_sigmai = (itemp-1)*nepdope*nk_irr_red*nbnd_red*2 + &
                                     (idope-1)*nk_irr_red*nbnd_red*2 + &
                                     (ik-1)*nbnd_red*2 + ibnd*2 - 1
                       if ((eimp_mode > 0) .and. (imode == 1)) then
                          WRITE (99998,REC=irec_sigmai)   sigmai_mode_all_ela_intra(itemp,idope,ibnd,ik)
                          WRITE (99998,REC=irec_sigmai+1) sigmai_mode_all_ela_inter(itemp,idope,ibnd,ik)
                          IF (alloy_pot) THEN
                          WRITE (77777,REC=irec_sigmai)   sigmai_mode_all_alloy_intra(itemp,idope,ibnd,ik)
                          WRITE (77777,REC=irec_sigmai+1) sigmai_mode_all_alloy_inter(itemp,idope,ibnd,ik)
                          ENDIF
                       endif

                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        CLOSE (99999)
        if (eimp_mode > 0) close(99998)
        IF (alloy_pot)  close (77777)
        !
     ENDIF
     !
     IF (bte .EQ. 0 .OR. bte .EQ. 1 .OR. bte .EQ. 3 .OR. bte .EQ. -1 .OR. bte .EQ. 19 .OR. bte .EQ. 18) THEN
        !
        OPEN (99999,FILE='BTE/META/electron_irr_red',FORM='unformatted',ACCESS='direct',RECL=4*DP,STATUS='replace')
        DO ik = 1, nk_irr_red
           DO ibnd = 1, nbndsub
              WRITE (99999,REC=(ik-1)*nbndsub+ibnd) etf_all(ibnd,ik), vel_all(1:3,ibnd,ik)
           ENDDO
        ENDDO
        CLOSE (99999)
        !
     ENDIF
     !
     IF ((bte .EQ. 0 .OR. bte .EQ. 1 .OR. bte .EQ. 3 .OR. bte .EQ. 18) .AND. (.NOT. phdrag)) THEN
        !
        OPEN (99999,FILE='BTE/META/phonon_ful',FORM='unformatted',ACCESS='direct',RECL=4*DP,STATUS='replace')
        DO iq = 1, nq_ful
           DO imode = 1, nmodes
              ! ERROR, not consider temperature and carrier concentration
              WRITE (99999,REC=(iq-1)*nmodes+imode) wf_ful(imode,1,1,iq), vph_ful(1:3,imode,iq)
           ENDDO
        ENDDO
        CLOSE (99999)
        !
     ENDIF
     !
     IF (phdrag) THEN
        !
        OPEN (99999,FILE='BTE/META/phonon_irr_red',FORM='unformatted',ACCESS='direct',RECL=4*DP,STATUS='replace')
        DO iq = 1, nq_irr_red
           DO imode = 1, nmodes
              WRITE (99999,REC=(iq-1)*nmodes+imode) wf_all(imode,iq), vph_all(1:3,imode,iq)
           ENDDO
        ENDDO
        CLOSE (99999)
        !
     ENDIF
     !
     IF (bte .EQ. 0 .OR. bte .EQ. 1 .OR. bte .EQ. 3) THEN
        !
        OPEN (99999,FILE='BTE/META/nscat',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='replace')
        DO ik = 1, nk_irr_red
           WRITE (99999,REC=ik) nscat_all(ik)
        ENDDO
        CLOSE (99999)
        !
     ENDIF
     ! 
  ENDIF
  !
  !
  IF ( ALLOCATED (sigmai_mode_all_abs) ) DEALLOCATE (sigmai_mode_all_abs) 
  IF ( ALLOCATED (sigmai_mode_all_emi) ) DEALLOCATE (sigmai_mode_all_emi) 
  IF ( ALLOCATED (sigmai_mode_all_inter) ) DEALLOCATE (sigmai_mode_all_inter) 
  IF ( ALLOCATED (sigmai_mode_all_intra) ) DEALLOCATE (sigmai_mode_all_intra) 
  if (eimp_mode > 0) then
     IF ( ALLOCATED (sigmai_mode_all_ela_intra) ) DEALLOCATE (sigmai_mode_all_ela_intra)
     IF ( ALLOCATED (sigmai_mode_all_ela_inter) ) DEALLOCATE (sigmai_mode_all_ela_inter)
     IF ( ALLOCATED (sigmai_mode_all_alloy_inter) ) DEALLOCATE (sigmai_mode_all_alloy_inter)
     IF ( ALLOCATED (sigmai_mode_all_alloy_intra) ) DEALLOCATE (sigmai_mode_all_alloy_intra)
  endif
  IF ( ALLOCATED (etf_all) )             DEALLOCATE (etf_all)  
  IF ( ALLOCATED (vel_all) )             DEALLOCATE (vel_all)  
  IF ( ALLOCATED (wf_ful) )              DEALLOCATE (wf_ful)  
  IF ( ALLOCATED (vph_ful) )             DEALLOCATE (vph_ful) 
  IF ( ALLOCATED (wf_all) )              DEALLOCATE (wf_all)  
  IF ( ALLOCATED (vph_all) )             DEALLOCATE (vph_all) 
  IF ( ALLOCATED (wf_irr) )              DEALLOCATE (wf_irr)  
  IF ( ALLOCATED (vph_irr) )             DEALLOCATE (vph_irr) 
  IF ( ALLOCATED (nscat_all) )           DEALLOCATE (nscat_all)  
  IF ( ALLOCATED (ef_epw) )              DEALLOCATE (ef_epw)  
  IF ( ALLOCATED (n_hole) )              DEALLOCATE (n_hole)  
  IF ( ALLOCATED (n_elec) )              DEALLOCATE (n_elec)  
  IF ( ALLOCATED (n_intr) )              DEALLOCATE (n_intr)  
  !
  !
  CALL mp_barrier (inter_pool_comm)
  !
  !
END SUBROUTINE meta_save



!-------------------------------------------------------------------------------
SUBROUTINE meta_load ()
!-------------------------------------------------------------------------------
!
! load the needed properties saved by meta_save
!
#INCLUDE "f_defs.h"
  USE kinds,         ONLY : DP
  USE phcom,         ONLY : nmodes
  USE pwcom,         ONLY : ef 
  USE epwcom,        ONLY : nbndsub, bte, nkf1, nkf2, nkf3, nqf1, nqf2, nqf3, neptemp, nepdope, &
                            efermi_read, fermi_energy, phdrag, eimp_mode, alloy_pot
  USE elph2,         ONLY : ibndmin, ibndmax, nbnd_red, ef_epw, dope_ef, &
                            etf_all, vel_all, wf_ful, vph_ful, sigmai_mode_all_abs, sigmai_mode_all_emi, &
                            n_hole, n_elec, n_intr, cbnd_emin, vbnd_emax, cfsthick, vfsthick, wf_irr, &
                            vph_irr, wf_all, vph_all, sigmai_mode_all_ela_intra, sigmai_mode_all_ela_inter, &
                            sigmai_mode_all_alloy_inter, sigmai_mode_all_alloy_intra, &
                            sigmai_mode_all_inter, sigmai_mode_all_intra
  USE bte_var
  USE constants_epw, ONLY : ryd2ev, ryd2mev, rydcm1, au2m, au2s
#ifdef __PARA
  USE mp,            ONLY : mp_bcast, mp_barrier
  USE io_global,     ONLY : ionode_id, stdout
  USE mp_global,     ONLY : my_pool_id, inter_pool_comm
#ENDIF
  !
  IMPLICIT NONE
  !
  INTEGER             :: ik, iq, ibnd, imode, ir, itemp, idope, irec_sigmai
  CHARACTER(LEN=256)  :: sigmai_mode_fmt
  CHARACTER(LEN=3)    :: itemp_num, idope_num
  !
  !
  WRITE (stdout,'(/5x,a)') 'Load files for BTE calculation'
  !
  !
  ! system parameters
  ALLOCATE (ef_epw(neptemp,nepdope)) 
  ALLOCATE (n_hole(neptemp,nepdope))  
  ALLOCATE (n_elec(neptemp,nepdope))   
  ALLOCATE (n_intr(neptemp,nepdope)) 
  ! 
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     OPEN (99999,FILE='BTE/META/parameters',FORM='unformatted',ACCESS='direct',RECL=(3*4)+(4*DP)+(4*neptemp*nepdope*DP),STATUS='old')
     READ (99999,REC=1) ibndmin, ibndmax, nbnd_red, vbnd_emax, cbnd_emin, vfsthick, cfsthick, &
                        ef_epw(1:neptemp,1:nepdope), n_hole(1:neptemp,1:nepdope), n_elec(1:neptemp,1:nepdope), n_intr(1:neptemp,1:nepdope)
     CLOSE (99999) 

     !
  ENDIF
  !
#ifdef __PARA
  CALL mp_bcast (ibndmin,ionode_id,inter_pool_comm)
  CALL mp_bcast (ibndmax,ionode_id,inter_pool_comm)
  CALL mp_bcast (nbnd_red,ionode_id,inter_pool_comm)
  CALL mp_bcast (vbnd_emax,ionode_id,inter_pool_comm)
  CALL mp_bcast (cbnd_emin,ionode_id,inter_pool_comm)
  CALL mp_bcast (vfsthick,ionode_id,inter_pool_comm)
  CALL mp_bcast (cfsthick,ionode_id,inter_pool_comm)
  CALL mp_bcast (ef_epw,ionode_id,inter_pool_comm)
  CALL mp_bcast (n_hole,ionode_id,inter_pool_comm)
  CALL mp_bcast (n_elec,ionode_id,inter_pool_comm)
  CALL mp_bcast (n_intr,ionode_id,inter_pool_comm)
#ENDIF
  !
  !
  ! sigmai_mode
  IF (bte .EQ. 0 .OR. bte .EQ. 1 .OR. bte .EQ. 10 .OR. bte .EQ. 19 .OR. bte .EQ. 18) THEN
     !
  ALLOCATE (sigmai_mode_all_abs(neptemp,nepdope,nmodes,nbnd_red,nk_irr_red))
  ALLOCATE (sigmai_mode_all_emi(neptemp,nepdope,nmodes,nbnd_red,nk_irr_red))
  ALLOCATE (sigmai_mode_all_inter(neptemp,nepdope,nmodes,nbnd_red,nk_irr_red))
  ALLOCATE (sigmai_mode_all_intra(neptemp,nepdope,nmodes,nbnd_red,nk_irr_red))
  sigmai_mode_all_abs = 0.0d0
  sigmai_mode_all_emi = 0.0d0
  sigmai_mode_all_inter = 0.0d0
  sigmai_mode_all_intra = 0.0d0
     !
  if (eimp_mode > 0) then
     ALLOCATE (sigmai_mode_all_ela_intra(neptemp,nepdope,nbnd_red,nk_irr_red))
     ALLOCATE (sigmai_mode_all_ela_inter(neptemp,nepdope,nbnd_red,nk_irr_red))
     sigmai_mode_all_ela_intra = 0.0d0
     sigmai_mode_all_ela_inter = 0.0d0
     IF (alloy_pot) THEN
      ALLOCATE (sigmai_mode_all_alloy_intra(neptemp,nepdope,nbnd_red,nk_irr_red))
      ALLOCATE (sigmai_mode_all_alloy_inter(neptemp,nepdope,nbnd_red,nk_irr_red))
      sigmai_mode_all_alloy_inter = 0.0d0
      sigmai_mode_all_alloy_intra = 0.0d0
     ENDIF
  endif
     !
     IF (my_pool_id .EQ. ionode_id) THEN
        !
        OPEN (99999,FILE='BTE/META/sigmai_mode',FORM='unformatted',ACCESS='direct',RECL=4*DP,STATUS='old')

        if (eimp_mode > 0) then
           OPEN (99998,FILE='BTE/META/sigmai_eimp_mode',FORM='unformatted',ACCESS='direct',RECL=DP,STATUS='old')
           IF (alloy_pot) &
            OPEN (77777,FILE='BTE/META/sigmai_alel_mode',FORM='unformatted',ACCESS='direct',RECL=DP,STATUS='old')
        endif

        DO itemp = 1, neptemp
           DO idope = 1, nepdope
              DO ik = 1, nk_irr_red
                 DO ibnd = 1, nbnd_red
                    DO imode = 1, nmodes
                       READ (99999,REC=(itemp-1)*nepdope*nk_irr_red*nbnd_red*nmodes+(idope-1)*nk_irr_red*nbnd_red*nmodes+(ik-1)*nbnd_red*nmodes+(ibnd-1)*nmodes+imode) &
                            sigmai_mode_all_abs(itemp,idope,imode,ibnd,ik), sigmai_mode_all_emi(itemp,idope,imode,ibnd,ik), &
sigmai_mode_all_intra(itemp,idope,imode,ibnd,ik), sigmai_mode_all_inter(itemp,idope,imode,ibnd,ik)

                       irec_sigmai = (itemp-1)*nepdope*nk_irr_red*nbnd_red*2 + &
                                     (idope-1)*nk_irr_red*nbnd_red*2 + &
                                     (ik-1)*nbnd_red*2 + ibnd*2 - 1

                       if ((eimp_mode > 0) .and. (imode == 1)) then
                          READ (99998,REC=irec_sigmai)   sigmai_mode_all_ela_intra(itemp,idope,ibnd,ik)
                          READ (99998,REC=irec_sigmai+1) sigmai_mode_all_ela_inter(itemp,idope,ibnd,ik)
                       endif
                       !
                       IF (alloy_pot) THEN
                         READ (77777,REC=irec_sigmai)   sigmai_mode_all_alloy_intra(itemp,idope,ibnd,ik)
                         READ (77777,REC=irec_sigmai+1) sigmai_mode_all_alloy_inter(itemp,idope,ibnd,ik)
                       ENDIF
                       !
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        CLOSE (99999)
        if (eimp_mode > 0) close(99998)
        IF (alloy_pot) CLOSE(77777)
        !
     ENDIF
     !
#ifdef __PARA
     CALL mp_bcast (sigmai_mode_all_abs,ionode_id,inter_pool_comm)
     CALL mp_bcast (sigmai_mode_all_emi,ionode_id,inter_pool_comm)
     CALL mp_bcast (sigmai_mode_all_inter,ionode_id,inter_pool_comm)
     CALL mp_bcast (sigmai_mode_all_intra,ionode_id,inter_pool_comm)
     if (eimp_mode > 0) then
        CALL mp_bcast (sigmai_mode_all_ela_intra,ionode_id,inter_pool_comm)
        CALL mp_bcast (sigmai_mode_all_ela_inter,ionode_id,inter_pool_comm)
        IF (alloy_pot)  CALL mp_bcast (sigmai_mode_all_alloy_intra,ionode_id,inter_pool_comm)
        IF (alloy_pot)  CALL mp_bcast (sigmai_mode_all_alloy_inter,ionode_id,inter_pool_comm)
     endif
#ENDIF
     !
     !
  ENDIF
  !
  !
  ! electron properties
  IF (bte .EQ. 0 .OR. bte .EQ. 1 .OR. bte .EQ. 10 .OR. bte .EQ. 3 .OR. bte .EQ. 30 .OR. bte .EQ. -1 .OR. bte .EQ. 19 .OR. bte .EQ. 18) THEN
     !
     ALLOCATE (etf_all(nbndsub,nk_irr_red))  
     ALLOCATE (vel_all(3,nbndsub,nk_irr_red)) 
     etf_all = 0.0d0
     vel_all = 0.0d0
     !
     IF (my_pool_id .EQ. ionode_id) THEN
        !
        OPEN (99999,FILE='BTE/META/electron_irr_red',FORM='unformatted',ACCESS='direct',RECL=4*DP,STATUS='old')
        DO ik = 1, nk_irr_red
           DO ibnd = 1, nbndsub
              READ (99999,REC=(ik-1)*nbndsub+ibnd) etf_all(ibnd,ik), vel_all(1:3,ibnd,ik)
           ENDDO
        ENDDO
        CLOSE (99999)
        !
     ENDIF
     !
#ifdef __PARA
     CALL mp_bcast (etf_all,ionode_id,inter_pool_comm)
     CALL mp_bcast (vel_all,ionode_id,inter_pool_comm)
#ENDIF
     !
  ENDIF
  !
  !
  ! index of scattering event
  IF (bte .EQ. 0 .OR. bte .EQ. 1 .OR. bte .EQ. 10 .OR. bte .EQ. 3 .OR. bte .EQ. 30) THEN
     !
     ALLOCATE (nscat_all(nk_irr_red))  
     nscat_all = 0
     !
     IF (my_pool_id .EQ. ionode_id) THEN
        !
        OPEN (99999,FILE='BTE/META/nscat',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='old')
        DO ik = 1, nk_irr_red
           READ (99999,REC=ik) nscat_all(ik)
        ENDDO
        CLOSE (99999)  
        !
     ENDIF
     !
#ifdef __PARA
     CALL mp_bcast (nscat_all,ionode_id,inter_pool_comm)
#ENDIF
     !
  ENDIF
  !
  !
  ! phonon properties
  IF ((bte .EQ. 0 .OR. bte .EQ. 1 .OR. bte .EQ. 10 .OR. bte .EQ. 3 .OR. bte .EQ. 30) .AND. (.NOT. phdrag)) THEN
     !
     ALLOCATE (wf_ful(nmodes,1,1,nq_ful))  
     ALLOCATE (vph_ful(3,nmodes,nq_ful))  
     wf_ful  = 0.0d0
     vph_ful = 0.0d0
     !
     IF (my_pool_id .EQ. ionode_id) THEN
        !
        OPEN (99999,FILE='BTE/META/phonon_ful',FORM='unformatted',ACCESS='direct',RECL=4*DP,STATUS='old')
        DO iq = 1, nq_ful
           DO imode = 1, nmodes
              ! ERROR, not consider temperature and carrier concentration
              READ (99999,REC=(iq-1)*nmodes+imode) wf_ful(imode,1,1,iq), vph_ful(1:3,imode,iq)
           ENDDO
        ENDDO
        CLOSE (99999)
        !
     ENDIF 
     !
#ifdef __PARA
     CALL mp_bcast (wf_ful,ionode_id,inter_pool_comm)
     CALL mp_bcast (vph_ful,ionode_id,inter_pool_comm)
#ENDIF
     !
  ENDIF
  !
  !
  ! phonon properties
  IF (phdrag) THEN
     !
     ALLOCATE (wf_all(nmodes,nq_irr_red))  
     ALLOCATE (vph_all(3,nmodes,nq_irr_red))  
     wf_all  = 0.0d0
     vph_all = 0.0d0
     !
     IF (my_pool_id .EQ. ionode_id) THEN
        !
        OPEN (99999,FILE='BTE/META/phonon_irr_red',FORM='unformatted',ACCESS='direct',RECL=4*DP,STATUS='old')
        DO iq = 1, nq_irr_red
           DO imode = 1, nmodes
              READ (99999,REC=(iq-1)*nmodes+imode) wf_all(imode,iq), vph_all(1:3,imode,iq)
           ENDDO
        ENDDO
        CLOSE (99999)
        !
     ENDIF 
     !
#ifdef __PARA
     CALL mp_bcast (wf_all,ionode_id,inter_pool_comm)
     CALL mp_bcast (vph_all,ionode_id,inter_pool_comm)
#ENDIF
     !
  ENDIF
  !
  !
  ! output and check
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     OPEN (99999,FILE='BTE/EPCHECK/parameters.dat',STATUS='replace')
     WRITE (99999,'(a,i15)')    'ibndmin   = ', ibndmin
     WRITE (99999,'(a,i15)')    'ibndmax   = ', ibndmax
     WRITE (99999,'(a,i15)')    'nbnd_red  = ', nbnd_red
     WRITE (99999,'(a,f15.4)')  'vbnd_emax = ', vbnd_emax*ryd2ev
     WRITE (99999,'(a,f15.4)')  'cbnd_emin = ', cbnd_emin*ryd2ev
     WRITE (99999,'(a,f15.4)')  'vfsthick  = ', vfsthick*ryd2ev
     WRITE (99999,'(a,f15.4)')  'cfsthick  = ', cfsthick*ryd2ev
     DO itemp = 1, neptemp
        DO idope = 1, nepdope
           WRITE (99999,'(a,f15.4)')  'ef        = ', ef_epw(itemp,idope)*ryd2ev
           WRITE (99999,'(a,es15.6)') 'n_hole    = ', n_hole(itemp,idope)
           WRITE (99999,'(a,es15.6)') 'n_elec    = ', n_elec(itemp,idope)
           WRITE (99999,'(a,es15.6)') 'n_intr    = ', n_intr(itemp,idope)
        ENDDO
     ENDDO
     CLOSE (99999)
     !
     IF (ALLOCATED(sigmai_mode_all_abs) .EQ. .TRUE. .AND. ALLOCATED(sigmai_mode_all_emi) .EQ. .TRUE.) THEN
        DO itemp = 1, neptemp
           DO idope = 1, nepdope
              !
              WRITE(itemp_num,'(i3)') itemp
              WRITE(idope_num,'(i3)') idope
              sigmai_mode_fmt = 'BTE/sigmai_mode_T'//TRIM(ADJUSTL(itemp_num))//'_N'//TRIM(ADJUSTL(idope_num))//'.dat'
              OPEN (99999,FILE='BTE/EPCHECK/sigmai_mode.dat',STATUS='replace')
              !
              DO ik = 1, nk_irr_red
                 DO ibnd = 1, nbnd_red 
                    DO imode = 1, nmodes  
                       WRITE (99999,'(i12,i5,i5,2f15.6)') ik, ibnd+ibndmin-1, imode, sigmai_mode_all_abs(itemp,idope,imode,ibnd,ik)*ryd2mev, &
                                                                                     sigmai_mode_all_emi(itemp,idope,imode,ibnd,ik)*ryd2mev
                    ENDDO
                 ENDDO
              ENDDO
              CLOSE (99999)
              !
           ENDDO
        ENDDO
     ENDIF
     !
     IF (ALLOCATED(etf_all) .EQ. .TRUE. .AND. ALLOCATED(vel_all) .EQ. .TRUE.) THEN
        OPEN (99999,FILE='BTE/EPCHECK/electron_irr_red.dat',STATUS='replace')
        DO ik = 1, nk_irr_red
           DO ibnd = 1, nbndsub
              WRITE (99999,'(i12,i5,f15.4,3f17.4)') ik, ibnd, etf_all(ibnd,ik)*ryd2ev, vel_all(1:3,ibnd,ik)*(au2m/au2s)
           ENDDO
        ENDDO
        CLOSE (99999)
     ENDIF
     !
     IF (ALLOCATED(wf_ful) .EQ. .TRUE. .AND. ALLOCATED(vph_ful) .EQ. .TRUE.) THEN
        OPEN (99999,FILE='BTE/EPCHECK/phonon_ful.dat',STATUS='replace')
        DO iq = 1, nq_ful
           DO imode = 1, nmodes  
              ! ERROR, not consider temperature and carrier concentration
              WRITE (99999,'(i12,i5,f15.4,3f17.4)') iq, imode, wf_ful(imode,1,1,iq)*rydcm1, vph_ful(1:3,imode,iq)*(au2m/au2s)
           ENDDO
        ENDDO
        CLOSE (99999)
     ENDIF
     !
     IF (ALLOCATED(wf_all) .EQ. .TRUE. .AND. ALLOCATED(vph_all) .EQ. .TRUE.) THEN
        OPEN (99999,FILE='BTE/EPCHECK/phonon_irr_red.dat',STATUS='replace')
        DO iq = 1, nq_irr_red
           DO imode = 1, nmodes  
              WRITE (99999,'(i12,i5,f15.4,3f17.4)') iq, imode, wf_all(imode,iq)*rydcm1, vph_all(1:3,imode,iq)*(au2m/au2s)
           ENDDO
        ENDDO
        CLOSE (99999)
     ENDIF
     !
     IF (ALLOCATED(nscat_all)) THEN
        OPEN (99999,FILE='BTE/EPCHECK/nscat.dat',STATUS='replace')
        DO ik = 1, nk_irr_red
           WRITE (99999,'(i12,i12)') ik, nscat_all(ik)
        ENDDO
        CLOSE (99999)
     ENDIF
     !
  ENDIF
  !
  !
  CALL mp_barrier (inter_pool_comm)
  !
  !
END SUBROUTINE meta_load



!-------------------------------------------------------------------------------
SUBROUTINE sigmai_load (ik_irr, ibnd, jbnd, imode, nscat, iq_ind, sigmai, itemp, idope)
!-------------------------------------------------------------------------------
!
#INCLUDE "f_defs.h"
  USE kinds,        ONLY : DP
  USE bte_var
  !
  IMPLICIT NONE
  !
  !
  INTEGER, INTENT(IN)        :: ik_irr, ibnd, jbnd, imode, itemp, idope
  INTEGER, INTENT(OUT)       :: nscat, iq_ind(nq_ful)
  REAL(KIND=DP), INTENT(OUT) :: sigmai(nq_ful)
  !
  CHARACTER(LEN=256)         :: file_ufmt
  CHARACTER(LEN=12)          :: tnph
  CHARACTER(LEN=10)          :: k_num, ibnd_num, jbnd_num, imode_num, itemp_num, idope_num
  !
  !
  WRITE(k_num,'(i10)') ik_irr
  WRITE(ibnd_num,'(i10)') ibnd
  WRITE(jbnd_num,'(i10)') jbnd
  WRITE(imode_num,'(i10)') imode
  WRITE(itemp_num,'(i10)') itemp
  WRITE(idope_num,'(i10)') idope
  !
  tnph = 'T'//TRIM(ADJUSTL(itemp_num))//'_N'//TRIM(ADJUSTL(idope_num))//'_ph'//TRIM(ADJUSTL(imode_num))
  !
  file_ufmt = 'BTE/SIGMAI/'//TRIM(ADJUSTL(tnph))//'/sigmai_'//TRIM(ADJUSTL(k_num))//'_'//TRIM(ADJUSTL(ibnd_num))//'_'//TRIM(ADJUSTL(jbnd_num))
  !
  OPEN (74639,FILE=file_ufmt,FORM='unformatted',ACCESS='direct',RECL=4+nq_ful*(4+DP),STATUS='old')
  READ (74639,REC=1) nscat, iq_ind(1:nscat), sigmai(1:nscat)
  CLOSE (74639)
  !
  !
END SUBROUTINE sigmai_load


!-------------------------------------------------------------------------------
SUBROUTINE gammai_load (iq_irr, ibnd, jbnd, imode, nscat, ik_ind, gammai, itemp, idope)
!-------------------------------------------------------------------------------
!
#INCLUDE "f_defs.h"
  USE kinds,        ONLY : DP
  USE bte_var
  !
  IMPLICIT NONE
  !
  !
  INTEGER, INTENT(IN)        :: iq_irr, ibnd, jbnd, imode, itemp, idope
  INTEGER, INTENT(OUT)       :: nscat, ik_ind(nk_ful)
  REAL(KIND=DP), INTENT(OUT) :: gammai(nk_ful)
  !
  CHARACTER(LEN=256)         :: file_ufmt
  CHARACTER(LEN=12)          :: tnph
  CHARACTER(LEN=10)          :: q_num, ibnd_num, jbnd_num, imode_num, itemp_num, idope_num
  !
  !
  WRITE(q_num,'(i10)') iq_irr
  WRITE(ibnd_num,'(i10)') ibnd
  WRITE(jbnd_num,'(i10)') jbnd
  WRITE(imode_num,'(i10)') imode
  WRITE(itemp_num,'(i10)') itemp
  WRITE(idope_num,'(i10)') idope
  !
  tnph = 'T'//TRIM(ADJUSTL(itemp_num))//'_N'//TRIM(ADJUSTL(idope_num))//'_ph'//TRIM(ADJUSTL(imode_num))
  !
  file_ufmt = 'BTE/GAMMAI/'//TRIM(ADJUSTL(tnph))//'/gammai_'//TRIM(ADJUSTL(q_num))//'_'//TRIM(ADJUSTL(ibnd_num))//'_'//TRIM(ADJUSTL(jbnd_num))
  !
  OPEN (73092,FILE=file_ufmt,FORM='unformatted',ACCESS='direct',RECL=4+nk_ful*(4+DP),STATUS='old')
  READ (73092,REC=1) nscat, ik_ind(1:nscat), gammai(1:nscat)
  CLOSE (73092)
  !
  !
END SUBROUTINE gammai_load



!-------------------------------------------------------------------------------
SUBROUTINE weight_load (ik_irr, ibnd, jbnd, imode, nscat, iq_ind, weight_abs, weight_emi)
!-------------------------------------------------------------------------------
!
#INCLUDE "f_defs.h"
  USE kinds,        ONLY : DP
  USE bte_var
  !
  IMPLICIT NONE
  !
  !
  INTEGER, INTENT(IN)        :: ik_irr, ibnd, jbnd, imode
  INTEGER, INTENT(OUT)       :: nscat, iq_ind(nq_ful)
  REAL(KIND=DP), INTENT(OUT) :: weight_abs(nq_ful), weight_emi(nq_ful)
  !
  CHARACTER(LEN=256)         :: file_ufmt
  CHARACTER(LEN=12)          :: tnph
  CHARACTER(LEN=8)           :: k_num, ibnd_num, jbnd_num, imode_num
  !
  !
  WRITE(k_num,'(i8)') ik_irr
  WRITE(ibnd_num,'(i8)') ibnd
  WRITE(jbnd_num,'(i8)') jbnd
  WRITE(imode_num,'(i8)') imode
  !
  tnph = 'T1_N1_ph'//TRIM(ADJUSTL(imode_num))
  !
  file_ufmt = 'BTE/WEIGHT/'//TRIM(ADJUSTL(tnph))//'/weight_'//TRIM(ADJUSTL(k_num))//'_'//TRIM(ADJUSTL(ibnd_num))//'_'//TRIM(ADJUSTL(jbnd_num))
  !
  OPEN (99999,FILE=file_ufmt,FORM='unformatted',ACCESS='direct',RECL=4+nq_ful*(4+2*DP),STATUS='old')
  READ (99999,REC=1) nscat, iq_ind(1:nscat), weight_abs(1:nscat), weight_emi(1:nscat)
  CLOSE (99999)
  !
  !
END SUBROUTINE weight_load



!-------------------------------------------------------------------------------
SUBROUTINE Fk_save (iter, itemp, idope)
!-------------------------------------------------------------------------------
#INCLUDE "f_defs.h"
  USE kinds,        ONLY : DP
  USE phcom,        ONLY : nmodes
  USE epwcom,       ONLY : nbndsub
  USE elph2,        ONLY : ibndmin, ibndmax, nbnd_red, F_k_ful
  USE bte_var,      ONLY : nk_ful_red
#ifdef __PARA
  USE mp,           ONLY : mp_barrier, mp_bcast
  USE io_global,    ONLY : ionode_id, stdout
  USE mp_global,    ONLY : my_pool_id, inter_pool_comm
#ENDIF
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iter, itemp, idope
  !
  INTEGER             :: ik, ibnd, imode, ir
  !
  !
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     OPEN (99999,FILE='BTE/META/F_k_ful_tmp',FORM='unformatted',ACCESS='direct',RECL=DP,STATUS='replace')
     WRITE (99999,REC=1) DBLE(iter)
     WRITE (99999,REC=2) DBLE(itemp)
     WRITE (99999,REC=3) DBLE(idope)
     DO ik = 1, nk_ful_red
        DO ibnd = 1, nbnd_red
           DO ir = 1, 3
              WRITE (99999,REC=(ik-1)*nbnd_red*3+(ibnd-1)*3+ir+3) F_k_ful(ir,ibnd,ik)
           ENDDO
        ENDDO
     ENDDO
     CLOSE (99999)
     !
     CALL SYSTEM ('cp -r BTE/META/F_k_ful_tmp BTE/META/F_k_ful')
     CALL SYSTEM ('rm BTE/META/F_k_ful_tmp')
     !
  ENDIF
  !
#ifdef __PARA
  CALL mp_barrier (inter_pool_comm)
#ENDIF
  !
  !
END SUBROUTINE Fk_save



!-------------------------------------------------------------------------------
SUBROUTINE Fk_load (iter, itemp, idope)
!-------------------------------------------------------------------------------
#INCLUDE "f_defs.h"
  USE kinds,        ONLY : DP
  USE phcom,        ONLY : nmodes
  USE epwcom,       ONLY : nbndsub
  USE elph2,        ONLY : ibndmin, ibndmax, nbnd_red, F_k_ful
  USE bte_var,        ONLY : nk_ful_red
#ifdef __PARA
  USE mp,           ONLY : mp_bcast, mp_barrier
  USE io_global,    ONLY : ionode_id, stdout
  USE mp_global,    ONLY : my_pool_id, inter_pool_comm
#ENDIF
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(OUT) :: iter, itemp, idope
  !
  INTEGER              :: ik, ibnd, imode, ir, cnt
  REAL(KIND=DP)        :: iter_tmp, temp_tmp, dope_tmp
  !
  !
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     OPEN (99999,FILE='BTE/META/F_k_ful',FORM='unformatted',ACCESS='direct',RECL=DP,STATUS='old')
     !
     cnt = 0
     READ (99999,REC=1) iter_tmp
     READ (99999,REC=2) temp_tmp
     READ (99999,REC=3) dope_tmp
     DO ik = 1, nk_ful_red
        DO ibnd = 1, nbnd_red
           DO ir = 1, 3
              READ (99999,REC=(ik-1)*nbnd_red*3+(ibnd-1)*3+ir+3) F_k_ful(ir,ibnd,ik)
              cnt = cnt + 1
           ENDDO
        ENDDO
     ENDDO
     CLOSE (99999)
     !
     IF (cnt .NE. nk_ful_red*nbnd_red*3) &
        CALL errore('Fk_load','length of metafile F_k_ful is incorrect, please delete BTE/META/F_k_ful manually',1)
     !
     iter = NINT(iter_tmp) + 1
     itemp = NINT(temp_tmp)
     idope = NINT(dope_tmp)
     !
  ENDIF
  !
#ifdef __PARA
     CALL mp_bcast (iter,ionode_id,inter_pool_comm)
     CALL mp_bcast (itemp,ionode_id,inter_pool_comm)
     CALL mp_bcast (idope,ionode_id,inter_pool_comm)
     CALL mp_bcast (F_k_ful,ionode_id,inter_pool_comm)
#ENDIF
  !
  CALL mp_barrier (inter_pool_comm)
  !
END SUBROUTINE Fk_load



!-------------------------------------------------------------------------------
SUBROUTINE ft_save (iter)
!-------------------------------------------------------------------------------
#INCLUDE "f_defs.h"
  USE kinds,        ONLY : DP
  USE phcom,        ONLY : nmodes
  USE epwcom,       ONLY : nbndsub, dt
  USE elph2,        ONLY : ibndmin, ibndmax, nbnd_red, f_t_ful
  USE bte_var
#ifdef __PARA
  USE mp,           ONLY : mp_barrier, mp_bcast
  USE io_global,    ONLY : ionode_id, stdout
  USE mp_global,    ONLY : my_pool_id, inter_pool_comm
#ENDIF
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iter
  INTEGER             :: ik, ibnd
  !
  !
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     OPEN (99999,FILE='BTE/META/f_t_ful_tmp',FORM='unformatted',ACCESS='direct',RECL=DP,STATUS='replace')
     WRITE (99999,REC=1) DBLE(iter)
     WRITE (99999,REC=2) dt
     DO ik = 1, nk_ful_red
        DO ibnd = 1, nbnd_red
           WRITE (99999,REC=(ik-1)*nbnd_red+ibnd+2) f_t_ful(ibnd,ik)
        ENDDO
     ENDDO
     CLOSE (99999)
     !
     CALL SYSTEM ('cp -r BTE/META/f_t_ful_tmp BTE/META/f_t_ful')
     CALL SYSTEM ('rm BTE/META/f_t_ful_tmp')
     !
  ENDIF
  !
#ifdef __PARA
  CALL mp_barrier (inter_pool_comm)
#ENDIF
  !
  !
END SUBROUTINE ft_save



!-------------------------------------------------------------------------------
SUBROUTINE ft_load (iter)
!-------------------------------------------------------------------------------
#INCLUDE "f_defs.h"
  USE kinds,        ONLY : DP
  USE phcom,        ONLY : nmodes
  USE epwcom,       ONLY : nbndsub, dt
  USE elph2,        ONLY : ibndmin, ibndmax, nbnd_red, f_0_ful
  USE bte_var
  USE constants_epw, ONLY : au2fs
#ifdef __PARA
  USE mp,           ONLY : mp_bcast, mp_barrier
  USE io_global,    ONLY : ionode_id, stdout
  USE mp_global,    ONLY : my_pool_id, inter_pool_comm
#ENDIF
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(OUT) :: iter
  !
  INTEGER              :: ik, ibnd, cnt
  REAL(KIND=DP)        :: iter_tmp, dt_tmp
  !
  !
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     OPEN (99999,FILE='BTE/META/f_t_ful',FORM='unformatted',ACCESS='direct',RECL=DP,STATUS='old')
     !
     cnt = 0
     READ (99999,REC=1) iter_tmp
     READ (99999,REC=2) dt_tmp
     !
     IF (ABS(dt-dt_tmp) .GT. 1.0d-6) THEN
         WRITE (stdout,'(a,f12.6,a)') 'dt in last run : ', dt_tmp*au2fs, ' fs'
         CALL errore('ft_load','dt is not inconsistent with taht in last run',1)
     ENDIF
     !
     DO ik = 1, nk_ful_red
        DO ibnd = 1, nbnd_red
           READ (99999,REC=(ik-1)*nbnd_red+ibnd+2) f_0_ful(ibnd,ik) ! load as f_0_ful
           cnt = cnt + 1
        ENDDO
     ENDDO
     CLOSE (99999)
     !
     IF (cnt .NE. nk_ful_red*nbnd_red) &
        CALL errore('ft_load','length of metafile f_t_ful is incorrect, please delete BTE/META/f_t_ful manually',1)
     !
     iter = NINT(iter_tmp) + 1
     !
  ENDIF
  !
#ifdef __PARA
     CALL mp_bcast (iter,ionode_id,inter_pool_comm)
     CALL mp_bcast (f_0_ful,ionode_id,inter_pool_comm)
     CALL mp_barrier (inter_pool_comm)
#ENDIF
  !
  !
END SUBROUTINE ft_load



!----------------------------------------------------------------------------
SUBROUTINE edos_load ()
!----------------------------------------------------------------------------
#INCLUDE "f_defs.h"
  USE kinds,         ONLY : DP
  USE lsda_mod,      ONLY : nspin
  USE elph2,         ONLY : edos, ndos, cbnd_emin, vbnd_emax
  USE constants_epw, ONLY : ryd2ev
#ifdef __PARA
  USE io_files,      ONLY : tmp_dir, prefix
  USE io_global,     ONLY : stdout, ionode_id
  USE mp_global,     ONLY : my_pool_id
#ENDIF
  !
  IMPLICIT NONE
  !
  REAL(KIND=DP), ALLOCATABLE :: DOSint(:)
  REAL(KIND=DP)              :: ndos_tmp, nkf1_tmp, nkf2_tmp, nkf3_tmp
  INTEGER                    :: n
  CHARACTER (LEN=256)        :: fildos
  LOGICAL                    :: exist_fildos
  REAL(KIND=DP)              :: egap, egap_dos, egap_min, egap_max
  !
  !
  IF (ALLOCATED(edos)) DEALLOCATE(edos)
  !
  ! read eDOS from file
  fildos = TRIM(tmp_dir)//TRIM(prefix)//'.edos'
  INQUIRE (FILE=fildos,EXIST=exist_fildos)
  IF (.NOT. exist_fildos) CALL errore('edos_load','the edos file does not exist',1)
  !
  OPEN (43234,FILE=fildos,FORM='unformatted',ACCESS='direct',RECL=DP,STATUS='old')
  !
  READ (43234,REC=1) ndos_tmp
  READ (43234,REC=2) nkf1_tmp
  READ (43234,REC=3) nkf2_tmp
  READ (43234,REC=4) nkf3_tmp
  !
  ndos = NINT(ndos_tmp)
  ALLOCATE (edos(3,ndos))
  ALLOCATE (DOSint(ndos))
  edos = 0.0d0
  DOSint = 0.0d0
  !
  IF (nspin .EQ. 1 .OR. nspin .EQ. 4) THEN
     DO n = 1, ndos
        READ (43234,REC=(3*(n-1)+1)+4) edos(1,n) ! [eV]
        READ (43234,REC=(3*(n-1)+2)+4) edos(2,n) ! [1/eV/cm^3]
        READ (43234,REC=(3*(n-1)+3)+4) DOSint(n) ! [1/cm^3]        
     ENDDO
  ELSE
     DO n = 1, ndos
        READ (43234,REC=(3*(n-1)+1)+4) edos(1,n) ! [eV]
        READ (43234,REC=(3*(n-1)+2)+4) edos(2,n) ! [1/eV/cm^3]
        READ (43234,REC=(3*(n-1)+3)+4) edos(3,n) ! [1/eV/cm^3]
        READ (43234,REC=(3*(n-1)+4)+4) DOSint(n) ! [1/cm^3] 
     ENDDO       
  ENDIF
  !
  CLOSE (43234)
  !
  !
  ! check the bandgap
  egap_min = +9999999999
  egap_max = -9999999999
  !
  DO n = 1, ndos
     IF (edos(2,n) .EQ. 0.0d0) THEN
        !
        IF (edos(1,n) .LT. egap_min) egap_min = edos(1,n)
        IF (edos(1,n) .GT. egap_max) egap_max = edos(1,n)
        !
     ENDIF
  ENDDO 
  !
  egap_dos = egap_max-egap_min
  egap = (cbnd_emin-vbnd_emax)*ryd2ev
  !
  WRITE (stdout,'(/5x,a,i4,a,i4,a,i4,a)') 'Load eDOS from file (', NINT(nkf1_tmp), '*', NINT(nkf2_tmp), '*', NINT(nkf3_tmp), ' k-mesh)'
  WRITE (stdout,'(5x,a,f9.4,a)')          'Bandgap of eDOS    = ', egap_dos, ' eV'
  !
  IF (ABS(egap-egap_dos) .GT. 1.0d-6) WRITE (stdout,'(5x,a)') 'WARNING: Bandgap mismatch between eDOS file and present calculation'
  !
  !
  ! export readable edos file
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     OPEN (43234,FILE='BTE/EPCHECK/edos.dat',STATUS='replace')
     !
     WRITE(43234,'(4x,a)') 'ndos'
     WRITE(43234,'(4x,i6)') ndos
     WRITE(43234,'(4x,a)') 'k-mesh'
     WRITE(43234,'(4x,3i5)') NINT(nkf1_tmp), NINT(nkf2_tmp), NINT(nkf3_tmp)
     !
     IF (nspin .EQ. 1 .OR. nspin .EQ. 4) THEN
        WRITE(43234,'("    E (eV)                   dos (1/eV/cm^3)          Int dos (1/cm^3)")')
     ELSE
        WRITE(43234,'("    E (eV)                   dosup (1/eV/cm^3)        dosdw (1/eV/cm^3)        Int dos (1/Bohr^3)")')
     END IF
     !
     DO n = 1, ndos
        !
        IF (nspin .EQ. 1 .OR. nspin .EQ. 4) THEN 
           WRITE (43234,'(es25.15,2es25.15)') edos(1,n), edos(2,n), DOSint(n)       
        ELSE
           WRITE (43234,'(es25.15,3es25.15)') edos(1,n), edos(2,n), edos(3,n), DOSint(n)   
        ENDIF 
        !
     ENDDO 
     !
     CLOSE (43234)
     !
  ENDIF
  !
  !
  DEALLOCATE (DOSint)
  !
END SUBROUTINE edos_load
