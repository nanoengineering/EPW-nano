!-----------------------------------------------------------------------
SUBROUTINE selfen_phon (iq0)
  !-----------------------------------------------------------------------
  !
#INCLUDE "f_defs.h"
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE io_files,      ONLY : tmp_dir, prefix
  USE spin_orb,      ONLY : lspinorb
  USE phcom,         ONLY : nmodes
  USE epwcom,        ONLY : nbndsub, fsthick, eptemp, degaussw, save_m_mat, epthick, phdrag, &
                            bte, efermi_read, fermi_energy, smearing, nptype, epdope, nepdope, neptemp, nepdope
  USE pwcom,         ONLY : nelec
  USE elph2,         ONLY : ibndmin, ibndmax, nbnd_red, cbnd_emin, vbnd_emax, cfsthick, vfsthick, &
                            etf, epf17, wf_irr, wf_all, gammai_mode_all, ef_epw
  USE constants_epw, ONLY : ryd2mev, ryd2ev, two, zero, pi
  USE tetrahedron
  USE bte_var         
#ifdef __PARA
  USE mp_global,     ONLY : my_pool_id, npool
#ENDIF
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)        :: iq0
  !
  REAL(KIND=DP), PARAMETER   :: eps = 1.0d-6
  ! local variable
  INTEGER                    :: ik, iq, ikk, ikq, ibnd0, jbnd0, ibnd, jbnd, imode, itemp, idope, iq_red, nscat_k, ik_exp(nk_ful)
  REAL(KIND=DP)              :: g2, xk(3), ekk, ekq, wq, wgq, ef0, wgkk, wgkq, inv_temp
  REAL(KIND=DP)              :: gammai(nk_ful), gammai_phd(nk_ful), weight(nk_ful), gammai_exp(nk_ful)
  COMPLEX(KIND=DP)           :: epf, epf_tmp(nbnd_red,nbnd_red,nmodes)
  LOGICAL                    :: within_range
  ! external fuction
  REAL(KIND=DP), EXTERNAL    :: wgauss, w0gauss
  ! gamma file
  CHARACTER(LEN=256)         :: file_ufmt
  CHARACTER(LEN=12)          :: tnph
  CHARACTER(LEN=10)          :: q_num, ibnd_num, jbnd_num, imode_num, itemp_num, idope_num
  ! tetrahedron
  REAL(KIND=DP)              :: ekq_tetra(nk_ful), wsp
  REAL(KIND=DP), ALLOCATABLE :: tetra_c_ph(:,:)
  INTEGER                    :: iqtetra, itetra_(2,30,nk_ful)
  ! epf17 file
  CHARACTER(LEN=256)         :: epf17_ufmt
  CHARACTER(LEN=3)           :: epf17_cpu
  ! 
  !
  ! if bte=2, iq0 belongs to nq_irr, while if phdrag=T, iq0 belongs to
  ! nq_irr_red
  IF (phdrag) THEN
     iq_red = iq0
     iq = rirr2irr_q(iq0) ! transform to nq_irr space
  ELSE ! bte=2
     iq_red = iq0 ! in the same space, for gammai_mode_all
     iq = iq0
  ENDIF
  !
  ! epf17 filename
  IF (save_m_mat) THEN
     !
#ifdef __PARA
     CALL set_ndnmbr (0,my_pool_id+1,1,npool,epf17_cpu)
     epf17_ufmt = TRIM(tmp_dir)//TRIM(prefix)//'.epf17_'//epf17_cpu
#ENDIF
     OPEN (17517,FILE=epf17_ufmt,FORM='unformatted',ACCESS='direct',RECL=nbnd_red*nbnd_red*nmodes*2*DP,STATUS='old') 
     !
  ENDIF
  !
  ! spin factor 2
  !IF (nspin .EQ. 1 .OR. nspin .EQ. 4) THEN
  IF (.NOT. lspinorb) THEN
     wsp = 2.0d0
  ELSE
     wsp = 1.0d0
  ENDIF
  !
  !
  ! initialization of tetrahedron
  IF (smearing .EQ. 'tetra') THEN
     itetra_ = itetra
     ALLOCATE (tetra_c_ph(4,ntetra))
  ENDIF
  !
  !
  DO imode = 1, nmodes
     !
     IF (phdrag) THEN
        wq = wf_all(imode,iq_red)
     ELSE ! bte=2
        wq = wf_irr(imode,iq)
     ENDIF
     !
     DO ibnd0 = 1, nbnd_red
        !
        DO jbnd0 = 1, nbnd_red
           !
           ibnd = ibnd0+ibndmin-1
           jbnd = jbnd0+ibndmin-1
           !
           IF (smearing .EQ. 'tetra') THEN   
              ! tetrahedron
              ! linewidth of phonon
              DO ik = 1, nk_ful
                 ekq_tetra(ik) = etf(jbnd,2*ik) - etf(ibnd,2*ik-1) 
              ENDDO
              CALL eigen_tet(ntetra,ekq_tetra,tetra_i,tetra_w,nk_ful)
              CALL weight_tet(nk_ful,ntetra,wq,tetra_i,tetra_w,tetra_c_ph,wkt)
              !
           ENDIF
           !
           DO itemp = 1, neptemp
              !
              inv_temp = 1.0d0/eptemp(itemp)
              !
              IF (phdrag) THEN
                 wgq = wgauss(-wq*inv_temp,-99)
                 wgq = wgq/(1.0d0-2.0d0*wgq)
              ENDIF
              !
              DO idope = 1, nepdope
                 !
                 IF (efermi_read) THEN
                    ef0 = fermi_energy
                 ELSE
                    ef0 = ef_epw(itemp,idope)
                 ENDIF
                 !
                 weight = 0.0d0
                 gammai = 0.0d0
                 gammai_phd = 0.0d0
                 !
                 DO ik = 1, nk_ful
                    !
                    ikq = 2 * ik
                    ikk = ikq - 1
                    !
                    xk(:) = xkf_ful(:,ik)
                    !
                    within_range = .FALSE.
                    !
                    IF ( (MAXVAL(ABS(xk(:))) .GT. eps) .AND. &
                         (etf(ibnd,ikk) .GE. vbnd_emax-vfsthick .AND. etf(ibnd,ikk) .LE. cbnd_emin+cfsthick) .AND. &
                         (etf(jbnd,ikq) .GE. vbnd_emax-vfsthick-epthick .AND. etf(jbnd,ikq) .LE. cbnd_emin+cfsthick+epthick) ) within_range = .TRUE.
                    !
                    IF (within_range) THEN
                       !     
                       ekk = etf(ibnd,ikk) - ef0
                       wgkk = wgauss(-ekk*inv_temp,-99) 
                       ! 
                       ekq = etf(jbnd,ikq) - ef0
                       wgkq = wgauss(-ekq*inv_temp,-99) 
                       !
                       IF (save_m_mat) THEN
                          READ (17517,REC=ik) epf_tmp(1:nbnd_red,1:nbnd_red,1:nmodes)
                          epf = epf_tmp(jbnd0,ibnd0,imode)
                       ELSE
                          ! ERROR, not consider temperature and carrier concentration
                          epf = epf17(ik,jbnd0,ibnd0,imode,1,1)
                       ENDIF
                       !
                       g2 = (ABS(epf)**2.0d0)/(2.0d0*wq)
                       !
                       ! tetrahedral smearing
                       IF (smearing .EQ. 'tetra') THEN
                          !
                          iqtetra = 1
                          DO WHILE (itetra_(1,iqtetra,ik) .NE. 0)
                             !
                             weight(ik) = weight(ik) + pi * g2 * tetra_c_ph(itetra_(2,iqtetra,ik),itetra_(1,iqtetra,ik))
                             !
                             iqtetra = iqtetra + 1
                             IF (iqtetra .GT. 30) CALL errore('selfen_phon','too many tetrahedron associated with one q-point',1)
                             !
                          ENDDO  
                          !
                          gammai(ik) = wsp * (wgkk-wgkq) * weight(ik)
                          IF (phdrag) gammai_phd(ik) = wsp * (wgkk*(1.0d0-wgkq)*wgq) * weight(ik)
                          !
                       ! lorentzian smearing
                       ELSEIF (smearing .EQ. 'lortz') THEN
                          !
                       ELSE
                          !
                          CALL errore('selfen_phon','wrong smearing flag',1)
                          !
                       ENDIF
                       !  
                    ENDIF ! within_range
                    !
                 ENDDO ! ik
                 !
                 gammai_mode_all(itemp,idope,imode,iq_red) = gammai_mode_all(itemp,idope,imode,iq_red) + SUM(gammai(:))
                 !
                 !
                 ! output binary file
                 IF (phdrag) THEN
                    !
                    nscat_k = 0
                    DO ik = 1, nk_ful
                       IF (gammai_phd(ik) .NE. 0.0d0) THEN
                          nscat_k = nscat_k + 1
                          ik_exp(nscat_k) = ik
                          gammai_exp(nscat_k) = gammai_phd(ik)
                       ENDIF
                    ENDDO
                    !
                    IF (nscat_k .EQ. 0) THEN
                       nscat_k = 1
                       ik_exp(nscat_k) = 1
                       gammai_exp(nscat_k) = 0.0d0
                    ENDIF
                    !
                    WRITE(q_num,'(i10)') iq_red
                    WRITE(ibnd_num,'(i10)') ibnd0
                    WRITE(jbnd_num,'(i10)') jbnd0
                    WRITE(imode_num,'(i10)') imode
                    WRITE(itemp_num,'(i10)') itemp
                    WRITE(idope_num,'(i10)') idope
                    !
                    tnph = 'T'//TRIM(ADJUSTL(itemp_num))//'_N'//TRIM(ADJUSTL(idope_num))//'_ph'//TRIM(ADJUSTL(imode_num))
                    !
                    file_ufmt = 'BTE/GAMMAI/'//TRIM(ADJUSTL(tnph))//'/gammai_'//TRIM(ADJUSTL(q_num))//'_'//TRIM(ADJUSTL(ibnd_num)) &
                                                                                                    //'_'//TRIM(ADJUSTL(jbnd_num))
                    OPEN (80300,FILE=file_ufmt,FORM='unformatted',ACCESS='direct',RECL=4+nscat_k*(4+DP),STATUS='replace')
                    WRITE (80300,REC=1) nscat_k, ik_exp(1:nscat_k), gammai_exp(1:nscat_k)
                    CLOSE (80300)
                    !
                 ENDIF
                 !
              ENDDO ! idope
              !
           ENDDO ! itemp
           ! 
        ENDDO ! jbnd
        !
     ENDDO ! ibnd
     !
  ENDDO ! imode
  !
  !
  IF (save_m_mat) CLOSE (17517)
  !
  IF (ALLOCATED(tetra_c_ph)) DEALLOCATE (tetra_c_ph)
  !
  RETURN
  !
  !
END SUBROUTINE selfen_phon
