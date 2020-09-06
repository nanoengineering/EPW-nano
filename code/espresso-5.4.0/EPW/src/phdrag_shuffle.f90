!---------------------------------------------------------------------------------
SUBROUTINE phdrag_shuffle (nqc, xqc)
!---------------------------------------------------------------------------------
#INCLUDE "f_defs.h"
  USE kinds,             ONLY : DP
  USE pwcom,            ONLY : nbnd, nks, nkstot, isk, et, xk, alat, nelec
  USE ions_base,         ONLY : nat
  USE phcom,             ONLY : nmodes
  USE cell_base,         ONLY : at, bg
  USE epwcom,            ONLY : nbndsub, lpolar, nbndskip, vme, eig_read, & 
                                nkf1, nkf2, nkf3, nqf1, nqf2, nqf3, bte, epdim, &
                                edos_read, asr_eph, smearing, ph_read, ephl_read, &
                                epdope, nepdope, save_m_mat, save_m_matw, shengbte_read, &
                                save_m_ph, save_t_el, epthick, phdrag, neptemp, nepdope, eptemp
  USE elph2,             ONLY : uf_ful, wf_ful, vph_ful, wf_all, uf_all, vph_all, wf_irr, uf_irr, vph_irr, &
                                nrr_k, nrr_q, cu, cuq, irvec, ndegen_k, ndegen_q, &
                                wslen, chw, chw_ks, cvmew, cdmew, rdw, epmatwp, epmatwe, epmatq, &
                                wf, etf, etf_k, etf_ks, xqf, xkf, wkf, wqf, &
                                dynq, nqtotf, nkqf, epf17, nkf, nqf, et_ks, &
                                ibndmin, ibndmax, lambda_all, dmec, dmef, vmef, &
                                gamma_all, nkqtotf, epsi, zstar, efnew, &
                                ! THL
                                nbnd_red, cbnd_emin, vbnd_emax, ef_m, delta_egap, outside_gap, cfsthick, vfsthick, vbnd_num, &
                                sigmai_mode_all_abs, sigmai_mode_all_emi, &
                                etf_ful, vel_ful, epmatwp_asr, eph_vogl, gammai_mode_all
  USE constants_epw, ONLY : ryd2ev, ryd2mev, rydcm1, one, two, twopi, pi, au2cm, au2ps, ryd2thz, kB
  USE bte_var
  USE bte_func
  USE tetrahedron     
#ifdef __PARA
  USE io_files,      ONLY : tmp_dir, prefix
  USE mp,            ONLY : mp_barrier, mp_bcast, mp_sum
  USE io_global,     ONLY : ionode_id, stdout, ionode
  USE mp_global,     ONLY : my_pool_id, nproc_pool, intra_image_comm, &
                            inter_pool_comm, me_pool, root_pool, intra_pool_comm, npool
  USE mp_world,      ONLY : mpime
#ENDIF
  !
  IMPLICIT NONE
  !
  complex(kind=DP), ALLOCATABLE :: &
    epmatwe_q (:,:,:,:),     &
    epmatwef (:,:,:,:)           ! e-p matrix  in el wannier - fine Bloch phonon grid
  complex(kind=DP), ALLOCATABLE :: &
    epmatf( :, :, :),           &! e-p matrix  in smooth Bloch basis, fine mesh
    cufkk ( :, :),              &! Rotation matrix, fine mesh, points k
    cufkq ( :, :),              &! the same, for points k+q
    uf    ( :, :),              &! Rotation matrix for phonons
    bmatf ( :, :)                ! overlap U_k+q U_k^\dagger in smooth Bloch basis, fine mesh
  integer, intent(in) :: &
    nqc                          ! number of qpoints in the coarse grid
  real(kind=DP), intent(in) :: &
    xqc (3, nqc)                 ! qpoint list, coarse mesh
  real(kind=DP), parameter :: eps = 0.01/ryd2mev
  !
  !
  integer :: ir, na, fermicount, nrec, indnew, indold, lrepmatw, ios, irq, lrepmatw_asr
  LOGICAL :: exst, file_exist
  real(kind=DP) :: xxq(3), xxk(3), xkk(3), xkq(3), xq_fbz(3), length
  ! 
  ! THL
  ! parallelization
  INTEGER                       :: ik_star, ik_stop, iq_star, iq_stop
  ! local variable
  INTEGER                       :: ik, ik0, ik_red, ik0_red, iq, iq_irr, iq_ful, iq_irr_red, iq0, io, itemp
  INTEGER                       :: ibnd, jbnd, ibnd0, jbnd0
  INTEGER                       :: imode, mu, nu, i, j, k
  INTEGER                       :: ikk, ikq, iqq
  REAL(KIND=DP)                 :: etf_ks_tmp(nbndsub), wmax_mode(nmodes), vel_tmp(3,nbndsub), vph_tmp(3,nmodes)
  COMPLEX(KIND=DP)              :: dmef_tmp(3,nbndsub,nbndsub)
  ! group velocity
  REAL(KIND=DP)                 :: vkq(3,nbndsub)
  ! polar interaction
  COMPLEX(KIND=DP)              :: uf_cart(nmodes,nmodes), ephl_vogl(nmodes)
  CHARACTER (LEN=256)           :: vogl_ufmt
  CHARACTER (LEN=3)             :: vogl_cpu
  CHARACTER (LEN=256)           :: filename_check
  INTEGER                       :: iq_vogl, vogl_num(npool)
  ! check within_range
  LOGICAL                       :: within_range
  ! save memory
  CHARACTER (LEN=256)           :: epf17_ufmt 
  CHARACTER (LEN=3)             :: epf17_cpu
  REAL(KIND=DP)                 :: wf_dyn(nmodes), vph_dyn(3,nmodes)
  COMPLEX(KIND=DP)              :: uf_dyn(nmodes,nmodes)
  CHARACTER (LEN=256)           :: dyn_ufmt, wf_ufmt, vph_ufmt, uf_ufmt
  CHARACTER (LEN=3)             :: dyn_cpu
  INTEGER                       :: iq_dyn, dyn_num(npool)
  ! save time
  CHARACTER (LEN=256)           :: ham_ufmt, etf_ufmt, cuf_ufmt
  CHARACTER (LEN=3)             :: ham_cpu
  INTEGER                       :: icpu, ik_ham, ham_num(npool)
  INTEGER                       :: ijk_k(3), ijk_q(3), ijk_kq(3), ikq_ham
  REAL(KIND=DP)                 :: etf_ham(nbndsub), etf_ks_ham(nbndsub)
  COMPLEX(KIND=DP)              :: cufkk_ham(nbndsub,nbndsub)
  ! date and time
  CHARACTER(LEN=8)              :: date_
  CHARACTER(LEN=10)             :: time_
  CHARACTER(LEN=5)              :: zone_
  INTEGER                       :: values_(8)
  ! time
  REAL(KIND=DP)                 :: t0=0.0d0, t0i, t0f, t1=0.0d0, t1i, t1f, t2=0.0d0, t2i, t2f, &
                                   t3=0.0d0, t3i, t3f, t4=0.0d0, t4i, t4f, t5=0.0d0, t5i, t5f, &
                                   t6=0.0d0, t6i, t6f, t7=0.0d0, t7i, t7f, t8=0.0d0, t8i, t8f, &
                                   t9=0.0d0, t9i, t9f, ta=0.0d0, tai, taf, tb=0.0d0, tbi, tbf, &
                                   tc=0.0d0, tci, tcf
  !
  COMPLEX(KIND=DP), ALLOCATABLE :: eph_vogl_tmp(:,:)
  ! For ShengBTE phv, Qian, Nov 2017
  INTEGER,ALLOCATABLE :: SBTE_q_ful(:),SBTE_q_1irr(:)
  REAL(kind=DP),ALLOCATABLE :: SBTE_q_crys_ful(:,:),S_fre_irr(:,:,:),S_sct_irr(:,:,:)
  CHARACTER(LEN=256)    :: S_sct_itemp_fmt, string_temp
  INTEGER :: iq_red, f2f, S_temp, nlines
  !
  WRITE (stdout,'(//5x,a)') '==================================================================='
  WRITE (stdout,'(5x,a)')   '                        Select out q points                        '
  WRITE (stdout,'(5x,a/)')  '==================================================================='
  !
  ! mapping ful -> irr for select q point
  WRITE (stdout,'(/5x,a)') 'Mapping wf, vph and uf from ful-BZ to irr-BZ'
  IF (ALLOCATED(wf_irr)) DEALLOCATE(wf_irr)
  IF (ALLOCATED(vph_irr)) DEALLOCATE(vph_irr)
  IF (ALLOCATED(uf_irr)) DEALLOCATE(uf_irr)
  !
  ALLOCATE (wf_irr(nmodes,nq_irr))
  ALLOCATE (vph_irr(3,nmodes,nq_irr))
  ALLOCATE (uf_irr(nmodes,nmodes,nq_irr))
  wf_irr = 0.0d0
  vph_irr = 0.0d0   
  uf_irr = (0.0d0,0.0d0)
  !
  DO iq = 1, nq_irr
     wf_irr(:,iq) = wf_ful(:,1,1,irr2ful_q(iq))
     vph_irr(:,:,iq) = vph_ful(:,:,irr2ful_q(iq))
     uf_irr(:,:,iq) = uf_ful(:,:,1,1,irr2ful_q(iq))
  ENDDO
  !
  WRITE (stdout,'(//5x,a)') '==================================================================='
  WRITE (stdout,'(5x,a)')   '                      Phonon-drag calculation                      '
  WRITE (stdout,'(5x,a/)')  '==================================================================='
  !
  IF (shengbte_read) THEN
   CALL read_shengbte()
   CALL phph_export()
  !
   filename_check = TRIM(prefix)//'.epcheck'
   INQUIRE (FILE=filename_check,EXIST=file_exist)
   IF (file_exist) CALL phcheck_scat (filename_check)
  ENDIF
  !
  DEALLOCATE (wf_ful)
  DEALLOCATE (vph_ful)
  DEALLOCATE (uf_ful)
  !
  !
  CALL reduce_index_ph (wf_irr)
  !
  !
  WRITE (stdout,'(/5x,a)') 'Mapping wf and vph from irr-BZ to red-irr-BZ'
  ALLOCATE (wf_all(nmodes,nq_irr_red))
  ALLOCATE (vph_all(3,nmodes,nq_irr_red))
  wf_all = 0.0d0
  vph_all = 0.0d0   
  !
  DO iq = 1, nq_irr_red
     wf_all(:,iq) = wf_irr(:,rirr2irr_q(iq))
     vph_all(:,:,iq) = vph_irr(:,:,rirr2irr_q(iq))
  ENDDO
  !
  DEALLOCATE (wf_irr)
  DEALLOCATE (vph_irr)
  !
  !
  IF (smearing .EQ. 'tetra') THEN
     !
     DEALLOCATE (tetra_i, tetra_c, tetra_w, wkt, itetra)
     !
     ntetra = 6*nkf1*nkf2*nkf3
     ALLOCATE (tetra_i(4,ntetra),tetra_c(4,ntetra),tetra_w(4,ntetra),wkt(nkf1*nkf2*nkf3),itetra(2,30,nkf1*nkf2*nkf3))
     CALL make_kp_reg (nkf1, nkf2, nkf3, wkt, ntetra, tetra_i, itetra)
     !
  ENDIF
  !
  !
  WRITE (stdout,'(//5x,a)') '==================================================================='
  WRITE (stdout,'(5x,a)')   '                Inverse Wannier interpolation (ph-e)               '
  WRITE (stdout,'(5x,a/)')  '==================================================================='
  !
  IF (save_m_mat) THEN
     WRITE (stdout,'(/5x,a/)')  'Save following variables to disk to reduce the using of internal memory :'
     WRITE (stdout,'(13x,a,f8.1,a)')  'epf17   ~', (DBLE(nk_ful)*1.0d-6)*DBLE(nbnd_red*nbnd_red*nmodes*2*DP), ' MB' 
  ENDIF
  !
  IF (save_m_matw) THEN
     WRITE (stdout,'(/5x,a)') 'epmatwe and epmatwp have been prepared in e-ph calculation'
  ENDIF
  !
  !
  IF (save_t_el) THEN
     WRITE (stdout,'(/5x,a)') 'Electron properties have been prepared in e-ph calculation'
  ELSE
     WRITE (stdout,'(/5x,a)') 'Electron part will be export automatically when phdrag = true'
     WRITE (stdout,'(5x,a/)') 'Read etf(k+q), etf_ks(k+q) and cufkk(k+q) from saved files to reduce the time consumption'
     WRITE (stdout,'(13x,a,f8.1,a)')  'etf     ~', (DBLE(nk_ful)*1.0d-6)*DBLE(nbndsub*DP), ' MB' 
     WRITE (stdout,'(13x,a,f8.1,a)')  'etf_ks  ~', (DBLE(nk_ful)*1.0d-6)*DBLE(nbndsub*DP), ' MB' 
     WRITE (stdout,'(13x,a,f8.1,a)')  'cufkk   ~', (DBLE(nk_ful)*1.0d-6)*DBLE(nbndsub*nbndsub*2*DP), ' MB' 
  ENDIF
  !
  !
  !
  IF (lpolar) THEN
     !
     WRITE (6,'(/5x,a/)') 'Dielectric tensor:'
     WRITE (6,'(13x,3f12.6)') ((epsi(i,j),j=1,3),i=1,3)
     !
     WRITE (6,'(/5x,a)') 'Born effective charges:'
     DO na = 1, nat
        WRITE (6,'(/13x,a,i3/)') 'Atom #', na
        WRITE (6,'(13x,3f12.6)') ((zstar(i,j,na),j=1,3),i=1,3)
     END DO
     !
     ! copy eph_vogl_ful -> eph_vogl_irr
     WRITE (stdout,'(/5x,a)') 'Mapping eph_vogl from full to irr BZ'
     ALLOCATE(eph_vogl_tmp(nmodes,nq_irr))
     DO iq = 1, nq_irr
        eph_vogl_tmp(:,iq) = eph_vogl(:,1,1,irr2ful_q(iq))
     ENDDO
     DEALLOCATE(eph_vogl)
     ALLOCATE(eph_vogl(nmodes,1,1,nq_irr))
     do iq = 1, nq_irr
        eph_vogl(:,1,1,iq) = eph_vogl_tmp(:,iq)
     enddo
     DEALLOCATE(eph_vogl_tmp)
     !
  ENDIF
  !
  !
  !
  CALL mp_barrier (inter_pool_comm)
     !
     !
  IF (.NOT. save_t_el) THEN
     !
     WRITE (stdout,'(/5x,a,f8.2,a)') 'Prepare electron properties :'
     WRITE (stdout,'(13x,a,f8.1,a)')  'etf     ~', (DBLE(nk_ful)*1.0d-6)*DBLE(nbndsub*DP), ' MB' 
     WRITE (stdout,'(13x,a,f8.1,a)')  'etf_ks  ~', (DBLE(nk_ful)*1.0d-6)*DBLE(nbndsub*DP), ' MB' 
     WRITE (stdout,'(13x,a,f8.1,a)')  'cufkk   ~', (DBLE(nk_ful)*1.0d-6)*DBLE(nbndsub*nbndsub*2*DP), ' MB' 
     !
     CALL CPU_TIME (tci)
     !
     CALL para_bounds (ik_star, ik_stop, nk_ful)
     !
     ham_num = 0
     ham_num(my_pool_id+1) = ik_stop-ik_star+1
     !
#ifdef __PARA
     CALL mp_barrier (inter_pool_comm)
     CALL mp_sum (ham_num,inter_pool_comm)
#ENDIF
     !
     !
#ifdef __PARA
     CALL set_ndnmbr (0,my_pool_id+1,1,npool,ham_cpu)
     ham_ufmt = TRIM(tmp_dir)//TRIM(prefix)//'.ham_'//ham_cpu
#ENDIF
     OPEN (40204,FILE=ham_ufmt,FORM='unformatted',ACCESS='direct',RECL=2*nbndsub*DP+2*nbndsub*nbndsub*DP,STATUS='replace') 
     !
     OPEN (1294,FILE='BTE/META/xkf_ful_cryst',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='old') 
     !
     etf_ham    = 0.0d0
     etf_ks_ham = 0.0d0
     cufkk_ham  = (0.0d0,0.0d0)
     !
     DO ik = ik_star, ik_stop
        !
        READ (1294,REC=ik) xkk(1:3)
        !
        IF (eig_read) THEN
           CALL hamwan2bloch (nbndsub, nrr_k, irvec, ndegen_k, xkk, cufkk_ham, etf_ks_ham, chw_ks)     
        ENDIF
        !
        CALL hamwan2bloch (nbndsub, nrr_k, irvec, ndegen_k, xkk, cufkk_ham, etf_ham, chw)        
        ! 
        DO ibnd = vbnd_num+1, nbndsub
           etf_ham(ibnd) = etf_ham(ibnd) + delta_egap
        ENDDO 
        !
        WRITE (40204,REC=ik-ik_star+1) etf_ham(1:nbndsub), etf_ks_ham(1:nbndsub), cufkk_ham(1:nbndsub,1:nbndsub)
        !
     ENDDO
     !
     CLOSE (1294)
     CLOSE (40204)
     !
     !
#ifdef __PARA
     CALL mp_barrier (inter_pool_comm)
#ENDIF
     !
     !
     IF (my_pool_id .EQ. ionode_id) THEN
        ! 
        etf_ufmt  = TRIM(tmp_dir)//TRIM(prefix)//'.ele'
        cuf_ufmt  = TRIM(tmp_dir)//TRIM(prefix)//'.elc'
        OPEN (99901,FILE=etf_ufmt,FORM='unformatted',ACCESS='direct',RECL=2*nbndsub*DP,STATUS='replace')
        OPEN (99902,FILE=cuf_ufmt,FORM='unformatted',ACCESS='direct',RECL=2*nbndsub*nbndsub*DP,STATUS='replace')
        !
        ik_ham = 1
        !
        DO icpu = 1, npool
           !
           WRITE(ham_cpu,'(i3)') icpu
           ham_ufmt = TRIM(tmp_dir)//TRIM(prefix)//'.ham_'//TRIM(ADJUSTL(ham_cpu))
           OPEN (40204,FILE=ham_ufmt,FORM='unformatted',ACCESS='direct',RECL=2*nbndsub*DP+2*nbndsub*nbndsub*DP,STATUS='old') 
           !
           DO ik = 1, ham_num(icpu)
              !
              READ (40204,REC=ik) etf_ham(1:nbndsub), etf_ks_ham(1:nbndsub), cufkk_ham(1:nbndsub,1:nbndsub)
              WRITE (99901,REC=ik_ham) etf_ham(1:nbndsub), etf_ks_ham(1:nbndsub)
              WRITE (99902,REC=ik_ham) cufkk_ham(1:nbndsub,1:nbndsub)       
              ik_ham = ik_ham + 1
              !
           ENDDO
           !
           CLOSE (40204)
           !
        ENDDO
        !
        CALL SYSTEM ('rm '//TRIM(tmp_dir)//TRIM(prefix)//'.ham_*')
        CLOSE (99901)
        CLOSE (99902)
        !
     ENDIF
     !
     mkq1 = nkf1/nqf1 ! must be 1 in phonon-drag calculation
     mkq2 = nkf2/nqf2
     mkq3 = nkf3/nqf3
     !
     CALL CPU_TIME (tcf)
     tc = tcf-tci
     WRITE (stdout,'(/13x,a,f8.2,a)') 'Out k | Out q | ham     W->B :', tc, 's'
     !
  ENDIF
  !
  !
#ifdef __PARA
  CALL mp_barrier (inter_pool_comm)
#ENDIF
  !
  !
  ! globle variable
  !
  IF (ALLOCATED(cufkk)) DEALLOCATE (cufkk)
   ALLOCATE (cufkk(nbndsub, nbndsub))
  IF (ALLOCATED(cufkq)) DEALLOCATE (cufkq)
   ALLOCATE (cufkq(nbndsub, nbndsub))
  IF (ALLOCATED(epmatf)) DEALLOCATE (epmatf)
   ALLOCATE (epmatf(nbndsub, nbndsub, nmodes))
  IF (ALLOCATED(bmatf)) DEALLOCATE (bmatf)
   ALLOCATE (bmatf(nbndsub, nbndsub))
  IF (ALLOCATED(etf)) DEALLOCATE (etf)
   ALLOCATE (etf(nbndsub, 2*nk_ful))
  IF (ALLOCATED(etf_ks)) DEALLOCATE (etf_ks)
   ALLOCATE (etf_ks(nbndsub, 2*nk_ful)) 
  IF (ALLOCATED(epmatwef)) DEALLOCATE (epmatwef)
   ALLOCATE (epmatwef(nbndsub, nbndsub, nrr_k, nmodes))
  IF (ALLOCATED(epf17)) DEALLOCATE (epf17)
   IF (.NOT. save_m_mat) ALLOCATE (epf17 (nk_ful, nbnd_red, nbnd_red, nmodes, 1, 1)) 
  !
  !
  !
  ! ======================================================
  !                         q loop
  ! ======================================================
  ! parallelization for nq_irr_red
  CALL para_bounds (iq_star, iq_stop, nq_irr_red)
  !
  ! BTE and selfen variable
  ALLOCATE (gammai_mode_all(neptemp,nepdope,nmodes,nq_irr_red)) ! here is nq_irr_red, different from bte=2
  gammai_mode_all = 0.0d0  
  !
  !
  DO iq_red = iq_star, iq_stop
     !
     CALL CPU_TIME (t0i)
     !
     CALL start_clock ( 'ep-interp' )
     !
     iq = rirr2irr_q(iq_red)
     !
     xxq(:) = xqf_irr(:,iq)     ! xqf_irr is assumed to be in crys coord
     !
     !
     ! --------------------------------------------------------------
     ! epmat : Wannier el and Wannier ph -> Wannier el and Bloch ph
     ! --------------------------------------------------------------
     !
     CALL CPU_TIME (t5i)
     !
     CALL ephwan2blochp (nmodes, xxq, irvec, ndegen_q, nrr_q, uf_irr(:,:,iq), epmatwef, nbndsub, nrr_k)
     !
     CALL CPU_TIME (t5f)
     t5 = t5f - t5i
     !
     !
     !
     t6 = 0.0d0
     t7 = 0.0d0
     t8 = 0.0d0
     t9 = 0.0d0
     !
     !
  IF (save_m_mat) THEN
     !
#ifdef __PARA
     CALL set_ndnmbr (0,my_pool_id+1,1,npool,epf17_cpu)
     epf17_ufmt = TRIM(tmp_dir)//TRIM(prefix)//'.epf17_'//epf17_cpu
#ENDIF
     OPEN (17017,FILE=epf17_ufmt,FORM='unformatted',ACCESS='direct',RECL=nbnd_red*nbnd_red*nmodes*2*DP,STATUS='replace') 
     !
  ELSE
     epf17 = (0.0d0,0.0d0)
     !
  ENDIF
     !
     !
     OPEN (12053,FILE='BTE/META/xkf_ful_cryst',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='old') 
     !
     etf_ufmt  = TRIM(tmp_dir)//TRIM(prefix)//'.ele'
     cuf_ufmt  = TRIM(tmp_dir)//TRIM(prefix)//'.elc'
     OPEN (40204,FILE=etf_ufmt,FORM='unformatted',ACCESS='direct',RECL=2*nbndsub*DP,STATUS='old')
     OPEN (40504,FILE=cuf_ufmt,FORM='unformatted',ACCESS='direct',RECL=2*nbndsub*nbndsub*DP,STATUS='old')
     ! 
     !
     DO ik = 1, nk_ful
        !
        ikq = 2 * ik
        ikk = ikq - 1  
        ! 
        READ (12053,REC=ik) xkk(1:3) ! xkf_ful is assumed to be in crys coord
        xkq = xkk + xxq
        !
        READ (40204,REC=ik) etf(1:nbndsub,ikk), etf_ks(1:nbndsub,ikk)
        !
        !
        ! ------------------------------------------------------        
        ! hamiltonian : Wannier -> Bloch 
        ! ------------------------------------------------------
        !
        CALL CPU_TIME (t6i)
        !
        IF (.NOT. save_t_el) THEN
           !
           CALL hamwan2bloch (nbndsub, nrr_k, irvec, ndegen_k, xkq, cufkq, etf(:,ikq), chw)
           !
           DO ibnd = vbnd_num+1, nbndsub
              etf(ibnd,ikq) = etf(ibnd,ikq) + delta_egap
           ENDDO
           !
        ELSE
           !
           ! transform k and q to crystal index
           ijk_k(1) = NINT(xkk(1)*DBLE(nkf1))
           ijk_k(2) = NINT(xkk(2)*DBLE(nkf2))
           ijk_k(3) = NINT(xkk(3)*DBLE(nkf3))
           ijk_q(1) = NINT(xxq(1)*DBLE(nqf1))
           ijk_q(2) = NINT(xxq(2)*DBLE(nqf2))
           ijk_q(3) = NINT(xxq(3)*DBLE(nqf3))
           !
           ijk_kq = ijk_fbz(ijk_k,ijk_q)
           ikq_ham = ijk2id(ijk_kq)
           !
           READ (40204,REC=ikq_ham) etf(1:nbndsub,ikq), etf_ks(1:nbndsub,ikq)
           READ (40504,REC=ikq_ham) cufkq(1:nbndsub,1:nbndsub)
           !
        ENDIF
        !
        CALL CPU_TIME (t6f)
        t6 = t6 + (t6f-t6i)
        !
        !
        ! check the scattering event is within_range
        within_range = .FALSE.
        !
        DO ibnd = ibndmin, ibndmax
           DO jbnd = ibndmin, ibndmax
              !
              IF ( (etf(ibnd,ikk) .GE. vbnd_emax-vfsthick .AND. etf(ibnd,ikk) .LE. cbnd_emin+cfsthick) .AND. &
                   (etf(jbnd,ikq) .GE. vbnd_emax-vfsthick-epthick .AND. etf(jbnd,ikq) .LE. cbnd_emin+cfsthick+epthick) ) within_range = .TRUE.
              !
           ENDDO
        ENDDO    
        ! 
        IF (within_range) THEN
           !
           ! --------------------------------------------------------------
           ! epmat : Wannier el and Bloch ph -> Bloch el and Bloch ph
           ! --------------------------------------------------------------
           !
           CALL CPU_TIME (t9i)
           !
           READ (40504,REC=ik) cufkk(1:nbndsub,1:nbndsub)
           !
           CALL ephwan2bloch (nbndsub, nrr_k, irvec, ndegen_k, epmatwef, xkk, cufkk, cufkq, epmatf, nmodes)
           !
           IF (lpolar) THEN
              !
              CALL compute_bmn_para2 (nbndsub, nkstot, cufkk, cufkq, bmatf)
              !
              IF ( (ABS(xxq(1)) .GT. eps) .OR. (ABS(xxq(2)) .GT. eps) .OR. (ABS(xxq(3)) .GT. eps) ) THEN
                 !
                 CALL cryst_to_cart (1, xxq, bg, 1)
                 !
                 DO ibnd = 1, nbndsub
                    DO jbnd = 1, nbndsub
                       !
                       CALL rgd_blk_epw3(uf_irr(:,:,iq), epmatf(ibnd,jbnd,:), 1, 1, eph_vogl(:,1,1,iq), nmodes, bmatf(ibnd,jbnd), +1.d0)
                       !
                    ENDDO
                 ENDDO
                 !
                 CALL cryst_to_cart (1, xxq, at, -1)
                 !
              ENDIF
              !
           ENDIF
           !
           CALL CPU_TIME (t9f)
           t9 = t9 + (t9f-t9i)
           !
           ! write epmatf to file / store in memory
           IF (save_m_mat) THEN
              !
              WRITE (17017,REC=ik) epmatf(ibndmin:ibndmax,ibndmin:ibndmax,1:nmodes)
              !
           ELSE
              !
              DO jbnd = ibndmin, ibndmax
                 DO ibnd = ibndmin, ibndmax
                    DO imode = 1, nmodes
                       !
                       ibnd0 = ibnd-ibndmin+1
                       jbnd0 = jbnd-ibndmin+1
                       epf17(ik,jbnd0,ibnd0,imode,1,1) = epmatf(jbnd,ibnd,imode)
                       !
                    ENDDO
                 ENDDO
              ENDDO
              !
           ENDIF
           !
        ENDIF ! within_range
        !
     ENDDO ! k loop
     !
     CLOSE (12053)
     !
     !
     ! ======================================================
     ! phonon self-energy
     ! ======================================================
     CALL CPU_TIME (tai)
     !
     CALL selfen_phon (iq_red)
     !
     CALL CPU_TIME (taf)
     ta = taf - tai
     !
     CALL CPU_TIME (t0f)
     t0 = t0f - t0i
     !
     CALL stop_clock ( 'ep-interp' )
     !
     !
     CALL DATE_AND_TIME (date_,time_,zone_,values_)
     WRITE (stdout,'(/5x,a,3f10.6,a,i4,a,i4,a)') 'k = (', xxq(1:3), '), No. (', iq_red-iq_star+1, '/', iq_stop-iq_star+1, ')'
     !WRITE (stdout,'(/13x,a,f8.2,a)') 'In  q | Out k | dyn     W->B :', t2, 's'
     !WRITE (stdout,'(13x,a,f8.2,a)')  'In  q | Out k | dme     W->B :', t3, 's'
     !WRITE (stdout,'(13x,a,f8.2,a)')  'In  q | Out k | vme     W->B :', t4, 's'
     WRITE (stdout,'(/13x,a,f8.2,a)')  'In  q | Out k | eph(ph) W->B :', t5, 's'
     WRITE (stdout,'(13x,a,f8.2,a)')  'In  q | In  k | ham     W->B :', t6, 's'
     !WRITE (stdout,'(13x,a,f8.2,a)')  'In  q | In  k | dme     W->B :', t7, 's'
     !WRITE (stdout,'(13x,a,f8.2,a)')  'In  q | In  k | vme     W->B :', t8, 's'
     WRITE (stdout,'(13x,a,f8.2,a)')  'In  q | In  k | eph(el) W->B :', t9, 's'
     WRITE (stdout,'(13x,a,f8.2,a)')  'In  q | Out k | self-energy  :', ta, 's'
     WRITE (stdout,'(13x,a,f8.2,a,i2,a,i2,a,i2)')  '                TOTAL        :', t0, 's | ', values_(5), ':', values_(6), ':', values_(7)
     !
     !
     IF (save_m_mat) CLOSE (17017)
     CLOSE (40204)
     CLOSE (40504)
     !
     !
  ENDDO ! q loop
  !
  !
  IF (save_m_mat) THEN
     IF (my_pool_id .EQ. ionode_id) CALL SYSTEM ('rm '//TRIM(tmp_dir)//TRIM(prefix)//'.epf17_*')
  ENDIF
  !
  !
#ifdef __PARA
  CALL mp_barrier (inter_pool_comm)
  CALL mp_sum (gammai_mode_all,inter_pool_comm)
#ENDIF
  !
  !
  ! Export
  CALL export_rate_ph ()
  !
  !
#ifdef __PARA
  CALL mp_barrier (inter_pool_comm)
#ENDIF
  !
  !
  !
END SUBROUTINE phdrag_shuffle
  !
  !
!---------------------------------
SUBROUTINE read_shengbte ()
!---------------------------------
#INCLUDE "f_defs.h"
  USE kinds,             ONLY : DP
  USE pwcom,            ONLY : nbnd, nks, nkstot, isk, et, xk, alat, nelec
  USE ions_base,         ONLY : nat
  USE phcom,             ONLY : nmodes
  USE cell_base,         ONLY : at, bg
  USE epwcom,            ONLY : nbndsub, lpolar, nbndskip, vme, eig_read, & 
                                nkf1, nkf2, nkf3, nqf1, nqf2, nqf3, bte, epdim, &
                                edos_read, asr_eph, smearing, ph_read, ephl_read, &
                                epdope, nepdope, save_m_mat, save_m_matw, &
                                save_m_ph, save_t_el, epthick, phdrag, neptemp, nepdope, eptemp
  USE elph2,             ONLY : uf_ful, wf_ful, vph_ful, wf_all, uf_all, vph_all, wf_irr, uf_irr, vph_irr, &
                                nrr_k, nrr_q, cu, cuq, irvec, ndegen_k, ndegen_q, &
                                wslen, chw, chw_ks, cvmew, cdmew, rdw, epmatwp, epmatwe, epmatq, &
                                wf, etf, etf_k, etf_ks, xqf, xkf, wkf, wqf, &
                                dynq, nqtotf, nkqf, epf17, nkf, nqf, et_ks, &
                                ibndmin, ibndmax, lambda_all, dmec, dmef, vmef, &
                                gamma_all, nkqtotf, epsi, zstar, efnew, &
                                ! THL
                                nbnd_red, cbnd_emin, vbnd_emax, ef_m, delta_egap, outside_gap, cfsthick, vfsthick, &
                                sigmai_mode_all_abs, sigmai_mode_all_emi, &
                                etf_ful, vel_ful, epmatwp_asr, eph_vogl, gammai_mode_all
  USE constants_epw, ONLY : ryd2ev, ryd2mev, rydcm1, one, two, twopi, pi, au2cm, au2ps, ryd2thz, kB
  USE bte_var
  USE bte_func
  USE tetrahedron     
#ifdef __PARA
  USE io_files,      ONLY : tmp_dir, prefix
  USE mp,            ONLY : mp_barrier, mp_bcast, mp_sum
  USE io_global,     ONLY : ionode_id, stdout, ionode
  USE mp_global,     ONLY : my_pool_id, nproc_pool, intra_image_comm, &
                            inter_pool_comm, me_pool, root_pool, intra_pool_comm, npool
  USE mp_world,      ONLY : mpime
#ENDIF
  !
  IMPLICIT NONE
  !
  INTEGER                       :: ik, ik0, ik_red, ik0_red, iq, iq_irr, iq_ful, iq_irr_red, iq0, iq_red
  INTEGER                       :: ibnd, jbnd, ibnd0, jbnd0
  INTEGER                       :: imode, mu, nu, i, j, k
  !
  ! For ShengBTE phv, Qian, Nov 2017
  INTEGER,ALLOCATABLE :: SBTE_q_ful(:),SBTE_q_1irr(:)
  REAL(kind=DP),ALLOCATABLE :: SBTE_q_crys_ful(:,:),S_fre_irr(:,:,:),S_sct_irr(:,:,:)
  CHARACTER(LEN=256)    :: S_sct_itemp_fmt, string_temp
  INTEGER :: f2f, S_temp, nlines, io, itemp
  !
  WRITE (stdout,'(/5x,a)') 'Reading phonon-phonon scattering rate from file BTE.w'
  !
  !
   ALLOCATE (ph_rate_ful(neptemp,nmodes,nq_ful))
   ph_rate_ful = 0.0d0
  
  ! read in qpoints_ful in ShengBTE and match its counterpart in nq_ful in epw
     ALLOCATE (SBTE_q_crys_ful(3,nq_ful))
     ALLOCATE (S_fre_irr(neptemp,nmodes,nq_irr))
     ALLOCATE (S_sct_irr(neptemp,nmodes,nq_irr))
     ALLOCATE (SBTE_q_ful(nq_ful))
     ALLOCATE (SBTE_q_1irr(nq_ful))
     !
     ! local
     SBTE_q_crys_ful = 0.0d0
     S_fre_irr = 0.0d0
     S_sct_irr = 0.0d0
     SBTE_q_ful = 0
     SBTE_q_1irr = 0
     !
     WRITE (stdout,'(/5x,a)') 'ShengBTE temp variables allocated'
     !
     !
     OPEN (22334,FILE='BTE.qpoints_full',STATUS='old')
     nlines = 0
     DO iq = 1, nq_ful
       READ (22334,*,iostat=io) SBTE_q_ful(iq), SBTE_q_1irr(iq), SBTE_q_crys_ful(1:3,iq)
       IF (io/=0) EXIT
       nlines = nlines + 1
     ENDDO
     ! check if lengthes of the 2 meshes matches each other
     IF (nlines .EQ. nq_ful) WRITE (stdout,'(/5x,a)') 'ShengBTE and EPW meshes matched, BTE.qpoints loaded'
     IF (nlines .NE. nq_ful) CALL errore ('meshes_not_matching','Error: qpoints_irr # in ShengBTE is different from EPW nq_irr',1)
     CLOSE (22334)
     !
     !
     DO itemp = 1, neptemp
     S_temp = NINT(eptemp(itemp)/kB)
     WRITE(string_temp, *) S_temp 
     S_sct_itemp_fmt = 'T'//TRIM(ADJUSTL(string_temp))//'K/BTE.w'
     OPEN (44556,FILE=S_sct_itemp_fmt,STATUS='old')
     !
      DO imode = 1, nmodes
     !
       DO iq = 1, nq_irr
     !
        READ (44556,*) S_fre_irr(itemp,imode,iq), S_sct_irr(itemp,imode,iq) ! ph scat_rate in ps-1
     !
        IF (iq .EQ. 500) WRITE (stdout,'(/5x,a,i4,a,i3,a,es17.10)') 'Check: Reading scat_rate, temp = ',S_temp,&
                              'K, imode =',imode, ', iq=500, scat_rate(THz) = ', S_sct_irr(itemp,imode,iq)
       ENDDO
      ENDDO
     CLOSE (44556)
     ENDDO
     WRITE (stdout,'(/5x,a)') 'BTE.w from ShengBTE loaded'

    
       DO iq = 1, nq_ful   ! qpoint_full are originally matched, sq = irr2ful_q(eq)
       !
                 DO imode = 1, nmodes !ShengBTE nmodes (w low to high) and EPW nmodes sequence should match
                 !
                 !
                  DO itemp = 1, neptemp
                 !
           ph_rate_ful(itemp,imode,iq) =  S_sct_irr(itemp,imode,SBTE_q_1irr(iq))  ! in ps-1
                 !
                 !
                  ENDDO ! neptemp
                 !
                 ENDDO  ! nmodes
                 !
                 ! ShengBTE includes no doping but temperature dependence in phonon lifetime 
              !
        !
        ENDDO ! iq  
     !
  !
  ph_rate_ful = ph_rate_ful * au2ps ! [1/ps]->[1/au_s]
  !
        DEALLOCATE (SBTE_q_1irr)
        DEALLOCATE (SBTE_q_ful)
        DEALLOCATE (SBTE_q_crys_ful)
        DEALLOCATE (S_fre_irr)
        DEALLOCATE (S_sct_irr)
        WRITE (stdout,'(/5x,a)') 'ShengBTE temp variables deallocated'
  !
  !
  ! ph_rate_ful with symmetry
  OPEN (24877,FILE='BTE/META/equiv_q',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='old')
  DO iq = 1, nq_ful
    DO itemp = 1, neptemp
     READ (24877,REC=iq) f2f 
     IF (iq .NE. f2f) ph_rate_ful(itemp,:,iq) = ph_rate_ful(itemp,:,f2f)
    ENDDO
  ENDDO
  CLOSE (24877) 
  !
  !
  IF (phdrag) THEN
  IF (ionode) THEN
        !
        DO itemp = 1, neptemp
        !
         S_temp = NINT(eptemp(itemp)/kB)
        !
         WRITE(string_temp, *) S_temp
        !
         CALL SYSTEM ('mkdir -p BTE/EPCHECK/T'//TRIM(ADJUSTL(string_temp))//'K') 
        !
         OPEN (12345,FILE='BTE/EPCHECK/T'//TRIM(ADJUSTL(string_temp))//'K/BTE_w_readcheck.dat',STATUS='replace')
        !
        ! Check ShengBTE input is read correctly
        !
          DO iq = 1, nq_irr
        !
  ! IF (.NOT. phdrag) WRITE (stdout,'(/5x,a,i10)') 'Debugging bte=2, nmodes = ', nmodes
        !
            DO nu = 1, nmodes
        !
  ! IF (.NOT. phdrag) WRITE (stdout,'(/5x,a,i10)') 'Debugging bte=2, imode = ', nu
        !
              WRITE (12345,'(2es16.6)') wf_irr(nu,iq)*ryd2thz, ph_rate_ful(itemp,nu,irr2ful_q(iq))/au2ps 
        !
            ENDDO
        !
          ENDDO
        !
         CLOSE (12345)
        ENDDO
  ENDIF
  !
  !
#ifdef __PARA
  CALL mp_barrier (inter_pool_comm)
#ENDIF
  !
  ENDIF
!---------------------------------
END SUBROUTINE read_shengbte
!---------------------------------

!---------------------------------
SUBROUTINE read_alloy_el ()
!---------------------------------
#INCLUDE "f_defs.h"
  USE kinds,             ONLY : DP
  USE pwcom,            ONLY : nbnd, nks, nkstot, isk, et, xk, alat, nelec
  USE ions_base,         ONLY : nat
  USE phcom,             ONLY : nmodes
  USE cell_base,         ONLY : at, bg
  USE epwcom,            ONLY : nbndsub, lpolar, nbndskip, vme, eig_read, & 
                                nkf1, nkf2, nkf3, nqf1, nqf2, nqf3, bte, epdim, &
                                edos_read, asr_eph, smearing, ph_read, ephl_read, &
                                epdope, nepdope, save_m_mat, save_m_matw, &
                                save_m_ph, save_t_el, epthick, phdrag, neptemp, nepdope, eptemp
  USE elph2,             ONLY : uf_ful, wf_ful, vph_ful, wf_all, uf_all, vph_all, wf_irr, uf_irr, vph_irr, &
                                nrr_k, nrr_q, cu, cuq, irvec, ndegen_k, ndegen_q, &
                                wslen, chw, chw_ks, cvmew, cdmew, rdw, epmatwp, epmatwe, epmatq, &
                                wf, etf, etf_k, etf_ks, xqf, xkf, wkf, wqf, &
                                dynq, nqtotf, nkqf, epf17, nkf, nqf, et_ks, &
                                ibndmin, ibndmax, lambda_all, dmec, dmef, vmef, &
                                gamma_all, nkqtotf, epsi, zstar, efnew, etf_all, &
                                ! THL
                                nbnd_red, cbnd_emin, vbnd_emax, ef_m, delta_egap, outside_gap, cfsthick, vfsthick, &
                                sigmai_mode_all_abs, sigmai_mode_all_emi, &
                                etf_ful, vel_ful, epmatwp_asr, eph_vogl, gammai_mode_all
  USE constants_epw, ONLY : ryd2ev, ryd2mev, rydcm1, one, two, twopi, pi, au2cm, au2ps, ryd2thz, kB
  USE bte_var
  USE bte_func
  USE tetrahedron     
#ifdef __PARA
  USE io_files,      ONLY : tmp_dir, prefix
  USE mp,            ONLY : mp_barrier, mp_bcast, mp_sum
  USE io_global,     ONLY : ionode_id, stdout, ionode
  USE mp_global,     ONLY : my_pool_id, nproc_pool, intra_image_comm, &
                            inter_pool_comm, me_pool, root_pool, intra_pool_comm, npool
  USE mp_world,      ONLY : mpime
#ENDIF
  !
  IMPLICIT NONE
  !
  INTEGER                       :: ik, ik0, ik_red, ik0_red, iq, iq_irr, iq_ful, iq_irr_red, iq0, iq_red
  INTEGER                       :: ibnd, jbnd, ibnd0, jbnd0
  INTEGER                       :: imode, mu, nu, i, j, k
  !
  ! For alloy electron scattering rate, Qian, June 2018
  INTEGER :: S_temp, io, itemp
  !
  WRITE (stdout,'(/5x,a)') 'Reading alloy-electron scattering rate from file alel_rate.txt'
  !
  !
   ALLOCATE (alelrate(4000,2))
   alelrate = 0.0d0
  
  ! read in alel_rate file 

     WRITE (stdout,'(/5x,a)') 'al_el temp variables allocated'
     !
     !
     OPEN (7655,FILE='alel_rate.txt',STATUS='old')
     DO iq = 1, 4000
       READ (7655,*,iostat=io) alelrate(iq,1), alelrate(iq,2)
     ENDDO
     CLOSE (7655) 
     !
     WRITE (stdout,'(/5x,a,i5)') 'al_el temp variables allocated, io read state = ', io
  !
    alelrate(:,2) = alelrate(:,2) * au2ps ! Thz = [1/ps]->[1/au_s]
  !
     IF (io .LE. 0) WRITE (stdout,'(/5x,a)') 'Alloy-electron scattering rate read'
  !
  IF (phdrag) THEN
  IF (ionode) THEN
        !
        !
        !
         OPEN (3345,FILE='BTE/EPCHECK/alelrate.dat',STATUS='replace')
        !
        ! Check if alel_rate input is read correctly
        !
          DO iq = 1, 4000
        !
              WRITE (3345,'(f12.6, es16.6)') alelrate(iq,1), alelrate(iq,2)/au2ps 
        !
          ENDDO
        !
         CLOSE (3345)
  ENDIF
  !
#ifdef __PARA
  CALL mp_barrier (inter_pool_comm)
#ENDIF
  !
  ENDIF
!---------------------------------
END SUBROUTINE read_alloy_el
!---------------------------------

!---------------------------------
SUBROUTINE phph_export ()
!---------------------------------
#INCLUDE "f_defs.h"
  USE kinds,             ONLY : DP
  USE pwcom,            ONLY : nbnd, nks, nkstot, isk, et, xk, alat, nelec
  USE ions_base,         ONLY : nat
  USE phcom,             ONLY : nmodes
  USE cell_base,         ONLY : at, bg
  USE epwcom,            ONLY : nbndsub, lpolar, nbndskip, vme, eig_read, & 
                                nkf1, nkf2, nkf3, nqf1, nqf2, nqf3, bte, epdim, &
                                edos_read, asr_eph, smearing, ph_read, ephl_read, &
                                epdope, nepdope, save_m_mat, save_m_matw, &
                                save_m_ph, save_t_el, epthick, phdrag, neptemp, nepdope, eptemp
  USE elph2,             ONLY : uf_ful, wf_ful, vph_ful, wf_all, uf_all, vph_all, wf_irr, uf_irr, vph_irr, &
                                nrr_k, nrr_q, cu, cuq, irvec, ndegen_k, ndegen_q, &
                                wslen, chw, chw_ks, cvmew, cdmew, rdw, epmatwp, epmatwe, epmatq, &
                                wf, etf, etf_k, etf_ks, xqf, xkf, wkf, wqf, &
                                dynq, nqtotf, nkqf, epf17, nkf, nqf, et_ks, &
                                ibndmin, ibndmax, lambda_all, dmec, dmef, vmef, &
                                gamma_all, nkqtotf, epsi, zstar, efnew, &
                                ! THL
                                nbnd_red, cbnd_emin, vbnd_emax, ef_m, delta_egap, outside_gap, cfsthick, vfsthick, &
                                sigmai_mode_all_abs, sigmai_mode_all_emi, &
                                etf_ful, vel_ful, epmatwp_asr, eph_vogl, gammai_mode_all
  USE constants_epw, ONLY : ryd2ev, ryd2mev, rydcm1, one, two, twopi, pi, au2cm, au2ps, ryd2thz, kB
  USE bte_var
  USE bte_func
  USE tetrahedron     
#ifdef __PARA
  USE io_files,      ONLY : tmp_dir, prefix
  USE mp,            ONLY : mp_barrier, mp_bcast, mp_sum
  USE io_global,     ONLY : ionode_id, stdout, ionode
  USE mp_global,     ONLY : my_pool_id, nproc_pool, intra_image_comm, &
                            inter_pool_comm, me_pool, root_pool, intra_pool_comm, npool
  USE mp_world,      ONLY : mpime
#ENDIF
  !
  IMPLICIT NONE
  !
  INTEGER                       :: ik, ik0, ik_red, ik0_red, iq, iq_irr, iq_ful, iq_irr_red, iq0, io, itemp
  INTEGER                       :: ibnd, jbnd, ibnd0, jbnd0
  INTEGER                       :: imode, mu, nu, i, j, k
  !
  ! For ShengBTE phv, Qian, Nov 2017
  !
  INTEGER,ALLOCATABLE :: SBTE_q_ful(:),SBTE_q_1irr(:)
  REAL(kind=DP),ALLOCATABLE :: SBTE_q_crys_ful(:,:),S_fre_irr(:,:,:),S_sct_irr(:,:,:)
  REAL(kind=DP) :: xxq(3), xxk(3), xkk(3), xkq(3), xq_fbz(3), length
  CHARACTER(LEN=256)    :: S_sct_itemp_fmt, string_temp
  INTEGER :: iq_red, S_temp, nlines
  !
  ! Output phonon-phonon scattering rate corresponding to format of total gammi output file
  !
  IF (ionode) THEN
     DO itemp = 1, neptemp
        S_temp = NINT(eptemp(itemp)/kB)
        WRITE(string_temp, *) S_temp
        ! 
        CALL SYSTEM ('mkdir -p BTE/EPCHECK/T'//TRIM(ADJUSTL(string_temp))//'K') 
        !
        OPEN (2224,FILE='BTE/EPCHECK/T'//TRIM(ADJUSTL(string_temp))//'K/BTE_w_qlength.dat',STATUS='replace')
        OPEN (7777,FILE='BTE/META/irr2ful_q',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='old')
        OPEN (8888,FILE='BTE/META/xqf_fbz_cart',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='old')
        ! 
         DO iq_irr = 1, nq_irr
           READ (7777,REC=iq_irr) iq_ful
           READ (8888,REC=iq_ful) xq_fbz(1:3)
           length = SQRT(DOT_PRODUCT(xq_fbz,xq_fbz))
           DO imode = 1, nmodes  
           IF (phdrag) THEN
            WRITE (2224,'(i8,f12.6,f11.4,es14.4)') iq_irr, length, wf_ful(imode,1,1,iq_ful)*rydcm1,&
                                                 ph_rate_ful(itemp,imode,iq_ful)/au2ps ! [cm-1] [THz]
           ELSE
            WRITE (2224,'(i8,f12.6,f11.4,es14.4)') iq_irr, length, wf_irr(imode,iq_irr)*rydcm1,&
                                                 ph_rate_ful(itemp,imode,iq_ful)/au2ps ! [cm-1] [THz]
           !
           ENDIF
           ENDDO
         ENDDO  
        !  
        !
        CLOSE (2224)
        CLOSE (7777)
        CLOSE (8888)        
        !
     ENDDO
  !
  ENDIF
  !
  !
#ifdef __PARA
  CALL mp_barrier (inter_pool_comm)
#ENDIF
!---------------------------------
END SUBROUTINE phph_export 
!---------------------------------


