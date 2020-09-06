  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                        
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE elphon_shuffle_wrap
  !-----------------------------------------------------------------------
  !
  ! Electron-phonon calculation with Wannier functions: load all phonon q's
  !
  ! 09/2009 This subroutine is the main driver of the electron-phonon 
  ! calculation. It first calculates the electron-phonon matrix elements
  ! on the coarse mesh and then passes the data off to ephwann_shuffle
  ! to perform the interpolation.
  !
  !-----------------------------------------------------------------------
  !
#include "f_defs.h"
  !
#ifdef __PARA
  USE mp_global,     ONLY : my_pool_id, inter_pool_comm, root_pool, &
                            intra_pool_comm,npool
  USE mp_world,      ONLY : mpime, root
  USE mp,            ONLY : mp_barrier, mp_bcast, mp_sum, mp_max
  USE io_global,     ONLY : ionode_id
#endif
  USE us,            ONLY : nqxq, dq, qrad
  USE gvect,         ONLY : gcutm
  USE gvecs,         ONLY : nls
  USE cellmd,        ONLY : cell_factor
  USE uspp_param,    ONLY : lmaxq, nbetam
  USE io_files,      ONLY : prefix, tmp_dir, iunigk, diropn, seqopn
  USE wavefunctions_module, ONLY: evc
  USE ions_base,     ONLY : nat, nsp, tau, ityp
  USE control_flags, ONLY : iverbosity
  USE io_global,     ONLY : stdout, ionode
  USE io_epw,        ONLY : iuepb, iueimpb
  USE kinds,         ONLY : DP
  USE pwcom,         ONLY : et, xk, ibrav, igk, nks, nbnd, nkstot, ngm, &
                            omega, g, npw, igk
  USE wvfct,         ONLY : npwx
  USE cell_base,     ONLY : at, bg
  USE symm_base,     ONLY : irt, s, nsym, ftau, sname, invs, &
                            s_axis_to_cart, sr, nrot, copy_sym, set_sym_bl, find_sym, inverse_s
  USE start_k,       ONLY : nk1, nk2, nk3
  USE phcom,         ONLY : dpsi, dvpsi, evq, nq1, nq3, nq2, iuwfc
  USE units_ph,      ONLY : lrwfc
  USE qpoint,        ONLY : igkq
  USE control_lr,    ONLY : lgamma 
  USE qpoint,        ONLY : xq
  USE modes,         ONLY : nmodes
  USE lr_symm_base,  ONLY : minus_q, rtau, gi, gimq, irotmq, nsymq, invsymq
  USE epwcom,        ONLY : epbread, epbwrite, epwread, phinterp, &
                            phonselfen, elecselfen, nbndsub, elinterp,    &
                            iswitch, kmaps, nest_fn, eig_read, &
                            band_plot, specfun, dvscf_dir, lpolar, bte, neptemp, nepdope, phdrag, &
                            eimp_mode, dielec, eimpbread, eimpbwrite, defectname, &
                            dvimpsr, epbjump, epbrestore, dvimpq, alloy_pot
  USE elph2,         ONLY : epmatq, dynq, sumr, et_all, xk_all, et_mb, et_ks, &
                            zstar, epsi, cu, cuq, lwin, lwinq, bmat, &
                            model_charge, alat_imp, scella, scellb, scellc, &
                            nfftmesh, rp, dvimp, eimpmatq, imp_meshmap, &
                            scell_vol, def_pos, bg_max, dvGr, dvimp_q, &
                            map_l2g, map_g2l, ng_max, evc_k, umat, umat_all, &
                            ig_max, dvimp_qc, bg_max_red, ng_max_red, &
                            map_kq2k
  USE constants_epw, ONLY : ryd2ev, ci, twopi, cone, czero
  USE fft_base,      ONLY : dfftp, dffts
  USE control_ph,    ONLY : u_from_file
  USE noncollin_module, ONLY : m_loc, npol, nspin_mag, noncolin
  use lsda_mod,      only : nspin
  use spin_orb,      only : domag
  USE iotk_module,   ONLY : iotk_open_read, iotk_scan_dat, iotk_free_unit,&
                            iotk_close_read
#ifdef __NAG
  USE f90_unix_io,    ONLY : flush
#endif
  implicit none

  !
  integer :: sym_smallq(48) 
  !
  real(kind=DP), allocatable :: xqc_irr(:,:), wqlist_irr(:), xqc(:,:), wqlist(:)
  ! the qpoints in the irr wedge
  ! the corresponding weigths
  ! the qpoints in the uniform mesh
  ! the corresponding weigths
  integer :: nqc_irr, nqc, max, nqxq_tmp, ibnd, ik, ios, &
             dummy1, dummy2, ik_start, ik_stop
  ! number of qpoints in the irreducible wedge
  ! number of qpoints on the uniform grid
  ! 
  ! symmetry-related variables
  !
  integer :: gmapsym(ngm,48)
  !  correspondence G -> S(G)
  complex(kind=DP) :: eigv (ngm, 48)
  ! e^{ iGv} for 1...nsym (v the fractional translation)
  complex(kind=DP) :: cz1( nmodes, nmodes), cz2(nmodes, nmodes)
  !  the eigenvectors for the first q in the star
  !  the rotated eigenvectors, for the current q in the star
  !
  integer :: nq, isq (48), imq 
  ! degeneracy of the star of q
  ! index of q in the star of a given sym.op.
  ! index of -q in the star of q (0 if not present)
  integer :: sym_sgq(48)
  ! the symmetries giving the q point iq in the star
  real(kind=DP) :: sxq (3, 48), et_tmp(nbnd, nkstot)
  ! list of vectors in the star of q
  integer :: i, j, iq, iq_irr, isym, &
     iq_first, jsym, ism1, nsq, ipol, jpol, ierr, iunpun
  real(kind=DP) xq0(3), aq(3), saq(3), raq(3), ft1, ft2, ft3
!  REAL(DP) :: w2(3*nat) ! dummy here
  logical :: sym(48),  eqvect_strict, nog, symmo, exst
  character (len=256) :: tempfile, dirname,filename
#ifdef __PARA
  character (len=3) :: filelab
#endif
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  !
  ! THL: BTE
  INTEGER          :: itemp, idope, imode
  CHARACTER(LEN=3) :: itemp_num, idope_num, imode_num
  LOGICAL          :: dir_exist
  !
  ! JZ: for e-imp
  integer :: nfftx_sc, nffty_sc, nfftz_sc, nfft_sc, ifft_star, &
             ifft_stop, finished_iq_irr, center_index, ind_r, &
             ind_g, irx, iry, irz, igx, igy, igz, imesh, ig, &
             i1, j1, k1, ig_global, ipooltmp, ig_nrmax, &
             nkk, nkk_abs, ipool, ind_gx, ind_gy, ind_gz, &
             iG1, iG2, iG3, iq1, iq2, iq3, nnq, ig_shift, jbnd
  real(kind=DP) :: temp(3), dvol, rp_(3), dielec0, &
                   rpl, rp_max, dv_tot_far, vcoul_far, &
                   rp_redm(3), xk_gamma(3), xk_x(3), &
                   dvimp_, r0(3), G0(3), dotGr, xq_cart(3), &
                   sc_un_ratio, xqG(3), norm0, zero_vect(3), &
                   dvimp_r, dvimp_i, xq1, xq2, xq3, xk_l(3)
  real(kind=DP) :: t1, t2, t3_1, t3_2, t3_3, t3_4, t3_5, t3_6, t3_7, t3_8, &
                   t12, t23, t34, t45, t56, t67, t78, t89
  real(kind=DP), allocatable :: vcoul_model(:), rp_red(:,:), rp_c(:,:)
  real(kind=DP), external :: V_coulomb
  complex(kind=DP) :: evctmp(nbnd)
  complex(kind=DP), allocatable :: eemat_check(:,:,:), evck_remap(:), evck(:,:)
  logical :: gamma_find, x_find, l_find
  COMPLEX(DP),EXTERNAL :: ZDOTC

  !
  !-------------------------------------------------------------------
  ! THL: prepare folders for BTE calculation
  !-------------------------------------------------------------------
  !
  !
  IF (bte .EQ. 0 .OR. bte .EQ. 1 .OR. bte .EQ. -1 .OR. bte .EQ. 19 .OR. bte .EQ. 18) THEN
     !
     CALL mp_barrier(inter_pool_comm)
     !
     IF (my_pool_id .EQ. ionode_id) THEN
        !
        INQUIRE(DIRECTORY='BTE',EXIST=dir_exist)
        IF (dir_exist) CALL SYSTEM ('rm -r BTE')
        !
        CALL SYSTEM ('mkdir BTE')
        CALL SYSTEM ('mkdir BTE/SIGMAI')
        CALL SYSTEM ('mkdir BTE/META')
        CALL SYSTEM ('mkdir BTE/EPCHECK')
        IF (phdrag) CALL SYSTEM ('mkdir BTE/GAMMAI')
        !
        DO itemp = 1, neptemp
           DO idope = 1, nepdope
              !
              WRITE(itemp_num,'(i3)') itemp
              WRITE(idope_num,'(i3)') idope
              CALL SYSTEM ( 'mkdir BTE/T'//TRIM(ADJUSTL(itemp_num))//'_N'//TRIM(ADJUSTL(idope_num)) )
              !
              DO imode = 1, nmodes
                 !
                 WRITE(imode_num,'(i3)') imode
                 CALL SYSTEM ( 'mkdir BTE/SIGMAI/T'//TRIM(ADJUSTL(itemp_num))//'_N'//TRIM(ADJUSTL(idope_num))//'_ph'//TRIM(ADJUSTL(imode_num)) )
                 IF (phdrag) CALL SYSTEM ( 'mkdir BTE/GAMMAI/T'//TRIM(ADJUSTL(itemp_num))//'_N'//TRIM(ADJUSTL(idope_num))//'_ph'//TRIM(ADJUSTL(imode_num)))
                 !
              ENDDO
              !
           ENDDO
        ENDDO
        !
     ENDIF 
     !
     CALL mp_barrier(inter_pool_comm)
     !
  ELSEIF (bte .EQ. 10) THEN
     !
     INQUIRE(DIRECTORY='BTE/META',EXIST=dir_exist)
     IF (.NOT. dir_exist) CALL errore('elphon_shuffle_wrap','There is no META directory',1)
     !
     INQUIRE(DIRECTORY='BTE/SIGMAI',EXIST=dir_exist)
     IF (.NOT. dir_exist) CALL errore('elphon_shuffle_wrap','There is no SIGMAI directory',1)
     !
     CALL mp_barrier(inter_pool_comm)
     !
     GOTO 96969
     !
  ELSEIF (bte .EQ. 2 .OR. bte .EQ. 29) THEN
     !
     IF (my_pool_id .EQ. ionode_id) THEN
        !
        INQUIRE(DIRECTORY='BTE',EXIST=dir_exist)
        IF (dir_exist) CALL SYSTEM ('rm -r BTE')
        !
        CALL SYSTEM ('mkdir BTE')
        CALL SYSTEM ('mkdir BTE/GAMMAI')
        CALL SYSTEM ('mkdir BTE/META')
        CALL SYSTEM ('mkdir BTE/EPCHECK')
        !
        DO itemp = 1, neptemp
           DO idope = 1, nepdope
              !
              WRITE(itemp_num,'(i3)') itemp
              WRITE(idope_num,'(i3)') idope
              CALL SYSTEM ( 'mkdir BTE/T'//TRIM(ADJUSTL(itemp_num))//'_N'//TRIM(ADJUSTL(idope_num)) )
              !
              DO imode = 1, nmodes
                 !
                 WRITE(imode_num,'(i3)') imode
                 CALL SYSTEM ( 'mkdir BTE/GAMMAI/T'//TRIM(ADJUSTL(itemp_num))//'_N'//TRIM(ADJUSTL(idope_num))//'_ph'//TRIM(ADJUSTL(imode_num)) )
                 !
              ENDDO
              !
           ENDDO
        ENDDO
        !
     ENDIF 
     !
     CALL mp_barrier(inter_pool_comm)
     !
  ELSEIF (bte .EQ. 3) THEN
     !
     IF (my_pool_id .EQ. ionode_id) THEN
        !
        INQUIRE(DIRECTORY='BTE',EXIST=dir_exist)
        IF (dir_exist) CALL SYSTEM ('rm -r BTE')
        !
        CALL SYSTEM ('mkdir BTE')
        CALL SYSTEM ('mkdir BTE/WEIGHT')
        CALL SYSTEM ('mkdir BTE/META')
        CALL SYSTEM ('mkdir BTE/EPCHECK')
        CALL SYSTEM ('mkdir BTE/T1_N1')
        !              
        DO imode = 1, nmodes
           !
           WRITE(imode_num,'(i3)') imode
           CALL SYSTEM ( 'mkdir BTE/WEIGHT/T1_N1_ph'//TRIM(ADJUSTL(imode_num)) )
           !
        ENDDO
        !
     ENDIF 
     !
     CALL mp_barrier(inter_pool_comm)
     !
  ELSEIF (bte .EQ. 30) THEN
     !
     INQUIRE(DIRECTORY='BTE/META',EXIST=dir_exist)
     IF (.NOT. dir_exist) CALL errore('elphon_shuffle_wrap','There is no META directory',1)
     !
     INQUIRE(DIRECTORY='BTE/WEIGHT',EXIST=dir_exist)
     IF (.NOT. dir_exist) CALL errore('elphon_shuffle_wrap','There is no WEIGHT directory',1)
     !
     CALL mp_barrier(inter_pool_comm)
     !
     GOTO 96969
     !
  ENDIF
  !
  !-------------------------------------------------------------------
  !
  !
  CALL start_clock ( 'elphon_wrap' )
  !
  IF ( elinterp .and. (.not.phinterp ) ) CALL errore &
        ('elphon_shuffle_wrap','elinterp requires phinterp' ,1)
  !
  ! READ qpoint list from stdin
  !
#ifdef __PARA
  IF (mpime.eq.ionode_id) &
#endif
  READ(5,*) nqc_irr
#ifdef __PARA
  CALL mp_bcast (nqc_irr, ionode_id, inter_pool_comm)
  CALL mp_bcast (nqc_irr, root_pool, intra_pool_comm)
#endif
  allocate ( xqc_irr(3,nqc_irr), wqlist_irr(nqc_irr) )
  allocate ( xqc(3,nq1*nq2*nq3), wqlist(nq1*nq2*nq3) )
  !  
#ifdef __PARA
  IF (mpime.eq.ionode_id) then
#endif
    DO iq = 1, nqc_irr
      READ (5,*) xqc_irr (:,iq), wqlist_irr (iq)
    ENDDO
#ifdef __PARA
  ENDIF
  CALL mp_bcast (xqc_irr, ionode_id, inter_pool_comm)
  CALL mp_bcast (xqc_irr, root_pool, intra_pool_comm)
#endif
  !
  ! fix for uspp
  max = nqxq
  DO iq = 1, nqc_irr
     nqxq_tmp = INT( ( (sqrt(gcutm) + sqrt(xqc_irr(1,iq)**2 + &
          xqc_irr(2,iq)**2 + xqc_irr(3,iq)**2) ) &
          / dq + 4) * cell_factor )
     IF (nqxq_tmp .gt. max)  max = nqxq_tmp
  ENDDO
  IF (max .gt. nqxq) then
     IF (allocated(qrad)) deallocate(qrad)
     allocate (qrad (max, nbetam*(nbetam+1)/2,lmaxq, nsp))
  ENDIF
  IF (nkstot .ne. nk1*nk2*nk3 ) &
       CALL errore('elphon_shuffle_wrap','nscf run inconsistent with epw input',1)  
  !
  ! READ in external electronic eigenvalues. e.g. GW 
  !
  IF ( .not. ALLOCATED(et_ks) ) ALLOCATE(et_ks(nbnd,nks))
  IF ( .not. ALLOCATED(et_mb) ) ALLOCATE(et_mb(nbnd,nks))
  et_ks(:,:) = 0.d0
  et_mb(:,:) = 0.d0
  IF (eig_read) then
#ifdef __PARA
  IF (mpime.eq.ionode_id) then
#endif
   WRITE (stdout,'(5x,a,i5,a,i5,a)') "Reading external electronic eigenvalues (", &
        nbnd, ",", nkstot,")"
   tempfile=trim(prefix)//'.eig'
   OPEN(1, file=tempfile, form='formatted', action='read', iostat=ios)
   IF (ios /= 0) CALL errore ('elphon_shuffle_wrap','error opening' // tempfile, 1)
    DO ik = 1, nkstot
       DO ibnd = 1, nbnd
          READ (1,*) dummy1, dummy2, et_tmp (ibnd,ik)
          IF (dummy1.ne.ibnd) CALL errore('elphon_shuffle_wrap', "Incorrect eigenvalue file", 1)
          IF (dummy2.ne.ik)   CALL errore('elphon_shuffle_wrap', "Incorrect eigenvalue file", 1)
       ENDDO
    ENDDO
    CLOSE(1)
    ! from eV to Ryd
    et_tmp = et_tmp / ryd2ev
#ifdef __PARA
    ENDIF
    CALL mp_bcast (et_tmp, ionode_id, inter_pool_comm)
    CALL mp_bcast (et_tmp, root_pool, intra_pool_comm)
#endif
    !
    CALL ckbounds(ik_start, ik_stop)
    et_ks(:,:)  = et(:,1:nks)
    et(:,1:nks) = et_tmp(:,ik_start:ik_stop)
    et_mb(:,:)  = et(:,1:nks)
 ENDIF
  !
  ! compute coarse grid dipole matrix elements.  Very fast 
  CALL compute_pmn_para
  !
  !  gather electronic eigenvalues for subsequent shuffle
  !  
  allocate ( xk_all( 3, nkstot) , et_all( nbnd, nkstot) )
  CALL poolgather (    3, nkstot, nks, xk(:,1:nks),      xk_all)
  CALL poolgather ( nbnd, nkstot, nks, et(1:nbnd,1:nks), et_all)
  !
  !
  IF (.not.kmaps) then
     CALL start_clock('kmaps')
     xq0(:)=0.d0
     CALL createkmap_pw2(xk_all,nkstot, xq0)
     CALL stop_clock('kmaps')
     CALL print_clock('kmaps')
  ELSE
     ! 
     ! 26/06/2012 RM
     ! if we do not have epmatq already on file then epbread=.false.
     ! .kgmap is used from disk and .kmap is regenerated for each q-point 
     ! 
     WRITE (stdout, '(/5x,a)')  'Using kmap and kgmap from disk'
  ENDIF
#ifdef __PARA
  CALL mp_barrier(inter_pool_comm)
#endif
  !  if we start with a gamma point calculation, ../PW/set_kplusq.f90
  !  is not active and the gmap has not been produced...
  !
  IF (lgamma) CALL errore &
    ('elphon_shuffle_wrap','tshuffle2 requires q!=0 starting nscf calculation',1)
  !
  ! NOTE: checked, sum of [nks] over pools gives nk1*nk2*nk3
!  write(*,*) ' nks =', nks
  !
  !  allocate dynamical matrix and ep matrix for all q's
  !
  allocate ( dynq (nmodes, nmodes, nq1*nq2*nq3), &
       epmatq (nbnd, nbnd, nks, nmodes, nq1*nq2*nq3), &
       eimpmatq (nbnd, nbnd, nks, nq1*nq2*nq3), &
       epsi(3,3), zstar(3,3,nat), bmat(nbnd, nbnd, nks, nq1*nq2*nq3), &
       cu ( nbnd, nbndsub, nks), cuq ( nbnd, nbndsub, nks), & 
       lwin ( nbnd, nks ), lwinq ( nbnd, nks ) )
  !
  epsi=0.d0
  zstar=0.d0
  !
  ! SP: The symmetries are now consistent with QE 5. This means that the order of the q in the star
  !     should be the same as in the .dyn files produced by QE 5.
  ! 
  !     First we start by setting up the lattice & crystal symm. as done in PH/init_representations.f90
  ! 
  ! ~~~~~~~~ setup bravais lattice symmetry ~~~~~~~~
  CALL set_sym_bl ( ) ! This should define the s matrix
  WRITE(stdout, '(5x,a,i3)') "Symmetries of bravais lattice: ", nrot
  !
  ! ~~~~~~~~ setup crystal symmetry ~~~~~~~~ 
  CALL find_sym ( nat, tau, ityp, dfftp%nr1,dfftp%nr2,dfftp%nr3, .false., m_loc )
  WRITE(stdout, '(5x,a,i3)') "Symmetries of crystal:         ", nsym
  !   
  ! The following loop is required to propertly set up the symmetry matrix s. 
  ! We here copy the calls made in PH/init_representations.f90 to have the same s as in QE 5.
  DO iq_irr = 1, nqc_irr
    xq = xqc_irr(:,iq_irr)
    CALL set_small_group_of_q(nsymq,invsymq,minus_q)
    CALL sgam_ph_new (at, bg, nsym, s, irt, tau, rtau, nat)
    CALL set_giq (xq,s,nsymq,nsym,irotmq,minus_q,gi,gimq)
  ENDDO
  ! 
  ! check dfft mesh   
  write(stdout,*) 'dfftp%nnr: ', dfftp%nnr
  write(stdout,*) 'dfftp%nr1..3: ', dfftp%nr1, dfftp%nr2, dfftp%nr3
  write(stdout,*) 'dffts%nnr: ', dffts%nnr
  write(stdout,*) 'dffts%nr1..3: ', dffts%nr1, dffts%nr2, dffts%nr3

  ! CV: if we read the .fmt files we don't need to read the .epb anymore
  !
  ! Read impurity potential
  !
  if ((.not. eimpbread .and. (eimp_mode == 5 .or. eimp_mode == 6)) .or. &
      (eimp_mode == 7 .or. eimp_mode == 8) ) then
     !
#ifdef __PARA
     IF (mpime .EQ. ionode_id) THEN
#ENDIF
        ! read in impurity potential
        !
        tempfile = 'dv_tot.' // trim(adjustl(defectname)) // '.dat'
        open(unit = 1000, file = tempfile)
        read(1000,*) model_charge
        read(1000,*) nfftx_sc, nffty_sc, nfftz_sc
        read(1000,*) alat_imp      ! this is the primitive lattice constant from DFT impurity calculation, in [aB]
                                   !    should be consistent with the value used
                                   !    here [alat]

        read(1000,*) def_pos(1:3)
        read(1000,*) scella(1:3)   ! in [alat_imp]
        read(1000,*) scellb(1:3)
        read(1000,*) scellc(1:3)
        read(1000,*) nfftmesh      ! here nfftmesh should be equal to nfft_sc

        nfft_sc = nfftx_sc * nffty_sc * nfftz_sc
        if (nfftmesh /= nfft_sc) &
           CALL errore('elphon_shuffle_wrap', 'nfftmesh not equal to nfft_sc',1)
        !
        allocate(rp(3,nfftmesh), dvimp(nfftmesh))
        do i = 1, nfftmesh
           ! rp in [alat_imp], dvimp in [Ryd]
           ! also NOTE: rp has been centered to the defect [rp = rp_true - r_defect]
           read(1000,*) rp(1:3,i), dvimp(i)
        enddo
        !
        close(1000)
#ifdef __PARA
     ENDIF
     !
     CALL mp_bcast (nfftmesh, ionode_id, inter_pool_comm)
     CALL mp_bcast (nfftmesh, root_pool, intra_pool_comm)
     if (.not.allocated(rp))    allocate(rp(3,nfftmesh))
     if (.not.allocated(dvimp)) allocate(dvimp(nfftmesh))
     CALL mp_bcast (rp, ionode_id, inter_pool_comm)
     CALL mp_bcast (rp, root_pool, intra_pool_comm)
     CALL mp_bcast (dvimp, ionode_id, inter_pool_comm)
     CALL mp_bcast (dvimp, root_pool, intra_pool_comm)
     !
     CALL mp_bcast (model_charge, ionode_id, inter_pool_comm)
     CALL mp_bcast (model_charge, root_pool, intra_pool_comm)
     CALL mp_bcast (def_pos, ionode_id, inter_pool_comm)
     CALL mp_bcast (def_pos, root_pool, intra_pool_comm)
     CALL mp_bcast (scella, ionode_id, inter_pool_comm)
     CALL mp_bcast (scella, root_pool, intra_pool_comm)
     CALL mp_bcast (scellb, ionode_id, inter_pool_comm)
     CALL mp_bcast (scellb, root_pool, intra_pool_comm)
     CALL mp_bcast (scellc, ionode_id, inter_pool_comm)
     CALL mp_bcast (scellc, root_pool, intra_pool_comm)
     CALL mp_bcast (alat_imp, ionode_id, inter_pool_comm)
     CALL mp_bcast (alat_imp, root_pool, intra_pool_comm)
     CALL mp_bcast (nfft_sc, ionode_id, inter_pool_comm)
     CALL mp_bcast (nfft_sc, root_pool, intra_pool_comm)
#ENDIF
     !
     call cross(scella, scellb, temp)
     scell_vol = abs(dot_product(scellc, temp)) * (alat_imp**3.d0)   ! in [aB]^3
     dvol = scell_vol/nfft_sc
     ! supercell / unit cell volume ratio
     sc_un_ratio = nint(scell_vol/omega)
     !
     !
     if (.not. dvimpsr) then
        !
        ! correct for long-range Coulomb part
        !
        allocate (vcoul_model(nfftmesh))
        vcoul_model = 0.d0
        !
        CALL para_bounds (ifft_star, ifft_stop, nfftmesh)
        !
        ! for even [lpolar = true], we still use dielec as epsi has not been read in
        dielec0 = dielec
        !
        do i = ifft_star, ifft_stop
           rp_ = rp(:,i)
           ! vcoul_model needs to be in [Ryd]
           vcoul_model(i) = V_coulomb(model_charge, scella, scellb, scellc, alat_imp, rp_, dielec0)
        enddo
        !
#ifdef __PARA
        CALL mp_sum (vcoul_model,inter_pool_comm)
#ENDIF
        !
        ! potential baseline
        !
        rp_max = 0.d0
        dv_tot_far = 0.d0
        vcoul_far = 0.d0
        center_index = 0
        do i = 1, nfftmesh
           rpl = (rp(1,i)**2.d0 + rp(2,i)**2.d0 + rp(3,i)**2.d0)**(0.5d0)
           if (rpl > rp_max) then
              rp_max = rpl
              dv_tot_far = dvimp(i)
              vcoul_far = vcoul_model(i)
           endif
           !
           if (rpl < 1d-5) center_index = i
        enddo

        if (center_index == 0) &
           CALL errore('elphon_shuffle_wrap', 'dv_tot.dat file does not have r=0 point',1)
        !
        ! check potential baseline correction
        !
        write(stdout,*)
        write(stdout,*) ' Potential far away (r, dv_tot, vcoul, eV):', rp_max, dv_tot_far*ryd2ev, vcoul_far*ryd2ev
        write(stdout,*) ' Potential baseline, eV:', (dv_tot_far - vcoul_far)*ryd2ev
        !
        ! we assume alloy potential is from neutral defect
        if (alloy_pot) then
           dvimp = dvimp - dv_tot_far
        else
           dvimp = dvimp - vcoul_model - (dv_tot_far - vcoul_far)
        endif

        ! to avoid arbitrary large number due to Coulomb correction
        dvimp(center_index) = 0.d0
        !
        write(stdout,*) ' [min, max] of corrected impurity potential, eV:', minval(dvimp)*ryd2ev, maxval(dvimp)*ryd2ev

        ! output short-range potential
        !
        tempfile = 'dv_tot_sr.' // trim(adjustl(defectname)) // '.dat'
        open(unit = 1002, file = tempfile)
        do i = 1, nfftmesh
           write(1002,'(3(f15.7,1x),es15.7)') rp(1:3,i), dvimp(i)
        enddo
        close(1002)
        !
     else
        !
        ! read short-range potential
        !
#ifdef __PARA
        IF (mpime .EQ. ionode_id) THEN
#ENDIF
           tempfile = 'dv_tot_sr.' // trim(adjustl(defectname)) // '.dat'
           open(unit = 1002, file = tempfile)
           do i = 1, nfftmesh
              read(1002,*) rp_, dvimp_
              dvimp(i) = dvimp_
           enddo
           close(1002)
#ifdef __PARA
        ENDIF
        CALL mp_bcast (dvimp, ionode_id, inter_pool_comm)
        CALL mp_bcast (dvimp, root_pool, intra_pool_comm)
#ENDIF
        !
     endif
     !
     !
!#ifdef __PARA
!     CALL mp_barrier(inter_pool_comm)
!#endif

!     stop
     !
     ! ---------
     ! mapping between supercell and primitive unit cell position, 
     ! for later electron-impurity scattering matrix calculation
     ! ---------
     !
     ! specifies the mapping from impurity mesh to crystal mesh on primitive unit cell
     allocate (imp_meshmap(nfftmesh, 4))
     allocate (rp_red(nfftmesh, 3))
     
     ! since rp is in unit of [alat_imp], which is primitive unit cell lattice
     ! constant, same as [alat], we will just directly use b1, b2, b3 to project
     ! on proper axes
     !
     do i = 1, nfftmesh
        !
        ! we use the true space vector to find corresponding reciprocal space
        ! index, so here we should add back the vector pointing to the defect
        rp_ = rp(:,i) + def_pos(:)
        !
        ! find the corresponding 
        rp_red(i,1) = rp_(1)*bg(1,1) + rp_(2)*bg(2,1) + rp_(3)*bg(3,1)
        rp_red(i,2) = rp_(1)*bg(1,2) + rp_(2)*bg(2,2) + rp_(3)*bg(3,2)
        rp_red(i,3) = rp_(1)*bg(1,3) + rp_(2)*bg(2,3) + rp_(3)*bg(3,3)
        !
        do while (rp_red(i,1) >= 1 - 1d-5) 
           rp_red(i,1) = rp_red(i,1) - 1.d0
        enddo
        do while (rp_red(i,1) <  0 - 1d-5)
           rp_red(i,1) = rp_red(i,1) + 1.d0
        enddo
        do while (rp_red(i,2) >= 1 - 1d-5) 
           rp_red(i,2) = rp_red(i,2) - 1.d0
        enddo
        do while (rp_red(i,2) <  0 - 1d-5)
           rp_red(i,2) = rp_red(i,2) + 1.d0
        enddo
        do while (rp_red(i,3) >= 1 - 1d-5) 
           rp_red(i,3) = rp_red(i,3) - 1.d0
        enddo
        do while (rp_red(i,3) <  0 - 1d-5)
           rp_red(i,3) = rp_red(i,3) + 1.d0
        enddo
        !
        rp_redm(1) = rp_red(i,1) * dffts%nr1
        rp_redm(2) = rp_red(i,2) * dffts%nr2
        rp_redm(3) = rp_red(i,3) * dffts%nr3

        ! check
        if ((abs(rp_redm(1)-nint(rp_redm(1)))>1d-4) .or. &
            (abs(rp_redm(2)-nint(rp_redm(2)))>1d-4) .or. &
            (abs(rp_redm(3)-nint(rp_redm(3)))>1d-4)) then
           !
           write(stdout,*) ' ERROR: Impurity potential mesh not consistent with FFT mesh, stop!'
           write(stdout,*) '      @', i
           write(stdout,*) '  rp_ =', rp_
           write(stdout,*) ' rp_red =', rp_red(i,1:3)
           write(stdout,*) ' rp_redm=', rp_redm(1:3)
           if (mpime == ionode_id) CALL errore('imp_meshmap', 'Inconsistency',1)
           !
        endif
        
        ! mapping
        imp_meshmap(i,1) = nint(rp_redm(1)) + 1
        imp_meshmap(i,2) = nint(rp_redm(2)) + 1
        imp_meshmap(i,3) = nint(rp_redm(3)) + 1
        imp_meshmap(i,4) = imp_meshmap(i,1) + (imp_meshmap(i,2)-1)*dffts%nr2 + &
                                              (imp_meshmap(i,3)-1)*dffts%nr2*dffts%nr3
        !
     enddo
     !
     if (mpime .eq. ionode_id) then
        !
        ! output mesh mapping
        !
        open(unit = 1004, file = 'imp_mesh_map.dat')
        do i = 1, nfftmesh
           write(1004,'(6(f15.7,1x),2x,3(i6,1x),i10)') rp(1:3,i), rp_red(i,1:3), imp_meshmap(i,1:4)
        enddo
        close(1004)

!        stop
        ! 
     endif
     !
     deallocate (rp_red)
     !
     ! we correct the space vector by adding back the vector pointing to the
     ! defect
     !
     allocate (rp_c(3,nfftmesh))
     do imesh = 1, nfftmesh
        rp_c(:,imesh) = rp(:,imesh) + def_pos(:)
     enddo
     !
  endif
  !
  write(*,*) ' pool#', my_pool_id, ' get here1'
  ! ========================================================
  ! set up exp(i*G*r) for interpolating wave function matrix
  ! ========================================================
  !
  if (eimp_mode == 5 .or. eimp_mode == 6 .or. eimp_mode == 7 .or. eimp_mode == 8) then
     !
     write(*,*) ' pool#', my_pool_id, ' get here1.0'
     bg_max = 5
     ng_max = (2*bg_max+1)**3
     bg_max_red = bg_max - 1
     ng_max_red = (2*bg_max_red + 1)**3

     write(*,*) ' pool#', my_pool_id, ' get here1.1'
     allocate (dvimp_q (nq1*nq2*nq3,2*bg_max+1,2*bg_max+1,2*bg_max+1))
     write(*,*) ' pool#', my_pool_id, ' get here1.2'

     ! for real space integration
!     allocate (dvGr (dffts%nr1*dffts%nr2*dffts%nr3,2*bg_max+1,2*bg_max+1,2*bg_max+1))        
!     write(*,*) ' pool#', my_pool_id, ' get here1.3'
!     ind_r = 0
!     do irz = 1, dffts%nr3
!        do iry = 1, dffts%nr2
!           do irx = 1, dffts%nr1
!              !
!              ind_r = ind_r + 1
!              r0 = (dble(irx-1)/dffts%nr1)*at(1:3,1) + (dble(iry-1)/dffts%nr2)*at(1:3,2) + &
!                   (dble(irz-1)/dffts%nr3)*at(1:3,3)
!
!              ind_g = 0
!              do igz = -bg_max, bg_max
!                 do igy = -bg_max, bg_max
!                    do igx = -bg_max, bg_max
!                     
!                       ind_g = ind_g + 1  
!                       G0 = igx*bg(1:3,1) + igy*bg(1:3,2) + igz*bg(1:3,3)
!                       dotGr = dot_product(G0,r0)
!                       dvGr(ind_r,igx+bg_max+1,igy+bg_max+1,igz+bg_max+1) = exp(ci*twopi*dotGr)
!
!                    enddo
!                 enddo
!              enddo
!              !
!           enddo
!        enddo
!     enddo
     !
     ! set up mapping between local G vectors and global index
     ! for G vectors, check [Modules/recvec_subs.f90]

     write(stdout,*) ' ngm = ', ngm
     allocate( map_l2g(ngm), map_g2l(dffts%nnr) )

     do ig = 1, ngm
        i1 = nint(dot_product(g(:,ig),at(:,1))) + dffts%nr1/2
        j1 = nint(dot_product(g(:,ig),at(:,2))) + dffts%nr2/2
        k1 = nint(dot_product(g(:,ig),at(:,3))) + dffts%nr3/2
        ig_global = k1*dffts%nr2*dffts%nr1 + j1*dffts%nr1 + i1
        map_l2g(ig) = ig_global
        map_g2l(ig_global) = ig
     enddo
     !
  endif
  !
  ! check spin parameters
  write(stdout,*) 'nspin = ', nspin
  write(stdout,*) 'npol = ', npol
  write(stdout,*) 'nspin_mag = ', nspin_mag
  write(stdout,*) 'domag = ', domag

#ifdef __PARA
   CALL mp_barrier(inter_pool_comm)
#endif
  !
  write(*,*) ' pool#', my_pool_id, ' get here2'
  !
  IF ((.not. epbread .and. .not. epwread) .or. &
      (.not. eimpbread .and. (eimp_mode == 5 .or. eimp_mode == 6)) .or. &
      (.not. dvimpq .and. (eimp_mode == 7 .or. eimp_mode  == 8)) ) THEN
    !
    nqc = 0
    !
    finished_iq_irr = 0
    if (epbrestore) then
       inquire(file = 'restore.epb', exist=exst)
       if (exst) then
          open(11110, file = 'restore.epb')
          read(11110,*) finished_iq_irr, nqc
          close(11110)
          !
          tempfile = trim(tmp_dir) // trim(prefix) // '.epb_temp' 
#ifdef __PARA
          CALL set_ndnmbr (0,my_pool_id+1,1,npool,filelab)
          tempfile = trim(tmp_dir) // trim(prefix) // '.epb_temp' // filelab
#endif
          OPEN  (11111, file = tempfile, form = 'unformatted')
          WRITE(stdout,'(/5x,"Reading temp epmatq on .epb_temp files"/)') 
          READ  (11111) epmatq
          CLOSE (11111)
          !
          tempfile = trim(tmp_dir) // trim(prefix) // '.eimpb_temp' 
#ifdef __PARA
          CALL set_ndnmbr (0,my_pool_id+1,1,npool,filelab)
          tempfile = trim(tmp_dir) // trim(prefix) // '.eimpb_temp' // filelab
#endif
          OPEN  (11112, file = tempfile, form = 'unformatted')
          WRITE(stdout,'(/5x,"Reading temp eimpmatq on .eimpb_temp files"/)') 
          READ  (11112) eimpmatq
          CLOSE (11112)
          !
       endif
    endif
    ! 
    ! In the loop over irr q-point, we need to read the pattern that
    ! correspond to the dvscf file computed with QE 5.
    !
    iq_first = 1
    DO iq_irr = 1, nqc_irr
      u_from_file = .TRUE.
      !tmp_dir_ph = './_ph0/'
      !
      call cpu_time(t1)
      !
      !  read the displacement patterns
      !
      IF (u_from_file) THEN
         ierr=0
         IF ( ionode ) THEN
         !
         ! ... look for an empty unit (only ionode needs it)
         !
            CALL iotk_free_unit( iunpun, ierr )
         !
         END IF
         dirname = TRIM(dvscf_dir) // TRIM( prefix ) // '.phsave'
         filename= TRIM( dirname ) // '/patterns.' // &
                   TRIM(int_to_char(iq_irr)) // '.xml'
         INQUIRE( FILE=TRIM(filename), EXIST=exst )
         IF (.NOT.exst) CALL errore( 'ph_restart_set_filename ', &
                   'cannot open file for reading or writing', ierr )
         CALL iotk_open_read( iunpun, FILE = TRIM( filename ), &
                                          BINARY = .FALSE., IERR = ierr )
         CALL read_modes(iunpun,iq_irr, ierr )
         IF (ierr /= 0) CALL errore('epw_setup', 'problem with modes file',1)
         IF (ionode) CALL iotk_close_read( iunpun )
      ENDIF
      !  
      WRITE(stdout,'(//5x,a)') repeat('=',67) 
      WRITE(stdout,'(5x,"irreducible q point # ",i4)') iq_irr
      WRITE(stdout,'(5x,a/)') repeat('=',67) 
      CALL flush(6)
      !
      xq = xqc_irr(:,iq_irr)
      !
      ! SP : The following is largely inspiered by PH/q2qstar.f90
      ! 
      ! ~~~~~~~~ setup small group of q symmetry ~~~~~~~~ 
      !
      minus_q = .true.
      sym = .false.
      sym(1:nsym) = .true.
      CALL smallg_q(xq, 0, at, bg, nsym, s, ftau, sym, minus_q) ! s is intent(in)
      !
      ! SP: Notice that the function copy_sym reshuffle the s matrix for each irr_q.  
      !     This is why we then need to call gmap_sym for each irr_q [see below]. 
      nsymq = copy_sym(nsym, sym)    
      !
      ! Recompute the inverses as the order of sym.ops. has changed
      CALL inverse_s ( )       
      call s_axis_to_cart ( )
      !
      ! This computes gi, gimq
      call set_giq (xq,s,nsymq,nsym,irotmq,minus_q,gi,gimq)
      WRITE(stdout, '(5x,a,i3)') "Symmetries of small group of q:", nsymq
      IF(minus_q) WRITE(stdout, '(10x,a)') "in addition sym. q -> -q+G:"
      ! 
      ! Finally this does some of the above again and also computes rtau...
      CALL sgam_ph_new(at, bg, nsym, s, irt, tau, rtau, nat)
      !
      ! ######################### star of q #########################
      ! 
      sym_smallq(:) = 0
      CALL star_q2(xq, at, bg, nsym, s, invs, nq, sxq, isq, imq, .true., sym_smallq )
      !
      ! The reason for xq instead of xq0 in the above is because xq is pass to QE through module  
      xq0 = xq

      ! here we jump if restoring previous calculation
      if (epbrestore .and. iq_irr <= finished_iq_irr) &
         cycle
      !
!      if (iq_irr == 3) exit

      !
      !  determine the G vector map S(G) -> G 
      !  SP: The mapping need to be done for each irr_q because the QE 5 symmetry routine
      !      reshuffle the s matrix for each irr_q [putting the sym of the small group of q first].
      !
      !  [I checked that gmapsym(gmapsym(ig,isym),invs(isym)) = ig]
      CALL gmap_sym ( nsym, s, ftau, gmapsym, eigv, invs)
      !
      !  Re-set the variables needed for the pattern representation
      !  and the symmetries of the small group of irr-q
      !  (from phq_setup.f90)
      !
      DO isym = 1, nsym
        sym (isym) = .true.
      ENDDO
      !
      CALL sgam_ph_new (at, bg, nsym, s, irt, tau, rtau, nat)
      !
#ifdef __PARA
     IF ( .not. allocated(sumr) ) allocate ( sumr(2,3,nat,3) )
     IF (mpime.eq.ionode_id) then
#endif
        CALL readmat_shuffle2 ( iq_irr, nqc_irr, nq, iq_first, sxq, imq,isq,&
                              invs, s, irt, rtau)
#ifdef __PARA
     ENDIF
     CALL mp_barrier(inter_pool_comm)
     CALL mp_barrier(intra_pool_comm)
     CALL mp_bcast (zstar, ionode_id, inter_pool_comm)
     CALL mp_bcast (zstar, root_pool, intra_pool_comm)
     CALL mp_bcast (epsi, ionode_id, inter_pool_comm)
     CALL mp_bcast (epsi, root_pool, intra_pool_comm)
     CALL mp_bcast (dynq, ionode_id, inter_pool_comm)
     CALL mp_bcast (dynq, root_pool, intra_pool_comm)
     CALL mp_bcast (sumr, ionode_id, inter_pool_comm)
     CALL mp_bcast (sumr, root_pool, intra_pool_comm)
#endif
      !
      ! now dynq is the cartesian dyn mat (NOT divided by the masses)
      !
      minus_q = (iswitch .gt. -3)  

      call cpu_time(t2)
      t12 = 0.d0
      t23 = 0.d0
      t34 = 0.d0
      t45 = 0.d0
      t56 = 0.d0
      t67 = 0.d0
      t78 = 0.d0
      t89 = 0.d0
      !
      !  loop over the q points of the star 
      !
      DO iq = 1, nq
        !
        call cpu_time(t3_1)

        ! SP: First the vlocq need to be initialized propertly with the first
        !     q in the star
        xq = xq0         
        CALL epw_init(.false.)
        !
        ! retrieve the q in the star
        xq = sxq(:,iq)                               
        !
        ! and populate the uniform grid
        nqc = nqc + 1
        xqc(:,nqc) = xq    ! this is on cartesian coordinate
        !
        IF (iq.eq.1) write(6,*)
        WRITE( stdout, 5) nqc, xq
        !
        !  prepare the gmap for the refolding
        !
        CALL createkmap ( xq )                      
        !
        IF (iverbosity.eq.1) then
          !
          !   description of symmetries
          !
          WRITE( stdout, '(36x,"s",24x,"frac. trans.")')
          CALL s_axis_to_cart() ! give sr (:,:, isym)
          DO isym = 1, nsym
            WRITE( stdout, '(/6x,"isym = ",i2,5x,a45/)') isym, sname(isym)
          !  CALL s_axis_to_cart (s(1,1,isym), sr, at, bg)
            IF (ftau(1,isym).ne.0.or.ftau(2,isym).ne.0.or.ftau(3,isym).ne.0) then
                ft1 = at(1,1)*ftau(1,isym)/dfftp%nr1 + at(1,2)*ftau(2,isym)/dfftp%nr2 + &
                      at(1,3)*ftau(3,isym)/dfftp%nr3
                ft2 = at(2,1)*ftau(1,isym)/dfftp%nr1 + at(2,2)*ftau(2,isym)/dfftp%nr2 + &
                      at(2,3)*ftau(3,isym)/dfftp%nr3
                ft3 = at(3,1)*ftau(1,isym)/dfftp%nr1 + at(3,2)*ftau(2,isym)/dfftp%nr2 + &
                      at(3,3)*ftau(3,isym)/dfftp%nr3
                WRITE( stdout, '(1x,"cryst.",3x,"s(",i2,") = (",3(i6,5x), &
                      &        " )    f =( ",f10.7," )")') &
                      isym, (s(1,ipol,isym),ipol=1,3), dble(ftau(1,isym))/dble(dfftp%nr1)
                WRITE( stdout, '(17x," (",3(i6,5x), " )       ( ",f10.7," )")') &
                            (s(2,ipol,isym),ipol=1,3), dble(ftau(2,isym))/dble(dfftp%nr2)
                WRITE( stdout, '(17x," (",3(i6,5x), " )       ( ",f10.7," )"/)') &
                            (s(3,ipol,isym),ipol=1,3), dble(ftau(3,isym))/dble(dfftp%nr3)
                WRITE( stdout, '(1x,"cart. ",3x,"s(",i2,") = (",3f11.7, &
                      &        " )    f =( ",f10.7," )")') &
                      isym, (sr(1,ipol,isym),ipol=1,3), ft1
                WRITE( stdout, '(17x," (",3f11.7, " )       ( ",f10.7," )")') &
                            (sr(2,ipol,isym),ipol=1,3), ft2
                WRITE( stdout, '(17x," (",3f11.7, " )       ( ",f10.7," )"/)') &
                            (sr(3,ipol,isym),ipol=1,3), ft3
            ELSE
                WRITE( stdout, '(1x,"cryst.",3x,"s(",i2,") = (",3(i6,5x), " )")') &
                                       isym,  (s (1, ipol, isym) , ipol = 1,3)
                WRITE( stdout, '(17x," (",3(i6,5x)," )")')  (s(2,ipol,isym), ipol=1,3)
                WRITE( stdout, '(17x," (",3(i6,5x)," )"/)') (s(3,ipol,isym), ipol=1,3)
                WRITE( stdout, '(1x,"cart. ",3x,"s(",i2,") = (",3f11.7," )")') &
                                              isym,  (sr (1, ipol,isym) , ipol = 1, 3)
                WRITE( stdout, '(17x," (",3f11.7," )")')  (sr (2, ipol,isym) , ipol = 1, 3)
                WRITE( stdout, '(17x," (",3f11.7," )"/)') (sr (3, ipol,isym) , ipol = 1, 3)
            ENDIF
          ENDDO
          !
        ENDIF
        !
        ! isq(isym)=iq means: when we apply symmetry isym to the originating q 
        ! of the star, we get the iq-th member of the star. There are as many 
        ! matches as the degeneracy of the star.
        !
        ! We now need to pick up the q in the small group of q* so that Sxq0+G=iq with G=0.
        ! If we choose another element in the small group
        ! the actual q-point may be Sq+G and we screw up the q-vector below to generate
        ! k+q from k and for the KB projectors
        !
        nsq = 0 ! nsq is the degeneracy of the small group for this iq in the star
        !
        DO jsym = 1, nsym
          IF ( isq(jsym) .eq. iq ) then
             nsq = nsq + 1
             sym_sgq(nsq) = jsym
          ENDIF
        ENDDO
        IF ( nsq*nq .ne. nsym ) CALL errore ('elphon_shuffle_wrap', 'wrong degeneracy', iq)
        ! 
        IF (iverbosity.eq.1) then
          !
          WRITE(stdout,*) 'iq, i, isym, nog, symmo'
          DO i = 1, nsq
            !
            isym = sym_sgq(i)
            ism1 = invs (isym)
            !
            !  check for G such that Sq = q* + G 
            ! 
            aq  = xq0
            saq = xq
            CALL cryst_to_cart (1, aq, at, -1)
            DO j = 1, 3
              raq (j) = s (j, 1, ism1) * aq (1) &
                      + s (j, 2, ism1) * aq (2) &
                      + s (j, 3, ism1) * aq (3)
            ENDDO
            CALL cryst_to_cart (1, saq, at, -1)
            nog = eqvect_strict( raq, saq) 
            !
            !  check whether the symmetry belongs to a symmorphic group
            !
            symmo = (ftau(1,isym).eq.0 .and. ftau(2,isym).eq.0 .and. ftau(3,isym).eq.0)
            !
            WRITE(stdout,'(3i5,2a)') iq, i, isym, nog, symmo
            !
          ENDDO  
          !
        ENDIF
        ! 
        ! SP: We now need to select one symmetry among the small group of q (i.e. that respect 
        !     Sq0+G=q ) that has G=0. There should always be such symmetry. 
        !     We enforce this for later easiness. 
        ! 
        aq = xq0
        saq = xq
        call cryst_to_cart (1, aq, at, - 1)
        CALL cryst_to_cart (1, saq, at, -1)
        !write(*,*)'xq0 ',aq
        !write(*,*)'xq ',saq
        DO jsym=1, nsq
          ism1 = invs (sym_sgq(jsym))
          raq = 0.d0
          DO ipol = 1, 3
             DO jpol = 1, 3
                raq (ipol) = raq (ipol) + s (ipol, jpol, ism1) * aq (jpol)
             ENDDO
          ENDDO
          nog = eqvect_strict( raq, saq)
          IF (nog) THEN ! This is the symmetry such that Sq=q
            isym = sym_sgq(jsym)
            EXIT
          ENDIF
          ! If we enter into that loop it means that we have not found 
          ! such symmetry within the small group of Q. 
          IF (jsym == nsq) THEN
            CALL errore( 'elphon_shuffle_wrap ', 'No sym. such that Sxq0=iq was found in the sgq !', 1 )
          ENDIF
        ENDDO

        !
        !DBSP
        !write(*,*)'isym ',isym 
        !END
        !
        !
        CALL loadumat ( nbnd, nbndsub, nks, nkstot, xq, cu, cuq, lwin, lwinq )

        call cpu_time(t3_2)
        !
        ! Calculate overlap U_k+q U_k^\dagger
        IF (lpolar) CALL compute_bmn_para3 ( nbnd, nbndsub, nks, cu, cuq,bmat(:,:,:,nqc) )
        !
        call cpu_time(t3_3)

        !   calculate the sandwiches
        !
        ! a more accurate way of doing is to symmetrize the matrix element w.r.t.
        ! the small group of the given q in the star. I'm not doint this here.
        ! (but I checked that even without symm the result of full zone and irr zone
        ! are equal to 5+ digits)
        ! For any volunteers, please write to giustino@civet.berkeley.edu
        !
!==============================================================

        if (.not. epbjump) &
           CALL elphon_shuffle ( iq_irr, nqc_irr, nqc, gmapsym, eigv, isym, invs, xq0, .false. )

        call cpu_time(t3_4)
!==============================================================
   
        ! calculate electron-impurity matrix
        if (eimp_mode == 5 .or. eimp_mode == 6) then
           !
           CALL eimpmat_shuffle ( iq_irr, nqc_irr, nqc, gmapsym, eigv, isym, invs, xq0, .false. )
           !
        endif
        !
        ! Fourier component of short range defect potential

        if (.not. dvimpq .and. (eimp_mode == 7 .or. eimp_mode == 8)) then

           xq_cart(1) = xq(1) !- NINT(xq(1))
           xq_cart(2) = xq(2) !- NINT(xq(2))
           xq_cart(3) = xq(3) !- NINT(xq(3))
           ! xq is already in cartesian coordinate
!           CALL cryst_to_cart ( 1, xq_cart, bg, 1 )   ! xq_cart in [2*pi/alat]

           do igz = -bg_max, bg_max
              ind_gz = igz+bg_max+1
              do igy = -bg_max, bg_max
                 ind_gy = igy+bg_max+1
                 do igx = -bg_max, bg_max
                    ind_gx = igx+bg_max+1

                    G0 = igx*bg(1:3,1) + igy*bg(1:3,2) + igz*bg(1:3,3)
                    xqG = xq_cart + G0

                    do imesh = 1, nfftmesh
                       dvimp_q(nqc,ind_gx,ind_gy,ind_gz) = dvimp_q(nqc,ind_gx,ind_gy,ind_gz) + &
                                            dvimp(imesh) * exp(-ci*twopi*dot_product(xqG,rp_c(:,imesh)))
                    enddo
                 enddo
              enddo
           enddo
           !
           ! normalization: dvimp(q+G) = 1/Omega_uc * integral dr^3 dvimp(r)*exp(-i*(q+G)*r)
           ! for discrete sum, dr^3/Omega_uc becomes sc_un_ratio / nfftmesh
           dvimp_q(nqc,:,:,:) = dvimp_q(nqc,:,:,:) * (sc_un_ratio/nfft_sc)

        endif

        call cpu_time(t3_5)
        !
        !  bring epmatq in the mode representation of iq_first, 
        !  and then in the cartesian representation of iq
        !
        CALL rotate_eigenm ( iq_first, iq, nqc, isym, nsym, s, invs, irt, &
           rtau, xq, isq, cz1, cz2 )

        call cpu_time(t3_6)
        !
        CALL rotate_epmat ( cz1, cz2, xq, nqc, lwin, lwinq )

        call cpu_time(t3_7)

  !DBSP
  !      write(*,*)'epmatq(:,:,2,:,nqc)',SUM(epmatq(:,:,2,:,nqc))
  !      write(*,*)'epmatq(:,:,2,:,nqc)**2',SUM((REAL(REAL(epmatq(:,:,2,:,nqc))))**2)+&
  !        SUM((REAL(AIMAG(epmatq(:,:,2,:,nqc))))**2)
  !END
        ! SP: Now we treat separately the case imq == 0
        IF (imq .eq. 0) then
          !
          ! SP: First the vlocq need to be initialized propertly with the first
          !     q in the star
          xq = -xq0
          CALL epw_init(.false.)
          !
          ! retrieve the q in the star
          xq = -sxq(:,iq)
          !
          ! and populate the uniform grid
          nqc = nqc + 1
          xqc(:,nqc) = xq
          !
          IF (iq.eq.1) write(stdout,*)
          WRITE( stdout, 5) nqc, xq
          !
          !  prepare the gmap for the refolding
          !
          CALL createkmap ( xq )
          !
          xq0 = -xq0
          !
          if (.not. epbjump) &
             CALL elphon_shuffle ( iq_irr, nqc_irr, nqc, gmapsym, eigv, isym, invs, xq0, .true. )
          !
          ! calculate electron-impurity matrix
          if (eimp_mode == 5 .or. eimp_mode == 6) then
             !
             CALL eimpmat_shuffle ( iq_irr, nqc_irr, nqc, gmapsym, eigv, isym, invs, xq0, .true. )
             !
          endif
          !
          ! Fourier component of short range defect potential

          if (.not. dvimpq .and. (eimp_mode == 7 .or. eimp_mode == 8)) then

             xq_cart(1) = xq(1) !- NINT(xq(1))
             xq_cart(2) = xq(2) !- NINT(xq(2))
             xq_cart(3) = xq(3) !- NINT(xq(3))
             ! NOTE: xq is already in cartesian coordinate
!             CALL cryst_to_cart ( 1, xq_cart, bg, 1 )   ! xq_cart in [2*pi/alat]

             do igz = -bg_max, bg_max
                ind_gz = igz+bg_max+1
                do igy = -bg_max, bg_max
                   ind_gy = igy+bg_max+1
                   do igx = -bg_max, bg_max
                      ind_gx = igx+bg_max+1

                      G0 = igx*bg(1:3,1) + igy*bg(1:3,2) + igz*bg(1:3,3)
                      xqG = xq_cart + G0

                      do imesh = 1, nfftmesh
                         dvimp_q(nqc,ind_gx,ind_gy,ind_gz) = dvimp_q(nqc,ind_gx,ind_gy,ind_gz) + &
                                              dvimp(imesh) * exp(-ci*twopi*dot_product(xqG,rp_c(:,imesh)))
                      enddo
                   enddo
                enddo
             enddo
             !
             ! normalization: dvimp(q+G) = 1/Omega_uc * integral dr^3 dvimp(r)*exp(-i*(q+G)*r)
             ! for discrete sum, dr^3/Omega_uc becomes sc_un_ratio / nfftmesh
             dvimp_q(nqc,:,:,:) = dvimp_q(nqc,:,:,:) * (sc_un_ratio/nfft_sc)

          endif
          !
          !  bring epmatq in the mode representation of iq_first, 
          !  and then in the cartesian representation of iq
          !
          CALL rotate_eigenm ( iq_first, iq, nqc, isym, nsym, s, invs, irt, &
             rtau, xq, isq, cz1, cz2 )
          !
          CALL rotate_epmat ( cz1, cz2, xq, nqc, lwin, lwinq )
          !
    !DBSP
    !      write(*,*)'epmatq(:,:,2,:,nqc)',SUM(epmatq(:,:,2,:,nqc))
    !      write(*,*)'epmatq(:,:,2,:,nqc)**2',SUM((REAL(REAL(epmatq(:,:,2,:,nqc))))**2)+&
    !        SUM((REAL(AIMAG(epmatq(:,:,2,:,nqc))))**2)
    !END
          xq0 = -xq0
        ENDIF ! end imq == 0  
        !
        call cpu_time(t3_8)
 
        t23 = t23 + (t3_2 - t3_1)
        t34 = t34 + (t3_3 - t3_2)
        t45 = t45 + (t3_4 - t3_3)
        t56 = t56 + (t3_5 - t3_4)
        t67 = t67 + (t3_6 - t3_5)
        t78 = t78 + (t3_7 - t3_6)
        t89 = t89 + (t3_8 - t3_7)
        !
      ENDDO ! iq loop
      !
      t12 = (t2 - t1)

      write(stdout,*) '-----------------------------'
      write(stdout,*) ' usage time, t12:', t12
      write(stdout,*) ' usage time, t23:', t23
      write(stdout,*) ' usage time, t34:', t34
      write(stdout,*) ' usage time, t45:', t45
      write(stdout,*) ' usage time, t56:', t56
      write(stdout,*) ' usage time, t67:', t67
      write(stdout,*) ' usage time, t78:', t78
      write(stdout,*) ' usage time, t89:', t89
      write(stdout,*) '-----------------------------'
      !
      iq_first = iq_first + nq
      if (imq .eq. 0) iq_first = iq_first + nq
      !
#ifdef __PARA
      CALL mp_barrier(inter_pool_comm)
#endif
      !
      if (epbrestore) then
         !
         open (11110, file = 'restore.epb', status = 'replace')
         write(11110,*) iq_irr, nqc
         close(11110)
         !
         tempfile = trim(tmp_dir) // trim(prefix) // '.epb_temp'
#ifdef __PARA
         CALL set_ndnmbr (0,my_pool_id+1,1,npool,filelab)
         tempfile = trim(tmp_dir) // trim(prefix) // '.epb_temp' // filelab
#endif
         OPEN  (11111, file = tempfile, form = 'unformatted', status = 'replace')
         WRITE(stdout,'(/5x,"Writing temp epmatq on .epb_temp files"/)') 
         WRITE (11111) epmatq
         CLOSE (11111)
         !
         tempfile = trim(tmp_dir) // trim(prefix) // '.eimpb_temp' 
#ifdef __PARA
         CALL set_ndnmbr (0,my_pool_id+1,1,npool,filelab)
         tempfile = trim(tmp_dir) // trim(prefix) // '.eimpb_temp' // filelab
#endif
         OPEN  (11112, file = tempfile, form = 'unformatted', status = 'replace')
         WRITE(stdout,'(/5x,"Writing temp eimpmatq on .eimpb_temp files"/)') 
         WRITE (11112) eimpmatq
         CLOSE (11112)
         !
      endif
      !
    ENDDO ! irr-q loop
    ! 
    ! write dvimp_q
    !
    if (.not.dvimpq .and. (eimp_mode == 7) .or. (eimp_mode == 8)) then
       !
       ! output short-range potential
       !
       tempfile = 'dvimp_q.' // trim(adjustl(defectname)) // '.dat'
       open(unit = 1003, file = tempfile)
       do iq = 1, nqc
          do igz = 1, 2*bg_max+1
             do igy = 1, 2*bg_max+1
                do igx = 1, 2*bg_max+1
                   write(1003,'(2(es20.12,1x))') dble(dvimp_q(iq,igx,igy,igz)), dimag(dvimp_q(iq,igx,igy,igz))
                enddo
             enddo
          enddo
       enddo
       close(1003)
       !
    endif
    !
    IF (nqc.ne.nq1*nq2*nq3) &
       CALL errore ('elphon_shuffle_wrap','nqc .ne. nq1*nq2*nq3',nqc)
    wqlist = dble(1)/dble(nqc)
    !
  ENDIF
  !
  !
  IF ( epbread .or. epbwrite ) THEN
    !
    ! write the e-ph matrix elements and other info in the Bloch representation
    ! (coarse mesh)
    ! in .epb files (one for each pool)
    !
    tempfile = trim(tmp_dir) // trim(prefix) // '.epb' 
#ifdef __PARA
     CALL set_ndnmbr (0,my_pool_id+1,1,npool,filelab)
     tempfile = trim(tmp_dir) // trim(prefix) // '.epb' // filelab
#endif
     !
     IF (epbread)  THEN
        inquire(file = tempfile, exist=exst)
        IF (.not. exst ) CALL errore( 'elphon_shuffle_wrap', 'epb files not found ', 1)
        OPEN  (iuepb, file = tempfile, form = 'unformatted')
        WRITE(stdout,'(/5x,"Reading epmatq from .epb files"/)') 
        READ  (iuepb) nqc, xqc, et, dynq, epmatq, zstar, epsi
        CLOSE (iuepb)
     ENDIF
     !
     IF (epbwrite) THEN
        OPEN  (iuepb, file = tempfile, form = 'unformatted')
        WRITE(stdout,'(/5x,"Writing epmatq on .epb files"/)') 
        WRITE (iuepb) nqc, xqc, et, dynq, epmatq, zstar, epsi
        CLOSE (iuepb)
     ENDIF
  ENDIF
  !
  !
  IF ( (eimpbread .or. eimpbwrite) .and. (eimp_mode == 5 .or. eimp_mode == 6) ) THEN
    !
    ! write the e-imp matrix elements in the Bloch representation
    ! (coarse mesh)
    ! in .epb files (one for each pool)
    !
    tempfile = trim(tmp_dir) // trim(prefix) // '.' // trim(adjustl(defectname)) // '.eimpb' 
#ifdef __PARA
     CALL set_ndnmbr (0,my_pool_id+1,1,npool,filelab)
     tempfile = trim(tmp_dir) // trim(prefix) // '.' // trim(adjustl(defectname)) // '.eimpb' // filelab
#endif
     !
     IF (eimpbread)  THEN
        inquire(file = tempfile, exist=exst)
        IF (.not. exst ) CALL errore( 'elphon_shuffle_wrap', 'eimpb files not found ', 1)
        OPEN  (iueimpb, file = tempfile, form = 'unformatted')
        WRITE(stdout,'(/5x,"Reading eimpmatq from .eimpb files"/)') 
        READ  (iueimpb) eimpmatq
        CLOSE (iueimpb)
     ENDIF
     !
     IF (eimpbwrite) THEN
        OPEN  (iueimpb, file = tempfile, form = 'unformatted')
        WRITE(stdout,'(/5x,"Writing eimpmatq on .eimpb files"/)') 
        WRITE (iueimpb) eimpmatq
        CLOSE (iueimpb)
     ENDIF
  ENDIF

  !
  if (eimp_mode == 7 .or. eimp_mode == 8) then
     !
     ! read dvimp_q (not centered)
     !
     if (dvimpq) then
#ifdef __PARA
     IF (mpime .EQ. ionode_id) THEN
#ENDIF
        tempfile = 'dvimp_q.' // trim(adjustl(defectname)) // '.dat'
        open(unit = 1003, file = tempfile)
        do iq = 1, nq1*nq2*nq3
           do igz = 1, 2*bg_max+1
              do igy = 1, 2*bg_max+1
                 do igx = 1, 2*bg_max+1
                    read(1003,*) dvimp_r, dvimp_i
                    dvimp_q(iq,igx,igy,igz) = dcmplx(dvimp_r,dvimp_i)
                 enddo
              enddo
           enddo
        enddo
        close(1003)
#ifdef __PARA
     ENDIF
     CALL mp_bcast (dvimp_q, ionode_id, inter_pool_comm)
     CALL mp_bcast (dvimp_q, root_pool, intra_pool_comm)
#ENDIF
     endif

     ! generate centered dvimp_qc
     !
     allocate (dvimp_qc (nq1*nq2*nq3,2*bg_max+1,2*bg_max+1,2*bg_max+1))
     dvimp_qc = 0.d0

     do iq = 1, nq1*nq2*nq3
        xq = xqc(:,iq)   ! xqc is in cartesian coordinate

        ! bring xq into crystal coordinate
        CALL cryst_to_cart (1, xq, at, -1)

        xq1 = xq(1)
        iG1 = 0
        do while (xq1 < 0)
           xq1 = xq1 + 1
           iG1 = iG1 - 1
        enddo
        do while (xq1 >= 1)
           xq1 = xq1 - 1
           iG1 = iG1 + 1
        enddo

        xq2 = xq(2)
        iG2 = 0
        do while (xq2 < 0)
           xq2 = xq2 + 1
           iG2 = iG2 - 1
        enddo
        do while (xq2 >= 1)
           xq2 = xq2 - 1
           iG2 = iG2 + 1
        enddo

        xq3 = xq(3)
        iG3 = 0
        do while (xq3 < 0)
           xq3 = xq3 + 1
           iG3 = iG3 - 1
        enddo
        do while (xq3 >= 1)
           xq3 = xq3 - 1
           iG3 = iG3 + 1
        enddo

        iq1 = nint(xq1*nq1)
        iq2 = nint(xq2*nq2)
        iq3 = nint(xq3*nq3)
        nnq = iq1*nq2*nq3 + iq2*nq3 + iq3 + 1

        do igz = -bg_max, bg_max
           ind_gz = igz+bg_max+1
           do igy = -bg_max, bg_max
              ind_gy = igy+bg_max+1
              do igx = -bg_max, bg_max
                 ind_gx = igx+bg_max+1

                 if ((ind_gx+iG1>0) .and. (ind_gy+iG2>0) .and. (ind_gz+iG3>0) .and. &
                     (ind_gx+iG1<=2*bg_max+1) .and. (ind_gy+iG2<=2*bg_max+1) .and. (ind_gz+iG3<=2*bg_max+1)) then
                    dvimp_qc(nnq,ind_gx+iG1,ind_gy+iG2,ind_gz+iG3) = dvimp_q(iq,ind_gx,ind_gy,ind_gz)

                    if (.false. .and. (igx==0) .and. (igy==0) .and. (igz==0)) then
                       write(stdout,*)
                       write(stdout,*) ' dvimp_qc test, iq#: ', iq
                       write(stdout,*) ' xq =', xq
                       write(stdout,*) ' iq1..3, nnq', iq1, iq2, iq3, nnq
                       write(stdout,*) ' ind_g', ind_gx, ind_gy, ind_gz
                       write(stdout,*) ' dvimp_qc:', dvimp_qc(nnq,ind_gx+iG1,ind_gy+iG2,ind_gz+iG3)
                       write(stdout,*) ' dvimp_q:', dvimp_q(iq,ind_gx,ind_gy,ind_gz)
                    endif
                 endif
              enddo
           enddo
        enddo
        !
     enddo

     do iq = 1, nq1*nq2*nq3
        if (dvimp_qc(iq,bg_max+1,bg_max+1,bg_max+1) == 0.d0) then
           write(stdout,*) ' zero at iq =', iq
           stop
        endif
     enddo

     !
     ! follow procedure in [elphel_shuffle2] to prepare the wavefunctions
     !
     ! close all sequential files in order to re-open them as direct access
     ! close all .wfc files in order to prepare shuffled read

     CLOSE (unit = iunigk, status = 'keep')
     CLOSE (unit = iuwfc,  status = 'keep')
     ipooltmp = my_pool_id+1
     !
     ! find largest global index of G across pools
     ig_max = 0
     DO ik = 1, nks
        CALL readigk (ipooltmp, ik, npw, igk)
        DO ig = 1, npw
           if (igk(ig)>ig_max) ig_max = igk(ig)
        ENDDO
        write(stdout,*) ' npw =', npw
     ENDDO
     write(*,*) ' before mp_max, pool# ',my_pool_id, 'ig_max =',ig_max
     CALL mp_max(ig_max,inter_pool_comm)
     write(*,*) ' after mp_max, pool# ',my_pool_id, 'ig_max =',ig_max
     !
     ! the first index (G) of evc_k is the global index of G
     if (.not.(allocated(evc_k))) allocate (evc_k(nbnd,npol,nks,ig_max))
     evc_k = (0.d0, 0.d0)
     !
     DO ik = 1, nks
        norm0 = 0
        CALL readwfc (ipooltmp, ik, evc)
        CALL readigk (ipooltmp, ik, npw, igk)

        CALL ktokpmq ( xk(:,ik),  zero_vect, +1, ipool, nkk, nkk_abs )
        IF ( .not. ALLOCATED (umat) )  ALLOCATE ( umat(nbnd,nbnd,nks) )
        umat(:,:,ik)  = umat_all(:,:, nkk_abs)
        !
        DO ig = 1, npw
           !
           DO ibnd = 1, nbnd
              evc_k (ibnd, 1, ik, igk(ig)) = evc (ig, ibnd)
              norm0 = norm0 + (abs(evc(ig,ibnd)))**2
              if (noncolin) &
                 evc_k (ibnd, 2, ik, igk(ig)) = evc (ig+npwx, ibnd)
           ENDDO
           if (nls(igk(ig))>ig_nrmax) ig_nrmax = nls(igk(ig))
           ! 
           ! gauge settings on wavefunctions
           CALL zgemm ('n', 'n', 1, nbnd, nbnd, cone, evc_k(:,1,ik,igk(ig)), & 
                       1, umat(:,:,ik), nbnd, czero, evctmp, 1)
           evc_k(:,1,ik,igk(ig)) = evctmp(:)
           if (noncolin) then
              CALL zgemm ('n', 'n', 1, nbnd, nbnd, cone, evc_k(:,2,ik,igk(ig)), & 
                          1, umat(:,:,ik), nbnd, czero, evctmp, 1)
              evc_k(:,2,ik,igk(ig)) = evctmp(:)
           endif
        ENDDO

!        write(stdout,*) ' normalization check for ik =',ik
!        write(stdout,*) '     xk    = ',xk(1:3,ik)
!        write(stdout,*) ' ----nbnd = ',nbnd
!        write(stdout,*) ' ----norm0 = ',norm0
     ENDDO

     CALL mp_max(ig_nrmax,inter_pool_comm)
     write(stdout,*) '   ig_max =',ig_max
     write(stdout,*) ' ig_nrmax =',ig_nrmax
     !
     !  restore original configuration of files
     !
     CALL seqopn (iunigk, 'igk', 'unformatted', exst)
     CALL diropn (iuwfc, 'wfc', lrwfc, exst)
     !
     ! set up mapping between local and global G points
     ind_g = 0
     allocate (map_kq2k (ig_max,ng_max_red))
     map_kq2k = 0
     do igz = -bg_max_red, bg_max_red
        do igy = -bg_max_red, bg_max_red
           do igx = -bg_max_red, bg_max_red
              ind_g = ind_g + 1
              ig_shift = igz*dffts%nr2*dffts%nr1 + igy*dffts%nr1 + igx

              do ig = 1, ig_max
                 if (((map_l2g(ig)-ig_shift) > 0) .and. ((map_l2g(ig)-ig_shift) < dffts%nnr)) then
                    if (map_g2l(map_l2g(ig)-ig_shift) <= ig_max) then
                       map_kq2k(ig,ind_g) = map_g2l(map_l2g(ig)-ig_shift)
                    endif
                 endif
              enddo
           enddo
        enddo
     enddo
     !
     if (.false.) then
        !
        ALLOCATE ( eemat_check( nbnd, nbnd, ng_max_red) )
        ALLOCATE ( evck ( ig_max, nbnd ) )
        ALLOCATE ( evck_remap ( ig_max ) )
        !
        xk_gamma = (/0.d0, 0.d0, 0.d0/)
        xk_x     = (/0.d0, 1.d0, 0.d0/)
        xk_l     = (/-0.5d0,0.5d0,-0.5d0/)

        do ik = 1, nks

           gamma_find = .false.
           if ((abs(xk(1,ik)-xk_gamma(1))<1d-5) .and. (abs(xk(2,ik)-xk_gamma(2))<1d-5) .and. &
               (abs(xk(3,ik)-xk_gamma(3))<1d-5)) then
              gamma_find = .true.
           endif

           l_find = .false.
           if ((abs(xk(1,ik)-xk_l(1))<1d-5) .and. (abs(xk(2,ik)-xk_l(2))<1d-5) .and. &
               (abs(xk(3,ik)-xk_l(3))<1d-5)) then
              l_find = .true.
              write(*,*) ' find L =', xk(:,ik), ' @ pool#', my_pool_id
           endif

           eemat_check = 0
   
           do ibnd = 1, nbnd
              !
              ind_g = 0
              do igz = -bg_max_red, bg_max_red
                 do igy = -bg_max_red, bg_max_red
                    do igx = -bg_max_red, bg_max_red
                       !
                       ind_g = ind_g + 1
                       ig_shift = igz*dffts%nr2*dffts%nr1 + igy*dffts%nr1 + igx
                       !
                       evck_remap = 0.d0
                       do ig = 1, ig_max
                          if ((map_kq2k(ig,ind_g) >= 1) .and. (map_kq2k(ig,ind_g) <= ig_max)) &
                             evck_remap(ig) = evc_k(ibnd,1,ik,map_kq2k(ig,ind_g))
                          evck(ig,1:nbnd) = evc_k(1:nbnd,1,ik,ig)
                       enddo
                       !
                       do jbnd = 1, nbnd
                          eemat_check(jbnd,ibnd,ind_g) = &   ! use reduced number of G points
                                 ZDOTC (ig_max, evck(1,jbnd), 1, evck_remap(1), 1)
                       enddo
                       !
                    enddo
                 enddo
              enddo
           enddo

           CALL mp_sum(eemat_check, intra_pool_comm)
           if (gamma_find) then

              ibnd = 10
              write(*,*) ' check_evc in elphon_shuffle_wrap:'
              write(*,*) ' xk =', xk(:,ik)
              write(*,*) ' evc(1..5):'
              write(*,*) 1, evc_k(ibnd,1,ik,1)
              write(*,*) 2, evc_k(ibnd,1,ik,2)
              write(*,*) 3, evc_k(ibnd,1,ik,3)
              write(*,*) 4, evc_k(ibnd,1,ik,4)
              write(*,*) 5, evc_k(ibnd,1,ik,5)

              ind_g = 0
              do igz = -bg_max_red, bg_max_red
                 do igy = -bg_max_red, bg_max_red
                    do igx = -bg_max_red, bg_max_red
                       ind_g = ind_g + 1
                       if ((igz==0) .and. (igy==0)) then
                          write(*,*) '    eemat_check =', igz,igy,igx, eemat_check(ibnd,ibnd,ind_g)
                       endif
                    enddo 
                 enddo
              enddo
           endif

        enddo
 
     endif
     !
  endif
  !
  !
#ifdef __PARA
  CALL mp_barrier(inter_pool_comm)
#endif
  !
  ! remove temporary restoring files

  IF ((.not. epbread .and. .not. epwread) .or. &
      (.not. eimpbread .and. (eimp_mode == 5 .or. eimp_mode == 6))) THEN
     !
     if (epbrestore .and. .true.) then
        !
        if (my_pool_id == ionode_id) then
           open (11110, file = 'restore.epb', status = 'old')
           close(11110, status = 'delete')
        endif
        !
        tempfile = trim(tmp_dir) // trim(prefix) // '.epb_temp'
#ifdef __PARA
        CALL set_ndnmbr (0,my_pool_id+1,1,npool,filelab)
        tempfile = trim(tmp_dir) // trim(prefix) // '.epb_temp' // filelab
#endif
        OPEN  (11111, file = tempfile, status = 'old')
        CLOSE (11111, status = 'delete')
        !
        tempfile = trim(tmp_dir) // trim(prefix) // '.eimpb_temp' 
#ifdef __PARA
        CALL set_ndnmbr (0,my_pool_id+1,1,npool,filelab)
        tempfile = trim(tmp_dir) // trim(prefix) // '.eimpb_temp' // filelab
#endif
        OPEN  (11112, file = tempfile, status = 'old')
        CLOSE (11112, status = 'delete')
        ! 
     endif 
     !
  ENDIF 
  !
  !
  ! test e-imp matrix element on coarse mesh
  !
  if (.false. .and. (eimp_mode == 5 .or. eimp_mode == 6)) then
     !
     xk_gamma = (/0.d0, 0.d0, 0.d0/)
     xk_x     = (/0.d0, 1.d0, 0.d0/)
     !
     do ik = 1, nks
        !
        gamma_find = .false.
        if ((abs(xk(1,ik)-xk_gamma(1))<1d-5) .and. (abs(xk(2,ik)-xk_gamma(2))<1d-5) .and. &
            (abs(xk(3,ik)-xk_gamma(3))<1d-5)) then
           gamma_find = .true.
        endif
        !
        x_find = .false.
        if ((abs(xk(1,ik)-xk_x(1))<1d-5) .and. (abs(xk(2,ik)-xk_x(2))<1d-5) .and. &
            (abs(xk(3,ik)-xk_x(3))<1d-5)) then
           x_find = .true.
        endif
        !
        if (gamma_find) then
           ibnd = 9
           write(*,*) ' eimp matrix at gamma with k2 = k1:'
           write(*,*) ' gamma_find =', xk(:,ik), ' @ pool#', my_pool_id
           write(*,*) '                      ', eimpmatq (ibnd, ibnd, ik, 1)
        endif
        !
        if (x_find) then
           ibnd = 10
           write(*,*) ' eimp matrix at X with k2 = k1:'
           write(*,*) ' x_find =', xk(:,ik), ' @ pool#', my_pool_id
           write(*,*) '                      ', eimpmatq (ibnd, ibnd, ik, 1)
        endif
     enddo
     !
  endif
  !
  !
  IF ( .not.epbread .and. epwread ) THEN
  !  CV: need dummy nqc, xqc for the ephwann_shuffle call
     nqc=1
     xqc=0.d0
     WRITE(stdout,'(/5x,"Do not need to read .epb files; read .fmt files"/)')
  !
  ENDIF
  !
  !
#ifdef __PARA
  CALL mp_barrier(inter_pool_comm)
#endif

!  stop
  !
  !   now dynq is the cartesian dyn mat ( NOT divided by the masses)
  !   and epmatq is the epmat in cartesian representation (rotation in elphon_shuffle)
  !
  !
  ! free up some memory
  !
  IF ( ASSOCIATED (evq)  ) NULLIFY    (evq)
  IF ( ALLOCATED  (evc)  ) DEALLOCATE (evc)
  IF ( ASSOCIATED (igkq) ) NULLIFY    (igkq)
  IF ( ALLOCATED  (igk)  ) DEALLOCATE (igk)
  IF ( ALLOCATED  (dvpsi)) DEALLOCATE (dvpsi)
  IF ( ALLOCATED  (dpsi) ) DEALLOCATE (dpsi)
  IF ( ALLOCATED  (sumr) ) DEALLOCATE (sumr)
  IF ( ALLOCATED (cu) )    DEALLOCATE (cu)
  IF ( ALLOCATED (cuq) )   DEALLOCATE (cuq)
  IF ( ALLOCATED (lwin) )  DEALLOCATE (lwin)
  IF ( ALLOCATED (lwinq) ) DEALLOCATE (lwinq)
  IF ( ALLOCATED (bmat) )  DEALLOCATE (bmat)
  ! 
  CALL stop_clock ( 'elphon_wrap' )
!DBSP
!  DO iq = 1, nqc
!    write(*,*) iq, xqc(:,iq)
!    write(*,*)'epmatq(:,:,2,:,iq)',SUM(epmatq(:,:,2,:,iq))
!    write(*,*)'epmatq(:,:,2,:,iq)**2',SUM((REAL(REAL(epmatq(:,:,2,:,iq))))**2)+&
!               SUM((REAL(AIMAG(epmatq(:,:,2,:,iq))))**2)
!  ENDDO
!END
  !
  ! the electron-phonon wannier interpolation
  !
  CALL ephwann_shuffle ( nqc, xqc )
  !
5 format (8x,"q(",i5," ) = (",3f12.7," )") 
  !
  !
  !-------------------------------------------------------------------
  ! THL: BTE calculation
  !-------------------------------------------------------------------
96969 CONTINUE
  !
#ifdef __PARA
  CALL mp_barrier (inter_pool_comm)
#ENDIF
  !
  WRITE (stdout,'(//5x,a)') '==================================================================='
  WRITE (stdout,'(5x,a)')   '                        Electron-phonon BTE                        '
  WRITE (stdout,'(5x,a/)')  '==================================================================='
  !
  ! start BTE
  IF (bte .EQ. 0 .OR. bte .EQ. -1) THEN
     !
     CALL driver_bte_rta ()
     !
  ELSEIF (bte .EQ. 1 .OR. bte .EQ. 10) THEN
     !
     CALL driver_bte_el ()
     !
  ELSEIF (bte .EQ. 3 .OR. bte .EQ. 30) THEN
     !
     CALL driver_tdbte_el ()
     !
  ELSEIF (bte .EQ. 19 .OR. bte .EQ. 18) THEN
     !
     CALL driver_bte_random ()
     !
  ENDIF
  !
  IF (ionode .AND. phdrag .AND. bte .EQ. 1) THEN
     CALL SYSTEM ('rm -r BTE/SIGMAI/ BTE/GAMMAI/')
  ELSE
     CALL SYSTEM ('rm -r BTE/SIGMAI/')
  ENDIF
  !
#ifdef __PARA
  CALL mp_barrier (inter_pool_comm)
#ENDIF
  !
  !-------------------------------------------------------------------
  !
  !
  RETURN
  END SUBROUTINE elphon_shuffle_wrap
  !
  !---------------------------------------------------------------------------
  SUBROUTINE irotate ( x, s, sx)
  !---------------------------------------------------------------------------
  !
  ! a simple symmetry operation in crystal coordinates ( s is integer!)
  !
  USE kinds, ONLY : DP
  implicit none
  real(kind=DP) x(3), sx(3)
  integer :: s(3,3),i
  !
  DO i = 1, 3
     sx (i) = dble(s (i, 1)) * x (1) &
            + dble(s (i, 2)) * x (2) &
            + dble(s (i, 3)) * x (3)
  ENDDO
  !
  RETURN
  end
  !---------------------------------------------------------------------------
  logical function eqvect_strict (x, y)
  !-----------------------------------------------------------------------
  !
  !   This function test if two tridimensional vectors are equal
  !
  USE kinds
  implicit none
  real(kind=DP) :: x (3), y (3)
  ! input: input vector
  ! input: second input vector
  real(kind=DP) :: accep
  ! acceptance PARAMETER
  PARAMETER (accep = 1.0d-5)
  !
  eqvect_strict = abs( x(1)-y(1) ).lt.accep .and. &
                  abs( x(2)-y(2) ).lt.accep .and. &
                  abs( x(3)-y(3) ).lt.accep
  END FUNCTION eqvect_strict
  !---------------------------------------------------------------------------
  SUBROUTINE read_modes(iunpun,current_iq, ierr )
  !---------------------------------------------------------------------------
  !
  ! This routine reads the displacement patterns.
  !
  USE modes,        ONLY : nirr, npert, u, name_rap_mode, num_rap_mode
  USE lr_symm_base, ONLY : minus_q, nsymq  
  USE iotk_module,  ONLY : iotk_index, iotk_scan_dat, iotk_scan_begin, &
                           iotk_scan_end
  USE io_global,    ONLY : ionode, ionode_id
  USE mp_images,    ONLY : intra_image_comm
  USE mp,           ONLY : mp_bcast

  IMPLICIT NONE

  INTEGER,          INTENT(IN) :: current_iq,iunpun
  INTEGER,          INTENT(OUT) :: ierr
  INTEGER :: imode0, imode, irr, ipert, iq 
  !
  ierr=0
  IF (ionode) THEN
     CALL iotk_scan_begin( iunpun, "IRREPS_INFO" )
     !
     CALL iotk_scan_dat(iunpun,"QPOINT_NUMBER",iq)
  ENDIF
  CALL mp_bcast( iq,  ionode_id, intra_image_comm )
  IF (iq /= current_iq) CALL errore('read_modes', &
            'problems with current_iq', 1 )

  IF (ionode) THEN

     CALL iotk_scan_dat(iunpun, "QPOINT_GROUP_RANK", nsymq)
     CALL iotk_scan_dat(iunpun, "MINUS_Q_SYM", minus_q)
     CALL iotk_scan_dat(iunpun, "NUMBER_IRR_REP", nirr)
     imode0=0
     DO irr=1,nirr
        CALL iotk_scan_begin( iunpun, "REPRESENTION"// &
                                   TRIM( iotk_index( irr ) ) )
        CALL iotk_scan_dat(iunpun,"NUMBER_OF_PERTURBATIONS", npert(irr))
        DO ipert=1,npert(irr)
           imode=imode0+ipert
           CALL iotk_scan_begin( iunpun, "PERTURBATION"// &
                                   TRIM( iotk_index( ipert ) ) )
           CALL iotk_scan_dat(iunpun,"SYMMETRY_TYPE_CODE", &
                                                      num_rap_mode(imode))
           CALL iotk_scan_dat(iunpun,"SYMMETRY_TYPE", name_rap_mode(imode))
           CALL iotk_scan_dat(iunpun,"DISPLACEMENT_PATTERN",u(:,imode))
           CALL iotk_scan_end( iunpun, "PERTURBATION"// &
                                   TRIM( iotk_index( ipert ) ) )
        ENDDO
        imode0=imode0+npert(irr)
        CALL iotk_scan_end( iunpun, "REPRESENTION"// &
                                   TRIM( iotk_index( irr ) ) )
     ENDDO
     !
     CALL iotk_scan_end( iunpun, "IRREPS_INFO" )
     !
  ENDIF

  CALL mp_bcast( nirr,  ionode_id, intra_image_comm )
  CALL mp_bcast( npert,  ionode_id, intra_image_comm )
  CALL mp_bcast( nsymq,  ionode_id, intra_image_comm )
  CALL mp_bcast( minus_q,  ionode_id, intra_image_comm )
  CALL mp_bcast( u,  ionode_id, intra_image_comm )
  CALL mp_bcast( name_rap_mode,  ionode_id, intra_image_comm )
  CALL mp_bcast( num_rap_mode,  ionode_id, intra_image_comm )

  RETURN
  END SUBROUTINE read_modes
