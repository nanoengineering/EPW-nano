  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !----------------------------------------------------------------------
SUBROUTINE ephwann_shuffle( nqc, xqc )
  !---------------------------------------------------------------------
  !
  !  Wannier interpolation of electron-phonon vertex
  !
  !  Scalar implementation   Feb 2006
  !  Parallel version        May 2006
  !  Disentenglement         Oct 2006
  !  Compact formalism       Dec 2006
  !  Phonon irreducible zone Mar 2007
  !
  !  No ultrasoft now
  !  No spin polarization
  !
  !  RM - add noncolin case
  !-----------------------------------------------------------------------
  !
  !
#include "f_defs.h"
  USE kinds,         ONLY : DP
  USE pwcom,         ONLY : nbnd, nks, nkstot, isk, &
                            et, xk, at, bg, ef, alat, nelec, omega
  USE spin_orb,      ONLY : lspinorb
  USE start_k,       ONLY : nk1, nk2, nk3
  USE ions_base,     ONLY : amass, ityp, nat
  USE phcom,         ONLY : nq1, nq2, nq3, nmodes, w2
  USE control_lr,    ONLY : lgamma
  USE epwcom,        ONLY : nbndsub, lrepmatf, fsthick, lretf, epwread,   &
                            epwwrite, ngaussw, degaussw, lpolar,          &
                            nbndskip, parallel_k, parallel_q, etf_mem,    &
                            nest_fn, a2f, indabs, &
                            epexst, vme, eig_read, ephwrite, mp_mesh_k,   & 
                            efermi_read, fermi_energy, specfun, band_plot, neptemp, &
                            ! THL
                            nkf1, nkf2, nkf3, nqf1, nqf2, nqf3, bte, epdim, &
                            edos_read, asr_eph, smearing, ph_read, ephl_read, &
                            epdope, nepdope, nptype, save_m_mat, save_m_matw, &
                            save_m_ph, save_t_el, vg_el, vg_ph, epthick, phdrag, &
                            shengbte_read, alloy_read, alloy_pot, &
                            eimp_mode, dielec, defectname, &
                            epcheck, nbnd_c, nbnd_v, xk_c, xk_v, filqf, &
                            dvimpsr, elop, L_D, imp_charge, eimp_ls_mode, &
                            screen_polar, eph_interp
  USE noncollin_module, ONLY : noncolin, npol
  USE constants_epw, ONLY : ryd2ev, ryd2mev, rydcm1, one, two, twopi, pi, au2cm, &
                            au2ps, ryd2thz, ci, bohr2ang
  USE control_flags, ONLY : iverbosity
  USE io_files,      ONLY : tmp_dir, prefix, diropn
  USE io_global,     ONLY : stdout, ionode
  USE io_epw,        ONLY : lambda_phself,linewidth_phself,linewidth_elself, iunepmatf, &
                            iuetf, iunepmatwe, iunepmatwp, iunepmatwp_asr
  USE elph2,         ONLY : nrr_k, nrr_q, cu, cuq, lwin, lwinq, irvec, ndegen_k, ndegen_q, &
                            wslen, chw, chw_ks, cvmew, cdmew, rdw, epmatwp, epmatwe, epmatq, &
                            wf, etf, etf_k, etf_ks, xqf, xkf, wkf, wqf, &
                            dynq, nqtotf, nkqf, epf17, nkf, nqf, et_ks, &
                            ibndmin, ibndmax, lambda_all, dmec, dmef, vmef, &
                            gamma_all, nkqtotf, epsi, zstar, efnew, &
                            ! THL
                            nbnd_red, cbnd_emin, vbnd_emax, ef_m, delta_egap, outside_gap, cfsthick, vfsthick, &
                            vbnd_num, &
                            sigmai_mode_all_abs, sigmai_mode_all_emi, sigmai_mode_all_alloy_inter, &
                            sigmai_mode_all_alloy_intra, &
                            sigmai_mode_all_intra, sigmai_mode_all_inter, &
                            sigmai_mode_all_ela_intra, sigmai_mode_all_ela_inter, &
                            etf_all, vel_all, uf_ful, wf_ful, vph_ful, &
                            wmax, epmatwp_asr, eph_vogl, &
                            gammai_mode_all, wf_irr, uf_irr, vph_irr, vkq, &
                            model_charge, alat_imp, scella, scellb, scellc, &
                            nfftmesh, rp, dvimp, eimpf_lr, eimpf_sr, alelpp, &
                            cbnd_emin_nxk, cbnd_emin_xk, def_pos, &
                            ! JWZ, for impurity scattering
                            eimpmatq, eimpmatwp, eimpmatwe, evc_k, evwan, &
                            evc_fmesh, map_kq2k, ig_max, bg_max, ng_max, &
                            map_l2g, map_g2l, dvimp_q, eemat_f, eimpf_full, &
                            bg_max_red, ng_max_red
  USE fft_base,      ONLY : dffts
  USE bte_var
  USE bte_func
  USE epw_explore
  USE tetrahedron
#ifdef __NAG
  USE f90_unix_io,   ONLY : flush
  USE,INTRINSIC :: f90_unix_file, ONLY:fstat, stat_t
#endif
#ifdef __PARA
  USE mp,            ONLY : mp_barrier, mp_bcast, mp_sum
  USE io_global,     ONLY : ionode_id
  USE mp_global,     ONLY : my_pool_id, nproc_pool, intra_image_comm, &
                            inter_pool_comm, me_pool, root_pool, intra_pool_comm, npool
  USE mp_world,      ONLY : mpime
#endif
  !
  implicit none
  !
#ifdef __NAG
  TYPE(stat_t) :: statb
#endif
#ifndef __NAG
  integer :: fstat,statb(13)
#endif
  !
  complex(kind=DP), ALLOCATABLE :: &
    epmatwe_q (:,:,:,:),     &
    epmatwef (:,:,:,:)           ! e-p matrix  in el wannier - fine Bloch phonon grid
  complex(kind=DP), ALLOCATABLE :: &
    epmatf( :, :, :, :, :),           &! e-p matrix  in smooth Bloch basis, fine mesh
    cufkk ( :, :),              &! Rotation matrix, fine mesh, points k
    cufkq ( :, :),              &! the same, for points k+q
    uf    ( :, :,:,:),              &! Rotation matrix for phonons
    bmatf ( :, :)                ! overlap U_k+q U_k^\dagger in smooth Bloch basis, fine mesh
  integer :: &
    nqc                          ! number of qpoints in the coarse grid
  real(kind=DP) :: &
    xqc (3, nqc)                 ! qpoint list, coarse mesh
  !
  integer :: ir, na, &
             fermicount, nrec, indnew, indold, lrepmatw, ios, irq, lrepmatw_asr, &
             nk_band, ik_band
  LOGICAL :: already_skipped, exst
  character (len=256) :: filint, tempfile
  character (len=30)  :: myfmt
  real(kind=DP) :: xxq(3), xxk(3), xkk(3), xkq(3), size_m, xxq_(3), rp_(3), &
                   temp(3), scell_vol, dvol, xq_cart(3), n_imp, dvimp_, &
                   xkk2(3), xxq_int(3)
  real(kind=DP), allocatable :: vcoul_model(:), eig_k(:), xxq_band(:,:)
  real(kind=DP), external :: efermig
  real(kind=DP), external :: efermig_seq
  real(kind=DP), external :: V_coulomb
  real(kind=DP), parameter :: eps = 0.01/ryd2mev
  ! 
  ! THL
  ! parallelization
  INTEGER                       :: ik_star, ik_stop, iq_star, iq_stop, &
                                   ifft_star, ifft_stop, imesh
  ! local variable
  INTEGER                       :: ik, ik0, ik_red, ik0_red, iq, iq0, idx1, idx2, idx3
  INTEGER                       :: ibnd, jbnd, ibnd0, jbnd0, mbnd
  INTEGER                       :: imode, mu, nu, i, j, k
  INTEGER                       :: ikk, ikq, iqq, ik_irr
  integer                       :: nfftx_sc, nffty_sc, nfftz_sc, nfft_sc
  REAL(KIND=DP)                 :: etf_ks_tmp(nbndsub), wmax_mode(nmodes), vel_tmp(3,nbndsub), vph_tmp(3,nmodes)
  COMPLEX(KIND=DP)              :: dmef_tmp(3,nbndsub,nbndsub)
  ! group velocity
  REAL(KIND=DP)                 :: velf(3,nbndsub), vq_(3,nmodes)
  ! BTE
  INTEGER                       :: nscat
  CHARACTER (LEN=3)             :: file_num
  INTEGER                       :: ifile
  ! polar interaction
  COMPLEX(KIND=DP)              :: uf_cart(nmodes,nmodes), ephl_vogl(nmodes,neptemp,nepdope)
  CHARACTER (LEN=256)           :: vogl_ufmt
  CHARACTER (LEN=3)             :: vogl_cpu
  INTEGER                       :: iq_vogl, vogl_num(npool)
  ! check within_range
  LOGICAL                       :: within_range
  ! save memory
  CHARACTER (LEN=256)           :: epf17_ufmt 
  CHARACTER (LEN=3)             :: epf17_cpu
  REAL(KIND=DP), allocatable    :: wf_dyn(:,:,:), vph_dyn(:,:,:,:), eig_q(:,:,:,:)
  COMPLEX(KIND=DP), allocatable :: uf_dyn(:,:,:,:), uf_(:,:,:,:)
  CHARACTER (LEN=256)           :: dyn_ufmt, wf_ufmt, vph_ufmt, uf_ufmt
  CHARACTER (LEN=3)             :: dyn_cpu
  INTEGER                       :: iq_dyn, dyn_num(npool)
  ! save memory (Wannier)
  CHARACTER (LEN=256)           :: matwe_ufmt
  ! save time
  CHARACTER (LEN=256)           :: ham_ufmt, etf_ufmt, cuf_ufmt
  CHARACTER (LEN=3)             :: ham_cpu
  INTEGER                       :: icpu, ik_ham, ham_num(npool)
  INTEGER                       :: ijk_k(3), ijk_q(3), ijk_kq(3), ikq_ham
  REAL(KIND=DP)                 :: etf_ham(nbndsub), etf_ks_ham(nbndsub)
  COMPLEX(KIND=DP)              :: cufkk_ham(nbndsub,nbndsub)
  ! random k-point calculation
  CHARACTER (LEN=256)           :: filename_check
  ! random k-point calculation
  CHARACTER (LEN=256)           :: filename_kmesh, filename_qmesh, coord_mesh
  ! phdrag
  CHARACTER (LEN=256)           :: filename_phd
  LOGICAL                       :: file_exist
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
                                   tc=0.0d0, tci, tcf, t9_1, t9_2, t9_3, t9_4, t9_5, t9_6
  ! impurity scattering
  integer :: ibnd_c, jbnd_c, ibnd_v, jbnd_v, center_index, ig, ind_g, &
             igx, igy, igz, ig_shift, &
             itemp, idope, ig_max_red, ikq_red, ik_red_, ikq_ful
  integer, allocatable :: map_mesh_ful(:), mapcheck_ful(:), &
                          map_mesh_rful(:), mapcheck_rful(:)
  REAL(KIND=DP)                 :: dielec0, rpl, rp_max, dv_tot_far, vcoul_far, &
                                   xql, xqG(3), G_cart(3), xqGl, norm_evcfk, &
                                   imp_charge_, xkq_(3), xkk_(3)
  complex(kind=DP), ALLOCATABLE :: eimpmatwe_q (:,:,:), eimpmatwef (:,:,:), &
                                   eimpmatf(:,:), evc_kw(:,:,:), dvimpqf_sr(:), &
                                   dvimpqf_lr(:,:), eimpmat_f(:,:,:,:), &
                                   evc_ks(:,:,:), evc_fk(:,:,:), evc_fkq(:,:,:)
  complex(kind=DP) :: eimpf_sr_, alelpp_, evc_remap(ig_max,npol), eimpmat_r(nbnd,nbnd)
  logical :: gamma_find, x_find, cond_find, vale_find, evcfkq_found, l_find
  character(len=256) :: filename
  COMPLEX(DP),EXTERNAL :: ZDOTC
  !
  !
  IF (nbndsub.ne.nbnd) &
       WRITE(stdout, '(/,14x,a,i4)' ) 'band disentanglement is used:  nbndsub = ', nbndsub
  !
  ALLOCATE ( cu ( nbnd, nbndsub, nks), & 
             cuq ( nbnd, nbndsub, nks), & 
             lwin ( nbnd, nks ), &
             lwinq ( nbnd, nks ), &
             irvec (3, 20*nk1*nk2*nk3), &
             ndegen_k (20*nk1*nk2*nk3), &
             ndegen_q (20*nq1*nq2*nq3), &
             wslen(20*nk1*nk2*nk3)      )
  !
  CALL start_clock ( 'ephwann' )
  !
  ! DBSP
  ! HERE loadkmesh
  !
  ! determine Wigner-Seitz points
  !
  CALL wigner_seitz2 (nk1, nk2, nk3, nq1, nq2, nq3, nrr_k, nrr_q, irvec, wslen, ndegen_k, ndegen_q)
  !
  ! 
  ! At this point, we will interpolate the Wannier rep to the Bloch rep 
  !
  IF ( epwread ) THEN
     !
     !  read all quantities in Wannier representation from file
     !  in parallel case all pools read the same file
     !
     CALL epw_read
     !
  ELSE !if not epwread (i.e. need to calculate fmt file)
     !
     xxq = 0.d0 
     CALL loadumat (nbnd, nbndsub, nks, nkstot, xxq, cu, cuq, lwin, lwinq)  
     !
     ! ------------------------------------------------------
     !   Bloch to Wannier transform
     ! ------------------------------------------------------
     !
     ALLOCATE (chw(nbndsub, nbndsub, nrr_k))
     ALLOCATE (chw_ks(nbndsub, nbndsub, nrr_k))
     ALLOCATE (cdmew(3, nbndsub, nbndsub, nrr_k))
     ALLOCATE (rdw(nmodes, nmodes, nrr_q))
     IF (vme) ALLOCATE(cvmew(3, nbndsub, nbndsub, nrr_k))
     IF (asr_eph) ALLOCATE (epmatwp_asr ( nbndsub, nbndsub, nrr_k, 3))
     !
     ALLOCATE(epmatwe_q(nbndsub, nbndsub, nrr_k, nmodes))
     IF (.NOT. save_m_matw) THEN
        ALLOCATE (epmatwe(nbndsub, nbndsub, nrr_k, nmodes, nqc)) 
        ALLOCATE (epmatwp(nbndsub, nbndsub, nrr_k, nmodes, nrr_q))
        epmatwe = (0.0d0,0.0d0)
        epmatwp = (0.0d0,0.0d0)
     ENDIF
     !
     if (eimp_mode == 5 .or. eimp_mode == 6) then
        ALLOCATE(eimpmatwe_q(nbndsub, nbndsub, nrr_k))
        IF (.NOT. save_m_matw) THEN
           ALLOCATE (eimpmatwe(nbndsub, nbndsub, nrr_k, nqc)) 
           ALLOCATE (eimpmatwp(nbndsub, nbndsub, nrr_k, nrr_q))
           eimpmatwe = (0.0d0,0.0d0)
           eimpmatwp = (0.0d0,0.0d0)
        ENDIF
     endif
     ! 
     if (eimp_mode == 7 .or. eimp_mode == 8) then
        ALLOCATE (evwan (nbndsub, npol, nrr_k , ig_max))
        ALLOCATE (evc_kw(nks, npol, ig_max))
        evwan = (0.d0,0.d0)
     endif
     !
     WRITE (stdout,'(/5x,a)')        'Check the use of memory:'  
     WRITE (stdout,'(/8x,a,i8)')     'nbndsub = ', nbndsub
     WRITE (stdout,'(8x,a,i8)')      'nmodes  = ', nmodes
     WRITE (stdout,'(8x,a,i8)')      'nrr_k   = ', nrr_k
     WRITE (stdout,'(8x,a,i8)')      'nqc     = ', nqc
     WRITE (stdout,'(8x,a,i8)')      'nrr_q   = ', nrr_q
     WRITE (stdout,'(8x,a,f9.1,a)')  'epmatwe ~ ', (DBLE(16*nbndsub*nbndsub*nmodes*nqc)*1.0d-6)*DBLE(nrr_k), ' MB'
     WRITE (stdout,'(8x,a,f9.1,a/)') 'epmatwp ~ ', (DBLE(16*nbndsub*nbndsub*nmodes*nrr_q)*1.0d-6)*DBLE(nrr_k), ' MB'
     !
     !
     ! Hamiltonian
     !
     CALL hambloch2wan &
          ( nbnd, nbndsub, nks, nkstot, lgamma, et, xk, cu, lwin, nrr_k, irvec, wslen, chw )
     !
     ! Kohn-Sham eigenvalues
     !
     IF (eig_read) THEN
       WRITE (6,'(5x,a)') "Interpolating MB and KS eigenvalues"
       CALL hambloch2wan (nbnd, nbndsub, nks, nkstot, lgamma, et_ks, xk, cu, lwin, nrr_k, irvec, wslen, chw_ks)
     ENDIF
     !
     ! Dipole
     !
     CALL dmebloch2wan (nbnd, nbndsub, nks, nkstot, dmec, xk, cu, nrr_k, irvec, wslen)
     !
     ! Dynamical Matrix 
     !
     CALL dynbloch2wan (nmodes, nqc, xqc, dynq, nrr_q, irvec, wslen)
     !
     ! Transform of position matrix elements
     ! PRB 74 195118  (2006)
     IF (vme) CALL vmebloch2wan (nbnd, nbndsub, nks, nks, nkstot, lgamma, xk, cu, nrr_k, irvec, wslen)
     !
     ! Electron-Phonon vertex (Bloch el and Bloch ph -> Wannier el and Bloch ph)
     !
     IF (save_m_matw .AND. my_pool_id .EQ. ionode_id) THEN
        !
        matwe_ufmt = TRIM(tmp_dir)//TRIM(prefix)//'.epmatwe_q'
        OPEN (81913,FILE=matwe_ufmt,FORM='unformatted',ACCESS='direct',RECL=2*nbndsub*nbndsub*nmodes*nrr_k*DP,STATUS='replace')
        !
     ENDIF
     !
     !
     if (eimp_mode == 7 .or. eimp_mode == 8) then
        !
        do ig = 1, ig_max
           !
           call evcbloch2wan (nbnd, nbndsub, nks, nkstot, xk, cu, evc_k(:,:,:,ig), &
                              nrr_k, irvec, evwan(:,:,:,ig), wslen, evc_kw(:,:,ig) )
           !
        enddo
        !
     endif
     !
     !
     DO iq = 1, nqc
       !
       xxq = xqc (:, iq)
       !
       ! we need the cu again for the k+q points, we generate the map here
       !
       CALL loadumat ( nbnd, nbndsub, nks, nkstot, xxq, cu, cuq, lwin, lwinq )
       !
       DO imode = 1, nmodes
         !
         CALL ephbloch2wane &
              (nbnd, nbndsub, nks, nkstot, lgamma, xk, cu, cuq, lwin, lwinq, epmatq(:,:,:,imode,iq), nrr_k, irvec, wslen, epmatwe_q(:,:,:,imode))
         !
       ENDDO 
       !
       if (eimp_mode == 5 .or. eimp_mode == 6) then
          !
          ! here we also transform the electron-impurity matrix element
          !
          CALL eimpbloch2wane &
               (nbnd, nbndsub, nks, nkstot, lgamma, xk, cu, cuq, lwin, lwinq, eimpmatq(:,:,:,iq), nrr_k, irvec, wslen, eimpmatwe_q(:,:,:))
          !
       endif
       !
       IF (save_m_matw) THEN
          IF (my_pool_id .EQ. ionode_id) WRITE (81913,REC=iq) epmatwe_q(1:nbndsub,1:nbndsub,1:nrr_k,1:nmodes)
       ELSE
          epmatwe(:,:,:,:,iq) = epmatwe_q(:,:,:,:)
          if (eimp_mode == 5 .or. eimp_mode == 6) &
             eimpmatwe(:,:,:,iq) = eimpmatwe_q(:,:,:)
       ENDIF
       !
     ENDDO
     !
     IF (save_m_matw .AND. my_pool_id .EQ. ionode_id) CLOSE (81913)
     !
#ifdef __PARA
     CALL mp_barrier(inter_pool_comm)
#endif
     !
     ! Electron-Phonon vertex (Wannier el and Bloch ph -> Wannier el and Wannier ph)
     ! Only master perform this task. Need to be parallelize in the future (SP)
     !
!potential ERROR
     ! Also, depending on [eimp_mode], also transform the electron-impurity matrix element
     ! Check: since only master performs this task, makes sure other nodes get
     ! correct values
     !
     IF (ionode) THEN
        !
        CALL ephbloch2wanp (nbndsub, nmodes, xqc, nqc, irvec, wslen, nrr_k, nrr_q)
        !
        if (eimp_mode == 5 .or. eimp_mode == 6) then
           !
           CALL eimpbloch2wanp (nbndsub, xqc, nqc, irvec, wslen, nrr_k, nrr_q)
           !
        endif
        !
     ENDIF
     !
#ifdef __PARA
     CALL mp_barrier(inter_pool_comm)
#endif
     !
     IF ( epwwrite ) THEN
        CALL epw_write 
        CALL epw_read
     ENDIF
     !
  ENDIF
  !
  !
  IF ( ALLOCATED (epmatwe) )     DEALLOCATE (epmatwe)
  IF ( ALLOCATED (epmatwe_q) )   DEALLOCATE (epmatwe_q)
  IF ( ALLOCATED (epmatq) )      DEALLOCATE (epmatq)
  IF ( ALLOCATED (cu) )          DEALLOCATE (cu)
  IF ( ALLOCATED (cuq) )         DEALLOCATE (cuq)
  IF ( ALLOCATED (lwin) )        DEALLOCATE (lwin)
  IF ( ALLOCATED (lwinq) )       DEALLOCATE (lwinq)
  !
  ! check
  if (epcheck .and. .false.) then
     write(stdout,*) ' check right after bloch-wan transform:'
     write(stdout,*) ' eimpmatwp at pool:', my_pool_id
     write(stdout,*) ' eimpmatwp:', eimpmatwp(10,10,1,1:5)
     if (my_pool_id == 4) then
        write(*,*) ' eimpmatwp at pool:', my_pool_id
        write(*,*) ' eimpmatwp:', eimpmatwp(10,10,1,1:5)
     endif
  endif
  !
  !----------------------------------------------------------------------------------------
  !
  WRITE (stdout,'(//5x,a)') '==================================================================='
  WRITE (stdout,'(5x,a)')   '                     Meshes and symmetry index                     '
  WRITE (stdout,'(5x,a/)')  '==================================================================='
  !
  CALL bz_index (nkf1,nkf2,nkf3,'k')
  CALL bz_index (nqf1,nqf2,nqf3,'q')
  !
  WRITE (stdout,'(/5x,a)') 'WARNING: check the file "xkf_irr_cryst/cart" is consistent with QE'
  !
  CALL kq_load ()
  !
  WRITE (stdout,'(/5x,a)') 'k and q points in full and irreducible BZ: '
  WRITE (stdout,'(/13x,a,i12)') 'nk_ful : ', nk_ful
  WRITE (stdout,'(13x,a,i12)')  'nk_irr : ', nk_irr
  WRITE (stdout,'(13x,a,i12)')  'nq_ful : ', nq_ful
  WRITE (stdout,'(13x,a,i12)')  'nq_irr : ', nq_irr
  !
  !
  !
  WRITE (stdout,'(//5x,a)') '==================================================================='
  WRITE (stdout,'(5x,a)')   '                            Fermi level                            '
  WRITE (stdout,'(5x,a/)')  '==================================================================='
  !
  ! compute energy on irr-mesh first
  CALL para_bounds (ik_star, ik_stop, nk_irr)
  !
  ALLOCATE (cufkk(nbndsub, nbndsub))
  ALLOCATE (etf(nbndsub,nk_irr))
  etf = 0.0d0
  !
  DO ik = ik_star, ik_stop
     !
     xxk(:) = xkf_irr(:,ik)
     CALL hamwan2bloch (nbndsub, nrr_k, irvec, ndegen_k, xxk, cufkk, etf(:,ik), chw)
     !
  ENDDO
  !
#ifdef __PARA
  CALL mp_barrier (inter_pool_comm)
  CALL mp_sum (etf,inter_pool_comm)
#ENDIF
  !
  ! etf will be used in fermiwindow and edos
  CALL fermiwindow ()
  !
  ! find out conduction band edge k-point; ERROR: currently only implements for
  ! n-type
  write(stdout,*)
  do ik = 1, nk_irr
     if (etf(vbnd_num+1,ik) == cbnd_emin) then
        write(stdout,'(a,1x,3(f15.7))') ' cbnd_emin: ', xkf_irr(:,ik)
        ik_irr = ik
     endif
  enddo
  !
  open (11111,FILE='BTE/META/xkf_ful_cryst',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='old') 
  !
  cbnd_emin_nxk = 0
  cbnd_emin_xk  = 0.d0
  write(stdout,*) 
  do ik = 1, nk_ful
     if (ful2irr(ik) == ik_irr) then
        read (11111, REC=ik) xkk(1:3)
        write(stdout,'(a,1x,3(f15.7))') ' degeneracy k: ', xkk(:)
        !
        cbnd_emin_nxk = cbnd_emin_nxk + 1
        cbnd_emin_xk(:,cbnd_emin_nxk) = xkk(:)  ! cbnd_emin_xk in crystal coord
     endif
  enddo

  !
  ! density of states
  IF (.NOT. edos_read) CALL edos_energy ()
  !
  ! fermi level location
  IF(.NOT. efermi_read) THEN
     CALL fermilocatioepdope ()
  ELSE
     CALL fermilocation_ef ()
  ENDIF
  !
#ifdef __PARA
  CALL mp_barrier (inter_pool_comm)
#ENDIF
  !
  !
  IF (bte .EQ. 2 .OR. bte .EQ. 29) GOTO 75302
  !
  !
  IF (bte .EQ. 19 .OR. bte .EQ. 18) THEN
     !
     WRITE (stdout,'(//5x,a)') '==================================================================='
     WRITE (stdout,'(5x,a)')   '                     Random k-point calculation                    '
     WRITE (stdout,'(5x,a/)')  '==================================================================='
     !
     WRITE (stdout,'(5x,a/)') 'WARNING: If nkf1~3 are chosen too small, ibndmin, ibndmax and ef_m may not be determined percisely'
     !
     ! check exisitence of file
     IF (bte .EQ. 19) THEN
        filename_kmesh = TRIM(prefix)//'.kmesh'
     ELSEIF (bte .EQ. 18) THEN
        filename_kmesh = TRIM(prefix)//'.channel'
     ENDIF
     INQUIRE (FILE=filename_kmesh,EXIST=file_exist)
     IF (.NOT. file_exist) CALL errore('ephwann_shuffle','cannot find .kmesh file (bte=19)',1)
     !
     OPEN (3029,FILE=filename_kmesh,STATUS='old')
     READ (3029,*) nk_kmesh, coord_mesh
     !
     ALLOCATE (xk_kmesh(3,nk_kmesh))
     xk_kmesh = 0.0d0
     !
     IF (bte .EQ. 18) THEN
        !
        ALLOCATE (band_ch(nk_kmesh))
        ALLOCATE (sigmai_ch(neptemp,nepdope,nmodes,nq_ful,nk_kmesh))
        band_ch  = 0
        sigmai_ch = 0.0d0
        !
     ENDIF
     !
     DO ik = 1, nk_kmesh
        IF (bte .EQ. 19) THEN
           READ (3029,*) xk_kmesh(1:3,ik)
        ELSEIF (bte .EQ. 18) THEN
           READ (3029,*) xk_kmesh(1:3,ik), band_ch(ik)
        ENDIF

        IF (coord_mesh .EQ. 'cart' .OR. coord_mesh .EQ. 'carte' .OR. coord_mesh .EQ. 'cartesian' .OR. &
            coord_mesh .EQ. 'Cart' .OR. coord_mesh .EQ. 'Carte' .OR. coord_mesh .EQ. 'Cartesian' .OR. &
            coord_mesh .EQ. '1') CALL cryst_to_cart (1, xk_kmesh(:,ik), at, -1)
     ENDDO 
     CLOSE (3029)
     !
     WRITE (stdout,'(5x,a,i5)') 'Number of k points in .kmesh file: ', nk_kmesh
     !
#ifdef __PARA
     CALL mp_barrier (inter_pool_comm)
#ENDIF
     !
  ENDIF
  !
  !
  !
  WRITE (stdout,'(//5x,a)') '==================================================================='
  WRITE (stdout,'(5x,a)')   '                        Select out k points                        '
  WRITE (stdout,'(5x,a/)')  '==================================================================='
  !
  IF (bte .EQ. 19 .OR. bte .EQ. 18) THEN
     !
     nk_ful = nk_kmesh
     nk_irr = nk_kmesh
     nk_irr_red = nk_kmesh
     nk_ful_red = nk_kmesh
     !
     WRITE (stdout,'(5x,a,i5)')  'Random k point, re-define nk_ful     = ', nk_kmesh
     WRITE (stdout,'(5x,a,i5)')  '                          nk_irr     = ', nk_kmesh
     WRITE (stdout,'(5x,a,i5)')  '                          nk_ful_red = ', nk_kmesh
     WRITE (stdout,'(5x,a,i5/)') '                          nk_irr_red = ', nk_kmesh
     WRITE (stdout,'(5x,a,f8.1,a/)')  'sigmai_mode ~ ', DBLE(nk_irr_red*neptemp*nepdope*nmodes*nbnd_red)*1.0d-6*(2*DP), ' MB' 
     !
  ELSE
     CALL reduce_index (etf)
     !
     WRITE (stdout,'(5x,a,f8.1,a/)')  'sigmai_mode ~ ', DBLE(nk_irr_red*neptemp*nepdope*nmodes*nbnd_red)*1.0d-6*(2*DP), ' MB' 
     !
     CALL kq_red_load ()
     !
  ENDIF
  !
  CALL mp_barrier (inter_pool_comm)
  !
  !
  !
  ! deallocate
  DEALLOCATE (etf)
  DEALLOCATE (cufkk)
  !
  !
  ! constant relaxation time approximation
  IF (bte .EQ. -1) GOTO 75301 
  !
  !
  !
  !----------------------------------------------------------------------------------------
  ! preparation of smearing method
  !----------------------------------------------------------------------------------------
  !
  ! tetrahedral smearing
  IF (smearing .EQ. 'tetra') THEN
     !
     IF (bte .EQ. 0 .OR. bte .EQ. 1 .OR. bte .EQ. 3 .OR. bte .EQ. 19 .OR. bte .EQ. 18) THEN
        !
        ntetra = 6*nqf1*nqf2*nqf3
        ALLOCATE (tetra_i(4,ntetra),tetra_c(4,ntetra),tetra_w(4,ntetra),wkt(nqf1*nqf2*nqf3),itetra(2,30,nqf1*nqf2*nqf3))
        ! need to check if in tetra generation the inherent k-point order is same 
        ! as the definition used in pw
        CALL make_kp_reg(nqf1, nqf2, nqf3, wkt, ntetra, tetra_i, itetra)
        !
     ENDIF
     !
  ENDIF
  !
  !
  !----------------------------------------------------------------------------------------
  !
  ! Wannier interpolation
  ! Only valid for electron self-energy currently and eft_mem = .true. and ephwrite = .false.
  !
  !----------------------------------------------------------------------------------------
  !
  ! inverse Wannier interpolation for electron and phonon
  WRITE (stdout,'(//5x,a)') '==================================================================='
  WRITE (stdout,'(5x,a)')   '                Inverse Wannier interpolation (e-ph)               '
  WRITE (stdout,'(5x,a/)')  '==================================================================='
  !
  IF (save_m_mat .EQ. .TRUE. .OR. save_m_matw .EQ. .TRUE. .OR. save_m_ph .EQ. .TRUE.) THEN
     !
     WRITE (stdout,'(5x,a/)')  'Save following variables to disk to reduce the using of internal memory :'
     !
     IF (save_m_mat) THEN
        WRITE (stdout,'(13x,a,f8.1,a)')  'epf17   ~', (DBLE(nq_ful)*1.0d-6)*DBLE(nbnd_red*nbnd_red*nmodes*2*DP), ' MB' 
     ENDIF
     !
     WRITE (stdout,'(13x,a,f8.1,a)')  'eimpf_full   ~', (DBLE(nq_ful)*1.0d-6)*DBLE(nbnd_red*nbnd_red*neptemp*nepdope*2*DP), ' MB'
     WRITE (stdout,'(13x,a,f8.1,a)')  'evc_fmesh   ~', (DBLE(nk_ful_red)*1.0d-6)*DBLE(nbnd_red*npol*ig_max_red*2*DP), ' MB'
     !
     IF (save_m_matw) THEN
        WRITE (stdout,'(13x,a,f8.1,a)')  'epmatwe ~', (DBLE(nrr_k*nqc)*1.0d-6)*DBLE(nbndsub*nbndsub*nmodes*2*DP), ' MB' 
        WRITE (stdout,'(13x,a,f8.1,a)')  'epmatwp ~', (DBLE(nrr_k*nrr_q)*1.0d-6)*DBLE(nbndsub*nbndsub*nmodes*2*DP), ' MB' 
     ENDIF
     !
     IF (save_m_ph) THEN
        WRITE (stdout,'(13x,a,f8.1,a)')  'uf_ful  ~', (DBLE(nq_ful)*1.0d-6)*DBLE(nmodes*nmodes*2*DP), ' MB' 
     ENDIF
     !
  ENDIF
  !
  IF (save_t_el .EQ. .TRUE.) THEN
     !
     WRITE (stdout,'(/5x,a/)') 'Save following variables to disk to reduce the time consumption :'
     !
     WRITE (stdout,'(13x,a,f8.1,a)')  'etf     ~', (DBLE(nk_ful)*1.0d-6)*DBLE(nbndsub*DP), ' MB' 
     WRITE (stdout,'(13x,a,f8.1,a)')  'etf_ks  ~', (DBLE(nk_ful)*1.0d-6)*DBLE(nbndsub*DP), ' MB' 
     WRITE (stdout,'(13x,a,f8.1,a)')  'cufkk   ~', (DBLE(nk_ful)*1.0d-6)*DBLE(nbndsub*nbndsub*2*DP), ' MB' 
     !
  ENDIF
  !
  !
  filename_check = TRIM(prefix)//'.epcheck'
  INQUIRE (FILE=filename_check,EXIST=file_exist)
  IF (file_exist) CALL ep_check (filename_check)
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
  ENDIF
  !
  ! globle variable
  ALLOCATE (epmatwef(nbndsub, nbndsub, nrr_q, nmodes))
  ALLOCATE (cufkk(nbndsub, nbndsub))
  ALLOCATE (cufkq(nbndsub, nbndsub))
  ALLOCATE (bmatf(nbndsub, nbndsub))
  ALLOCATE (epmatf(nbndsub, nbndsub, nmodes, neptemp, nepdope))
  ! if (screen_polar) then
  ALLOCATE (uf(nmodes, nmodes, neptemp, nepdope))
  ALLOCATE (etf(nbndsub, 2*nq_ful))
  ALLOCATE (etf_ks(nbndsub, 2*nq_ful)) 
  IF (.NOT. save_m_mat) ALLOCATE (epf17 (nq_ful, nbnd_red, nbnd_red, nmodes, neptemp, nepdope)) 
  !ALLOCATE (dmef(3, nbndsub, nbndsub, 2*nq_ful))
  IF (vme) ALLOCATE (vmef(3, nbndsub, nbndsub, 2*nq_ful))
  !
  ! JWZ, for impurity scattering
  if (eimp_mode == 2 .or. eimp_mode == 4 .or. eimp_mode == 6 .or. eimp_mode == 8) &
     ALLOCATE (vkq(3, nbndsub, 2*nq_ful))
  !
  if (eimp_mode > 0) then
     !
     if (eimp_mode < 7) then
        ALLOCATE ( eimpf_lr (nq_ful, neptemp, nepdope) )
        ALLOCATE ( eimpf_sr (nq_ful, nbnd_red, nbnd_red) )
        eimpf_lr = 0.d0
        eimpf_sr = 0.d0
     endif
     !
     if (eimp_mode == 5 .or. eimp_mode == 6) then
        ALLOCATE (eimpmatwef(nbndsub, nbndsub, nrr_q))
        ALLOCATE (eimpmatf(nbndsub, nbndsub))
        eimpmatwef = 0.d0
        eimpmatf = 0.d0

        IF (alloy_pot) THEN
           allocate ( alelpp (nq_ful, nbnd_red, nbnd_red) )
           alelpp = 0.d0
        ENDIF
     endif
     !
     if (eimp_mode == 7 .or. eimp_mode == 8) then

!        ALLOCATE (eemat_f (nbnd_red, nbnd_red, ng_max_red, nq_ful)) 
        ALLOCATE (eemat_f (nbnd_red, nbnd_red, npol, ng_max_red)) 

        ALLOCATE (dvimpqf_sr (ng_max_red))
        ALLOCATE (dvimpqf_lr (neptemp, nepdope))
        ALLOCATE (eimpmat_f (nbnd_red, nbnd_red, neptemp, nepdope))
        dvimpqf_sr = 0.d0
        dvimpqf_lr = 0.d0
        eimpmat_f = 0.d0

        ALLOCATE (eimpf_full (nbnd_red, nbnd_red, neptemp, nepdope, nq_ful))
        eimpf_full = 0.d0
     endif
     !
     if (eimp_mode == 7 .or. eimp_mode == 8) then
        !
        ! collect wavefunctions on fine mesh
        if (ig_max <= 1000) then
           ig_max_red = ig_max
        else
           ig_max_red = 1000
        endif
        ALLOCATE (evc_ks  (ig_max_red, npol, nbnd_red))
!        allocate (evc_fk  (ig_max, nbndsub))
!        allocate (evc_fkq (ig_max, nbndsub))
        ALLOCATE (evc_fk  (ig_max_red, npol, nbnd_red))
        ALLOCATE (evc_fkq (ig_max_red, npol, nbnd_red))
        ALLOCATE (map_mesh_rful (nkf1*nkf2*nkf3))   ! NOTE: only works when nkf is integer multiples of nqf
        ALLOCATE (mapcheck_rful (nkf1*nkf2*nkf3))   ! NOTE: only works when nkf is interger multiples of nqf
        ALLOCATE (map_mesh_ful (nkf1*nkf2*nkf3))   ! NOTE: only works when nkf is integer multiples of  nqf
        ALLOCATE (mapcheck_ful (nkf1*nkf2*nkf3))   ! NOTE: only works when nkf is integer multiples nqf
        map_mesh_rful = 0
        mapcheck_rful = 0
        map_mesh_ful = 0
        mapcheck_ful = 0
        !
        if (.not. epcheck) then
           ALLOCATE (evc_fmesh (ig_max_red, npol, nbnd_red, nk_ful_red))
        else
           ALLOCATE (evc_fmesh (ig_max_red, npol, nbnd_red, nk_ful))
        endif
        evc_fmesh = (0.d0, 0.d0)

        ! read ful_red mesh and mapping
        if (allocated(xkf_ful))   deallocate(xkf_ful)
        if (allocated(rful2ful))  deallocate(rful2ful)
        if (allocated(ful2rful))  deallocate(ful2rful)
        if (allocated(rful2rirr)) deallocate(rful2rirr)
        ALLOCATE(xkf_ful(3,nk_ful))
        ALLOCATE(rful2ful(nk_ful_red))
        ALLOCATE(ful2rful(nk_ful))
        ALLOCATE(rful2rirr(nk_ful_red))
        OPEN (1200,FILE='BTE/META/xkf_ful_cryst',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='old')
        do ik = 1, nk_ful
           READ (1200,REC=ik) xkf_ful(1:3,ik)
        enddo
        close (1200)
        OPEN (1201,FILE='BTE/META/rful2ful',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='old')
        DO ik = 1, nk_ful_red
           READ (1201,REC=ik) rful2ful(ik)
        ENDDO
        CLOSE (1201)
        OPEN (1202,FILE='BTE/META/ful2rful',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='old')
        DO ik = 1, nk_ful
           READ (1202,REC=ik) ful2rful(ik)
        ENDDO
        CLOSE (1202)
        !
        OPEN (1203,FILE='BTE/META/rful2rirr',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='old')
        DO ik = 1, nk_ful_red
           READ (1203,REC=ik) rful2rirr(ik)
        ENDDO
        CLOSE (1203)
        !
        !
        if (.not. epcheck) then
           !
           CALL para_bounds (ik_star, ik_stop, nk_ful_red)
           DO ik_red = ik_star, ik_stop
              !
              ik = rful2ful(ik_red)       ! ik is the original index of k point in full BZ
              xkk(:) = xkf_ful(:,ik)     ! xkf_ful is assumed to be in crys coord
   
              CALL hamwan2bloch (nbndsub, nrr_k, irvec, ndegen_k, xkk, cufkk, etf_ks_tmp, chw)
              CALL evcwan2bloch (xkk, irvec, ndegen_k, cufkk, nbndsub, nrr_k, evc_ks, ig_max_red)
              evc_fmesh(:,:,:,ik_red) = evc_ks(:,:,:)
              !
              do while (xkk(1) < 0)
                 xkk(1) = xkk(1) + 1
              enddo
              do while (xkk(1) >= 1)
                 xkk(1) = xkk(1) - 1
              enddo
              do while (xkk(2) < 0)
                 xkk(2) = xkk(2) + 1
              enddo
              do while (xkk(2) >= 1)
                 xkk(2) = xkk(2) - 1
              enddo
              do while (xkk(3) < 0)
                 xkk(3) = xkk(3) + 1
              enddo
              do while (xkk(3) >= 1)
                 xkk(3) = xkk(3) - 1
              enddo
              ijk_k(1) = nint(xkk(1)*nkf1)
              ijk_k(2) = nint(xkk(2)*nkf2)
              ijk_k(3) = nint(xkk(3)*nkf3)
              map_mesh_rful(ijk_k(1)*nkf2*nkf3+ijk_k(2)*nkf3+ijk_k(3)+1) = ik_red
              mapcheck_rful(ijk_k(1)*nkf2*nkf3+ijk_k(2)*nkf3+ijk_k(3)+1) = 1
              !
           ENDDO
#ifdef __PARA
           CALL mp_barrier (inter_pool_comm)
           CALL mp_sum (evc_fmesh,inter_pool_comm)
           CALL mp_sum (map_mesh_rful,inter_pool_comm)
           CALL mp_sum (mapcheck_rful,inter_pool_comm)
#ENDIF
           do ik = 1, nkf1*nkf2*nkf3
              if (mapcheck_rful(ik) > 1) then
                 write(stdout,*) ' ERROR: mapping to overlapped k point. Stop'
                 stop
              endif
           enddo
           !
        endif

        CALL para_bounds (ik_star, ik_stop, nk_ful)
        DO ik_red = ik_star, ik_stop
           !
           ik = ik_red
           xkk(:) = xkf_ful(:,ik)     ! xkf_irr is assumed to be in crys coord

           if (epcheck) then
              CALL hamwan2bloch (nbndsub, nrr_k, irvec, ndegen_k, xkk, cufkk, etf_ks_tmp, chw)
              CALL evcwan2bloch (xkk, irvec, ndegen_k, cufkk, nbndsub, nrr_k, evc_ks, ig_max_red)
              evc_fmesh(:,:,:,ik_red) = evc_ks(:,:,:)
           endif
           !
           do while (xkk(1) < 0)
              xkk(1) = xkk(1) + 1
           enddo
           do while (xkk(1) >= 1)
              xkk(1) = xkk(1) - 1
           enddo
           do while (xkk(2) < 0)
              xkk(2) = xkk(2) + 1
           enddo
           do while (xkk(2) >= 1)
              xkk(2) = xkk(2) - 1
           enddo
           do while (xkk(3) < 0)
              xkk(3) = xkk(3) + 1
           enddo
           do while (xkk(3) >= 1)
              xkk(3) = xkk(3) - 1
           enddo
           ijk_k(1) = nint(xkk(1)*nkf1)
           ijk_k(2) = nint(xkk(2)*nkf2)
           ijk_k(3) = nint(xkk(3)*nkf3)
           map_mesh_ful(ijk_k(1)*nkf2*nkf3+ijk_k(2)*nkf3+ijk_k(3)+1) = ik_red
           mapcheck_ful(ijk_k(1)*nkf2*nkf3+ijk_k(2)*nkf3+ijk_k(3)+1) = 1
           !
        ENDDO
#ifdef __PARA
        CALL mp_barrier (inter_pool_comm)
        if (epcheck) CALL mp_sum (evc_fmesh,inter_pool_comm)
        CALL mp_sum (map_mesh_ful,inter_pool_comm)
        CALL mp_sum (mapcheck_ful,inter_pool_comm)
#ENDIF
        do ik = 1, nkf1*nkf2*nkf3
           if (mapcheck_ful(ik) > 1) then
              write(stdout,*) ' ERROR: mapping to overlapped k point. Stop'
              stop
           endif
        enddo
        !
     endif
  endif
  !


  !
  ! ======================================================
  ! dynamical matrix : Wannier -> Bloch
  ! ======================================================
  !
  ! Prepare phonon eigenvector of each iq first for saving computation
  ! parralelized version
  if (screen_polar) then
     ALLOCATE (wf_ful(nmodes,neptemp,nepdope,nq_ful))
  else
     ALLOCATE (wf_ful(nmodes,1,1,nq_ful))
  endif
  ALLOCATE (vph_ful(3,nmodes,nq_ful))
  wf_ful  = 0.0d0
  vph_ful = 0.0d0
  IF (.NOT. save_m_ph) THEN
     if (screen_polar) then
        ALLOCATE (uf_ful(nmodes,nmodes,neptemp,nepdope,nq_ful))
     else
        ALLOCATE (uf_ful(nmodes,nmodes,1,1,nq_ful))
     endif
     uf_ful = (0.0d0,0.0d0)
  ENDIF
  ! 
  WRITE (stdout,'(/5x,a,f8.2,a)') 'Prepare phonon properties :'
  !
  IF (ph_read) THEN
     !
     CALL phonon_read
     !
  ELSE
     !
     CALL CPU_TIME (t1i)
     !
     CALL para_bounds (iq_star, iq_stop, nq_ful)
     !
     dyn_num = 0
     dyn_num(my_pool_id+1) = iq_stop-iq_star+1
     !
#ifdef __PARA
     CALL mp_barrier (inter_pool_comm)
     CALL mp_sum (dyn_num,inter_pool_comm)
#ENDIF
     !
     !
#ifdef __PARA
     CALL set_ndnmbr (0,my_pool_id+1,1,npool,dyn_cpu)
     dyn_ufmt = TRIM(tmp_dir)//TRIM(prefix)//'.dyn_'//dyn_cpu
#ENDIF
     if (screen_polar) then
        OPEN (41214,FILE=dyn_ufmt,FORM='unformatted',ACCESS='direct',RECL=(4*nmodes+2*nmodes*nmodes)*neptemp*nepdope*DP,STATUS='replace') 
     else
        OPEN (41214,FILE=dyn_ufmt,FORM='unformatted',ACCESS='direct',RECL=(4*nmodes+2*nmodes*nmodes)*DP,STATUS='replace') 
     endif
     !
     if (screen_polar) then
        ALLOCATE (wf_dyn(nmodes, neptemp, nepdope))
        ALLOCATE (vph_dyn(3, nmodes, neptemp, nepdope))
        ALLOCATE (uf_dyn(nmodes, nmodes, neptemp, nepdope))
     else
        ALLOCATE (wf_dyn(nmodes, 1, 1))
        ALLOCATE (vph_dyn(3, nmodes, 1, 1))
        ALLOCATE (uf_dyn(nmodes, nmodes, 1, 1))
     endif
     !
     wf_dyn  = 0.0d0
     vph_dyn = 0.0d0
     uf_dyn  = (0.0d0,0.0d0)
     !
     DO iq = iq_star, iq_stop
        !
        xxq = xqf_ful(:,iq)
        !
        if (screen_polar) then
           CALL dynwan2bloch_s (nmodes, nrr_q, irvec, ndegen_q, xxq, uf_dyn, wf_dyn, vph_dyn, neptemp, nepdope)
        else
           CALL dynwan2bloch (nmodes, nrr_q, irvec, ndegen_q, xxq, uf_dyn, wf_dyn, vph_dyn)
        endif
        !
        IF (epdim .EQ. 2) vph_dyn(3,:,:,:) = 0.0d0
        !
        DO nu = 1, nmodes
           !
           if (screen_polar) then
              do itemp = 1, neptemp
                 do idope = 1, nepdope
                    IF (wf_dyn(nu,itemp,idope) .GT. 0.0d0) THEN
                       wf_dyn(nu,itemp,idope) =  SQRT(ABS(wf_dyn(nu,itemp,idope)))
                    ELSE
                       wf_dyn(nu,itemp,idope) = -SQRT(ABS(wf_dyn(nu,itemp,idope)))
                    ENDIF
                 enddo
              enddo
           else
              IF (wf_dyn(nu,1,1) .GT. 0.0d0) THEN
                 wf_dyn(nu,1,1) =  SQRT(ABS(wf_dyn(nu,1,1)))
              ELSE
                 wf_dyn(nu,1,1) = -SQRT(ABS(wf_dyn(nu,1,1)))
              ENDIF
           endif
           !
           DO mu = 1, nmodes
              na = (mu-1)/3 + 1
              if (screen_polar) then
                 uf_dyn(mu,nu,:,:) = uf_dyn(mu,nu,:,:) / SQRT(amass(ityp(na)))
              else
                 uf_dyn(mu,nu,1,1) = uf_dyn(mu,nu,1,1) / SQRT(amass(ityp(na)))
              endif
           ENDDO
           !
        ENDDO 
        !
        IF (vg_ph .EQ. 'linear') THEN
           !
           if (screen_polar) then
              CALL vgwan2bloch_q (0.0d0, xxq, vph_dyn, neptemp, nepdope)
           else
              CALL vgwan2bloch_q (0.0d0, xxq, vph_dyn, 1, 1)
           endif
           IF (epdim .EQ. 2) vph_dyn(3,:,:,:) = 0.0d0
           !
        ENDIF
        !
        if (screen_polar) then
           WRITE (41214,REC=iq-iq_star+1) wf_dyn(1:nmodes,1:neptemp,1:nepdope), vph_dyn(1:3,1:nmodes,1:neptemp,1:nepdope), &
                                          uf_dyn(1:nmodes,1:nmodes,1:neptemp,1:nepdope)
        else
           WRITE (41214,REC=iq-iq_star+1) wf_dyn(1:nmodes,1,1), vph_dyn(1:3,1:nmodes,1,1), &
                                          uf_dyn(1:nmodes,1:nmodes,1,1)
        endif
        !
     ENDDO
     !
     CLOSE (41214)
     !
     !
     CALL mp_barrier (inter_pool_comm)
     !
     !
     IF (my_pool_id .EQ. ionode_id) THEN
        !
        wf_ufmt  = TRIM(tmp_dir)//TRIM(prefix)//'.phw'
        uf_ufmt  = TRIM(tmp_dir)//TRIM(prefix)//'.phu'
        vph_ufmt = TRIM(tmp_dir)//TRIM(prefix)//'.phv'  

        if (screen_polar) then
           OPEN (99901,FILE=wf_ufmt,FORM='unformatted',ACCESS='direct',RECL=nmodes*neptemp*nepdope*DP,STATUS='replace')
           OPEN (99902,FILE=uf_ufmt,FORM='unformatted',ACCESS='direct',RECL=2*nmodes*nmodes*neptemp*nepdope*DP,STATUS='replace')
           OPEN (99903,FILE=vph_ufmt,FORM='unformatted',ACCESS='direct',RECL=3*nmodes*neptemp*nepdope*DP,STATUS='replace')
        else
           OPEN (99901,FILE=wf_ufmt,FORM='unformatted',ACCESS='direct',RECL=nmodes*DP,STATUS='replace')
           OPEN (99902,FILE=uf_ufmt,FORM='unformatted',ACCESS='direct',RECL=2*nmodes*nmodes*DP,STATUS='replace')
           OPEN (99903,FILE=vph_ufmt,FORM='unformatted',ACCESS='direct',RECL=3*nmodes*DP,STATUS='replace')
        endif
        !
        iq_dyn = 1
        !
        DO icpu = 1, npool
           !
           WRITE(dyn_cpu,'(i3)') icpu
           dyn_ufmt = TRIM(tmp_dir)//TRIM(prefix)//'.dyn_'//TRIM(ADJUSTL(dyn_cpu))
           if (screen_polar) then
              OPEN (41214,FILE=dyn_ufmt,FORM='unformatted',ACCESS='direct',RECL=(4*nmodes+2*nmodes*nmodes)*neptemp*nepdope*DP,STATUS='old') 
           else
              OPEN (41214,FILE=dyn_ufmt,FORM='unformatted',ACCESS='direct',RECL=(4*nmodes+2*nmodes*nmodes)*DP,STATUS='old') 
           endif
           !
           DO iq = 1, dyn_num(icpu)
              !
              if (screen_polar) then
                 READ (41214,REC=iq) wf_dyn(1:nmodes,1:neptemp,1:nepdope), vph_dyn(1:3,1:nmodes,1:neptemp,1:nepdope), &
                                     uf_dyn(1:nmodes,1:nmodes,1:neptemp,1:nepdope)
                 WRITE (99901,REC=iq_dyn) wf_dyn(1:nmodes,1:neptemp,1:nepdope)
                 WRITE (99902,REC=iq_dyn) uf_dyn(1:nmodes,1:nmodes,1:neptemp,1:nepdope)
                 WRITE (99903,REC=iq_dyn) vph_dyn(1:3,1:nmodes,1:neptemp,1:nepdope)
              else
                 READ (41214,REC=iq) wf_dyn(1:nmodes,1,1), vph_dyn(1:3,1:nmodes,1,1), uf_dyn(1:nmodes,1:nmodes,1,1)
                 WRITE (99901,REC=iq_dyn) wf_dyn(1:nmodes,1,1)
                 WRITE (99902,REC=iq_dyn) uf_dyn(1:nmodes,1:nmodes,1,1)
                 WRITE (99903,REC=iq_dyn) vph_dyn(1:3,1:nmodes,1,1)
              endif
              iq_dyn = iq_dyn + 1
              !
           ENDDO
           !
           CLOSE (41214)
           !
        ENDDO
        !
        CALL SYSTEM ('rm '//TRIM(tmp_dir)//TRIM(prefix)//'*.dyn_*')
        CLOSE (99901)
        CLOSE (99902)
        CLOSE (99903)
        !
     ENDIF
     !
     !
     CALL mp_barrier (inter_pool_comm)
     !
     CALL phonon_read
     !
     !
     CALL CPU_TIME (t1f)
     t1 = t1f-t1i
     WRITE (stdout,'(/13x,a,f8.2,a)') 'Out k | Out q | dyn     W->B :', t1, 's'
     !
  ENDIF
  !
  ! find maximun frequency
  !
  ! for screen_polar case, we only use one temperature, 
  ! carrier concentration combination to find maximum frequency
  !
  wmax_mode(1:nmodes) = MAXVAL(wf_ful(1:nmodes,1,1,:))
  wmax = MAXVAL(wmax_mode(:))
  !
  !
  ! ======================================================
  ! eph-vertex : Wannier -> Bloch
  ! ======================================================
  !
  IF (lpolar .and. elop) THEN
     !
     ALLOCATE(eph_vogl(nmodes,neptemp,nepdope,nq_ful))
     eph_vogl = 0.0d0
     !
     WRITE (stdout,'(/5x,a,f8.2,a)') 'Prepare long-range part of e-ph matrix :'
     !
     IF (ephl_read) THEN
        !
        CALL eph_long_read
        !
     ELSE
        !
        CALL CPU_TIME (tbi)
        !

        CALL para_bounds (iq_star, iq_stop, nq_ful)
        !
        vogl_num = 0
        vogl_num(my_pool_id+1) = iq_stop-iq_star+1
        !
#ifdef __PARA
        CALL mp_barrier (inter_pool_comm)
        CALL mp_sum (vogl_num,inter_pool_comm)
#ENDIF
        !
        !
#ifdef __PARA
        CALL set_ndnmbr (0,my_pool_id+1,1,npool,vogl_cpu)
        vogl_ufmt = TRIM(tmp_dir)//TRIM(prefix)//'.ephl_'//vogl_cpu
#ENDIF
        OPEN (41914,FILE=vogl_ufmt,FORM='unformatted',ACCESS='direct',RECL=2*nmodes*neptemp*nepdope*DP,STATUS='replace') 
        !
        ephl_vogl  = (0.0d0,0.0d0)
        !
        !
        uf_cart  = 0.0d0
        DO imode = 1, nmodes
           uf_cart(imode,imode) = 1.0d0
        ENDDO
        !
        DO iq = iq_star, iq_stop
           !
           xxq = xqf_ful (:,iq)
           CALL cryst_to_cart (1, xxq, bg, 1) 
           !
           IF (MAXVAL(ABS(xxq)) .NE. 0.0d0) THEN
              CALL polar_eph (uf_cart, xxq, ephl_vogl, neptemp, nepdope)
           ELSE
              ephl_vogl = 0.0d0
           ENDIF
           !
           WRITE (41914,REC=iq-iq_star+1) ephl_vogl(1:nmodes,1:neptemp,1:nepdope)
           !
        ENDDO
        !
        CLOSE (41914)
        !
        !
        CALL mp_barrier (inter_pool_comm)
        !
        !
        IF (my_pool_id .EQ. ionode_id) THEN
           !
           vogl_ufmt  = TRIM(tmp_dir)//TRIM(prefix)//'.ephl'
           OPEN (99901,FILE=vogl_ufmt,FORM='unformatted',ACCESS='direct',RECL=2*nmodes*neptemp*nepdope*DP,STATUS='replace')
           !
           iq_vogl = 1
           !
           DO icpu = 1, npool
              !
              WRITE(vogl_cpu,'(i3)') icpu
              vogl_ufmt = TRIM(tmp_dir)//TRIM(prefix)//'.ephl_'//TRIM(ADJUSTL(vogl_cpu))
              OPEN (41914,FILE=vogl_ufmt,FORM='unformatted',ACCESS='direct',RECL=2*nmodes*neptemp*nepdope*DP,STATUS='old') 
              !
              DO iq = 1, vogl_num(icpu)
                 !
                 READ (41914,REC=iq) ephl_vogl(1:nmodes,1:neptemp,1:nepdope)
                 WRITE (99901,REC=iq_vogl) ephl_vogl(1:nmodes,1:neptemp,1:nepdope)
                 iq_vogl = iq_vogl + 1
                 !
              ENDDO
              !
              CLOSE (41914)
              !
           ENDDO
           !
           CALL SYSTEM ('rm '//TRIM(tmp_dir)//TRIM(prefix)//'.ephl_*')
           CLOSE (99901)
           !
        ENDIF
        !
        !
        CALL mp_barrier (inter_pool_comm)
        !
        CALL eph_long_read
        !
        CALL CPU_TIME (tbf)
        tb = tbf - tbi
        WRITE (stdout,'(/13x,a,f8.2,a)') 'Out k | Out q | eph-L   W->B :', tb, 's'
        !
     ENDIF
     !
  ENDIF
  !
  !
  ! ======================================================
  ! ham : Wannier -> Bloch
  ! ======================================================
  !
  IF (save_t_el) THEN
     !
     WRITE (stdout,'(/5x,a,f8.2,a)') 'Prepare electron properties :'
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
        !DO ibnd = 1, nbndsub
        !   IF (etf_ham(ibnd) .GT. ef_m) etf_ham(ibnd) = etf_ham(ibnd) + delta_egap
        !ENDDO
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
     CALL mp_barrier (inter_pool_comm)
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
     mkq1 = nkf1/nqf1
     mkq2 = nkf2/nqf2
     mkq3 = nkf3/nqf3
     !
     CALL CPU_TIME (tcf)
     tc = tcf-tci
     WRITE (stdout,'(/13x,a,f8.2,a)') 'Out k | Out q | ham     W->B :', tc, 's'
     !
     !
     CALL mp_barrier (inter_pool_comm)
     !
  ENDIF
  !
  ! ======================================================
  !                   impurity potential
  ! ======================================================
  !
  !
  if (eimp_mode > 2 .and. eimp_mode < 5) then
     !
     CALL CPU_TIME (tci)
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
                                   !    should be consistent with the value used here [alat]

        read(1000,*) def_pos(1:3)
        read(1000,*) scella(1:3)   ! in [alat_imp]
        read(1000,*) scellb(1:3)
        read(1000,*) scellc(1:3)
        read(1000,*) nfftmesh

        nfft_sc = nfftx_sc * nffty_sc * nfftz_sc
        !
        ALLOCATE(rp(3,nfftmesh), dvimp(nfftmesh))
        do i = 1, nfftmesh
           read(1000,*) rp(1:3,i), dvimp(i)   ! rp in [alat_imp], dvimp in [Ryd]
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
        if (.not. lpolar) then
           dielec0 = dielec
        else
           dielec0 = (epsi(1,1) + epsi(2,2) + epsi(3,3)) / 3.d0
        endif
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
        !
        ! to avoid arbitrary large number due to Coulomb correction
        dvimp(center_index) = 0.d0
        !
        write(stdout,*) ' [min, max] of corrected impurity potential, eV:', minval(dvimp)*ryd2ev, maxval(dvimp)*ryd2ev
        !
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
     ! normalization over space, for coupling matrix integration
     !
     dvimp = dvimp * dvol / omega
     ! 
     ! now calculate short-range coupling matrix for impurity scattering
     !
     CALL para_bounds (iq_star, iq_stop, nq_ful)
     !
     do iq = iq_star, iq_stop
        !
        xxq(:) = xqf_ful(:,iq)
        !
        ! !IMPORTANT NOTE! this brings xq to first-BZ, but it is only centered around
        ! Gamma, not yet fully symmetric
        xq_cart(1) = xxq(1) - NINT(xxq(1))
        xq_cart(2) = xxq(2) - NINT(xxq(2))
        xq_cart(3) = xxq(3) - NINT(xxq(3))
        !
        CALL cryst_to_cart ( 1, xq_cart, bg, 1 )   ! xq_cart in [2*pi/alat]
        !
        ! NOTE: xq_cart in [2*pi/alat], while rp in [alat_imp]
        !       according to above, [alat] = [alat_imp]
        !
        if (.not. alloy_pot) then
           eimpf_sr_ = 0.d0
           do imesh = 1, nfftmesh
              ! * dvol / omega is added in after reading dvimp
              ! better to check the plane wave approximation here
              eimpf_sr_ = eimpf_sr_ + exp(-ci*twopi*dot_product(xq_cart,rp(:,imesh))) * &
                                           dvimp(imesh)
           enddo
           eimpf_sr(iq,:,:) = eimpf_sr_
        else
           alelpp_ = 0.d0
           do imesh = 1, nfftmesh
              alelpp_ = alelpp_ + exp(-ci*twopi*dot_product(xq_cart,rp(:,imesh))) * &
                                           dvimp(imesh)
           enddo
           alelpp(iq,:,:) = alelpp_
        endif
     enddo
     !
     eimpf_sr = real(eimpf_sr)
     IF (alloy_pot) alelpp = real(alelpp)
     !
#ifdef __PARA
     CALL mp_sum (eimpf_sr,inter_pool_comm)
     IF (alloy_pot)  CALL mp_sum (alelpp,inter_pool_comm)
#ENDIF
     !
     CALL CPU_TIME (tcf)
     tc = tcf - tci
     WRITE (stdout,'(/5x,a,f8.2,a)') 'Prepare impurity potential:', tc, 's'
     !
  endif
  !
  !
  !
  if (eimp_mode == 7 .or. eimp_mode == 8) then
     !
     if (.not. lpolar) then
        dielec0 = dielec
     else
        dielec0 = (epsi(1,1) + epsi(2,2) + epsi(3,3)) / 3.d0
     endif
     !
  endif
  !
  !
  ! -------------------- ep/eimp check --------------------
  !
  filename = trim(adjustl(prefix)) // '.eph_fine_mesh_cond.txt'
  open(unit = 30401, file = filename)
  filename = trim(adjustl(prefix)) // '.eph_fine_mesh_vale.txt'
  open(unit = 30403, file = filename)
  write(30401,*) nqf, nbnd_c
  write(30403,*) nqf, nbnd_v
  !
  if (eimp_mode == 5 .or. eimp_mode == 6) then
     !
     filename = trim(adjustl(prefix)) // '.eimp_fine_mesh_cond.txt'
     open(unit = 30405, file = filename)
     filename = trim(adjustl(prefix)) // '.eimp_fine_mesh_vale.txt'
     open(unit = 30406, file = filename)
     write(30405,*) nqf, nbnd_c
     write(30406,*) nqf, nbnd_v
     !   
  endif

  ! ==============================
  ! Check electron band structure
  ! ==============================
  if ((epcheck) .and. (my_pool_id .EQ. ionode_id)) then
     !
     open(unit=10101, file='cz.dat')
     ! here we output the eigenvalues on the fine mesh for checking
     allocate (eig_k(nbndsub))
     open (2105, file='band_k.in')
     open (2106, file='band_energy_fine.txt')
     read (2105,*) nk_band
     read (2105,*) nu
     write(2106,*) nk_band*nu+1, nbndsub
     !

     do mu = 1, nu
        read(2105,*) xkk, xkk2
        if (mu==1) then
           xxq = xkk
           CALL hamwan2bloch &
                ( nbndsub, nrr_k, irvec, ndegen_k, xxq, cufkk, eig_k, chw)
           !
           write (2106,'(f12.8,f12.8,f12.8,48f10.4)') &
                xxq(1:3), (eig_k(ibnd)*ryd2ev, ibnd=1,nbndsub)
        endif
        do iq = 1, nk_band
           xxq = xkk + ((xkk2-xkk)*iq)/nk_band
           CALL hamwan2bloch &
                ( nbndsub, nrr_k, irvec, ndegen_k, xxq, cufkk, eig_k, chw)
           !
           write (2106,'(f12.8,f12.8,f12.8,48f10.4)') &
                xxq(1:3), (eig_k(ibnd)*ryd2ev, ibnd=1,nbndsub)
        enddo
     enddo
     !
     deallocate (eig_k)
     close (2105)
     close (2106)
     close(10101)
  endif

  ! ==============================
  ! Check phonon dispersion
  ! ==============================
  !
  IF ((epcheck) .and. (my_pool_id .EQ. ionode_id)) THEN
     !
     ! output phonon frequency for checking

     open(2105, file='band_q.in')
     filename = trim(adjustl(prefix)) // '.phonon_freq_fine.txt'
     open(2106, file=filename)
     read(2105,*) nk_band
     read(2105,*) nu
     if (.not. screen_polar) then
        write(2106,*) nk_band*nu+1, nmodes
     else
        write(2106,*) nk_band*nu+1, nmodes, neptemp, nepdope
     endif
     !
     if (.not. screen_polar) then
        allocate (eig_q(nmodes, 1, 1, nk_band*nu+1))
        allocate (uf_(nmodes, nmodes, 1, 1))
     else
        allocate (eig_q(nmodes, neptemp, nepdope, nk_band*nu+1))
        allocate (uf_(nmodes, nmodes, neptemp, nepdope))
     endif
     allocate (xxq_band(3, nk_band*nu+1))
     uf_ = (0.0d0,0.0d0)
     !
     ik_band = 0
     do mu = 1, nu
!        write(stdout,*) ' direction:',mu
        read(2105,*) xkk, xkk2
        if (mu==1) then
           ik_band = ik_band + 1
           xxq_band(:,ik_band) = xkk
           if (eph_interp) then
!              call ph_freq_interp (nmodes, xxq, eig_k)
              if (screen_polar) then
                 write(stdout,*) ' eph_interp cannot be used with screen_polar'
                 stop
              endif
!CHECK: add this subroutine
!              call ph_freq_interp_tetra (nmodes, xxq_band(:,ik_band), eig_q(:,1,1,1))
              write (2106,'(f12.8,f12.8,f12.8,48f10.4)') &
                   xxq(1:3), (abs(eig_q(ibnd,1,1,1))*rydcm1, ibnd=1,nmodes)
           else
              if (screen_polar) then
                 CALL dynwan2bloch_s &
                      ( nmodes, nrr_q, irvec, ndegen_q, xxq_band(:,ik_band), uf_, eig_q(:,:,:,ik_band), vq_,neptemp,nepdope)
              else
                 CALL dynwan2bloch &
                      ( nmodes, nrr_q, irvec, ndegen_q, xxq_band(:,ik_band), uf_, eig_q(:,1,1,ik_band), vq_)
              endif
           endif
        endif
        !
        do iq = 1, nk_band
           ik_band = ik_band + 1
           xxq_band(:,ik_band) = xkk + ((xkk2-xkk)*iq)/nk_band
!           write(stdout,*) ' iq:',iq, '-------------------------------'
!           write(stdout,*) ' xxq =',xxq(1:3)
           if (eph_interp) then
!              call ph_freq_interp (nmodes, xxq, eig_k)
!CHECK: add this subroutine
!              call ph_freq_interp_tetra (nmodes, xxq_band(:,ik_band), eig_q(:,1,1,1))
              write (2106,'(f12.8,f12.8,f12.8,48f10.4)') &
                   xxq(1:3), (abs(eig_q(ibnd,1,1,1))*rydcm1, ibnd=1,nmodes)
           else
              if (screen_polar) then
                 CALL dynwan2bloch_s &
                      ( nmodes, nrr_q, irvec, ndegen_q, xxq_band(:,ik_band), uf_, eig_q(:,:,:,ik_band), vq_,neptemp,nepdope)
              else
                 CALL dynwan2bloch &
                      ( nmodes, nrr_q, irvec, ndegen_q, xxq_band(:,ik_band), uf_, eig_q(:,1,1,ik_band), vq_)
              endif
           endif
        enddo
     enddo
     !
     do itemp = 1, neptemp
        do idope = 1, nepdope
           if (.not. screen_polar) then
              if ((itemp > 1) .or. (idope > 1)) cycle
           endif
           !
           ik_band = 0
           do mu = 1, nu
              if (mu==1) then
                 ik_band = ik_band + 1
                 write (2106,'(f12.8,f12.8,f12.8,48f10.4)') &
                      xxq_band(1:3,ik_band), (sqrt(abs(eig_q(ibnd,itemp,idope,ik_band)))*rydcm1, ibnd=1,nmodes)
              endif
              do iq = 1, nk_band
                 ik_band = ik_band + 1
                 write (2106,'(f12.8,f12.8,f12.8,48f10.4)') &
                      xxq_band(1:3,ik_band), (sqrt(abs(eig_q(ibnd,itemp,idope,ik_band)))*rydcm1, ibnd=1,nmodes)
              enddo
           enddo
        enddo
     enddo
     !
     deallocate (eig_q)
     deallocate (uf_)
     close(2105)
     close(2106)
  ENDIF

  !
  ! ======================================================
  !                         k loop
  ! ======================================================
  !
  ! parallelization for nk_irr_red
  CALL para_bounds (ik_star, ik_stop, nk_irr_red)
  !
  ! BTE and selfen variable
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
  IF (bte .EQ. 0 .OR. bte .EQ. 1 .OR. bte .EQ. 3) THEN
     ALLOCATE (nscat_all(nk_irr_red))
     nscat_all = 0 
  ENDIF
  !
  !
  DO ik_red = ik_star, ik_stop
     !
     CALL CPU_TIME (t0i)
     !
     CALL start_clock ( 'ep-interp' )
     !
     IF (bte .EQ. 19 .OR. bte .EQ. 18) THEN
        ik = ik_red
        xkk(:) = xk_kmesh(:,ik)
     ELSE
        ik = rirr2irr(ik_red)       ! ik is the original index of k point in irreducible BZ
        xkk(:) = xkf_irr(:,ik)     ! xkf_irr is assumed to be in crys coord
     ENDIF
!     xkk = -xkk
     !
     write(*,*) 'my pool id #', my_pool_id, 'xk =', xkk
     !
     ! find special k-points
     gamma_find = .false.
     if (epcheck) then
        if ((abs(xkk(1))==0).and.(abs(xkk(2))==0).and.(abs(xkk(3))==0)) then
           gamma_find = .true.
           write(*,*) ' finding gamma-point at #pool:',my_pool_id
           write(*,*) ' xk =',xkk(1:3)
        endif
     endif
     !
     x_find = .false.
     if (epcheck) then
        if (((abs(xkk(1))==0.0).and.(abs(xkk(2))==0.5).and.(abs(xkk(3))==0.5)) .or. &
            ((abs(xkk(1))==0.5).and.(abs(xkk(2))==0.0).and.(abs(xkk(3))==0.5)) .or. &
            ((abs(xkk(1))==0.5).and.(abs(xkk(2))==0.5).and.(abs(xkk(3))==0.0))) then
           x_find = .true.
           write(*,*) ' finding x-point at #pool:',my_pool_id
           write(*,*) ' xk =',xkk(1:3)
        endif
     endif
     !
     l_find = .false.
     if (epcheck) then
        if (((abs(xkk(1))==0.0).and.(abs(xkk(2))==0.0).and.(abs(xkk(3))==0.5)) .or. &
            ((abs(xkk(1))==0.0).and.(abs(xkk(2))==0.5).and.(abs(xkk(3))==0.0)) .or. &
            ((abs(xkk(1))==0.5).and.(abs(xkk(2))==0.0).and.(abs(xkk(3))==0.0))) then
           l_find = .true.
           write(*,*) ' finding l-point at #pool:',my_pool_id
           write(*,*) ' xk =',xkk(1:3)
        endif
     endif

     cond_find = .false.
     if ((gamma_find .and. (xk_c=='G')) .or. (x_find .and. (xk_c=='X'))) then
        cond_find = .true.
        write(*,*) ' cond_find at #pool:',my_pool_id
!        if ((epcheck).and.(filqf/='')) then
!           filename = trim(adjustl(prefix)) // '.eph_interpolated_fine_mesh_cond.txt'
!           open(unit = 30402, file = filename)
!           write(30402,'(2(2x,i5))') nqf, nbnd_c
!           cond_output = .true.
!        endif
     endif

     vale_find = .false.
     if ((gamma_find .and. (xk_v=='G')) .or. (x_find .and. (xk_v=='X')) .or. &
         (l_find .and. (xk_v=='L'))) then
        vale_find = .true.
        write(*,*) ' vale_find at #pool:',my_pool_id
!        if ((epcheck).and.(filqf/='')) then
!           filename = trim(adjustl(prefix)) // '.eph_interpolated_fine_mesh_vale.txt'
!           open(unit = 30404, file = filename)
!           write(30404,'(2(2x,i5))') nqf, nbnd_v
!           vale_output = .true.
!        endif
     endif

     if ((epcheck) .and. (.not. (cond_find .or. vale_find))) then
        cycle
     endif

     !
     !
     ! ======================================================
     ! hamiltonian : Wannier -> Bloch
     ! ======================================================
     CALL CPU_TIME (t2i)
     !
     IF (eig_read) THEN
        CALL hamwan2bloch (nbndsub, nrr_k, irvec, ndegen_k, xkk, cufkk, etf_ks(:,1), chw_ks)     
     ENDIF
     !
     CALL hamwan2bloch (nbndsub, nrr_k, irvec, ndegen_k, xkk, cufkk, etf(:,1), chw)        
     ! 
     !DO ibnd = 1, nbndsub
     !   IF (etf(ibnd,1) .GT. ef_m) etf(ibnd,1) = etf(ibnd,1) + delta_egap
     !ENDDO
        DO ibnd = vbnd_num+1, nbndsub
           etf(ibnd,1) = etf(ibnd,1) + delta_egap
        ENDDO  
     !
     CALL CPU_TIME (t2f)
     t2 = t2f - t2i
     !
     !
     ! ======================================================
     !  dipole: Wannier -> Bloch
     ! ======================================================
     CALL CPU_TIME (t3i)
     !
     IF (eimp_mode == 2 .or. eimp_mode == 4 .or. eimp_mode == 6 .or. eimp_mode == 8) THEN
        !
        IF (vg_el .EQ. 'matrix') THEN
           !
           CALL dmewan2bloch (nbndsub, nrr_k, irvec, ndegen_k, xkk, cufkk, dmef_tmp, etf(:,1), etf_ks(:,1))
           !
           DO ibnd = 1, nbndsub
              vkq(:,ibnd,1) = 2.0d0*REAL(dmef_tmp(:,ibnd,ibnd))
           ENDDO
           !
        ELSEIF (vg_el .EQ. 'linear') THEN
           xxq_ = 0.d0
!           CALL velwan2bloch (xkk, xxq_, etf(:,1), chw, nrr_k, irvec, ndegen_k, velf)

           CALL vgwan2bloch_k (xkk, 0.0d0, velf)
           IF (epdim .EQ. 2) velf(3,:) = 0.0d0
           !
           DO ibnd = 1, nbndsub
              vkq(:,ibnd,1) = velf(:,ibnd)
           ENDDO
           !
        ENDIF
     ENDIF
     !
     CALL CPU_TIME (t3f)
     t3 = t3f - t3i
     !
     ! wavefunction at k

     if (eimp_mode == 7 .or. eimp_mode == 8) then
        !
!        call evcwan2bloch (xkk, irvec, ndegen_k, cufkk, nbndsub, nrr_k, evc_fk, ig_max_red)

        xkk_ = xkk
        do while (xkk_(1) < 0)
           xkk_(1) = xkk_(1) + 1
        enddo
        do while (xkk_(1) >= 1)
           xkk_(1) = xkk_(1) - 1
        enddo
        do while (xkk_(2) < 0)
           xkk_(2) = xkk_(2) + 1
        enddo
        do while (xkk_(2) >= 1)
           xkk_(2) = xkk_(2) - 1
        enddo
        do while (xkk_(3) < 0)
           xkk_(3) = xkk_(3) + 1
        enddo
        do while (xkk_(3) >= 1)
           xkk_(3) = xkk_(3) - 1
        enddo
        ijk_k(1) = nint(xkk_(1)*nkf1)
        ijk_k(2) = nint(xkk_(2)*nkf2)
        ijk_k(3) = nint(xkk_(3)*nkf3)
        if (.not. epcheck) then
           ik_red_ = map_mesh_rful(ijk_k(1)*nkf2*nkf3+ijk_k(2)*nkf3+ijk_k(3)+1)
        else
           ik_red_ = map_mesh_ful(ijk_k(1)*nkf2*nkf3+ijk_k(2)*nkf3+ijk_k(3)+1)
        endif

        if (ik_red_ /= 0) then
           evc_fk(:,:,:) = evc_fmesh(:,:,:,ik_red_)
        else
           CALL errore('ephwann_shuffle','evc_fk not found',1)
        endif
        if (.not. epcheck) then
           if (rful2rirr(ik_red_) /= ik_red) &
              CALL errore('ephwann_shuffle','rful2rirr(ik_ful_red) /= ik_irr_red',1)
        else
           if (rful2rirr(ful2rful(ik_red_)) /= ik_red) &
              CALL errore('ephwann_shuffle','rful2rirr(ful2rful(ik_ful)) /= ik_irr_red',1)
        endif
        !
        if (.false. .and. (gamma_find .or. x_find .or. l_find)) then

           ibnd = 10
           write(*,*)
           write(*,*) ' xk =',xkk(1:3)
           write(*,*) ' check evc_fk:'
           write(*,*) 1, evc_fk(1,1,ibnd)
           write(*,*) 2, evc_fk(2,1,ibnd)
           write(*,*) 3, evc_fk(3,1,ibnd)
           write(*,*) 4, evc_fk(4,1,ibnd)
           write(*,*) 5, evc_fk(5,1,ibnd)

        endif
        !
     endif
     !
     ! ======================================================
     !  velocity: Wannier -> Bloch
     ! ======================================================
     CALL CPU_TIME (t4i)
     !
     IF (vme) THEN
        IF (eig_read) THEN
           CALL vmewan2bloch (nbndsub, nrr_k, irvec, ndegen_k, xkk, cufkk, vmef(:,:,:,1), etf(:,1), etf_ks(:,1), chw_ks)
        ELSE
           CALL vmewan2bloch (nbndsub, nrr_k, irvec, ndegen_k, xkk, cufkk, vmef(:,:,:,1), etf(:,1), etf_ks(:,1), chw)
        ENDIF
     ENDIF 
     CALL CPU_TIME (t4f)
     t4 = t4f - t4i
     !
     !
     ! ======================================================
     ! epmat : Wannier el and Wannier ph -> Bloch el and Wannier ph
     ! ======================================================
     CALL CPU_TIME (t5i)
     !
     CALL ephwan2bloche (nmodes, xkk, irvec, ndegen_k, nrr_q, cufkk, epmatwef, nbndsub, nrr_k)
     !
     if (eimp_mode == 5 .or. eimp_mode == 6) then
        !
        CALL eimpwan2bloche (xkk, irvec, ndegen_k, nrr_q, cufkk, eimpmatwef, nbndsub, nrr_k)
        !
     endif
     !
     CALL CPU_TIME (t5f)
     t5 = t5f - t5i
     !
     !
     ! ======================================================
     !                         q loop
     ! ======================================================
     !
     t6 = 0.0d0
     t7 = 0.0d0
     t8 = 0.0d0
     t9 = 0.0d0
     !
     !
     IF (save_m_mat) THEN
!ERROR: save_m_mat for [epmatf] and [epf17] are wrong (recl length not correct)
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
        if (eimp_mode == 5 .or. eimp_mode == 6) then
           eimpf_sr = (0.0d0,0.0d0)
        endif
        !
        if (eimp_mode == 7 .or. eimp_mode == 8) then
           eimpf_full = (0.0d0,0.0d0)
        endif
        !
     ENDIF
     !
     IF (save_m_ph) THEN
        !
        uf_ufmt  = TRIM(tmp_dir)//TRIM(prefix)//'.phu'
        OPEN (17817,FILE=uf_ufmt,FORM='unformatted',ACCESS='direct',RECL=2*nmodes*nmodes*DP,STATUS='old')
        !
     ENDIF
     !
     IF (save_t_el) THEN
        !
        etf_ufmt  = TRIM(tmp_dir)//TRIM(prefix)//'.ele'
        cuf_ufmt  = TRIM(tmp_dir)//TRIM(prefix)//'.elc'
        OPEN (40204,FILE=etf_ufmt,FORM='unformatted',ACCESS='direct',RECL=2*nbndsub*DP,STATUS='old')
        OPEN (40504,FILE=cuf_ufmt,FORM='unformatted',ACCESS='direct',RECL=2*nbndsub*nbndsub*DP,STATUS='old')
        !
     ENDIF
     !
     !
     DO iq = 1, nq_ful
        !
        ikq = 2 * iq
        ikk = ikq - 1  
        ! 
        xxq(:) = xqf_ful(:,iq) ! xqf_ful is assumed to be in crys coord
!        xxq = -xxq
        xkq = xkk + xxq
        !
        if (.not.screen_polar) then
           IF (save_m_ph) THEN
              READ (17817,REC=iq) uf(1:nmodes,1:nmodes,1,1)
           ELSE
              uf(:,:,1,1) = uf_ful(:,:,1,1,iq)
           ENDIF
           do itemp = 1, neptemp
              do idope = 1, nepdope
                 uf(:,:,itemp,idope) = uf(:,:,1,1)
              enddo
           enddo
        else
           IF (save_m_ph) THEN
              READ (17817,REC=iq) uf(1:nmodes,1:nmodes,1:neptemp,1:nepdope)
           ELSE
              uf(:,:,:,:) = uf_ful(:,:,:,:,iq)
           ENDIF
        endif
        !
        IF (eig_read) etf_ks (:,ikk) = etf_ks(:,1)
        etf (:,ikk) = etf(:,1)
        IF (eimp_mode == 2 .or. eimp_mode == 4 .or. eimp_mode == 6 .or. eimp_mode == 8) then
           !dmef(:,:,:,ikk) = dmef(:,:,:,1)
           vkq(:,:,ikk) = vkq(:,:,1)
        ENDIF
        IF (vme) vmef(:,:,:,ikk) = vmef(:,:,:,1) 
        !
!        if (iq > 1) cycle
        !
        ! ======================================================
        ! hamiltonian : Wannier -> Bloch
        ! ======================================================
        CALL CPU_TIME (t6i)
        !
        IF (.NOT. save_t_el) THEN
           !
           IF (eig_read) THEN    
              CALL hamwan2bloch (nbndsub, nrr_k, irvec, ndegen_k, xkq, cufkq, etf_ks (:,ikq), chw_ks)
           ENDIF
           !   
           CALL hamwan2bloch (nbndsub, nrr_k, irvec, ndegen_k, xkq, cufkq, etf (:,ikq), chw)
           !
           !DO ibnd = 1, nbndsub
           !   IF (etf(ibnd,ikq) .GT. ef_m) etf(ibnd,ikq) = etf(ibnd,ikq) + delta_egap
           !ENDDO
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
        ! ======================================================
        !  dipole: Wannier -> Bloch
        ! ======================================================
        CALL CPU_TIME (t7i)
        !
        IF (eimp_mode == 2 .or. eimp_mode == 4 .or. eimp_mode == 6 .or. eimp_mode == 8) THEN
           !
           IF (vg_el .EQ. 'matrix') THEN
              !
              CALL dmewan2bloch (nbndsub, nrr_k, irvec, ndegen_k, xkq, cufkq, dmef_tmp, etf(:,ikq), etf_ks(:,ikq))
              !
              DO ibnd = 1, nbndsub
                 vkq(:,ibnd,ikq) = 2.0d0*REAL(dmef_tmp(:,ibnd,ibnd))
              ENDDO
              !
           ELSEIF (vg_el .EQ. 'linear') THEN
              xxq_ = 0.d0
!              write(stdout,*) ' vg_el: ', vg_el
!              CALL velwan2bloch (xkq, xxq_, etf(:,ikq), chw, nrr_k, irvec, ndegen_k, velf)

              CALL vgwan2bloch_k (xkq, 0.0d0, velf)
              IF (epdim .EQ. 2) velf(3,:) = 0.0d0
              !
              DO ibnd = 1, nbndsub
                 vkq(:,ibnd,ikq) = velf(:,ibnd)
              ENDDO
              !
           ENDIF
        ENDIF
        !
        CALL CPU_TIME (t7f)
        t7 = t7 + (t7f-t7i)
        !
        !
        ! ======================================================
        !  velocity: Wannier -> Bloch
        ! ======================================================
        CALL CPU_TIME (t8i)
        !
        IF (vme) THEN
           IF (eig_read) THEN
              CALL vmewan2bloch (nbndsub, nrr_k, irvec, ndegen_k, xkq, cufkq, vmef (:,:,:, ikq), etf(:,ikq), etf_ks(:,ikq), chw_ks)
           ELSE
              CALL vmewan2bloch (nbndsub, nrr_k, irvec, ndegen_k, xkq, cufkq, vmef (:,:,:, ikq), etf(:,ikq), etf_ks(:,ikq), chw)
           ENDIF
        ENDIF 
        !
        CALL CPU_TIME (t8f)
        t8 = t8 + (t8f-t8i)
        !
        !
        ! check the scattering event is within_range
        within_range = .FALSE.
        !
        DO ibnd = ibndmin, ibndmax
           DO jbnd = ibndmin, ibndmax
              !
              IF (nptype .EQ. 'n') THEN
                 IF ( (etf(ibnd,ikk) .GE. cbnd_emin .AND. etf(ibnd,ikk) .LE. cbnd_emin+cfsthick+epthick) .and. &
                      (etf(jbnd,ikq) .GE. cbnd_emin .AND. etf(jbnd,ikq) .LE. cbnd_emin+cfsthick+epthick) ) within_range = .TRUE.
              ELSEIF (nptype .EQ. 'p') THEN
                 IF ( (etf(ibnd,ikk) .GE. vbnd_emax-vfsthick-epthick .AND. etf(ibnd,ikk) .LE. vbnd_emax) .and. &
                      (etf(jbnd,ikq) .GE. vbnd_emax-vfsthick-epthick .AND. etf(jbnd,ikq) .LE. vbnd_emax) ) within_range = .TRUE.
              ELSE
                 IF ( (etf(ibnd,ikk) .GE. vbnd_emax-vfsthick .AND. etf(ibnd,ikk) .LE. cbnd_emin+cfsthick) .AND. &
                      (etf(jbnd,ikq) .GE. vbnd_emax-vfsthick-epthick .AND. etf(jbnd,ikq) .LE. cbnd_emin+cfsthick+epthick) ) within_range = .TRUE.
              ENDIF

              !
           ENDDO
        ENDDO    
        ! 
        ! when ep-check, always consider every e-ph coupling
        if (epcheck) within_range = .true.
        !
        ! make sure the corresponding irreducible point of xkq is also within range
        if (eimp_mode == 7 .or. eimp_mode == 8) then
           !
           xkq_ = xkq
           do while (xkq_(1) < 0)
              xkq_(1) = xkq_(1) + 1
           enddo
           do while (xkq_(1) >= 1)
              xkq_(1) = xkq_(1) - 1
           enddo
           do while (xkq_(2) < 0)
              xkq_(2) = xkq_(2) + 1
           enddo
           do while (xkq_(2) >= 1)
              xkq_(2) = xkq_(2) - 1
           enddo
           do while (xkq_(3) < 0)
              xkq_(3) = xkq_(3) + 1
           enddo
           do while (xkq_(3) >= 1)
              xkq_(3) = xkq_(3) - 1
           enddo
           ijk_kq(1) = nint(xkq_(1)*nkf1)
           ijk_kq(2) = nint(xkq_(2)*nkf2)
           ijk_kq(3) = nint(xkq_(3)*nkf3)
           ikq_ful = map_mesh_ful(ijk_kq(1)*nkf2*nkf3+ijk_kq(2)*nkf3+ijk_kq(3)+1)

           if (.not. epcheck) then
              if (ful2rful(ikq_ful) == 0) then
!                 write(*,*) ' energy of kq point ignored:', etf(ibndmin:ibndmax,ikq)*ryd2ev
                 within_range = .false.
              endif
           endif
           !
        endif
        !
        CALL CPU_TIME (t9_1)

        IF (within_range) THEN
           !
           ! wavefunction at k+q

           if (eimp_mode == 7 .or. eimp_mode == 8) then
              !
!              call evcwan2bloch (xkq, irvec, ndegen_k, cufkq, nbndsub, nrr_k, evc_fkq, ig_max_red)

              if (.not. epcheck) then
                 ikq_red = map_mesh_rful(ijk_kq(1)*nkf2*nkf3+ijk_kq(2)*nkf3+ijk_kq(3)+1)
              else
                 ikq_red = map_mesh_ful(ijk_kq(1)*nkf2*nkf3+ijk_kq(2)*nkf3+ijk_kq(3)+1)
              endif

              if (ikq_red /= 0) then
                 evc_fkq(:,:,:) = evc_fmesh(:,:,:,ikq_red)
!                 write(stdout,*) ' ikq_red =', ikq_red
!                 write(stdout,*) ' xkq =', xkq
!                 write(stdout,*) ' xk(ikq_red) =', xkf_ful(:,ikq_red)
!                 write(stdout,*) ' ijk_kq =', ijk_kq(1:3)
!                 write(stdout,*) ' map_mesh_ful =', map_mesh_ful(1:3)
              else
                 write(*,*) ' ERROR. index of kq point:', ijk_kq(1)*nkf2*nkf3+ijk_kq(2)*nkf3+ijk_kq(3)+1
                 write(*,*) ' kq point:', xkq_(1:3)
                 write(*,*) ' energy of kq point:', etf(ibndmin:ibndmax,ikq)
                 CALL errore('ephwann_shuffle','evc_fkq not found but within range',1)
              endif
              !
           endif
           !
           CALL CPU_TIME (t9_2)
           ! ======================================================
           ! epmat : Bloch el and Wannier ph -> Bloch el and Bloch ph
           ! ======================================================
           CALL CPU_TIME (t9i)
           !
!CHECK: need to add e_ph_interp subroutine
           if (eph_interp) then
              !
           else
              if (.not. screen_polar) then
                 CALL ephwan2bloch2 (nbndsub, nrr_q, irvec, ndegen_q, epmatwef, xxq, uf, cufkk, cufkq, epmatf(:,:,:,1,1), nmodes,1,1)
              else
                 CALL ephwan2bloch2 (nbndsub, nrr_q, irvec, ndegen_q, epmatwef, xxq, uf, cufkk, cufkq, epmatf(:,:,:,:,:), nmodes,neptemp,nepdope)
              endif
           endif

           !  when tranferring from Wanneir to Bloch
           !  for the case of not using sreen_polar we only change epmatf(:,:,:,1,1)
           !  here we need to copy
           !
           if (eph_interp .or. (.not. screen_polar)) then
              do itemp = 1, neptemp
                 do idope = 1, nepdope
                    if ((itemp /= 1) .or. (idope /= 1)) &
                       epmatf(:,:,:,itemp,idope) = epmatf(:,:,:,1,1)
                 enddo
              enddo
           endif

           ! 
           CALL CPU_TIME (t9_3)

           IF (lpolar .and. elop) THEN
              !
              CALL compute_bmn_para2 (nbndsub, nkstot, cufkk, cufkq, bmatf)
              !
              IF ( (abs(xxq(1)).gt.eps) .or. (abs(xxq(2)).gt.eps) .or. (abs(xxq(3)).gt.eps) ) THEN
                 CALL cryst_to_cart (1, xxq, bg, 1)
                 DO ibnd = 1, nbndsub
                    DO jbnd = 1, nbndsub
                       !CALL rgd_blk_epw2(nq1, nq2, nq3, xxq, uf, epmatf(ibnd,jbnd,:), &
                       !                  nmodes, epsi, zstar, bmatf(ibnd,jbnd), +1.d0)

                       if (.not. screen_polar) then
                          CALL rgd_blk_epw3(uf, epmatf(ibnd,jbnd,:,:,:), neptemp, nepdope, eph_vogl(:,:,:,iq), &
                                            nmodes, bmatf(ibnd,jbnd), +1.d0)
                       else
                          CALL rgd_blk_epw3(uf, epmatf(ibnd,jbnd,:,:,:), neptemp, nepdope, eph_vogl(:,:,:,iq), &
                                            nmodes, bmatf(ibnd,jbnd), +1.d0)
                       endif
                    ENDDO
                 ENDDO
                 CALL cryst_to_cart (1, xxq, at, -1)
              ENDIF
              !
           ENDIF
           !
           if (eimp_mode == 5 .or. eimp_mode == 6) then
              !
              CALL eimpwan2bloch2 (nbndsub, nrr_q, irvec, ndegen_q, eimpmatwef, xxq, cufkk, cufkq, eimpmatf)
              ! 
           endif
           !
           if (eimp_mode == 7 .or. eimp_mode == 8) then
              !
              ! Fourier component of electron wavefunction products
              !
              eemat_f = 0

              do ibnd = ibndmin, ibndmax
                 ibnd0 = ibnd-ibndmin+1
                 !
                 ind_g = 0
                 do igz = -bg_max_red, bg_max_red
                    do igy = -bg_max_red, bg_max_red
                       do igx = -bg_max_red, bg_max_red
                          !
                          ind_g = ind_g + 1
                          ig_shift = igz*dffts%nr2*dffts%nr1 + igy*dffts%nr1 + igx
                          !
                          evc_remap = 0.d0
                          do ig = 1, ig_max_red
                             if ((map_kq2k(ig,ind_g) >= 1) .and. (map_kq2k(ig,ind_g) <= ig_max_red)) then
                                evc_remap(ig,1) = evc_fk(map_kq2k(ig,ind_g),1,ibnd0)
                                if (noncolin) &
                                   evc_remap(ig,2) = evc_fk(map_kq2k(ig,ind_g),2,ibnd0)
                             endif
                          enddo
                          !
                          do jbnd = ibndmin, ibndmax
                             jbnd0 = jbnd-ibndmin+1
                             eemat_f(jbnd0,ibnd0,1,ind_g) = &   ! use reduced number of G points
                                    ZDOTC (ig_max_red, evc_fkq(1:ig_max_red,1,jbnd0), 1, evc_remap(1:ig_max_red,1), 1)
                             if (noncolin) then
                                eemat_f(jbnd0,ibnd0,2,ind_g) = &   ! use reduced number of G points
                                       ZDOTC (ig_max_red, evc_fkq(1:ig_max_red,2,jbnd0), 1, evc_remap(1:ig_max_red,2), 1)
                             endif
                          enddo
                          !
                       enddo
                    enddo
                 enddo
              enddo
              !
              CALL CPU_TIME (t9_4)

              ! interpolate Fourier component of impurity potential
              ! here the q-point used should correspond to the (k+q)-(k) used by wavefunctions
              xxq_int = xxq
              xxq_int(1) = xxq_int(1) - nint(xxq_int(1))
              xxq_int(2) = xxq_int(2) - nint(xxq_int(2))
              xxq_int(3) = xxq_int(3) - nint(xxq_int(3))
              !
              if (eimp_ls_mode == 0 .or. eimp_ls_mode == 2) then
                 call eimp_interp_tetra(xxq_int, bg_max_red, ng_max_red, dvimpqf_sr)
              else
                 dvimpqf_sr = 0.d0
              endif

              CALL CPU_TIME (t9_5)
              !
              ! calculate electron-impurity scattering matrix - short range plus long range
              !
              ! [IMPORTANT NOTE] this brings xq to first-BZ, but it is only centered around Gamma, not yet fully symmetric
              xq_cart(1) = xxq(1) - NINT(xxq(1))
              xq_cart(2) = xxq(2) - NINT(xxq(2))
              xq_cart(3) = xxq(3) - NINT(xxq(3))
              CALL cryst_to_cart ( 1, xq_cart, bg, 1 )   ! xq_cart in [2*pi/alat], [alat] being primitive unit cell lattice constant
 
              ! here the mapping should be established from xxq to the coarse mesh
!              xq_cart(1:3) = xxq(1:3)
!              CALL cryst_to_cart ( 1, xq_cart, bg, 1 )   ! xq_cart in [2*pi/alat]
              !
              ind_g = 0
              eimpmat_f = (0.d0, 0.d0)
              if (eimp_ls_mode < 2) then
                 imp_charge_ = imp_charge
              else
                 imp_charge_ = 0.d0
              endif
              do igz = -bg_max_red, bg_max_red
                 do igy = -bg_max_red, bg_max_red
                    do igx = -bg_max_red, bg_max_red
                       !
                       ind_g = ind_g + 1
                       G_cart = igx*bg(1:3,1) + igy*bg(1:3,2) + igz*bg(1:3,3)
                       xqG = xq_cart + G_cart
                       xqGl  = sqrt(xqG(1)**2 + xqG(2)**2 + xqG(3)**2) * (2.d0*pi/alat)
                       !
                       do idope = 1, nepdope
                          do itemp = 1, neptemp
                             ! note that the factor of 2*pi/hbar (hbar=1 in AU) is added later,
                             ! specifically, pi is added in [selfen_elec.f90], while 2 is considered in later
                             ! transport calculations
                             dvimpqf_lr(itemp,idope) = - imp_charge_*(2.d0*(4.d0*pi)/omega/dielec0) / &
                                         (xqGl**2.d0 + 1.d0/(L_D(itemp,idope)**2.d0))
 
                             eimpmat_f(:,:,itemp,idope) = eimpmat_f(:,:,itemp,idope) + &
                                (dvimpqf_sr(ind_g) + dvimpqf_lr(itemp,idope)) * eemat_f(:,:,1,ind_g)
                             if (noncolin) then
                                eimpmat_f(:,:,itemp,idope) = eimpmat_f(:,:,itemp,idope) + &
                                   (dvimpqf_sr(ind_g) + dvimpqf_lr(itemp,idope)) * eemat_f(:,:,2,ind_g)
                             endif
                          enddo
                       enddo
                    enddo
                 enddo
              enddo
              !
           endif
           !
!           if ((epcheck) .and. (filqf /= '')) then
           if (epcheck) then
              xq_cart = xxq
              call cryst_to_cart (1, xq_cart, bg, 1)

              IF (.NOT. lspinorb) THEN
                 mbnd = NINT(nelec)/2 + 1
              ELSE
                 mbnd = NINT(nelec) + 1
              ENDIF
              !
              if (cond_find) then
                 !
                 write(*,*) ' ------------------------------------------------ '
                 write(*,*) ' print e-ph matrix at #pool:',my_pool_id, iq
                 write(*,*) ' check xq_cart :', xq_cart

                 ! test impurity matrix element
                 eimpmat_r = (0.d0, 0.d0)
                 do jbnd = ibndmin, ibndmax
                    jbnd0 = jbnd-ibndmin+1
                    do ibnd = ibndmin, ibndmax
                       ibnd0 = ibnd-ibndmin+1
                       !
                       ind_g = 0
                       do igz = -bg_max_red, bg_max_red
                          do igy = -bg_max_red, bg_max_red
                             do igx = -bg_max_red, bg_max_red
                                !
                                ind_g = ind_g + 1
                                eimpmat_r (jbnd,ibnd) = eimpmat_r (jbnd,ibnd) + &
!                                        dvimp_q(1,igx+bg_max+1,igy+bg_max+1,igz+bg_max+1) * eemat_f(jbnd0,ibnd0,ind_g,iq)
                                        dvimpqf_sr(ind_g) * eemat_f(jbnd0,ibnd0,1,ind_g)
                             enddo
                          enddo
                       enddo
                       !
!                       do ind_g = 1, ng_max_red
!                          eimpmat_r (jbnd,ibnd) = eimpmat_r (jbnd,ibnd) + &
!                                  dvimpqf_sr(ind_g) * eemat_f(jbnd0,ibnd0,ind_g)
!                       enddo
                    enddo
                 enddo
           
                 write(*,*) ' check evc_fk:'
                 norm_evcfk = 0.d0
                 do ig = 1, ig_max_red
                    norm_evcfk = norm_evcfk + abs(evc_fk(ig,1,mbnd-ibndmin+1))**2.d0
                 enddo
                 write(*,*) ' norm_evcfk =', norm_evcfk

                 write(*,*)
                 write(*,*) ' ig_max_red =', ig_max_red
                 write(*,*) ' check eemat_f at G=0:'
                 ind_g = 0
                 do igz = -bg_max_red, bg_max_red
                    do igy = -bg_max_red, bg_max_red
                       do igx = -bg_max_red, bg_max_red
                          !
                          ind_g = ind_g + 1
                          if ((igy == 0) .and. (igz == 0)) then
                             if (cond_find) then
                                write(*,*) igz, igy, igx, eemat_f(mbnd-ibndmin+1,mbnd-ibndmin+1,1,ind_g)
                             elseif (vale_find) then
                                write(*,*) igz, igy, igx, eemat_f(mbnd-ibndmin,mbnd-ibndmin,1,ind_g)
                             endif
                          endif
                       enddo
                    enddo
                 enddo
                 write(*,*) ' check dvimp_q at G=0:', dvimp_q(1,bg_max+1,bg_max+1,bg_max+1)
                 write(*,*) ' check dvimpqf_sr at G=0:', dvimpqf_sr(bg_max_red*((2*bg_max_red+1)**2)+bg_max_red*(2*bg_max_red+1)+bg_max_red+1)
      
                 write(*,*)
                 write(*,*) ' xxq = ', xxq(:)
                 !write(*,*) ' eimpmat_r at CBM:', eimpmat_r(mbnd,mbnd)
                 write(*,*) ' mbnd =', mbnd
                 if (cond_find) then
                    write(*,*) ' eimpmat_r at CBM:', eimpmat_r(mbnd,mbnd)
                 elseif (vale_find) then
                    write(*,*) ' eimpmat_r at VBM:', eimpmat_r(mbnd-3,mbnd-3), &
                                                     eimpmat_r(mbnd-2,mbnd-2), eimpmat_r(mbnd-1,mbnd-1)
                 endif
                 !
              endif
              !
              if (cond_find) then
                 do ibnd_c = mbnd, mbnd+nbnd_c-1
                    do jbnd_c = mbnd, mbnd+nbnd_c-1
                       write(30401,'(3(1x,i7),5(1x,f15.7),100(1x,E15.7))') iq, ibnd_c, jbnd_c, xq_cart, &
                                          etf(ibnd_c,2*iq-1), etf(jbnd_c,2*iq), &
                                          wf_ful(1:nmodes,1,1,iq), &
                                          abs(epmatf(ibnd_c,jbnd_c,1:nmodes,1,1))
                       if (eimp_mode == 5 .or. eimp_mode == 6) then
                          write(30405,'(3(1x,i7),5(1x,f15.7),100(1x,E15.7))') iq, ibnd_c, jbnd_c, xq_cart, &
                                          etf(ibnd_c,2*iq-1), etf(jbnd_c,2*iq), &
                                          real(eimpmatf(ibnd_c,jbnd_c)), aimag(eimpmatf(ibnd_c,jbnd_c))
                       endif
                    enddo
                 enddo
              endif
              !
              if (vale_find) then
                 do ibnd_v = mbnd-nbnd_v, mbnd-1
                    do jbnd_v = mbnd-nbnd_v, mbnd-1
                       write(30403,'(3(1x,i7),5(1x,f15.7),100(1x,E14.7))') iq, ibnd_v, jbnd_v, xq_cart, &
                                          etf(ibnd_v,2*iq-1), etf(jbnd_v,2*iq), &
                                          wf_ful(1:nmodes,1,1,iq), &
                                          abs(epmatf(ibnd_v,jbnd_v,1:nmodes,1,1))
                       if (eimp_mode == 5 .or. eimp_mode == 6) then
                          write(30405,'(3(1x,i7),5(1x,f15.7),100(1x,E15.7))') iq, ibnd_v, jbnd_v, xq_cart, &
                                          etf(ibnd_v,2*iq-1), etf(jbnd_v,2*iq), &
                                          real(eimpmatf(ibnd_v,jbnd_v)), aimag(eimpmatf(ibnd_v,jbnd_v))
                       endif
                    enddo
                 enddo
              endif
           endif

           !    
           ! ignore some phonon mode in scattering
           IF (ignph) THEN
              !
              DO imode = 1, ignph_num
                 do itemp = 1, neptemp
                    do idope = 1, nepdope
                       epmatf(:,:,ignph_mode(imode),itemp,idope) = (0.0d0,0.0d0)
                    enddo
                 enddo
              ENDDO
              !
           ENDIF
           ! 
           ! write epmatf to file / store in memory
           IF (save_m_mat) THEN
              !
              !ERROR: writing epmatf is not correct
              WRITE (17017,REC=iq) epmatf(ibndmin:ibndmax,ibndmin:ibndmax,1:nmodes,1,1)
              !
           ELSE
              !
              DO jbnd = ibndmin, ibndmax
                 DO ibnd = ibndmin, ibndmax
                    ibnd0 = ibnd-ibndmin+1
                    jbnd0 = jbnd-ibndmin+1
                    !
                    DO imode = 1, nmodes
                       epf17(iq,jbnd0,ibnd0,imode,:,:) = epmatf(jbnd,ibnd,imode,:,:)
                    ENDDO
                    !
                    if (eimp_mode == 5 .or. eimp_mode == 6) then
                       if (alloy_pot) then
                          alelpp(iq,jbnd0,ibnd0) = eimpmatf(jbnd,ibnd)
                       else
                          eimpf_sr(iq,jbnd0,ibnd0) = eimpmatf(jbnd,ibnd)
                       endif
                    endif
                    !
                    if (eimp_mode == 7 .or. eimp_mode == 8) then
                       eimpf_full(:,:,:,:,iq) = eimpmat_f(:,:,:,:)
                    endif
                    !
                 ENDDO
              ENDDO
              !
           ENDIF
           !
           CALL CPU_TIME (t9f)
           t9 = t9 + (t9f-t9i)
           ! 
        ENDIF ! within_range
        !
        CALL CPU_TIME (t9_6)

        if (within_range .and. .false.) then
           write(stdout,*) ' time test, iq#:', iq
           write(stdout,*) ' total: t9_1 -> t9_6:', t9_6 - t9_1
           write(stdout,*) ' t9_1 -> t9_2:', t9_2 - t9_1
           write(stdout,*) ' t9_2 -> t9_3:', t9_3 - t9_2
           write(stdout,*) ' t9_3 -> t9_4:', t9_4 - t9_3
           write(stdout,*) ' t9_4 -> t9_5:', t9_5 - t9_4
           write(stdout,*) ' t9_5 -> t9_6:', t9_6 - t9_5
        endif
        !
     ENDDO ! q loop 
     !
     ! ======================================================
     ! compute short-range coupling matrix
     ! ======================================================
     !     
     if (eimp_mode > 2) then
     endif
     !
     ! ======================================================
     ! electron self-energy
     ! ======================================================
     CALL CPU_TIME (tai)
     !
     CALL selfen_elec (ik_red, ik_star, xkk, nscat)
     !
     IF (bte .EQ. 0 .OR. bte .EQ. 1 .OR. bte .EQ. 3) THEN
        nscat_all(ik_red) = nscat
     ENDIF
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
     WRITE (stdout,'(/5x,a,3f10.6,a,i4,a,i4,a)') 'k = (', xkk(1:3), '), No. (', ik_red-ik_star+1, '/', ik_stop-ik_star+1, ')'
     WRITE (stdout,'(/13x,a,f8.2,a)') 'In  k | Out q | ham     W->B :', t2, 's'
     WRITE (stdout,'(13x,a,f8.2,a)')  'In  k | Out q | dme     W->B :', t3, 's'
     WRITE (stdout,'(13x,a,f8.2,a)')  'In  k | Out q | vme     W->B :', t4, 's'
     WRITE (stdout,'(13x,a,f8.2,a)')  'In  k | Out q | eph(el) W->B :', t5, 's'
     WRITE (stdout,'(13x,a,f8.2,a)')  'In  k | In  q | ham     W->B :', t6, 's'
     WRITE (stdout,'(13x,a,f8.2,a)')  'In  k | In  q | dme     W->B :', t7, 's'
     WRITE (stdout,'(13x,a,f8.2,a)')  'In  k | In  q | vme     W->B :', t8, 's'
     WRITE (stdout,'(13x,a,f8.2,a)')  'In  k | In  q | eph(ph) W->B :', t9, 's'
     WRITE (stdout,'(13x,a,f8.2,a)')  'In  k | Out q | self-energy  :', ta, 's'
     WRITE (stdout,'(13x,a,f8.2,a,i2,a,i2,a,i2)')  '                TOTAL        :', t0, 's | ', values_(5), ':', values_(6), ':', values_(7)
     !
     !
     IF (save_m_mat) CLOSE (17017)
     IF (save_m_ph)  CLOSE (17817)
     IF (save_t_el) THEN
        CLOSE (40204)
        CLOSE (40504)
     ENDIF
     !
  ENDDO ! k loop 
  ! 
  !
  if ((epcheck) .and. (filqf/='')) then
     close(30401)
     close(30403)
     if (eimp_mode == 5 .or. eimp_mode == 6) then
        close(30405)
        close(30406)
     endif
  endif
  !
  !
#ifdef __PARA
  CALL mp_barrier (inter_pool_comm)
  CALL mp_sum (sigmai_mode_all_abs,inter_pool_comm)
  CALL mp_sum (sigmai_mode_all_emi,inter_pool_comm)
  CALL mp_sum (sigmai_mode_all_inter,inter_pool_comm)
  CALL mp_sum (sigmai_mode_all_intra,inter_pool_comm)

  if (eimp_mode > 0) CALL mp_sum (sigmai_mode_all_ela_intra,inter_pool_comm)
  if (eimp_mode > 0) CALL mp_sum (sigmai_mode_all_ela_inter,inter_pool_comm)
  if (alloy_pot)  CALL mp_sum (sigmai_mode_all_alloy_inter,inter_pool_comm)
  if (alloy_pot)  CALL mp_sum (sigmai_mode_all_alloy_intra,inter_pool_comm)

  IF (bte .EQ. 0 .OR. bte .EQ. 1 .OR. bte .EQ. 3) CALL mp_sum (nscat_all,inter_pool_comm)
  IF (bte .EQ. 18) CALL mp_sum (sigmai_ch,inter_pool_comm)
#ENDIF
  !
  IF (save_m_mat) THEN
     IF (my_pool_id .EQ. ionode_id) CALL SYSTEM ('rm '//TRIM(tmp_dir)//TRIM(prefix)//'.epf17_*')
  ENDIF
  !
  IF (save_t_el .and. .NOT. phdrag) THEN
     IF (my_pool_id .EQ. ionode_id) CALL SYSTEM ('rm '//TRIM(tmp_dir)//TRIM(prefix)//'.ele')
     IF (my_pool_id .EQ. ionode_id) CALL SYSTEM ('rm '//TRIM(tmp_dir)//TRIM(prefix)//'.elc')
  ENDIF
  !
  !
  !
75301 CONTINUE
  !
  ! compute electron and phonon energy and group velocity 
  ! in irreducible BZ on fine mesh will be used in BTE calculations
  IF (bte .EQ. 0 .OR. bte .EQ. 1 .OR. bte .EQ. 3 .OR. bte .EQ. -1 .OR. bte .EQ. 19 .OR. bte .EQ. 18) THEN
     !
     ! electron
     CALL para_bounds (ik_star, ik_stop, nk_irr_red)
     !
     IF (.NOT. ALLOCATED(cufkk)) ALLOCATE (cufkk(nbndsub, nbndsub))
     ALLOCATE (etf_all(nbndsub, nk_irr_red))
     ALLOCATE (vel_all(3, nbndsub, nk_irr_red))
     etf_all = 0.0d0
     vel_all = 0.0d0
     !
     DO ik_red = ik_star, ik_stop
        !
        IF (bte .EQ. 19 .OR. bte .EQ. 18) THEN
           ik = ik_red
           xkk(:) = xk_kmesh(:,ik)
        ELSE
           ik = rirr2irr(ik_red)       ! ik is the original index of k point in irreducible BZ
           xkk(:) = xkf_irr(:,ik)     ! xkf_irr is assumed to be in crys coord
        ENDIF
        !
        IF (eig_read) CALL hamwan2bloch (nbndsub, nrr_k, irvec, ndegen_k, xkk, cufkk, etf_ks_tmp, chw_ks)     
        CALL hamwan2bloch (nbndsub, nrr_k, irvec, ndegen_k, xkk, cufkk, etf_all(:,ik_red), chw)        
        !
        !DO ibnd = 1, nbndsub
        !   IF (etf_all(ibnd,ik_red) .GT. ef_m) etf_all(ibnd,ik_red) = etf_all(ibnd,ik_red) + delta_egap
        !ENDDO
        DO ibnd = vbnd_num+1, nbndsub
           etf_all(ibnd,ik_red) = etf_all(ibnd,ik_red) + delta_egap
        ENDDO
        !        
  
        IF (vg_el .EQ. 'linear') THEN
           !
           ! compute electron velocity using linear interpolation
           CALL vgwan2bloch_k (xkk, 0.0d0, vel_all(:,:,ik_red))
           IF (epdim .EQ. 2) vel_all(3,:,ik_red) = 0.0d0
           !
        ELSEIF (vg_el .EQ. 'matrix') THEN
           !
           ! compute electron velocity using electric dipole
           CALL dmewan2bloch (nbndsub, nrr_k, irvec, ndegen_k, xkk, cufkk, dmef_tmp, etf_all(:,ik_red), etf_ks_tmp)
           DO ibnd = 1, nbndsub
              vel_all(:,ibnd,ik_red) = 2.0d0*REAL(dmef_tmp(:,ibnd,ibnd))
              IF (epdim .EQ. 2) vel_all(3,ibnd,ik_red) = 0.0d0
           ENDDO
           !
        ENDIF
        !
     ENDDO
     !
#ifdef __PARA
     CALL mp_barrier (inter_pool_comm)
     CALL mp_sum (etf_all,inter_pool_comm)
     CALL mp_sum (vel_all,inter_pool_comm)
#ENDIF
     !
     DEALLOCATE (cufkk)
     !
  ENDIF
  !
  !
  IF (phdrag) CALL phdrag_shuffle (nqc, xqc)
  IF (alloy_read) CALL read_alloy_el ()
  !
  !
  CALL meta_save ()
  !
  !
  !
  !
  !
  !
  !
  !
75302 CONTINUE ! for bte=2
  !
  IF (bte .EQ. 2 .OR. bte .EQ. 29) THEN
     !
     !
     IF (bte .EQ. 29) THEN
        !
        WRITE (stdout,'(//5x,a)') '==================================================================='
        WRITE (stdout,'(5x,a)')   '                     Random q-point calculation                    '
        WRITE (stdout,'(5x,a/)')  '==================================================================='
        !
        ! check exisitence of file
        filename_qmesh = TRIM(prefix)//'.qmesh'
        INQUIRE (FILE=filename_qmesh,EXIST=file_exist)
        IF (.NOT. file_exist) CALL errore('ephwann_shuffle','cannot find .qmesh file (bte=29)',1)
        !
        OPEN (3029,FILE=filename_qmesh,STATUS='old')
        READ (3029,*) nq_qmesh, coord_mesh
        !
        ALLOCATE (xq_qmesh(3,nq_qmesh))
        xq_qmesh = 0.0d0
        !
        DO iq = 1, nq_qmesh
           READ (3029,*) xq_qmesh(1:3,iq)
           IF (coord_mesh .EQ. 'cart' .OR. coord_mesh .EQ. 'carte' .OR. coord_mesh .EQ. 'cartesian' .OR. &
               coord_mesh .EQ. 'Cart' .OR. coord_mesh .EQ. 'Carte' .OR. coord_mesh .EQ. 'Cartesian' .OR. &
               coord_mesh .EQ. '1') CALL cryst_to_cart (1, xq_qmesh(:,iq), at, -1)
        ENDDO 
        CLOSE (3029)
        !
        WRITE (stdout,'(5x,a,i5)') 'Number of q points in .qmesh file: ', nq_qmesh
        !
        !
        nq_ful = nq_qmesh
        nq_irr = nq_qmesh
        !
        !
#ifdef __PARA
        CALL mp_barrier (inter_pool_comm)
#ENDIF
        !
     ENDIF
     !
     !
     IF (smearing .EQ. 'tetra') THEN
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
     IF (save_m_mat .EQ. .TRUE. .OR. save_m_matw .EQ. .TRUE.) THEN
        !
        WRITE (stdout,'(5x,a/)')  'Save following variables to disk to reduce the using of internal memory :'
        !
        IF (save_m_mat) THEN
           WRITE (stdout,'(13x,a,f8.1,a)')  'epf17   ~', (DBLE(nk_ful)*1.0d-6)*DBLE(nbnd_red*nbnd_red*nmodes*2*DP), ' MB' 
        ENDIF
        !
        IF (save_m_matw) THEN
           WRITE (stdout,'(13x,a,f8.1,a)')  'epmatwe ~', (DBLE(nrr_k*nqc)*1.0d-6)*DBLE(nbndsub*nbndsub*nmodes*2*DP), ' MB' 
           WRITE (stdout,'(13x,a,f8.1,a)')  'epmatwp ~', (DBLE(nrr_k*nrr_q)*1.0d-6)*DBLE(nbndsub*nbndsub*nmodes*2*DP), ' MB' 
        ENDIF
        !
     ENDIF
     !
     !
     WRITE (stdout,'(/5x,a)') 'Electron part will be export automatically when bte = 2'
     WRITE (stdout,'(5x,a/)') 'Read etf(k+q), etf_ks(k+q) and cufkk(k+q) from saved files to reduce the time consumption'
     WRITE (stdout,'(13x,a,f8.1,a)')  'etf     ~', (DBLE(nk_ful)*1.0d-6)*DBLE(nbndsub*DP), ' MB' 
     WRITE (stdout,'(13x,a,f8.1,a)')  'etf_ks  ~', (DBLE(nk_ful)*1.0d-6)*DBLE(nbndsub*DP), ' MB' 
     WRITE (stdout,'(13x,a,f8.1,a)')  'cufkk   ~', (DBLE(nk_ful)*1.0d-6)*DBLE(nbndsub*nbndsub*2*DP), ' MB' 
     !
     !
     filename_check = TRIM(prefix)//'.epcheck'
     INQUIRE (FILE=filename_check,EXIST=file_exist)
     IF (file_exist) CALL ep_check (filename_check)
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
     ENDIF
     !
     !
     ! ======================================================
     ! eph-vertex : Wannier -> Bloch (out k loop)
     ! ======================================================
     !
     IF (lpolar) THEN
        !
        ALLOCATE(eph_vogl(nmodes,neptemp,nepdope,nq_irr))
        eph_vogl = 0.0d0
        !
        WRITE (stdout,'(/5x,a,f8.2,a)') 'Prepare long-range part of e-ph matrix :'
        !
        CALL CPU_TIME (tbi)
        !
        CALL para_bounds (iq_star, iq_stop, nq_irr)
        !
        uf_cart  = 0.0d0
        DO imode = 1, nmodes
           uf_cart(imode,imode) = 1.0d0
        ENDDO
        !
        DO iq = iq_star, iq_stop
           !
           IF (bte .EQ. 29) THEN
              xxq(:) = xq_qmesh(:,iq)
           ELSE
              xxq(:) = xqf_irr(:,iq)     ! xqf_irr is assumed to be in crys coord
           ENDIF
           CALL cryst_to_cart (1, xxq, bg, 1) 
           !
           IF (MAXVAL(ABS(xxq)) .NE. 0.0d0) THEN
              CALL polar_eph (uf_cart, xxq, eph_vogl(:,:,:,iq), neptemp, nepdope)
           ELSE
              eph_vogl(:,:,:,iq) = 0.0d0
           ENDIF
           !
        ENDDO
        !
#ifdef __PARA
        CALL mp_barrier (inter_pool_comm)
        CALL mp_sum (eph_vogl,inter_pool_comm)
#ENDIF
        !
        CALL CPU_TIME (tbf)
        tb = tbf - tbi
        WRITE (stdout,'(/13x,a,f8.2,a)') 'Out k | Out q | eph-L   W->B :', tb, 's'
        !
     ENDIF
     !
     !
     ! ======================================================
     ! hamiltonian : Wannier -> Bloch (out k loop)
     ! ======================================================
     !
     WRITE (stdout,'(/5x,a,f8.2,a)') 'Prepare electron properties :'
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
        !DO ibnd = 1, nbndsub
        !   IF (etf_ham(ibnd) .GT. ef_m) etf_ham(ibnd) = etf_ham(ibnd) + delta_egap
        !ENDDO
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
     CALL mp_barrier (inter_pool_comm)
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
     mkq1 = nkf1/nqf1
     mkq2 = nkf2/nqf2
     mkq3 = nkf3/nqf3
     !
     CALL CPU_TIME (tcf)
     tc = tcf-tci
     WRITE (stdout,'(/13x,a,f8.2,a)') 'Out k | Out q | ham     W->B :', tc, 's'
     !
     CALL mp_barrier (inter_pool_comm)
     !
     !
     !
     ! globle variable
     !
     IF (ALLOCATED(cufkk)) DEALLOCATE (cufkk)
      ALLOCATE (cufkk(nbndsub, nbndsub))
     IF (ALLOCATED(cufkq)) DEALLOCATE (cufkq)
      ALLOCATE (cufkq(nbndsub, nbndsub))
     IF (ALLOCATED(epmatf)) DEALLOCATE (epmatf)
      ALLOCATE (epmatf(nbndsub, nbndsub, nmodes, neptemp, nepdope))
     IF (ALLOCATED(bmatf)) DEALLOCATE (bmatf)
      ALLOCATE (bmatf(nbndsub, nbndsub))
     IF (ALLOCATED(etf)) DEALLOCATE (etf)
      ALLOCATE (etf(nbndsub, 2*nk_ful))
     IF (ALLOCATED(etf_ks)) DEALLOCATE (etf_ks)
      ALLOCATE (etf_ks(nbndsub, 2*nk_ful)) 
     IF (ALLOCATED(epmatwef)) DEALLOCATE (epmatwef)
      ALLOCATE (epmatwef(nbndsub, nbndsub, nrr_k, nmodes))
     IF (ALLOCATED(epf17)) DEALLOCATE (epf17)
      IF (.NOT. save_m_mat) ALLOCATE (epf17 (nk_ful, nbnd_red, nbnd_red, nmodes, neptemp, nepdope)) 
     !
     !
     ! ======================================================
     !                         q loop
     ! ======================================================
     ! parallelization for nq_irr
     CALL para_bounds (iq_star, iq_stop, nq_irr)
     !
     ! BTE and selfen variable
     ALLOCATE (gammai_mode_all(neptemp,nepdope,nmodes,nq_irr))
     gammai_mode_all = 0.0d0
     !
     !
     ALLOCATE (wf_irr(nmodes,nq_irr))
     ALLOCATE (uf_irr(nmodes,nmodes,nq_irr))
     ALLOCATE (vph_irr(3,nmodes,nq_irr))
     wf_irr = 0.0d0
     uf_irr = (0.0d0,0.0d0)
     vph_irr = 0.0d0     
     !
     !
     DO iq = iq_star, iq_stop
        !
        CALL CPU_TIME (t0i)
        !
        CALL start_clock ( 'ep-interp' )
        !
        !
        IF (bte .EQ. 29) THEN
           xxq(:) = xq_qmesh(:,iq)
        ELSE
           xxq(:) = xqf_irr(:,iq)     ! xqf_irr is assumed to be in crys coord
        ENDIF
        !
        !
        ! ------------------------------------------------------
        ! dynamical matrix : Wannier -> Bloch
        ! ------------------------------------------------------
        !
        CALL CPU_TIME (t2i)
        !
        CALL dynwan2bloch (nmodes, nrr_q, irvec, ndegen_q, xxq, uf_irr(:,:,iq), wf_irr(:,iq), vph_irr(:,:,iq))
        IF (epdim .EQ. 2) vph_irr(3,:,iq) = 0.0d0
        !
        !
        DO nu = 1, nmodes
           !
           IF (wf_irr(nu,iq) .GT. 0.0d0) THEN
              wf_irr(nu,iq) = SQRT(ABS(wf_irr(nu,iq)))
           ELSE
              wf_irr(nu,iq) = -SQRT(ABS(wf_irr(nu,iq)))
           ENDIF
           !
           DO mu = 1, nmodes
              na = (mu - 1) / 3 + 1
              uf_irr(mu,nu,iq) = uf_irr(mu,nu,iq)/SQRT(amass(ityp(na)))
           ENDDO
           !
        ENDDO
        !
        IF (vg_ph .EQ. 'linear') THEN
           !
           CALL vgwan2bloch_q (0.0d0, xxq, vph_irr(:,:,iq))
           IF (epdim .EQ. 2) vph_irr(3,:,iq) = 0.0d0
           !
        ENDIF
        !
        CALL CPU_TIME (t2f)
        t2 = t2f - t2i
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
        OPEN (12943,FILE='BTE/META/xkf_ful_cryst',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='old') 
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
           READ (12943,REC=ik) xkk(1:3) ! xkf_ful is assumed to be in crys coord
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
              !DO ibnd = 1, nbndsub
              !   IF (etf(ibnd,ikq) .GT. ef_m) etf(ibnd,ikq) = etf(ibnd,ikq) + delta_egap
              !ENDDO
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
                          CALL rgd_blk_epw3(uf_irr(:,:,iq), epmatf(ibnd,jbnd,:,:,:), neptemp, nepdope, eph_vogl(:,:,:,iq), nmodes, bmatf(ibnd,jbnd), +1.d0)
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
                 ! writing epmatf is wrong
                 WRITE (17017,REC=ik) epmatf(ibndmin:ibndmax,ibndmin:ibndmax,1:nmodes,1,1)
                 !
              ELSE
                 !
                 DO jbnd = ibndmin, ibndmax
                    DO ibnd = ibndmin, ibndmax
                       DO imode = 1, nmodes
                          !
                          ibnd0 = ibnd-ibndmin+1
                          jbnd0 = jbnd-ibndmin+1
                          epf17(ik,jbnd0,ibnd0,imode,:,:) = epmatf(jbnd,ibnd,imode,:,:)
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
        CLOSE (12943)
        !
        !
        ! ======================================================
        ! phonon self-energy
        ! ======================================================
        CALL CPU_TIME (tai)
        !
        CALL selfen_phon (iq)
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
        WRITE (stdout,'(/5x,a,3f10.6,a,i4,a,i4,a)') 'k = (', xxq(1:3), '), No. (', iq-iq_star+1, '/', iq_stop-iq_star+1, ')'
        WRITE (stdout,'(/13x,a,f8.2,a)') 'In  q | Out k | dyn     W->B :', t2, 's'
        !WRITE (stdout,'(13x,a,f8.2,a)')  'In  q | Out k | dme     W->B :', t3, 's'
        !WRITE (stdout,'(13x,a,f8.2,a)')  'In  q | Out k | vme     W->B :', t4, 's'
        WRITE (stdout,'(13x,a,f8.2,a)')  'In  q | Out k | eph(ph) W->B :', t5, 's'
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
     CALL mp_sum (wf_irr,inter_pool_comm)
     CALL mp_sum (uf_irr,inter_pool_comm)
     CALL mp_sum (vph_irr,inter_pool_comm)
     CALL mp_sum (gammai_mode_all,inter_pool_comm)
#ENDIF
     !
     !
     ! Export
     CALL export_rate_ph ()
 !    IF (shengbte_read) CALL read_shengbte ()
 !    IF (shengbte_read) CALL phph_export ()
 !    IF ((shengbte_read .EQ. .TRUE.) .AND. (file_exist)) CALL phcheck_scat (filename_check)
     IF (bte .EQ. 2) CALL export_ShengBTE ()
     !
     CALL mp_barrier (inter_pool_comm)
     !
  ENDIF ! bte=2
  !
  !
  !
  CALL deallocate_shuffle ()
  !
  CALL stop_clock ( 'ephwann' )
  !
END SUBROUTINE ephwann_shuffle
!
!-------------------------------------------
SUBROUTINE epw_write
!-------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE epwcom,    ONLY : nbndsub, vme, eig_read
  USE ions_base, ONLY : nat
  USE pwcom,     ONLY : ef
  USE elph2,     ONLY : nrr_k, nrr_q, chw, rdw, cdmew, cvmew, chw_ks, zstar, epsi, epmatwp
  USE phcom,     ONLY : nmodes  
  USE io_global, ONLY : stdout
  USE io_files,  ONLY : tmp_dir, prefix
#ifdef __PARA
  USE mp,        ONLY : mp_barrier
  USE mp_global, ONLY : my_pool_id, inter_pool_comm
  USE io_global, ONLY : ionode_id
#ENDIF
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256) :: file_var, file_chw, file_cdmew, file_cvmew, file_chw_ks, file_rdw
  !
  WRITE(stdout,'(/5x,"Writing Hamiltonian, Dynamical matrix and EP vertex in Wann rep to file"/)')
  !
  IF (my_pool_id .EQ. ionode_id) THEN  
     !
     file_var = TRIM(tmp_dir)//TRIM(prefix)//'.var'
     OPEN (7701,FILE=file_var,FORM='unformatted',ACCESS='direct',RECL=4*4+DP+9*nat*DP+9*DP,STATUS='replace')
     WRITE (7701,REC=1) nbndsub, nrr_k, nmodes, nrr_q, ef, zstar(1:3,1:3,1:nat), epsi(1:3,1:3)
     CLOSE (7701)
     !
     file_chw = TRIM(tmp_dir)//TRIM(prefix)//'.chw'
     OPEN (7702,FILE=file_chw,FORM='unformatted',ACCESS='direct',RECL=2*nbndsub*nbndsub*nrr_k*DP,STATUS='replace')
     WRITE (7702,REC=1) chw(1:nbndsub,1:nbndsub,1:nrr_k)
     CLOSE (7702)
     !
     file_cdmew = TRIM(tmp_dir)//TRIM(prefix)//'.cdmew'
     OPEN (7703,FILE=file_cdmew,FORM='unformatted',ACCESS='direct',RECL=2*3*nbndsub*nbndsub*nrr_k*DP,STATUS='replace')
     WRITE (7703,REC=1) cdmew(1:3,1:nbndsub,1:nbndsub,1:nrr_k)
     CLOSE (7703)
     !
     IF (vme) THEN
        file_cvmew = TRIM(tmp_dir)//TRIM(prefix)//'.cvmew'
        OPEN (7704,FILE=file_cvmew,FORM='unformatted',ACCESS='direct',RECL=2*3*nbndsub*nbndsub*nrr_k*DP,STATUS='replace')
        WRITE (7704,REC=1) cvmew(1:3,1:nbndsub,1:nbndsub,1:nrr_k)
        CLOSE (7704)
     ENDIF
     !
     IF (eig_read) THEN
        file_chw_ks = TRIM(tmp_dir)//TRIM(prefix)//'.chw_ks'
        OPEN (7705,FILE=file_chw_ks,FORM='unformatted',ACCESS='direct',RECL=2*nbndsub*nbndsub*nrr_k*DP,STATUS='replace')
        WRITE (7705,REC=1) chw_ks(1:nbndsub,1:nbndsub,1:nrr_k)
        CLOSE (7705)
     ENDIF
     !
     file_rdw = TRIM(tmp_dir)//TRIM(prefix)//'.rdw'
     OPEN (7706,FILE=file_rdw,FORM='unformatted',ACCESS='direct',RECL=2*nmodes*nmodes*nrr_q*DP,STATUS='replace')
     WRITE (7706,REC=1) rdw(1:nmodes,1:nmodes,1:nrr_q) 
     CLOSE (7706)
     !
  ENDIF
  !
#ifdef __PARA
  CALL mp_barrier (inter_pool_comm)
#ENDIF     
  !
!---------------------------------
END SUBROUTINE epw_write
!---------------------------------
!---------------------------------
SUBROUTINE epw_read()
!---------------------------------
  USE kinds,     ONLY : DP
  USE epwcom,    ONLY : nbndsub, vme, eig_read, save_m_matw, eimp_mode
  USE pwcom,     ONLY : ef
  USE elph2,     ONLY : nrr_k, nrr_q, chw, rdw, epmatwp, cdmew, cvmew, chw_ks, zstar, epsi, irvec, &
                        eimpmatwp
  USE phcom,     ONLY : nmodes  
  USE ions_base, ONLY : nat
  USE io_global, ONLY : stdout
  USE io_files,  ONLY : tmp_dir, prefix
#ifdef __NAG
  USE f90_unix_io, ONLY : flush
#endif
#ifdef __PARA
  USE io_global, ONLY : ionode_id
  USE mp,        ONLY : mp_barrier, mp_bcast
  USE mp_global, ONLY : my_pool_id, intra_pool_comm, inter_pool_comm, root_pool
#endif
  !
  implicit none
  !
  INTEGER            :: irq
  CHARACTER(LEN=256) :: file_var, file_chw, file_cdmew, file_cvmew, file_chw_ks, file_rdw, filename_q, filename_k
  !
  WRITE(stdout,'(/5x,"Reading Hamiltonian, Dynamical matrix and EP vertex in Wann rep from file"/)')
  !
  CALL flush(6)
  !
  IF (my_pool_id .EQ. ionode_id) THEN
    !
    file_var = TRIM(tmp_dir)//TRIM(prefix)//'.var'
    OPEN (7701,FILE=file_var,FORM='unformatted',ACCESS='direct',RECL=4*4+DP+9*nat*DP+9*DP,STATUS='old')
    READ (7701,REC=1) nbndsub, nrr_k, nmodes, nrr_q, ef, zstar(1:3,1:3,1:nat), epsi(1:3,1:3)
    CLOSE (7701)
    ! 
  ENDIF
  !
#ifdef __PARA
  CALL mp_bcast (ef, ionode_id, inter_pool_comm)
  CALL mp_bcast (ef, root_pool, intra_pool_comm)
  CALL mp_bcast (nbndsub, ionode_id, inter_pool_comm)
  CALL mp_bcast (nbndsub, root_pool, intra_pool_comm)
  CALL mp_bcast (nrr_k, ionode_id, inter_pool_comm)
  CALL mp_bcast (nrr_k, root_pool, intra_pool_comm)
  CALL mp_bcast (nmodes, ionode_id, inter_pool_comm)
  CALL mp_bcast (nmodes, root_pool, intra_pool_comm)
  CALL mp_bcast (nrr_q, ionode_id, inter_pool_comm)
  CALL mp_bcast (nrr_q, root_pool, intra_pool_comm)
  CALL mp_bcast (zstar, ionode_id, inter_pool_comm)
  CALL mp_bcast (zstar, root_pool, intra_pool_comm)
  CALL mp_bcast (epsi, ionode_id, inter_pool_comm)
  CALL mp_bcast (epsi, root_pool, intra_pool_comm)
  CALL mp_barrier (inter_pool_comm)
#ENDIF
  !
  IF (.NOT. ALLOCATED(chw))                     ALLOCATE (chw(nbndsub,nbndsub,nrr_k))
  IF (.NOT. ALLOCATED(cdmew))                   ALLOCATE (cdmew(3,nbndsub,nbndsub,nrr_k))
  IF (vme .AND. (.NOT. ALLOCATED(cvmew)))       ALLOCATE (cvmew(3,nbndsub,nbndsub,nrr_k))
  IF (eig_read .AND. (.NOT. ALLOCATED(chw_ks))) ALLOCATE (chw_ks(nbndsub,nbndsub,nrr_k))
  IF (.NOT. ALLOCATED(rdw))                     ALLOCATE (rdw(nmodes,nmodes,nrr_q ))
  !
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     file_chw = TRIM(tmp_dir)//TRIM(prefix)//'.chw'
     OPEN (7702,FILE=file_chw,FORM='unformatted',ACCESS='direct',RECL=2*nbndsub*nbndsub*nrr_k*DP,STATUS='old')
     READ (7702,REC=1) chw(1:nbndsub,1:nbndsub,1:nrr_k)
     CLOSE (7702)
     !
     file_cdmew = TRIM(tmp_dir)//TRIM(prefix)//'.cdmew'
     OPEN (7703,FILE=file_cdmew,FORM='unformatted',ACCESS='direct',RECL=2*3*nbndsub*nbndsub*nrr_k*DP,STATUS='old')
     READ (7703,REC=1) cdmew(1:3,1:nbndsub,1:nbndsub,1:nrr_k)
     CLOSE (7703)
     !
     IF (vme) THEN
        file_cvmew = TRIM(tmp_dir)//TRIM(prefix)//'.cvmew'
        OPEN (7704,FILE=file_cvmew,FORM='unformatted',ACCESS='direct',RECL=2*3*nbndsub*nbndsub*nrr_k*DP,STATUS='old')
        READ (7704,REC=1) cvmew(1:3,1:nbndsub,1:nbndsub,1:nrr_k)
        CLOSE (7704)
     ENDIF
     !
     IF (eig_read) THEN
        file_chw_ks = TRIM(tmp_dir)//TRIM(prefix)//'.chw_ks'
        OPEN (7705,FILE=file_chw_ks,FORM='unformatted',ACCESS='direct',RECL=2*nbndsub*nbndsub*nrr_k*DP,STATUS='old')
        READ (7705,REC=1) chw_ks(1:nbndsub,1:nbndsub,1:nrr_k)
        CLOSE (7705)
     ENDIF
     !
     file_rdw = TRIM(tmp_dir)//TRIM(prefix)//'.rdw'
     OPEN (7706,FILE=file_rdw,FORM='unformatted',ACCESS='direct',RECL=2*nmodes*nmodes*nrr_q*DP,STATUS='old')
     READ (7706,REC=1) rdw(1:nmodes,1:nmodes,1:nrr_q) 
     CLOSE (7706)
     !
  ENDIF
  !
#ifdef __PARA
  CALL mp_bcast (chw, ionode_id, inter_pool_comm)
  CALL mp_bcast (chw, root_pool, intra_pool_comm)
  IF (eig_read) CALL mp_bcast (chw_ks, ionode_id, inter_pool_comm)
  IF (eig_read) CALL mp_bcast (chw_ks, root_pool, intra_pool_comm)
  CALL mp_bcast (rdw, ionode_id, inter_pool_comm)
  CALL mp_bcast (rdw, root_pool, intra_pool_comm)
  CALL mp_bcast (cdmew, ionode_id, inter_pool_comm)
  CALL mp_bcast (cdmew, root_pool, intra_pool_comm)
  IF (vme) CALL mp_bcast (cvmew, ionode_id, inter_pool_comm)
  IF (vme) CALL mp_bcast (cvmew, root_pool, intra_pool_comm)
  CALL mp_barrier (inter_pool_comm)
#ENDIF
  !
  !
  IF (.NOT. save_m_matw) THEN
     !
     IF (.NOT. ALLOCATED(epmatwp)) ALLOCATE (epmatwp(nbndsub,nbndsub,nrr_k,nmodes,nrr_q))
     !
     IF (my_pool_id .EQ. ionode_id) THEN
        !
        ! epmatwp(nbndsub, nbndsub, nrr_k, nmodes, nrr_q))
        ! epmatwp_q(nbndsub, nbndsub, nrr_k, nmodes)
        ! epmatwp_k(nbndsub, nbndsub, nmodes, nrr_q)
        !
        filename_q = TRIM(tmp_dir)//TRIM(prefix)//'.epmatwp_q'
        OPEN (91915,FILE=filename_q,FORM='unformatted',ACCESS='direct',RECL=2*nbndsub*nbndsub*nrr_k*nmodes*DP,STATUS='old')
        !
        DO irq = 1, nrr_q
           READ (91915,REC=irq) epmatwp(1:nbndsub,1:nbndsub,1:nrr_k,1:nmodes,irq)
        ENDDO
        !
        CLOSE (91915)
        !
     ENDIF
     !
#ifdef __PARA
    CALL mp_bcast (epmatwp, ionode_id, inter_pool_comm)
    CALL mp_bcast (epmatwp, root_pool, intra_pool_comm)
    CALL mp_barrier(inter_pool_comm)
#ENDIF
  !
  ENDIF
  !
  ! read electron-impurity coupling matrix
  !
  if (eimp_mode == 5 .or. eimp_mode == 6) then
     !
     if (.not. allocated(eimpmatwp)) ALLOCATE (eimpmatwp(nbndsub, nbndsub, nrr_k, nrr_q))
     !
     if (my_pool_id .eq. ionode_id) then
        !
        filename_q = TRIM(tmp_dir)//TRIM(prefix)//'.eimpmatwp_q'
        OPEN (91916,FILE=filename_q,FORM='unformatted',ACCESS='direct',RECL=2*nbndsub*nbndsub*nrr_k*DP,STATUS='old')
        !
        DO irq = 1, nrr_q
           READ (91916,REC=irq) eimpmatwp(1:nbndsub,1:nbndsub,1:nrr_k,irq)
        ENDDO
        !
        CLOSE (91916)
        !
     endif
     !
#ifdef __PARA
    CALL mp_bcast (eimpmatwp, ionode_id, inter_pool_comm)
    CALL mp_bcast (eimpmatwp, root_pool, intra_pool_comm)
    CALL mp_barrier(inter_pool_comm)
#ENDIF
     !
  endif
  !
  ! copy rdw to ifc
  CALL ifc_load (nrr_q, irvec)
  ! 
  !
#ifdef __PARA
  CALL mp_barrier(inter_pool_comm)
#ENDIF
  !
  WRITE(stdout,'(/5x,"Finished reading Wann rep data from file"/)')
  !
!---------------------------------
END SUBROUTINE epw_read
!---------------------------------



SUBROUTINE phonon_read ()
  !
#INCLUDE "f_defs.h"
  !
  USE kinds,     ONLY : DP
  USE phcom,     ONLY : nmodes
  USE epwcom,    ONLY : nqf1, nqf2, nqf3, smearing, bte, save_m_ph, ph_read, &
                        screen_polar, neptemp, nepdope
  USE io_files,  ONLY : tmp_dir, prefix
  USE elph2,     ONLY : uf_ful, wf_ful, vph_ful
  USE bte_var,     ONLY : nq_ful
#ifdef __PARA
  USE mp,        ONLY : mp_barrier, mp_bcast
  USE io_global, ONLY : ionode_id, stdout
  USE mp_global, ONLY : my_pool_id, inter_pool_comm
#ENDIF
  !
  IMPLICIT NONE
  !
  INTEGER             :: iq, imode, jmode, ir
  CHARACTER(LEN=256)  :: filename1, filename2, filename3
  LOGICAL             :: file_exist
  REAL(KIND=DP)       :: vph_temp(3,nmodes,neptemp,nepdope)
  !
!ERROR, vph_ful is not allocated to [neptemp,nepdope] when considering carrier screeening effect
  !
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     filename1 = TRIM(tmp_dir)//TRIM(prefix)//'.phw'
     filename2 = TRIM(tmp_dir)//TRIM(prefix)//'.phu'
     filename3 = TRIM(tmp_dir)//TRIM(prefix)//'.phv'  
     !
     ! check exisitence of file
     INQUIRE (FILE=filename1,EXIST=file_exist)
     IF (.NOT. file_exist) CALL errore('ph_read','cannot find .phw file (phonon frequency)',1)
     INQUIRE (FILE=filename2,EXIST=file_exist)
     IF (.NOT. file_exist) CALL errore('ph_read','cannot find .phu file (phonon eigenmode)',1)
     INQUIRE (FILE=filename3,EXIST=file_exist)
     IF (.NOT. file_exist) CALL errore('ph_read','cannot find .phv file (phonon velocity)',1)
     !
     if (screen_polar) then
        OPEN (99901,FILE=filename1,FORM='unformatted',ACCESS='direct',RECL=nmodes*neptemp*nepdope*DP,STATUS='old')
        OPEN (99902,FILE=filename2,FORM='unformatted',ACCESS='direct',RECL=2*nmodes*nmodes*neptemp*nepdope*DP,STATUS='old')
        OPEN (99903,FILE=filename3,FORM='unformatted',ACCESS='direct',RECL=3*nmodes*neptemp*nepdope*DP,STATUS='old')
     else
        OPEN (99901,FILE=filename1,FORM='unformatted',ACCESS='direct',RECL=nmodes*DP,STATUS='old')
        OPEN (99902,FILE=filename2,FORM='unformatted',ACCESS='direct',RECL=2*nmodes*nmodes*DP,STATUS='old')
        OPEN (99903,FILE=filename3,FORM='unformatted',ACCESS='direct',RECL=3*nmodes*DP,STATUS='old')
     endif
     !
     ! frequency
     if (screen_polar) then
        DO iq = 1, nq_ful
           READ (99901,REC=iq) wf_ful(1:nmodes,1:neptemp,1:nepdope,iq)
        ENDDO 
     else
        DO iq = 1, nq_ful
           READ (99901,REC=iq) wf_ful(1:nmodes,1,1,iq)
        ENDDO 
     endif
     !
     IF (.NOT. save_m_ph) THEN
        !
        ! eigenvector
        if (screen_polar) then
           DO iq = 1, nq_ful
              READ (99902,REC=iq) uf_ful(1:nmodes,1:nmodes,1:neptemp,1:nepdope,iq)
           ENDDO
        else
           DO iq = 1, nq_ful
              READ (99902,REC=iq) uf_ful(1:nmodes,1:nmodes,1,1,iq)
           ENDDO
        endif
        !
     ENDIF
     !
     ! velocity
     if (screen_polar) then
        DO iq = 1, nq_ful
           READ (99903,REC=iq) vph_temp(1:3,1:nmodes,1:neptemp,1:nepdope)
           vph_ful(1:3,1:nmodes,iq) = vph_temp(1:3,1:nmodes,1,1)
        ENDDO
     else
        DO iq = 1, nq_ful
           READ (99903,REC=iq) vph_ful(1:3,1:nmodes,iq)
        ENDDO
     endif
     !
     CLOSE (99901)
     CLOSE (99902)
     CLOSE (99903)
     !
     !
     IF (ph_read) WRITE (stdout,'(/13x,a,i3,a,i3,a,i3,a)') 'Load phonon properties from files (', &
                                                           nqf1, '*', nqf2, '*', nqf3, ' q-mesh)'
     !
  ENDIF
  !
#ifdef __PARA
  CALL mp_bcast (wf_ful,ionode_id,inter_pool_comm)
  CALL mp_bcast (vph_ful,ionode_id,inter_pool_comm)
  IF (.NOT. save_m_ph) CALL mp_bcast (uf_ful,ionode_id,inter_pool_comm)
  CALL mp_barrier (inter_pool_comm)
#ENDIF
  !
END SUBROUTINE phonon_read



SUBROUTINE eph_long_read
  !
#INCLUDE "f_defs.h"
  !
  USE kinds,     ONLY : DP
  USE phcom,     ONLY : nmodes
  USE epwcom,    ONLY : nqf1, nqf2, nqf3, ephl_read, neptemp, nepdope
  USE io_files,  ONLY : tmp_dir, prefix
  USE elph2,     ONLY : eph_vogl
  USE bte_var,     ONLY : nq_ful
#ifdef __PARA
  USE mp,        ONLY : mp_barrier, mp_bcast
  USE io_global, ONLY : ionode_id, stdout
  USE mp_global, ONLY : my_pool_id, inter_pool_comm
#ENDIF
  !
  IMPLICIT NONE
  !
  INTEGER             :: iq, imode, jmode, ir
  CHARACTER(LEN=256)  :: filename
  LOGICAL             :: file_exist
  !
  !
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     filename = TRIM(tmp_dir)//TRIM(prefix)//'.ephl' 
     !
     ! check exisitence of file
     INQUIRE (FILE=filename,EXIST=file_exist)
     IF (.NOT. file_exist) CALL errore('ephl_read','cannot find .ephl file (ephmat_L)',1)
     !
     OPEN (99901,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=2*nmodes*neptemp*nepdope*DP,STATUS='old')
     DO iq = 1, nq_ful
        READ (99901,REC=iq) eph_vogl(1:nmodes,1:neptemp,1:nepdope,iq)
     ENDDO
     CLOSE (99901)
     !
     IF (ephl_read) WRITE (stdout,'(/13x,a,i3,a,i3,a,i3,a)') 'Load long-range part of e-ph matrix from files (', &
                                                             nqf1, '*', nqf2, '*', nqf3, ' q-mesh)'
     !
  ENDIF
  !
#ifdef __PARA
  CALL mp_bcast (eph_vogl,ionode_id,inter_pool_comm)
  CALL mp_barrier (inter_pool_comm)
#ENDIF
  !
END SUBROUTINE eph_long_read



SUBROUTINE deallocate_shuffle ()
  !
  use epwcom, only: eimp_mode
  USE elph2, ONLY : sigmai_mode_all_abs, sigmai_mode_all_emi, &
                    sigmai_mode_all_ela_intra, sigmai_mode_all_ela_inter, sigmai_mode_all_alloy_inter, &
                    sigmai_mode_all_intra, sigmai_mode_all_inter, sigmai_mode_all_intra, &
                    sigmai_mode_all_alloy_intra, &
                    etf_all, vel_all, wf_irr, vph_irr, wf_all, vph_all, uf_all, &
                    etf_ful, vel_ful, wf_ful, vph_ful, uf_ful, uf_irr, &
                    epf17, etf, etf_ks, wf, dmef, vmef, eph_vogl, eimpmatq, eimpmatwp, eimpmatwe
  USE bte_var
#ifdef __PARA
  USE mp,        ONLY : mp_barrier
  USE io_global, ONLY : ionode_id, stdout
  USE mp_global, ONLY : my_pool_id, inter_pool_comm
#ENDIF
  !
  IMPLICIT NONE
  !
  IF (ALLOCATED (epf17))               DEALLOCATE (epf17)
  IF (ALLOCATED (etf))                 DEALLOCATE (etf)
  IF (ALLOCATED (etf_ks))              DEALLOCATE (etf_ks)  
  IF (ALLOCATED (wf))                  DEALLOCATE (wf)
  IF (ALLOCATED (dmef))                DEALLOCATE (dmef)
  IF (ALLOCATED (vmef))                DEALLOCATE (vmef)
  IF (ALLOCATED (etf_all))             DEALLOCATE (etf_all)
  IF (ALLOCATED (vel_all))             DEALLOCATE (vel_all)
  IF (ALLOCATED (wf_irr))              DEALLOCATE (wf_irr)
  IF (ALLOCATED (vph_irr))             DEALLOCATE (vph_irr)
  IF (ALLOCATED (uf_irr))              DEALLOCATE (uf_irr)
  IF (ALLOCATED (wf_all))              DEALLOCATE (wf_all)
  IF (ALLOCATED (vph_all))             DEALLOCATE (vph_all)
  IF (ALLOCATED (uf_all))              DEALLOCATE (uf_all)
  IF (ALLOCATED (etf_ful))             DEALLOCATE (etf_ful)
  IF (ALLOCATED (vel_ful))             DEALLOCATE (vel_ful)
  IF (ALLOCATED (uf_ful))              DEALLOCATE (uf_ful)
  IF (ALLOCATED (eph_vogl))            DEALLOCATE (eph_vogl)
  !
  IF (ALLOCATED (sigmai_mode_all_abs)) DEALLOCATE (sigmai_mode_all_abs)
  IF (ALLOCATED (sigmai_mode_all_emi)) DEALLOCATE (sigmai_mode_all_emi)
  if (eimp_mode > 0) then
     IF (ALLOCATED (sigmai_mode_all_ela_intra)) DEALLOCATE (sigmai_mode_all_ela_intra)
     IF (ALLOCATED (sigmai_mode_all_ela_inter)) DEALLOCATE (sigmai_mode_all_ela_inter)
     IF (ALLOCATED (sigmai_mode_all_alloy_intra)) DEALLOCATE (sigmai_mode_all_alloy_intra)
     IF (ALLOCATED (sigmai_mode_all_alloy_inter)) DEALLOCATE (sigmai_mode_all_alloy_inter)
  endif
  IF (ALLOCATED (eimpmatq)) DEALLOCATE (eimpmatq)
  IF (ALLOCATED (eimpmatwp)) DEALLOCATE (eimpmatwp)
  IF (ALLOCATED (eimpmatwe)) DEALLOCATE (eimpmatwe)
  IF (ALLOCATED (sigmai_mode_all_intra)) DEALLOCATE (sigmai_mode_all_intra)
  IF (ALLOCATED (sigmai_mode_all_inter)) DEALLOCATE (sigmai_mode_all_inter)
  IF (ALLOCATED (wf_ful))              DEALLOCATE (wf_ful)
  IF (ALLOCATED (vph_ful))             DEALLOCATE (vph_ful)
  IF (ALLOCATED (nscat_all))           DEALLOCATE (nscat_all)
  !
  IF (ALLOCATED (symmat_lat))          DEALLOCATE (symmat_lat)
  !
  !  
  CALL mp_barrier (inter_pool_comm)
  !
END SUBROUTINE deallocate_shuffle


FUNCTION V_coulomb (icharge, a1, a2, a3, alat, rp, dielec)

   USE kinds,         ONLY : DP
   use constants_epw, only : pi, au2m, ryd2ev
   implicit none
   !
   real(kind=DP) :: V_coulomb
   real(kind=DP) :: icharge, q0, e0, a1(3), a2(3), a3(3), alat, rp(3), dielec, &
              omega, omega_, temp(3), b1(3), b2(3), b3(3), igamma, &
              dv_R, ix_, iy_, iz_, rx, ry, rz, rdiff, &
              gx, gy, gz, g2, gdotr, epsi0, dv_C
   complex(kind=DP) :: dv_G, one
   integer :: nR, nG, ix, iy, iz
   !
   if ((rp(1)==0).and.(rp(2)==0).and.(rp(3)==0)) then
      V_coulomb = 0.d0
      return
   endif
   !
   ! constants
   ! electron carries negative charge
   q0 = sqrt(2.d0)
   e0 = -q0
   epsi0 = 1/4.d0/pi
   !
   one = dcmplx(0.d0, 1.d0)
   !
   call cross(a1, a2, temp)
   omega = abs(dot_product(a3, temp))

   call cross(a2, a3, temp)
   omega_ = abs(dot_product(a1, temp))
   b1 = 2*pi*temp/omega_
   call cross(a3, a1, temp)
   omega_ = abs(dot_product(a2, temp))
   b2 = 2*pi*temp/omega_
   call cross(a1, a2, temp)
   omega_ = abs(dot_product(a3, temp))
   b3 = 2*pi*temp/omega_

   ! convergence parameter
   igamma = 2
   nR = 21
   nG = 21
   !
   ! eWald summation over real space
   dv_R = 0.d0
   do ix = 1, nR
      ix_ = ix - (nR+1)/2.d0
      do iy = 1, nR
         iy_ = iy - (nR+1)/2.d0
         do iz = 1, nR
            iz_ = iz - (nR+1)/2.d0

            rx = ix_*a1(1) + iy_*a2(1) + iz_*a3(1) - rp(1)
            ry = ix_*a1(2) + iy_*a2(2) + iz_*a3(2) - rp(2)
            rz = ix_*a1(3) + iy_*a2(3) + iz_*a3(3) - rp(3)
            rdiff = sqrt(rx**2 + ry**2 + rz**2)
            dv_R = dv_R + erfc(igamma*rdiff)/rdiff
         enddo
      enddo
   enddo
   !
   ! eWald summation over reciprocal space
   dv_G = 0.d0
   do ix = 1, nG
      ix_ = ix - (nG+1)/2.d0
      do iy = 1, nG
         iy_ = iy - (nG+1)/2.d0
         do iz = 1, nG
            iz_ = iz - (nG+1)/2.d0

            gx = ix_*b1(1) + iy_*b2(1) + iz_*b3(1)
            gy = ix_*b1(2) + iy_*b2(2) + iz_*b3(2)
            gz = ix_*b1(3) + iy_*b2(3) + iz_*b3(3)
            g2 = gx**2 + gy**2 + gz**2
            gdotr = gx*rp(1) + gy*rp(2) + gz*rp(3)
            if (g2 > 0) &
               dv_G = dv_G + (4*pi/omega)*exp(one*gdotr - g2/igamma/igamma/4.d0)/g2
         enddo
      enddo
   enddo
   !
   ! constant residue
   dv_C = -pi/omega/igamma/igamma
   !
   ! should be in unit of [Ryd]
   V_coulomb = (((icharge*q0)*e0)/(dielec*epsi0))*(dv_R + dv_G + dv_C)/4/pi/alat
   !
   RETURN
   !
END FUNCTION V_coulomb


subroutine cross (a, b, c)
   USE kinds,         ONLY : DP
   implicit none
   real(kind=DP) :: a(3), b(3), c(3)
   !
   c(1) = a(2)*b(3) - a(3)*b(2)
   c(2) = a(3)*b(1) - a(1)*b(3)
   c(3) = a(1)*b(2) - a(2)*b(1)
   !
   return
end subroutine cross


! =========================================================================
SUBROUTINE eimp_interp_tetra (xq, bg_max_red, ng_max_red, dvimp_qf)
! =========================================================================
   !
   USE kinds,         ONLY : DP
   USE pwcom,         ONLY : bg
   USE phcom,         ONLY : nq1, nq2, nq3
   use elph2,         only : dvimp_qc, bg_max, ng_max
   use io_global,     only : stdout
 
   implicit none
 
   integer, intent(in) :: bg_max_red, ng_max_red
   real(kind=DP), intent(in)   :: xq(3)
   complex(kind=DP), intent(out) :: dvimp_qf(ng_max_red)
 
   integer :: nq000(4), nq100(4), nq010(4), nq001(4), &
              nq101(4), nq011(4), nq110(4), nq111(4), &
              iG1, iG1p, iG2, iG2p, iG3, iG3p, &
              iq1, iq1p, iq2, iq2p, iq3, iq3p, &
              igx, igy, igz, ind_g, ind_gx, ind_gy, ind_gz
   real(kind=DP) :: xq1, xq2, xq3, xq000(3), xq100(3), xq010(3), xq001(3), &
                    xq101(3), xq011(3), xq110(3), xq111(3), xq_v1(3), xq_v2(3), &
                    xq_v3(3), xq_v4(3), xqint(3), interpx, interpy, interpz
   complex(kind=DP) :: dvimpqf_int, dvimpq_v1, dvimpq_v2, dvimpq_v3, &
                       dvimpq_v4

   dvimp_qf = 0.d0
   iG1 = 0
   iG2 = 0
   iG3 = 0

   xq1 = xq(1)
   do while (xq1 < 0)
      xq1 = xq1 + 1
      iG1 = iG1 - 1
   enddo
   do while (xq1 >= 1)
      xq1 = xq1 - 1
      iG1 = iG1 + 1
   enddo

   xq2 = xq(2)
   do while (xq2 < 0)
      xq2 = xq2 + 1
      iG2 = iG2 - 1
   enddo
   do while (xq2 >= 1)
      xq2 = xq2 - 1
      iG2 = iG2 + 1
   enddo

   xq3 = xq(3)
   do while (xq3 < 0)
      xq3 = xq3 + 1
      iG3 = iG3 - 1
   enddo
   do while (xq3 >= 1)
      xq3 = xq3 - 1
      iG3 = iG3 + 1
   enddo

   iG1p = iG1
   iq1 = int(xq1*nq1)
   interpx = xq1*nq1 - iq1
   iq1p = iq1 + 1
   if (iq1p == nq1) then
      iq1p = 0
      iG1p = iG1p + 1
   endif

   iG2p = iG2
   iq2 = int(xq2*nq2)
   interpy = xq2*nq2 - iq2
   iq2p = iq2 + 1
   if (iq2p == nq2) then
      iq2p = 0
      iG2p = iG2p + 1
   endif

   iG3p = iG3
   iq3 = int(xq3*nq3)
   interpz = xq3*nq3 - iq3
   iq3p = iq3 + 1
   if (iq3p == nq3) then
      iq3p = 0
      iG3p = iG3p + 1
   endif

   nq000(1:4) = (/iq1*nq2*nq3 + iq2*nq3 + iq3 + 1, iG1, iG2, iG3/)
   nq100(1:4) = (/iq1p*nq2*nq3+ iq2*nq3 + iq3 + 1, iG1p,iG2, iG3/)
   nq010(1:4) = (/iq1*nq2*nq3 + iq2p*nq3+ iq3 + 1, iG1, iG2p,iG3/)
   nq001(1:4) = (/iq1*nq2*nq3 + iq2*nq3 + iq3p+ 1, iG1, iG2, iG3p/)
   nq101(1:4) = (/iq1p*nq2*nq3+ iq2*nq3 + iq3p+ 1, iG1p,iG2, iG3p/)
   nq011(1:4) = (/iq1*nq2*nq3 + iq2p*nq3+ iq3p+ 1, iG1, iG2p,iG3p/)
   nq110(1:4) = (/iq1p*nq2*nq3+ iq2p*nq3+ iq3 + 1, iG1p,iG2p,iG3/)
   nq111(1:4) = (/iq1p*nq2*nq3+ iq2p*nq3+ iq3p+ 1, iG1p,iG2p,iG3p/)

!   write(stdout,*) ' xq =', xq
!   write(stdout,*) ' iG =', iG1,iG2,iG3
!   write(stdout,*) ' iq =', iq1,iq2,iq3
!   write(stdout,*) ' interp =', interpx,interpy,interpz
!   write(stdout,*) ' nq000 =', nq000

   xq000 = 0.d0
   xq100 = bg(1:3,1)
   xq010 = bg(1:3,2)
   xq001 = bg(1:3,3)
   xq110 = bg(1:3,1) + bg(1:3,2)
   xq101 = bg(1:3,1) + bg(1:3,3)
   xq011 = bg(1:3,2) + bg(1:3,3)
   xq111 = bg(1:3,1) + bg(1:3,2) + bg(1:3,3)
   xqint = interpx*bg(1:3,1) + interpy*bg(1:3,2) + interpz*bg(1:3,3)
   CALL cryst_to_cart ( 1, xq000, bg, 1 )
   CALL cryst_to_cart ( 1, xq100, bg, 1 )
   CALL cryst_to_cart ( 1, xq010, bg, 1 )
   CALL cryst_to_cart ( 1, xq001, bg, 1 )
   CALL cryst_to_cart ( 1, xq110, bg, 1 )
   CALL cryst_to_cart ( 1, xq101, bg, 1 )
   CALL cryst_to_cart ( 1, xq011, bg, 1 )
   CALL cryst_to_cart ( 1, xq111, bg, 1 )
   CALL cryst_to_cart ( 1, xqint, bg, 1 )

   ! interpolate the Fourier component of defect potential

   ind_g = 0
   do igz = -bg_max_red, bg_max_red
      ind_gz = igz+bg_max+1
      do igy = -bg_max_red, bg_max_red
         ind_gy = igy+bg_max+1
         do igx = -bg_max_red, bg_max_red
            ind_gx = igx+bg_max+1
            !
            ind_g = ind_g + 1

            ! find the tetrahedra inside the cube that encloses the q-point
 
            if (interpx>=interpy) then
               if (interpz>=interpx) then
                  dvimpq_v1 = dvimp_qc(nq001(1),nq001(2)+ind_gx,nq001(3)+ind_gy,nq001(4)+ind_gz)
                  dvimpq_v2 = dvimp_qc(nq101(1),nq101(2)+ind_gx,nq101(3)+ind_gy,nq101(4)+ind_gz)
                  dvimpq_v3 = dvimp_qc(nq000(1),nq000(2)+ind_gx,nq000(3)+ind_gy,nq000(4)+ind_gz)
                  dvimpq_v4 = dvimp_qc(nq111(1),nq111(2)+ind_gx,nq111(3)+ind_gy,nq111(4)+ind_gz)
                  xq_v1 = xq001
                  xq_v2 = xq101
                  xq_v3 = xq000
                  xq_v4 = xq111
               else
                  if (interpz>=interpy) then
                     dvimpq_v1 = dvimp_qc(nq000(1),nq000(2)+ind_gx,nq000(3)+ind_gy,nq000(4)+ind_gz)
                     dvimpq_v2 = dvimp_qc(nq101(1),nq101(2)+ind_gx,nq101(3)+ind_gy,nq101(4)+ind_gz)
                     dvimpq_v3 = dvimp_qc(nq100(1),nq100(2)+ind_gx,nq100(3)+ind_gy,nq100(4)+ind_gz)
                     dvimpq_v4 = dvimp_qc(nq111(1),nq111(2)+ind_gx,nq111(3)+ind_gy,nq111(4)+ind_gz)
                     xq_v1 = xq000
                     xq_v2 = xq101
                     xq_v3 = xq100
                     xq_v4 = xq111
                  else
                     dvimpq_v1 = dvimp_qc(nq000(1),nq000(2)+ind_gx,nq000(3)+ind_gy,nq000(4)+ind_gz)
                     dvimpq_v2 = dvimp_qc(nq110(1),nq110(2)+ind_gx,nq110(3)+ind_gy,nq110(4)+ind_gz)
                     dvimpq_v3 = dvimp_qc(nq100(1),nq100(2)+ind_gx,nq100(3)+ind_gy,nq100(4)+ind_gz)
                     dvimpq_v4 = dvimp_qc(nq111(1),nq111(2)+ind_gx,nq111(3)+ind_gy,nq111(4)+ind_gz)
                     xq_v1 = xq000
                     xq_v2 = xq110
                     xq_v3 = xq100
                     xq_v4 = xq111
                  endif
               endif
            else
               if (interpz<=interpx) then
                  dvimpq_v1 = dvimp_qc(nq010(1),nq010(2)+ind_gx,nq010(3)+ind_gy,nq010(4)+ind_gz)
                  dvimpq_v2 = dvimp_qc(nq111(1),nq111(2)+ind_gx,nq111(3)+ind_gy,nq111(4)+ind_gz)
                  dvimpq_v3 = dvimp_qc(nq110(1),nq110(2)+ind_gx,nq110(3)+ind_gy,nq110(4)+ind_gz)
                  dvimpq_v4 = dvimp_qc(nq000(1),nq000(2)+ind_gx,nq000(3)+ind_gy,nq000(4)+ind_gz)
                  xq_v1 = xq010
                  xq_v2 = xq111
                  xq_v3 = xq110
                  xq_v4 = xq000
               else
                  if (interpz<=interpy) then
                     dvimpq_v1 = dvimp_qc(nq011(1),nq011(2)+ind_gx,nq011(3)+ind_gy,nq011(4)+ind_gz)
                     dvimpq_v2 = dvimp_qc(nq111(1),nq111(2)+ind_gx,nq111(3)+ind_gy,nq111(4)+ind_gz)
                     dvimpq_v3 = dvimp_qc(nq010(1),nq010(2)+ind_gx,nq010(3)+ind_gy,nq010(4)+ind_gz)
                     dvimpq_v4 = dvimp_qc(nq000(1),nq000(2)+ind_gx,nq000(3)+ind_gy,nq000(4)+ind_gz)
                     xq_v1 = xq011
                     xq_v2 = xq111
                     xq_v3 = xq010
                     xq_v4 = xq000
                  else
                     dvimpq_v1 = dvimp_qc(nq001(1),nq001(2)+ind_gx,nq001(3)+ind_gy,nq001(4)+ind_gz)
                     dvimpq_v2 = dvimp_qc(nq011(1),nq011(2)+ind_gx,nq011(3)+ind_gy,nq011(4)+ind_gz)
                     dvimpq_v3 = dvimp_qc(nq111(1),nq111(2)+ind_gx,nq111(3)+ind_gy,nq111(4)+ind_gz)
                     dvimpq_v4 = dvimp_qc(nq000(1),nq000(2)+ind_gx,nq000(3)+ind_gy,nq000(4)+ind_gz)
                     xq_v1 = xq001
                     xq_v2 = xq011
                     xq_v3 = xq111
                     xq_v4 = xq000
                  endif
               endif
            endif
            !
            ! check
!            if ((igx == 0) .and. (igy == 0) .and. (igz == 0)) then
!               write(stdout,*)
!               write(stdout,*) ' inside eimp_interp, test:'
!               write(stdout,*) ' xq_int', xqint
!               write(stdout,*) ' xq_v1..4', xq_v1, xq_v2, xq_v3, xq_v4
!               write(stdout,*) ' dvimpq_v1..4', dvimpq_v1, dvimpq_v2, dvimpq_v3, dvimpq_v4
!               write(stdout,*) ' dvimp_qc000:', dvimp_qc(nq000(1),nq000(2)+ind_gx,nq000(3)+ind_gy,nq000(4)+ind_gz)
!            endif

            call tetra_interpl(xqint, xq_v1, xq_v2, xq_v3, xq_v4, &
                               dvimpqf_int, dvimpq_v1, dvimpq_v2, dvimpq_v3, dvimpq_v4)
 
            dvimp_qf(ind_g) = dvimpqf_int
            !
         enddo
      enddo
   enddo
   !
END SUBROUTINE eimp_interp_tetra


subroutine tetra_interpl(x0,x1,x2,x3,x4,f0,f1,f2,f3,f4)
  USE kinds,         ONLY : DP
  use io_global,     only : stdout

  implicit none
  integer     :: info, ipiv(4)
  real(DP)    :: x0(3), x1(3), x2(3), x3(3), x4(3)
  complex(DP) :: f0, f1, f2, f3, f4
  real(DP)    :: eta_trans(4,4), eta_list(4), x0_list(4)

  x0_list(1) = 1.d0
  x0_list(2:4) = x0(1:3)

  eta_trans(1:4,1) = (/1.d0, x1(1), x1(2), x1(3)/)
  eta_trans(1:4,2) = (/1.d0, x2(1), x2(2), x2(3)/)
  eta_trans(1:4,3) = (/1.d0, x3(1), x3(2), x3(3)/)
  eta_trans(1:4,4) = (/1.d0, x4(1), x4(2), x4(3)/)

  CALL DGESV(4, 1, eta_trans, 4, ipiv, x0_list, 4, info)

  eta_list = x0_list
  f0 = eta_list(1)*f1 + eta_list(2)*f2 + eta_list(3)*f3 + eta_list(4)*f4
  !
end subroutine tetra_interpl

