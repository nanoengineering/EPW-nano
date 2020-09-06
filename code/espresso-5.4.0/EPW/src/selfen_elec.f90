!-----------------------------------------------------------------------
SUBROUTINE selfen_elec (ik_red, ik_star, xkk, nscat)
!-----------------------------------------------------------------------
  !
  !  Compute the imaginary part of the electron self energy due to electron-
  !  phonon interaction in the Migdal approximation. This corresponds to 
  !  the electron linewidth (half width). The phonon frequency is taken into
  !  account in the energy selection rule.
  !
  !  Use matrix elements, electronic eigenvalues and phonon frequencies
  !  from ep-wannier interpolation
  !
  !-----------------------------------------------------------------------
#INCLUDE "f_defs.h"
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE ions_base,     ONLY : nat
  USE io_files,      ONLY : tmp_dir, prefix
  USE phcom,         ONLY : nmodes
  USE epwcom,        ONLY : nbndsub, fsthick, eptemp, degaussw, &
                            bte, efermi_read, fermi_energy, smearing, &
                            nptype, epdope, nepdope, save_m_mat, neptemp, nepdope, &
                            eimp_mode, L_D, dielec, lpolar, imp_charge, defectdens, &
                            eimp_sr, screen_polar, frac_type, alloy_pot
  USE pwcom,         ONLY : ef, at, bg, alat, omega
  USE elph2,         ONLY : ibndmin, ibndmax, nbnd_red, xqf, cbnd_emin, vbnd_emax, cfsthick, vfsthick, &
                            etf, epf17, wf_ful, wqf, sigmai_mode_all_abs, sigmai_mode_all_emi, wmax, &
                            sigmai_mode_all_inter, sigmai_mode_all_intra, &
                            epsi, epsil_dyn, ef_epw, sigmai_mode_all_ela_intra, vkq, &
                            nfftmesh, rp, dvimp, eimpf_lr, eimpf_sr, sigmai_mode_all_ela_inter, &
                            eimpf_full, sigmai_mode_all_alloy_inter, sigmai_mode_all_alloy_intra, alelpp
  USE constants_epw, ONLY : pi, twopi, ryd2ev, au2cm, bohr2ang, ci, ryd2thz
  USE tetrahedron
  USE bte_var         
#ifdef __PARA
  USE mp,            ONLY : mp_sum
  USE mp_global,     ONLY : my_pool_id, npool, inter_pool_comm
#ENDIF
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)        :: ik_red, ik_star
  REAL(KIND=DP), INTENT(IN)  :: xkk(3)
  INTEGER, INTENT(OUT)       :: nscat
  !
  REAL(KIND=DP), PARAMETER   :: eps = 1.0d-6
  ! local variable
  INTEGER                    :: ikk, ikq, iq, ibnd0, jbnd0, ibnd, jbnd, imode, itemp, idope, &
                                nscat_q, iq_exp(nq_ful), imesh, iq_star, iq_stop
  REAL(KIND=DP)              :: g2, xq(3), ekk, ekq, wq, ef0, wgq, wgkq, inv_temp, &
                                xxq(3), xkq(3), xkq_cart(3), xq_cart(3), xxq_m(3), &
                                xkl, xkql, xql, xkdotxkq, n_imp, &
                                g2_eimp, dielec0, n_dop, vkdotvkq, vkl, vkql, &
                                cosa(nbndsub,nbndsub,nq_ful), g2_alel
  REAL(KIND=DP)              :: sigmai_abs(nq_ful), sigmai_emi(nq_ful), sigmai_exp(nq_ful), &
                                sigmai_intra(nq_ful), sigmai_inter(nq_ful), &
                                weight_abs(nq_ful), weight_emi(nq_ful), weight_abs_exp(nq_ful), weight_emi_exp(nq_ful), &
                                sigmai_ela_intra(nq_ful), sigmai_ela_inter(nq_ful), &
                                weight_ela(nq_ful), weight_alel(nq_ful), sigmai_alloy_intra(nq_ful), &
                                sigmai_alloy_inter(nq_ful)
  COMPLEX(KIND=DP)           :: epf, epf_tmp(nbnd_red,nbnd_red,nmodes)
  LOGICAL                    :: within_range
  ! external fuction
  REAL(KIND=DP), EXTERNAL    :: wgauss, w0gauss
  ! tetrahedron
  REAL(KIND=DP)              :: ekq_tetra(nq_ful)
  REAL(KIND=DP), ALLOCATABLE :: tetra_c_abs(:,:,:,:), tetra_c_emi(:,:,:,:), &
                                tetra_c_ela(:,:)
  INTEGER                    :: iktetra, itetra_(2,30,nq_ful)
  ! 
  REAL(KIND=DP)              :: gauss_1, gauss_2
  ! sigma file
  CHARACTER(LEN=256)         :: file_ufmt
  CHARACTER(LEN=12)          :: tnph
  CHARACTER(LEN=10)          :: k_num, ibnd_num, jbnd_num, imode_num, itemp_num, idope_num
  ! epf17 file
  CHARACTER(LEN=256)         :: epf17_ufmt
  CHARACTER(LEN=3)           :: epf17_cpu
  !
  integer                    :: ik_v, ikq_v
  real(kind=DP)              :: xk_v(3), xkq_v(3)
  integer                    :: valley_map(2,nq_ful)
  integer, external          :: kq_label
  complex(kind=DP)           :: eimpf_tmp, alelpp_tmp
  complex(kind=DP), allocatable :: eimpf_lr_valley(:,:,:), eimpf_sr_valley(:,:,:), alelpp_valley(:,:,:)
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
  !
  ! initialization of tetrahedron
  IF (smearing .EQ. 'tetra') THEN
     itetra_ = itetra
     if (screen_polar) then
        ALLOCATE (tetra_c_abs(4,ntetra,neptemp,nepdope))
        ALLOCATE (tetra_c_emi(4,ntetra,neptemp,nepdope))
     else
        ALLOCATE (tetra_c_abs(4,ntetra,1,1))
        ALLOCATE (tetra_c_emi(4,ntetra,1,1))
     endif
     if (eimp_mode > 0) then
        ALLOCATE (tetra_c_ela(4,ntetra))
     endif
  ENDIF
  !
  ! initialization of e-imp scattering - angle dependence
  !
  if ((eimp_mode == 1) .or. (eimp_mode == 3) .or. (eimp_mode == 5) .or. (eimp_mode == 7)) then
     !
     cosa(:,:,:) = 0.d0
     !
  elseif ((eimp_mode == 2) .or. (eimp_mode == 4) .or. (eimp_mode == 6) .or. (eimp_mode == 8)) then
     !
     DO iq = 1, nq_ful
        !
        ikq = 2 * iq

        ! consider the fact that elastic scattering is strong for near-angle
        ! case
        do ibnd = 1, nbndsub
           do jbnd = 1, nbndsub
              vkl  = sqrt(vkq(1,ibnd,1)**2 + vkq(2,ibnd,1)**2 + vkq(3,ibnd,1)**2)
              vkql = sqrt(vkq(1,jbnd,ikq)**2 + vkq(2,jbnd,ikq)**2 + vkq(3,jbnd,ikq)**2)
              vkdotvkq = vkq(1,ibnd,1)*vkq(1,jbnd,ikq) + &
                         vkq(2,ibnd,1)*vkq(2,jbnd,ikq) + &
                         vkq(3,ibnd,1)*vkq(3,jbnd,ikq)

              if ((vkl > 1.d-20) .and. (vkql > 1.d-20)) then
                 cosa(ibnd,jbnd,iq) = vkdotvkq / vkl / vkql
              else
                 cosa(ibnd,jbnd,iq) = 0.d0
              endif
           enddo
        enddo
     enddo
     !
  else
     !
     cosa(:,:,:) = 0.d0
     !
  endif
  !
  ! long-range impurity scattering coupling matrix, only do for the first time
  !
  if ((eimp_mode > 0) .and. (eimp_mode < 7) .and. (ik_red == ik_star)) then
     !
     if (.not. lpolar) then
        dielec0 = dielec
     else
        dielec0 = (epsi(1,1) + epsi(2,2) + epsi(3,3)) / 3.d0
     endif
     !
     DO iq = 1, nq_ful
        !
        xxq(:) = xqf_ful(:,iq)
        !
        ! !IMPORTANT NOTE! this brings xq to first-BZ, but it is only centered around
        ! Gamma, not yet fully symmetric
        xq_cart(1) = xxq(1) - NINT(xxq(1))
        xq_cart(2) = xxq(2) - NINT(xxq(2))
        xq_cart(3) = xxq(3) - NINT(xxq(3))
        !
        CALL cryst_to_cart ( 1, xq_cart, bg, 1 )   ! xq_cart in [2*pi/alat],
                                                   ! [alat] being primitive unit cell lattice constant
        xql  = sqrt(xq_cart(1)**2 + xq_cart(2)**2 + xq_cart(3)**2) * &
               (2.d0*pi/alat)
        !
        do idope = 1, nepdope
           do itemp = 1, neptemp
              ! note that the factor of 2*pi/hbar (hbar=1 in AU) is added later,
              ! specifically, pi is added below, while 2 is considered in later
              ! transport calculations
              eimpf_lr(iq,itemp,idope) = - imp_charge*(2.d0*(4.d0*pi)/omega/dielec0) / &
                          (xql**2.d0 + 1.d0/(L_D(itemp,idope)**2.d0))
           enddo
        enddo
     enddo
     !
  endif
  !
  ! differentiate between intra- and inter- processes,
  ! currently only applied to impurity scatterings
  !
!ERROR: currently only works for conduction band
  ! for eimp_mode = 3 or 4, and eimp_mode = 5 or 6 with eimp_sr = false

  call valley_degen(xkk, xk_v, ik_v)

  write(stdout,'(a,1x,3(f15.7))') ' test valley (xkk): ', xkk
  write(stdout,'(a,1x,3(f15.7),1x,i5)') '      near ',xk_v, ik_v
  !
  do iq = 1, nq_ful
     !
     xxq(:) = xqf_ful(:,iq) ! xqf_ful is assumed to be in crys coord
     xkq = xkk + xxq
     !
     call valley_degen(xkq, xkq_v, ikq_v)
     !
     ! modified momentum transfer according to the nearest valley
     xxq_m = (xkq - xkq_v) - (xkk - xk_v)
     !
     if (ikq_v == ik_v) then
        valley_map(1,iq) = 0   ! intra-valley
     else
        valley_map(1,iq) = 1   ! inter-valley
     endif
     valley_map(2,iq) = kq_label(xxq_m)
   
!     if (my_pool_id == 3 ) then
!        write(*,*) ' find negative valley_map at pool #', my_pool_id
!        write(*,*) iq
!        write(*,*) xkk, xkq
!        write(*,*) xk_v, xkq_v
!        write(*,*) xxq, xxq_m
!     endif
     !
!     if (ik_red == ik_star) then
!        ikq = kq_label(xkq)
!        write(stdout,*)
!        write(stdout,*) ' qmesh map and valley test: '
!        write(stdout,'(i5,a,1x,3(f15.7))') iq, ' - xkq: ', xkq
!        write(stdout,'(i5,a,1x,3(f15.7))') ikq, ' - mapped: ', xqf_ful(:,ikq)
!        write(stdout,'(a,1x,3(f15.7),a,i5)') '    near ', xkq_v, ' at #', ikq_v
!     endif
     !
  enddo
  !
  ! add short-range coupling matrix, if needed, and remove long range for
  ! intervaley processes
  ! depending on whether it is intra- or inter- process
  !
  ! ==========================================================================
  ! When there are additional defects other than dopants (when [defectdens] /= 0),
  ! eimpf_lr_valley and eimpf_sr_valley are meant to use n_dop and n_imp density, respectively
  ! ==========================================================================

  if (eimp_mode > 0) then
     !
     if (eimp_mode < 7) then
        allocate ( eimpf_lr_valley (nq_ful,neptemp,nepdope) )
        allocate ( eimpf_sr_valley (nq_ful,nbnd_red,nbnd_red) )
        allocate ( alelpp_valley (nq_ful,nbnd_red,nbnd_red) )
        eimpf_lr_valley = 0.d0
        eimpf_sr_valley = 0.d0
        eimpf_sr_valley = 0.d0
     endif
     !
     if (eimp_mode == 7 .or. eimp_mode == 8) then
        !
        ! we do not allocate valley arrays but directly use eimpf_full
        !
     elseif (eimp_mode == 5 .or. eimp_mode == 6) then
        !
! NOTE: as short range we now use Wannier transform, so adding long range has
! to be careful
        if (.not. eimp_sr) then
           do iq = 1, nq_ful
              if (valley_map(1,iq) == 0) then
                 ! only assign intravalley, intervalley is zero
                 eimpf_lr_valley(iq,:,:) = eimpf_lr(iq,:,:)
              endif
           enddo
        endif
        eimpf_sr_valley = eimpf_sr
        IF (alloy_pot) alelpp_valley = alelpp
        !
     elseif (eimp_mode == 3 .or. eimp_mode == 4) then
        !
        do iq = 1, nq_ful
           if (valley_map(1,iq) == 0) then
              ! intravalley
              if (.not. eimp_sr) eimpf_lr_valley(iq,:,:) = eimpf_lr(iq,:,:)
              eimpf_sr_valley(iq,:,:) = eimpf_sr(iq,:,:)
           elseif (valley_map(1,iq) == 1) then
              ! intervalley
              eimpf_sr_valley(iq,:,:) = eimpf_sr(valley_map(2,iq),:,:)
           endif
        enddo
        !
     else
        !
        do iq = 1, nq_ful
           if (valley_map(1,iq) == 0) then
              ! intravalley
              eimpf_lr_valley(iq,:,:) = eimpf_lr(iq,:,:)
           endif
        enddo
        !
     endif
     !
  endif
  !
  !
  DO ibnd0 = 1, nbnd_red
     !
     DO jbnd0 = 1, nbnd_red
        !
        ibnd = ibnd0+ibndmin-1
        jbnd = jbnd0+ibndmin-1
        !
        if (smearing .EQ. 'tetra') then
           if (eimp_mode > 0) then
              ! linewidth due to impurity scattering
              DO iq = 1, nq_ful
                 ekq_tetra(iq) = etf(jbnd,2*iq)
              ENDDO
              CALL eigen_tet(ntetra,ekq_tetra,tetra_i,tetra_w,nq_ful)
              CALL weight_tet(nq_ful,ntetra,etf(ibnd,1),tetra_i,tetra_w,tetra_c_ela,wkt)
           endif
        endif
        !
        DO imode = 1, nmodes
           !
           IF (smearing .EQ. 'tetra') THEN   
              !
              if (screen_polar) then
                 do itemp = 1, neptemp
                    do idope = 1, nepdope
                       ! tetrahedron
                       ! linewidth due to absorption of phonon
                       DO iq = 1, nq_ful
                          ekq_tetra(iq) = -wf_ful(imode,itemp,idope,iq) + etf(jbnd,2*iq)
                       ENDDO
                       CALL eigen_tet(ntetra,ekq_tetra,tetra_i,tetra_w,nq_ful)
                       CALL weight_tet(nq_ful,ntetra,etf(ibnd,1),tetra_i,tetra_w,tetra_c_abs(:,:,itemp,idope),wkt)
                       !
                       ! linewidth due to emission of phonon
                       DO iq = 1, nq_ful
                          ekq_tetra(iq) = wf_ful(imode,itemp,idope,iq) + etf(jbnd,2*iq)
                       ENDDO
                       CALL eigen_tet(ntetra,ekq_tetra,tetra_i,tetra_w,nq_ful)
                       CALL weight_tet(nq_ful,ntetra,etf(ibnd,1),tetra_i,tetra_w,tetra_c_emi(:,:,itemp,idope),wkt)
                       !
                    enddo
                 enddo
              else
                 ! tetrahedron
                 ! linewidth due to absorption of phonon
                 DO iq = 1, nq_ful
                    ekq_tetra(iq) = -wf_ful(imode,1,1,iq) + etf(jbnd,2*iq)
                 ENDDO
                 CALL eigen_tet(ntetra,ekq_tetra,tetra_i,tetra_w,nq_ful)
                 CALL weight_tet(nq_ful,ntetra,etf(ibnd,1),tetra_i,tetra_w,tetra_c_abs(:,:,1,1),wkt)
                 !
                 ! linewidth due to emission of phonon
                 DO iq = 1, nq_ful
                    ekq_tetra(iq) = wf_ful(imode,1,1,iq) + etf(jbnd,2*iq)
                 ENDDO
                 CALL eigen_tet(ntetra,ekq_tetra,tetra_i,tetra_w,nq_ful)
                 CALL weight_tet(nq_ful,ntetra,etf(ibnd,1),tetra_i,tetra_w,tetra_c_emi(:,:,1,1),wkt)
                 !
              endif
           ENDIF
           !
           DO itemp = 1, neptemp
              !
              inv_temp = 1.0d0/eptemp(itemp)
              !
              DO idope = 1, nepdope
                 !
                 n_dop = abs(epdope(idope)) * ((bohr2ang*1.0d-8)**3.d0) * omega
                 n_imp = abs(defectdens) * ((bohr2ang*1.0d-8)**3.d0) * omega
                 !
                 nscat = 0 
                 !
                 IF (efermi_read) THEN
                    ef0 = fermi_energy
                 ELSE
                    ef0 = ef_epw(itemp,idope)
                 ENDIF
                 !
                 weight_abs = 0.0d0
                 weight_emi = 0.0d0
                 sigmai_abs = 0.0d0
                 sigmai_emi = 0.0d0
                 sigmai_inter = 0.0d0
                 sigmai_intra = 0.0d0
                 if (eimp_mode > 0) then
                    weight_ela = 0.0d0
                    sigmai_ela_intra = 0.0d0
                    sigmai_ela_inter = 0.0d0
                    IF (alloy_pot) THEN
                       weight_alel= 0.0d0
                       sigmai_alloy_inter = 0.0d0
                       sigmai_alloy_intra = 0.0d0
                    ENDIF
                 endif
                 !
                 DO iq = 1, nq_ful
                    !
                    ikq = 2 * iq
                    ikk = ikq - 1
                    !
                    xq(:) = xqf_ful(:,iq)
                    !
                    within_range = .FALSE.
                    !
                    IF ( (MAXVAL(ABS(xq(:))) .GT. eps) .AND. &
                         (etf(ibnd,ikk) .GE. vbnd_emax-vfsthick .AND. etf(ibnd,ikk) .LE. cbnd_emin+cfsthick) .AND. &
                         (etf(jbnd,ikq) .GE. vbnd_emax-vfsthick .AND. etf(jbnd,ikq) .LE. cbnd_emin+cfsthick) ) within_range = .TRUE.
                    !
                    IF (within_range) THEN
                       !     
                       ! the energy of the electron at k (relative to Ef)
                       ekk = etf(ibnd,ikk) - ef0
                       ! the phonon frequency and Bose occupation
                       if (screen_polar) then
                          wq = wf_ful(imode,itemp,idope,iq)
                       else
                          wq = wf_ful(imode,1,1,iq)
                       endif
                       wgq = wgauss(-wq*inv_temp,-99)
                       wgq = wgq/(1.0d0-2.0d0*wgq)
                       ! the fermi occupation for k+q
                       ekq = etf(jbnd,ikq) - ef0
                       wgkq = wgauss(-ekq*inv_temp,-99)  
                       !
                       ! NOTE: epf17 takes into account the screening effect
                       ! when [screen_polar] is true
                       !
                       IF (save_m_mat) THEN
                          READ (17517,REC=iq) epf_tmp(1:nbnd_red,1:nbnd_red,1:nmodes)
                          epf = epf_tmp(jbnd0,ibnd0,imode)
                       ELSE
                          epf = epf17(iq,jbnd0,ibnd0,imode,itemp,idope)
                       ENDIF
                       !
                       if (wq /= 0.d0) then
                          g2 = (ABS(epf)**2.0d0)/(2.0d0*wq)
                       else
                          g2 = 0.d0
                       endif

                       if (wq /= 0.d0) then
                          if (eimp_mode > 0) then
                             if (defectdens /= 0) then
                                g2_eimp = n_dop * (abs(eimpf_lr_valley(iq,itemp,idope))**2.d0) + &
                                          n_imp * (abs(eimpf_sr_valley(iq,jbnd0,ibnd0))**2.d0)
                             else
                                if (alloy_pot) then
                                   alelpp_tmp = alelpp_valley(iq,jbnd0,ibnd0)  ! for 2x2x2 supercel binary Si-fcc lattice only
                                   g2_alel = ABS(alelpp_tmp)**2
                                   eimpf_tmp = eimpf_lr_valley(iq,itemp,idope)
                                   g2_eimp = n_dop * (abs(eimpf_tmp)**2.d0)
                                elseif (eimp_mode < 7) then
                                   eimpf_tmp = eimpf_lr_valley(iq,itemp,idope) + &
                                               eimpf_sr_valley(iq,jbnd0,ibnd0)
                                   g2_eimp = n_dop * (abs(eimpf_tmp)**2.d0)
                                elseif (eimp_mode == 7 .or. eimp_mode == 8) then
                                  
                                   eimpf_tmp = eimpf_full(jbnd0,ibnd0,itemp,idope,iq)
                                   g2_eimp = n_dop * (abs(eimpf_tmp)**2.d0)
                                endif
                             endif
                          else
                             g2_eimp = 0.d0
                          endif
                       else
                          g2_eimp = 0.d0
                       endif
                       !
                       ! tetrahedral smearing
                       IF (smearing .EQ. 'tetra') THEN
                          !
                          iktetra = 1
                          DO WHILE (itetra_(1,iktetra,iq) .NE. 0)
                             !
                             if (screen_polar) then
                                weight_abs(iq) = weight_abs(iq) + pi * g2 * tetra_c_abs(itetra_(2,iktetra,iq),itetra_(1,iktetra,iq),itemp,idope)
                                weight_emi(iq) = weight_emi(iq) + pi * g2 * tetra_c_emi(itetra_(2,iktetra,iq),itetra_(1,iktetra,iq),itemp,idope)
                             else
                                weight_abs(iq) = weight_abs(iq) + pi * g2 * tetra_c_abs(itetra_(2,iktetra,iq),itetra_(1,iktetra,iq),1,1)
                                weight_emi(iq) = weight_emi(iq) + pi * g2 * tetra_c_emi(itetra_(2,iktetra,iq),itetra_(1,iktetra,iq),1,1)
                             endif

                             if ((eimp_mode > 0) .and. (imode == 1)) then
                                weight_ela(iq) = weight_ela(iq) + pi * g2_eimp * tetra_c_ela(itetra_(2,iktetra,iq),itetra_(1,iktetra,iq))
                             endif
                             !
                             if (alloy_pot .and. (imode == 1)) then
                                weight_alel(iq) = weight_alel(iq) + pi * g2_alel * tetra_c_ela(itetra_(2,iktetra,iq),itetra_(1,iktetra,iq))
                             endif
                             !
                             iktetra = iktetra + 1
                             IF (iktetra .GT. 30) CALL errore('selfen_elec','too many tetrahedron associated with one k-point',1)
                             !
                          ENDDO  
                          !
                          sigmai_abs(iq) = (wgkq+wgq) * weight_abs(iq)
                          sigmai_emi(iq) = (1.0d0-wgkq+wgq) * weight_emi(iq)
                          if ((eimp_mode > 0)) then

                             if (valley_map(1,iq) == 0) then
                                ! NOTE: no (1 - f_k+q) factor
                                if (imode == 1)  then
                                   sigmai_ela_intra(iq) = (1.d0 - cosa(ibnd,jbnd,iq)) * weight_ela(iq)
                                end if
                                sigmai_intra(iq) = sigmai_abs(iq) + sigmai_emi(iq)
                                IF ((alloy_pot) .and. (imode == 1)) sigmai_alloy_intra(iq) = 64*weight_alel(iq)*frac_type(1)*frac_type(2)
                             elseif (valley_map(1,iq) == 1) then
                                if (imode == 1) sigmai_ela_inter(iq) = (1.d0 - cosa(ibnd,jbnd,iq)) * weight_ela(iq)
                                sigmai_inter(iq) = sigmai_abs(iq) + sigmai_emi(iq)
                                IF ((alloy_pot) .and. (imode == 1))  sigmai_alloy_inter(iq) = 64*weight_alel(iq)*frac_type(1)*frac_type(2)
                             endif
                          endif
                          !
                       ! lorentzian smearing
                       ELSEIF (smearing .EQ. 'lortz') THEN
                          !
                          weight_abs(iq) = pi * g2 * (1.0d0/nq_ful) * w0gauss((ekk-ekq+wq)/degaussw,0) / degaussw
                          weight_emi(iq) = pi * g2 * (1.0d0/nq_ful) * w0gauss((ekk-ekq-wq)/degaussw,0) / degaussw 
                          sigmai_abs(iq) = (wgkq+wgq) * weight_abs(iq)
                          sigmai_emi(iq) = (1.0d0-wgkq+wgq) * weight_emi(iq)

       if ((eimp_mode > 0) .and. (imode == 1))  weight_ela(iq) = weight_ela(iq) + pi * g2_eimp * (1.0d0/nq_ful) * w0gauss((ekk-ekq)/degaussw,0) / degaussw

       if (alloy_pot .and. (imode == 1))  weight_alel(iq) = weight_alel(iq) + pi * g2_alel * (1.0d0/nq_ful) * w0gauss((ekk-ekq)/degaussw,0) / degaussw

                          if ((eimp_mode > 0)) then
                             if (valley_map(1,iq) == 0) then
                                ! NOTE: no (1 - f_k+q) factor
            if (imode == 1)  sigmai_ela_intra(iq) = (1.d0 - cosa(ibnd,jbnd,iq)) * weight_ela(iq)
            sigmai_intra(iq) = sigmai_abs(iq) + sigmai_emi(iq)
            IF ((alloy_pot) .and. (imode == 1)) sigmai_alloy_intra(iq) = 64*weight_alel(iq)*frac_type(1)*frac_type(2)

                             elseif (valley_map(1,iq) == 1) then
            if (imode == 1) sigmai_ela_inter(iq) = (1.d0 - cosa(ibnd,jbnd,iq)) * weight_ela(iq)
            sigmai_inter(iq) = sigmai_abs(iq) + sigmai_emi(iq)
            IF ((alloy_pot) .and. (imode == 1))  sigmai_alloy_inter(iq) = 64*weight_alel(iq)*frac_type(1)*frac_type(2)

                             endif
                          endif
                          !
                       ELSE
                          !
                          CALL errore('selfen_elec','wrong smearing flag',1)
                          !
                       ENDIF
                       !
                       IF (bte .EQ. 18) THEN
                          !
                          ! only check electron-phonon coupling scattering rate
                          ! here
                          IF (ibnd .EQ. band_ch(ik_red)) sigmai_ch(itemp,idope,imode,iq,ik_red) = sigmai_ch(itemp,idope,imode,iq,ik_red) + &
                                                                                                  sigmai_abs(iq) + sigmai_emi(iq)
                          !
                       ENDIF
                       !
                       nscat = nscat + 1
                       !
                    ENDIF ! within_range
                    !
                 ENDDO ! iq
                 !
                 !
                 sigmai_mode_all_abs(itemp,idope,imode,ibnd0,ik_red) = sigmai_mode_all_abs(itemp,idope,imode,ibnd0,ik_red) + SUM(sigmai_abs(:))
                 sigmai_mode_all_emi(itemp,idope,imode,ibnd0,ik_red) = sigmai_mode_all_emi(itemp,idope,imode,ibnd0,ik_red) + SUM(sigmai_emi(:))

                sigmai_mode_all_inter(itemp,idope,imode,ibnd0,ik_red) = sigmai_mode_all_inter(itemp,idope,imode,ibnd0,ik_red) + SUM(sigmai_inter(:))
                 sigmai_mode_all_intra(itemp,idope,imode,ibnd0,ik_red) = sigmai_mode_all_intra(itemp,idope,imode,ibnd0,ik_red) + SUM(sigmai_intra(:))

!ERROR, inter-valley scattering from short-range potential seems to be still
!much weaker than intra-valley process

                 if ((eimp_mode > 0) .and. (imode == 1)) then
                    ! NOTE: differentia intra from inter
                    sigmai_mode_all_ela_intra(itemp,idope,ibnd0,ik_red)       = sigmai_mode_all_ela_intra(itemp,idope,ibnd0,ik_red)       + SUM(sigmai_ela_intra(:))
                    sigmai_mode_all_ela_inter(itemp,idope,ibnd0,ik_red)       = sigmai_mode_all_ela_inter(itemp,idope,ibnd0,ik_red)       + SUM(sigmai_ela_inter(:))
                 !  al-el scattering rate
                 IF (alloy_pot)  THEN
 sigmai_mode_all_alloy_inter(itemp,idope,ibnd0,ik_red)  = sigmai_mode_all_alloy_inter(itemp,idope,ibnd0,ik_red)   + SUM(sigmai_alloy_inter(:))
 sigmai_mode_all_alloy_intra(itemp,idope,ibnd0,ik_red)  = sigmai_mode_all_alloy_intra(itemp,idope,ibnd0,ik_red)   + SUM(sigmai_alloy_intra(:))
                 ENDIF
                 endif
                 !
                 !
                 ! outpur binary file
                 IF (bte .EQ. 1) THEN
                    !
                    nscat_q = 0
                    DO iq = 1, nq_ful
                       IF (sigmai_abs(iq)+sigmai_emi(iq) .NE. 0.0d0) THEN
                          nscat_q = nscat_q + 1
                          iq_exp(nscat_q) = iq
                          sigmai_exp(nscat_q) = sigmai_abs(iq)+sigmai_emi(iq)
                       ENDIF
                    ENDDO
                    !
                    IF (nscat_q .EQ. 0) THEN
                       nscat_q = 1
                       iq_exp(nscat_q) = 1
                       sigmai_exp(nscat_q) = 0.0d0
                    ENDIF
                    !
                    WRITE(k_num,'(i10)') ik_red
                    WRITE(ibnd_num,'(i10)') ibnd0
                    WRITE(jbnd_num,'(i10)') jbnd0
                    WRITE(imode_num,'(i10)') imode
                    WRITE(itemp_num,'(i10)') itemp
                    WRITE(idope_num,'(i10)') idope
                    !
                    tnph = 'T'//TRIM(ADJUSTL(itemp_num))//'_N'//TRIM(ADJUSTL(idope_num))//'_ph'//TRIM(ADJUSTL(imode_num))
                    !
                    file_ufmt = 'BTE/SIGMAI/'//TRIM(ADJUSTL(tnph))//'/sigmai_'//TRIM(ADJUSTL(k_num))//'_'//TRIM(ADJUSTL(ibnd_num)) &
                                                                                                    //'_'//TRIM(ADJUSTL(jbnd_num))
                    OPEN (80700,FILE=file_ufmt,FORM='unformatted',ACCESS='direct',RECL=4+nscat_q*(4+DP),STATUS='replace')
                    WRITE (80700,REC=1) nscat_q, iq_exp(1:nscat_q), sigmai_exp(1:nscat_q)
                    CLOSE (80700)
                    !
                 ELSEIF (bte .EQ. 3) THEN
                    !
                    nscat_q = 0
                    DO iq = 1, nq_ful
                       IF (weight_abs(iq)+weight_emi(iq) .NE. 0.0d0) THEN
                          nscat_q = nscat_q + 1
                          iq_exp(nscat_q) = iq
                          weight_abs_exp(nscat_q) = weight_abs(iq)
                          weight_emi_exp(nscat_q) = weight_emi(iq)
                       ENDIF
                    ENDDO
                    !
                    IF (nscat_q .EQ. 0) THEN
                       nscat_q = 1
                       iq_exp(nscat_q) = 1
                       weight_abs_exp(nscat_q) = 0.0d0
                       weight_emi_exp(nscat_q) = 0.0d0
                    ENDIF
                    !
                    WRITE(k_num,'(i10)') ik_red
                    WRITE(ibnd_num,'(i10)') ibnd0
                    WRITE(jbnd_num,'(i10)') jbnd0
                    WRITE(imode_num,'(i10)') imode
                    !
                    tnph = 'T1_N1_ph'//TRIM(ADJUSTL(imode_num))
                    !
                    file_ufmt = 'BTE/WEIGHT/'//TRIM(ADJUSTL(tnph))//'/weight_'//TRIM(ADJUSTL(k_num))//'_'//TRIM(ADJUSTL(ibnd_num)) &
                                                                                                    //'_'//TRIM(ADJUSTL(jbnd_num))
                    OPEN (80700,FILE=file_ufmt,FORM='unformatted',ACCESS='direct',RECL=4+nscat_q*(4+2*DP),STATUS='replace')
                    WRITE (80700,REC=1) nscat_q, iq_exp(1:nscat_q), weight_abs_exp(1:nscat_q), weight_emi_exp(1:nscat_q)
                    CLOSE (80700)
                    !
                 ENDIF
                 !
              ENDDO ! idope
              !
           ENDDO ! itemp
           !
        ENDDO ! imode
        !
     ENDDO ! jbnd0
     !
  ENDDO ! ibnd0
  !
  ! check
  if (.false. .and. (xkk(1) == 0) .and. (xkk(2) == 0.25) .and. (xkk(3) == 0.25)) then

     write(*,*)
     write(*,*) ' -- check inside selfen_elec_k, @ pool#', my_pool_id
     write(*,*) ' band energy =', etf(ibndmin,1)*ryd2ev
     write(*,*) ' sigmai_mode_all_abs =', (0.5d0/sigmai_mode_all_abs(1,13,imode,1,ik_red), imode=1,9)
     write(*,*) ' sigmai_mode_all_emi =', (0.5d0/sigmai_mode_all_emi(1,13,imode,1,ik_red), imode=1,9)
     if (eimp_mode > 0) &
        write(*,*) ' sigmai_mode_all_ela_intra =', 0.5d0/sigmai_mode_all_ela_intra(1,13,1,ik_red)
     if (alloy_pot) &
        write(*,*) ' sigmai_mode_all_alloy_intra =', 0.5d0/sigmai_mode_all_alloy_intra(1,13,1,ik_red)   
  endif

  !
  IF (save_m_mat) CLOSE (17517)
  !
  IF (ALLOCATED(tetra_c_abs)) DEALLOCATE (tetra_c_abs)
  IF (ALLOCATED(tetra_c_emi)) DEALLOCATE (tetra_c_emi)
  if (eimp_mode > 0) then
     IF (ALLOCATED(tetra_c_ela)) DEALLOCATE (tetra_c_ela)
  endif
  !
  RETURN
  !
END SUBROUTINE selfen_elec


! -----------------------------------------------------------------------------
function kq_label(xq)
! -----------------------------------------------------------------------------
!NOTE: works only for nkf = nqf
   USE kinds,         ONLY : DP
   use io_global,     only : stdout
   use epwcom,        only : nqf1, nqf2, nqf3
   !
   implicit none
   !
   real(kind=DP) :: xq(3)
   integer  :: kq_label
   
   integer  :: i, j, k

   ! WRITE (8888,REC=nk) DBLE(i-1)/DBLE(nkf1), DBLE(j-1)/DBLE(nkf2), DBLE(k-1)/DBLE(nkf3)
   i = mod(nint(xq(1)*nqf1), nqf1)
   j = mod(nint(xq(2)*nqf2), nqf2)
   k = mod(nint(xq(3)*nqf3), nqf3)
   if (i < 0) i = i + nqf1
   if (j < 0) j = j + nqf2
   if (k < 0) k = k + nqf3

   kq_label = i*nqf2*nqf3 + j*nqf3 + k + 1
   !
!   write(stdout,*) xq(1:3)
!   write(stdout,*) i, j, k, nqf1, nqf2, nqf3
   !
end function kq_label


! -----------------------------------------------------------------------------
subroutine valley_degen(xk, xk_v, ik_v)
! -----------------------------------------------------------------------------
   use kinds,       only : DP
   use pwcom,       only : bg
   USE io_files,      ONLY : prefix
   USE epwcom,        ONLY : alloy_pot
   use elph2,       only : cbnd_emin_nxk, cbnd_emin_xk
   USE io_global,   ONLY :stdout
   !
   implicit none
   !
   real(kind=DP)  :: xk(3), xk_v(3)
   integer        :: ik_v

   integer        :: ik, dk1, dk2, dk3
   real(kind=DP)  :: dxk(3), l_dxk, minl_dxk, dxk2(3), &
                     minl_xk
   !
 IF ((TRIM(prefix) .eq. 'Si') .or. (TRIM(prefix) .eq. 'Ge')) THEN   ! only for Qian's SiGe alloys

IF (SQRT((xk(1)-0.5d0)**2 + (xk(2)-0.5d0)**2+ (xk(3)-0.5d0)**2) .LT. 0.0625) THEN
xk_v = (/0.5,0.5,0.5/)
ik_v =1
ENDIF
IF (SQRT((xk(1)-0.5d0)**2 + xk(2)**2+ xk(3)**2) .LT. 0.0625) THEN
xk_v = (/0.5,0.0,0.0/)
ik_v =2
ENDIF
IF (SQRT(xk(1)**2 + (xk(2)-0.5d0)**2+ xk(3)**2) .LT. 0.0625) THEN
xk_v = (/0.0,0.5,0.0/)
ik_v =3
ENDIF
IF (SQRT(xk(1)**2 + xk(2)**2+ (xk(3)-0.5d0)**2) .LT. 0.0625) THEN
xk_v = (/0.0,0.0,0.5/)
ik_v =4
ENDIF
IF (SQRT((xk(1)-0.5d0)**2 + (xk(2)-1.0d0)**2+ (xk(3)-1.0d0)**2) .LT. 0.0625) THEN
xk_v = (/0.5,1.0,1.0/)
ik_v =2
ENDIF
IF (SQRT((xk(1)-1.0d0)**2 + (xk(2)-0.5d0)**2+ (xk(3)-1.0d0)**2) .LT. 0.0625) THEN
xk_v = (/1.0,0.5,1.0/)
ik_v =3
ENDIF
IF (SQRT((xk(1)-1.0d0)**2 + (xk(2)-1.0d0)**2+ (xk(3)-0.5d0)**2) .LT. 0.0625) THEN
xk_v = (/1.0,1.0,0.5/)
ik_v =4
ENDIF
IF (SQRT((xk(1)-0.5d0)**2 + (xk(2)-1.0d0)**2+ xk(3)**2) .LT. 0.0625)  THEN
xk_v = (/0.5,1.0,0.0/)
ik_v =2
ENDIF
IF (SQRT((xk(1)-1.0d0)**2 + (xk(2)-0.5d0)**2+ xk(3)**2) .LT. 0.0625)  THEN
xk_v = (/1.0,0.5,0.0/)
ik_v =3
ENDIF
IF (SQRT((xk(1)-1.0d0)**2 + xk(2)**2+ (xk(3)-0.5d0)**2) .LT. 0.0625)  THEN
xk_v = (/1.0,0.0,0.5/)
ik_v =4
ENDIF
IF (SQRT((xk(1)-0.5d0)**2 + xk(2)**2+ (xk(3)-1.0d0)**2) .LT. 0.0625)  THEN
xk_v = (/0.5,0.0,1.0/)
ik_v =2
ENDIF
IF (SQRT(xk(1)**2 + (xk(2)-0.5d0)**2+ (xk(3)-1.0d0)**2) .LT. 0.0625)  THEN
xk_v = (/0.0,0.5,1.0/)
ik_v =3
ENDIF
IF (SQRT(xk(1)**2 + (xk(2)-1.0d0)**2+ (xk(3)-0.5d0)**2) .LT. 0.0625) THEN
xk_v = (/0.0,1.0,0.5/)
ik_v =4
ENDIF
 ! crystal: L 0.5 0.5 0.5, 0/1 0/1 0.5, 0/1 0.5 0/1, 0.5 0/1 0/1

IF (SQRT((xk(1)-0.425d0)**2 + (xk(2)-0.425d0)**2+ xk(3)**2) .LT. 0.0625)  THEN
xk_v = (/0.425,0.425,0.0/)
ik_v =5
ENDIF
IF (SQRT((xk(1)-0.425d0)**2 + (xk(3)-0.425d0)**2+ xk(2)**2) .LT. 0.0625)  THEN
xk_v = (/0.425,0.0,0.425/)
ik_v =6
ENDIF
IF (SQRT((xk(3)-0.425d0)**2 + (xk(2)-0.425d0)**2+ xk(1)**2) .LT. 0.0625)  THEN
xk_v = (/0.0,0.425,0.425/)
ik_v =7
ENDIF
IF (SQRT((xk(1)-0.425d0)**2 + (xk(2)-0.425d0)**2+ (xk(3)-1.0d0)**2) .LT. 0.0625) THEN
xk_v = (/0.425,0.425,1.0/)
ik_v =5
ENDIF
IF (SQRT((xk(1)-0.425d0)**2 + (xk(3)-0.425d0)**2+ (xk(2)-1.0d0)**2) .LT. 0.0625) THEN
xk_v = (/0.425,1.0,0.425/)
ik_v =6
ENDIF
IF (SQRT((xk(3)-0.425d0)**2 + (xk(2)-0.425d0)**2+ (xk(1)-1.0d0)**2) .LT. 0.0625) THEN
xk_v = (/1.0,0.425,0.425/)
ik_v =7
ENDIF
IF (SQRT((xk(1)-0.575d0)**2 + (xk(2)-0.575d0)**2+ xk(3)**2) .LT. 0.0625)  THEN
xk_v = (/0.575,0.575,0.0/)
ik_v =8
ENDIF
IF (SQRT((xk(1)-0.575d0)**2 + (xk(3)-0.575d0)**2+ xk(2)**2) .LT. 0.0625)  THEN
xk_v = (/0.575,0.0,0.575/)
ik_v =9
ENDIF
IF (SQRT((xk(3)-0.575d0)**2 + (xk(2)-0.575d0)**2+ xk(1)**2) .LT. 0.0625)  THEN
xk_v = (/0.0,0.575,0.575/)
ik_v =10
ENDIF
IF (SQRT((xk(1)-0.575d0)**2 + (xk(2)-0.575d0)**2+ (xk(3)-1.0d0)**2) .LT. 0.0625) THEN
xk_v = (/0.575,0.575,1.0/)
ik_v =8
ENDIF
IF (SQRT((xk(1)-0.575d0)**2 + (xk(3)-0.575d0)**2+ (xk(2)-1.0d0)**2) .LT. 0.0625) THEN
xk_v = (/0.575,1.0,0.575/)
ik_v =9
ENDIF
IF (SQRT((xk(3)-0.575d0)**2 + (xk(2)-0.575d0)**2+ (xk(1)-1.0d0)**2) .LT. 0.0625) THEN
xk_v = (/1.0,0.575,0.575/)
ik_v =10
ENDIF
  ! Si crystal: near X 0/1 0.425 0.425 or 0/1 0.575 0.575

IF (SQRT(xk(1)**2+ xk(2)**2 + xk(3)**2) .LT. 0.0625)  THEN
xk_v = (/0.0,0.0,0.0/)
ik_v =11
ENDIF
IF (SQRT((xk(1)-1.0d0)**2+ (xk(2)-1.0d0)**2 + (xk(3)-1.0d0)**2) .LT. 0.0625) THEN
xk_v = (/1.0,1.0,1.0/)
ik_v =11
ENDIF
 !gamma

 ELSE
   !
   minl_xk = 1d8
   do ik = 1, cbnd_emin_nxk
      dxk(:) = xk(:) - cbnd_emin_xk(:,ik)
    ! WRITE (*,'(3x,a,3f8.2)')  'cbnd_emin_xk(:,ik) =', cbnd_emin_xk(1:3,ik)
    ! WRITE (*,'(3x,a,i4,a,3f8.2)')  'ik =', ik, ', dxk =', dxk(1:3)
      dxk(1) = mod(dxk(1), 1.0d0)
      dxk(2) = mod(dxk(2), 1.0d0)
      dxk(3) = mod(dxk(3), 1.0d0)

      minl_dxk = 1d8
      do dk1 = -2, 2
         do dk2 = -2, 2
            do dk3 = -2, 2
               dxk2 = dxk + (/dk1, dk2, dk3/)
               CALL cryst_to_cart ( 1, dxk2, bg, 1 )   ! xq_cart in [2*pi/alat],
               l_dxk = sqrt(dxk2(1)**2.d0 + dxk2(2)**2.d0 + dxk2(3)**2.d0)
               if (l_dxk < minl_dxk) then
                  minl_dxk = l_dxk
               endif
            enddo
         enddo
      enddo

      if (minl_dxk < minl_xk) then
         minl_xk = minl_dxk
         ik_v = ik
      endif
   enddo
   !
   xk_v = cbnd_emin_xk(:,ik_v)
   ENDIF
   !
end subroutine valley_degen
