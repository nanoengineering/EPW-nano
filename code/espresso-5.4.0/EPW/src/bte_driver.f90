!----------------------------------------------------------------------------
SUBROUTINE driver_bte_el ()
!----------------------------------------------------------------------------
  USE kinds,      ONLY : DP
  USE cell_base,  ONLY : at, bg
  USE phcom,      ONLY : nmodes
  USE epwcom,     ONLY : nbndsub, bte, nkf1, nkf2, nkf3, nqf1, nqf2, nqf3, neptemp, nepdope, mixing, phdrag, &
                         eimp_mode, bte_o, smearing, alloy_pot, eptemp, epdope
  USE elph2,      ONLY : ibndmin, ibndmax, nbnd_red, &
                         sigmai_mode_all_abs, sigmai_mode_all_emi, vel_ful, &
                         F_k_ful, F_k_ful_tmp, tau0_k_ful, mfp_k_ful, mfp_q_ful, wf_all, vph_ful, &
                         sigmai_mode_all_ela_intra, etf_all, vel_all, ef_epw, sigmai_mode_all_ela_inter, &
                         sigmai_mode_all_alloy_inter, sigmai_mode_all_alloy_intra, &
                         sigmai_mode_all_inter, sigmai_mode_all_intra
  USE bte_var
  USE epw_explore
  USE para_thl
  USE constants_epw, ONLY : ryd2thz, au2ps, au2m, au2s, au2nm, ryd2ev
#ifdef __PARA
  USE io_global,  ONLY : ionode_id, stdout
  USE mp,         ONLY : mp_barrier, mp_sum
  USE mp_world,   ONLY : world_comm, mpime
  USE mp_global,  ONLY : my_pool_id, inter_pool_comm
  use cell_base,  only : omega
  use io_files,   only : prefix
#ENDIF
  !
  IMPLICIT NONE
  !
  !
  REAL(KIND=DP), PARAMETER   :: conv_thr = 1.0d-3
  !
  REAL(KIND=DP), ALLOCATABLE :: F_k_pol(:,:,:), N_q_pol(:,:,:)
  REAL(KIND=DP)              :: transpt(3,3,24), mob_0(3), mob_1(3), spd_0(3), spd_1(3), &
                                eptemp0, vel_k(3), vel0, tau0, tau_abs(nmodes), &
                                tau_emi(nmodes), tau_eimp, tau_eimp_intra, tau_eimp_inter, tau_alloy, &
                                scat0, scat_abs, scat_emi, scat_inter, scat_intra, &
                                scat_eimp, scat_eimp_intra, scat_eimp_inter, scat_alloy, scat_alloy_inter, &
                                scat_alloy_intra
  INTEGER                    :: ik, ik0, ibnd, imode, ir, itemp, idope, iq, iq0, &
                                ik_ful_red, ik_irr, ik_ful, ibnd0, nxk
  INTEGER                    :: nk_pol, ik_pol, ik_star, ik_stop, nq_pol, iq_pol, iq_star, iq_stop
  LOGICAL                    :: converge = .FALSE., file_exist
  INTEGER                    :: iter = 0, iter_rec = 0, temp_rec = 1, dope_rec = 1
  ! date and time
  CHARACTER(LEN=8)           :: date_
  CHARACTER(LEN=10)          :: time_
  CHARACTER(LEN=5)           :: zone_
  INTEGER                    :: values_(8)
  REAL(KIND=DP)              :: t0, t1, ta, tb, tab, ta_phd, tb_phd, tab_phd
  character(len=256) :: tempfilek_fbz, filint_k_fbz, chr_deg, chr_temp, chr_ef, &
                        tempfilek_ibz, filint_k_ibz, file_scat_rate

  !
  !
  mkq1 = nkf1/nqf1
  mkq2 = nkf2/nqf2
  mkq3 = nkf3/nqf3
  !
  ! load the needed file
  CALL kq_load ()
  CALL kq_red_load ()
  CALL meta_load ()
  !
  IF (channel .EQ. .TRUE.) CALL export_channel ()
  !
  ! generate new k list in red-ful-BZ according to the number of scattering events
  CALL cpu_index ()
  !
  ! rotate and copy group velocity from irr- to red-ful-BZ
  CALL rotate_vel ()
  !
  ! initialization of paralellization
  CALL para_bounds (ik_star, ik_stop, nk_ful_red)
  nk_pol = ik_stop-ik_star+1
  !
  IF (phdrag) THEN
     !
     CALL rotate_vph ()
     !
     CALL para_bounds (iq_star, iq_stop, nq_ful_red)
     nq_pol = iq_stop-iq_star+1
     !
  ENDIF
  !
  !
  ! load F_k_ful
  IF (bte .EQ. 10) THEN
     !
     ALLOCATE (F_k_ful(3,nbnd_red,nk_ful_red), F_k_ful_tmp(3,nbnd_red,nk_ful_red))
     F_k_ful = 0.0d0
     F_k_ful_tmp = 0.0d0
     INQUIRE(FILE='BTE/META/F_k_ful',EXIST=file_exist)
     IF (file_exist) CALL Fk_load (iter_rec, temp_rec, dope_rec)
     !
  ENDIF
  !
  if (bte_o) then
     tempfilek_ibz = trim(prefix) // '.k_ibz.'
     tempfilek_fbz = trim(prefix) // '.k_fbz.'
  endif
  !
  WRITE (stdout,'(/5x,a)') 'Iterative scheme'
  !
  DO itemp = temp_rec, neptemp
     !
     eptemp0 = eptemp(itemp) / 0.000086173423d0 * ryd2ev
     write(chr_temp,'(f6.1)') eptemp0
     chr_temp = trim(adjustl(chr_temp))
     !
     DO idope = dope_rec, nepdope
        !
        write(chr_ef,'(ES10.2)') epdope(idope)
        chr_ef = trim(adjustl(chr_ef))
        !
        if (bte_o) then
           filint_k_ibz = trim(adjustl(tempfilek_ibz)) // trim(adjustl(smearing)) &
                         // '_T' // trim(adjustl(chr_temp)) &
                         // '_' // trim(adjustl(chr_ef))
           file_scat_rate = 'Scat_rate_' // trim(adjustl(smearing)) &
                         // '_T' // trim(adjustl(chr_temp)) &
                         // '_' // trim(adjustl(chr_ef))
        endif
        ! 
        iter = iter_rec
        converge = .FALSE.
        !
        ALLOCATE (F_k_pol(3,nbnd_red,nk_pol))
        ALLOCATE (tau0_k_ful(nbnd_red,nk_ful_red))
        ALLOCATE (mfp_k_ful(nbnd_red,nk_ful_red))
        F_k_pol    = 0.0d0
        tau0_k_ful = 0.0d0
        mfp_k_ful  = 0.0d0
        !
        IF (.NOT. ALLOCATED(F_k_ful)) THEN
           ALLOCATE (F_k_ful(3,nbnd_red,nk_ful_red))
           F_k_ful = 0.0d0
        ENDIF
        !
        IF (.NOT. ALLOCATED(F_k_ful_tmp)) THEN
           ALLOCATE (F_k_ful_tmp(3,nbnd_red,nk_ful_red))
           F_k_ful_tmp = 0.0d0
        ENDIF
        !
        ! compute phonon-limited electron relaxation time in red-ful-BZ
        CALL bte_tau0 (itemp, idope)
        !
        !
        IF (phdrag) THEN
           !
           ALLOCATE (mfp_q_ful(nmodes,nq_ful_red))  
           ALLOCATE (int_N_q_ful(3,neptemp,nmodes,nq_ful_red))
           ALLOCATE (N_q_ful(3,neptemp,nmodes,nq_ful_red))
           ALLOCATE (dN_q_ful(3,neptemp,nmodes,nq_ful_red))
           ALLOCATE (tau0_q_ful(neptemp,nmodes,nq_ful_red))
           N_q_ful = 0.0d0
           int_N_q_ful = 0.0d0
           dN_q_ful = 0.0d0
           tau0_q_ful = 0.0d0
           !
           CALL bte_tau0_ph (itemp, idope)
        ENDIF
        !
        !
        ! iteration
        WRITE (stdout,'(/8x,a,i2,a,i2,a/)') '[T', itemp, ' | N', idope, ']'
        !
        DO WHILE (.NOT. converge)
           !
           CALL CPU_TIME (t0)
           !
           !
           ! phonon drag
           IF (phdrag) THEN
              !
              N_q_ful = 0.0d0
              dN_q_ful = 0.0d0
           int_N_q_ful = 0.0d0
              !
              tab_phd = 0.0d0
              DO iq0 = iq_star, iq_stop
                 !
                 CALL CPU_TIME (ta_phd)
                 CALL iter_bte_ph (iq0, itemp, idope, iter)
                 CALL CPU_TIME (tb_phd)
                 tab_phd = tab_phd + (tb_phd-ta_phd)
                 !
                 IF (MOD(iq0-iq_star+1,1000) .EQ. 0 .OR. MOD(iq0-iq_star+1,nq_pol) .EQ. 0) THEN
                    WRITE (stdout,'(16x,a,i6,a,i6,a,f7.1,a,i2,a,i2,a,i2)') 'phdrag  |  q = (', iq0-iq_star+1, '/', nq_pol, &
                    ') completed | Time : ', tab_phd, ' s | ', values_(5), ':', values_(6), ':', values_(7)
                    tab_phd = 0.0d0
                 ENDIF
                 !
              ENDDO
              !
#ifdef __PARA
  CALL mp_barrier (world_comm)
  CALL mp_sum (N_q_ful,inter_pool_comm)
  CALL mp_sum (dN_q_ful,inter_pool_comm)
  CALL mp_sum (int_N_q_ful,inter_pool_comm)
#ENDIF
 !  open(12159,file='N.dat',status='replace')
 !  do iq=1,nq_ful_red
 !  do imode=1,nmodes
 !    write(12159,'(4es15.5)') wf_all(imode,rful2rirr_q(iq))*ryd2thz, tau0_q_ful(imode,iq)*au2ps, &
 !                             SQRT(DOT_PRODUCT(vph_ful(:,imode,iq),vph_ful(:,imode,iq)))*(au2m/au2s), &
 !                             SQRT(DOT_PRODUCT(dN_q_ful(:,imode,iq),dN_q_ful(:,imode,iq)))*au2nm
 ! enddo
 ! enddo
 !STOP
              !
           ENDIF
           !
           !
           tab = 0.0d0
           DO ik0 = ik_star, ik_stop
              !
              ! ik0 is the id of kpoint red-ful-BZ ranging from ik_star to ik_stop (totally, from 1 to nk_ful_red)
              !
              ik_pol = ik0-ik_star+1
              !
              CALL CPU_TIME (ta)
              CALL iter_bte_el (ik0, ik_pol, nk_pol, F_k_pol, itemp, idope, iter)
              CALL CPU_TIME (tb)
              tab = tab + (tb-ta)
              !
              IF (MOD(ik0-ik_star+1,1000) .EQ. 0 .OR. MOD(ik0-ik_star+1,nk_pol) .EQ. 0) THEN
                 CALL DATE_AND_TIME (date_,time_,zone_,values_)
                 WRITE (stdout,'(13x,a,i3,a,i6,a,i6,a,f7.1,a,i2,a,i2,a,i2)') 'Iter. ', iter, '  |  k = (', ik0-ik_star+1, '/', nk_pol, &
                 ') completed | Time : ', tab, ' s | ', values_(5), ':', values_(6), ':', values_(7)
                 tab = 0.0d0
              ENDIF
              !
           ENDDO
           !
#ifdef __PARA
  CALL mp_barrier (world_comm)
#ENDIF
           !
           CALL gether_thl (F_k_ful_tmp, F_k_pol, ik_star, ik_stop, seq2nscat)
           !
           IF (iter .EQ. 0) THEN
              F_k_ful = F_k_ful_tmp !RTA
           ELSE
              F_k_ful = mixing*(F_k_ful_tmp-F_k_ful) + F_k_ful 
           ENDIF
           !
#ifdef __PARA
  CALL mp_barrier (world_comm)
#ENDIF
           !
       ! output to file
        IF (iter .GE. 0) THEN
        if (bte_o .and. (mpime == ionode_id)) then
           OPEN (unit = 30000, file = filint_k_ibz)
           OPEN (unit = 4444, file = file_scat_rate)
!           OPEN (unit = 30001, file = filint_k_fbz)

           write (30000,*) trim(adjustl(prefix))
           write (30000,*) omega
           write (30000,*) ef_epw(itemp,idope) * ryd2ev
           write (30000,*) nk_irr_red, nbnd_red, eimp_mode, nk_ful
           do ik_irr = 1, nk_irr_red
              !
              nxk = 0
              do ik_ful = 1, nk_ful_red
                 if (rful2rirr(ik_ful) == ik_irr) then
                    ik_ful_red = ik_ful
                    nxk = nxk + 1
                 endif
              enddo
              DO ibnd = 1, nbnd_red
                 ibnd0 = ibnd + ibndmin - 1

                 vel_k = vel_ful(:,ibnd,ik_ful_red)
!                 vel0 = sqrt(vel_all(1,ibnd0,ik_irr)**2 + vel_all(2,ibnd0,ik_irr)**2 + &
!                             vel_all(3,ibnd0,ik_irr)**2)
                 vel0 = sqrt(vel_k(1)**2 + vel_k(2)**2 + vel_k(3)**2)

                 scat0 = 2*(sum(sigmai_mode_all_abs(itemp,idope,1:nmodes,ibnd,ik_irr)) + &
                        sum(sigmai_mode_all_emi(itemp,idope,1:nmodes,ibnd,ik_irr)))
                 scat_inter = 2*sum(sigmai_mode_all_inter(itemp,idope,1:nmodes,ibnd,ik_irr))
                 scat_intra = 2*sum(sigmai_mode_all_intra(itemp,idope,1:nmodes,ibnd,ik_irr))
                 scat_abs = 2*sum(sigmai_mode_all_abs(itemp,idope,1:nmodes,ibnd,ik_irr))
                 scat_emi = 2*sum(sigmai_mode_all_emi(itemp,idope,1:nmodes,ibnd,ik_irr))

                 tau0 = sum(sigmai_mode_all_abs(itemp,idope,1:nmodes,ibnd,ik_irr)) + &
                        sum(sigmai_mode_all_emi(itemp,idope,1:nmodes,ibnd,ik_irr))
                 tau_abs = sigmai_mode_all_abs(itemp,idope,1:nmodes,ibnd,ik_irr)
                 tau_emi = sigmai_mode_all_emi(itemp,idope,1:nmodes,ibnd,ik_irr)

                 if (eimp_mode > 0) then
                    scat_eimp_intra = 2*sigmai_mode_all_ela_intra(itemp,idope,ibnd,ik_irr)
                    scat_eimp_inter = 2*sigmai_mode_all_ela_inter(itemp,idope,ibnd,ik_irr)
                    scat_eimp = scat_eimp_intra + scat_eimp_inter
                    IF (alloy_pot) THEN
                    scat_alloy = 2*(sigmai_mode_all_alloy_inter(itemp,idope,ibnd,ik_irr) + sigmai_mode_all_alloy_intra(itemp,idope,ibnd,ik_irr))
                    scat_alloy_inter = 2*sigmai_mode_all_alloy_inter(itemp,idope,ibnd,ik_irr)
                    scat_alloy_intra = 2*sigmai_mode_all_alloy_intra(itemp,idope,ibnd,ik_irr)
                    ENDIF

                    tau_eimp_intra = sigmai_mode_all_ela_intra(itemp,idope,ibnd,ik_irr)
                    tau_eimp_inter = sigmai_mode_all_ela_inter(itemp,idope,ibnd,ik_irr)
                    tau_eimp = tau_eimp_intra + tau_eimp_inter
                    IF (alloy_pot) tau_alloy = sigmai_mode_all_alloy_inter(itemp,idope,ibnd,ik_irr) + sigmai_mode_all_alloy_intra(itemp,idope,ibnd,ik_irr) 
                 endif

                 do imode = 1, nmodes
                    if (tau_abs(imode)/=0.0d0) tau_abs(imode) = 0.5d0 / tau_abs(imode)
                    if (tau_emi(imode)/=0.0d0) tau_emi(imode) = 0.5d0 / tau_emi(imode)
                 enddo

                 if (eimp_mode == 0) then
                    if (tau0/=0.0d0) tau0 = 0.5d0 / tau0
                 else
                    tau0 = tau0 + tau_eimp
                    if (tau0/=0.0d0) tau0 = 0.5d0 / tau0

                    if (tau_eimp/=0.0d0) tau_eimp = 0.5d0 / tau_eimp
                    if (tau_eimp_intra/=0.0d0) tau_eimp_intra = 0.5d0 / tau_eimp_intra
                    if (tau_eimp_inter/=0.0d0) tau_eimp_inter = 0.5d0 / tau_eimp_inter
                    IF ((alloy_pot) .and. (tau_alloy/=0.0d0)) tau_alloy = 0.5d0/tau_alloy
                 endif

                 if (eimp_mode == 0) then
                    write(4444, 200) xkf_irr(1:3,rirr2irr(ik_irr)), etf_all(ibnd0,ik_irr) * ryd2ev, &
                 2.0670687d4*scat0, 2.0670687d4*scat_intra, 2.0670687d4*scat_inter, &
                 2.0670687d4*scat_abs, 2.0670687d4*scat_emi

                    write(30000, 100) ik_irr, ibnd0, xkf_irr(1:3,rirr2irr(ik_irr)), nxk, &
                               etf_all(ibnd0,ik_irr) * ryd2ev, &
                               vel_k(1:3), tau0, tau_abs(1:nmodes), tau_emi(1:nmodes)

                 elseif (eimp_mode > 0) then

                    IF (alloy_pot) THEN
                 write(4444, 200) xkf_irr(1:3,rirr2irr(ik_irr)), etf_all(ibnd0,ik_irr) * ryd2ev, &
                 2.0670687d4*scat0, 2.0670687d4*scat_intra, 2.0670687d4*scat_inter, &
                 2.0670687d4*scat_abs, 2.0670687d4*scat_emi, &
                 2.0670687d4*scat_eimp, 2.0670687d4*scat_eimp_intra, 2.0670687d4*scat_eimp_inter, &
                 2.0670687d4*scat_alloy, 2.0670687d4*scat_alloy_intra, 2.0670687d4*scat_alloy_inter

                     write(30000, 100) ik_irr, ibnd0, xkf_irr(1:3,rirr2irr(ik_irr)), nxk, &
                               etf_all(ibnd0,ik_irr) * ryd2ev, &
                               vel_k(1:3), tau0, tau_abs(1:nmodes), tau_emi(1:nmodes), &
                               tau_eimp, tau_eimp_intra, tau_eimp_inter, tau_alloy
                    ELSE
                     write(30000, 100) ik_irr, ibnd0, xkf_irr(1:3,rirr2irr(ik_irr)), nxk, &
                               etf_all(ibnd0,ik_irr) * ryd2ev, &
                               vel_k(1:3), tau0, tau_abs(1:nmodes), tau_emi(1:nmodes), &
                               tau_eimp, tau_eimp_intra, tau_eimp_inter   

                       write(4444, 200) xkf_irr(1:3,rirr2irr(ik_irr)), etf_all(ibnd0,ik_irr) * ryd2ev, &
                 2.0670687d4*scat0, 2.0670687d4*scat_intra, 2.0670687d4*scat_inter, &
                 2.0670687d4*scat_abs, 2.0670687d4*scat_emi, &
                 2.0670687d4*scat_eimp, 2.0670687d4*scat_eimp_intra, 2.0670687d4*scat_eimp_inter
                   
                    ENDIF
                 endif
              ENDDO
           ENDDO
           close(30000)
           close(4444)
        endif
        ENDIF
        !
           !
           ! calculate transport properties
           CALL bte_transpt (transpt, itemp, idope)
           !
           ! check mobility is converged
           IF (iter-iter_rec .EQ. 0) THEN
              !
              DO ir = 1, 3
                 mob_0(ir) = transpt(ir,ir,1)
              ENDDO
              !
              IF (phdrag) THEN
                 !
                 DO ir = 1, 3
                    spd_0(ir) = transpt(ir,ir,9)
                 ENDDO
                 !
              ENDIF
              !
           ELSE
              !
              DO ir = 1, 3
                 mob_1(ir) = transpt(ir,ir,1)
              ENDDO
              IF (MAXVAL(ABS((mob_1(:)-mob_0(:))/mob_0(:))) .LT. conv_thr) converge = .TRUE.
              mob_0 = mob_1
              !
              IF (phdrag) THEN
                 !
                 DO ir = 1, 3
                    spd_1(ir) = transpt(ir,ir,9)
                 ENDDO
                 IF (MAXVAL(ABS((spd_1(:)-spd_0(:))/spd_0(:))) .LT. conv_thr .AND. converge .EQ. .TRUE.) THEN
                    converge = .TRUE.
                 ELSE
                    converge = .FALSE.
                 ENDIF
                 spd_0 = spd_1
                 !
              ENDIF
              !
           ENDIF
           !
           CALL CPU_TIME (t1)
           !
           IF (iter .GE. 1) THEN 
             converge = .TRUE.  ! Forced (might be false) convergence, check the output yourself carefully
           ENDIF
           !
           ! output
           CALL export_mfp (iter, converge, itemp, idope)
           CALL export_rate (iter, converge, itemp, idope)
           CALL export_boltzep (iter, converge, itemp, idope)
           CALL export_cumu (iter, converge, itemp, idope)
           CALL export_fprime (iter, converge, itemp, idope)
           CALL export_onsager (iter, converge, itemp, idope)
           CALL export_result (transpt, iter, converge, itemp, idope, t1-t0)
           !
           ! save the F_k of current iteration as a metafile
           IF (.NOT. converge) THEN
              CALL Fk_save (iter,itemp,idope)
              iter = iter + 1
           ENDIF
           !
           IF (iter .GT. 1) THEN 
                EXIT ! 1 time iteration at maximum
           ENDIF
           !
        ENDDO ! converge 
        !
#ifdef __PARA
  CALL mp_barrier (world_comm)
#ENDIF
        !
        DEALLOCATE (F_k_ful)
        DEALLOCATE (F_k_ful_tmp)
        DEALLOCATE (F_k_pol)
        DEALLOCATE (mfp_k_ful)
        DEALLOCATE (tau0_k_ful)
        IF (phdrag) THEN
           DEALLOCATE (mfp_q_ful)
           DEALLOCATE (N_q_ful)
           DEALLOCATE (int_N_q_ful)
           DEALLOCATE (dN_q_ful)
           DEALLOCATE (tau0_q_ful)
        ENDIF
        !
        iter_rec = 0
        !
     ENDDO ! idope
     !
     dope_rec = 1
     !
  ENDDO ! itemp
  !
  IF (my_pool_id .EQ. ionode_id) CALL SYSTEM ('rm BTE/META/F_k_ful')
! 
100 format(5x,i6,2x,i5,3(1x,f16.10),2x,i4,35(2x,g14.7))
200 format(3(1x,f16.10),2x,35(2x,g14.7))
  !
END SUBROUTINE driver_bte_el



!----------------------------------------------------------------------------
SUBROUTINE driver_bte_rta ()
!----------------------------------------------------------------------------
#INCLUDE "f_defs.h"
  USE kinds,      ONLY : DP
  USE phcom,      ONLY : nmodes
  USE epwcom,     ONLY : nbndsub, bte, nkf1, nkf2, nkf3, nqf1, nqf2, nqf3, neptemp, nepdope, &
                         eptemp, epdope, relax_time, phdrag, eimp_mode, bte_o, smearing, alloy_pot
  USE elph2,      ONLY : ibndmin, ibndmax, nbnd_red, &
                         sigmai_mode_all_abs, sigmai_mode_all_emi, vel_ful, &
                         sigmai_mode_all_inter, sigmai_mode_all_intra, & 
                         F_k_ful, tau0_k_ful, mfp_k_ful, mfp_q_ful, vph_ful, ef_epw, &
                         sigmai_mode_all_ela_inter, sigmai_mode_all_ela_intra, &
                         sigmai_mode_all_alloy_inter, sigmai_mode_all_alloy_intra, etf_all, vel_all
  USE bte_var
  USE epw_explore
  USE constants_epw, ONLY : kB, au2fs, ryd2ev
#ifdef __PARA
  USE io_global,  ONLY : ionode_id, stdout
  USE mp,         ONLY : mp_barrier, mp_sum
  USE mp_world,   ONLY : world_comm, mpime
  USE mp_global,  ONLY : my_pool_id, inter_pool_comm
  use cell_base,  only : omega
  use io_files,   only : prefix
#ENDIF
  !
  IMPLICIT NONE
  !
  !
  REAL(KIND=DP)     :: transpt(3,3,24), &  !zjw
                   eptemp0, vel_k(3), vel0, tau0, tau_abs(nmodes), &
                   tau_emi(nmodes), tau_eimp, tau_eimp_intra, &
                   tau_eimp_inter, tau_alloy, &
                   scat0, scat_abs, scat_emi, scat_inter, scat_intra, &
                   scat_eimp, scat_eimp_intra, scat_eimp_inter, scat_alloy, scat_alloy_inter, &
                   scat_alloy_intra
  INTEGER           :: ik_ful_red, ibnd, imode, itemp, idope, &
                       ik_irr, ik_ful, ibnd0, nxk
  INTEGER           :: iq0, iq_star, iq_stop, nq_pol
  INTEGER           :: iter = 0
  LOGICAL           :: converge = .TRUE.
  ! zjw
  character(len=256) :: tempfilek_fbz, filint_k_fbz, chr_deg, chr_temp, chr_ef, &
                        tempfilek_ibz, filint_k_ibz, file_scat_rate

  ! date and time
  CHARACTER(LEN=8)  :: date_
  CHARACTER(LEN=10) :: time_
  CHARACTER(LEN=5)  :: zone_
  INTEGER           :: values_(8)
  REAL(KIND=DP)     :: t0, t1, ta, tb, tab, ta_phd, tb_phd, tab_phd
  !
  !
  ! bte =  0, RTA
  ! bte = -1, cRTA
  !
  ! load the needed file
  CALL kq_load ()
  CALL kq_red_load ()
  CALL meta_load ()
  !
  IF (channel .EQ. .TRUE.) CALL export_channel ()
  !
  ! rotate and copy group velocity from irr- to red-ful-BZ
  CALL rotate_vel ()
  !
  IF (phdrag) THEN
     !
     mkq1 = nkf1/nqf1 ! must be 1 if phdrag=T
     mkq2 = nkf2/nqf2
     mkq3 = nkf3/nqf3
     !
     CALL rotate_vph ()
     !
     CALL para_bounds (iq_star, iq_stop, nq_ful_red)
     nq_pol = iq_stop-iq_star+1
     !
  ENDIF
  !
  !
  IF (bte .EQ. 0) THEN
     WRITE (stdout,'(/5x,a)') 'Relaxation time approximation'
  ELSEIF (bte .EQ. -1) THEN
     WRITE (stdout,'(/5x,a,f10.3,a)') 'Constant relaxation time approximation (', relax_time, ' fs)'
  ENDIF
  !
  if (bte_o) then
     tempfilek_ibz = trim(prefix) // '.k_ibz.'
     tempfilek_fbz = trim(prefix) // '.k_fbz.'
  endif
  ! zjw
  !
  DO itemp = 1, neptemp
     ! zjw
     eptemp0 = eptemp(itemp) / 0.000086173423d0 * ryd2ev
     write(chr_temp,'(f6.1)') eptemp0
     chr_temp = trim(adjustl(chr_temp))
     !
     DO idope = 1, nepdope
        ! zjw
        write(chr_ef,'(ES10.2)') epdope(idope)
        chr_ef = trim(adjustl(chr_ef))
        !
        if (bte_o) then
           filint_k_ibz = trim(adjustl(tempfilek_ibz)) // trim(adjustl(smearing)) &
                         // '_T' // trim(adjustl(chr_temp)) &
                         // '_' // trim(adjustl(chr_ef))
           file_scat_rate = 'Scat_' // trim(adjustl(smearing)) &
                         // '_T' // trim(adjustl(chr_temp)) &
                         // '_' // trim(adjustl(chr_ef))
        endif
        !
        !
        ! allocate
        ALLOCATE (F_k_ful(3,nbnd_red,nk_ful_red))
        ALLOCATE (tau0_k_ful(nbnd_red,nk_ful_red))
        ALLOCATE (mfp_k_ful(nbnd_red,nk_ful_red))
        F_k_ful    = 0.0d0
        tau0_k_ful = 0.0d0
        mfp_k_ful  = 0.0d0
        !
        ! compute phonon-limited electron relaxation time
        IF (bte .EQ. 0) THEN
           CALL bte_tau0 (itemp,idope)
        ELSEIF (bte .EQ. -1) THEN
           tau0_k_ful = relax_time / au2fs
        ENDIF
        !
        !
        IF (phdrag) THEN
           !
           ALLOCATE (int_N_q_ful(3,neptemp,nmodes,nq_ful_red))
           ALLOCATE (N_q_ful(3,neptemp,nmodes,nq_ful_red))
           ALLOCATE (dN_q_ful(3,neptemp,nmodes,nq_ful_red))
           ALLOCATE (tau0_q_ful(neptemp,nmodes,nq_ful_red))
           ALLOCATE (mfp_q_ful(nmodes,nq_ful_red))
           N_q_ful = 0.0d0
       int_N_q_ful = 0.0d0
           dN_q_ful = 0.0d0
           tau0_q_ful = 0.0d0
           !
           CALL bte_tau0_ph (itemp, idope)
        ENDIF
        !
        !
        CALL CPU_TIME (t0)
        !
        ! compute electron mean free path       
        WRITE (stdout,'(/8x,a,i2,a,i2,a/)') '[T', itemp, ' | N', idope, ']'
        !
        DO ik_ful_red = 1, nk_ful_red
           DO ibnd = 1, nbnd_red
              F_k_ful(:,ibnd,ik_ful_red) = vel_ful(:,ibnd,ik_ful_red) * tau0_k_ful(ibnd,ik_ful_red)
           ENDDO
        ENDDO
        !
        ! zjw, output to file
        if (bte_o .and. (mpime == ionode_id)) then
           OPEN (unit = 30000, file = filint_k_ibz)
           OPEN (unit = 4444, file = file_scat_rate)
!           OPEN (unit = 30001, file = filint_k_fbz)

           write (30000,*) trim(adjustl(prefix))
           write (30000,*) omega
           write (30000,*) ef_epw(itemp,idope) * ryd2ev
           write (30000,*) nk_irr_red, nbnd_red, eimp_mode, nk_ful
           do ik_irr = 1, nk_irr_red
              !
              nxk = 0
              do ik_ful = 1, nk_ful_red
                 if (rful2rirr(ik_ful) == ik_irr) then
                    ik_ful_red = ik_ful
                    nxk = nxk + 1
                 endif
              enddo
              DO ibnd = 1, nbnd_red
                 ibnd0 = ibnd + ibndmin - 1
                 scat0 = 2*(sum(sigmai_mode_all_abs(itemp,idope,1:nmodes,ibnd,ik_irr)) + &
                        sum(sigmai_mode_all_emi(itemp,idope,1:nmodes,ibnd,ik_irr)))
                 scat_inter = 2*sum(sigmai_mode_all_inter(itemp,idope,1:nmodes,ibnd,ik_irr))
                 scat_intra = 2*sum(sigmai_mode_all_intra(itemp,idope,1:nmodes,ibnd,ik_irr))
                 scat_abs = 2*sum(sigmai_mode_all_abs(itemp,idope,1:nmodes,ibnd,ik_irr))
                 scat_emi = 2*sum(sigmai_mode_all_emi(itemp,idope,1:nmodes,ibnd,ik_irr))

                 tau0 = sum(sigmai_mode_all_abs(itemp,idope,1:nmodes,ibnd,ik_irr)) + &
                        sum(sigmai_mode_all_emi(itemp,idope,1:nmodes,ibnd,ik_irr))
                 tau_abs = sigmai_mode_all_abs(itemp,idope,1:nmodes,ibnd,ik_irr)
                 tau_emi = sigmai_mode_all_emi(itemp,idope,1:nmodes,ibnd,ik_irr)

                 if (eimp_mode > 0) then
                    scat_eimp_intra = 2*sigmai_mode_all_ela_intra(itemp,idope,ibnd,ik_irr)
                    scat_eimp_inter = 2*sigmai_mode_all_ela_inter(itemp,idope,ibnd,ik_irr)
                    scat_eimp = scat_eimp_intra + scat_eimp_inter
                    IF (alloy_pot) THEN
                    scat_alloy = 2*(sigmai_mode_all_alloy_inter(itemp,idope,ibnd,ik_irr) + sigmai_mode_all_alloy_intra(itemp,idope,ibnd,ik_irr))
                    scat_alloy_inter = 2*sigmai_mode_all_alloy_inter(itemp,idope,ibnd,ik_irr)
                    scat_alloy_intra = 2*sigmai_mode_all_alloy_intra(itemp,idope,ibnd,ik_irr)
                    ENDIF

                    tau_eimp_intra = sigmai_mode_all_ela_intra(itemp,idope,ibnd,ik_irr)
                    tau_eimp_inter = sigmai_mode_all_ela_inter(itemp,idope,ibnd,ik_irr)
                    tau_eimp = tau_eimp_intra + tau_eimp_inter
                    IF (alloy_pot) tau_alloy = sigmai_mode_all_alloy_inter(itemp,idope,ibnd,ik_irr) + sigmai_mode_all_alloy_intra(itemp,idope,ibnd,ik_irr) 
                 endif

                 do imode = 1, nmodes
                    if (tau_abs(imode)/=0.0d0) tau_abs(imode) = 0.5d0 / tau_abs(imode)
                    if (tau_emi(imode)/=0.0d0) tau_emi(imode) = 0.5d0 / tau_emi(imode)
                 enddo

                 if (eimp_mode == 0) then
                    if (tau0/=0.0d0) tau0 = 0.5d0 / tau0
                 else
                    tau0 = tau0 + tau_eimp
                    if (tau0/=0.0d0) tau0 = 0.5d0 / tau0

                    if (tau_eimp/=0.0d0) tau_eimp = 0.5d0 / tau_eimp
                    if (tau_eimp_intra/=0.0d0) tau_eimp_intra = 0.5d0 / tau_eimp_intra
                    if (tau_eimp_inter/=0.0d0) tau_eimp_inter = 0.5d0 / tau_eimp_inter
                    IF ((alloy_pot) .and. (tau_alloy/=0.0d0)) tau_alloy = 0.5d0/tau_alloy
                 endif

                 if (eimp_mode == 0) then
                    write(4444, 200) xkf_irr(1:3,rirr2irr(ik_irr)), etf_all(ibnd0,ik_irr) * ryd2ev, &
                 2.0670687d4*scat0, 2.0670687d4*scat_intra, 2.0670687d4*scat_inter, &
                 2.0670687d4*scat_abs, 2.0670687d4*scat_emi

                    write(30000, 100) ik_irr, ibnd0, xkf_irr(1:3,rirr2irr(ik_irr)), nxk, &
                               etf_all(ibnd0,ik_irr) * ryd2ev, &
                               vel_k(1:3), tau0, tau_abs(1:nmodes), tau_emi(1:nmodes)

                 elseif (eimp_mode > 0) then

                    IF (alloy_pot) THEN
                 
                 write(4444, 200) xkf_irr(1:3,rirr2irr(ik_irr)), etf_all(ibnd0,ik_irr) * ryd2ev, &
                 2.0670687d4*scat0, 2.0670687d4*scat_intra, 2.0670687d4*scat_inter, &
                 2.0670687d4*scat_abs, 2.0670687d4*scat_emi, &
                 2.0670687d4*scat_eimp, 2.0670687d4*scat_eimp_intra, 2.0670687d4*scat_eimp_inter, 2.0670687d4*scat_alloy, 2.0670687d4*scat_alloy_intra, 2.0670687d4*scat_alloy_inter

                  write(30000, 100) ik_irr, ibnd0, xkf_irr(1:3,rirr2irr(ik_irr)), nxk, &
                               etf_all(ibnd0,ik_irr) * ryd2ev, &
                               vel_k(1:3), tau0, tau_abs(1:nmodes), tau_emi(1:nmodes), &
                               tau_eimp, tau_eimp_intra, tau_eimp_inter, tau_alloy
                    ELSE
                     write(30000, 100) ik_irr, ibnd0, xkf_irr(1:3,rirr2irr(ik_irr)), nxk, &
                               etf_all(ibnd0,ik_irr) * ryd2ev, &
                               vel_k(1:3), tau0, tau_abs(1:nmodes), tau_emi(1:nmodes), &
                               tau_eimp, tau_eimp_intra, tau_eimp_inter   

                     write(4444, 200) xkf_irr(1:3,rirr2irr(ik_irr)), etf_all(ibnd0,ik_irr) * ryd2ev, &
                 2.0670687d4*scat0,2.0670687d4*scat_intra, 2.0670687d4*scat_inter, &
                 2.0670687d4*scat_abs, 2.0670687d4*scat_emi, &                  
                 2.0670687d4*scat_eimp, 2.0670687d4*scat_eimp_intra, 2.0670687d4*scat_eimp_inter
                
                    ENDIF
                 endif
              ENDDO
           ENDDO
           close(30000)
           close(4444)
        endif
        !
        ! phonon drag
        IF (phdrag) THEN
           !
           N_q_ful = 0.0d0
           dN_q_ful = 0.0d0
        int_N_q_ful = 0.0d0
           !
           tab_phd = 0.0d0
           DO iq0 = iq_star, iq_stop
              !
              CALL CPU_TIME (ta_phd)
              CALL iter_bte_ph (iq0, itemp, idope, iter)
              CALL CPU_TIME (tb_phd)
              tab_phd = tab_phd + (tb_phd-ta_phd)
              !
              IF (MOD(iq0-iq_star+1,1000) .EQ. 0 .OR. MOD(iq0-iq_star+1,nq_pol) .EQ. 0) THEN
                 CALL DATE_AND_TIME (date_,time_,zone_,values_)
                 WRITE (stdout,'(16x,a,i6,a,i6,a,f7.1,a,i2,a,i2,a,i2)') 'phdrag  |  q = (', iq0-iq_star+1, '/', nq_pol, &
                 ') completed | Time : ', tab_phd, ' s | ', values_(5), ':', values_(6), ':', values_(7)
                 tab_phd = 0.0d0
              ENDIF
              !
           ENDDO
           !
#ifdef __PARA
  CALL mp_barrier (world_comm)
  CALL mp_sum (int_N_q_ful,inter_pool_comm)
  CALL mp_sum (N_q_ful,inter_pool_comm)
  CALL mp_sum (dN_q_ful,inter_pool_comm)
#ENDIF
           !
        ENDIF
        !
        ! compute electron transport properties
        CALL bte_transpt (transpt, itemp, idope)
        !
        CALL CPU_TIME (t1)  
        !
        ! output
        CALL export_mfp (iter, converge, itemp, idope)
        IF (bte .EQ. 0) CALL export_rate (iter, converge, itemp, idope)
        IF (bte .EQ. 0) CALL export_boltzep (iter, converge, itemp, idope)
        IF (bte .EQ. 0) CALL export_cumu (iter, converge, itemp, idope)
        IF (bte .EQ. 0) CALL export_fprime (iter, converge, itemp, idope)
        IF (bte .EQ. 0) CALL export_onsager (iter, converge, itemp, idope)
        CALL export_result (transpt, iter, converge, itemp, idope, t1-t0)
        !
        !
        DEALLOCATE (F_k_ful)
        DEALLOCATE (mfp_k_ful)
        DEALLOCATE (tau0_k_ful)
        IF (phdrag) THEN
           DEALLOCATE (mfp_q_ful)
           DEALLOCATE (N_q_ful)
           DEALLOCATE (int_N_q_ful)
           DEALLOCATE (dN_q_ful)
           DEALLOCATE (tau0_q_ful)
        ENDIF
        !
     ENDDO
     !
  ENDDO
!
100 format(5x,i6,2x,i5,3(1x,f16.10),2x,i4,35(2x,g14.7))
200 format(3(1x,f16.10),2x,35(2x,g14.7))
  !
END SUBROUTINE driver_bte_rta



!----------------------------------------------------------------------------
SUBROUTINE driver_tdbte_el ()
!----------------------------------------------------------------------------
  USE kinds,      ONLY : DP
  USE cell_base,  ONLY : at, bg
  USE phcom,      ONLY : nmodes
  USE epwcom,     ONLY : nbndsub, bte, nkf1, nkf2, nkf3, nqf1, nqf2, nqf3, run, phdrag, alloy_pot
  USE elph2,      ONLY : ibndmin, ibndmax, nbnd_red, &
                         vel_ful, f_t_ful, f_0_ful
  USE bte_var
  USE para_thl
#ifdef __PARA
  USE io_global,  ONLY : ionode_id, stdout
  USE mp,         ONLY : mp_barrier
  USE mp_global,  ONLY : my_pool_id, inter_pool_comm
#ENDIF
  !
  IMPLICIT NONE
  !
  REAL(KIND=DP)              :: transpt(3,3,24)
  REAL(KIND=DP), ALLOCATABLE :: f_t_pol(:,:)
  INTEGER                    :: ik, ik0, ibnd, imode
  INTEGER                    :: nk_pol, ik_pol, ik_star, ik_stop
  LOGICAL                    :: converge = .FALSE., file_exist
  INTEGER                    :: iter = 0, iter_rec = 0
  ! date and time
  CHARACTER(LEN=8)           :: date_
  CHARACTER(LEN=10)          :: time_
  CHARACTER(LEN=5)           :: zone_
  INTEGER                    :: values_(8)
  REAL(KIND=DP)              :: t0, t1, ta, tb, tab
  !
  !
  mkq1 = nkf1/nqf1
  mkq2 = nkf2/nqf2
  mkq3 = nkf3/nqf3
  !
  ! load the needed file
  CALL kq_load ()
  CALL kq_red_load ()
  CALL meta_load ()
  !
  ! generate new k list
  CALL cpu_index ()
  !
  ! rotate and copy group velocity from irreducible to full BZ
  CALL rotate_vel ()
  !
  ! initialization of paralellization
  CALL para_bounds (ik_star, ik_stop, nk_ful_red)
  nk_pol = ik_stop-ik_star+1
  !
  ALLOCATE (f_t_pol(nbnd_red,nk_pol))
  ALLOCATE (f_t_ful(nbnd_red,nk_ful_red))
  ALLOCATE (f_0_ful(nbnd_red,nk_ful_red))
  f_t_ful = 0.0d0
  f_0_ful = 0.0d0
  f_t_pol = 0.0d0
  !
  ! compute the initial (excited) electron distribution function
  CALL tdbte_f_0 ()
  !
  ! load f_0_ful if BTE/META/f_t_ful exists
  IF (bte .EQ. 30) THEN
     INQUIRE(FILE='BTE/META/f_t_ful',EXIST=file_exist)
     IF (file_exist) CALL ft_load (iter_rec)
     iter = iter_rec
  ENDIF
  !
  ! td-bte
  WRITE (stdout,'(/5x,a/)') 'Start calculation...'
  !
  DO WHILE (.NOT. converge)
     !
     CALL CPU_TIME (t0)
     !
     DO ik0 = ik_star, ik_stop
        !
        ! ik0 is the id of original kpoint in full BZ ranging from ik_star to ik_stop (totally, from 1 to nk_ful)
        ! ik_pol = 1 ~ nk_pol, used in storing F_k_pol
        ik_pol = ik0-ik_star+1
        !
        CALL CPU_TIME (ta)
        CALL iter_tdbte_el (ik0, ik_pol, nk_pol, f_t_pol)
        CALL CPU_TIME (tb)
        tab = tab + (tb-ta)
        !
        IF (MOD(ik0-ik_star+1,1000) .EQ. 0 .OR. MOD(ik0-ik_star+1,nk_pol) .EQ. 0) THEN
           CALL DATE_AND_TIME (date_,time_,zone_,values_)
           WRITE (stdout,'(13x,a,i3,a,i6,a,i6,a,f8.1,a,i2,a,i2,a,i2)') 'Timestep ', iter, '  |  k = (', ik0-ik_star+1, '/', nk_pol, &
           ') completed  |  Time : ', tab, ' s', values_(5), ':', values_(6), ':', values_(7)
           tab = 0.0d0
        ENDIF
        !
     ENDDO
     !
     CALL mp_barrier (inter_pool_comm)
     !
     CALL gether_thl (f_t_ful, f_t_pol, ik_star, ik_stop, seq2nscat)
     !
     ! check simulation step
     IF (iter .EQ. run) converge = .TRUE.
     !
     CALL CPU_TIME (t1)
     !
     ! output
     CALL export_ft (iter, converge, 1, 1)
     CALL export_result (transpt, iter, converge, 1, 1, t1-t0)
     !
     ! save the f_t of current iteration as a metafile
     IF (.NOT. converge) THEN
        !
        CALL ft_save (iter)
        f_0_ful = f_t_ful
        f_t_ful = 0.0d0
        !
     ELSE
        IF (my_pool_id .EQ. ionode_id) CALL SYSTEM ('rm BTE/META/f_t_ful') 
     ENDIF
     !
     iter = iter + 1
     !
  ENDDO ! converge 
  !
  !
  DEALLOCATE (f_0_ful)
  DEALLOCATE (f_t_ful)
  DEALLOCATE (f_t_pol)
  !
END SUBROUTINE driver_tdbte_el



!----------------------------------------------------------------------------
SUBROUTINE driver_bte_random ()
!----------------------------------------------------------------------------
#INCLUDE "f_defs.h"
  USE kinds,      ONLY : DP
  USE phcom,      ONLY : nmodes
  USE epwcom,     ONLY : nbndsub, bte, nkf1, nkf2, nkf3, neptemp, nepdope, eptemp, epdope
  USE elph2,      ONLY : ibndmin, ibndmax, nbnd_red, &
                         sigmai_mode_all_abs, sigmai_mode_all_emi, vel_ful, vel_all, &
                         F_k_ful, tau0_k_ful, mfp_k_ful
  USE bte_var
  USE constants_epw, ONLY : kB, au2fs
#ifdef __PARA
  USE io_global,  ONLY : ionode_id, stdout
  USE mp,         ONLY : mp_barrier
  USE mp_global,  ONLY : my_pool_id, inter_pool_comm
#ENDIF
  !
  IMPLICIT NONE
  !
  !
  REAL(KIND=DP) :: transpt(3,3,24)
  INTEGER       :: ik, ik_irr_red, ibnd, ibnd0, imode, itemp, idope
  REAL(KIND=DP) :: t0, t1
  INTEGER       :: iter = 0
  LOGICAL       :: converge = .TRUE.
  !
  ! load the needed file
  ! xk_kmesh and nk_kmesh are always exist
  nk_ful = nk_kmesh
  nk_irr = nk_kmesh
  nk_irr_red = nk_kmesh
  nk_ful_red = nk_kmesh
  !
  CALL meta_load ()
  !
  ! get vel_ful from vel_all
  ALLOCATE (vel_ful(3,nbnd_red,nk_ful_red))
  vel_ful = 0.0d0
  !
  DO ik = 1, nk_ful_red
     DO ibnd = 1, nbnd_red
        ibnd0 = ibnd+ibndmin-1
        vel_ful(:,ibnd,ik) = vel_all(:,ibnd0,ik)
     ENDDO
  ENDDO
  !
  CALL mp_barrier (inter_pool_comm)
  !
  WRITE (stdout,'(/5x,a)') 'Random k-point calculation'
  !
  DO itemp = 1, neptemp
     !
     DO idope = 1, nepdope
        !
        ! allocate
        ALLOCATE (F_k_ful(3,nbnd_red,nk_ful_red))
        ALLOCATE (tau0_k_ful(nbnd_red,nk_ful_red))
        ALLOCATE (mfp_k_ful(nbnd_red,nk_ful_red))
        F_k_ful    = 0.0d0
        tau0_k_ful = 0.0d0
        mfp_k_ful  = 0.0d0
        !
        CALL CPU_TIME (t0)
        !
        ! compute phonon-limited electron relaxation time
        CALL bte_tau0 (itemp,idope)
        !
        ! compute electron mean free path       
        WRITE (stdout,'(/8x,a,i2,a,i2,a/)') '[T', itemp, ' | N', idope, ']'
        !
        DO ik = 1, nk_ful_red
           DO ibnd = 1, nbnd_red
              F_k_ful(:,ibnd,ik) = vel_ful(:,ibnd,ik) * tau0_k_ful(ibnd,ik)
           ENDDO
        ENDDO
        !
        ! compute electron transport properties
        !CALL bte_transpt (transpt, itemp, idope)
        !
        CALL CPU_TIME (t1)  
        !
        ! output
        IF (bte .EQ. 18) CALL export_channel (itemp, idope)
        CALL export_mfp (iter, converge, itemp, idope)
        CALL export_rate (iter, converge, itemp, idope)
        !CALL export_cumu (iter, converge, itemp, idope)
        !CALL export_onsager (iter, converge, itemp, idope)
        !CALL export_result (transpt, iter, converge, itemp, idope, t1-t0)
        !
        !
        DEALLOCATE (F_k_ful)
        DEALLOCATE (mfp_k_ful)
        DEALLOCATE (tau0_k_ful)
        !
     ENDDO
     !
  ENDDO
  !
END SUBROUTINE driver_bte_random



!----------------------------------------------------------------------------
SUBROUTINE bte_tau0 (itemp, idope)
!----------------------------------------------------------------------------
#INCLUDE "f_defs.h"
  USE kinds,  ONLY : DP
  USE phcom,  ONLY : nmodes
  USE epwcom, ONLY : bte, nbndsub, nkf1, nkf2, nkf3, alloy_read, eimp_mode, alloy_pot
  USE constants_epw, ONLY : ryd2ev
  USE elph2,  ONLY : nbnd_red, sigmai_mode_all_abs, sigmai_mode_all_emi, tau0_k_ful, ibndmin, etf_all, &
                     sigmai_mode_all_ela_intra, sigmai_mode_all_ela_inter, sigmai_mode_all_alloy_inter, &
                     sigmai_mode_all_alloy_intra  !zjw
  USE bte_var
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: itemp, idope
  REAL(KIND=DP)       :: sigmai, energy
  INTEGER             :: ik_ful_red, ibnd, ibnd0, imode, ik_irr_red, i
  !
  !
  DO ik_ful_red = 1, nk_ful_red
     !
     IF (bte .EQ. 19 .OR. bte .EQ. 18) THEN
        ik_irr_red = ik_ful_red
     ELSE
        ik_irr_red = rful2rirr(ik_ful_red)
     ENDIF
     !
  ! etf_all(ibnd0,ik)*ryd2ev
!
     DO ibnd = 1, nbnd_red
        ibnd0 = ibnd+ibndmin-1
        !
        sigmai = 0.0d0
        !
        sigmai =  SUM(sigmai_mode_all_abs(itemp,idope,1:nmodes,ibnd,ik_irr_red)) + &
                  SUM(sigmai_mode_all_emi(itemp,idope,1:nmodes,ibnd,ik_irr_red)) 
        ! zjw 
        if ((eimp_mode > 0) .and. (.not. alloy_pot)) then
           sigmai = sigmai + sigmai_mode_all_ela_intra(itemp,idope,ibnd,ik_irr_red) + &
                             sigmai_mode_all_ela_inter(itemp,idope,ibnd,ik_irr_red)
        endif
        ! al-el method 1: calculated from input perturbed potential dv_tot
        !
        If (alloy_pot)  THEN
                  sigmai = sigmai + sigmai_mode_all_alloy_inter(itemp,idope,ibnd,ik_irr_red) + &
                                sigmai_mode_all_alloy_intra(itemp,idope,ibnd,ik_irr_red) + &
                                sigmai_mode_all_ela_intra(itemp,idope,ibnd,ik_irr_red) + &
                                sigmai_mode_all_ela_inter(itemp,idope,ibnd,ik_irr_red)
        ENDIF
        !
        ! al-el method 2: use input al-el scattering rate
        IF (alloy_read)  THEN
            energy = 0.0d0
            !
            energy = etf_all(ibnd0,ik_irr_red)*ryd2ev
            !
            i = 1
         IF (energy .GE. alelrate(1,1) .AND. energy .LE. alelrate(4000,1)) THEN
            i = int((energy - alelrate(1,1))/1.0d-4)    ! alel_rate.txt energy range 0.4eV with 4000 data points
            !
            IF (i .EQ. 4000) THEN
       sigmai = sigmai + 0.5d0 * alelrate(i,2)
            !
            ELSE
       sigmai = sigmai + 0.5d0 * alelrate(i,2)*(energy - alelrate(i,1))/1.0d-4 + 0.5d0 * alelrate(i+1,2)*(alelrate(i+1,1)-energy)/1.0d-4
            !
            ENDIF  ! i
         ENDIF  ! energy range
        ENDIF   ! alloy_read
        !
        IF (sigmai .NE. 0.0d0) THEN
           tau0_k_ful(ibnd,ik_ful_red) = 0.5d0 / sigmai
        ELSE
           tau0_k_ful(ibnd,ik_ful_red) = 0.0d0
        ENDIF
        !
     ENDDO
     !
  ENDDO
  !
  !
END SUBROUTINE bte_tau0


!----------------------------------------------------------------------------
SUBROUTINE bte_tau0_ph (itemp, idope)
!----------------------------------------------------------------------------
#INCLUDE "f_defs.h"
  USE kinds,  ONLY : DP
  USE phcom,  ONLY : nmodes
  USE epwcom, ONLY : bte
  USE elph2,  ONLY : gammai_mode_all, vph_ful
  USE constants_epw, ONLY : au2ps, ryd2thz
  USE bte_var
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: itemp, idope
  REAL(KIND=DP)       :: gammai
  INTEGER             :: iq_ful_red, imode, iq_irr_red, iq_ful
  !
  !
  DO iq_ful_red = 1, nq_ful_red
     !
     iq_irr_red = rful2rirr_q(iq_ful_red)
     iq_ful     = rful2ful_q(iq_ful_red)
     !
     DO imode = 1, nmodes
        !
        !
        gammai = 0.0d0
        !
        gammai = ph_rate_ful(itemp,imode,iq_ful) + 2.0d0*gammai_mode_all(itemp,idope,imode,iq_irr_red)
        !
        IF (gammai .NE. 0.0d0) THEN
           tau0_q_ful(itemp,imode,iq_ful_red) = 1.0d0 / gammai
        ELSE
           tau0_q_ful(itemp,imode,iq_ful_red) = 0.0d0
        ENDIF
        !
     ENDDO
     !
  ENDDO
  !
  !
END SUBROUTINE bte_tau0_ph




!----------------------------------------------------------------------------
SUBROUTINE tdbte_f_0
!----------------------------------------------------------------------------
#INCLUDE "f_defs.h"
  USE kinds,     ONLY : DP
  USE phcom,     ONLY : nmodes
  USE pwcom,     ONLY : ef 
  USE lsda_mod,  ONLY : nspin
  USE epwcom,    ONLY : nbndsub, eptemp, nkf1, nkf2, nkf3
  USE elph2,     ONLY : nbnd_red, f_0_ful, etf_all, ibndmin, ibndmax, ef_epw
  USE bte_var
#ifdef __PARA
  USE io_files,  ONLY : tmp_dir, prefix
  USE io_global, ONLY : stdout
#ENDIF
  !
  IMPLICIT NONE
  !
  INTEGER                 :: ik_ful_red, ik_irr_red, ibnd, ibnd0
  REAL(KIND=DP)           :: eptemp0, ef0, ekk
  REAL(KIND=DP), EXTERNAL :: wgauss
  !
  !
  ef0 = ef_epw(1,1)
  eptemp0 = eptemp(1)!*3.3333333333333333333d0 ! 1000K

  DO ik_ful_red = 1, nk_ful_red
     !
     ik_irr_red = rful2rirr(ik_ful_red)
     !
     DO ibnd = 1, nbnd_red
        !
        ibnd0 = ibnd+ibndmin-1
        !
        ekk = etf_all(ibnd0,ik_irr_red) - ef0
        !
        f_0_ful(ibnd,ik_ful_red) = wgauss(-ekk/eptemp0,-99)
        !
     ENDDO
     !
  ENDDO
  !
  !
END SUBROUTINE tdbte_f_0
