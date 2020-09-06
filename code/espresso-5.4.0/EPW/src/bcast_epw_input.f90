  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                        
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the ionode_id, world_comm directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! Code adapted from PH/bcast_ph_input - Quantum-ESPRESSO group
  ! 09/2009 Very little of this subroutine in necessary.  Many 
  ! excess variables
  !
  !-----------------------------------------------------------------------
  SUBROUTINE bcast_ph_input
  !-----------------------------------------------------------------------
  !
  !     In this routine the first processor sends the input to all
  !     the other processors
  !
  !
#ifdef __PARA
#include "f_defs.h"

  USE phcom,         ONLY : zue, trans, tr2_ph, recover, nmix_ph, niter_ph, &
                            lnscf, ldisp, fildvscf, fildrho, epsil, alpha_mix 
  USE epwcom,        ONLY : epexst, epbwrite, ep_coupling, eminabs, emaxabs, &
                            elinterp, eliashberg, elecselfen, eig_read, &
                            efermi_read, dvscf_dir, deltaeabs, delta_smear, &
                            delta_qsmear, degaussw, degaussq, conv_thr_raxis, &
                            conv_thr_racon, conv_thr_iaxis, broyden_ndim, &
                            broyden_beta, band_plot, a2f, lacon, &
                            kmaps, kerwrite, kerread, indabs, imag_read, &
                            gap_edge, fsthick, filukq, filukk, filqf, filkf, &
                            filelph, fileig, fildvscf0, fila2f, fermi_energy, &
                            etf_mem, epwwrite, epwread, eptemp, epstrict, &
                            eps_acustic, ephwrite, epbread, nsiter, nqstep, &
                            nqsmear, nqf3, nqf2, nqf1, nkf3, nkf2, nkf1, &
                            ngaussw, nest_fn, neptemp, nbndsub, nbndskip, &
                            muc, mp_mesh_q, mp_mesh_k, max_memlt, lunif, &
                            lreal, lpolar, lpade, liso, limag, laniso, &
                            specfun, selfen_type, &
                            rand_q, rand_nq, rand_nk, rand_k, pwc, phonselfen, &
                            phinterp, parallel_q, parallel_k, &
                            nw_specfun, nw, nswi, nswfc, nswc, nstemp, nsmear, &
                            wsfc, wscut, write_wfn, wmin_specfun, wmin, &
                            wmax_specfun, wmax, wepexst, wannierize, &
                            vme, twophoton, tshuffle2, tshuffle, tphases, &
                            tempsmin, tempsmax, temps, delta_approx, &
                            ! THL
                            egap_rbm, edos_read, epdope, nepdope, smearing, epdim, &
                            bte, asr_eph, ph_read, ephl_read, dt, run, efield, gradt, mixing, &
                            nptype, save_m_mat, save_m_matw, save_m_ph, save_t_el, nkfdos1, nkfdos2, nkfdos3, &
                            vg_el, vg_ph, relax_time, fscheck, ifc_read, &
                            epthick, phdrag, phwmax, phkmax, alloy_read, &
                            shengbte_read, alloy_pot, frac_type, &
                            ! JWZ
                            eimp_mode, bte_o, dielec, imp_charge, eimpbread, &
                            eimpbwrite, defectname, dvimpsr, defectdens, &
                            epcheck, nbnd_c, nbnd_v, xk_c, xk_v, epbjump, &
                            epbrestore, eimp_sr, elop, screen_polar, dvimpq, &
                            eimp_ls_mode, eph_interp

  USE elph2,         ONLY : elph 
  USE mp,            ONLY : mp_bcast
  USE mp_world,      ONLY : world_comm
  USE io_files,      ONLY : prefix, tmp_dir
  USE qpoint,        ONLY : xq
  USE control_lr,    ONLY : lgamma
  USE io_global,     ONLY : ionode_id
  USE control_flags, ONLY : iverbosity
  USE ions_base,     ONLY : amass
  USE printout_base, ONLY : title   ! title of the run


  implicit none
  !
  ! THL
  ! real
  CALL mp_bcast (egap_rbm, ionode_id, world_comm)
  CALL mp_bcast (epdope, ionode_id, world_comm)
  CALL mp_bcast (dt, ionode_id, world_comm)
  CALL mp_bcast (efield, ionode_id, world_comm)
  CALL mp_bcast (gradt, ionode_id, world_comm)
  CALL mp_bcast (relax_time, ionode_id, world_comm)
  CALL mp_bcast (mixing, ionode_id, world_comm)
  CALL mp_bcast (dielec, ionode_id, world_comm)
  CALL mp_bcast (imp_charge, ionode_id, world_comm)
  CALL mp_bcast (phwmax, ionode_id, world_comm)
  CALL mp_bcast (phkmax, ionode_id, world_comm)
  CALL mp_bcast (epthick, ionode_id, world_comm)
  ! integer
  CALL mp_bcast (bte, ionode_id, world_comm)
  CALL mp_bcast (run, ionode_id, world_comm)
  CALL mp_bcast (nkfdos1, ionode_id, world_comm)
  CALL mp_bcast (nkfdos2, ionode_id, world_comm)
  CALL mp_bcast (nkfdos3, ionode_id, world_comm)
  CALL mp_bcast (epdim, ionode_id, world_comm)
  CALL mp_bcast (nepdope, ionode_id, world_comm)
  CALL mp_bcast (frac_type, ionode_id, world_comm)
  ! logical
  CALL mp_bcast (edos_read, ionode_id, world_comm)
  CALL mp_bcast (asr_eph, ionode_id, world_comm)
  CALL mp_bcast (ph_read, ionode_id, world_comm)
  CALL mp_bcast (ephl_read, ionode_id, world_comm)
  CALL mp_bcast (ifc_read, ionode_id, world_comm)
  CALL mp_bcast (save_m_mat, ionode_id, world_comm)
  CALL mp_bcast (save_m_matw, ionode_id, world_comm)
  CALL mp_bcast (save_m_ph, ionode_id, world_comm)
  CALL mp_bcast (save_t_el, ionode_id, world_comm)
  CALL mp_bcast (fscheck, ionode_id, world_comm)
  CALL mp_bcast (bte_o, ionode_id, world_comm)
  CALL mp_bcast (elop, ionode_id, world_comm)
  CALL mp_bcast (phdrag, ionode_id, world_comm)
  CALL mp_bcast (shengbte_read, ionode_id, world_comm)
  CALL mp_bcast (alloy_read, ionode_id, world_comm)
  CALL mp_bcast (alloy_pot, ionode_id, world_comm)
  ! char
  CALL mp_bcast (nptype, ionode_id, world_comm)
  CALL mp_bcast (smearing, ionode_id, world_comm)
  CALL mp_bcast (vg_el, ionode_id, world_comm)
  CALL mp_bcast (vg_ph, ionode_id, world_comm)
  !
  ! logicals
  !
  CALL mp_bcast (lgamma, ionode_id, world_comm)
  CALL mp_bcast (epsil, ionode_id, world_comm)
  CALL mp_bcast (trans, ionode_id, world_comm)
  CALL mp_bcast (zue, ionode_id, world_comm)
  CALL mp_bcast (elph, ionode_id, world_comm)
  CALL mp_bcast (lnscf, ionode_id, world_comm)
  CALL mp_bcast (ldisp, ionode_id, world_comm)
  CALL mp_bcast (tshuffle, ionode_id, world_comm)  ! 
  CALL mp_bcast (tshuffle2, ionode_id, world_comm) !
  CALL mp_bcast (elecselfen, ionode_id, world_comm)!
  CALL mp_bcast (phonselfen, ionode_id, world_comm)!
  CALL mp_bcast (ephwrite, ionode_id, world_comm)! RM
  CALL mp_bcast (band_plot, ionode_id, world_comm)! RM
  CALL mp_bcast (vme, ionode_id, world_comm)!
  CALL mp_bcast (recover, ionode_id, world_comm)!
  CALL mp_bcast (epbread, ionode_id, world_comm)   !
  CALL mp_bcast (epbwrite, ionode_id, world_comm)  !
  CALL mp_bcast (phinterp, ionode_id, world_comm)  !
  CALL mp_bcast (elinterp, ionode_id, world_comm)  !
  CALL mp_bcast (tphases, ionode_id, world_comm)   !
  CALL mp_bcast (epstrict, ionode_id, world_comm)  !
  CALL mp_bcast (fsthick, ionode_id, world_comm)   !
  CALL mp_bcast (eptemp, ionode_id, world_comm)    !
  CALL mp_bcast (wmin, ionode_id, world_comm)      !
  CALL mp_bcast (wmax, ionode_id, world_comm)      !
  CALL mp_bcast (epwread, ionode_id, world_comm)   !
  CALL mp_bcast (epwwrite, ionode_id, world_comm)  !
  CALL mp_bcast (specfun, ionode_id, world_comm)   !
  CALL mp_bcast (wannierize, ionode_id, world_comm)! JN
  CALL mp_bcast (write_wfn, ionode_id, world_comm) ! 
  CALL mp_bcast (kmaps, ionode_id, world_comm) ! 
  CALL mp_bcast (nest_fn, ionode_id, world_comm) ! 
  CALL mp_bcast (indabs, ionode_id, world_comm) ! 
  CALL mp_bcast (twophoton, ionode_id, world_comm) ! 
  CALL mp_bcast (eig_read, ionode_id, world_comm) ! 
  CALL mp_bcast (parallel_k, ionode_id, world_comm) 
  CALL mp_bcast (parallel_q, ionode_id, world_comm)
  CALL mp_bcast (a2f, ionode_id, world_comm)
  CALL mp_bcast (etf_mem, ionode_id, world_comm)
  CALL mp_bcast (rand_q, ionode_id, world_comm)
  CALL mp_bcast (rand_k, ionode_id, world_comm)
  CALL mp_bcast (mp_mesh_q, ionode_id, world_comm)
  CALL mp_bcast (mp_mesh_k, ionode_id, world_comm)
  CALL mp_bcast (wepexst, ionode_id, world_comm)
  CALL mp_bcast (epexst, ionode_id, world_comm)
  CALL mp_bcast (lreal, ionode_id, world_comm)     ! RM
  CALL mp_bcast (limag, ionode_id, world_comm)     !
  CALL mp_bcast (lpade, ionode_id, world_comm)     !  
  CALL mp_bcast (lacon, ionode_id, world_comm)     !
  CALL mp_bcast (liso, ionode_id, world_comm)     !
  CALL mp_bcast (laniso, ionode_id, world_comm)     !
  CALL mp_bcast (lpolar, ionode_id, world_comm)     !
  CALL mp_bcast (lunif, ionode_id, world_comm)     !
  CALL mp_bcast (kerwrite, ionode_id, world_comm)     !
  CALL mp_bcast (kerread, ionode_id, world_comm)     !
  CALL mp_bcast (imag_read, ionode_id, world_comm ) !
  CALL mp_bcast (eliashberg, ionode_id, world_comm ) !
  CALL mp_bcast (ep_coupling, ionode_id, world_comm ) !
  CALL mp_bcast (efermi_read, ionode_id, world_comm)
  CALL mp_bcast (wmin_specfun, ionode_id, world_comm)      !
  CALL mp_bcast (wmax_specfun, ionode_id, world_comm)      !
  CALL mp_bcast (delta_approx, ionode_id, world_comm)      !
  CALL mp_bcast (eimpbread, ionode_id, world_comm)  !
  CALL mp_bcast (eimpbwrite, ionode_id, world_comm)  !
  CALL mp_bcast (dvimpsr, ionode_id, world_comm)  !
  CALL mp_bcast (dvimpq, ionode_id, world_comm)  !
  CALL mp_bcast (epbjump, ionode_id, world_comm)  !
  CALL mp_bcast (epbrestore, ionode_id, world_comm)  !
  CALL mp_bcast (eimp_sr, ionode_id, world_comm)  !
  CALL mp_bcast (screen_polar, ionode_id, world_comm)  !
  CALL mp_bcast (eph_interp, ionode_id, world_comm)  !
  !
  ! integers
  !
  CALL mp_bcast (niter_ph, ionode_id, world_comm)
  CALL mp_bcast (nmix_ph, ionode_id, world_comm)
  CALL mp_bcast (iverbosity, ionode_id, world_comm)
  CALL mp_bcast (ngaussw, ionode_id, world_comm)     ! FG
  CALL mp_bcast (nw, ionode_id, world_comm)          ! 
  CALL mp_bcast (selfen_type, ionode_id, world_comm) ! 
  CALL mp_bcast (nbndsub, ionode_id, world_comm)     ! 
  CALL mp_bcast (nbndskip, ionode_id, world_comm)    ! 
  CALL mp_bcast (nsmear, ionode_id, world_comm)      ! 
  CALL mp_bcast (rand_nq, ionode_id, world_comm)     ! 
  CALL mp_bcast (rand_nk, ionode_id, world_comm)     ! 
  CALL mp_bcast (nkf1, ionode_id, world_comm)
  CALL mp_bcast (nkf2, ionode_id, world_comm)
  CALL mp_bcast (nkf3, ionode_id, world_comm)
  CALL mp_bcast (nqf1, ionode_id, world_comm)
  CALL mp_bcast (nqf2, ionode_id, world_comm)
  CALL mp_bcast (nqf3, ionode_id, world_comm)
  CALL mp_bcast (nqsmear, ionode_id, world_comm )    ! 
  CALL mp_bcast (nqstep, ionode_id, world_comm)      ! 
  CALL mp_bcast (neptemp, ionode_id, world_comm)     !
  CALL mp_bcast (nswfc, ionode_id, world_comm )      ! 
  CALL mp_bcast (nswc, ionode_id, world_comm )       !
  CALL mp_bcast (nswi, ionode_id, world_comm )       !
  CALL mp_bcast (broyden_ndim, ionode_id, world_comm)!
  CALL mp_bcast (nstemp, ionode_id, world_comm )     !
  CALL mp_bcast (nsiter, ionode_id, world_comm )     !
  CALL mp_bcast (nw_specfun, ionode_id, world_comm)  !
  CALL mp_bcast (eimp_mode, ionode_id, world_comm)  !
  CALL mp_bcast (eimp_ls_mode, ionode_id, world_comm)  !
  !
  ! real*8
  !
  CALL mp_bcast (tr2_ph, ionode_id, world_comm)
  CALL mp_bcast (amass, ionode_id, world_comm)
  CALL mp_bcast (alpha_mix, ionode_id, world_comm)
  CALL mp_bcast (xq, ionode_id, world_comm)
  CALL mp_bcast (degaussw, ionode_id, world_comm)  ! FG
  CALL mp_bcast (delta_smear, ionode_id, world_comm)    ! 
  CALL mp_bcast (eminabs, ionode_id, world_comm)    ! 
  CALL mp_bcast (emaxabs, ionode_id, world_comm)    ! 
  CALL mp_bcast (deltaeabs, ionode_id, world_comm)    ! 
  CALL mp_bcast (eps_acustic, ionode_id, world_comm)     ! RM
  CALL mp_bcast (degaussq, ionode_id, world_comm)        !
  CALL mp_bcast (delta_qsmear, ionode_id, world_comm)    ! 
  CALL mp_bcast (pwc, ionode_id, world_comm )            !
  CALL mp_bcast (wsfc, ionode_id, world_comm )           !
  CALL mp_bcast (wscut, ionode_id, world_comm )          !
  CALL mp_bcast (broyden_beta, ionode_id, world_comm )   !
  CALL mp_bcast (tempsmin, ionode_id, world_comm )       !
  CALL mp_bcast (tempsmax, ionode_id, world_comm )       !
  CALL mp_bcast (temps, ionode_id, world_comm )       !
  CALL mp_bcast (conv_thr_raxis, ionode_id, world_comm ) !
  CALL mp_bcast (conv_thr_iaxis, ionode_id, world_comm ) !
  CALL mp_bcast (conv_thr_racon, ionode_id, world_comm ) !
  CALL mp_bcast (gap_edge, ionode_id, world_comm ) !
  CALL mp_bcast (muc, ionode_id, world_comm )            !
  CALL mp_bcast (max_memlt, ionode_id, world_comm)       !
  CALL mp_bcast (fermi_energy, ionode_id, world_comm)    !
  !
  ! characters
  !
  CALL mp_bcast (title, ionode_id, world_comm)
  CALL mp_bcast (filelph, ionode_id, world_comm)
  CALL mp_bcast (fildvscf, ionode_id, world_comm)
  CALL mp_bcast (fildrho, ionode_id, world_comm)
  CALL mp_bcast (tmp_dir, ionode_id, world_comm)
  CALL mp_bcast (prefix, ionode_id, world_comm)
  !
  CALL mp_bcast (filkf, ionode_id, world_comm)     ! FG
  CALL mp_bcast (filqf, ionode_id, world_comm)     ! FG
  CALL mp_bcast (filukk, ionode_id, world_comm)    ! FG
  CALL mp_bcast (filukq, ionode_id, world_comm)    ! FG
  CALL mp_bcast (fileig, ionode_id, world_comm)    ! FG
  CALL mp_bcast (fildvscf0, ionode_id, world_comm) !
  CALL mp_bcast (dvscf_dir, ionode_id, world_comm)
  CALL mp_bcast (fila2f, ionode_id, world_comm)     ! RM

  CALL mp_bcast (defectname, ionode_id, world_comm)
  CALL mp_bcast (defectdens, ionode_id, world_comm)
  !
  ! ep/eimp check, JWZ
  CALL mp_bcast (epcheck, ionode_id, world_comm)
  CALL mp_bcast (nbnd_c,  ionode_id, world_comm)
  CALL mp_bcast (nbnd_v,  ionode_id, world_comm)
  CALL mp_bcast (xk_c,    ionode_id, world_comm)
  CALL mp_bcast (xk_v,    ionode_id, world_comm)
  
  !
#endif
  !
END SUBROUTINE bcast_ph_input
!
!-----------------------------------------------------------------------
SUBROUTINE bcast_ph_input1
  !-----------------------------------------------------------------------
  !
#ifdef __PARA
#include "f_defs.h"

  USE pwcom
  USE phcom
  USE mp,         ONLY: mp_bcast
  USE mp_world,   ONLY : world_comm
  USE io_global,  ONLY : ionode_id
  implicit none

  !
  ! integers
  !
  CALL mp_bcast (nat_todo, ionode_id, world_comm)
  IF (nat_todo.gt.0) THEN
     CALL mp_bcast (atomo, ionode_id, world_comm)
  ENDIF
#endif
  !  
END SUBROUTINE bcast_ph_input1
