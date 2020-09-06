!-----------------------------------------------------------------------
SUBROUTINE fermiwindow ()
!-----------------------------------------------------------------------
  !
#INCLUDE "f_defs.h"
  !
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE spin_orb,      ONLY : lspinorb
  USE lsda_mod,      ONLY : nspin
  USE elph2,         ONLY : etf, ibndmin, ibndmax, ebndmin, ebndmax, nbnd_red, cbnd_emin, vbnd_emax, ef_m, delta_egap, &
                            vbnd_num, cbnd_num, cfsthick, vfsthick, nelec_red
  USE epwcom,        ONLY : fsthick, nbndsub, egap_rbm, bte, nptype, nbndskip
  USE pwcom,         ONLY : ef, nelec
  USE constants_epw, ONLY : ryd2ev
  USE bte_var,         ONLY : nk_irr
#ifdef __PARA
  USE mp,            ONLY : mp_max, mp_min, mp_barrier
  USE mp_global,     ONLY : inter_pool_comm
#ENDIF
  !
  IMPLICIT NONE
  !
  INTEGER       :: ik, ibnd
#ifdef __PARA
  REAL(KIND=DP) :: tmp
#ENDIF
  ! local variable
  REAL(KIND=DP) :: egap_dft
  ! parallelization
  INTEGER       :: ik_star, ik_stop
  !
  !
  CALL para_bounds (ik_star,ik_stop,nk_irr)
  !
  ! check nbndskip and spin
  IF (nbndskip .NE. 0) THEN
     IF (.NOT. lspinorb) THEN
        nelec_red = NINT(nelec) - 2*nbndskip
     ELSE
        nelec_red = NINT(nelec) - nbndskip
     ENDIF
  ELSE
     nelec_red = nelec
  ENDIF
  !

  !IF (nspin .EQ. 1 .OR. nspin .EQ. 4) THEN
  IF (.NOT. lspinorb) THEN
     vbnd_num = NINT(dble(nelec_red))/2
  ELSE
     vbnd_num = nelec_red
  ENDIF
  !
  cbnd_num  = nbndsub - vbnd_num
  cbnd_emin = +1.0d+8
  vbnd_emax = -1.0d+8
  !
  DO ik = ik_star, ik_stop
     DO ibnd = 1, vbnd_num
        IF (etf(ibnd,ik) .GT. vbnd_emax) vbnd_emax = etf(ibnd,ik)
     ENDDO
  ENDDO
  !
  DO ik = ik_star, ik_stop
     DO ibnd = vbnd_num+1, nbndsub
        IF (etf(ibnd,ik) .LT. cbnd_emin) cbnd_emin = etf(ibnd,ik)
     ENDDO
  ENDDO
  !
  !
#ifdef __PARA
  CALL mp_barrier (inter_pool_comm)
  CALL mp_min(cbnd_emin,inter_pool_comm)
  CALL mp_max(vbnd_emax,inter_pool_comm)
#endif
  !
  egap_dft = cbnd_emin - vbnd_emax
  ef_m = (vbnd_emax+cbnd_emin)/2.0d0
  ef = ef_m
  !
  WRITE(stdout,'(/5x,a,f9.4,a)')       'DFT/GW bandgap     = ', egap_dft*ryd2ev, ' eV' 
  WRITE(stdout,'(5x,a,f9.4,a,f9.4,a)') 'Bandgap            = ', vbnd_emax*ryd2ev, ' eV to ', cbnd_emin*ryd2ev, ' eV'
  WRITE(stdout,'(5x,a,f9.4,a)')        'Fermi level        = ', ef*ryd2ev, ' eV (center of gap)'
  !
  !
  IF (egap_rbm .NE. 0.0d0) THEN
     delta_egap = (egap_rbm/ryd2ev) - egap_dft ! plus=shift upward; mins=shift downward of conduction bands
  ELSE
     delta_egap = 0.0d0
  ENDIF
  !
  cbnd_emin = cbnd_emin + delta_egap
  !
  CALL mp_barrier (inter_pool_comm)
  !
  !
  ! shift the electron energy
  ! do not loop from ik_star to ik_stop due to only part of etf will be shifted
  !DO ik = 1, nk_irr
  !   DO ibnd = 1, nbndsub
  !      IF (etf(ibnd,ik) .GT. ef_m) etf(ibnd,ik) = etf(ibnd,ik) + delta_egap
  !   ENDDO
  !ENDDO
  DO ik = 1, nk_irr
     DO ibnd = vbnd_num+1, nbndsub
        etf(ibnd,ik) = etf(ibnd,ik) + delta_egap
     ENDDO
  ENDDO
  !
  !ef_m not changed
  ef = ef_m + 0.5d0*delta_egap
  !
  CALL mp_barrier (inter_pool_comm)
  !
  !
  IF (egap_rbm .NE. 0.0d0) THEN
     WRITE(stdout,'(/5x,a,f9.4,a)')    'RBM bandgap        = ', egap_rbm, ' eV'  
  ELSE
     WRITE(stdout,'(/5x,a)')           'RBM bandgap        =    N/A'
  ENDIF
  WRITE(stdout,'(5x,a,f9.4,a,f9.4,a)') 'New bandgap        = ', vbnd_emax*ryd2ev, ' eV to ', cbnd_emin*ryd2ev, ' eV'
  WRITE(stdout,'(5x,a,f9.4,a)')        'New Fermi level    = ', ef*ryd2ev, ' eV (center of gap)' 
  !
  !
  ibndmin = 100000
  ibndmax = 0
  ebndmin =  1.d8
  ebndmax = -1.d8
  !
  cfsthick = fsthick
  vfsthick = fsthick
  !
  DO ik = ik_star, ik_stop
     DO ibnd = 1, nbndsub
        !
        IF (nptype .EQ. 'n') THEN
           if ((etf(ibnd,ik)>=cbnd_emin).and.(etf(ibnd,ik)<=cbnd_emin+cfsthick)) then
              ibndmax = max(ibnd,ibndmax)
              ibndmin = min(ibnd,ibndmin)
              ebndmax = max(etf(ibnd,ik),ebndmax)
              ebndmin = min(etf(ibnd,ik),ebndmin)
           endif
        ELSEIF (nptype .EQ. 'p') THEN 
           if ((etf(ibnd,ik)<=vbnd_emax).and.(etf(ibnd,ik)>=vbnd_emax-vfsthick)) then
              ibndmax = max(ibnd,ibndmax)
              ibndmin = min(ibnd,ibndmin)
              ebndmax = max(etf(ibnd,ik),ebndmax)
              ebndmin = min(etf(ibnd,ik),ebndmin)
           endif
        ELSE
           IF (ABS(etf(ibnd,ik)-cbnd_emin) .LE. cfsthick) THEN
              ibndmax = max(ibnd,ibndmax)
              ebndmax = max(etf(ibnd,ik),ebndmax)
           ENDIF
           !
           IF (ABS(etf(ibnd,ik)-vbnd_emax) .LE. vfsthick) THEN
              ibndmin = min(ibnd,ibndmin)
              ebndmin = min(etf(ibnd,ik),ebndmin)
           ENDIF
        ENDIF
        !
     ENDDO
  ENDDO
  !
#ifdef __PARA
  CALL mp_barrier (inter_pool_comm)
  !
  tmp = dble (ibndmin)
  CALL mp_min(tmp,inter_pool_comm)
  ibndmin = nint (tmp)
  CALL mp_min(ebndmin,inter_pool_comm)
  !
  tmp = dble (ibndmax)
  CALL mp_max(tmp, inter_pool_comm)
  ibndmax = nint (tmp)
  CALL mp_max(ebndmax,inter_pool_comm)
#endif 
  !
  IF (bte .EQ. 3) THEN
     !
     IF (nptype .EQ. 'n') THEN
        ibndmin = vbnd_num+1
     ELSEIF (nptype .EQ. 'p') THEN
        ibndmax = vbnd_num
     ENDIF
     !
  ENDIF
  !
  nbnd_red = ibndmax - ibndmin + 1
  !
  !
  WRITE(stdout,'(/5x,a,i4)')           'Number of skipped bands = ', nbndskip
  WRITE(stdout,'(/5x,a,f9.4)')         'Number of electron = ', nelec
  WRITE(stdout,'(5x,a,i4)')            'Number of EPW electrons = ', nelec_red
  WRITE(stdout,'(5x,a,i4,a,i4)')       'Number of v/c band = ', vbnd_num, '          / ', cbnd_num
  WRITE(stdout,'(5x,a,i4,a,i4)')       'Band window        = ', ibndmin, '         to ', ibndmax
  WRITE(stdout,'(5x,a,f9.4,a,f9.4,a)') 'v/c fsthick        = ', vfsthick*ryd2ev, '     / ', cfsthick*ryd2ev,  ' eV'
  WRITE(stdout,'(5x,a,f9.4,a,f9.4,a)') 'Energy window      = ', ebndmin*ryd2ev, ' eV to ', ebndmax*ryd2ev, ' eV'
  !
  END SUBROUTINE fermiwindow
  !
  !
  !
!-----------------------------------------------------------------------
SUBROUTINE fermiwindow_reassign (ef0, ibndmin, ibndmax, ebndmin, ebndmax, vfsthick, cfsthick)
!-----------------------------------------------------------------------
#INCLUDE "f_defs.h"
  !
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE spin_orb,      ONLY : lspinorb
  USE lsda_mod,      ONLY : nspin
  USE elph2,         ONLY : etf, cbnd_emin, vbnd_emax, vbnd_num, cbnd_num, &
                            cfs => cfsthick, vfs => vfsthick
  USE epwcom,        ONLY : fsthick, nbndsub, bte, nptype
  USE pwcom,         ONLY : ef, nelec
  USE constants_epw, ONLY : ryd2ev
  USE bte_var,       ONLY : nk_irr
#ifdef __PARA
  USE mp,            ONLY : mp_max, mp_min, mp_barrier
  USE mp_global,     ONLY : inter_pool_comm
#ENDIF
  !
  IMPLICIT NONE
  !
  ! local variable
  REAL(KIND=DP), INTENT(IN)  :: ef0
  REAL(KIND=DP), INTENT(OUT) :: ebndmin, ebndmax, cfsthick, vfsthick
  INTEGER, INTENT(OUT)       :: ibndmin, ibndmax
  INTEGER                    :: ik, ibnd
  ! parallelization
  INTEGER                    :: ik_star, ik_stop
  !
#ifdef __PARA
  REAL(KIND=DP)              :: tmp
#ENDIF
  !
  !
  CALL para_bounds (ik_star,ik_stop,nk_irr)
  !
  ibndmin = 100000
  ibndmax = 0
  ebndmin =  1.d8
  ebndmax = -1.d8
  !
  IF (ef0 .GT. cbnd_emin) THEN
     cfsthick = cfs + (ef0-cbnd_emin)
     vfsthick = vfs
  ELSEIF (ef0 .LT. vbnd_emax) THEN
     cfsthick = cfs
     vfsthick = vfs + (vbnd_emax-ef0)
  ENDIF
  !
  DO ik = ik_star, ik_stop
     DO ibnd = 1, nbndsub
        !
        IF (nptype .EQ. 'n') THEN
           if ((etf(ibnd,ik)>=cbnd_emin).and.(etf(ibnd,ik)<=cbnd_emin+cfsthick)) then
              ibndmax = max(ibnd,ibndmax)
              ibndmin = min(ibnd,ibndmin)
              ebndmax = max(etf(ibnd,ik),ebndmax)
              ebndmin = min(etf(ibnd,ik),ebndmin)
           endif
        ELSEIF (nptype .EQ. 'p') THEN 
           if ((etf(ibnd,ik)<=vbnd_emax).and.(etf(ibnd,ik)>=vbnd_emax-vfsthick)) then
              ibndmax = max(ibnd,ibndmax)
              ibndmin = min(ibnd,ibndmin)
              ebndmax = max(etf(ibnd,ik),ebndmax)
              ebndmin = min(etf(ibnd,ik),ebndmin)
           endif
        ELSE
           IF (ABS(etf(ibnd,ik)-cbnd_emin) .LE. cfsthick) THEN
              ibndmax = max(ibnd,ibndmax)
              ebndmax = max(etf(ibnd,ik),ebndmax)
           ENDIF
           !
           IF (ABS(etf(ibnd,ik)-vbnd_emax) .LE. vfsthick) THEN
              ibndmin = min(ibnd,ibndmin)
              ebndmin = min(etf(ibnd,ik),ebndmin)
           ENDIF
        ENDIF
        !
     ENDDO
  ENDDO
  !
#ifdef __PARA
  CALL mp_barrier (inter_pool_comm)
  !
  tmp = dble (ibndmin)
  CALL mp_min(tmp,inter_pool_comm)
  ibndmin = nint (tmp)
  CALL mp_min(ebndmin,inter_pool_comm)
  !
  tmp = dble (ibndmax)
  CALL mp_max(tmp, inter_pool_comm)
  ibndmax = nint (tmp)
  CALL mp_max(ebndmax,inter_pool_comm)
#endif 
  !
  IF (bte .EQ. 3) THEN
     !
     IF (nptype .EQ. 'n') THEN
        ibndmin = vbnd_num+1
     ELSEIF (nptype .EQ. 'p') THEN
        ibndmax = vbnd_num
     ENDIF
     !
  ENDIF
  !
  CALL mp_barrier (inter_pool_comm)
  !
END SUBROUTINE fermiwindow_reassign
