!-----------------------------------------------------------------------
SUBROUTINE fermilocatioepdope ()
!-----------------------------------------------------------------------
!
!  find the fermi-level location with specified doping-level
! 
!-----------------------------------------------------------------------
  !
#INCLUDE "f_defs.h"
  !
  USE kinds,         ONLY : DP
  USE pwcom,         ONLY : ef
  USE lsda_mod,      ONLY : nspin
  USE epwcom,        ONLY : eptemp, epdope, neptemp, nepdope, fscheck, L_D, eimp_mode, &
                            lpolar, dielec, screen_polar
  USE elph2,         ONLY : n_elec, n_hole, n_intr, cbnd_emin, vbnd_emax, outside_gap, &
                            ndos, edos, ef_epw, ibndmin, ibndmax, ebndmin, ebndmax, cfsthick, vfsthick, nbnd_red, &
                            epsi
  USE constants_epw, ONLY : ryd2ev, bohr2ang, kB, pi
#ifdef __PARA
  USE mp,            ONLY : mp_barrier
  USE mp_global,     ONLY : inter_pool_comm, my_pool_id
  USE io_files,      ONLY : tmp_dir, prefix
  USE io_global,     ONLY : stdout
#ENDIF
  !
  IMPLICIT NONE
  !
  REAL(KIND=DP), PARAMETER   :: ef_extent = 2.0d0/ryd2ev, eps = 1.0d-6
  REAL(KIND=DP)              :: E_cen, ef_min, ef_max, ef_pl, ef_pr, &
                                ef_tmp, ef_new, temp, doping, non_r, n_plus, n_mins, &
                                LD2inv, dielec0, deltaE, ef0, E, d_wgk, &
                                dos_sp1, dos_sp2, dos_tot
  REAL(KIND=DP), EXTERNAL    :: char_dens
  INTEGER                    :: itemp, idope, i
  ! reassign
  INTEGER                    :: ibnd_0(neptemp,nepdope), ibnd_1(neptemp,nepdope)
  REAL(KIND=DP)              :: ebnd_0(neptemp,nepdope), ebnd_1(neptemp,nepdope), cf_0(neptemp,nepdope), vf_0(neptemp,nepdope)

  REAL(kind=DP), external      :: de_wgauss
  !
  ! ERROR: assign epsi to Silicon's values if [lpolar] = .false.
  if (.not. lpolar) then
     epsi(1,1) = dielec
     epsi(2,2) = dielec
     epsi(3,3) = dielec
  endif
  dielec0 = (epsi(1,1)+epsi(2,2)+epsi(3,3)) / 3.d0
  !
  ! read eDOS from file
  CALL edos_load ()
  !
  E_cen = 0.5d0 * (cbnd_emin + vbnd_emax)
  ef_min = (E_cen - (0.5d0 * (cbnd_emin - vbnd_emax) + ef_extent)) * ryd2ev ! initial lower bound in searching fermi level
  ef_max = (E_cen + (0.5d0 * (cbnd_emin - vbnd_emax) + ef_extent)) * ryd2ev ! initial upper bound in searching fermi level
  !
  ALLOCATE (ef_epw(neptemp,nepdope))
  ALLOCATE (n_elec(neptemp,nepdope))
  ALLOCATE (n_hole(neptemp,nepdope))
  ALLOCATE (n_intr(neptemp,nepdope))
  ef_epw = 0.0d0
  n_elec  = 0.0d0
  n_hole  = 0.0d0
  n_intr  = 0.0d0
  !
  ! compute the fermi energy
  DO itemp = 1, neptemp
     !
     temp = eptemp(itemp)*ryd2ev ! in atomic units
     ! 
     DO idope = 1, nepdope
        !
        doping = epdope(idope)
        !
        ef_pl = ef_min
        ef_pr = ef_max
        !
        n_plus = char_dens (E_cen*ryd2ev, E_cen*ryd2ev, temp, 0.0d0, +1)
        n_mins = char_dens (E_cen*ryd2ev, E_cen*ryd2ev, temp, 0.0d0, -1)
        n_intr(itemp,idope) = SQRT((n_plus+n_mins)*(n_plus-n_mins)/4.0d0)
        !
        !
        IF ( (char_dens (ef_min, E_cen*ryd2ev, temp, doping, 0) .LT. 0.0d0) .OR. &
             (char_dens (ef_max, E_cen*ryd2ev, temp, doping, 0) .GT. 0.0d0) ) THEN
           !
           CALL errore('fermilocation', "fermi level starting range not including the specified doping level",1)
           !
        ENDIF
        !
        DO WHILE (ABS(ef_pr-ef_pl) .GE. eps) 
           ef_tmp = (ef_pr + ef_pl) / 2.0d0
           IF (char_dens (ef_tmp, E_cen*ryd2ev, temp, doping, 0) .LT. 0.0d0) THEN
              ef_pr = ef_tmp
           ELSE
              ef_pl = ef_tmp
           ENDIF
        ENDDO
        !          
        ef_new = (ef_pl + ef_pr) / 2.0d0  
        ef_epw(itemp,idope) = ef_new/ryd2ev ! change to [Ry]
        !
        n_plus = char_dens (ef_new, E_cen*ryd2ev, temp, doping, +1)
        n_mins = char_dens (ef_new, E_cen*ryd2ev, temp, doping, -1)
        n_hole(itemp,idope) = (n_plus + n_mins) / 2.0d0
        n_elec(itemp,idope) = (n_plus - n_mins) / 2.0d0
        !
        !
        WRITE (stdout,'(/5x,a,i2,a,i2,a)')  '[T', itemp, ' | N', idope, ']'
        WRITE (stdout,'(5x,a,f8.1,a)')      'Temperature     = ', eptemp(itemp)/kB, ' K'
        WRITE (stdout,'(5x,a,es15.6,a)')    'Doping          = ', epdope(idope), ' cm^-3'
        WRITE (stdout,'(5x,a,f9.4,a)')      'New Fermi level = ', ef_epw(itemp,idope)*ryd2ev
        WRITE (stdout,'(5x,a,es15.6,a)')    'n_hole          = ', n_hole(itemp,idope), ' cm^-3'
        WRITE (stdout,'(5x,a,es15.6,a)')    'n_elec          = ', n_elec(itemp,idope), ' cm^-3'
        WRITE (stdout,'(5x,a,es15.6,a)')    'n_intr          = ', n_intr(itemp,idope), ' cm^-3'
        !
        CALL mp_barrier (inter_pool_comm)
        !
     ENDDO
     !
  ENDDO
  !
  !
  CALL mp_barrier (inter_pool_comm)
  !
  ! if eimp_mode > 0 or screen_polar = true, calculate the Debye screening length here
  L_D = 0.d0

  if (eimp_mode > 0 .or. screen_polar) then
     ! dos(1,n) in eV, dos(2/3,n) in 1/eV/cm^3
     deltaE = 9.876543210d-6 / ryd2ev
     !
     if (dielec0 == 0) &
        call errore ('fermilocation', 'epsi is zero, cannot perform Debye length calculation', dielec0)
     !
     do itemp = 1, neptemp
        temp = eptemp(itemp) ! in atomic units [Ryd]
        !
        do idope = 1, nepdope
           doping = epdope(idope)
           ef0 = ef_epw(itemp,idope)
           !
           LD2inv = 0
           i = 1
           E = edos(1,1)/ryd2ev + deltaE
           DO WHILE (E .LT. edos(1,ndos)/ryd2ev)
              ! get dos_tot(E) by using interpolation
              ! dos_sp2!=0 if spin up and spin down are distinguishable
              dos_sp1 = edos(2,i) + ( (edos(2,i+1)-edos(2,i)) / (edos(1,i+1)-edos(1,i)) ) * (E*ryd2ev-edos(1,i))
              dos_sp2 = edos(3,i) + ( (edos(3,i+1)-edos(3,i)) / (edos(1,i+1)-edos(1,i)) ) * (E*ryd2ev-edos(1,i))
              ! change to atomic unit, 1/Ry/aB^3
              dos_tot = (dos_sp1 + dos_sp2) * ryd2ev * ((bohr2ang*1.0d-8)**3.d0)
              !
              if (doping < 0) then
                 IF (E .GE. E_cen) THEN
                    d_wgk = de_wgauss( -(E-ef0)/temp, temp)
                    LD2inv = LD2inv + dos_tot * d_wgk * deltaE
                 ENDIF
              else
                 IF (E .LE. E_cen) THEN
                    d_wgk = de_wgauss( (E-ef0)/temp, temp)
                    LD2inv = LD2inv + dos_tot * d_wgk * deltaE
                 ENDIF
              endif
              !
              E = E + deltaE
              !
              IF (E .GT. edos(1,i+1)/ryd2ev) i = i + 1
              !
           ENDDO
           !
           LD2inv = LD2inv * 2.d0 * (4.d0 * pi) / dielec0
           if (LD2inv == 0) &
              call errore ('fermilocation', 'Debye length turns out to be zero..', LD2inv)

           L_D(itemp,idope) = sqrt(1.d0 / LD2inv)    ! in atomic unit [Bohr radius]

           write(stdout,*) ' (temp, doping):', eptemp(itemp)/6.33362d-6, doping
           write(stdout,*) ' --------- L_D =', L_D(itemp,idope)*0.529177d0, 'Ang'
           !
        enddo
     enddo
  endif
  !
  ! recompute the fermiwindow
  ibnd_0 = 0
  ibnd_1 = 0
  ebnd_0 = 0.0d0
  ebnd_1 = 0.0d0
  cf_0   = 0.0d0
  vf_0   = 0.0d0
  outside_gap = .FALSE.
  !
  !
  IF (fscheck .EQ. .TRUE.) THEN
     !
     WRITE (stdout,'(/5x,a/)') 'Compute Fermiwindow for all temperatures and concentrations : ' 
     WRITE (stdout,'(9x,a)')   '|  temp  |     dope     | ibndmin | ibndmax | ebndmin | ebndmax | vfsthick | cfsthick | ef_loc |' 
     WRITE (stdout,'(9x,a)')   '------------------------------------------------------------------------------------------------' 
     !
     DO itemp = 1, neptemp
        !
        DO idope = 1, nepdope
           !
           IF (ef_epw(itemp,idope) .GT. cbnd_emin) THEN
              !
              outside_gap = .TRUE.
              CALL fermiwindow_reassign (ef_epw(itemp,idope), ibnd_0(itemp,idope), ibnd_1(itemp,idope), ebnd_0(itemp,idope), ebnd_1(itemp,idope), &
                                         vf_0(itemp,idope), cf_0(itemp,idope))
              !
              WRITE (stdout,'(5x,f11.1,es16.6,i10,i10,f10.4,f10.4,f11.4,f11.4,a)') &
                     eptemp(itemp)/kB, epdope(idope), ibnd_0(itemp,idope), ibnd_1(itemp,idope), ebnd_0(itemp,idope)*ryd2ev, ebnd_1(itemp,idope)*ryd2ev, &
                     vf_0(itemp,idope)*ryd2ev, cf_0(itemp,idope)*ryd2ev, '       CB'
              !
           ELSEIF (ef_epw(itemp,idope) .LT. vbnd_emax) THEN
              !
              outside_gap = .TRUE.
              CALL fermiwindow_reassign (ef_epw(itemp,idope), ibnd_0(itemp,idope), ibnd_1(itemp,idope), ebnd_0(itemp,idope), ebnd_1(itemp,idope), &
                                         vf_0(itemp,idope), cf_0(itemp,idope))
              !
              WRITE (stdout,'(5x,f11.1,es16.6,i10,i10,f10.4,f10.4,f11.4,f11.4,a)') &
                     eptemp(itemp)/kB, epdope(idope), ibnd_0(itemp,idope), ibnd_1(itemp,idope), ebnd_0(itemp,idope)*ryd2ev, ebnd_1(itemp,idope)*ryd2ev, &
                     vf_0(itemp,idope)*ryd2ev, cf_0(itemp,idope)*ryd2ev, '       VB'
              !
           ELSE
              !
              ibnd_0(itemp,idope) = ibndmin
              ibnd_1(itemp,idope) = ibndmax
              ebnd_0(itemp,idope) = ebndmin
              ebnd_1(itemp,idope) = ebndmax
              vf_0(itemp,idope)   = vfsthick
              cf_0(itemp,idope)   = cfsthick
              !
              WRITE (stdout,'(5x,f11.1,es16.6,i10,i10,f10.4,f10.4,f11.4,f11.4,a)') &
                     eptemp(itemp)/kB, epdope(idope), ibnd_0(itemp,idope), ibnd_1(itemp,idope), ebnd_0(itemp,idope)*ryd2ev, ebnd_1(itemp,idope)*ryd2ev, &
                     vf_0(itemp,idope)*ryd2ev, cf_0(itemp,idope)*ryd2ev, '      GAP'
              !
           ENDIF
           !
        ENDDO
        !
     ENDDO
     !
     ibndmin  = MINVAL(ibnd_0)
     ibndmax  = MAXVAL(ibnd_1)
     ebndmin  = MINVAL(ebnd_0)
     ebndmax  = MAXVAL(ebnd_1)
     vfsthick = MAXVAL(vf_0)
     cfsthick = MAXVAL(cf_0)
     !
  ELSE
     !
     WRITE (stdout,'(/5x,a/)') 'Do not compute Fermiwindow for all temperatures and concentrations' 
     !
  ENDIF
  !
  nbnd_red = ibndmax - ibndmin + 1
  !
  WRITE(stdout,'(/5x,a,i4,a,i4)')      'New band window    = ', ibndmin, '         to ', ibndmax
  WRITE(stdout,'(5x,a,f9.4,a,f9.4,a)') 'New v/c fsthick    = ', vfsthick*ryd2ev, '     / ', cfsthick*ryd2ev,  ' eV'
  WRITE(stdout,'(5x,a,f9.4,a,f9.4,a)') 'New energy window  = ', ebndmin*ryd2ev, ' eV to ', ebndmax*ryd2ev, ' eV'
  !
  CALL mp_barrier (inter_pool_comm)
  DEALLOCATE (edos)
  !
END SUBROUTINE fermilocatioepdope



!-----------------------------------------------------------------------
SUBROUTINE fermilocation_ef ()
!-----------------------------------------------------------------------
!
!  find the doping-level with specified fermi-level location 
! 
!-----------------------------------------------------------------------
  !
#INCLUDE "f_defs.h"
  !
  USE kinds,         ONLY : DP
  USE pwcom,         ONLY : ef
  USE lsda_mod,      ONLY : nspin
  USE epwcom,        ONLY : eptemp, neptemp, fermi_energy
  USE elph2,         ONLY : n_elec, n_hole, n_intr, cbnd_emin, vbnd_emax, outside_gap, dope_ef, &
                            ndos, edos, ef_epw, ibndmin, ibndmax, ebndmin, ebndmax, cfsthick, vfsthick, nbnd_red
  USE constants_epw, ONLY : ryd2ev, bohr2ang, kB
#ifdef __PARA
  USE mp,            ONLY : mp_barrier
  USE mp_global,     ONLY : inter_pool_comm
  USE io_files,      ONLY : tmp_dir, prefix
  USE io_global,     ONLY : stdout
#ENDIF
  !
  IMPLICIT NONE
  !
  REAL(KIND=DP), PARAMETER   :: ef_extent = 1.0d-6/ryd2ev
  REAL(KIND=DP)              :: E_cen, ef_min, ef_max, &
                                ef_tmp, ef_new, temp, non_r, n_plus, n_mins, dope(neptemp)
  REAL(KIND=DP), EXTERNAL    :: char_dens
  INTEGER                    :: itemp
  ! reassign
  INTEGER                    :: ibnd_0, ibnd_1
  REAL(KIND=DP)              :: ebnd_0, ebnd_1, cf_0, vf_0
  !
  !
  ! read eDOS from file
!  CALL edos_load ()
!  !
!  ALLOCATE (ef_epw(neptemp))
!  ALLOCATE (dope_ef(neptemp))
!  ALLOCATE (n_elec(neptemp))
!  ALLOCATE (n_hole(neptemp))
!  ALLOCATE (n_intr(neptemp))
!  ef_epw = 0.0d0
!  dope_ef = 0.0d0
!  n_elec  = 0.0d0
!  n_hole  = 0.0d0
!  n_intr  = 0.0d0
!  !
!  ! compute the fermi energy
!  WRITE (stdout,'(/5x,a/)') 'Check the concentration of all temperatures : ' 
!  WRITE (stdout,'(9x,a)')   '   temp  |     Dope     |     Hole     |     Elec     |     Intr   '
!  WRITE (stdout,'(9x,a)')   '-------------------------------------------------------------------' 
!  DO itemp = 1, neptemp
!     !
!     temp = eptemp(itemp)*ryd2ev ! in atomic units
!     ! 
!     E_cen = fermi_energy
!     !
!     n_plus = char_dens (E_cen*ryd2ev, E_cen*ryd2ev, temp, 0.0d0, +1)
!     n_mins = char_dens (E_cen*ryd2ev, E_cen*ryd2ev, temp, 0.0d0, -1)
!     !
!     n_intr(itemp) = SQRT((n_plus+n_mins)*(n_plus-n_mins)/4.0d0)
!     n_hole(itemp) = (n_plus + n_mins) / 2.0d0
!     n_elec(itemp) = (n_plus - n_mins) / 2.0d0
!     dope_ef(itemp) = n_hole(itemp) - n_elec(itemp)
!     !
!     WRITE (stdout,'(9x,f8.1,4es15.6)') eptemp(itemp)/kB, dope_ef(itemp), n_hole(itemp), n_elec(itemp), n_intr(itemp)
!     !
!  ENDDO
!  !
!  !
!  ! recompute the fermiwindow
!  ibnd_0 = 0
!  ibnd_1 = 0
!  ebnd_0 = 0.0d0
!  ebnd_1 = 0.0d0
!  cf_0   = 0.0d0
!  vf_0   = 0.0d0
!  outside_gap = .FALSE.
!  !
!  IF (fermi_energy .GT. cbnd_emin) THEN
!     !
!     WRITE (stdout,'(/5x,a)') 'Fermi level is at conduction band'
!     outside_gap = .TRUE.
!     CALL fermiwindow_reassign (fermi_energy, ibnd_0, ibnd_1, ebnd_0, ebnd_1, vf_0, cf_0)
!     !
!     ibndmin  = ibnd_0
!     ibndmax  = ibnd_1
!     ebndmin  = ebnd_0
!     ebndmax  = ebnd_1
!     vfsthick = vf_0
!     cfsthick = cf_0
!     nbnd_red = ibndmax - ibndmin + 1
!     !
!  ELSEIF (fermi_energy .LT. vbnd_emax) THEN
!     !
!     WRITE (stdout,'(/5x,a)') 'Fermi level is at valence band'
!     outside_gap = .TRUE.
!     CALL fermiwindow_reassign (fermi_energy, ibnd_0, ibnd_1, ebnd_0, ebnd_1, vf_0, cf_0)
!     !
!     ibndmin  = ibnd_0
!     ibndmax  = ibnd_1
!     ebndmin  = ebnd_0
!     ebndmax  = ebnd_1
!     vfsthick = vf_0
!     cfsthick = cf_0
!     nbnd_red = ibndmax - ibndmin + 1
!     !
!  ELSE
!     !
!     WRITE (stdout,'(/5x,a)') 'Fermi level is inside bandgap'
!     !
!  ENDIF
!  !
!  WRITE(stdout,'(/5x,a,f9.4,a)')       'Set Fermi level    = ', fermi_energy*ryd2ev, ' eV'
!  WRITE(stdout,'(5x,a,i4,a,i4)')       'New band window    = ', ibndmin, '         to ', ibndmax
!  WRITE(stdout,'(5x,a,f9.4,a,f9.4,a)') 'New v/c fsthick    = ', vfsthick*ryd2ev, '     / ', cfsthick*ryd2ev,  ' eV'
!  WRITE(stdout,'(5x,a,f9.4,a,f9.4,a)') 'New energy window  = ', ebndmin*ryd2ev, ' eV to ', ebndmax*ryd2ev, ' eV'
!  !
!  CALL mp_barrier (inter_pool_comm)
!  DEALLOCATE (edos)
  !
END SUBROUTINE fermilocation_ef



!-----------------------------------------------------------------------
FUNCTION char_dens (ef0, e_cen, temp, ndoping, plus)
!-----------------------------------------------------------------------
!
!  calculate the carrier density or charge density based on the
!     input DOS data and fermi-level, can also consider dopant level
!     plus =  0 : net charge density = n_hole - n_elec - NA + ND
!     plus = -1 : charge density  = nhole - nelec
!     plus =  1 : carrier density = nhole + nelec
!
!-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE elph2,      ONLY : ndos, edos
  USE io_global,  ONLY : stdout
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)       :: plus
  REAL(KIND=DP), INTENT(IN) :: ef0, e_cen, ndoping, temp
  ! local variable
  REAL(KIND=DP), PARAMETER  :: deltaE = 9.876543210d-6
  INTEGER                   :: i
  REAL(KIND=DP)             :: E, edos_tot, edos_sp1, edos_sp2, wgk, nelec, nhole
  REAL(KIND=DP)             :: char_dens, wgauss
  !
  !
  nelec = 0.d0
  nhole = 0.d0
  !
  !
  ! count electrons and holes
  !
  i = 1
  E = edos(1,1) + deltaE
  DO WHILE (E .LT. edos(1,ndos))
     ! get edos_tot(E) by using interpolation
     ! edos_sp2!=0 if spin up and spin down are distinguishable
     edos_sp1 = edos(2,i) + ( (edos(2,i+1)-edos(2,i)) / (edos(1,i+1)-edos(1,i)) ) * (E-edos(1,i))
     edos_sp2 = edos(3,i) + ( (edos(3,i+1)-edos(3,i)) / (edos(1,i+1)-edos(1,i)) ) * (E-edos(1,i))
     edos_tot = edos_sp1 + edos_sp2
     !
     IF (E .GE. e_cen) THEN
        wgk = wgauss( -(E-ef0)/temp, -99)
        nelec = nelec + edos_tot * wgk * deltaE
     ELSE
        wgk = wgauss( (E-ef0)/temp, -99)
        nhole = nhole + edos_tot * wgk * deltaE
     ENDIF
     !
     E = E + deltaE
     !
     IF (E .GT. edos(1,i+1)) i = i + 1
     !
  ENDDO
  !
  !
  ! here different dopant models can be applied (X)
  IF (plus .EQ. 0) THEN
     char_dens = (nhole - nelec) - ndoping
  ELSE
     char_dens = nhole + plus * nelec
  ENDIF
  !
  !
END FUNCTION char_dens


!-----------------------------------------------------------------------
  FUNCTION de_wgauss (x, temp)
  !-----------------------------------------------------------------------
  !
  !     this function computes the derivative of the fermi-dirac distribution
  !     at the point x. essentially it's calculating:
  !
  !     -partial(f)/partial(E) = (1/kB.T)*exp((E-Ef)/(kB.T))/(1+exp(..))^2
  !
  USE kinds, ONLY : DP
  use io_global, only : stdout
  !
  implicit none
  real(DP) :: de_wgauss, x, temp
  ! output: the value of the function
  ! input: the argument of the function
  !
  !    the local variables
  !
  real(DP), parameter :: ryd2ev = 13.6058d0
  real(DP), parameter :: maxarg = 200.d0
  real(DP) :: tempk
  ! maximum value for the argument of the exponential

  ! change unit of temp [Ry] -> [K]
  !tempk = (temp * ryd2ev) / 0.000086173423d0

  IF (x.gt.maxarg) THEN
     de_wgauss = 0.d0
  ELSEIF (x.lt.-maxarg) THEN
     de_wgauss = 0.d0
  ELSE
     de_wgauss = (exp(-x) / temp) / ((1.0d0 + exp(-x)) ** 2)    ! [1/Ry]
  ENDIF

  RETURN
  !
  END FUNCTION de_wgauss

