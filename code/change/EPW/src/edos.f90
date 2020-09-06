#INCLUDE "f_defs.h"
!
!--------------------------------------------------------------------------
MODULE edos
!-------------------------------------------------------------------------- 
  USE kinds, ONLY :  DP
  !
  SAVE
  !
  REAL(KIND=DP), ALLOCATABLE ::  &
          xkf_ful_dos(:,:)      ,&
          xkf_irr_dos(:,:)      ,&
          wkf_irr_dos(:)        ,&
          symmat_lat_dos(:,:,:)
  !
  INTEGER, ALLOCATABLE       ::  &
          ful2irr_dos(:)        ,&
          equiv_dos(:)
  !
  INTEGER                    ::  &
          nk_ful_dos            ,&
          nk_irr_dos
  !                 
END MODULE edos
!
!
!--------------------------------------------------------------------
SUBROUTINE edos_energy ()
  !--------------------------------------------------------------------
  !
  ! Calculates the Density of States (DOS),
  ! to do: separated into up and down components for LSDA
  !
  ! use tetrahedron method, adapted from dos.f90 and dos_t.f90
  !
  !
  USE kinds,         ONLY : DP
  USE io_files,      ONLY : tmp_dir, prefix
  USE io_global,     ONLY : stdout
  USE spin_orb,      ONLY : lspinorb
  USE lsda_mod,      ONLY : nspin
  USE pwcom,         ONLY : omega
  USE constants_epw, ONLY : ryd2ev, bohr2ang
  USE epwcom,        ONLY : nbndsub, nkf1, nkf2, nkf3, nkfdos1, nkfdos2, nkfdos3
  USE elph2,         ONLY : vbnd_emax, cbnd_emin, ef_m, delta_egap, nrr_k, irvec, &
                            vbnd_num, ndegen_k, chw
  USE tetrahedron,   ONLY : ntetra_dose, tetra_dose
  USE edos
#ifdef __PARA
  USE io_global,     ONLY : ionode_id, stdout
  USE mp,            ONLY : mp_barrier, mp_sum
  USE mp_global,     ONLY : my_pool_id, inter_pool_comm
#ENDIF
  !
  IMPLICIT NONE
  !
  REAL(KIND=DP), PARAMETER      :: E_extent = 3.5d0/ryd2ev, deltaE = 1.0d-4/ryd2ev
  !
  REAL(KIND=DP)                 :: DOSofE(2), DOSint, E_min, E_max, E_cen
  INTEGER                       :: n, ik, ik0, ibnd
  ! output edos file
  CHARACTER(LEN=256)            :: fildos
  LOGICAL                       :: minus_q, od
  ! parallelization
  INTEGER                       :: ndos_ful, ndos_pol, ndos_star, ndos_stop
  REAL(KIND=DP), ALLOCATABLE    :: dos_ful(:,:)
  ! computational time
  REAL(KIND=DP)                 :: t0, t1
  ! band energy
  INTEGER                       :: nk_pol, ik_star, ik_stop
  REAL(KIND=DP)                 :: xxk(3)
  REAL(KIND=DP), ALLOCATABLE    :: et(:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: cufkk(:,:)
  !
  REAL(KIND=DP)                 :: t2i, t2f, t2=0.0d0 
  !
  !
  WRITE (stdout,'(/5x,a/)') 'Start computing eDOS: '
  !
  IF (nkfdos1 .EQ. 0) nkfdos1 = nkf1
  IF (nkfdos2 .EQ. 0) nkfdos2 = nkf2
  IF (nkfdos3 .EQ. 0) nkfdos3 = nkf3
  !
  ! get ful2all and xkf_irr_dos used in edos
  CALL edos_mesh (nkfdos1, nkfdos2, nkfdos3)
  !
  !
  ! calculate the band energy used in edos
  CALL para_bounds (ik_star, ik_stop, nk_irr_dos)
  !
  ALLOCATE (et(nbndsub, nk_irr_dos))
  ALLOCATE (cufkk(nbndsub,nbndsub))
  et     = 0.0d0
  cufkk  = (0.0d0,0.0d0)
  !
  DO ik = ik_star, ik_stop
     !
     xxk = xkf_irr_dos(:,ik)
     !
     CALL hamwan2bloch (nbndsub, nrr_k, irvec, ndegen_k, xxk, cufkk, et(:,ik), chw)
     !
     !DO ibnd = 1, nbndsub
     !   IF (et(ibnd,ik) .GT. ef_m) et(ibnd,ik) = et(ibnd,ik) + delta_egap
     !ENDDO
     DO ibnd = vbnd_num+1, nbndsub
        et(ibnd,ik) = et(ibnd,ik) + delta_egap
     ENDDO 
     ! 
  ENDDO
  !
#ifdef __PARA
  CALL mp_barrier (inter_pool_comm)
  CALL mp_sum (et,inter_pool_comm)
#ENDIF
  !
  !
  ! generate tetrahedra grid
  ntetra_dose = 6*nk_ful_dos
  ALLOCATE (tetra_dose(4,ntetra_dose))
  !
  CALL tetra_generate (nkfdos1, nkfdos2, nkfdos3, ntetra_dose, tetra_dose)
  !
  ! find band extrema (in Ryd)
  ! bandgap = cbnd_emin-vbnd_emax
  E_cen = 0.5d0 * (cbnd_emin + vbnd_emax)
  E_min  = E_cen - (0.5d0 * (cbnd_emin - vbnd_emax) + E_extent) 
  E_max  = E_cen + (0.5d0 * (cbnd_emin - vbnd_emax) + E_extent) 
  ndos_ful = NINT ( (E_max-E_min)/deltaE + 0.500001d0 )
  !
  ! density-of-state calculation
  !
  CALL CPU_TIME (t0)
  !
  CALL para_bounds (ndos_star, ndos_stop, ndos_ful)
  !
  ! dos(:,1) = energy; dos(:,1) = dos_up; dos(:,1) = dos_down; 
  ALLOCATE (dos_ful(3,ndos_ful))
  dos_ful = 0.0d0
  !
  DO n = ndos_star, ndos_stop
     !
     CALL CPU_TIME (t2i)
     !
     dos_ful(1,n) = E_min + DBLE(n-1)*deltaE
     CALL dos_tetra (et, nspin, nbndsub, nk_irr_dos, ntetra_dose, tetra_dose, dos_ful(1,n), DOSofE)
     !
     IF (nspin .EQ. 1 .OR. nspin .EQ. 4) THEN 
        dos_ful(2,n) = DOSofE(1)
     ELSE
        dos_ful(2,n) = DOSofE(1)
        dos_ful(3,n) = DOSofE(2)
     ENDIF  
     !
     CALL CPU_TIME (t2f)
     t2 = t2 + (t2f-t2i)
     IF (MOD(n-ndos_star+1,100) .EQ. 0 .OR. MOD(n-ndos_star+1,ndos_stop-ndos_star+1) .EQ. 0) THEN
        WRITE (stdout,'(9x,a,i5,a,i5,a,f7.1,a)') 'ndos = (', n-ndos_star+1, '/', ndos_stop-ndos_star+1, ') completed | Time : ', t2, ' s'
        t2 = 0.0d0 
     ENDIF
     !
  ENDDO
  !
#ifdef __PARA
  CALL mp_barrier (inter_pool_comm)
  CALL mp_sum (dos_ful,inter_pool_comm)
#ENDIF
  !
  dos_ful(1,:) = dos_ful(1,:) * ryd2ev                                      ! change to eV
  dos_ful(2,:) = dos_ful(2,:) / ryd2ev / (omega*((bohr2ang*1.0d-8)**3.0d0)) ! change to 1/eV/cm^3
  dos_ful(3,:) = dos_ful(3,:) / ryd2ev / (omega*((bohr2ang*1.0d-8)**3.0d0)) ! change to 1/eV/cm^3
  !
  !
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     fildos = TRIM(tmp_dir)//TRIM(prefix)//'.edos'
     OPEN (43234,FILE=fildos,FORM='unformatted',ACCESS='direct',RECL=DP,STATUS='replace')
     !
     WRITE(43234,REC=1) DBLE(ndos_ful)
     WRITE(43234,REC=2) DBLE(nkfdos1)
     WRITE(43234,REC=3) DBLE(nkfdos2) 
     WRITE(43234,REC=4) DBLE(nkfdos3)
     !
     DOSint = 0.0d0
     DO n = 1, ndos_ful
        IF (nspin .EQ. 1 .OR. nspin .EQ. 4) THEN 
           DOSint = DOSint + dos_ful(2,n) * (deltaE*ryd2ev)
           WRITE (43234,REC=(3*(n-1)+1)+4) dos_ful(1,n)
           WRITE (43234,REC=(3*(n-1)+2)+4) dos_ful(2,n)
           WRITE (43234,REC=(3*(n-1)+3)+4) DOSint        
        ELSE
           DOSint = DOSint + (dos_ful(2,n)+dos_ful(3,n)) * (deltaE*ryd2ev)
           WRITE (43234,REC=(3*(n-1)+1)+4) dos_ful(1,n)
           WRITE (43234,REC=(3*(n-1)+2)+4) dos_ful(2,n)
           WRITE (43234,REC=(3*(n-1)+3)+4) dos_ful(3,n)
           WRITE (43234,REC=(3*(n-1)+4)+4) DOSint  
        ENDIF 
     ENDDO 
     !
     CLOSE (43234)
     !
  ENDIF
  !
  !
  CALL CPU_TIME (t1)
  !
  WRITE (stdout,'(/5x,a,i4,a,i4,a,i4,a,f10.1,a)') 'Time spent on computing eDOS (', nkfdos1, '*', nkfdos2, '*', nkfdos3, ' k-mesh): ', &
                                                   t1-t0, ' s'
  !
  !
  DEALLOCATE (dos_ful)
  DEALLOCATE (tetra_dose)
  DEALLOCATE (et)
  DEALLOCATE (cufkk)
  DEALLOCATE (xkf_ful_dos)
  DEALLOCATE (xkf_irr_dos)
  DEALLOCATE (wkf_irr_dos) 
  DEALLOCATE (symmat_lat_dos)
  DEALLOCATE (equiv_dos)  
  DEALLOCATE (ful2irr_dos)
  !
  CALL mp_barrier (inter_pool_comm)
  !
  !
END SUBROUTINE edos_energy

!--------------------------------------------------------------------
subroutine dos_tetra (et, nspin, nbnd, nks, ntetra, tetra, e, dost)
  !------------------------------------------------------------------
  !
  USE kinds, only : DP
  implicit none
  integer :: nspin, nbnd, nks, ntetra, tetra (4, ntetra)

  real(DP) :: et (nbnd, nks), e, dost (2)
  integer :: itetra (4), nk, ns, nt, ibnd, i

  real(DP) :: etetra (4), e1, e2, e3, e4
  integer :: nspin0

  if (nspin==4) then
     nspin0=1
  else 
     nspin0=nspin
  endif
  do ns = 1, nspin0
     dost (ns) = 0.d0
     !
     ! nk is used to select k-points with up (ns=1) or down (ns=2) spin
     !
     if (ns.eq.1) then
        nk = 0
     else
        nk = nks / 2
     endif
     do nt = 1, ntetra
        do ibnd = 1, nbnd
           ! these are the energies at the vertexes of the nt-th tetrahedron
           do i = 1, 4
              etetra (i) = et (ibnd, tetra (i, nt) + nk)
           enddo
           itetra (1) = 0
           call hpsort (4, etetra, itetra)
           e1 = etetra (1)
           e2 = etetra (2)
           e3 = etetra (3)
           e4 = etetra (4)
           if (e.lt.e4.and.e.ge.e3) then
              dost (ns) = dost (ns) + 1.d0 / ntetra * (3.0d0 * (e4 - e) **2 / &
                   (e4 - e1) / (e4 - e2) / (e4 - e3) )
           elseif (e.lt.e3.and.e.ge.e2) then
              dost (ns) = dost (ns) + 1.d0 / ntetra / (e3 - e1) / (e4 - e1) &
                   * (3.0d0 * (e2 - e1) + 6.0d0 * (e-e2) - 3.0d0 * (e3 - e1 + e4 - e2) &
                   / (e3 - e2) / (e4 - e2) * (e-e2) **2)
           elseif (e.lt.e2.and.e.gt.e1) then
              dost (ns) = dost (ns) + 1.d0 / ntetra * 3.0d0 * (e-e1) **2 / &
                   (e2 - e1) / (e3 - e1) / (e4 - e1)
           endif
        enddo
     enddo
     ! add correct spin normalization : 2 for LDA, 1 for LSDA or
     ! noncollinear calculations 
     if ( nspin == 1 ) dost (ns) = dost (ns) * 2.d0
  enddo
  return
end subroutine dos_tetra


!-----------------------------------------------------------------------
SUBROUTINE tetra_generate (nk1, nk2, nk3, ntetra, tetra)
  !-----------------------------------------------------------------------
  !
  ! Tetrahedron method according to P. E. Bloechl et al, PRB49, 16223 (1994)
  !
  USE kinds, ONLY : DP
  USE edos,  ONLY : ful2irr_dos, nk_irr_dos
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)  :: nk1, nk2, nk3, ntetra
  INTEGER, INTENT(OUT) :: tetra(4,ntetra)
  ! local variable
  INTEGER              :: i, j, k, n, ip1, jp1, kp1, &
                          n1, n2, n3, n4, n5, n6, n7, n8
  !
  ! Need the map of index from fine mesh to coarse mesh (ful2irr_dos) and
  ! the number of k points in coarse mesh (nk_irr_dos)
  ! time-reversal symmetry is already considered in ful2irr_dos
  !
  DO i = 1, nk1
     DO j = 1, nk2
        DO k = 1, nk3
           !  n1-n8 are the indices of k-point 1-8 forming a cube
           ip1 = MOD(i,nk1)+1
           jp1 = MOD(j,nk2)+1
           kp1 = MOD(k,nk3)+1
           n1 = (  k-1) + (  j-1)*nk3 + (  i-1)*nk2*nk3 + 1
           n2 = (  k-1) + (  j-1)*nk3 + (ip1-1)*nk2*nk3 + 1
           n3 = (  k-1) + (jp1-1)*nk3 + (  i-1)*nk2*nk3 + 1
           n4 = (  k-1) + (jp1-1)*nk3 + (ip1-1)*nk2*nk3 + 1
           n5 = (kp1-1) + (  j-1)*nk3 + (  i-1)*nk2*nk3 + 1
           n6 = (kp1-1) + (  j-1)*nk3 + (ip1-1)*nk2*nk3 + 1
           n7 = (kp1-1) + (jp1-1)*nk3 + (  i-1)*nk2*nk3 + 1
           n8 = (kp1-1) + (jp1-1)*nk3 + (ip1-1)*nk2*nk3 + 1
           !  there are 6 tetrahedra per cube (and nk1*nk2*nk3 cubes)
           n  = 6 * ( (k-1) + (j-1)*nk3 + (i-1)*nk3*nk2 )

           tetra (1,n+1) = ful2irr_dos(n1)
           tetra (2,n+1) = ful2irr_dos(n2)
           tetra (3,n+1) = ful2irr_dos(n3)
           tetra (4,n+1) = ful2irr_dos(n6)

           tetra (1,n+2) = ful2irr_dos(n2)
           tetra (2,n+2) = ful2irr_dos(n3)
           tetra (3,n+2) = ful2irr_dos(n4)
           tetra (4,n+2) = ful2irr_dos(n6)

           tetra (1,n+3) = ful2irr_dos(n1)
           tetra (2,n+3) = ful2irr_dos(n3)
           tetra (3,n+3) = ful2irr_dos(n5)
           tetra (4,n+3) = ful2irr_dos(n6)

           tetra (1,n+4) = ful2irr_dos(n3)
           tetra (2,n+4) = ful2irr_dos(n4)
           tetra (3,n+4) = ful2irr_dos(n6)
           tetra (4,n+4) = ful2irr_dos(n8)

           tetra (1,n+5) = ful2irr_dos(n3)
           tetra (2,n+5) = ful2irr_dos(n6)
           tetra (3,n+5) = ful2irr_dos(n7)
           tetra (4,n+5) = ful2irr_dos(n8)

           tetra (1,n+6) = ful2irr_dos(n3)
           tetra (2,n+6) = ful2irr_dos(n5)
           tetra (3,n+6) = ful2irr_dos(n6)
           tetra (4,n+6) = ful2irr_dos(n7)
        ENDDO
     ENDDO
  ENDDO
  !
  ! check
  DO n = 1, ntetra
     DO i = 1, 4
        IF (tetra(i,n) .LT. 1 .OR. tetra(i,n) .GT. nk_irr_dos) &
           CALL errore ('tetra_generate','something wrong',n)
     ENDDO
  ENDDO
  !
  RETURN
  !
END SUBROUTINE tetra_generate


!----------------------------------------------------------------------------
SUBROUTINE edos_mesh (nkf1, nkf2, nkf3)
!----------------------------------------------------------------------------
#INCLUDE "f_defs.h"
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : at, bg
  USE symm_base, ONLY : set_sym_bl, s, t_rev, time_reversal, nrot
  USE edos
#ifdef __PARA
  USE mp,        ONLY : mp_barrier, mp_sum, mp_min
  USE io_global, ONLY : ionode_id, stdout
  USE mp_global, ONLY : my_pool_id, inter_pool_comm
#ENDIF
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)        :: nkf1, nkf2, nkf3
  ! 
  ! generate full and irreducible kpoints
  REAL(KIND=DP), PARAMETER   :: eps = 1.0d-5
  REAL(KIND=DP)              :: xk_rot(3), xx, yy, zz, wk_tot
  LOGICAL                    :: in_the_list
  ! full index to irreducible index
  REAL(KIND=DP), ALLOCATABLE :: wkf_ful_dos(:)
  INTEGER                    :: mini
  ! local variables
  INTEGER                    :: nk, nr, n, i, j, k, cnt, n_rot
  INTEGER                    :: ik, ibnd, ir, isym, is, it, ik0, ik_star, ik_stop
  !
#ifdef __PARA
  REAL(KIND=DP)              :: tmp
#ENDIF 
  !
  !
  !------------------------------------------------
  ! 0. clean the variables first
  !------------------------------------------------
  IF (ALLOCATED(xkf_ful_dos))     DEALLOCATE (xkf_ful_dos)
  IF (ALLOCATED(xkf_irr_dos))     DEALLOCATE (xkf_irr_dos)
  IF (ALLOCATED(wkf_irr_dos))     DEALLOCATE (wkf_irr_dos) 
  IF (ALLOCATED(symmat_lat_dos))  DEALLOCATE (symmat_lat_dos)
  IF (ALLOCATED(equiv_dos))       DEALLOCATE (equiv_dos) 
  IF (ALLOCATED(ful2irr_dos))     DEALLOCATE (ful2irr_dos) 
  !
  !------------------------------------------------
  ! 1. generate xkf_ful_dos
  !------------------------------------------------
  nk_ful_dos = nkf1*nkf2*nkf3
  !
  ALLOCATE (xkf_ful_dos(3,nk_ful_dos))
  ALLOCATE (equiv_dos(nk_ful_dos))    
  ALLOCATE (wkf_ful_dos(nk_ful_dos))    
  !
  DO i = 1, nkf1
     DO j = 1, nkf2
        DO k = 1, nkf3
           nk = (i-1)*nkf2*nkf3 + (j-1)*nkf3 + (k-1) + 1
           !
           xkf_ful_dos(1,nk) = DBLE(i-1)/nkf1
           xkf_ful_dos(2,nk) = DBLE(j-1)/nkf2
           xkf_ful_dos(3,nk) = DBLE(k-1)/nkf3
           !
        ENDDO
     ENDDO
  ENDDO
  !
  DO nk = 1, nk_ful_dos
     equiv_dos(nk) = nk
  ENDDO
  !
  !
  !------------------------------------------------
  ! 2. get s, t_rev, time_reversal and nrot
  !------------------------------------------------
  CALL set_sym_bl ()
  !
  ALLOCATE (symmat_lat_dos(3,3,48))
  symmat_lat_dos = s
  n_rot = nrot
  !
  !
  !------------------------------------------------
  ! 3. find xkf_irr_dos and wkf_irr_dos
  !------------------------------------------------
  DO nk = 1, nk_ful_dos
     IF (equiv_dos(nk) .EQ. nk) THEN
        wkf_ful_dos(nk) = 1.0d0
        DO nr = 1, n_rot
           DO i = 1, 3
              xk_rot(i) = DBLE(s(i,1,nr))*xkf_ful_dos(1,nk) + DBLE(s(i,2,nr))*xkf_ful_dos(2,nk) + DBLE(s(i,3,nr))*xkf_ful_dos(3,nk)
              xk_rot(i) = xk_rot(i) - NINT(xk_rot(i))
           ENDDO
           IF (t_rev(nr) .EQ. 1) xk_rot = -xk_rot
           xx = xk_rot(1)*nkf1
           yy = xk_rot(2)*nkf2
           zz = xk_rot(3)*nkf3
           in_the_list = (ABS(xx-NINT(xx)) .LE. eps) .AND. (ABS(yy-NINT(yy)) .LE. eps) .AND. (ABS(zz-NINT(zz)) .LE. eps) 
           IF (in_the_list) THEN
              i = MOD ( NINT(xk_rot(1)*nkf1 + 2*nkf1), nkf1 ) + 1
              j = MOD ( NINT(xk_rot(2)*nkf2 + 2*nkf2), nkf2 ) + 1
              k = MOD ( NINT(xk_rot(3)*nkf3 + 2*nkf3), nkf3 ) + 1
              n = (i-1)*nkf2*nkf3 + (j-1)*nkf3 + (k-1) + 1
              IF ((n .GT. nk) .AND. (equiv_dos(n) .EQ. n)) then
                 equiv_dos(n) = nk
                 wkf_ful_dos(nk) = wkf_ful_dos(nk) + 1.0d0
              ELSE
                 IF ((equiv_dos(n) .NE. nk) .OR. (n .LT. nk)) CALL errore ('edos_mesh', 'something wrong in the checking algorithm',1)
              ENDIF
           ENDIF
           IF (time_reversal) THEN
              xx = -xk_rot(1)*nkf1
              yy = -xk_rot(2)*nkf2
              zz = -xk_rot(3)*nkf3
              in_the_list = (ABS(xx-NINT(xx)) .LE. eps) .AND. (ABS(yy-NINT(yy)) .LE. eps) .AND. (ABS(zz-NINT(zz)) .LE. eps)  
              IF (in_the_list) THEN
                 i = MOD ( NINT(-xk_rot(1)*nkf1 + 2*nkf1), nkf1 ) + 1
                 j = MOD ( NINT(-xk_rot(2)*nkf2 + 2*nkf2), nkf2 ) + 1
                 k = MOD ( NINT(-xk_rot(3)*nkf3 + 2*nkf3), nkf3 ) + 1
                 n = (i-1)*nkf2*nkf3 + (j-1)*nkf3 + (k-1) + 1
                 IF ((n .GT. nk) .AND. (equiv_dos(n) .EQ. n)) then
                    equiv_dos(n) = nk
                    wkf_ful_dos(nk) = wkf_ful_dos(nk) + 1.0d0
                 ELSE
                    IF ((equiv_dos(n) .NE. nk) .OR. (n .LT. nk)) CALL errore ('edos_mesh', 'something wrong in the checking algorithm',1)
                 ENDIF
              ENDIF
           ENDIF
        ENDDO
     ENDIF
  ENDDO
  !
  nk_irr_dos = 0
  DO nk = 1, nk_ful_dos
     IF (equiv_dos(nk) .EQ. nk) nk_irr_dos = nk_irr_dos + 1
  ENDDO
  !
  ALLOCATE (xkf_irr_dos(3,nk_irr_dos))
  ALLOCATE (wkf_irr_dos(nk_irr_dos)) 
  !
  n = 0
  wk_tot = 0.0d0
  DO nk = 1, nk_ful_dos
     IF (equiv_dos(nk) .EQ. nk) THEN
        n = n + 1
        wkf_irr_dos(n) = wkf_ful_dos(nk)
        wk_tot = wk_tot + wkf_irr_dos(n)
        DO i = 1, 3
           xkf_irr_dos(i,n) = xkf_ful_dos(i,nk)-NINT(xkf_ful_dos(i,nk))
        ENDDO
     ENDIF
  ENDDO
  !
  !  normalize weights to one
  DO nk=1,nk_irr_dos
     wkf_irr_dos(nk) = wkf_irr_dos(nk) / wk_tot
  ENDDO
  !
  !
  !------------------------------------------------
  ! 4. find ful2irr_dos
  !------------------------------------------------
  !
  ALLOCATE (ful2irr_dos(nk_ful_dos))
  ful2irr_dos = 0
  !
  ik = 0
  DO nk = 1, nk_ful_dos
     IF (equiv_dos(nk) .EQ. nk) THEN
        ik = ik + 1
        ful2irr_dos(nk) = ik
     ENDIF
  ENDDO
  !
  DO nk = 1, nk_ful_dos
     ful2irr_dos(nk) = ful2irr_dos(equiv_dos(nk))
  ENDDO
  !
  !
  DEALLOCATE (wkf_ful_dos)  
  !
  !
  CALL mp_barrier(inter_pool_comm)
  !
  !
END SUBROUTINE edos_mesh
