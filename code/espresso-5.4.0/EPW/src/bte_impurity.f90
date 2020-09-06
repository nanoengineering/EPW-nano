!-----------------------------------------------------------------------
SUBROUTINE impurity_factor (iq, fac_imp)
!-----------------------------------------------------------------------
#INCLUDE "f_defs.h"
  USE kinds,         ONLY : DP
  USE epwcom,        ONLY : z_imp, epdope
  USE pwcom,         ONLY : at, bg, alat, omega
  USE elph2,         ONLY : epsi, epsil_dyn
  USE constants_epw, ONLY : pi, twopi, au2cm
  USE bte_var            
#ifdef __PARA
  USE io_global,     ONLY : stdout
  USE mp_global,     ONLY : my_pool_id, npool
#ENDIF
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)        :: iq
  REAL(KIND=DP), INTENT(OUT) :: fac_imp
  ! local variable
  REAL(KIND=DP)              :: xxq(3), doping
  REAL(KIND=DP)              :: coef, qeqe, olap
  !
  ! 
  ! previous coefficient
  doping = epdope(1)
  coef = 128.0d0*pi*pi*pi*z_imp*z_imp*(ABS(doping)*au2cm*au2cm*au2cm) / omega
  !
  ! qeq * epsil_dyn
  xxq(:) = xqf_ful(:,iq)
  !
  CALL cryst_to_cart (1, xxq, bg, 1) 
  !
  qeqe = ( xxq(1) * (epsi(1,1)*xxq(1)+epsi(1,2)*xxq(2)+epsi(1,3)*xxq(3) ) + &
           xxq(2) * (epsi(2,1)*xxq(1)+epsi(2,2)*xxq(2)+epsi(2,3)*xxq(3) ) + &
           xxq(3) * (epsi(3,1)*xxq(1)+epsi(3,2)*xxq(2)+epsi(3,3)*xxq(3) ) ) * twopi/alat
  !
  qeqe = qeqe * epsil_dyn(ful2irr_q(iq)) 
  !
  ! overlap factor
  olap = 1.0d0 ! for N-process
  !
  !
  fac_imp = coef*olap/qeqe/qeqe
  !
  !
END SUBROUTINE impurity_factor



!-----------------------------------------------------------------------
SUBROUTINE epsilon_dynamic (etf)
!-----------------------------------------------------------------------
  !
  !  Compute the actual Lindhard dielectric constant screening due to
  !  the electrons through the RPA Lindhaed potential function
  !
  !-----------------------------------------------------------------------
#INCLUDE "f_defs.h"
  USE kinds,         ONLY : DP
  USE pwcom,         ONLY : at, bg, ef, alat, omega
  USE phcom,         ONLY : nmodes
  USE epwcom,        ONLY : nbndsub, efermi_read, fermi_energy, eptemp, degaussw, &
                            nkf1, nkf2, nkf3, nqf1, nqf2, nqf3
  USE elph2,         ONLY : epsil_dyn, epsi, ef_epw
  USE constants_epw, ONLY : twopi, fpi, ryd2ev
  USE bte_func
  USE bte_var             
#ifdef __PARA
  USE mp,            ONLY : mp_bcast, mp_sum, mp_barrier
  USE mp_global,     ONLY : inter_pool_comm, my_pool_id
  USE io_global,     ONLY : ionode_id, stdout
#ENDIF
  !
  IMPLICIT NONE
  !
  REAL(KIND=DP), INTENT(IN)  :: etf(nbndsub,nk_irr)
  REAL(KIND=DP)              :: wgk(nbndsub,nk_irr), qeq, uq, lindh
  INTEGER                    :: ik, iq, ikq, ibnd, jbnd, isym, ik_star, ik_stop, iq_star, iq_stop, iq_map, itemp
  REAL(KIND=DP)              :: eptemp0, ef0, ekk, ekq, fkk, fkq
  INTEGER, DIMENSION(3)      :: ijk_k, ijk_q, ijk_kq, ijk_q_rot, ijk_q_fbz
  REAL(KIND=DP)              :: q(3), q_rot(3), q_fbz(3), non_r
  ! external fuction
  REAL(KIND=DP), EXTERNAL    :: wgauss
  REAL(KIND=DP)              :: t0, t1, tt
  !
  !
  WRITE(stdout,'(5x,a/)') 'Start computing epsilon(q,w=0)'
  !
  tt = 0.0d0
  !
  ! initialization of system
  eptemp0 = eptemp(1)
  IF (efermi_read) THEN
     ef0 = fermi_energy
  ELSE
     ef0 = ef_epw(1,1)
  ENDIF
  !
  !
  ! FD distribution in irr-BZ
  CALL para_bounds (ik_star, ik_stop, nk_irr)
  wgk = 0.0d0
  !
  DO ik = ik_star, ik_stop
     DO ibnd = 1, nbndsub
        wgk(ibnd,ik) = wgauss(-(etf(ibnd,ik)-ef0)/eptemp0,-99)  
     ENDDO
  ENDDO
  !
#ifdef __PARA
  CALL mp_barrier (inter_pool_comm)
  CALL mp_sum (wgk,inter_pool_comm)
#ENDIF
  !
  !
  ! lindhard function
  ALLOCATE (epsil_dyn(nq_irr))
  epsil_dyn = 0.0d0
  !
  CALL para_bounds (iq_star, iq_stop, nq_irr)
  !
  OPEN (1111,FILE='BTE/META/xqf_irr_cryst',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='old')
  OPEN (2222,FILE='BTE/META/irr2ful_sym',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='old')  
  ! 
  DO iq = iq_star, iq_stop
     !
     CALL CPU_TIME(t0)
     !
     READ (1111,REC=iq) q(1:3)
     !
     CALL cryst_to_cart (1, q, bg, 1) 
     !
     qeq = ( q(1) * (epsi(1,1)*q(1)+epsi(1,2)*q(2)+epsi(1,3)*q(3) ) + &
             q(2) * (epsi(2,1)*q(1)+epsi(2,2)*q(2)+epsi(2,3)*q(3) ) + &
             q(3) * (epsi(3,1)*q(1)+epsi(3,2)*q(2)+epsi(3,3)*q(3) ) ) * twopi/alat
     !
     CALL cryst_to_cart (1, q, at, -1) 
     !
     IF (qeq .NE. 0.0d0) THEN
        uq = fpi*2.0d0/qeq/omega
     ELSE
        uq = 0.0d0
        GOTO 987
     ENDIF
     !
     lindh = 0.0d0
     !
     DO ik = 1, nk_ful
        !
        READ (2222,REC=ik) isym
        !
        q_rot(:) = q(1)*symmat_lat(:,1,isym) + q(2)*symmat_lat(:,2,isym) + q(3)*symmat_lat(:,3,isym)
        !
        ijk_q_rot(1) = NINT(q_rot(1)*DBLE(nqf1))
        ijk_q_rot(2) = NINT(q_rot(2)*DBLE(nqf2))
        ijk_q_rot(3) = NINT(q_rot(3)*DBLE(nqf3))
        ! 
        ijk_k =  id2ijk(ik)
        ijk_kq = ijk_fbz(ijk_k,ijk_q_rot)
        !
        ikq = ijk2id(ijk_kq)
        !
        DO ibnd = 1, nbndsub
           !
           ekk = etf(ibnd,ful2irr(ik)) ! - ef0
           fkk = wgk(ibnd,ful2irr(ik))
           !
           DO jbnd = 1, nbndsub
              !
              ekq = etf(jbnd,ful2irr(ikq)) ! - ef0
              !
              IF (ekk .NE. ekq) THEN
                 !               
                 fkq = wgk(jbnd,ful2irr(ikq))
                 !
                 lindh = lindh + (fkq-fkk)/(ekq-ekk)
                 !
              ENDIF
              ! 
           ENDDO
           !
        ENDDO
        !
     ENDDO
     !
987 CONTINUE
     epsil_dyn(iq) = 1.0d0 - uq * lindh / DBLE(nk_ful)
     !
     CALL CPU_TIME(t1)
     tt = tt + (t1-t0)
     IF (MOD(iq,200) .EQ. 0 .OR. MOD(iq,iq_stop-iq_star+1) .EQ. 0) THEN
         WRITE(stdout,'(13x,a,i6,a,i6,a,f7.1,a)') '(', iq, '/', iq_stop-iq_star+1,') completed, ', tt, ' s'
         tt = 0.0d0
     ENDIF
     !
  ENDDO
  !
  CLOSE (1111)
  CLOSE (2222)
  !
#ifdef __PARA
  CALL mp_barrier (inter_pool_comm)
  CALL mp_sum (epsil_dyn,inter_pool_comm)
#ENDIF
  !
  ! export
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     OPEN (99999,FILE='BTE/META/epsil_dyn',FORM='unformatted',ACCESS='direct',RECL=DP,STATUS='replace')
     DO iq = 1, nq_irr
        WRITE (99999,REC=iq) epsil_dyn(iq)
     ENDDO
     CLOSE (99999)
     !
     !
     OPEN (77777,FILE='BTE/META/irr2ful_q',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='old')
     OPEN (88888,FILE='BTE/META/xqf_fbz_cart',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='old')
     OPEN (99999,file='BTE/EPCHECK/epsil_dyn.dat',STATUS='replace')
     DO iq = 1, nq_irr
        !
        READ (77777,REC=iq) iq_map
        READ (88888,REC=iq_map) q_fbz(1:3)
        !
        WRITE (99999,'(2es13.4)') SQRT(DOT_PRODUCT(q_fbz,q_fbz)), epsil_dyn(iq)
        !
     ENDDO
     CLOSE (77777)
     CLOSE (88888)
     CLOSE (99999)
     !
  ENDIF
  !
  !
END SUBROUTINE epsilon_dynamic
