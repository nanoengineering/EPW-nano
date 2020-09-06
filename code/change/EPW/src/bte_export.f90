!----------------------------------------------------------------------------
SUBROUTINE export_cumu (iter, converge, itemp, idope)
!----------------------------------------------------------------------------
! 
! export cumulative transport properties
!
! mfp (nm), cod_x, etc_x, seb_x, pof_x ... (y and z directions)
!
! export cumulative phonon drag (contribution by phonons) properties, added by Qian, Jan 2018 
!
#INCLUDE "f_defs.h"
  USE kinds,         ONLY : DP
  USE phcom,         ONLY : nmodes
  USE epwcom,        ONLY : nbndsub, epdim, phdrag
  USE elph2,         ONLY : gammai_mode_all, ibndmin, ibndmax, nbnd_red, wf_all, etf_all, F_k_ful, &
                            sigmai_mode_all_abs, sigmai_mode_all_emi, mfp_k_ful, mfp_q_ful
  USE constants_epw, ONLY : au2nm, ryd2thz, twopi, ryd2ev
  USE bte_var
  USE bte_func
#ifdef __PARA
  USE mp,            ONLY : mp_bcast, mp_sum
  USE io_global,     ONLY : ionode_id
  USE mp,            ONLY : mp_barrier
  USE mp_global,     ONLY : my_pool_id, inter_pool_comm
#ENDIF
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)        :: iter, itemp, idope
  LOGICAL, INTENT(IN)        :: converge
  INTEGER                    :: ip_star, ip_stop
  REAL(KIND=DP)              :: L11(3,3), L12(3,3), L21(3,3), L22(3,3), &
                                L11_tot(3,3), inv_L11(3,3), inv_L11_tot(3,3), &
                                L11_k(3,3), L12_k(3,3), L21_k(3,3), L22_k(3,3), &
                                L11_upper(3,3), L12_upper(3,3), L21_upper(3,3), L22_upper(3,3), &
                                L11_lower(3,3), L12_lower(3,3), L21_lower(3,3), L22_lower(3,3)
  !
  REAL(KIND=DP)              :: L33(3,3), L44(3,3), L33_q(3,3), L44_q(3,3), wqq, L33_m(3,3,nmodes), L44_m(3,3,nmodes),&
                                length, xq_fbz(3), xk_fbz(3), L44_e(3,3), L44_k(3,3), ekk, sigmai, xful_cry(3)
  REAL(KIND=DP), ALLOCATABLE :: mfp(:), mfp_ph(:), freq_ph(:)
  REAL(KIND=DP)              :: mfp_max, mfp_min, mfp_del, mfp_ph_max, mfp_ph_min, mfp_ph_del, freq_max, freq_del
  INTEGER                    :: mfp_num = 1000
  !
  REAL(KIND=DP), ALLOCATABLE :: cod(:,:,:), etc(:,:,:), seb(:,:,:), pof(:,:,:), spd(:,:,:), ltc(:,:,:), &
                                cod_upper(:,:,:), etc_upper(:,:,:), seb_upper(:,:,:), pof_upper(:,:,:), &
                                cod_lower(:,:,:), etc_lower(:,:,:), seb_lower(:,:,:), pof_lower(:,:,:), &
                                spd_f(:,:,:), ltc_f(:,:,:),spd_f_m(:,:,:,:), ltc_f_m(:,:,:,:), &
                                spd_m(:,:,:,:), ltc_m(:,:,:,:), spd_e(:,:,:), spd_x(:,:,:), spd_l(:,:,:), &
                                seb_x(:,:,:), seb_l(:,:,:), seb_g(:,:,:), spd_g(:,:,:)
  !
  INTEGER                    :: ik, ik_ful, ik_irr_red, ik_ful_red, iq, iq_ful, iq_irr_red, ibnd, ibnd0, imode,&
                                i, j, ip, ir, uorl
  CHARACTER(LEN=256)         :: onsager_ufmt, phdrag_ufmt, phdragcumu_MFP_fmt, cumu_fmt, cumu_upper_fmt, &
                                cumu_lower_fmt, phdragcumu_freq_fmt, rate_RTA_3Dphdrag_fmt, phdragcumu_k_MFP_fmt, &
                                rate_k_RTA_3Dphdrag_fmt, phdrag_k_ufmt
  CHARACTER(LEN=12)          :: txnx
  CHARACTER(LEN=3)           :: itemp_num, idope_num
  !
  !
  IF (converge) THEN
     !
     WRITE(itemp_num,'(i3)') itemp
     WRITE(idope_num,'(i3)') idope
     txnx = 'T'//TRIM(ADJUSTL(itemp_num))//'_N'//TRIM(ADJUSTL(idope_num))
     onsager_ufmt = 'BTE/META/onsager_'//TRIM(ADJUSTL(txnx))
     !
     phdrag_ufmt = 'BTE/META/phdrag_q_'//TRIM(ADJUSTL(txnx))
     phdrag_k_ufmt = 'BTE/META/phdrag_k_'//TRIM(ADJUSTL(txnx))
     !
     !
     !cumu_fmt = 'BTE/'//TRIM(ADJUSTL(txnx))//'/cumu_'//TRIM(ADJUSTL(txnx))//'.dat'
     !cumu_upper_fmt = 'BTE/'//TRIM(ADJUSTL(txnx))//'/cumu_upper_'//TRIM(ADJUSTL(txnx))//'.dat'
     !cumu_lower_fmt = 'BTE/'//TRIM(ADJUSTL(txnx))//'/cumu_lower_'//TRIM(ADJUSTL(txnx))//'.dat'
     phdragcumu_MFP_fmt = 'BTE/'//TRIM(ADJUSTL(txnx))//'/phdrag_MFPcumu.dat'
     phdragcumu_k_MFP_fmt = 'BTE/'//TRIM(ADJUSTL(txnx))//'/phdrag_k_MFPcumu.dat'
     phdragcumu_freq_fmt = 'BTE/'//TRIM(ADJUSTL(txnx))//'/phdrag_freqcumu.dat'
     cumu_fmt = 'BTE/'//TRIM(ADJUSTL(txnx))//'/cumu.dat'
     cumu_upper_fmt = 'BTE/'//TRIM(ADJUSTL(txnx))//'/cumu_upper.dat'
     cumu_lower_fmt = 'BTE/'//TRIM(ADJUSTL(txnx))//'/cumu_lower.dat'
     !
     ! 1. find the maximum and minimun of mfp, freq
     mfp_min = +1.0d+30
     mfp_max = -1.0d+30
     mfp_ph_min = +1.0d+30
     mfp_ph_max = -1.0d+30
     freq_max  = 0.0d0
     !
     DO ik = 1, nk_ful_red
        DO ibnd = 1, nbnd_red
           !
           IF (mfp_k_ful(ibnd,ik) .LT. mfp_min .AND. mfp_k_ful(ibnd,ik) .GE. 0.0d0) mfp_min = mfp_k_ful(ibnd,ik)
           IF (mfp_k_ful(ibnd,ik) .GT. mfp_max) mfp_max = mfp_k_ful(ibnd,ik)
           !
        ENDDO
     ENDDO
     !
     mfp_del = (mfp_max-mfp_min)/DBLE(mfp_num-1)
     !
     IF (phdrag) THEN 
      DO iq = 1, nq_ful_red
        iq_irr_red = rful2rirr_q(iq)
        DO imode = 1, nmodes
           !
           IF (mfp_q_ful(imode,iq) .LT. mfp_ph_min .AND. mfp_q_ful(imode,iq) .GE. 0.0d0) mfp_ph_min = mfp_q_ful(imode,iq)
           IF (mfp_q_ful(imode,iq) .GT. mfp_ph_max) mfp_ph_max = mfp_q_ful(imode,iq)
           !
           IF (wf_all(imode,iq_irr_red) .GT. freq_max) freq_max = wf_all(imode,iq_irr_red)
           !
        ENDDO
      ENDDO
      !
      mfp_ph_del = (mfp_ph_max-mfp_ph_min)/DBLE(mfp_num-1)
      freq_del = (freq_max)/DBLE(mfp_num-1)
      !
     ENDIF
     !
     ! 2. compute L11_tot and inv_L11_tot
     IF (my_pool_id .EQ. ionode_id) THEN
        !
        L11_tot     = 0.0d0
        inv_L11_tot = 0.0d0
        !
        OPEN (7777,FILE=onsager_ufmt,FORM='unformatted',ACCESS='direct',RECL=4+9*4*DP)
        !
        DO ik = 1, nk_ful_red
           DO ibnd = 1, nbnd_red
              !
              READ (7777,REC=(ik-1)*nbnd_red+ibnd) L11_k(1:3,1:3), L12_k(1:3,1:3), L21_k(1:3,1:3), L22_k(1:3,1:3), uorl
              L11_tot = L11_tot + L11_k
              !
           ENDDO
        ENDDO
        CLOSE (7777)
        !
        IF (epdim .EQ. 2) THEN
           L11_tot(3,3) = 1.0d0
        ENDIF
        !
        inv_L11_tot = imat(L11_tot)
        !
     ENDIF
     !
#ifdef __PARA
     CALL mp_bcast (inv_L11_tot,ionode_id,inter_pool_comm)
     CALL mp_barrier (inter_pool_comm)
#ENDIF
     !
     !
     ! 3. compute cumulative properties
     CALL para_bounds (ip_star, ip_stop, mfp_num)
     !
     ALLOCATE (mfp(mfp_num))
     IF (phdrag) ALLOCATE (mfp_ph(mfp_num))
     IF (phdrag) ALLOCATE (freq_ph(mfp_num))
     ALLOCATE (cod(3,3,mfp_num))
     ALLOCATE (etc(3,3,mfp_num))
     ALLOCATE (seb(3,3,mfp_num))
     ALLOCATE (seb_x(3,3,mfp_num))
     ALLOCATE (seb_l(3,3,mfp_num))
     ALLOCATE (seb_g(3,3,mfp_num))
     ALLOCATE (pof(3,3,mfp_num))
     IF (phdrag) ALLOCATE (spd_e(3,3,mfp_num))
     IF (phdrag) ALLOCATE (spd_x(3,3,mfp_num))
     IF (phdrag) ALLOCATE (spd_l(3,3,mfp_num))
     IF (phdrag) ALLOCATE (spd_g(3,3,mfp_num))
     IF (phdrag) ALLOCATE (spd(3,3,mfp_num))
     IF (phdrag) ALLOCATE (ltc(3,3,mfp_num))
     IF (phdrag) ALLOCATE (spd_f(3,3,mfp_num))
     IF (phdrag) ALLOCATE (ltc_f(3,3,mfp_num))
     IF (phdrag) ALLOCATE (spd_m(3,3,nmodes,mfp_num))
     IF (phdrag) ALLOCATE (ltc_m(3,3,nmodes,mfp_num))
     IF (phdrag) ALLOCATE (spd_f_m(3,3,nmodes,mfp_num))
     IF (phdrag) ALLOCATE (ltc_f_m(3,3,nmodes,mfp_num))
     ALLOCATE (cod_upper(3,3,mfp_num))
     ALLOCATE (etc_upper(3,3,mfp_num))
     ALLOCATE (seb_upper(3,3,mfp_num))
     ALLOCATE (pof_upper(3,3,mfp_num))
     ALLOCATE (cod_lower(3,3,mfp_num))
     ALLOCATE (etc_lower(3,3,mfp_num))
     ALLOCATE (seb_lower(3,3,mfp_num))
     ALLOCATE (pof_lower(3,3,mfp_num))
     mfp = 0.0d0
     IF (phdrag) mfp_ph = 0.0d0
     IF (phdrag) freq_ph = 0.0d0
     cod = 0.0d0
     etc = 0.0d0
     seb = 0.0d0
     seb_x = 0.0d0
     seb_l = 0.0d0
     seb_g = 0.0d0
     pof = 0.0d0
     IF (phdrag) ltc = 0.0d0
     IF (phdrag) spd = 0.0d0
     IF (phdrag) spd_e = 0.0d0
     IF (phdrag) spd_x = 0.0d0
     IF (phdrag) spd_l = 0.0d0
     IF (phdrag) spd_g = 0.0d0
     IF (phdrag) ltc_f = 0.0d0
     IF (phdrag) spd_f = 0.0d0
     IF (phdrag) ltc_m = 0.0d0
     IF (phdrag) spd_m = 0.0d0
     IF (phdrag) ltc_f_m = 0.0d0
     IF (phdrag) spd_f_m = 0.0d0
     cod_upper = 0.0d0
     etc_upper = 0.0d0
     seb_upper = 0.0d0
     pof_upper = 0.0d0
     cod_lower = 0.0d0
     etc_lower = 0.0d0
     seb_lower = 0.0d0
     pof_lower = 0.0d0
     !
     OPEN (8888,FILE=onsager_ufmt,FORM='unformatted',ACCESS='direct',RECL=4+9*4*DP,STATUS='old')
     IF (phdrag) OPEN (33333,FILE=phdrag_k_ufmt,FORM='unformatted',ACCESS='direct',RECL=(1+9)*DP,STATUS='old')
     OPEN (7777,FILE='BTE/META/rful2ful',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='old')
     OPEN (6666,FILE='BTE/META/xkf_ful_cryst',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='old')
     !
     DO ip = ip_star, ip_stop
        !
        mfp(ip) = mfp_min + DBLE(ip-1)*mfp_del
        !
        L11   = 0.0d0
        L12   = 0.0d0
        L21   = 0.0d0
        L22   = 0.0d0
        IF (phdrag) L44_e = 0.0d0
        L11_upper = 0.0d0
        L12_upper = 0.0d0
        L21_upper = 0.0d0
        L22_upper = 0.0d0
        L11_lower = 0.0d0
        L12_lower = 0.0d0
        L21_lower = 0.0d0
        L22_lower = 0.0d0
        !
        DO ik = 1, nk_ful_red
           !
           READ (7777,REC=ik) ik_ful
           READ (6666,REC=ik_ful) xful_cry(1:3)
           !
           DO ibnd = 1, nbnd_red
              !
              IF (mfp_k_ful(ibnd,ik) .LE. mfp(ip)) THEN
                 READ (8888,REC=(ik-1)*nbnd_red+ibnd) L11_k(1:3,1:3), L12_k(1:3,1:3), L21_k(1:3,1:3), L22_k(1:3,1:3), uorl
                 IF (phdrag) READ (33333,REC=(ik-1)*nbnd_red+ibnd)  ekk, L44_k(1:3,1:3)
              ELSE
                 !
                 L11_k = 0.0d0
                 L12_k = 0.0d0
                 L21_k = 0.0d0
                 L22_k = 0.0d0   
                 IF (phdrag) L44_k = 0.0d0
                 !  
              ENDIF
              !
              L11 = L11 + L11_k
              L12 = L12 + L12_k
              L21 = L21 + L21_k
              L22 = L22 + L22_k
              !
              ! seperate valleys' contributions
        IF ((SQRT((xful_cry(1)-0.5d0)**2 + (xful_cry(2)-0.5d0)**2+ (xful_cry(3)-0.5d0)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-0.5d0)**2 + xful_cry(2)**2+ xful_cry(3)**2) .LT. 0.0625) .OR. (SQRT(xful_cry(1)**2 + (xful_cry(2)-0.5d0)**2+ xful_cry(3)**2) .LT. 0.0625) .OR. (SQRT(xful_cry(1)**2 + xful_cry(2)**2+ (xful_cry(3)-0.5d0)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-0.5d0)**2 + (xful_cry(2)-1.0d0)**2+ (xful_cry(3)-1.0d0)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-1.0d0)**2 + (xful_cry(2)-0.5d0)**2+ (xful_cry(3)-1.0d0)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-1.0d0)**2 + (xful_cry(2)-1.0d0)**2+ (xful_cry(3)-0.5d0)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-0.5d0)**2 + (xful_cry(2)-1.0d0)**2+ xful_cry(3)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-1.0d0)**2 + (xful_cry(2)-0.5d0)**2+ xful_cry(3)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-1.0d0)**2 + xful_cry(2)**2+ (xful_cry(3)-0.5d0)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-0.5d0)**2 + xful_cry(2)**2+ (xful_cry(3)-1.0d0)**2) .LT. 0.0625) .OR. (SQRT(xful_cry(1)**2 + (xful_cry(2)-0.5d0)**2+ (xful_cry(3)-1.0d0)**2) .LT. 0.0625) .OR. (SQRT(xful_cry(1)**2 + (xful_cry(2)-1.0d0)**2+ (xful_cry(3)-0.5d0)**2) .LT. 0.0625)) THEN  ! crystal: L 0.5 0.5 0.5, 0/1 0/1 0.5, 0/1 0.5 0/1, 0.5 0/1 0/1
        !
        seb_l(:,:,ip) =  seb_l(:,:,ip) + L12_k
        !
        ELSEIF ((SQRT((xful_cry(1)-0.425d0)**2 + (xful_cry(2)-0.425d0)**2+ xful_cry(3)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-0.425d0)**2 + (xful_cry(3)-0.425d0)**2+ xful_cry(2)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(3)-0.425d0)**2 + (xful_cry(2)-0.425d0)**2+ xful_cry(1)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-0.425d0)**2 + (xful_cry(2)-0.425d0)**2+ (xful_cry(3)-1.0d0)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-0.425d0)**2 + (xful_cry(3)-0.425d0)**2+ (xful_cry(2)-1.0d0)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(3)-0.425d0)**2 + (xful_cry(2)-0.425d0)**2+ (xful_cry(1)-1.0d0)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-0.575d0)**2 + (xful_cry(2)-0.575d0)**2+ xful_cry(3)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-0.575d0)**2 + (xful_cry(3)-0.575d0)**2+ xful_cry(2)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(3)-0.575d0)**2 + (xful_cry(2)-0.575d0)**2+ xful_cry(1)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-0.575d0)**2 + (xful_cry(2)-0.575d0)**2+ (xful_cry(3)-1.0d0)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-0.575d0)**2 + (xful_cry(3)-0.575d0)**2+ (xful_cry(2)-1.0d0)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(3)-0.575d0)**2 + (xful_cry(2)-0.575d0)**2+ (xful_cry(1)-1.0d0)**2) .LT. 0.0625)) THEN   ! Si crystal: near X 0/1 0.425 0.425 or 0/1 0.575 0.575
        !
        seb_x(:,:,ip) =  seb_x(:,:,ip) + L12_k
        !
        ELSEIF ((SQRT(xful_cry(1)**2+ xful_cry(2)**2 + xful_cry(3)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-1.0d0)**2+ (xful_cry(2)-1.0d0)**2 + (xful_cry(3)-1.0d0)**2) .LT. 0.0625)) THEN  !gamma
        !
        seb_g(:,:,ip) =  seb_g(:,:,ip) + L12_k
        !
        ENDIF ! select coordinates
        !
        !
              !
              IF (phdrag) THEN 
                L44_e = L44_e + L44_k
              !
              ! seperate valleys' contributions
        IF ((SQRT((xful_cry(1)-0.5d0)**2 + (xful_cry(2)-0.5d0)**2+ (xful_cry(3)-0.5d0)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-0.5d0)**2 + xful_cry(2)**2+ xful_cry(3)**2) .LT. 0.0625) .OR. (SQRT(xful_cry(1)**2 + (xful_cry(2)-0.5d0)**2+ xful_cry(3)**2) .LT. 0.0625) .OR. (SQRT(xful_cry(1)**2 + xful_cry(2)**2+ (xful_cry(3)-0.5d0)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-0.5d0)**2 + (xful_cry(2)-1.0d0)**2+ (xful_cry(3)-1.0d0)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-1.0d0)**2 + (xful_cry(2)-0.5d0)**2+ (xful_cry(3)-1.0d0)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-1.0d0)**2 + (xful_cry(2)-1.0d0)**2+ (xful_cry(3)-0.5d0)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-0.5d0)**2 + (xful_cry(2)-1.0d0)**2+ xful_cry(3)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-1.0d0)**2 + (xful_cry(2)-0.5d0)**2+ xful_cry(3)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-1.0d0)**2 + xful_cry(2)**2+ (xful_cry(3)-0.5d0)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-0.5d0)**2 + xful_cry(2)**2+ (xful_cry(3)-1.0d0)**2) .LT. 0.0625) .OR. (SQRT(xful_cry(1)**2 + (xful_cry(2)-0.5d0)**2+ (xful_cry(3)-1.0d0)**2) .LT. 0.0625) .OR. (SQRT(xful_cry(1)**2 + (xful_cry(2)-1.0d0)**2+ (xful_cry(3)-0.5d0)**2) .LT. 0.0625)) THEN  ! crystal: L 0.5 0.5 0.5, 0/1 0/1 0.5, 0/1 0.5 0/1, 0.5 0/1 0/1
        !
        spd_l(:,:,ip) =  spd_l(:,:,ip) + L44_k
        !
        ELSEIF ((SQRT((xful_cry(1)-0.425d0)**2 + (xful_cry(2)-0.425d0)**2+ xful_cry(3)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-0.425d0)**2 + (xful_cry(3)-0.425d0)**2+ xful_cry(2)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(3)-0.425d0)**2 + (xful_cry(2)-0.425d0)**2+ xful_cry(1)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-0.425d0)**2 + (xful_cry(2)-0.425d0)**2+ (xful_cry(3)-1.0d0)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-0.425d0)**2 + (xful_cry(3)-0.425d0)**2+ (xful_cry(2)-1.0d0)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(3)-0.425d0)**2 + (xful_cry(2)-0.425d0)**2+ (xful_cry(1)-1.0d0)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-0.575d0)**2 + (xful_cry(2)-0.575d0)**2+ xful_cry(3)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-0.575d0)**2 + (xful_cry(3)-0.575d0)**2+ xful_cry(2)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(3)-0.575d0)**2 + (xful_cry(2)-0.575d0)**2+ xful_cry(1)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-0.575d0)**2 + (xful_cry(2)-0.575d0)**2+ (xful_cry(3)-1.0d0)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-0.575d0)**2 + (xful_cry(3)-0.575d0)**2+ (xful_cry(2)-1.0d0)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(3)-0.575d0)**2 + (xful_cry(2)-0.575d0)**2+ (xful_cry(1)-1.0d0)**2) .LT. 0.0625)) THEN   ! Si crystal: near X 0/1 0.425 0.425 or 0/1 0.575 0.575
        !
        spd_x(:,:,ip) =  spd_x(:,:,ip) + L44_k
        !
        ELSEIF ((SQRT(xful_cry(1)**2+ xful_cry(2)**2 + xful_cry(3)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-1.0d0)**2+ (xful_cry(2)-1.0d0)**2 + (xful_cry(3)-1.0d0)**2) .LT. 0.0625)) THEN  !gamma
        !
        !
        spd_g(:,:,ip) =  spd_g(:,:,ip) + L44_k
        !
        ENDIF ! select coordinates
        !
              ENDIF ! if phonon drag
              !
              IF (uorl .EQ. +1) THEN
                 L11_upper = L11_upper + L11_k
                 L12_upper = L12_upper + L12_k
                 L21_upper = L21_upper + L21_k
                 L22_upper = L22_upper + L22_k
              ENDIF
              !
              IF (uorl .EQ. -1) THEN
                 L11_lower = L11_lower + L11_k
                 L12_lower = L12_lower + L12_k
                 L21_lower = L21_lower + L21_k
                 L22_lower = L22_lower + L22_k
              ENDIF
              !
           ENDDO !ibnd
        ENDDO !nk_ful_red
        !
        IF (epdim .EQ. 2) THEN
           !
           L11(3,3)   = 1.0d0
           L12(3,3)   = 1.0d0
           L21(3,3)   = 1.0d0
           L22(3,3)   = 1.0d0
           IF (phdrag) L44_e(3,3)   = 1.0d0
           L11_upper(3,3) = 1.0d0
           L12_upper(3,3) = 1.0d0
           L21_upper(3,3) = 1.0d0
           L22_upper(3,3) = 1.0d0
           L11_lower(3,3) = 1.0d0
           L12_lower(3,3) = 1.0d0
           L21_lower(3,3) = 1.0d0
           L22_lower(3,3) = 1.0d0
           !
        ENDIF
        !
        IF (FindDet(L11,3) .EQ. 0.0d0) THEN
        !
         cod(:,:,ip) = 0.0d0                       
        !
         cod_upper(:,:,ip) = 0                       
        !
         cod_lower(:,:,ip) = 0                        
        !
        ELSE
         inv_L11 = imat(L11)
         cod(:,:,ip) = L11                         
         etc(:,:,ip) = L22 - MATMUL(MATMUL(L21,inv_L11_tot),L12)  
         seb(:,:,ip) = MATMUL(inv_L11_tot,L12)* 1.0d+6        ! [uV/K]
         seb_x(:,:,ip) = MATMUL(inv_L11_tot,seb_x(:,:,ip))* 1.0d+6        ! [uV/K]
         seb_l(:,:,ip) = MATMUL(inv_L11_tot,seb_l(:,:,ip))* 1.0d+6        ! [uV/K]
         seb_g(:,:,ip) = MATMUL(inv_L11_tot,seb_g(:,:,ip))* 1.0d+6        ! [uV/K]
         pof(:,:,ip) = MATMUL(MATMUL(seb(:,:,ip),seb(:,:,ip)),cod(:,:,ip))
         IF (phdrag) spd_e(:,:,ip) = MATMUL(inv_L11_tot,L44_e) * 1.0d+6        ! [uV/K]
         IF (phdrag) spd_x(:,:,ip) = MATMUL(inv_L11_tot,spd_x(:,:,ip)) * 1.0d+6        ! [uV/K]
	 IF (phdrag) spd_l(:,:,ip) = MATMUL(inv_L11_tot,spd_l(:,:,ip)) * 1.0d+6        ! [uV/K]
	 IF (phdrag) spd_g(:,:,ip) = MATMUL(inv_L11_tot,spd_g(:,:,ip)) * 1.0d+6        ! [uV/K]
        !
         cod_upper(:,:,ip) = L11_upper                         
         etc_upper(:,:,ip) = L22_upper - MATMUL(MATMUL(L21_upper,inv_L11_tot),L12_upper)  
         seb_upper(:,:,ip) = MATMUL(inv_L11_tot,L12_upper)
         pof_upper(:,:,ip) = MATMUL(MATMUL(seb_upper(:,:,ip),seb_upper(:,:,ip)),cod_upper(:,:,ip))
        !
         cod_lower(:,:,ip) = L11_lower                         
         etc_lower(:,:,ip) = L22_lower - MATMUL(MATMUL(L21_lower,inv_L11_tot),L12_lower)  
         seb_lower(:,:,ip) = MATMUL(inv_L11_tot,L12_lower)
         pof_lower(:,:,ip) = MATMUL(MATMUL(seb_lower(:,:,ip),seb_lower(:,:,ip)),cod_lower(:,:,ip))
        !
        ENDIF
     ENDDO  ! ip
     !
#ifdef __PARA
     CALL mp_barrier (inter_pool_comm)
     CALL mp_sum (mfp,inter_pool_comm)
     CALL mp_sum (cod,inter_pool_comm)
     CALL mp_sum (etc,inter_pool_comm)
     CALL mp_sum (seb,inter_pool_comm)
     CALL mp_sum (seb_x,inter_pool_comm)
     CALL mp_sum (seb_l,inter_pool_comm)
     CALL mp_sum (seb_g,inter_pool_comm)
     CALL mp_sum (pof,inter_pool_comm)
     IF (phdrag) CALL mp_sum (spd_e,inter_pool_comm)
     IF (phdrag) CALL mp_sum (spd_x,inter_pool_comm)
     IF (phdrag) CALL mp_sum (spd_l,inter_pool_comm)
     IF (phdrag) CALL mp_sum (spd_g,inter_pool_comm)
     CALL mp_sum (cod_upper,inter_pool_comm)
     CALL mp_sum (etc_upper,inter_pool_comm)
     CALL mp_sum (seb_upper,inter_pool_comm)
     CALL mp_sum (pof_upper,inter_pool_comm)
     CALL mp_sum (cod_lower,inter_pool_comm)
     CALL mp_sum (etc_lower,inter_pool_comm)
     CALL mp_sum (seb_lower,inter_pool_comm)
     CALL mp_sum (pof_lower,inter_pool_comm)
#ENDIF
     !
     CLOSE (8888)
     CLOSE (7777)
     CLOSE (6666)
     IF (phdrag) CLOSE (33333)
     !
     !
     IF (phdrag) THEN
      !
      ! cumulation on phonon MFP
      !
      CALL para_bounds (ip_star, ip_stop, mfp_num)
      !
      OPEN (1221,FILE=phdrag_ufmt,FORM='unformatted',ACCESS='direct',RECL=(1+2*9)*DP,STATUS='old')
      !
      DO ip = ip_star, ip_stop
        !
        mfp_ph(ip) = mfp_ph_min + DBLE(ip-1)*mfp_ph_del
        !
        L33 = 0.0d0
        L44 = 0.0d0
        L33_m = 0.0d0
        L44_m = 0.0d0         
        !
        DO imode = 1, nmodes
           DO iq = 1, nq_ful_red
              !
              IF (mfp_q_ful(imode,iq) .LE. mfp_ph(ip)) THEN
                 READ (1221,REC=(iq-1)*nmodes+imode) wqq, L44_q(1:3,1:3), L33_q(1:3,1:3)
              ELSE
                 !
                 L33_q = 0.0d0
                 L44_q = 0.0d0
                 !  
              ENDIF
                 L33 = L33 + L33_q
                 L44 = L44 + L44_q
                 L33_m(:,:,imode) = L33_m(:,:,imode) + L33_q
                 L44_m(:,:,imode) = L44_m(:,:,imode) + L44_q
            ENDDO ! iq
            !
        ENDDO  ! imode      
           !
           IF (epdim .EQ. 2) THEN
           !
           L33(3,3)   = 1.0d0
           L44(3,3)   = 1.0d0
           L33_m(3,3,:)   = 1.0d0
           L44_m(3,3,:)   = 1.0d0
           !
           ENDIF
           !
           ltc(:,:,ip) = L33                                     ! [W/m/K]
           spd(:,:,ip) = MATMUL(inv_L11_tot,L44) * 1.0d+6        ! [uV/K]
           DO imode = 1, nmodes
            ltc_m(:,:,imode,ip) = L33_m(:,:,imode)                                     ! [W/m/K]
            spd_m(:,:,imode,ip) = MATMUL(inv_L11_tot,L44_m(:,:,imode)) * 1.0d+6        ! [uV/K]
           ENDDO  
           !
     ENDDO ! MFP(ip)
     !
#ifdef __PARA
     CALL mp_barrier (inter_pool_comm)
     CALL mp_sum (mfp_ph,inter_pool_comm)
     CALL mp_sum (ltc,inter_pool_comm)
     CALL mp_sum (spd,inter_pool_comm)
     CALL mp_sum (ltc_m,inter_pool_comm)
     CALL mp_sum (spd_m,inter_pool_comm)
#ENDIF 
     !
     CLOSE (1221)
     !
     CALL para_bounds (ip_star, ip_stop, mfp_num)
     ! cumulation on phonon frequency
     !
     OPEN (2332,FILE=phdrag_ufmt,FORM='unformatted',ACCESS='direct',RECL=(1+2*9)*DP,STATUS='old')
     !
      DO ip = ip_star, ip_stop
        !
        freq_ph(ip) = DBLE(ip-1)*freq_del
        !
        L33 = 0.0d0
        L44 = 0.0d0
        L33_m = 0.0d0
        L44_m = 0.0d0 
        !
        DO imode = 1, nmodes
           DO iq = 1, nq_ful_red
           iq_irr_red = rful2rirr_q(iq)
              !
              IF (wf_all(imode,iq_irr_red) .LE. freq_ph(ip)) THEN
                 READ (2332,REC=(iq-1)*nmodes+imode) wqq, L44_q(1:3,1:3), L33_q(1:3,1:3)
              ELSE
                 !
                 L33_q = 0.0d0
                 L44_q = 0.0d0
                 !  
              ENDIF
                 L33 = L33 + L33_q
                 L44 = L44 + L44_q
                 L33_m(:,:,imode) = L33_m(:,:,imode) + L33_q
                 L44_m(:,:,imode) = L44_m(:,:,imode) + L44_q
           ENDDO  ! iq
        ENDDO  ! imode     
           !
           IF (epdim .EQ. 2) THEN
           !
           L33(3,3)   = 1.0d0
           L44(3,3)   = 1.0d0
           L33_m(3,3,:)   = 1.0d0
           L44_m(3,3,:)   = 1.0d0
           !
           ENDIF
           !
           ltc_f(:,:,ip) = L33                                     ! [W/m/K]
           spd_f(:,:,ip) = MATMUL(inv_L11_tot,L44) * 1.0d+6        ! [uV/K]
           !
           DO imode = 1, nmodes
            ltc_f_m(:,:,imode,ip) = L33_m(:,:,imode)                                     ! [W/m/K]
            spd_f_m(:,:,imode,ip) = MATMUL(inv_L11_tot,L44_m(:,:,imode)) * 1.0d+6        ! [uV/K]
           ENDDO 
     ENDDO ! freq_ph(ip)
     !
#ifdef __PARA
     CALL mp_barrier (inter_pool_comm)
     CALL mp_sum (freq_ph,inter_pool_comm)
     CALL mp_sum (ltc_f,inter_pool_comm)
     CALL mp_sum (spd_f,inter_pool_comm)
     CALL mp_sum (ltc_f_m,inter_pool_comm)
     CALL mp_sum (spd_f_m,inter_pool_comm)
#ENDIF   
     !
     CLOSE (2332) 
     !
     ENDIF  ! phdrag
     !
     ! 4. normalization
     !
     ! 5. export the cumulative cod, etc, seb and pof
     IF (my_pool_id .EQ. ionode_id) THEN
        !
        OPEN (9999,FILE=cumu_fmt,STATUS='replace')
        IF (phdrag) OPEN (7223,FILE=phdragcumu_MFP_fmt,STATUS='replace')
        IF (phdrag)  OPEN (8992,FILE=phdragcumu_freq_fmt,STATUS='replace')
        IF (phdrag) OPEN (4455,FILE=phdragcumu_k_MFP_fmt,STATUS='replace')
        OPEN (9191,FILE=cumu_upper_fmt,STATUS='replace')
        OPEN (9393,FILE=cumu_lower_fmt,STATUS='replace')
        OPEN (2224,FILE='BTE/'//TRIM(ADJUSTL(txnx))//'/Seebeck_k_MFPcum_XLG.dat',STATUS='replace')
     !
        DO ip = 1, mfp_num
           !
           IF (epdim .EQ. 3) THEN
               !
               WRITE (9999,'(es14.4,12es13.5)') mfp(ip)*au2nm, cod(1,1,ip), etc(1,1,ip), seb(1,1,ip), pof(1,1,ip)*1.0d+4, &
                                                               cod(2,2,ip), etc(2,2,ip), seb(2,2,ip), pof(2,2,ip)*1.0d+4, &
                                                               cod(3,3,ip), etc(3,3,ip), seb(3,3,ip), pof(3,3,ip)*1.0d+4
               !
               WRITE (9191,'(es14.4,6es13.5)') mfp(ip)*au2nm, cod_upper(1,1,ip), seb_upper(1,1,ip)*1.0d+6, &
                                                              cod_upper(2,2,ip), seb_upper(2,2,ip)*1.0d+6, &
                                                              cod_upper(3,3,ip), seb_upper(3,3,ip)*1.0d+6
               !
               WRITE (9393,'(es14.4,6es13.5)') mfp(ip)*au2nm, cod_lower(1,1,ip), seb_lower(1,1,ip)*1.0d+6, &
                                                              cod_lower(2,2,ip), seb_lower(2,2,ip)*1.0d+6, &
                                                              cod_lower(3,3,ip), seb_lower(3,3,ip)*1.0d+6
               !
               WRITE (2224,'(es14.4,4es13.5)') mfp(ip)*au2nm, (seb(1,1,ip)+seb(2,2,ip)+seb(3,3,ip))/3.0d0, &
               (seb_x(1,1,ip)+seb_x(2,2,ip)+seb_x(3,3,ip))/3.0d0,(seb_l(1,1,ip)+seb_l(2,2,ip)+seb_l(3,3,ip))/3.0d0, &
               (seb_g(1,1,ip)+seb_g(2,2,ip)+seb_g(3,3,ip))/3.0d0
               !
               IF (phdrag) WRITE (4455,'(es14.4,4es13.5)') mfp(ip)*au2nm, &
               (spd_e(1,1,ip)+spd_e(2,2,ip)+spd_e(3,3,ip))/3.0d0, (spd_x(1,1,ip)+spd_x(2,2,ip)+spd_x(3,3,ip))/3.0d0, &
               (spd_l(1,1,ip)+spd_l(2,2,ip)+spd_l(3,3,ip))/3.0d0, (spd_g(1,1,ip)+spd_g(2,2,ip)+spd_g(3,3,ip))/3.0d0
               IF (phdrag) WRITE (7223,'(es14.4,18es13.5)') mfp_ph(ip)*au2nm, spd(1,1,ip), ltc(1,1,ip),spd(2,2,ip), ltc(2,2,ip),&
                                                                 spd(3,3,ip), ltc(3,3,ip),spd_m(1,1,1:6,ip),&
                                                                 ltc_m(1,1,1:6,ip)
               IF (phdrag) WRITE (8992,'(es14.4,18es13.5)') freq_ph(ip)*ryd2thz, spd_f(1,1,ip), ltc_f(1,1,ip), spd_f(2,2,ip), &
                                                                 ltc_f(2,2,ip), spd_f(3,3,ip), ltc_f(3,3,ip), &
                                                                 spd_f_m(1,1,1:6,ip), ltc_f_m(1,1,1:6,ip)
               !
           ELSEIF (epdim .EQ. 2) THEN
               !
               WRITE (9999,'(es14.4,8es13.5)') mfp(ip)*au2nm, cod(1,1,ip), etc(1,1,ip), seb(1,1,ip), pof(1,1,ip)*1.0d+4, &
                                                              cod(2,2,ip), etc(2,2,ip), seb(2,2,ip), pof(2,2,ip)*1.0d+4
               !
               WRITE (9191,'(es14.4,4es13.5)') mfp(ip)*au2nm, cod_upper(1,1,ip), seb_upper(1,1,ip)*1.0d+6, &
                                                              cod_upper(2,2,ip), seb_upper(2,2,ip)*1.0d+6
               !
               WRITE (9393,'(es14.4,4es13.5)') mfp(ip)*au2nm, cod_lower(1,1,ip), seb_lower(1,1,ip)*1.0d+6, &
                                                              cod_lower(2,2,ip), seb_lower(2,2,ip)*1.0d+6
               IF (phdrag) WRITE (7223,'(es14.4,16es13.5)') mfp_ph(ip)*au2nm, spd(1,1,ip), ltc(1,1,ip),spd(2,2,ip), ltc(2,2,ip),&
                                                                 spd_m(1,1,1:6,ip), ltc_m(1,1,1:6,ip)
               IF (phdrag) WRITE (8992,'(es14.4,16es13.5)') freq_ph(ip)*ryd2thz, spd_f(1,1,ip), ltc_f(1,1,ip), spd_f(2,2,ip), &
                                                                 ltc_f(2,2,ip), spd_f_m(1,1,1:6,ip), &
                                                                 ltc_f_m(1,1,1:6,ip)              
               !
           ENDIF
           !
        ENDDO
        CLOSE (9999)
        CLOSE (9191)
        CLOSE (9393)
        CLOSE (2224)
        IF (phdrag) CLOSE (7223)
        IF (phdrag) CLOSE (8992)
        IF (phdrag) CLOSE (4455)
        !
        ! output each q phdrag contribution, el-ph scattering rate, MFP
           !
           !
           rate_RTA_3Dphdrag_fmt = 'BTE/'//TRIM(ADJUSTL(txnx))//'/rful_q_phdrag_3D.dat'
           !
           !
           IF (phdrag) THEN
           !
           OPEN (2332,FILE=phdrag_ufmt,FORM='unformatted',ACCESS='direct',RECL=(1+2*9)*DP,STATUS='old')
           OPEN (8788,FILE='BTE/META/xqf_fbz_cart',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='old')
           OPEN (9899,FILE=rate_RTA_3Dphdrag_fmt,STATUS='replace')
           
           DO iq = 1, nq_ful_red
            iq_ful = rful2ful_q(iq)
            READ (8788,REC=iq_ful) xq_fbz(1:3)
           !
            length = SQRT(DOT_PRODUCT(xq_fbz,xq_fbz))
           !
            DO imode = 1, nmodes
             READ (2332,REC=(iq-1)*nmodes+imode) wqq, L44_q(1:3,1:3), L33_q(1:3,1:3)
             L44_q = MATMUL(inv_L11_tot,L44_q)*1.0d+6  ! = spd_q 
           !
             WRITE (9899,'(i8,f12.6,f11.4,3es14.4)') iq_irr_red, length, wqq*ryd2thz, mfp_q_ful(imode,iq)*au2nm, &
                   				     (L44_q(1,1)+L44_q(2,2)+L44_q(3,3))/3.0d0, L33_q(1,1) 
           !
           ! [THz]
           ! L33_q                                     ! [W/m/K]
           ! MATMUL(inv_L11_tot,L44_q) * 1.0d+6        ! [uV/K]
           !
            ENDDO  ! nmodes
           ENDDO  ! nq_ful_red	   
           !
           CLOSE (2332) 
           CLOSE (8788)
           CLOSE (9899)
           !
           ! output each k phdrag contribution, el-ph scattering rate, MFP
           !
           !
           rate_k_RTA_3Dphdrag_fmt = 'BTE/'//TRIM(ADJUSTL(txnx))//'/rful_k_phdrag_3D.dat'
           !
           !
           OPEN (44444,FILE=phdrag_k_ufmt,FORM='unformatted',ACCESS='direct',RECL=(1+9)*DP,STATUS='old')
           OPEN (99999,FILE='BTE/META/xkf_fbz_cart',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='old')
           OPEN (66666,FILE=rate_k_RTA_3Dphdrag_fmt,STATUS='replace')
           
           DO ik_ful_red = 1, nk_ful_red
            ik_ful = rful2ful(ik_ful_red)
            ! ik = ful2rful(ik_ful) ! nk_ful_red
            READ (99999,REC=ik_ful) xk_fbz(1:3)
           !
            length = SQRT(DOT_PRODUCT(xk_fbz,xk_fbz))
           !
            DO ibnd = 1, nbnd_red
           !
           !sigmai = 2.0d0 * ( SUM(sigmai_mode_all_abs(itemp,idope,1:nmodes,ibnd,ik_irr_red)) + SUM(sigmai_mode_all_emi(itemp,idope,1:nmodes,ibnd,ik_irr_red)) )
           !
               READ (44444,REC=(ik_ful_red-1)*nbnd_red+ibnd) ekk, L44_k(1:3,1:3)
               L44_k = MATMUL(inv_L11_tot,L44_k)*1.0d+6  ! = spd_k
           !
               WRITE (66666,'(i8,f12.6,f11.4,2es14.4)') ik_ful_red, length, ekk*ryd2ev, &
                                                        mfp_k_ful(ibnd,ik_ful_red)*au2nm, L44_k(1,1)
           !
           ! [THz] [Thz]
           ! Output L44_k = (inv_L11,L44_k) * 1.0d+6        ! [uV/K]
           !
            ENDDO  ! nbnd_red
           ENDDO  ! nk_irr_red	   
           !

           CLOSE (44444) 
           CLOSE (99999)
           CLOSE (66666)
           !
           ENDIF !phdrag
     ENDIF ! mypool
     !
     !
     DEALLOCATE (mfp)
     IF (phdrag) DEALLOCATE (mfp_ph)
     IF (phdrag) DEALLOCATE (freq_ph)
     DEALLOCATE (cod)
     DEALLOCATE (etc)
     DEALLOCATE (seb)
     DEALLOCATE (seb_x)
     DEALLOCATE (seb_l)
     DEALLOCATE (seb_g)
     DEALLOCATE (pof)
     IF (phdrag) DEALLOCATE (ltc)
     IF (phdrag) DEALLOCATE (spd)
     IF (phdrag) DEALLOCATE (spd_e)
     IF (phdrag) DEALLOCATE (spd_x)
     IF (phdrag) DEALLOCATE (spd_l)
     IF (phdrag) DEALLOCATE (spd_g)
     IF (phdrag) DEALLOCATE (ltc_f)
     IF (phdrag) DEALLOCATE (spd_f)
     IF (phdrag) DEALLOCATE (ltc_m)
     IF (phdrag) DEALLOCATE (spd_m)
     IF (phdrag) DEALLOCATE (ltc_f_m)
     IF (phdrag) DEALLOCATE (spd_f_m)
     DEALLOCATE (cod_upper)
     DEALLOCATE (etc_upper)
     DEALLOCATE (seb_upper)
     DEALLOCATE (pof_upper)
     DEALLOCATE (cod_lower)
     DEALLOCATE (etc_lower)
     DEALLOCATE (seb_lower)
     DEALLOCATE (pof_lower)
     !
  ENDIF 
  !
END SUBROUTINE export_cumu



!----------------------------------------------------------------------------
SUBROUTINE export_rate (iter, converge, itemp, idope)
!----------------------------------------------------------------------------
! 
! export el-ph electron scattering rate of each mode under RTA, and al-el scattering rate
!
! ik (irr), ibnd (red), energy [eV], scattering rate [THz]
!
#INCLUDE "f_defs.h"
  USE kinds,         ONLY : DP
  USE phcom,         ONLY : nmodes
  USE epwcom,        ONLY : nbndsub, bte, alloy_pot
  USE elph2,         ONLY : ibndmin, ibndmax, nbnd_red, &
                            etf_all, sigmai_mode_all_abs, sigmai_mode_all_emi, sigmai_mode_all_alloy_inter, &
                            sigmai_mode_all_alloy_intra, mfp_k_ful, vel_ful
  USE bte_var
  USE constants_epw, ONLY : ryd2ev, ryd2thz, twopi, au2m, au2s
#ifdef __PARA
  USE io_global,     ONLY : ionode_id
  USE mp_global,     ONLY : my_pool_id
#ENDIF
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iter, itemp, idope
  LOGICAL, INTENT(IN) :: converge
  REAL(KIND=DP)       :: sigmai, v_abs
  INTEGER             :: ik, ik_ful, ik_ful_red, ik_irr, ibnd, ibnd0, imode
  CHARACTER(LEN=256)  :: rate_RTA_fmt, rate_mode_RTA_fmt, rate_EFF_fmt, velo_fmt
  CHARACTER(LEN=12)   :: txnx
  CHARACTER(LEN=3)    :: itemp_num, idope_num
  !
  !
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     IF (iter .EQ. 0) THEN
        !
        WRITE(itemp_num,'(i3)') itemp
        WRITE(idope_num,'(i3)') idope
        txnx = 'T'//TRIM(ADJUSTL(itemp_num))//'_N'//TRIM(ADJUSTL(idope_num))
        !rate_RTA_fmt = 'BTE/'//TRIM(ADJUSTL(txnx))//'/rate_RTA_'//TRIM(ADJUSTL(txnx))//'.dat'
        !rate_mode_RTA_fmt = 'BTE/'//TRIM(ADJUSTL(txnx))//'/rate_mode_RTA_'//TRIM(ADJUSTL(txnx))//'.dat'
        rate_RTA_fmt = 'BTE/'//TRIM(ADJUSTL(txnx))//'/rate_RTA.dat'
        rate_mode_RTA_fmt = 'BTE/'//TRIM(ADJUSTL(txnx))//'/rate_mode_RTA.dat'
        velo_fmt = 'BTE/'//TRIM(ADJUSTL(txnx))//'/velo.dat'
        OPEN (8888,FILE=rate_RTA_fmt,STATUS='replace')
        OPEN (9999,FILE=rate_mode_RTA_fmt,STATUS='replace')
        OPEN (7777,FILE=velo_fmt,STATUS='replace')
        !
        DO ik = 1, nk_irr_red
           !
           IF (bte .EQ. 19 .OR. bte .EQ. 18) THEN
              ik_ful = ik
              ik_ful_red = ik
           ELSE
              ik_ful = irr2ful(rirr2irr(ik))
              ik_ful_red = ful2rful(ik_ful)
           ENDIF           
           !
           DO ibnd = 1, nbnd_red
              !
              ibnd0 = ibnd+ibndmin-1
              !
              IF (alloy_pot) THEN
              sigmai = 2.0d0 * ( SUM(sigmai_mode_all_abs(itemp,idope,1:nmodes,ibnd,ik)) + SUM(sigmai_mode_all_emi(itemp,idope,1:nmodes,ibnd,ik)) + sigmai_mode_all_alloy_inter(itemp,idope,ibnd,ik) + sigmai_mode_all_alloy_intra(itemp,idope,ibnd,ik))
              ELSE
              sigmai = 2.0d0 * ( SUM(sigmai_mode_all_abs(itemp,idope,1:nmodes,ibnd,ik)) + SUM(sigmai_mode_all_emi(itemp,idope,1:nmodes,ibnd,ik)))
              ENDIF
              !
              IF (sigmai .NE. 0.0d0) WRITE (8888,'(i12,i5,f10.4,es14.4)') ik_ful, ibnd0, etf_all(ibnd0,ik)*ryd2ev, sigmai*ryd2thz*twopi   
              !
              IF (sigmai .NE. 0.0d0) THEN
               IF (alloy_pot) THEN
WRITE (9999,'(i12,i5,f10.4,96es14.4)') ik_ful, ibnd0, etf_all(ibnd0,ik)*ryd2ev, &
2.0d0*sigmai_mode_all_abs(itemp,idope,1:nmodes,ibnd,ik)*ryd2thz*twopi, & ! ryd2thz*twopi = 1/au2ps
2.0d0*sigmai_mode_all_emi(itemp,idope,1:nmodes,ibnd,ik)*ryd2thz*twopi, &    ! 2*sigmai is in units of [1/au_s]
2.0d0*sigmai_mode_all_alloy_intra(itemp,idope,ibnd,ik)*ryd2thz*twopi, &
2.0d0*sigmai_mode_all_alloy_inter(itemp,idope,ibnd,ik)*ryd2thz*twopi
               ELSE
WRITE (9999,'(i12,i5,f10.4,96es14.4)') ik_ful, ibnd0, etf_all(ibnd0,ik)*ryd2ev, &
2.0d0*sigmai_mode_all_abs(itemp,idope,1:nmodes,ibnd,ik)*ryd2thz*twopi, & ! ryd2thz*twopi = 1/au2ps
2.0d0*sigmai_mode_all_emi(itemp,idope,1:nmodes,ibnd,ik)*ryd2thz*twopi   ! 2*sigmai is in units of [1/au_s]
               ENDIF
              ENDIF
              !
              IF (sigmai .NE. 0.0d0) WRITE (7777,'(i12,i5,f10.4,3es16.6)') ik_ful, ibnd0, etf_all(ibnd0,ik)*ryd2ev, &
                                     vel_ful(1:3,ibnd,ik_ful_red)*(au2m/au2s)
        
              !
           ENDDO
        ENDDO
        CLOSE (8888)
        CLOSE (9999)
        CLOSE (7777)
        !
     ELSEIF (iter .GT. 0 .AND. converge .EQ. .TRUE.) THEN
        !
        WRITE(itemp_num,'(i3)') itemp
        WRITE(idope_num,'(i3)') idope
        txnx = 'T'//TRIM(ADJUSTL(itemp_num))//'_N'//TRIM(ADJUSTL(idope_num))
        !rate_EFF_fmt = 'BTE/'//TRIM(ADJUSTL(txnx))//'/rate_EFF_'//TRIM(ADJUSTL(txnx))//'.dat'
        rate_EFF_fmt = 'BTE/'//TRIM(ADJUSTL(txnx))//'/rate_EFF.dat'
        OPEN (9999,FILE=rate_EFF_fmt,STATUS='replace')
        !
        DO ik = 1, nk_ful_red
           !
           ik_ful = rful2ful(ik)
           ik_irr = rful2rirr(ik)
           !
           DO ibnd = 1, nbnd_red
              !
              ibnd0 = ibnd+ibndmin-1
              !
              v_abs = SQRT(DOT_PRODUCT(vel_ful(:,ibnd,ik),vel_ful(:,ibnd,ik)))
              !
              IF (mfp_k_ful(ibnd,ik) .NE. 0.0d0) WRITE (9999,'(i12,i5,f10.4,es14.4)') ik_ful, ibnd0, etf_all(ibnd0,ik_irr)*ryd2ev, &
                                                       (v_abs/mfp_k_ful(ibnd,ik))*ryd2thz*twopi                                 
              !
           ENDDO
        ENDDO
        !
        CLOSE (9999)
        !
     ENDIF
     !
  ENDIF
  !
  !
END SUBROUTINE export_rate



!----------------------------------------------------------------------------
SUBROUTINE export_mfp (iter, converge, itemp, idope)
!----------------------------------------------------------------------------
! 
! export electron MFP
! 
! ik (irr), ibnd (red), energy [eV], MFP [nm]
!
#INCLUDE "f_defs.h"
  USE kinds,         ONLY : DP
  USE phcom,         ONLY : nmodes
  USE epwcom,        ONLY : nbndsub, bte
  USE elph2,         ONLY : ibndmin, ibndmax, nbnd_red, etf_all, mfp_k_ful
  USE bte_var
  USE constants_epw, ONLY : ryd2ev, au2nm
#ifdef __PARA
  USE io_global,     ONLY : ionode_id
  USE mp_global,     ONLY : my_pool_id
#ENDIF
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iter, itemp, idope
  LOGICAL, INTENT(IN) :: converge
  REAL(KIND=DP)       :: v_abs
  INTEGER             :: ik, ik_ful, ik_irr, ibnd, ibnd0, imode
  CHARACTER(LEN=256)  :: mfp_RTA_fmt, mfp_EFF_fmt
  CHARACTER(LEN=12)   :: txnx
  CHARACTER(LEN=3)    :: itemp_num, idope_num
  !
  !
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     IF ( (iter .EQ. 0) .OR. (iter .GT. 0 .AND. converge .EQ. .TRUE.) ) THEN
        !
        IF (iter .EQ. 0) THEN
           !
           WRITE(itemp_num,'(i3)') itemp
           WRITE(idope_num,'(i3)') idope
           txnx = 'T'//TRIM(ADJUSTL(itemp_num))//'_N'//TRIM(ADJUSTL(idope_num))
           !mfp_RTA_fmt = 'BTE/'//TRIM(ADJUSTL(txnx))//'/mfp_RTA_'//TRIM(ADJUSTL(txnx))//'.dat'
           mfp_RTA_fmt = 'BTE/'//TRIM(ADJUSTL(txnx))//'/mfp_RTA.dat'
           OPEN (9999,FILE=mfp_RTA_fmt,STATUS='replace')
           !
        ENDIF
        !
        IF (iter .GT. 0 .AND. converge .EQ. .TRUE.) THEN
           !
           WRITE(itemp_num,'(i3)') itemp
           WRITE(idope_num,'(i3)') idope
           txnx = 'T'//TRIM(ADJUSTL(itemp_num))//'_N'//TRIM(ADJUSTL(idope_num))
           !mfp_EFF_fmt = 'BTE/'//TRIM(ADJUSTL(txnx))//'/mfp_EFF_'//TRIM(ADJUSTL(txnx))//'.dat'
           mfp_EFF_fmt = 'BTE/'//TRIM(ADJUSTL(txnx))//'/mfp_EFF.dat'
           OPEN (9999,FILE=mfp_EFF_fmt,STATUS='replace')
           !
        ENDIF
        !
        DO ik = 1, nk_ful_red
           !
           IF (bte .EQ. 19 .OR. bte .EQ. 18) THEN
              ik_ful = ik
              ik_irr = ik
           ELSE
              ik_ful = rful2ful(ik)
              ik_irr = rful2rirr(ik)
           ENDIF  
           !
           DO ibnd = 1, nbnd_red
              !
              ibnd0 = ibnd+ibndmin-1
              !
              IF (mfp_k_ful(ibnd,ik) .NE. 0.0d0) WRITE (9999,'(i12,i5,f10.4,es14.4)') ik_ful, ibnd0, etf_all(ibnd0,ik_irr)*ryd2ev, &
                                                                                      mfp_k_ful(ibnd,ik)*au2nm                                  
              !
           ENDDO
        ENDDO
        !
        CLOSE (9999)
        !
     ENDIF
     !
  ENDIF
  !
  !
END SUBROUTINE export_mfp



!----------------------------------------------------------------------------
SUBROUTINE export_boltzep (iter, converge, itemp, idope)
!----------------------------------------------------------------------------
! 
! export electron info
! 
! ik (irr), ibnd (red), energy [eV], scattering rate[THz], velocity[m/s]
!
#INCLUDE "f_defs.h"
  USE kinds,         ONLY : DP
  USE cell_base,     ONLY : omega
  USE phcom,         ONLY : nmodes
  USE spin_orb,      ONLY : lspinorb
  USE epwcom,        ONLY : nkf1, nkf2, nkf3, nbndsub, bte, eptemp
  USE elph2,         ONLY : ibndmin, ibndmax, nbnd_red, etf_all, vel_ful, mfp_k_ful, ef_epw
  USE bte_var
  USE constants_epw, ONLY : ryd2ev, au2m, au2s, ryd2thz, twopi, kB
#ifdef __PARA
  USE io_global,     ONLY : ionode_id
  USE mp_global,     ONLY : my_pool_id
#ENDIF
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iter, itemp, idope
  LOGICAL, INTENT(IN) :: converge
  REAL(KIND=DP)       :: temp, ef0, vol, ekk, vkk(3), sigmai
  INTEGER             :: ik, ik_ful, ik_ful_red, ik_irr_red, ibnd, ibnd0
  CHARACTER(LEN=256)  :: boltzep_fmt
  CHARACTER(LEN=12)   :: txnx
  CHARACTER(LEN=3)    :: itemp_num, idope_num
  !
  !
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     IF (converge .EQ. .TRUE.) THEN
        !
        WRITE(itemp_num,'(i3)') itemp
        WRITE(idope_num,'(i3)') idope
        txnx = 'T'//TRIM(ADJUSTL(itemp_num))//'_N'//TRIM(ADJUSTL(idope_num))
        boltzep_fmt = 'BTE/'//TRIM(ADJUSTL(txnx))//'/BoltzEP_elph.dat'
        OPEN (9999,FILE=boltzep_fmt,STATUS='replace')
        !
        temp = eptemp(itemp)/kB
        vol = omega*au2m*au2m*au2m
        ef0  = ef_epw(itemp,idope) * ryd2ev
        !
        WRITE (9999,'(3i5,i12,i5,3es15.6,l4)') nkf1, nkf2, nkf3, nk_ful_red, nbnd_red, temp, vol, ef0, lspinorb
        !
        DO ik_ful_red = 1, nk_ful_red
           !
           ik_irr_red = rful2rirr(ik_ful_red)
           ik_ful     = rful2ful(ik_ful_red)
           !
           DO ibnd = 1, nbnd_red
              !
              ibnd0 = ibnd+ibndmin-1
              !
              ekk = etf_all(ibnd0,ik_irr_red)
              vkk(:) = vel_ful(:,ibnd,ik_ful_red)
              IF (mfp_k_ful(ibnd,ik_ful_red) .NE. 0.0d0) THEN
                 sigmai = (SQRT(DOT_PRODUCT(vkk,vkk)) / mfp_k_ful(ibnd,ik_ful_red))
              ELSE
                 sigmai = 0.0d0
              ENDIF
              !
              WRITE (9999,'(i12,i5,5es24.15)') ik_ful, ibnd0, ekk*ryd2ev, sigmai*ryd2thz*twopi, vkk(1:3)*au2m/au2s ! [eV], [THz], [m/s]
              !
           ENDDO
           !
        ENDDO
        !
        CLOSE (9999)
        !
     ENDIF
     !
  ENDIF
  !
  !
END SUBROUTINE export_boltzep



!-------------------------------------------------------------------------------
SUBROUTINE export_fprime (iter, converge, itemp, idope)
!-------------------------------------------------------------------------------
! 
! export electron distribution function
! 
! ik (irr), ibnd (red), xkf_fbz [cart], energy [eV], f0 [none], f' [m/V]
!
#INCLUDE "f_defs.h"
  USE kinds,         ONLY : DP
  USE pwcom,         ONLY : ef 
  USE phcom,         ONLY : nmodes
  USE epwcom,        ONLY : nbndsub, eptemp, efermi_read, fermi_energy, bte
  USE elph2,         ONLY : ibndmin, ibndmax, nbnd_red, etf_all, F_k_ful, ef_epw
  USE bte_var
  USE constants_epw, ONLY : ryd2ev, au2m, au2V
#ifdef __PARA
  USE io_global,     ONLY : ionode_id
  USE mp_global,     ONLY : my_pool_id
#ENDIF
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)     :: iter, itemp, idope
  LOGICAL, INTENT(IN)     :: converge
  REAL(KIND=DP)           :: f0, fprime(3), eptemp0, ef0, xk_fbz(3), non_r
  INTEGER                 :: ik, ik_ful, ik_irr, ibnd, ibnd0, imode
  CHARACTER(LEN=256)      :: fprime_RTA_fmt, fprime_EFF_fmt
  CHARACTER(LEN=12)       :: txnx
  CHARACTER(LEN=3)        :: itemp_num, idope_num
  !
  REAL(KIND=DP), EXTERNAL :: wgauss
  !
  !
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     IF ( (iter .EQ. 0) .OR. (iter .GT. 0 .AND. converge .EQ. .TRUE.) ) THEN
        !
        OPEN (8888,FILE='BTE/META/xkf_fbz_cart',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='old')
        !
        IF (iter .EQ. 0) THEN
           !
           WRITE(itemp_num,'(i3)') itemp
           WRITE(idope_num,'(i3)') idope
           txnx = 'T'//TRIM(ADJUSTL(itemp_num))//'_N'//TRIM(ADJUSTL(idope_num))
           !fprime_RTA_fmt = 'BTE/'//TRIM(ADJUSTL(txnx))//'/fprime_RTA_'//TRIM(ADJUSTL(txnx))//'.dat'
           fprime_RTA_fmt = 'BTE/'//TRIM(ADJUSTL(txnx))//'/fprime_RTA.dat'
           OPEN (9999,FILE=fprime_RTA_fmt,STATUS='replace')
           !
        ENDIF
        !
        IF (iter .GT. 0 .AND. converge .EQ. .TRUE.) THEN
           !
           WRITE(itemp_num,'(i3)') itemp
           WRITE(idope_num,'(i3)') idope
           txnx = 'T'//TRIM(ADJUSTL(itemp_num))//'_N'//TRIM(ADJUSTL(idope_num))
           !fprime_EFF_fmt = 'BTE/'//TRIM(ADJUSTL(txnx))//'/fprime_EFF_'//TRIM(ADJUSTL(txnx))//'.dat'
           fprime_EFF_fmt = 'BTE/'//TRIM(ADJUSTL(txnx))//'/fprime_EFF.dat'
           OPEN (9999,FILE=fprime_EFF_fmt,STATUS='replace')
           !
        ENDIF
        !
        eptemp0 = eptemp(itemp)  ! in units of Ry, essentially kBT
        !
        IF (efermi_read) THEN
           ef0 = fermi_energy
        ELSE
           ef0 = ef_epw(itemp,idope)
        ENDIF
        !
        DO ik = 1, nk_ful_red
           !
           ik_ful = rful2ful(ik)
           ik_irr = rful2rirr(ik)
           !
           DO ibnd = 1, nbnd_red
              !
              ibnd0 = ibnd+ibndmin-1
              !
              f0 = wgauss(-(etf_all(ibnd0,ik_irr)-ef0)/eptemp0,-99) 
              fprime(:) = (1.0d0/eptemp0)*f0*(1.0d0-f0)*F_k_ful(:,ibnd,ik)*SQRT(2.0d0)*(au2m/au2V) ! ~[m/V] normalized by electric field
              !
              IF (MAXVAL(ABS(fprime(:))) .NE. 0.0d0) THEN
                 READ (8888,REC=ik_ful) xk_fbz(1:3)
                 WRITE (9999,'(i12,i5,3f12.6,f10.4,4es14.4)') ik_ful, ibnd0, xk_fbz(1:3), etf_all(ibnd0,ik_irr)*ryd2ev, f0, fprime(1:3)
              ENDIF
              !
           ENDDO 
        ENDDO 
        !
        CLOSE (8888)
        CLOSE (9999)
        !
     ENDIF
     !
  ENDIF
  !
  !
END SUBROUTINE export_fprime



!-------------------------------------------------------------------------------
SUBROUTINE export_onsager (iter, converge, itemp, idope)
!-------------------------------------------------------------------------------
!
! export electron transport properties and Onsager coefficients
!
! L11 [1/Ohm/m], L12 [A/m/K], L21 [A/m], L22 [W/m/K]
!
#INCLUDE "f_defs.h"
  USE kinds,        ONLY : DP
  USE epwcom,       ONLY : nbndsub
  USE elph2,        ONLY : nbnd_red
  USE bte_var
#ifdef __PARA
  USE io_global,    ONLY : ionode_id
  USE mp_global,    ONLY : me_pool, my_pool_id
#ENDIF
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iter, itemp, idope
  LOGICAL, INTENT(IN) :: converge
  REAL(KIND=DP)       :: L11(3,3), L12(3,3), L21(3,3), L22(3,3), &
                         L11_k(3,3), L12_k(3,3), L21_k(3,3), L22_k(3,3)
  INTEGER             :: ik, ibnd, i, j, uorl
  CHARACTER(LEN=256)  :: onsager_ufmt, onsager_fmt
  CHARACTER(LEN=12)   :: txnx
  CHARACTER(LEN=3)    :: itemp_num, idope_num
  !
  !
  IF (converge .EQ. .TRUE. .AND. my_pool_id .EQ. ionode_id) THEN
     !
     WRITE(itemp_num,'(i3)') itemp
     WRITE(idope_num,'(i3)') idope
     txnx = 'T'//TRIM(ADJUSTL(itemp_num))//'_N'//TRIM(ADJUSTL(idope_num))
     onsager_ufmt = 'BTE/META/onsager_'//TRIM(ADJUSTL(txnx))
     !onsager_fmt = 'BTE/'//TRIM(ADJUSTL(txnx))//'/onsager_'//TRIM(ADJUSTL(txnx))//'.dat'
     onsager_fmt = 'BTE/'//TRIM(ADJUSTL(txnx))//'/onsager.dat'
     !
     L11   = 0.0d0
     L12   = 0.0d0
     L21   = 0.0d0
     L22   = 0.0d0
     !
     OPEN (8888,FILE=onsager_ufmt,FORM='unformatted',ACCESS='direct',RECL=4+9*4*DP,STATUS='old')
     DO ik = 1, nk_ful_red
        DO ibnd = 1, nbnd_red
           !
           READ (8888,REC=(ik-1)*nbnd_red+ibnd) L11_k(1:3,1:3), L12_k(1:3,1:3), L21_k(1:3,1:3), L22_k(1:3,1:3), uorl
           !
           L11 = L11 + L11_k
           L12 = L12 + L12_k
           L21 = L21 + L21_k
           L22 = L22 + L22_k
           !
        ENDDO
     ENDDO
     CLOSE (8888)
     !
     OPEN (9999,FILE=onsager_fmt,STATUS='replace')
     WRITE (9999,'(3x,a)') 'L11'
     WRITE (9999,'(3es15.6)') L11(1,1), L11(1,2), L11(1,3)
     WRITE (9999,'(3es15.6)') L11(2,1), L11(2,2), L11(2,3)
     WRITE (9999,'(3es15.6)') L11(3,1), L11(3,2), L11(3,3)
     WRITE (9999,'(3x,a)') 'L12'
     WRITE (9999,'(3es15.6)') L12(1,1), L12(1,2), L12(1,3)
     WRITE (9999,'(3es15.6)') L12(2,1), L12(2,2), L12(2,3)
     WRITE (9999,'(3es15.6)') L12(3,1), L12(3,2), L12(3,3)
     WRITE (9999,'(3x,a)') 'L21'
     WRITE (9999,'(3es15.6)') L21(1,1), L21(1,2), L21(1,3)
     WRITE (9999,'(3es15.6)') L21(2,1), L21(2,2), L21(2,3)
     WRITE (9999,'(3es15.6)') L21(3,1), L21(3,2), L21(3,3)
     WRITE (9999,'(3x,a)') 'L22'
     WRITE (9999,'(3es15.6)') L22(1,1), L22(1,2), L22(1,3)
     WRITE (9999,'(3es15.6)') L22(2,1), L22(2,2), L22(2,3)
     WRITE (9999,'(3es15.6)') L22(3,1), L22(3,2), L22(3,3)
     CLOSE (9999)
     !
  ENDIF
  !
  !
END SUBROUTINE export_onsager



!----------------------------------------------------------------------------
SUBROUTINE export_rate_ph ()
!----------------------------------------------------------------------------
! 
! export phonon scattering rate of each mode under RTA
!
! iq_irr, frequency 1~nmodes [cm^-1], scattering rate 1~nmodes [THz]
!
#INCLUDE "f_defs.h"
  USE kinds,         ONLY : DP
  USE phcom,         ONLY : nmodes
  USE epwcom,        ONLY : bte, neptemp, nepdope, nqf1, nqf2, nqf3, phdrag
  USE elph2,         ONLY : gammai_mode_all, wf_irr, vph_irr, wf_all
  USE bte_var
  USE bte_func
  USE constants_epw, ONLY : ryd2ev, rydcm1, ryd2thz, twopi, au2m, au2s
#ifdef __PARA
  USE io_global,     ONLY : ionode_id, stdout, ionode
  USE mp,            ONLY : mp_barrier, mp_bcast, mp_sum
  USE mp_global,     ONLY : my_pool_id, nproc_pool, intra_image_comm, &
                            inter_pool_comm, me_pool, root_pool, intra_pool_comm, npool
#ENDIF
  !
  IMPLICIT NONE
  !
  INTEGER             :: iq, iq_irr, imode, itemp, idope, iq_ful, nq_tot
  REAL(KIND=DP)       :: xq_fbz(3), length, wf(nmodes)
  
  CHARACTER(LEN=256)  :: rate_RTA_fmt, rate_RTA_3D_fmt
  CHARACTER(LEN=12)   :: txnx
  CHARACTER(LEN=3)    :: itemp_num, idope_num
  !
  !
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     DO itemp = 1, neptemp
        DO idope = 1, nepdope
           !
           WRITE(itemp_num,'(i3)') itemp
           WRITE(idope_num,'(i3)') idope
           txnx = 'T'//TRIM(ADJUSTL(itemp_num))//'_N'//TRIM(ADJUSTL(idope_num))
           rate_RTA_fmt = 'BTE/'//TRIM(ADJUSTL(txnx))//'/rate_RTA_phe.dat'
           rate_RTA_3D_fmt = 'BTE/'//TRIM(ADJUSTL(txnx))//'/rate_RTA_phe_3D.dat'
           !
           OPEN (1111,FILE=rate_RTA_fmt,STATUS='replace')

           IF (phdrag) THEN
              nq_tot = nq_irr_red
           ELSE
              nq_tot = nq_irr
           ENDIF
           !
           DO iq = 1, nq_tot
              !
              IF (phdrag) THEN
                 wf(:) = wf_all(:,iq)
              ELSE
                 wf(:) = wf_irr(:,iq)
              ENDIF
              !
              WRITE (1111,'(i8,96es14.4)') iq, wf(1:nmodes)*rydcm1, 2.0d0*gammai_mode_all(itemp,idope,1:nmodes,iq)*ryd2thz*twopi            
           ENDDO
           CLOSE (1111)
           !
           !
              OPEN (7777,FILE='BTE/META/irr2ful_q',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='old')
              OPEN (8888,FILE='BTE/META/xqf_fbz_cart',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='old')
              OPEN (9999,FILE=rate_RTA_3D_fmt,STATUS='replace')
                 !
              DO iq = 1, nq_tot
                 !
                 IF (phdrag) THEN
                   iq_irr = rirr2irr_q(iq)
                   READ (7777,REC=iq_irr) iq_ful
                 ELSE
                   READ (7777,REC=iq) iq_ful
                 ENDIF
                 !
                 READ (8888,REC=iq_ful) xq_fbz(1:3)
                 length = SQRT(DOT_PRODUCT(xq_fbz,xq_fbz))
                 !
                 DO imode = 1, nmodes
                 !
                 IF (.NOT. phdrag) THEN
                   WRITE (9999,'(i8,f12.6,f11.4,es14.4)') iq, length, wf_irr(imode,iq)*rydcm1, 2.0d0*gammai_mode_all(itemp,idope,imode,iq)*ryd2thz*twopi 
                 ELSE
                   WRITE (9999,'(i8,f12.6,f11.4,es14.4)') iq, length, wf_all(imode,iq)*rydcm1, 2.0d0*gammai_mode_all(itemp,idope,imode,iq)*ryd2thz*twopi
                 ! [cm-1] [THz]
                 ENDIF
                 !    
                 ENDDO
                 !            
              ENDDO
              CLOSE (7777)
              CLOSE (8888)
              CLOSE (9999)
              !
              !

        ENDDO
     ENDDO
     !
  ENDIF
  !
     CALL mp_barrier (inter_pool_comm)
  !
END SUBROUTINE export_rate_ph



!----------------------------------------------------------------------------
SUBROUTINE export_ShengBTE ()
!----------------------------------------------------------------------------
! 
! export tau used in ShengBTE
! xq_ful(1:3) [cryst], scattering rate(1:nmodes) [THz,ps^-1]
!
! export el-ph scattering rate, Jan, 2018 
! in the format of BTE.qpoints, BTE.v, BTE.w_final, BTE.omega, BTE.cumulative_kappa_scalar
!
#INCLUDE "f_defs.h"
  USE kinds,         ONLY : DP
  USE phcom,         ONLY : nmodes
  USE epwcom,        ONLY : bte, eptemp, epdope, neptemp, nepdope, nqf1, nqf2, nqf3
  USE elph2,         ONLY : gammai_mode_all, wf_irr, vph_irr
  USE bte_var
  USE bte_func
  USE constants_epw, ONLY : ryd2ev, rydcm1, ryd2thz, twopi, au2m, au2s, kelvin2eV
#ifdef __PARA
  USE io_global,     ONLY : ionode_id, stdout, ionode
  USE mp,            ONLY : mp_barrier, mp_bcast, mp_sum
  USE mp_global,     ONLY : my_pool_id, nproc_pool, intra_image_comm, &
                            inter_pool_comm, me_pool, root_pool, intra_pool_comm, npool
#ENDIF
  !
  IMPLICIT NONE
  !
  INTEGER             :: iq, iq_ful, imode, itemp, idope, i, j, k
  REAL(KIND=DP)       :: rate_ful(nmodes), xq_irr(3), wf_weight
  CHARACTER(LEN=256)  :: ShengBTE_fmt, bte_w_final, bte_v, bte_qpoints
  CHARACTER(LEN=12)   :: txnx
  CHARACTER(LEN=3)    :: itemp_num, idope_num
  !
  !
  IF (my_pool_id .EQ. ionode_id) THEN
     !
      bte_qpoints = 'BTE/epw-BTE.qpoints'
      bte_v = 'BTE/epw-BTE.v'
      OPEN (2345,FILE=bte_v,STATUS='replace') 
      OPEN (3456,FILE=bte_qpoints,STATUS='replace')   
      OPEN (3333,FILE='BTE/META/xqf_irr_cryst',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='old')
      OPEN (1155,FILE='BTE/META/wqf_irr',FORM='unformatted',ACCESS='direct',RECL=DP,STATUS='old')
      OPEN (7777,FILE='BTE/META/irr2ful_q',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='old')
     !
      DO imode = 1, nmodes
         DO iq = 1, nq_irr     
            WRITE(2345,'(3e20.10)') vph_irr(:,imode,iq)*(au2m/au2s)/1000.0d0   
         ENDDO       
      ENDDO
     !
     DO iq = 1, nq_irr     
      !    
      READ (7777,REC=iq) iq_ful
      !
      READ (3333,REC=iq) xq_irr(1), xq_irr(2), xq_irr(3)
      READ (1155,REC=iq) wf_weight
      !
      WRITE (3456,'(i9,x,i9,x,i9,x,3(e20.10,1x))') iq, iq_ful, NINT(nqf1*nqf2*nqf3*wf_weight), xq_irr(1:3)
      !
      !
     ENDDO
      CLOSE (7777)
      CLOSE (2345)
      CLOSE (3456)
      CLOSE (3333)
      CLOSE (1155)
     !
     DO itemp = 1, neptemp
        DO idope = 1, nepdope
           !
           WRITE(itemp_num,'(i3)') itemp
           WRITE(idope_num,'(i3)') idope
           txnx = 'T'//TRIM(ADJUSTL(itemp_num))//'_N'//TRIM(ADJUSTL(idope_num))
           ShengBTE_fmt = 'BTE/'//TRIM(ADJUSTL(txnx))//'/ShengBTE.phe'
           bte_w_final = 'BTE/'//TRIM(ADJUSTL(txnx))//'/epw-BTE.w_final'
           !
           ! ShengBTE phonon_electron=.true.
           OPEN (2222,FILE=ShengBTE_fmt,STATUS='replace') 
           !
           WRITE (2222,'(3i4,f10.4,es14.4)') nqf1, nqf2, nqf3, eptemp(itemp)*ryd2ev/kelvin2eV, epdope(idope)
           !
           DO i = 1, nqf1
              DO j = 1, nqf2
                 DO k = 1, nqf3
                    !
                    iq = (i-1)*nqf2*nqf3 + (j-1)*nqf3 + k
                    !
                    DO imode = 1, nmodes
                       !
                       IF (gammai_mode_all(itemp,idope,imode,ful2irr_q(iq)) .GT. 0.0d0) THEN
                          rate_ful(imode) = 2.0d0*gammai_mode_all(itemp,idope,imode,ful2irr_q(iq))*ryd2thz*twopi ! THz, ps^-1
                       ELSE
                          rate_ful(imode) = 2.0d0*ABS(gammai_mode_all(itemp,idope,imode,ful2irr_q(iq)))*ryd2thz*twopi ! THz, ps^-1
                       ENDIF
                       !
                    ENDDO    
                    !
                    WRITE (2222,'(96es18.10)') DBLE(i-1)/DBLE(nqf1), DBLE(j-1)/DBLE(nqf2), DBLE(k-1)/DBLE(nqf3), rate_ful(1:nmodes)
                    !
                 ENDDO 
              ENDDO
           ENDDO  
           !
           CLOSE (2222)
           !
           OPEN (1234,FILE=bte_w_final,STATUS='replace') 
           !
           DO imode = 1, nmodes
             DO iq = 1, nq_irr
             !
              WRITE(1234,'(2e20.10)') wf_irr(imode,iq)*ryd2thz*twopi, 2.0d0*ABS(gammai_mode_all(itemp,idope,imode,iq))*ryd2thz*twopi       
             ENDDO
             !
           ENDDO
           !
           CLOSE (1234)
           !
        ENDDO
     ENDDO
     !
  ENDIF
  !
     CALL mp_barrier (inter_pool_comm)
  !
END SUBROUTINE export_ShengBTE



!----------------------------------------------------------------------------
SUBROUTINE export_channel (itemp, idope)
!----------------------------------------------------------------------------
! 
! export 
!
! xkf_fbz(1:3), length (2pi/a), phonon frequency (THz), scattering rate (THz)
!
#INCLUDE "f_defs.h"
  USE kinds,         ONLY : DP
  USE phcom,         ONLY : nmodes
  USE constants_epw, ONLY : ryd2thz, rydcm1, twopi
  USE bte_var
  USE bte_func
#ifdef __PARA
  USE mp,            ONLY : mp_bcast, mp_sum
  USE io_global,     ONLY : ionode_id
  USE mp,            ONLY : mp_barrier
  USE mp_global,     ONLY : my_pool_id, inter_pool_comm
#ENDIF
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)     :: itemp, idope
  REAL(KIND=DP)           :: wf, xk_fbz(3), leng(nq_ful), leng_d, leng_u, leng_del, non_r(3), scatt_ch(nmodes), scatt_ch_cumu(nmodes)
  INTEGER                 :: leng_num=500
  INTEGER                 :: ik, iq, imode, il
  CHARACTER(LEN=256)      :: channel_fmt, channel_dist_fmt
  CHARACTER(LEN=12)       :: txnx
  CHARACTER(LEN=3)        :: itemp_num, idope_num, ik_num, ibnd_num
  !
  !
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     WRITE(itemp_num,'(i3)') itemp
     WRITE(idope_num,'(i3)') idope
     txnx = 'T'//TRIM(ADJUSTL(itemp_num))//'_N'//TRIM(ADJUSTL(idope_num))
     !
     OPEN (2222,FILE='BTE/META/xqf_fbz_cart',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='old')
     OPEN (3333,FILE='BTE/META/phonon_ful',FORM='unformatted',ACCESS='direct',RECL=4*DP,STATUS='old')
     !
     DO ik = 1, nk_kmesh
        !
        WRITE(ik_num,'(i3)') ik
        WRITE(ibnd_num,'(i3)') band_ch(ik)
        channel_fmt = 'BTE/'//TRIM(ADJUSTL(txnx))//'/channel_k-'//TRIM(ADJUSTL(ik_num))//'_band-'//TRIM(ADJUSTL(ibnd_num))//'.dat'
        !
        leng = 0.0d0
        OPEN (9999,FILE=channel_fmt,STATUS='replace')
        DO iq = 1, nq_ful
           !
           READ (2222,REC=iq) xk_fbz(1:3)
           leng(iq) = SQRT(DOT_PRODUCT(xk_fbz,xk_fbz))
           !
           DO imode = 1, nmodes
              !
              READ (3333,REC=(iq-1)*nmodes+imode) wf, non_r(1:3)
              !
              WRITE (9999,'(3es14.4)') leng(iq), wf*rydcm1, 2.0d0*sigmai_ch(itemp,idope,imode,iq,ik)*ryd2thz*twopi
              !
           ENDDO
           !
        ENDDO
        CLOSE (9999)
        !
        !
        leng_del = MAXVAL(leng(:))/DBLE(leng_num)
        leng_d = 0.0d0
        leng_u = leng_del
        !
        channel_dist_fmt = 'BTE/'//TRIM(ADJUSTL(txnx))//'/channel_dist_k-'//TRIM(ADJUSTL(ik_num))//'_band-'//TRIM(ADJUSTL(ibnd_num))//'.dat'
        OPEN (8888,FILE=channel_dist_fmt,STATUS='replace')
        DO il = 1, leng_num
           !
           scatt_ch = 0.0d0
           scatt_ch_cumu = 0.0d0
           !
           DO iq = 1, nq_ful         
              DO imode = 1, nmodes
                 !
                 IF (leng(iq) .GT. leng_d .AND. leng(iq) .LE. leng_u) scatt_ch(imode) = scatt_ch(imode) + 2.0d0*sigmai_ch(itemp,idope,imode,iq,ik)*ryd2thz*twopi
                 IF (leng(iq) .LE. leng_u) scatt_ch_cumu(imode) = scatt_ch_cumu(imode) + 2.0d0*sigmai_ch(itemp,idope,imode,iq,ik)*ryd2thz*twopi

                 !
              ENDDO
           ENDDO
           !
           WRITE (8888,'(100es14.4)') leng_d, scatt_ch(1:nmodes), scatt_ch_cumu(1:nmodes)
           !
           leng_d = leng_u
           leng_u = leng_u + leng_del
           !
        ENDDO
        CLOSE (8888)
        !
     ENDDO
     !
     CLOSE (2222)
     CLOSE (3333)
     !
  ENDIF
  !
  CALL mp_barrier (inter_pool_comm)
  !
END SUBROUTINE export_channel



!-------------------------------------------------------------------------------
SUBROUTINE export_ft (iter, converge, itemp, idope)
!-------------------------------------------------------------------------------
! 
! export electron distribution function
! 
! ik (irr), ibnd (red), xkf_fbz [cart], energy [eV], f(t) [none]
!
#INCLUDE "f_defs.h"
  USE kinds,         ONLY : DP
  USE pwcom,         ONLY : ef 
  USE phcom,         ONLY : nmodes
  USE epwcom,        ONLY : nbndsub, dt
  USE elph2,         ONLY : ibndmin, ibndmax, nbnd_red, cbnd_emin, vbnd_emax, &
                            etf_all, f_t_ful, f_0_ful
  USE bte_var
  USE constants_epw, ONLY : ryd2ev, au2fs
#ifdef __PARA
  USE io_global,     ONLY : ionode_id
  USE mp_global,     ONLY : my_pool_id
#ENDIF
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iter, itemp, idope
  LOGICAL, INTENT(IN) :: converge
  INTEGER             :: ik, ik_ful, ik_irr, ibnd, ibnd0, imode
  CHARACTER(LEN=256)  :: filename_0, filename_t, file_num
  CHARACTER(LEN=12)   :: txnx
  CHARACTER(LEN=3)    :: itemp_num, idope_num
  REAL(KIND=DP)       :: xk_fbz(3)
  !
  !
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     OPEN (8888,FILE='BTE/META/xkf_fbz_cart',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='old')
     !
     WRITE(itemp_num,'(i3)') itemp
     WRITE(idope_num,'(i3)') idope
     txnx = 'T'//TRIM(ADJUSTL(itemp_num))//'_N'//TRIM(ADJUSTL(idope_num))
  !   bte_fmt = 'BTE/'//TRIM(ADJUSTL(txnx))//'/bte_'//TRIM(ADJUSTL(txnx))//'.out'
     !
     IF (iter .EQ. 0) THEN
        filename_0 = 'BTE/'//TRIM(ADJUSTL(txnx))//'/ft_0.dat'
        OPEN (9998,FILE=filename_0,STATUS='replace')
     ENDIF
     !
     WRITE (file_num,*) iter+1
     filename_t = 'BTE/'//TRIM(ADJUSTL(txnx))//'/ft_'//TRIM(ADJUSTL(file_num))//'.dat'
     OPEN (9999,FILE=filename_t,STATUS='replace')
     !
     DO ik = 1, nk_ful_red
        !
        ik_ful = rful2ful(ik)
        ik_irr = rful2rirr(ik)
        !
        DO ibnd = 1, nbnd_red
           !
           ibnd0 = ibnd+ibndmin-1
           !
            IF (f_0_ful(ibnd,ik) .NE. 0.0d0) THEN
              !
              READ (8888,REC=ik_ful) xk_fbz(1:3)
              !
              IF (iter .EQ. 0) &
              WRITE (9998,'(i12,i5,3f12.6,f10.4,es14.4)') ik_ful, ibnd0, xk_fbz(1:3), etf_all(ibnd0,ik_irr)*ryd2ev, f_0_ful(ibnd,ik)
              !
              WRITE (9999,'(i12,i5,3f12.6,f10.4,es14.4)') ik_ful, ibnd0, xk_fbz(1:3), etf_all(ibnd0,ik_irr)*ryd2ev, f_t_ful(ibnd,ik)
              !
           ENDIF
           !
        ENDDO 
     ENDDO 
     !
     CLOSE (8888)
     IF (iter .EQ. 0) CLOSE (9998)
     CLOSE (9999)
     !
     !
  ENDIF
  !
  !
END SUBROUTINE export_ft



!-------------------------------------------------------------------------------
SUBROUTINE export_result (transpt, iter, converge, itemp, idope, time)
!-------------------------------------------------------------------------------
!
#INCLUDE "f_defs.h"
  USE kinds,         ONLY : DP
  USE phcom,         ONLY : nmodes
  USE pwcom,         ONLY : ef 
  USE spin_orb,      ONLY : lspinorb
  USE epwcom,        ONLY : nkf1, nkf2, nkf3, nqf1, nqf2, nqf3, epdim, phdrag, neptemp, nepdope, &
                            nbndsub, eptemp, efermi_read, fsthick, epthick, fermi_energy, epdope, egap_rbm, mixing, &
                            lpolar, eig_read, asr_eph, nptype, bte, dt, run, smearing, eimp_mode, &
                            elop, screen_polar, eimp_ls_mode, eph_interp
  USE elph2,         ONLY : ibndmin, ibndmax, n_hole, n_elec, n_intr, cbnd_emin, vbnd_emax, vfsthick, cfsthick, ef_epw, dope_ef
  USE bte_var
  USE constants_epw, ONLY : kB, ryd2ev, au2fs
#ifdef __PARA
  USE io_global,     ONLY : ionode_id
  USE mp_global,     ONLY : my_pool_id
#ENDIF
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)       :: iter, itemp, idope
  REAL(KIND=DP), INTENT(IN) :: transpt(3,3,24), time
  LOGICAL, INTENT(IN)       :: converge
  CHARACTER(LEN=256)        :: bte_fmt
  CHARACTER(LEN=12)         :: txnx
  CHARACTER(LEN=3)          :: itemp_num, idope_num
  REAL(KIND=DP)             :: ef0, doping, temp, mob0, cod0, seb0, etc0, pf0, &
                               pft0
  ! date and time
  CHARACTER(LEN=8)          :: date_
  CHARACTER(LEN=10)         :: time_
  CHARACTER(LEN=5)          :: zone_
  INTEGER                   :: values_(8)
  !
  !
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     CALL DATE_AND_TIME (date_,time_,zone_,values_)
     !
     WRITE(itemp_num,'(i3)') itemp
     WRITE(idope_num,'(i3)') idope
     txnx = 'T'//TRIM(ADJUSTL(itemp_num))//'_N'//TRIM(ADJUSTL(idope_num))
     !bte_fmt = 'BTE/'//TRIM(ADJUSTL(txnx))//'/bte_'//TRIM(ADJUSTL(txnx))//'.out'
     bte_fmt = 'BTE/'//TRIM(ADJUSTL(txnx))//'/bte.out'
     !
     IF (iter .EQ. 0) THEN
        !
        IF (efermi_read) THEN
           ef0 = fermi_energy
           doping = dope_ef(itemp)
        ELSE
           ef0    = ef_epw(itemp,idope)
           temp   = eptemp(itemp)/kB
           doping = epdope(idope)
        ENDIF
        !
        OPEN (100000,FILE=bte_fmt,STATUS='replace')
        !      
        WRITE (100000,'(7x,a)')                     ''                                                        
        WRITE (100000,'(7x,a)')                     '~B(                 .BNz.          zNNN+                        '
        WRITE (100000,'(7x,a)')                     '~NNNNNNNNNNNNNNNNNNNNNNN+          zNh+                         '
        WRITE (100000,'(7x,a)')                     '~Nz  .         ~     Nz    ~NNs~   zN~                          '
        WRITE (100000,'(7x,a)')                     '~Nz  ~N-      .NNND  Nz     NN     zN~                          '
        WRITE (100000,'(7x,a)')                     '~Nz    Nh     zN+    Nz     NN     zN~                          '
        WRITE (100000,'(7x,a)')                     '~Nz    (Nz   =N~     Nz     NN     zN~                          '
        WRITE (100000,'(7x,a)')                     '~Nz.+<~~Dz~~zh~+NNB( Nz     NN     zN~                          '
        WRITE (100000,'(7x,a)')                     '~Nz  -    ~Nz        Nz     NN     zN~                          '
        WRITE (100000,'(7x,a)')                     '~Nz       ~Nz        Nz     NN     zN~                          '
        WRITE (100000,'(7x,a)')                     '~Nz  +NNB-~Nz  ~NNN< Nz     NN     zN~                          '
        WRITE (100000,'(7x,a)')                     '~Nz  +N+  ~Nz  ~Nz   Nz     NN     zN~                          '
        WRITE (100000,'(7x,a)')                     '~Nz  +N<  ~Nz  ~Nz   Nz     NN     zN~                          '
        WRITE (100000,'(7x,a)')                     '~Nz  +N+  ~Nz  ~Nz   Nz     NN     zN~                          '
        WRITE (100000,'(7x,a)')                     '~Nz  <N+  ~Nz  ~Nz   Nz     NN     zN~                          '
        WRITE (100000,'(7x,a)')                     '~Nz  +N+  ~Nz  ~Nz   Nz    ~NN     zN~                          '
        WRITE (100000,'(7x,a)')                     '~Nz  zNNNNNNNNNNNz   Nz    .s+     zN~                          '
        WRITE (100000,'(7x,a)')                     '~Nz             +(   Nz            zN~                          '
        WRITE (100000,'(7x,a)')                     '~Nz                  Nz            zN~                          '
        WRITE (100000,'(7x,a)')                     '+Nz           ~+zNNNNNz     <Bzz<+zNN~                          '
        WRITE (100000,'(7x,a)')                     '+Nz                hN+           sNB<                           '
        WRITE (100000,'(7x,a)')                     '' 
        WRITE (100000,'(7x,a)')                     '' 
        WRITE (100000,'(7x,a)')                     'MIT NANOENGINEERING GROUP                                       '
        WRITE (100000,'(7x,a)')                     'TE-HUAN LIU and JIAWEI ZHOU                                     '
       !WRITE (100000,'(7x,a)')                     'ELECTRON-PHONON BOLTZMANN TRANSPORT EQUATION SOLVER (v2.0)      ' ! 05/25/2016 
       !WRITE (100000,'(7x,a)')                     'ELECTRON-PHONON BOLTZMANN TRANSPORT EQUATION SOLVER (v2.1)      ' ! 09/13/2016
        WRITE (100000,'(7x,a)')                     'ELECTRON-PHONON BOLTZMANN TRANSPORT EQUATION SOLVER (v2.2)      ' ! 01/05/2017
        WRITE (100000,'(7x,a)')                     ''    
        WRITE (100000,'(7x,a,i2,a,i2,a,i4,a,i2,a,i2,a,i2)') 'Program starts on ', values_(2), '/', values_(3), '/', values_(1), ' at ', &
                                                                                  values_(5), ':', values_(6), ':', values_(7)
        WRITE (100000,'(7x,a)')                     '' 
        WRITE (100000,'(7x,a)')                     '======================================================='                                                        
        WRITE (100000,'(7x,a)')                     ''        
        WRITE (100000,'(7x,a,i4,a,i4,a,i4,a,i12,a)') 'k-mesh: ', nkf1, ' * ', nkf2, ' * ', nkf3, '  = ', nkf1*nkf2*nkf3, ' points'
        WRITE (100000,'(7x,a,i4,a,i4,a,i4,a,i12,a)') 'q-mesh: ', nqf1, ' * ', nqf2, ' * ', nqf3, '  = ', nqf1*nqf2*nqf3, ' points'
        WRITE (100000,'(a)') ''
        WRITE (100000,'(7x,a,i5)')                  'Dimension   = ', epdim
        WRITE (100000,'(7x,a,i5,a,i3)')             'Band (e)    = ', ibndmin, '  to', ibndmax
        WRITE (100000,'(7x,a,i5)')                  'Mode (ph)   = ', nmodes
        WRITE (100000,'(7x,a,f10.4,a,f10.4,a)')     'Bandgap     = ', vbnd_emax*ryd2ev, '  to', cbnd_emin*ryd2ev, ' [eV]'
        WRITE (100000,'(7x,a,f10.4,a,f10.4,a)')     'fsthick     = ', vfsthick*ryd2ev, '  / ', cfsthick*ryd2ev, ' [eV]'
        WRITE (100000,'(7x,a,f10.4,a)')             'epthick     = ', epthick*ryd2ev, ' [eV]'
        WRITE (100000,'(a)') ''
        WRITE (100000,'(7x,a,f7.1,a)')              'Temperature = ', temp, ' [K]'
        WRITE (100000,'(7x,a,es14.4,a)')            'Doping      = ', doping, ' [cm^-3]'
        WRITE (100000,'(7x,a,f10.4,a)')             'Fermi level = ', ef0 * ryd2ev, ' [eV]'
        WRITE (100000,'(7x,a,es14.4,a)')            'n_hole      = ', n_hole(itemp,idope), ' [cm^-3]'
        WRITE (100000,'(7x,a,es14.4,a)')            'n_elec      = ', n_elec(itemp,idope), ' [cm^-3]'
        WRITE (100000,'(7x,a,es14.4,a)')            'n_intr      = ', n_intr(itemp,idope), ' [cm^-3]'
        WRITE (100000,'(a)') ''
        WRITE (100000,'(7x,a,a5)')                  'Smearing                 =    ', smearing
        WRITE (100000,'(7x,a,f8.4)')                'Mixing ratio             = ', mixing
        WRITE (100000,'(7x,a,l8)')                  'GW band structure        = ', eig_read
        WRITE (100000,'(7x,a,f8.4,a)')              'RBM bandgap              = ', egap_rbm, ' [eV]'
        WRITE (100000,'(7x,a,l8)')                  'Polar interaction        = ', lpolar
        WRITE (100000,'(7x,a,l8)')                  'Spin-orbital interaction = ', lspinorb
        WRITE (100000,'(7x,a,l8)')                  'ASR to ep-matrix (eph)   = ', asr_eph
        WRITE (100000,'(7x,a,a5/)')                 'n/p-type semiconductor   =        ', nptype
        WRITE (100000,'(7x,a,2x,i6,a)')             'e-imp scattering mode    = ', eimp_mode
        WRITE (100000,'(7x,a,2x,i6,a)')             '    (0: none; 1: Brooks-Herring relaxation time; 2: B-H momentum relaxation time)'
        WRITE (100000,'(7x,a,2x,i6,a)')             'e-imp-ls mode    = ', eimp_ls_mode
        WRITE (100000,'(7x,a,2x,i6,a)')             '    (0: both long and short range; 1: only long range; '
        WRITE (100000,'(7x,a,2x,i6,a)')             '     2: only short range;          3: none )'
        !
        IF (bte .EQ. 19 .OR. bte .EQ. 18) THEN
           !
           WRITE (100000,'(7x,a)')                  '-------------------------------------------------------' 
           WRITE (100000,'(7x,a)')                  '               Random k-point calculation              '
           WRITE (100000,'(7x,a)')                  '' 
           WRITE (100000,'(7x,a,i5)')               'Number of random k point: ', nk_kmesh
           WRITE (100000,'(7x,a)')                  '' 
           !
        ENDIF
        !
        IF (bte .EQ. 3 .OR. bte .EQ. 30) THEN
           !
           WRITE (100000,'(7x,a)')                  '-------------------------------------------------------------------------' 
           WRITE (100000,'(7x,a)')                  '                          Time-dependent ep-BTE                          '
           WRITE (100000,'(7x,a)')                  '' 
           WRITE (100000,'(7x,a,f6.2,a)')           'Timestep = ', dt*au2fs, ' [fs]'
           WRITE (100000,'(7x,a)')                  '' 
           !
        ENDIF
        !
     ELSE
        !
        OPEN (100000,FILE=bte_fmt,POSITION='append')
        !
     ENDIF
     !
     IF (bte .EQ. 0 .OR. bte .EQ. 1 .OR. bte .EQ. 10 .OR. bte .EQ. -1 .OR. bte .EQ. 19 .OR. bte .EQ. 18) THEN
        !
        WRITE (100000,'(7x,a)')                     '-------------------------------------------------------' 
        WRITE (100000,'(7x,a,i5,a,f8.1,a)')         'ITERATION:', iter, '                     CPU TIME:', time, ' s'
        WRITE (100000,'(7x,a)')                     '-------------------------------------------------------'
        WRITE (100000,'(7x,a,3es16.6)')             '  MOB |', transpt(1,1,1), transpt(1,2,1), transpt(1,3,1)  
        WRITE (100000,'(7x,a,3es16.6)')             '      |', transpt(2,1,1), transpt(2,2,1), transpt(2,3,1)  
        WRITE (100000,'(7x,a,3es16.6)')             '      |', transpt(3,1,1), transpt(3,2,1), transpt(3,3,1) 
        WRITE (100000,'(7x,a)')                     '-------------------------------------------------------'
        WRITE (100000,'(7x,a,3es16.6)')             'MOB_e |', transpt(1,1,16), transpt(1,2,16), transpt(1,3,16)  
        WRITE (100000,'(7x,a,3es16.6)')             '      |', transpt(2,1,16), transpt(2,2,16), transpt(2,3,16)  
        WRITE (100000,'(7x,a,3es16.6)')             '      |', transpt(3,1,16), transpt(3,2,16), transpt(3,3,16)  
        WRITE (100000,'(7x,a)')                     '-------------------------------------------------------'
        WRITE (100000,'(7x,a,3es16.6)')             'MOB_h |', transpt(1,1,17), transpt(1,2,17), transpt(1,3,17)  
        WRITE (100000,'(7x,a,3es16.6)')             '      |', transpt(2,1,17), transpt(2,2,17), transpt(2,3,17)  
        WRITE (100000,'(7x,a,3es16.6)')             '      |', transpt(3,1,17), transpt(3,2,17), transpt(3,3,17)   
        WRITE (100000,'(7x,a)')                     '      -------------------------------------------------'
        WRITE (100000,'(7x,a,3es16.6)')             '  COD |', transpt(1,1,2), transpt(1,2,2), transpt(1,3,2)  
        WRITE (100000,'(7x,a,3es16.6)')             '      |', transpt(2,1,2), transpt(2,2,2), transpt(2,3,2)  
        WRITE (100000,'(7x,a,3es16.6)')             '      |', transpt(3,1,2), transpt(3,2,2), transpt(3,3,2) 
        WRITE (100000,'(7x,a)')                     '      -------------------------------------------------'
        WRITE (100000,'(7x,a,3es16.6)')             'COD_e |', transpt(1,1,18), transpt(1,2,18), transpt(1,3,18)  
        WRITE (100000,'(7x,a,3es16.6)')             '      |', transpt(2,1,18), transpt(2,2,18), transpt(2,3,18)  
        WRITE (100000,'(7x,a,3es16.6)')             '      |', transpt(3,1,18), transpt(3,2,18), transpt(3,3,18) 
        WRITE (100000,'(7x,a)')                     '      -------------------------------------------------'
        WRITE (100000,'(7x,a,3es16.6)')             'COD_h |', transpt(1,1,19), transpt(1,2,19), transpt(1,3,19)  
        WRITE (100000,'(7x,a,3es16.6)')             '      |', transpt(2,1,19), transpt(2,2,19), transpt(2,3,19)  
        WRITE (100000,'(7x,a,3es16.6)')             '      |', transpt(3,1,19), transpt(3,2,19), transpt(3,3,19) 
        WRITE (100000,'(7x,a)')                     '      -------------------------------------------------'
        WRITE (100000,'(7x,a,3es16.6)')             '  ETC |', transpt(1,1,3), transpt(1,2,3), transpt(1,3,3)  
        WRITE (100000,'(7x,a,3es16.6)')             '      |', transpt(2,1,3), transpt(2,2,3), transpt(2,3,3)  
        WRITE (100000,'(7x,a,3es16.6)')             '      |', transpt(3,1,3), transpt(3,2,3), transpt(3,3,3)
        WRITE (100000,'(7x,a)')                     '      -------------------------------------------------'
        WRITE (100000,'(7x,a,3es16.6)')             '  SEB |', transpt(1,1,4), transpt(1,2,4), transpt(1,3,4)  
        WRITE (100000,'(7x,a,3es16.6)')             '      |', transpt(2,1,4), transpt(2,2,4), transpt(2,3,4)  
        WRITE (100000,'(7x,a,3es16.6)')             '      |', transpt(3,1,4), transpt(3,2,4), transpt(3,3,4) 
        WRITE (100000,'(7x,a)')                     '      -------------------------------------------------'
        WRITE (100000,'(7x,a,3es16.6)')             'SEB_e |', transpt(1,1,20), transpt(1,2,20), transpt(1,3,20)  
        WRITE (100000,'(7x,a,3es16.6)')             '      |', transpt(2,1,20), transpt(2,2,20), transpt(2,3,20)  
        WRITE (100000,'(7x,a,3es16.6)')             '      |', transpt(3,1,20), transpt(3,2,20), transpt(3,3,20) 
        WRITE (100000,'(7x,a)')                     '      -------------------------------------------------'
        WRITE (100000,'(7x,a,3es16.6)')             'SEB_h |', transpt(1,1,21), transpt(1,2,21), transpt(1,3,21)  
        WRITE (100000,'(7x,a,3es16.6)')             '      |', transpt(2,1,21), transpt(2,2,21), transpt(2,3,21)  
        WRITE (100000,'(7x,a,3es16.6)')             '      |', transpt(3,1,21), transpt(3,2,21), transpt(3,3,21) 
        WRITE (100000,'(7x,a)')                     '      -------------------------------------------------'
        WRITE (100000,'(7x,a,3es16.6)')             '  PEL |', transpt(1,1,5), transpt(1,2,5), transpt(1,3,5)  
        WRITE (100000,'(7x,a,3es16.6)')             '      |', transpt(2,1,5), transpt(2,2,5), transpt(2,3,5)  
        WRITE (100000,'(7x,a,3es16.6)')             '      |', transpt(3,1,5), transpt(3,2,5), transpt(3,3,5) 
        WRITE (100000,'(7x,a)')                     '      -------------------------------------------------'
        WRITE (100000,'(7x,a,3es16.6)')             '  POF |', transpt(1,1,6), transpt(1,2,6), transpt(1,3,6)  
        WRITE (100000,'(7x,a,3es16.6)')             '      |', transpt(2,1,6), transpt(2,2,6), transpt(2,3,6)  
        WRITE (100000,'(7x,a,3es16.6)')             '      |', transpt(3,1,6), transpt(3,2,6), transpt(3,3,6) 
        WRITE (100000,'(7x,a)')                     '      -------------------------------------------------'
        WRITE (100000,'(7x,a,3es16.6)')             '  LOZ |', transpt(1,1,7), transpt(1,2,7), transpt(1,3,7)  
        WRITE (100000,'(7x,a,3es16.6)')             '      |', transpt(2,1,7), transpt(2,2,7), transpt(2,3,7)  
        WRITE (100000,'(7x,a,3es16.6)')             '      |', transpt(3,1,7), transpt(3,2,7), transpt(3,3,7) 
        WRITE (100000,'(7x,a)')                     '      -------------------------------------------------'
        IF (phdrag) WRITE (100000,'(7x,a,3es16.6)') '  LTC |', transpt(1,1,8), transpt(1,2,8), transpt(1,3,8)  
        IF (phdrag) WRITE (100000,'(7x,a,3es16.6)') '      |', transpt(2,1,8), transpt(2,2,8), transpt(2,3,8)  
        IF (phdrag) WRITE (100000,'(7x,a,3es16.6)') '      |', transpt(3,1,8), transpt(3,2,8), transpt(3,3,8) 
        IF (phdrag) WRITE (100000,'(7x,a)')         '      -------------------------------------------------'
        IF (phdrag) WRITE (100000,'(7x,a,3es16.6)') ' SLTC |', transpt(1,1,9), transpt(1,2,9), transpt(1,3,9)  
        IF (phdrag) WRITE (100000,'(7x,a,3es16.6)') '      |', transpt(2,1,9), transpt(2,2,9), transpt(2,3,9)  
        IF (phdrag) WRITE (100000,'(7x,a,3es16.6)') '      |', transpt(3,1,9), transpt(3,2,9), transpt(3,3,9)
        WRITE (100000,'(7x,a)')                     '      -------------------------------------------------'
        WRITE (100000,'(7x,a,3es16.6)')             '  BTC |', transpt(1,1,10), transpt(1,2,10), transpt(1,3,10)
        WRITE (100000,'(7x,a,3es16.6)')             '      |', transpt(2,1,10), transpt(2,2,10), transpt(2,3,10)
        WRITE (100000,'(7x,a,3es16.6)')             '      |', transpt(3,1,10), transpt(3,2,10), transpt(3,3,10) 
        IF (phdrag) WRITE (100000,'(7x,a)')         '      -------------------------------------------------'
        IF (phdrag) WRITE (100000,'(7x,a,3es16.6)') '  SPD |', transpt(1,1,11), transpt(1,2,11), transpt(1,3,11)  
        IF (phdrag) WRITE (100000,'(7x,a,3es16.6)') '      |', transpt(2,1,11), transpt(2,2,11), transpt(2,3,11)  
        IF (phdrag) WRITE (100000,'(7x,a,3es16.6)') '      |', transpt(3,1,11), transpt(3,2,11), transpt(3,3,11) 
        IF (phdrag) WRITE (100000,'(7x,a/)')        '      -------------------------------------------------' 
        IF (phdrag) WRITE (100000,'(7x,a,3es16.6)') '  SPK |', transpt(1,1,12), transpt(1,2,12), transpt(1,3,12)  
        IF (phdrag) WRITE (100000,'(7x,a,3es16.6)') '      |', transpt(2,1,12), transpt(2,2,12), transpt(2,3,12)  
        IF (phdrag) WRITE (100000,'(7x,a,3es16.6)') '      |', transpt(3,1,12), transpt(3,2,12), transpt(3,3,12) 
        IF (phdrag) WRITE (100000,'(7x,a/)')        '      -------------------------------------------------'
        IF (phdrag) WRITE (100000,'(7x,a,3es16.6)') '  SPX |', transpt(1,1,13), transpt(1,2,13), transpt(1,3,13)  
        IF (phdrag) WRITE (100000,'(7x,a,3es16.6)') '      |', transpt(2,1,13), transpt(2,2,13), transpt(2,3,13)  
        IF (phdrag) WRITE (100000,'(7x,a,3es16.6)') '      |', transpt(3,1,13), transpt(3,2,13), transpt(3,3,13) 
        IF (phdrag) WRITE (100000,'(7x,a/)')        '      -------------------------------------------------' 
        IF (phdrag) WRITE (100000,'(7x,a,3es16.6)') '  SPL |', transpt(1,1,14), transpt(1,2,14), transpt(1,3,14)  
        IF (phdrag) WRITE (100000,'(7x,a,3es16.6)') '      |', transpt(2,1,14), transpt(2,2,14), transpt(2,3,14)  
        IF (phdrag) WRITE (100000,'(7x,a,3es16.6)') '      |', transpt(3,1,14), transpt(3,2,14), transpt(3,3,14) 
        IF (converge) WRITE (100000,'(7x,a/)')      '-------------------------------------------------------' 
        IF (phdrag) WRITE (100000,'(7x,a,3es16.6)') '  SPG |', transpt(1,1,15), transpt(1,2,15), transpt(1,3,15)  
        IF (phdrag) WRITE (100000,'(7x,a,3es16.6)') '      |', transpt(2,1,15), transpt(2,2,15), transpt(2,3,15)  
        IF (phdrag) WRITE (100000,'(7x,a,3es16.6)') '      |', transpt(3,1,15), transpt(3,2,15), transpt(3,3,15) 
        IF (converge) WRITE (100000,'(7x,a/)')      '-------------------------------------------------------' 
        IF (converge) WRITE (100000,'(7x,a)')               'MOB = Electron mobility [cm^2/V/s]'
        IF (converge) WRITE (100000,'(7x,a)')               'COD = Electron conductivity [1/Ohm/m]'
        IF (converge) WRITE (100000,'(7x,a)')               'ETC = Electron thermal conductivity [W/m/K]'
        IF (converge) WRITE (100000,'(7x,a)')               'SEB = Seebeck coefficient [uV/K]'
        IF (converge) WRITE (100000,'(7x,a)')               'PEL = Peltier coefficient [uV]'
        IF (converge) WRITE (100000,'(7x,a)')               'POF = Power factor [uW/cm/K^2]'
        IF (converge) WRITE (100000,'(7x,a)')               'LOZ = Lorentz number [W*Ohm/K^2]'
        IF (converge .AND. phdrag) WRITE (100000,'(7x,a)')  'LTC = Lattice thermal conductivity [W/m/K]'
        IF (converge .AND. phdrag) WRITE (100000,'(7x,a)')  'SLTC = Lattice thermal conductivity using ph lifetime from ShengBTE [W/m/K]'
        IF (converge) WRITE (100000,'(7x,a)')               'BTC = Bipolar thermal conductivity [W/m/K]'
        IF (converge .AND. phdrag) WRITE (100000,'(7x,a/)') 'SPD = Phonon-drag Seebeck coefficient [uV/K]'
        IF (converge .AND. phdrag) WRITE (100000,'(7x,a/)') 'SPK ~ SPD but sum_k(sum_q) [uV/K]'
        IF (converge .AND. phdrag) WRITE (100000,'(7x,a/)') 'SPX ~ SPD but sum_k(sum_q) from near X valley [uV/K]'
        IF (converge .AND. phdrag) WRITE (100000,'(7x,a/)') 'SPL ~ SPD but sum_k(sum_q) from L valley [uV/K]'
        IF (converge .AND. phdrag) WRITE (100000,'(7x,a/)') 'SPG ~ SPD but sum_k(sum_q) from Gamma walley [uV/K]'
        !
     ELSEIF (bte .EQ. 3 .OR. bte .EQ. 30) THEN
        !
        WRITE (100000,'(7x,a,i10,a,f12.2,a,f8.1,a)') 'ITERATION:', iter+1, '  |  BTE TIME:', DBLE(iter+1)*dt*au2fs, ' fs  |  CPU TIME:', time, ' s'
        IF (converge) WRITE (100000,'(7x,a/)')       '-------------------------------------------------------------------------' 
        !
     ENDIF
     !
     ! ============================== write a more compact version
     ! ====================== JW
     !
     if (iter .eq. 0) then
        OPEN (100000,FILE='bte.out',POSITION='append')
        if ((itemp == 1) .and. (idope == 1)) then
           WRITE (100000,'(7x,a)')                     ''
           WRITE (100000,'(7x,a)')                     'MIT NANOENGINEERING GROUP'
           WRITE (100000,'(7x,a)')                     'T.H. LIU and J. ZHOU'
           WRITE (100000,'(7x,a)')                     'ELECTRON-PHONON BOLTZMANN TRANSPORT EQUATION SOLVER (v5.0)     '
           WRITE (100000,'(7x,a)')                     ''                      
           WRITE (100000,'(7x,a)')                     ''

           WRITE (100000,'(7x,a,i4,a,i4,a,i4,a,i12,a)') 'k-mesh: ', nkf1, ' * ', nkf2, ' * ', nkf3, '  = ', nkf1*nkf2*nkf3, ' points'
           WRITE (100000,'(7x,a,i4,a,i4,a,i4,a,i12,a)') 'q-mesh: ', nqf1, ' * ', nqf2, ' * ', nqf3, '  = ', nqf1*nqf2*nqf3, ' points'
           WRITE (100000,'(a)') ''
           WRITE (100000,'(7x,a,i5)')                  'Dimension   = ', epdim
           WRITE (100000,'(7x,a,i5,a,i3)')             'Band (e)    = ', ibndmin, '  to', ibndmax
           WRITE (100000,'(7x,a,i5)')                  'Mode (ph)   = ', nmodes
           WRITE (100000,'(7x,a,f10.4,a,f10.4,a)')     'Bandgap     = ', vbnd_emax*ryd2ev, '  to', cbnd_emin*ryd2ev, ' [eV]'
           WRITE (100000,'(7x,a,f10.4,a,f10.4,a)')     'fsthick     = ', vfsthick*ryd2ev, '  / ', cfsthick*ryd2ev, ' [eV]'
           WRITE (100000,'(a)') ''
           WRITE (100000,'(7x,a,es14.4,a,f7.1,a)')     'n_intr      = ', n_intr(1,1), ' [cm^-3] at T =', temp, ' [K]'
           WRITE (100000,'(a)') ''
           WRITE (100000,'(7x,a,a5)')                  'Smearing                 =    ', smearing
           WRITE (100000,'(7x,a,f8.4)')                'Mixing ratio             = ', mixing
           WRITE (100000,'(7x,a,l8)')                  'GW band structure        = ', eig_read
           WRITE (100000,'(7x,a,f8.4,a)')              'RBM bandgap              = ', egap_rbm, ' [eV]'
           WRITE (100000,'(7x,a,l8)')                  'Spin-orbital interaction = ', lspinorb
           WRITE (100000,'(a)') ''
           WRITE (100000,'(7x,a,l8)')                  'Polar interaction        = ', lpolar
           WRITE (100000,'(7x,a,l8)')                  'Polar e-LO phonon considered?  = ', elop
           WRITE (100000,'(7x,a,l8)')                  'Screening of polar scattering considered?  = ', screen_polar
           WRITE (100000,'(7x,a)') ''
           WRITE (100000,'(7x,a,l8)')                  'Use direct interpolation for ep-matrix?   : ', eph_interp
           WRITE (100000,'(7x,a,l8)')                  'ASR to ep-matrix (eph)   = ', asr_eph
           WRITE (100000,'(a)') ''
           WRITE (100000,'(7x,a,a5/)')                 'n/p-type semiconductor   =        ', nptype
           WRITE (100000,'(7x,a,2x,i6,a)')             'e-imp scattering mode    = ', eimp_mode
           WRITE (100000,'(7x,a,2x,i6,a)')             '    (0: none; 1: Brooks-Herring relaxation time; 2: B-H momentum relaxation time)'
           WRITE (100000,'(7x,a,2x,i6,a)')             'e-imp-ls mode    = ', eimp_ls_mode
           WRITE (100000,'(7x,a,2x,i6,a)')             '    (0: both long and short range; 1: only long range; '
           WRITE (100000,'(7x,a,2x,i6,a)')             '     2: only short range;          3: none )'

           WRITE (100000,'(7x,a)') ''
           write (100000,'()')
           WRITE (100000,'(7x,a)') '========================================================================================================================='
           WRITE (100000,'(7x,a)')                     '   T(K)  Fermi level Doping(cm^-3)  n_carr     n_charge    Mobility      cond      Seebeck        PF         PF*T    thermal cond'
        endif

        mob0 = (transpt(1,1,1)+transpt(2,2,1)+transpt(3,3,1))/3.d0
        cod0 = (transpt(1,1,2)+transpt(2,2,2)+transpt(3,3,2))/3.d0
        etc0 = (transpt(1,1,3)+transpt(2,2,3)+transpt(3,3,3))/3.d0
        seb0 = (transpt(1,1,4)+transpt(2,2,4)+transpt(3,3,4))/3.d0
        pf0 = cod0 * ((seb0/1d6)**2.d0) * 1d4
        pft0 = (pf0 / 1d4) * (eptemp(itemp) / kB)
        WRITE (100000,'(7x,f7.1,2x,f10.4,1x,15(2x,es10.2))') eptemp(itemp)/kB, ef_epw(itemp,idope)*ryd2ev, epdope(idope), &
                                                    n_elec(itemp,idope)+n_hole(itemp,idope), n_hole(itemp,idope)-n_elec(itemp,idope), &
                                                    mob0, cod0, seb0, pf0, pft0, etc0

        if ((itemp == neptemp) .and. (idope == nepdope)) then
           WRITE (100000,'(7x,a)') '-------------------------------------------------------'
           WRITE (100000,'(7x,a)')                        ''
           WRITE (100000,'(7x,a)')                        'MOB = Electron mobility [cm^2/V/s]'
           WRITE (100000,'(7x,a)')                        'COD = Electron conductivity [1/Ohm/m]'
           WRITE (100000,'(7x,a)')                        'ETC = Electron thermal conductivity [W/m/K]'
           WRITE (100000,'(7x,a)')                        'SEB = Seebeck coefficient [uV/K]'
           WRITE (100000,'(7x,a)')                        'PEL = Peltier coefficient [uV]'
           WRITE (100000,'(7x,a)')                        ''
        endif
        !
        close(100000)
     endif
     !
  ENDIF
  !
  !
END SUBROUTINE export_result 
