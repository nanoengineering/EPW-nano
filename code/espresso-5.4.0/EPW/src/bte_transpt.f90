!---------------------------------------------------------------------------------
SUBROUTINE bte_transpt (transpt, itemp, idope)
!---------------------------------------------------------------------------------
#INCLUDE "f_defs.h"
  USE kinds,         ONLY : DP
  USE pwcom,         ONLY : ef 
  USE phcom,         ONLY : nmodes
  USE cell_base,     ONLY : at, bg, omega
  USE spin_orb,      ONLY : lspinorb
  USE lsda_mod,      ONLY : nspin
  USE epwcom,        ONLY : bte, nbndsub, eptemp, efermi_read, fermi_energy, nptype, epdope, epdim, neptemp, phdrag
  USE elph2,         ONLY : ibndmin, ibndmax, nbnd_red, n_elec, n_hole, etf_all, vel_ful, ef_epw, vph_ful, wf_all,  &
                            F_k_ful, mfp_k_ful, mfp_q_ful, cbnd_emin, vbnd_emax, vbnd_num, cbnd_num
  USE bte_var
  USE bte_func
  USE constants_epw, ONLY : kB, au2Ohm, au2Amp, au2j, au2m, au2cm, au2s, twopi
#ifdef __PARA
  USE io_global,     ONLY : ionode_id
  USE mp_global,     ONLY : me_pool, inter_pool_comm, my_pool_id
#ENDIF
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)        :: itemp, idope
  REAL(KIND=DP), INTENT(OUT) :: transpt(3,3,24)
  ! locol vairable
  INTEGER                    :: ik_ful_red, ik_ful, ik_irr_red, ibnd, ibnd0, imode, ir, i, j, uorl, &
                                io, iq_ful_red, iq_irr_red
  REAL(KIND=DP)              :: wsp, temp, ef0, h_carr, e_carr, f0, df0, ekk, vkk(3,1), n0, dn0, wqq, vqq(3,1), &
                                vtau(3,1), vtau_q(3,1), int_vtau_q(3,1),vtau_q_phd(3,1)
  REAL(KIND=DP)              :: factor(15)  
  ! transport property
  REAL(KIND=DP)              :: L11(3,3), L11_e(3,3), L11_h(3,3), L12(3,3), L12_e(3,3), L12_h(3,3), &
                                L21(3,3), L21_e(3,3), L21_h(3,3), L22(3,3), L33(3,3), L3i(3,3), L44(3,3), inv_L11(3,3), &
                                inv_L11_h(3,3), inv_L11_e(3,3), L11_k(3,3),L11_e_k(3,3), L11_h_k(3,3), L12_k(3,3), &
                                L12_e_k(3,3), L12_h_k(3,3), L21_k(3,3), L21_e_k(3,3), L21_h_k(3,3), L22_k(3,3), &
                                L44_k(3,3), L33_q(3,3), L44_q(3,3), &
                                L44_e(3,3), L44_x(3,3), L44_l(3,3), L44_g(3,3), xful_cry(3), L3i_q(3,3)
  REAL(KIND=DP)              :: mob(3,3), cod(3,3), etc(3,3), seb(3,3), pel(3,3), pof(3,3), loz(3,3), ltc(3,3), sltc(3,3), &
                                btc(3,3),spd(3,3), spd_k(3,3), spd_x(3,3), spd_l(3,3), spd_g(3,3), phdrag_k_mat(3,3), &
                                mob_e(3,3), mob_h(3,3), cod_e(3,3), cod_h(3,3), seb_e(3,3), seb_h(3,3)
  !
  CHARACTER(LEN=256)         :: onsager_ufmt, phdrag_ufmt, phdrag_k_ufmt, tnpe, phdragsumk_ufmt
  CHARACTER(LEN=10)           :: itemp_num, idope_num, ibnd_num
  REAL(KIND=DP), EXTERNAL    :: wgauss
  !
  !
  ! initialization
  temp = eptemp(itemp) ! in units of Ry
  !
  IF (efermi_read) THEN
     ef0 = fermi_energy
  ELSE
     ef0 = ef_epw(itemp,idope)
  ENDIF
  !
  h_carr = ABS(n_hole(itemp,idope)) ! in units of cm^3
  IF (h_carr .LT. 1.0d-50) h_carr = 1.0d-50
  e_carr = ABS(n_elec(itemp,idope)) ! in units of cm^3
  IF (e_carr .LT. 1.0d-50) e_carr = 1.0d-50
  !
  ! check spin
  !IF (nspin .EQ. 1 .OR. nspin .EQ. 4) THEN
  IF (.NOT. lspinorb) THEN
     wsp = 2.0d0
  ELSE
     wsp = 1.0d0
  ENDIF
  !
  L11       = 0.0d0
  L11_e     = 0.0d0
  L11_h     = 0.0d0
  L12       = 0.0d0
  L12_e     = 0.0d0
  L12_h     = 0.0d0
  L21       = 0.0d0
  L21_e     = 0.0d0
  L21_h     = 0.0d0
  L22       = 0.0d0
  L33       = 0.0d0
  L3i       = 0.0d0
  L44       = 0.0d0
  L44_e     = 0.0d0
  L44_x     = 0.0d0
  L44_l     = 0.0d0
  L44_g     = 0.0d0
  L11_k     = 0.0d0
  L11_e_k   = 0.0d0
  L11_h_k   = 0.0d0
  L12_k     = 0.0d0
  L12_e_k   = 0.0d0
  L12_h_k   = 0.0d0
  L21_k     = 0.0d0
  L21_e_k   = 0.0d0
  L21_h_k   = 0.0d0
  L22_k     = 0.0d0
  L33_q     = 0.0d0
  L3i_q     = 0.0d0
  L44_q     = 0.0d0
  L44_k     = 0.0d0
  mob       = 0.0d0
  mob_e     = 0.0d0
  mob_h     = 0.0d0
  cod       = 0.0d0
  cod_e     = 0.0d0
  cod_h     = 0.0d0
  etc       = 0.0d0
  seb       = 0.0d0
  pel       = 0.0d0
  loz       = 0.0d0
  btc       = 0.0d0
  ltc       = 0.0d0
  sltc       = 0.0d0
  spd       = 0.0d0
  transpt   = 0.0d0
  mfp_k_ful = 0.0d0
  mfp_q_ful = 0.0d0
  !
  factor(1) = -1.0d0*(wsp/DBLE(nk_ful))*(2.0d0/omega) / au2Ohm / au2m          
  factor(2) = (wsp/DBLE(nk_ful))*(SQRT(2.0d0)/omega/temp) * kB * au2Amp / au2m
  factor(3) = (wsp/DBLE(nk_ful))*(SQRT(2.0d0)/omega) * au2Amp / au2m             
  factor(4) = -1.0d0*(wsp/DBLE(nk_ful))*(1.0d0/omega/temp) * kB * au2j / au2s / au2m ! au2m 05252017
  factor(5) = 1.0/omega/DBLE(nq_ful) * kB * au2j / au2s / au2m
  factor(6) = (2.0d0/DBLE(nq_ful))*(SQRT(2.0d0)/omega/temp/temp) * kB * au2Amp / au2m ! wsp and nk_ful have been considered in gammai_phd 
  !
  ! only last iteration will be preserved
  IF (my_pool_id .EQ. ionode_id) THEN
      !
      WRITE(itemp_num,'(i3)') itemp
      WRITE(idope_num,'(i3)') idope
      onsager_ufmt = 'BTE/META/onsager_T'//TRIM(ADJUSTL(itemp_num))//'_N'//TRIM(ADJUSTL(idope_num))
      OPEN (21212,FILE=onsager_ufmt,FORM='unformatted',ACCESS='direct',RECL=4+9*4*DP,STATUS='replace')
      phdrag_ufmt = 'BTE/META/phdrag_q_T'//TRIM(ADJUSTL(itemp_num))//'_N'//TRIM(ADJUSTL(idope_num))
      OPEN (22222,FILE=phdrag_ufmt,FORM='unformatted',ACCESS='direct',RECL=(1+2*9)*DP,STATUS='replace')
      phdrag_k_ufmt = 'BTE/META/phdrag_k_T'//TRIM(ADJUSTL(itemp_num))//'_N'//TRIM(ADJUSTL(idope_num))
      OPEN (33333,FILE=phdrag_k_ufmt,FORM='unformatted',ACCESS='direct',RECL=(1+9)*DP,STATUS='replace')
      !
  ENDIF
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
     DO ibnd = 1, nbnd_red
        !
        ibnd0 = ibnd+ibndmin-1
        !
        !
        ekk       = etf_all(ibnd0,ik_irr_red) - ef0
        vkk(:,1)  = vel_ful(:,ibnd,ik_ful_red)
        vtau(:,1) = F_k_ful(:,ibnd,ik_ful_red)
        !
        f0 = wgauss(-ekk/temp,-99)
        df0 = (-1.0d0/temp)*f0*(1.0d0-f0) ! 1.0 -> -1.0 05252017
        !
        ! Onsager coefficients
        IF (ibnd0 .LE. vbnd_num) THEN  !  to include semi-metal
           !
           L11_h_k(1:3,1:3) = factor(1) * df0 * MATMUL(vkk,TRANSPOSE(vtau)) ! f0_h = 1-f0_e, and df0_h = -df0_e
           L11_h = L11_h + L11_h_k
           !
           L12_h_k(1:3,1:3) = factor(2) * df0 * MATMUL(vkk,TRANSPOSE(vtau)) * ekk
           L12_h = L12_h + L12_h_k
           !
           L21_h_k(1:3,1:3) = factor(3) * df0 * MATMUL(vkk,TRANSPOSE(vtau)) * ekk
           L21_h = L21_h + L21_h_k
           !
         ELSE
           !
           L11_e_k(1:3,1:3) = factor(1) * df0 * MATMUL(vkk,TRANSPOSE(vtau))
           L11_e = L11_e + L11_e_k
           !
           L12_e_k(1:3,1:3) = factor(2) * df0 * MATMUL(vkk,TRANSPOSE(vtau)) * ekk
           L12_e = L12_e + L12_e_k
           !
           L21_e_k(1:3,1:3) = factor(3) * df0 * MATMUL(vkk,TRANSPOSE(vtau)) * ekk
           L21_e = L21_e + L21_e_k
           !
        ENDIF
           !
        L11_k(1:3,1:3) = factor(1) * df0 * MATMUL(vkk,TRANSPOSE(vtau))
        L12_k(1:3,1:3) = factor(2) * df0 * MATMUL(vkk,TRANSPOSE(vtau)) * ekk
        L21_k(1:3,1:3) = factor(3) * df0 * MATMUL(vkk,TRANSPOSE(vtau)) * ekk
        L22_k(1:3,1:3) = factor(4) * df0 * MATMUL(vkk,TRANSPOSE(vtau)) * ekk * ekk
        !
        L11 = L11 + L11_k
        L12 = L12 + L12_k
        L21 = L21 + L21_k
        L22 = L22 + L22_k
        !  
        ! mean free path of each k
        IF (MAXVAL(ABS(vkk(:,1))) .NE. 0.0d0) THEN
           mfp_k_ful(ibnd,ik_ful_red) = DOT_PRODUCT(vtau(:,1),vkk(:,1)) / SQRT(DOT_PRODUCT(vkk(:,1),vkk(:,1)))
        ELSE
           mfp_k_ful(ibnd,ik_ful_red) = 0.0d0
        ENDIF
        !
        ! save Onsager coefficients of each k
        IF (my_pool_id .EQ. ionode_id) THEN 
           !
           IF (ekk .GE. 0.0d0) THEN
              uorl = +1 ! upper Fermi level
           ELSE
              uorl = -1 ! lower Fermi level
           ENDIF
           !
           WRITE (21212,REC=(ik_ful_red-1)*nbnd_red+ibnd) L11_k(1:3,1:3), L12_k(1:3,1:3), L21_k(1:3,1:3), L22_k(1:3,1:3), uorl
           !
        ENDIF
        !
     ENDDO
  ENDDO
  !
  !
  IF (phdrag) THEN
     !
     WRITE(itemp_num,'(i10)') itemp
     !
     WRITE(idope_num,'(i10)') idope
     !
     OPEN (9999,FILE='BTE/META/xkf_ful_cryst',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='old')
     !
     DO ik_ful_red = 1, nk_ful_red
        !
        ik_ful = rful2ful(ik_ful_red)
        ik_irr_red = rful2rirr(ik_ful_red)
        READ (9999,REC=ik_ful) xful_cry(1:3)
        !
       DO ibnd = 1, nbnd_red
        !
        ibnd0 = ibnd+ibndmin-1
        WRITE(ibnd_num,'(i10)') ibnd0
        !
        tnpe = 'T'//TRIM(ADJUSTL(itemp_num))//'_N'//TRIM(ADJUSTL(idope_num))//'_i'//TRIM(ADJUSTL(ibnd_num))
        !
        phdragsumk_ufmt = 'BTE/META/phdragsumk_'//TRIM(ADJUSTL(tnpe))
        !
        OPEN (18888,FILE=phdragsumk_ufmt,FORM='unformatted',ACCESS='direct',RECL=9*DP,STATUS='old')
        !
        READ (18888,REC=ik_ful,IOSTAT=io) phdrag_k_mat(1,:),phdrag_k_mat(2,:),phdrag_k_mat(3,:)
        !
        IF (io .NE. 0) phdrag_k_mat = 0.0d0 ! some ik_ful might not have corresponding nscat written
        CLOSE (18888)
        !
        L44_k(1:3,1:3) = factor(6) * phdrag_k_mat(1:3,1:3)
        !
        L44_e =  L44_e + L44_k
        !
        IF ((SQRT((xful_cry(1)-0.5d0)**2 + (xful_cry(2)-0.5d0)**2+ (xful_cry(3)-0.5d0)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-0.5d0)**2 + xful_cry(2)**2+ xful_cry(3)**2) .LT. 0.0625) .OR. (SQRT(xful_cry(1)**2 + (xful_cry(2)-0.5d0)**2+ xful_cry(3)**2) .LT. 0.0625) .OR. (SQRT(xful_cry(1)**2 + xful_cry(2)**2+ (xful_cry(3)-0.5d0)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-0.5d0)**2 + (xful_cry(2)-1.0d0)**2+ (xful_cry(3)-1.0d0)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-1.0d0)**2 + (xful_cry(2)-0.5d0)**2+ (xful_cry(3)-1.0d0)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-1.0d0)**2 + (xful_cry(2)-1.0d0)**2+ (xful_cry(3)-0.5d0)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-0.5d0)**2 + (xful_cry(2)-1.0d0)**2+ xful_cry(3)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-1.0d0)**2 + (xful_cry(2)-0.5d0)**2+ xful_cry(3)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-1.0d0)**2 + xful_cry(2)**2+ (xful_cry(3)-0.5d0)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-0.5d0)**2 + xful_cry(2)**2+ (xful_cry(3)-1.0d0)**2) .LT. 0.0625) .OR. (SQRT(xful_cry(1)**2 + (xful_cry(2)-0.5d0)**2+ (xful_cry(3)-1.0d0)**2) .LT. 0.0625) .OR. (SQRT(xful_cry(1)**2 + (xful_cry(2)-1.0d0)**2+ (xful_cry(3)-0.5d0)**2) .LT. 0.0625)) THEN  ! crystal: L 0.5 0.5 0.5, 0/1 0/1 0.5, 0/1 0.5 0/1, 0.5 0/1 0/1
        !
        L44_l =  L44_l + L44_k
        !
        ELSEIF ((SQRT((xful_cry(1)-0.425d0)**2 + (xful_cry(2)-0.425d0)**2+ xful_cry(3)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-0.425d0)**2 + (xful_cry(3)-0.425d0)**2+ xful_cry(2)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(3)-0.425d0)**2 + (xful_cry(2)-0.425d0)**2+ xful_cry(1)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-0.425d0)**2 + (xful_cry(2)-0.425d0)**2+ (xful_cry(3)-1.0d0)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-0.425d0)**2 + (xful_cry(3)-0.425d0)**2+ (xful_cry(2)-1.0d0)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(3)-0.425d0)**2 + (xful_cry(2)-0.425d0)**2+ (xful_cry(1)-1.0d0)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-0.575d0)**2 + (xful_cry(2)-0.575d0)**2+ xful_cry(3)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-0.575d0)**2 + (xful_cry(3)-0.575d0)**2+ xful_cry(2)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(3)-0.575d0)**2 + (xful_cry(2)-0.575d0)**2+ xful_cry(1)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-0.575d0)**2 + (xful_cry(2)-0.575d0)**2+ (xful_cry(3)-1.0d0)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-0.575d0)**2 + (xful_cry(3)-0.575d0)**2+ (xful_cry(2)-1.0d0)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(3)-0.575d0)**2 + (xful_cry(2)-0.575d0)**2+ (xful_cry(1)-1.0d0)**2) .LT. 0.0625)) THEN   ! Si crystal: near X 0/1 0.425 0.425 or 0/1 0.575 0.575
        !
        L44_x =  L44_x + L44_k
        !
        ELSEIF ((SQRT(xful_cry(1)**2+ xful_cry(2)**2 + xful_cry(3)**2) .LT. 0.0625) .OR. (SQRT((xful_cry(1)-1.0d0)**2+ (xful_cry(2)-1.0d0)**2 + (xful_cry(3)-1.0d0)**2) .LT. 0.0625)) THEN  !gamma
        !
        L44_g =  L44_g + L44_k
        !
        ENDIF
        !
        ! save phdrag components of each k
        IF (my_pool_id .EQ. ionode_id) THEN 
           !
           WRITE (33333,REC=(ik_ful_red-1)*nbnd_red+ibnd)  etf_all(ibnd0,ik_irr_red), L44_k(1:3,1:3)
        !
        ENDIF
        !
       ENDDO
     ENDDO
     CLOSE (9999)
     !
     DO iq_ful_red = 1, nq_ful_red
        !
        iq_irr_red = rful2rirr_q(iq_ful_red)
        !
        DO imode = 1, nmodes
           !
           wqq             = wf_all(imode,iq_irr_red)
           vqq(:,1)        = vph_ful(:,imode,iq_ful_red)
           int_vtau_q(:,1) = int_N_q_ful(:,itemp,imode,iq_ful_red)
           vtau_q(:,1)     = N_q_ful(:,itemp,imode,iq_ful_red)
           vtau_q_phd(:,1) = dN_q_ful(:,itemp,imode,iq_ful_red)
           !
           IF (wqq .NE. 0.0d0) THEN
             n0 = wgauss(-wqq/temp,-99)
             n0 = n0/(1.0d0-2.0d0*n0)
             dn0 = (wqq/temp/temp)*n0*(n0+1.0d0) 
           ELSE
             n0  = 0.0d0
             dn0 = 0.0d0
           ENDIF
           !   
           L33_q(1:3,1:3) = factor(5) * dn0 * MATMUL(vqq,TRANSPOSE(vtau_q)) * wqq
           L3i_q(1:3,1:3) = factor(5) * dn0 * MATMUL(vqq,TRANSPOSE(int_vtau_q)) * wqq
           L44_q(1:3,1:3) = factor(6) * MATMUL(vqq,TRANSPOSE(vtau_q_phd)) * wqq
           !
           L33 = L33 + L33_q
           L3i = L3i + L3i_q
           L44 = L44 + L44_q
           !
           ! mean free path of each q
        IF (MAXVAL(ABS(vqq(:,1))) .NE. 0.0d0) THEN
           mfp_q_ful(imode,iq_ful_red) = DOT_PRODUCT(vtau_q(:,1),vqq(:,1)) / SQRT(DOT_PRODUCT(vqq(:,1),vqq(:,1)))
        ELSE
           mfp_q_ful(imode,iq_ful_red) = 0.0d0
        ENDIF
        !
        ! save phdrag components of each q
        IF (my_pool_id .EQ. ionode_id) THEN 
           !
           !
           WRITE (22222,REC=(iq_ful_red-1)*nmodes+imode)  wf_all(imode,iq_irr_red), L44_q(1:3,1:3), L33_q(1:3,1:3)
        !
        ENDIF
        !
       ENDDO ! imodes
    ENDDO ! iq_ful_red
     !
  ENDIF
  !
  !
  IF (epdim .EQ. 2) THEN
     !
     L11(3,3)   = 1.0d0
     L11_h(3,3) = 1.0d0
     L11_e(3,3) = 1.0d0
     L12(3,3)   = 1.0d0
     L12_h(3,3) = 1.0d0
     L12_e(3,3) = 1.0d0
     L21(3,3)   = 1.0d0
     L21_h(3,3) = 1.0d0
     L21_e(3,3) = 1.0d0
     L22(3,3) = 1.0d0
     L33(3,3) = 1.0d0
     L3i(3,3) = 1.0d0
     L44(3,3) = 1.0d0
     !
  ENDIF
  !
  ! compute electron transport properties
  inv_L11 = imat(L11)
  !
  IF (nptype .NE. 'p' .AND. nptype .NE. 'n') THEN
     !
    IF ((L11_h(1,1)*L11_h(2,2)*L11_h(3,3) + L11_h(1,2)*L11_h(2,3)*L11_h(3,1) + L11_h(1,3)*L11_h(3,2)*L11_h(2,1) - L11_h(1,3)*L11_h(2,2)*L11_h(3,1) - L11_h(2,3)*L11_h(3,2)*L11_h(1,1) - L11_h(3,3)*L11_h(1,2)*L11_h(2,1)) .NE. 0.0d0) THEN
     inv_L11_h = imat(L11_h)
    ELSE
     inv_L11_h = 0.0d0
    ENDIF
    IF ((L11_e(1,1)*L11_e(2,2)*L11_e(3,3) + L11_e(1,2)*L11_e(2,3)*L11_e(3,1) + L11_e(1,3)*L11_e(3,2)*L11_e(2,1) - L11_e(1,3)*L11_e(2,2)*L11_e(3,1) - L11_e(2,3)*L11_e(3,2)*L11_e(1,1) - L11_e(3,3)*L11_e(1,2)*L11_e(2,1)) .NE. 0.0d0) THEN
     inv_L11_e = imat(L11_e)
    ELSE
     inv_L11_e = 0.0d0
    ENDIF
     !
  ELSEIF (nptype .EQ. 'p') THEN
     !L33
     inv_L11_h = imat(L11_h)
     inv_L11_e = 0.0d0
     !
  ELSEIF (nptype .EQ. 'n') THEN
     !
     inv_L11_h = 0.0d0
     inv_L11_e = imat(L11_e)
     !
  ENDIF
  !
  cod_h = L11_h                                      ! [1/Ohm/m]
  cod_e = L11_e 
  cod = L11    
  !
  IF (h_carr .NE. 0.0d0) THEN
     mob_h = L11_h/(h_carr) / 1.60217662089d-19 / 100.0d0      ! [cm^2/V/s]
  ELSE
     mob_h = 0.0d0
  ENDIF
  !
  IF (e_carr .NE. 0.0d0) THEN
     mob_e = L11_e/(e_carr) / 1.60217662089d-19 / 100.0d0
  ELSE
     mob_e = 0.0d0
  ENDIF
  !
  mob   = (h_carr*mob_h+e_carr*mob_e)/(h_carr+e_carr) !
 ! mob = (L11/ndoping/1.602176565d-19) / 100.0d0 ! [cm^2/V/s]   ! -1.0 -> 1.0 05252017
 ! cod = L11                                     ! [1/Ohm/m]    ! -1.0 -> 1.0 05252017
  etc = L22 - MATMUL(MATMUL(L21,inv_L11),L12)   ! [W/m/K]
  seb = MATMUL(inv_L11,L12) * 1.0d+6            ! [uV/K]
  seb_e = MATMUL(inv_L11_e,L12_e) * 1.0d+6                ! [uV/K]
  seb_h = MATMUL(inv_L11_h,L12_h) * 1.0d+6                ! [uV/K]
  pel = MATMUL(L21,inv_L11) * 1.0d+6            ! [uV]
  pof = MATMUL(MATMUL(seb,seb),cod) * 1.0d-8    ! [uW/cm/K^2]
  loz = MATMUL(etc,inv_L11)/temp * kB           ! [W*Ohm/K^2]
  sltc = L3i                                     ! [W/m/K]
  ltc = L33
  btc = MATMUL(MATMUL(MATMUL(MATMUL(L21_e,inv_L11_e),L11_h)-L21_h,inv_L11),L12_e) + &
        MATMUL(MATMUL(MATMUL(MATMUL(L21_h,inv_L11_h),L11_e)-L21_e,inv_L11),L12_h) ! [W/m/K]
  spd = MATMUL(inv_L11,L44) * 1.0d+6            ! [uV/K]
  spd_k = MATMUL(inv_L11,L44_e) * 1.0d+6            ! [uV/K]
  spd_x = MATMUL(inv_L11,L44_x) * 1.0d+6            ! [uV/K]
  spd_l = MATMUL(inv_L11,L44_l) * 1.0d+6            ! [uV/K]
  spd_g = MATMUL(inv_L11,L44_g) * 1.0d+6            ! [uV/K]
  !
  transpt(:,:,1) = mob(:,:)
  transpt(:,:,16)  = mob_e(:,:)
  transpt(:,:,17)  = mob_h(:,:)
  transpt(:,:,2) = cod(:,:)
  transpt(:,:,18)  = cod_e(:,:)
  transpt(:,:,19)  = cod_h(:,:)
  transpt(:,:,3) = etc(:,:)
  transpt(:,:,4) = seb(:,:)
  transpt(:,:,20) = seb_e(:,:)
  transpt(:,:,21) = seb_h(:,:)
  transpt(:,:,5) = pel(:,:)
  transpt(:,:,6) = pof(:,:)
  transpt(:,:,7) = loz(:,:)
  transpt(:,:,8) = ltc(:,:)
  transpt(:,:,9) = sltc(:,:)
  transpt(:,:,10) = btc(:,:)
  transpt(:,:,11) = spd(:,:)
  transpt(:,:,12) = spd_k(:,:)
  transpt(:,:,13) = spd_x(:,:)
  transpt(:,:,14) = spd_l(:,:)
  transpt(:,:,15) = spd_g(:,:)
  !
  IF (epdim .EQ. 2) transpt(3,3,:) = 0.0d0
  !
  IF (my_pool_id .EQ. ionode_id) CLOSE (21212)
  IF (my_pool_id .EQ. ionode_id) CLOSE (22222)
  IF (my_pool_id .EQ. ionode_id) CLOSE (33333)
  !
END SUBROUTINE bte_transpt
