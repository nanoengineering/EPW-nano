!----------------------------------------------------------------------------
SUBROUTINE iter_bte_el (ik0, ik_pol, nk_pol, F_k_pol, itemp, idope, iter)
!----------------------------------------------------------------------------
#INCLUDE "f_defs.h"
  USE kinds,        ONLY : DP
  USE cell_base,    ONLY : at
  USE phcom,        ONLY : nmodes
  USE epwcom,       ONLY : nkf1, nkf2, nkf3, nqf1, nqf2, nqf3, fermi_energy, efermi_read, phdrag
  USE elph2,        ONLY : ibndmin, ibndmax, nbnd_red, tau0_k_ful, F_k_ful, vel_ful, etf_all, ef_epw, wf_all
  USE bte_var
  USE bte_func
#ifdef __PARA
  USE io_global,    ONLY : ionode_id
  USE mp,           ONLY : mp_barrier
  USE mp_global,    ONLY : my_pool_id
#ENDIF
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)        :: ik0, ik_pol, nk_pol, itemp, idope, iter
  REAL(KIND=DP), INTENT(OUT) :: F_k_pol(3,nbnd_red,nk_pol)
  !
  INTEGER                    :: iq_id(nq_ful)
  REAL(KIND=DP)              :: sigmai(nq_ful)
  REAL(KIND=DP)              :: Delta(3,nbnd_red), ef0
  ! id and index
  INTEGER                    :: ik, ik_irr_red, ik_ful, ikq_ful_red, iq_ful_red, iq_irr_red, ikq_irr_red, iq, iq_irr, ikq, &
                                ibnd, imode, jbnd, ibnd0, jbnd0, nscat, ns, isym
  INTEGER                    :: ijk_k(3), ijk_q(3), ijk_kq(3), ijk_q_rot(3), ijk_q_fbz(3)
  REAL(KIND=DP)              :: xq(3), xq_rot(3), F_tmp_kq, e_tmp_kq
  !
  !
  IF (efermi_read) THEN
     ef0 = fermi_energy
  ELSE
     ef0 = ef_epw(itemp,idope)
  ENDIF
  !
  !
  ! ik0 is the kpoint in red-ful-BZ ranging from ik_star to ik_stop (totally, from 1 to nk_ful_red)
  ! ik_pol = 1 ~ nk_pol, used in F_k_pol
  ! ik is the kpoint in red-ful-BZ sequenced by the number of scattering events
  ik = seq2nscat(ik0)
  ik_ful = rful2ful(ik)
  ik_irr_red = rful2rirr(ik)
  !
  ijk_k= id2ijk(ik_ful)
  !
  !
  OPEN (1111,FILE='BTE/META/irr2ful_sym',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='old')
  READ (1111,REC=ik_ful) isym
  CLOSE (1111)
  !
  !
  ! iteration
  F_k_pol(:,:,ik_pol) = 0.0d0 
  Delta = 0.0d0
  IF (iter .EQ. 0) GOTO 17120
  !
  DO ibnd = 1, nbnd_red
     DO jbnd = 1, nbnd_red 
        !
        jbnd0 = jbnd + ibndmin - 1
        !
        DO imode = 1, nmodes
           !
           sigmai = 0.0d0
           iq_id  = 0
           !
           CALL sigmai_load (ik_irr_red, ibnd, jbnd, imode, nscat, iq_id, sigmai, itemp, idope)
           !
           DO ns = 1, nscat 
              !
              iq = iq_id(ns)
              !
              ! crystal coordinate of q
              xq(:) = xqf_ful(:,iq)
              !
              ! 1. symmat_lat (S) is the rotation matrix used to rotate the coordinate of k from k_irr to k_ful: k_ful=S*k_irr.
              ! 2. sigmai is the electron self-energy on kq coordinate in which k is in the irreducible BZ.
              ! 3. In this subroutine, k is already rotated from k_ful to k_irr, that is, k_irr=P*k_ful where P=S^-1.
              ! 4. Therefore we have to rotate q with an inverse angle, which is q_rot=P^-1*q=S*q.
              ! 5. We can find k_irr+q_rot is equivalent to k_ful+q
              xq_rot(:) = xq(1)*symmat_lat(:,1,isym) + xq(2)*symmat_lat(:,2,isym) + xq(3)*symmat_lat(:,3,isym)
              !
              ! transform q to crystal index
              ijk_q_rot(1) = NINT(xq_rot(1)*DBLE(nqf1))
              ijk_q_rot(2) = NINT(xq_rot(2)*DBLE(nqf2))
              ijk_q_rot(3) = NINT(xq_rot(3)*DBLE(nqf3))
              !
              ! fold crystal index of kq in FBZ
              ijk_kq = ijk_fbz(ijk_k,ijk_q_rot)
              !
              ! find the id of kq
              ikq = ijk2id(ijk_kq) ! ijk2id
              !
              ikq_ful_red = ful2rful(ikq) ! ful-BZ to red-ful-BZ
              !
              IF (ikq_ful_red .NE. 0) THEN
                 !
                 ikq_irr_red = rful2rirr(ikq_ful_red)
                 Delta(:,ibnd) = Delta(:,ibnd) + ( sigmai(ns) * (etf_all(jbnd0,ikq_irr_red)-ef0) * F_k_ful(:,jbnd,ikq_ful_red) )
                 !
              ENDIF
              !
           ENDDO ! ns
           !
        ENDDO ! mode
     ENDDO !jbnd
  ENDDO ! ibnd
  !
  !
17120 CONTINUE ! if iter=0
  !
  DO ibnd = 1, nbnd_red
     !
     ibnd0 = ibnd+ibndmin-1
     !
     ! F_k_pol is in units of v*tau
     F_k_pol(:,ibnd,ik_pol) = tau0_k_ful(ibnd,ik) * ( vel_ful(:,ibnd,ik) + ( 2.0d0*Delta(:,ibnd) / (etf_all(ibnd0,ik_irr_red)-ef0) )  )
     !
     ! impose symmetry on every kpoint to prevent from drift
     CALL symmetrizer (F_k_pol(:,ibnd,ik_pol),ik_ful,'k')
     !
  ENDDO 
  !
  !
END SUBROUTINE iter_bte_el


!----------------------------------------------------------------------------
SUBROUTINE iter_bte_ph (iq0, itemp, idope, iter)
!----------------------------------------------------------------------------
#INCLUDE "f_defs.h"
  USE kinds,        ONLY : DP
  USE cell_base,    ONLY : at
  USE phcom,        ONLY : nmodes
  USE epwcom,       ONLY : nkf1, nkf2, nkf3, nqf1, nqf2, nqf3, fermi_energy, efermi_read, phdrag
  USE elph2,        ONLY : ibndmin, ibndmax, nbnd_red, F_k_ful, vel_ful, etf_all, ef_epw, wf_all, vph_ful
  USE bte_var
  USE bte_func
#ifdef __PARA
  USE io_global,    ONLY : ionode_id, stdout
  USE mp,           ONLY : mp_barrier
  USE mp_world,     ONLY : world_comm
  USE mp_global,    ONLY : my_pool_id
#ENDIF
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)        :: iq0, itemp, idope, iter
  !
  INTEGER                    :: ik_id(nk_ful)
  REAL(KIND=DP)              :: gammai(nk_ful)
  REAL(KIND=DP)              :: Delta(3,nmodes), ef0
  ! id and index
  INTEGER                    :: ik, iq_irr_red, iq_ful, ik_ful, ikq_ful_red, ikq_irr_red, iq_ful_red, ik_irr_red, ik_ful_red, iq_irr, &
                                iq, ikq, ibnd, imode, jbnd, ibnd0, jbnd0, nscat, ns, isym, isym_q
  INTEGER                    :: ijk_k(3), ijk_q(3), ijk_kq(3), ijk_k_rot(3), ijk_q_rot(3), ijk_k_fbz(3)
  REAL(KIND=DP)              :: xk(3), xq(3), xk_rot(3), xq_rot(3), F_tmp_k(3), F_tmp_kq(3), e_tmp_k, e_tmp_kq, &
                                phdrag_k_term(3,1), vqq(3,1), phdrag_k_mat(3,3), phdrag_k_mat_tmp(3,3)
  !
  CHARACTER(LEN=256)         :: file_ufmt
  CHARACTER(LEN=12)          :: tnpe
  CHARACTER(LEN=10)          :: ibnd_num, jbnd_num, itemp_num, idope_num
  LOGICAL                    :: exst
  INTEGER                    :: io
  !
  IF (efermi_read) THEN
     ef0 = fermi_energy
  ELSE
     ef0 = ef_epw(itemp,idope)
  ENDIF
  !
  !
  ! iq0 belongs to nq_ful_red
  iq_ful_red = iq0
  iq_ful     = rful2ful_q(iq_ful_red)
  iq_irr_red = rful2rirr_q(iq_ful_red)
  iq_irr     = rirr2irr_q(iq_irr_red)
  !
  ijk_q= id2ijk_q(iq_ful)
  !
  OPEN (1111,FILE='BTE/META/irr2ful_sym_q',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='old')
  READ (1111,REC=iq_ful) isym_q
  CLOSE (1111)
  !
  ! iteration
  Delta = 0.0d0
!  IF (iter .EQ. 0) GOTO 17126
  !
  DO imode = 1, nmodes
     !
     vqq(:,1) = vph_ful(:,imode,iq_ful_red)
     !
     DO ibnd = 1, nbnd_red
        !
        ibnd0 = ibnd + ibndmin - 1
        !
        WRITE(ibnd_num,'(i10)') ibnd0
        WRITE(itemp_num,'(i10)') itemp
        WRITE(idope_num,'(i10)') idope
        tnpe = 'T'//TRIM(ADJUSTL(itemp_num))//'_N'//TRIM(ADJUSTL(idope_num))//'_i'//TRIM(ADJUSTL(ibnd_num))
        !
        file_ufmt = 'BTE/META/phdragsumk_'//TRIM(ADJUSTL(tnpe))
        !
        DO jbnd = 1, nbnd_red 
           !
           jbnd0 = jbnd + ibndmin - 1
           !
           gammai = 0.0d0
           ik_id  = 0
           !
           CALL gammai_load (iq_irr_red, ibnd, jbnd, imode, nscat, ik_id, gammai, itemp, idope)
           !
           !
           DO ns = 1, nscat 
              !
              ! ik belongs to nk_ful
              ik = ik_id(ns) 
              !
              xk(:) = xkf_ful(:,ik)
              !
              xk_rot(:) = xk(1)*symmat_lat(:,1,isym_q) + xk(2)*symmat_lat(:,2,isym_q) + xk(3)*symmat_lat(:,3,isym_q)
              !
              ijk_k_rot(1) = NINT(xk_rot(1)*DBLE(nkf1))
              ijk_k_rot(2) = NINT(xk_rot(2)*DBLE(nkf2))
              ijk_k_rot(3) = NINT(xk_rot(3)*DBLE(nkf3))
              !
              ijk_kq = ijk_fbz(ijk_q,ijk_k_rot)
              !
              ikq = ijk2id(ijk_kq)
              !
              !
              ik_ful_red  = ful2rful(ik)
              IF (ik_ful_red .NE. 0) THEN
                 !
                 !ik_irr_red  = rful2rirr(ik_ful_red)
                 !e_tmp_k = etf_all(ibnd0,ik_irr_red) - ef0
                 F_tmp_k(:) = F_k_ful(:,ibnd,ik_ful_red)
                 !
              ELSE
                 !
                 F_tmp_k(:) = 0.0d0
                 !e_tmp_k = 0.0d0
                 !
              ENDIF
              !
              !
              ikq_ful_red = ful2rful(ikq)
              IF (ikq_ful_red .NE. 0) THEN
                 !
                 !ikq_irr_red = rful2rirr(ikq_ful_red)
                 !e_tmp_kq = etf_all(jbnd0,ikq_irr_red)-ef0
                 F_tmp_kq(:) = F_k_ful(:,jbnd,ikq_ful_red)
                 !
              ELSE
                 !
                 F_tmp_kq(:) = 0.0d0
                 !e_tmp_kq = 0.0d0
                 !
              ENDIF
              !
              !Delta(:,imode) = Delta(:,imode) + ( gammai(ns) * (e_tmp_k*F_tmp_k(:) - e_tmp_kq*F_tmp_kq(:)) )
              Delta(:,imode) = Delta(:,imode) + ( gammai(ns) * (F_tmp_k(:) - F_tmp_kq(:)) )
              !
              phdrag_k_term(:,1) = tau0_q_ful(itemp,imode,iq_ful_red)*( gammai(ns) * (F_tmp_k(:) - F_tmp_kq(:)) )
              phdrag_k_mat(1:3,1:3) = 2.0d0*MATMUL(vqq,TRANSPOSE(phdrag_k_term)) * wf_all(imode,iq_irr_red) 
              !
              INQUIRE (FILE=file_ufmt,EXIST=exst)
              phdrag_k_mat_tmp = 0.0d0
              !
              IF (exst) THEN
                OPEN (18888,FILE=file_ufmt,FORM='unformatted',ACCESS='direct',RECL=9*DP,STATUS='old')
                READ (18888,REC=ik,IOSTAT=io) phdrag_k_mat_tmp(1,:),phdrag_k_mat_tmp(2,:),phdrag_k_mat_tmp(3,:)
                IF ( io .NE. 0) phdrag_k_mat_tmp = 0.0d0
                phdrag_k_mat(:,:) = phdrag_k_mat(:,:) + phdrag_k_mat_tmp(:,:)
                WRITE (18888,REC=ik) phdrag_k_mat(1,:),phdrag_k_mat(2,:),phdrag_k_mat(3,:)
              !
              ELSE
                OPEN (18888,FILE=file_ufmt,FORM='unformatted',ACCESS='direct',RECL=9*DP,STATUS='replace')
                WRITE (18888,REC=ik) phdrag_k_mat(1,:),phdrag_k_mat(2,:),phdrag_k_mat(3,:)
              ENDIF
              !
              CLOSE (18888)

              !
           ENDDO ! ns
           !
        ENDDO ! jmode
        !
     ENDDO ! ibnd
     !
  ENDDO ! imode
  !
  !
!17126 CONTINUE ! if iter=0
  !
  DO imode = 1, nmodes
     !
     N_q_ful(:,itemp,imode,iq_ful_red) = tau0_q_ful(itemp,imode,iq_ful_red) * vph_ful(:,imode,iq_ful_red)
     !dN_q_ful(:,imode,iq_ful_red) = tau0_q_ful(imode,iq_ful_red) * ( 2.0d0*Delta(:,imode)/wf_all(imode,iq_irr_red) )
     dN_q_ful(:,itemp,imode,iq_ful_red) = tau0_q_ful(itemp,imode,iq_ful_red) * ( 2.0d0*Delta(:,imode) )
     !
     IF (ph_rate_ful(itemp,imode,iq_ful) .NE. 0.0d0) THEN
     int_N_q_ful(:,itemp,imode,iq_ful_red) = 1.0d0 * vph_ful(:,imode,iq_ful_red) / ph_rate_ful(itemp,imode,iq_ful)
     ELSE
     int_N_q_ful(:,itemp,imode,iq_ful_red) = 0.0d0
     ENDIF  
     !
     CALL symmetrizer (N_q_ful(:,itemp,imode,iq_ful_red),iq_ful,'q')
     CALL symmetrizer (dN_q_ful(:,itemp,imode,iq_ful_red),iq_ful,'q')
     !
  ENDDO 
  !
  !
END SUBROUTINE iter_bte_ph





!----------------------------------------------------------------------------
SUBROUTINE iter_tdbte_el (ik0, ik_pol, nk_pol, f_t_pol)
!----------------------------------------------------------------------------
#INCLUDE "f_defs.h"
  USE kinds,        ONLY : DP
  USE cell_base,    ONLY : at
  USE phcom,        ONLY : nmodes
  USE epwcom,       ONLY : nkf1, nkf2, nkf3, nqf1, nqf2, nqf3, eptemp, dt, efield, gradt, phdrag
  USE elph2,        ONLY : ibndmin, ibndmax, nbnd_red, &
                           f_t_ful, f_0_ful, vel_ful, wf_ful, wf_irr
  USE bte_var
  USE bte_func
#ifdef __PARA
  USE io_global,    ONLY : ionode_id
  USE mp,           ONLY : mp_barrier
  USE mp_global,    ONLY : my_pool_id
#ENDIF
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)        :: ik0, ik_pol, nk_pol
  REAL(KIND=DP), INTENT(OUT) :: f_t_pol(nbnd_red,nk_pol)
  !
  INTEGER                    :: iq_id(nq_ful)
  REAL(KIND=DP)              :: weight_abs(nq_ful), weight_emi(nq_ful)
  REAL(KIND=DP)              :: Delta(nbnd_red,nk_pol)
  ! id and index
  INTEGER                    :: ik, ik_irr_red, ik_ful, ikq_ful_red, iq, ikq, ibnd, imode, jbnd, nscat, ns, isym
  INTEGER                    :: ijk_k(3), ijk_q(3), ijk_kq(3), ijk_q_rot(3), ijk_q_fbz(3)
  REAL(KIND=DP)              :: xq(3), xq_rot(3)
  !
  REAL(KIND=DP)              :: f_k, f_kq, n_q, f_abs, f_emi, inv_eptemp0
  REAL(KIND=DP), EXTERNAL    :: wgauss
  !
  inv_eptemp0 = 1.0d0/eptemp(1)
  !
  ! ik0 is the kpoint in red-ful-BZ ranging from ik_star to ik_stop (totally, from 1 to nk_ful_red)
  ! ik_pol = 1 ~ nk_pol, used in f_t_pol
  ! ik is the kpoint in red-ful-BZ sequenced by the number of scattering events
  ik = seq2nscat(ik0)
  ik_ful = rful2ful(ik)
  !
  OPEN (1111,FILE='BTE/META/irr2ful_sym',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='old')
  READ (1111,REC=ik_ful) isym
  CLOSE (1111)
  !
  ! ik_irr_red is the kpoint in red-irr-BZ
  ik_irr_red = rful2rirr(ik)
  !
  f_t_pol(:,ik_pol) = 0.0d0 
  !
  ! ijk index of ik in ful-BZ
  ijk_k= id2ijk(ik_ful)
  !
  ! time evolution
  Delta = 0.0d0
  DO ibnd = 1, nbnd_red
     DO jbnd = 1, nbnd_red 
        DO imode = 1, nmodes
           !
           weight_abs = 0.0d0
           weight_emi = 0.0d0
           iq_id  = 0
           !
           CALL weight_load (ik_irr_red, ibnd, jbnd, imode, nscat, iq_id, weight_abs, weight_emi)
           !
           DO ns = 1, nscat 
              !
              iq = iq_id(ns)
              !
              ! crystal coordinate of q
              xq(:) = xqf_ful(:,iq)
              !
              ! 1. symmat_lat (S) is the rotation matrix used to rotate the coordinate of k from k_irr to k_ful: k_ful=S*k_irr.
              ! 2. sigmai is the electron self-energy on kq coordinate in which k is in the irreducible BZ.
              ! 3. In this subroutine, k is already rotated from k_ful to k_irr, that is, k_irr=P*k_ful where P=S^-1.
              ! 4. Therefore we have to rotate q with an inverse angle, which is q_rot=P^-1*q=S*q.
              ! 5. We can find k_irr+q_rot is equivalent to k_ful+q
              xq_rot(:) = xq(1)*symmat_lat(:,1,isym) + xq(2)*symmat_lat(:,2,isym) + xq(3)*symmat_lat(:,3,isym)
              !
              ! transform q to crystal index
              ijk_q_rot(1) = NINT(xq_rot(1)*DBLE(nqf1))
              ijk_q_rot(2) = NINT(xq_rot(2)*DBLE(nqf2))
              ijk_q_rot(3) = NINT(xq_rot(3)*DBLE(nqf3))
              !
              ! fold crystal index of kq in FBZ
              ! k and q-mesh should be the same
              ijk_kq = ijk_fbz(ijk_k,ijk_q_rot)
              !
              ! find the id of kq
              ikq = ijk2id(ijk_kq)
              !
              f_k  = f_0_ful(ibnd,ik)
              !
              f_kq = 0.0d0
              ikq_ful_red = ful2rful(ikq)
              IF (ikq_ful_red .NE. 0) f_kq = f_0_ful(jbnd,ikq_ful_red)
              !
              n_q = n_q/(1.0d0-2.0d0*n_q)    
              !
              f_abs = (1.0d0+n_q)*f_kq - (n_q+f_kq)*f_k ! (1.0d0-f_k)*(1.0d0+n_q)*f_kq - f_k*n_q*(1.0d0-f_kq)
              f_emi = (f_k+n_q)*f_kq - (1.0d0+n_q)*f_k  ! (1.0d0-f_k)*n_q*f_kq - f_k*(1.0d0+n_q)*(1.0d0-f_kq)
              !
              Delta(ibnd,ik_pol) = Delta(ibnd,ik_pol) + ( f_abs*weight_abs(ns) + f_emi*weight_emi(ns) )


             ! write (6,'(i6,7es15.6)') ns, f_k, f_kq, n_q, f_abs, f_emi, f_abs*weight_abs(ns), f_emi*weight_emi(ns)

              !
           ENDDO ! ns
           !
        ENDDO ! mode
     ENDDO !jbnd
  ENDDO ! ibnd
  !
  DO ibnd = 1, nbnd_red
     f_t_pol(ibnd,ik_pol) = f_0_ful(ibnd,ik) + dt*2.0d0*Delta(ibnd,ik_pol)
  ENDDO 
  !
  !
END SUBROUTINE iter_tdbte_el
