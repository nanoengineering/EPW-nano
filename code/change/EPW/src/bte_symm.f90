!----------------------------------------------------------------------------
SUBROUTINE bz_index (nkf1, nkf2, nkf3, flag)
!----------------------------------------------------------------------------
#INCLUDE "f_defs.h"
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : at, bg
  USE symm_base, ONLY : set_sym_bl, s, t_rev, time_reversal, nrot
  USE bte_var
  USE bte_func
#ifdef __PARA
  USE mp,        ONLY : mp_barrier, mp_sum, mp_min, mp_bcast
  USE io_global, ONLY : ionode_id, stdout
  USE mp_global, ONLY : my_pool_id, inter_pool_comm
  USE mp_world,      ONLY : world_comm
#ENDIF
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)          :: nkf1, nkf2, nkf3
  CHARACTER(LEN=*), INTENT(IN) :: flag
  ! flag = 'k' or 'q'
  ! 
  ! generate full and irreducible kpoints
  REAL(KIND=DP), PARAMETER   :: eps = 1.0d-5
  REAL(KIND=DP)              :: xk_ful(3), xk_irr(3), xk_rot(3), xk(3), xx, yy, zz, wk
  LOGICAL                    :: in_the_list, find_kpt
  ! full index to irreducible index
  INTEGER, ALLOCATABLE       :: equiv(:)
  INTEGER                    :: mini
  ! local variables
  REAL(KIND=DP), ALLOCATABLE :: wkf_ful(:)
  REAL(KIND=DP)              :: wk_tot
  INTEGER                    :: nk, nr, n, i, j, k, cnt
  REAL(KIND=DP)              :: t0, t1, t2, t3
  !
  INTEGER                    :: ik, ik0, ibnd, ir, isym, is, it, ik_star, ik_stop, nk_pol
  REAL(KIND=DP)              :: inv_at(3,3), rot_tmp(3,3), at_rot(3,3), inv_rot(3,3,48), xk_tmp(3), symmat(3,3)
  ! check xkf_irr
  INTEGER                    :: nk_tmp
  ! xkf_fbz
  INTEGER                    :: ig1, ig2, ig3, ic1, ic2, ic3
  REAL(KIND=DP)              :: xk_g(3), xk_fbz(3), leng, leng_new
  !
#ifdef __PARA
  REAL(KIND=DP)              :: tmp
#ENDIF 
  !
  !
  !
  IF (flag .NE. 'k' .AND. flag .NE. 'q') CALL errore ('bz_index','Wrong input flag',1)
  !
  WRITE (stdout,'(5x,a,a1,a/)') 'Prepare the information for ', flag, '-mesh'
  CALL CPU_TIME (t0)
  !
  !------------------------------------------------
  ! 1. nk_ful
  !------------------------------------------------
  !
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     nk_ful = nkf1*nkf2*nkf3
     !
     IF (flag .EQ. 'k') OPEN (9999,FILE='BTE/META/nk_ful',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='replace')
     IF (flag .EQ. 'q') OPEN (9999,FILE='BTE/META/nq_ful',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='replace')
     WRITE (9999,REC=1) nk_ful
     CLOSE (9999)
     !
  ENDIF
  !
#ifdef __PARA
  CALL mp_bcast (nk_ful,ionode_id,inter_pool_comm)
  CALL mp_barrier (inter_pool_comm)
#ENDIF
  !
  !
  !------------------------------------------------
  ! 2. xkf_ful
  !------------------------------------------------
  !
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     IF (flag .EQ. 'k') OPEN (8888,FILE='BTE/META/xkf_ful_cryst',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='replace')
     IF (flag .EQ. 'q') OPEN (8888,FILE='BTE/META/xqf_ful_cryst',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='replace')
     !
     DO i = 1, nkf1
        DO j = 1, nkf2
           DO k = 1, nkf3
              !
              nk = (i-1)*nkf2*nkf3 + (j-1)*nkf3 + k
              !
              WRITE (8888,REC=nk) DBLE(i-1)/DBLE(nkf1), DBLE(j-1)/DBLE(nkf2), DBLE(k-1)/DBLE(nkf3)
              !
           ENDDO
        ENDDO
     ENDDO
     !
     CLOSE (8888)
     !
  ENDIF
  !
#ifdef __PARA
  CALL mp_barrier (inter_pool_comm)
#ENDIF
  !
  !
  !------------------------------------------------
  ! 3. symmat_lat, n_rot (t_rev, time_reversal)
  !------------------------------------------------
  ALLOCATE (symmat_lat(3,3,48))
  symmat_lat = 0.0d0
  !
  CALL set_sym_bl ()
  !
  symmat_lat = DBLE(s)
  n_rot = nrot
  !
  OPEN (9999,FILE='BTE/META/symmat_lat',FORM='unformatted',ACCESS='direct',RECL=9*DP,STATUS='replace')
  DO isym = 1, n_rot
     WRITE (9999,REC=isym) symmat_lat(1,1:3,isym), symmat_lat(2,1:3,isym), symmat_lat(3,1:3,isym)
  ENDDO
  CLOSE (9999)
  !
  OPEN (9999,FILE='BTE/META/n_rot',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='replace')
  WRITE (9999,REC=1) n_rot
  CLOSE (9999)
  !
#ifdef __PARA
  CALL mp_barrier (inter_pool_comm)
#ENDIF
  !
  !
  !------------------------------------------------
  ! 3. nk_irr, xkf_irr, wkf_irr
  !------------------------------------------------
  !
  ALLOCATE (equiv(nk_ful))      
  equiv = 0
  !
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     ALLOCATE (wkf_ful(nk_ful)) 
     wkf_ful = 0.0
     !
     DO nk = 1, nk_ful
        equiv(nk) = nk
     ENDDO
     !
     IF (flag .EQ. 'k') OPEN (1111,FILE='BTE/META/xkf_ful_cryst',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='old')
     IF (flag .EQ. 'q') OPEN (1111,FILE='BTE/META/xqf_ful_cryst',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='old')
     !
     DO nk = 1, nk_ful
        !
        READ (1111,REC=nk) xk_ful(1:3)
        !
        IF (equiv(nk) .EQ. nk) THEN
           wkf_ful(nk) = 1.0d0
           DO nr = 1, n_rot
              DO i = 1, 3
                 xk_rot(i) = DBLE(s(i,1,nr))*xk_ful(1) + DBLE(s(i,2,nr))*xk_ful(2) + DBLE(s(i,3,nr))*xk_ful(3)
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
                 IF ((n .GT. nk) .AND. (equiv(n) .EQ. n)) then
                    equiv(n) = nk
                    wkf_ful(nk) = wkf_ful(nk) + 1.0d0
                 ELSE
                    IF ((equiv(n) .NE. nk) .OR. (n .LT. nk)) CALL errore ('bz_index', 'something wrong in the checking algorithm',1)
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
                    IF ((n .GT. nk) .AND. (equiv(n) .EQ. n)) then
                       equiv(n) = nk
                       wkf_ful(nk) = wkf_ful(nk) + 1.0d0
                    ELSE
                       IF ((equiv(n) .NE. nk) .OR. (n .LT. nk)) CALL errore ('bz_index', 'something wrong in the checking algorithm',2)
                    ENDIF
                 ENDIF
              ENDIF
           ENDDO
        ENDIF
     ENDDO
     CLOSE (1111)
     !
     !
     nk_irr = 0
     DO nk = 1, nk_ful
        IF (equiv(nk) .EQ. nk) nk_irr = nk_irr + 1
     ENDDO
     !
     IF (flag .EQ. 'k') OPEN (9999,FILE='BTE/META/nk_irr',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='replace')
     IF (flag .EQ. 'q') OPEN (9999,FILE='BTE/META/nq_irr',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='replace')
     WRITE (9999,REC=1) nk_irr
     CLOSE (9999)
     !
     IF (flag .EQ. 'k') OPEN (1111,FILE='BTE/META/xkf_ful_cryst',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='old')
     IF (flag .EQ. 'q') OPEN (1111,FILE='BTE/META/xqf_ful_cryst',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='old')
     IF (flag .EQ. 'k') OPEN (9999,FILE='BTE/META/xkf_irr_cryst',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='replace')
     IF (flag .EQ. 'q') OPEN (9999,FILE='BTE/META/xqf_irr_cryst',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='replace')
     !
     !
     ALLOCATE (wkf_irr(nk_irr)) 
     wkf_irr = 0.0d0
     !
     n = 0
     wk_tot = 0.0d0
     DO nk = 1, nk_ful
        IF (equiv(nk) .EQ. nk) THEN
           n = n + 1
           wkf_irr(n) = wkf_ful(nk)
           wk_tot = wk_tot + wkf_irr(n)
           ! 
           READ (1111,REC=nk) xk_ful(1:3)
           WRITE (9999,REC=n) xk_ful(1)-NINT(xk_ful(1)), xk_ful(2)-NINT(xk_ful(2)), xk_ful(3)-NINT(xk_ful(3))
           !
        ENDIF
     ENDDO
     CLOSE (1111)
     CLOSE (9999)
     !
     !  normalize weights to one
     IF (flag .EQ. 'k') OPEN (9999,FILE='BTE/META/wkf_irr',FORM='unformatted',ACCESS='direct',RECL=DP,STATUS='replace')
     IF (flag .EQ. 'q') OPEN (9999,FILE='BTE/META/wqf_irr',FORM='unformatted',ACCESS='direct',RECL=DP,STATUS='replace')
     DO nk = 1, nk_irr
        WRITE (9999,REC=nk) wkf_irr(nk) / wk_tot
     ENDDO
     CLOSE (9999)
     !
     !
     DEALLOCATE (wkf_irr)
     DEALLOCATE (wkf_ful)
     !
  ENDIF
  !
#ifdef __PARA
  CALL mp_bcast (nk_irr,ionode_id,inter_pool_comm)
  CALL mp_bcast (equiv,ionode_id,inter_pool_comm)
  CALL mp_barrier (inter_pool_comm)
#ENDIF
  !
  !
  !------------------------------------------------
  ! 4. ful2irr
  !------------------------------------------------
  ! explaination:
  !   xkf_ful(:,ik) and xkf_irr(:,ful2irr(ik)) are degenerate kpoint
  !   "ik" is in full BZ
  ! note:
  !   equiv points the degeneracies in full list
  !   ful2irr points the degeneracies in irreducible list
  !
  ALLOCATE (ful2irr(nk_ful))
  ful2irr = 0
  !
  ik = 0
  DO nk = 1, nk_ful
     IF (equiv(nk) .EQ. nk) THEN
        ik = ik + 1
        ful2irr(nk) = ik
     ENDIF
  ENDDO
  !
  DO nk = 1, nk_ful
     ful2irr(nk) = ful2irr(equiv(nk))
  ENDDO
  !
  ! export for checking
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     IF (flag .EQ. 'k') OPEN (1217,FILE='BTE/EPCHECK/equiv_k.dat',STATUS='replace')
     IF (flag .EQ. 'q') OPEN (1217,FILE='BTE/EPCHECK/equiv_q.dat',STATUS='replace')
     !
     DO nk = 1, nk_ful
        WRITE (1217,'(3i12)') nk, equiv(nk), ful2irr(nk)
     ENDDO
     CLOSE (1217)
     !
  ENDIF
  !
  !
#ifdef __PARA
  CALL mp_barrier (inter_pool_comm)
#ENDIF
  !
  !
  ! check, very very time-comsuming
!  DO n = 1, nk_irr
!     DO nk = 1, nk_ful
!        IF (ful2irr(nk) .EQ. n) GOTO 201
!     ENDDO
!     CALL errore ('bz_index','cannot remap the index from fine mesh to coarse mesh',n)
!201  CONTINUE
!  ENDDO
  !
  !
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     IF (flag .EQ. 'k') OPEN (8888,FILE='BTE/META/equiv',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='replace')
     IF (flag .EQ. 'q') OPEN (8888,FILE='BTE/META/equiv_q',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='replace')
     IF (flag .EQ. 'k') OPEN (9999,FILE='BTE/META/ful2irr',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='replace')
     IF (flag .EQ. 'q') OPEN (9999,FILE='BTE/META/ful2irr_q',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='replace')
     DO ik = 1, nk_ful
        WRITE (8888,REC=ik) equiv(ik)
        WRITE (9999,REC=ik) ful2irr(ik)
     ENDDO
     CLOSE (9999)
     !
  ENDIF
  !
  !
#ifdef __PARA
  CALL mp_barrier (inter_pool_comm)
#ENDIF
  !
  !
  !------------------------------------------------
  ! 5. irr2ful
  !------------------------------------------------
  !
  IF (my_pool_id .EQ. ionode_id) THEN
     ! 
     IF (flag .EQ. 'k') OPEN (9999,FILE='BTE/META/irr2ful',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='replace')
     IF (flag .EQ. 'q') OPEN (9999,FILE='BTE/META/irr2ful_q',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='replace')
     !
     ik = 0
     DO nk = 1, nk_ful
        IF (equiv(nk) .EQ. nk) THEN
           ik = ik + 1
           WRITE (9999,REC=ik) equiv(nk)
        ENDIF
     ENDDO    
     !
     CLOSE (9999)
     !
  ENDIF
  !
  DEALLOCATE (equiv)   
  !
#ifdef __PARA
  CALL mp_barrier (inter_pool_comm)
#ENDIF
  !
  !
  !------------------------------------------------
  ! 6. irr2ful_sym
  !------------------------------------------------
  ! explaination:
  !   If you got "246 63" in ful2irr.dat and "246 10" in irr2ful_sym.dat, this means the kpoint 246 in full BZ 
  !   can be obtained from kpoint 63 in irreducible BZ transformed by 10th transform matrix
  ! example:
  !   xkf_ful(:,ik) = xkf_irr(:,ful2irr(ik)) * s(:,:,10)
  !   "ik" is in full BZ and "xkf" is in crystal coordinate
  !
  IF (flag .EQ. 'k') OPEN (1111,FILE='BTE/META/xkf_ful_cryst',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='old')
  IF (flag .EQ. 'q') OPEN (1111,FILE='BTE/META/xqf_ful_cryst',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='old')
  IF (flag .EQ. 'k') OPEN (2222,FILE='BTE/META/xkf_irr_cryst',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='old')
  IF (flag .EQ. 'q') OPEN (2222,FILE='BTE/META/xqf_irr_cryst',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='old')
  !  
  ALLOCATE (irr2ful_sym(nk_ful))
  irr2ful_sym = 0
  !
  CALL para_bounds (ik_star, ik_stop, nk_ful)
  !
  DO ik = ik_star, ik_stop
     !
     READ (1111,REC=ik) xk_ful(1), xk_ful(2), xk_ful(3)
     READ (2222,REC=ful2irr(ik)) xk_irr(1), xk_irr(2), xk_irr(3)
     !
     IF (equiv_k(xk_ful,xk_irr)) THEN
        ! s(:,:,1) should be unit matrix
        irr2ful_sym(ik) = 1
        GOTO 404
     ENDIF
     !
     DO isym = 1, n_rot
        !
        xk_tmp(:) = xk_irr(1)*DBLE(s(:,1,isym)) + xk_irr(2)*DBLE(s(:,2,isym)) + xk_irr(3)*DBLE(s(:,3,isym))
        !
        IF (equiv_k(xk_ful,xk_tmp)) irr2ful_sym(ik) = isym
        IF (irr2ful_sym(ik) .GT. 1) EXIT
        !
     ENDDO
     !
     IF (irr2ful_sym(ik) .EQ. 0) CALL errore ('bz_index','Cannot find the rotation matrix',ik)
     !
404  CONTINUE
  ENDDO
  !
#ifdef __PARA
  CALL mp_barrier (inter_pool_comm)
  CALL mp_sum (irr2ful_sym,inter_pool_comm)
#ENDIF
  !
  !
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     IF (flag .EQ. 'k') OPEN (9999,FILE='BTE/META/irr2ful_sym',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='replace')
     IF (flag .EQ. 'q') OPEN (9999,FILE='BTE/META/irr2ful_sym_q',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='replace')
     DO ik = 1, nk_ful
        WRITE (9999,REC=ik) irr2ful_sym(ik)
     ENDDO
     CLOSE (9999)
     !
  ENDIF
  !
  !
  CLOSE (1111)
  CLOSE (2222)
  !
  DEALLOCATE (irr2ful_sym)
  !
#ifdef __PARA
  CALL mp_barrier (inter_pool_comm)
#ENDIF
  !
  !
  !
  !------------------------------------------------
  ! 7. symmat_kpt and ignore_kpt
  !------------------------------------------------
  !
  ! symmat_kpt is in crystal coordinate
  !
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     IF (flag .EQ. 'k') OPEN (1111,FILE='BTE/META/xkf_ful_cryst',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='old')
     IF (flag .EQ. 'q') OPEN (1111,FILE='BTE/META/xqf_ful_cryst',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='old')
     IF (flag .EQ. 'k') OPEN (8888,FILE='BTE/META/symmat_kpt',FORM='unformatted',ACCESS='direct',RECL=9*DP,STATUS='replace')
     IF (flag .EQ. 'q') OPEN (8888,FILE='BTE/META/symmat_qpt',FORM='unformatted',ACCESS='direct',RECL=9*DP,STATUS='replace')
     IF (flag .EQ. 'k') OPEN (9999,FILE='BTE/META/ignore_kpt',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='replace')
     IF (flag .EQ. 'q') OPEN (9999,FILE='BTE/META/ignore_qpt',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='replace')
     !
     DO ik = 1, nk_ful
        !
        READ (1111,REC=ik) xk_ful(1:3)
        cnt = 0 
        symmat = 0.0d0
        !
        DO isym = 1, n_rot
           !
           xk_tmp(:) = xk_ful(1)*DBLE(s(:,1,isym)) + xk_ful(2)*DBLE(s(:,2,isym)) + xk_ful(3)*DBLE(s(:,3,isym))
           !
           IF (equiv_k(xk_ful,xk_tmp)) THEN
              !
              symmat = symmat + DBLE(s(:,:,isym))
              cnt = cnt + 1
              !
           ENDIF
           !
        ENDDO
        !
        IF (cnt .GE. 1) THEN
           WRITE (8888,REC=ik) symmat(1,1:3)/DBLE(cnt), symmat(2,1:3)/DBLE(cnt), symmat(3,1:3)/DBLE(cnt)
        ELSE
           CALL errore ('bz_index','Cannot find symmat_kpt matrix',ik)
        ENDIF
        !
        IF ( .NOT. ( ABS(symmat(1,1)-1.0d0) .LT. eps .AND. &
                     ABS(symmat(1,2)-0.0d0) .LT. eps .AND. &
                     ABS(symmat(1,3)-0.0d0) .LT. eps .AND. &
                     ABS(symmat(2,1)-0.0d0) .LT. eps .AND. &
                     ABS(symmat(2,2)-1.0d0) .LT. eps .AND. &
                     ABS(symmat(2,3)-0.0d0) .LT. eps .AND. &
                     ABS(symmat(3,1)-0.0d0) .LT. eps .AND. &
                     ABS(symmat(3,2)-0.0d0) .LT. eps .AND. &
                     ABS(symmat(3,3)-1.0d0) .LT. eps ) ) THEN
           WRITE (9999,REC=ik) 1
        ELSE
           WRITE (9999,REC=ik) 0
        ENDIF
        !
     ENDDO
     !
     CLOSE (1111)
     CLOSE (8888)
     CLOSE (9999)
     !
  ENDIF
  !
#ifdef __PARA
  CALL mp_barrier (inter_pool_comm)
#ENDIF
  !
  !
  !------------------------------------------------
  ! 8. xkf_fbz
  !------------------------------------------------
  !
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     IF (flag .EQ. 'k') OPEN (2222,FILE='BTE/META/xkf_ful_cryst',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='old')
     IF (flag .EQ. 'q') OPEN (2222,FILE='BTE/META/xqf_ful_cryst',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='old')
     IF (flag .EQ. 'k') OPEN (9999,FILE='BTE/META/xkf_fbz_cart',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='replace')
     IF (flag .EQ. 'q') OPEN (9999,FILE='BTE/META/xqf_fbz_cart',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='replace')
     !
     DO ik = 1, nk_ful
        !
        READ (2222,REC=ik) xk(1:3)
        !
        xk_fbz = 0.0d0
        leng = +9999999999.9d0
        !
        DO ig1 = -1, 1
           DO ig2 = -1, 1
              DO ig3 = -1, 1
                 !
                 xk_g(1) = xk(1) + DBLE(ig1)
                 xk_g(2) = xk(2) + DBLE(ig2)
                 xk_g(3) = xk(3) + DBLE(ig3)
                 !
                 CALL cryst_to_cart (1,xk_g,bg,1)
                 !  
                 leng_new = DOT_PRODUCT(xk_g,xk_g)
                 !
                 IF (leng_new .LT. leng) THEN
                    !
                    leng = leng_new
                    xk_fbz(:) = xk_g(:) ! cart
                    !
                 ENDIF
                 !    
              ENDDO
           ENDDO
        ENDDO
        !
        WRITE (9999,REC=ik) xk_fbz(1:3)    
        !
     ENDDO
     !
     CLOSE (2222)
     CLOSE (9999)
     !
  ENDIF
  !
#ifdef __PARA
  CALL mp_barrier (inter_pool_comm)
#ENDIF
  !
  !
  !------------------------------------------------
  ! 9. check and done
  !------------------------------------------------
  !
  IF (flag .EQ. 'k') THEN
     !
     OPEN (1111,FILE='BTE/META/xkf_irr_cryst',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='old')
     OPEN (1155,FILE='BTE/META/wkf_irr',FORM='unformatted',ACCESS='direct',RECL=DP,STATUS='old')
     OPEN (2222,FILE='BTE/EPCHECK/xkf_irr_cryst.dat',STATUS='replace')
     OPEN (3333,FILE='BTE/EPCHECK/xkf_irr_cart.dat',STATUS='replace')
     !
     DO ik = 1, nk_irr
        !
        READ (1111,REC=ik) xk(1:3)
        READ (1155,REC=ik) wk
        !
        WRITE (2222,'(4f17.12)') xk(1:3), wk
        CALL cryst_to_cart (1,xk(:),bg,1)
        WRITE (3333,'(4f17.12)') xk(1:3), wk
        !
     ENDDO
     !
     CLOSE (1111)
     CLOSE (1155)
     CLOSE (2222)
     CLOSE (3333)
     !
  ENDIF
  !
  CALL CPU_TIME (t1)
  WRITE (stdout,'(13x,a,f8.2,a/)') 'Time : ', (t1-t0), ' s'
  !
  DEALLOCATE (symmat_lat) 
  DEALLOCATE (ful2irr)
  !
#ifdef __PARA
  CALL mp_barrier (inter_pool_comm)
#ENDIF
  !
  !
END SUBROUTINE bz_index



!---------------------------------------------------------------------------------
SUBROUTINE reduce_index (etf)
!---------------------------------------------------------------------------------
#INCLUDE "f_defs.h"
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : at, bg
  USE epwcom,    ONLY : nptype, nbndsub
  USE elph2,     ONLY : nbnd_red, cbnd_emin, vbnd_emax, ibndmin, ibndmax, &
                        ef_m, delta_egap, outside_gap, cfsthick, vfsthick
  USE bte_var
#ifdef __PARA
  USE mp,        ONLY : mp_barrier, mp_sum, mp_min, mp_bcast
  USE io_global, ONLY : ionode_id, stdout
  USE mp_global, ONLY : my_pool_id, inter_pool_comm
#ENDIF
  !
  IMPLICIT NONE
  !
  REAL(KIND=DP), INTENT(IN) :: etf(nbndsub,nk_irr)
  REAL(KIND=DP)             :: xk(3), t0, t1
  INTEGER                   :: ik, ik0, ik_red, ibnd
  LOGICAL                   :: within_range
  INTEGER, ALLOCATABLE      :: rirr2irr_tmp(:), rful2ful_tmp(:)
  !
  !
  CALL CPU_TIME (t0)
  !
  ALLOCATE (rirr2irr_tmp(nk_irr))
  rirr2irr_tmp = 0
  !
  nk_irr_red = 0
  DO ik = 1, nk_irr
     !
     within_range = .FALSE.
     !
     DO ibnd = ibndmin, ibndmax
        !
        IF (nptype .EQ. 'n') THEN
           IF (etf(ibnd,ik) .GE. cbnd_emin .AND. etf(ibnd,ik) .LE. cbnd_emin+cfsthick) within_range = .TRUE.
        ELSEIF (nptype .EQ. 'p') THEN 
           IF (etf(ibnd,ik) .GE. vbnd_emax-vfsthick .AND. etf(ibnd,ik) .LE. vbnd_emax) within_range = .TRUE.
        ELSE
           IF (etf(ibnd,ik) .GE. vbnd_emax-vfsthick .AND. etf(ibnd,ik) .LE. cbnd_emin+cfsthick) within_range = .TRUE.
        ENDIF
        !
     ENDDO
     !
     IF (within_range) THEN
        !
        nk_irr_red = nk_irr_red + 1
        rirr2irr_tmp(nk_irr_red) = ik
        !
     ENDIF
     !
  ENDDO
  !
  ALLOCATE (rirr2irr(nk_irr_red))
  rirr2irr(1:nk_irr_red) = rirr2irr_tmp(1:nk_irr_red)
  !
  !
  ! mapping to ful-BZ
  ALLOCATE(rful2ful_tmp(nk_ful))
  rful2ful_tmp = 0
  !
  nk_ful_red = 0
  DO ik_red = 1, nk_irr_red
     !
     DO ik = 1, nk_ful
        !
        IF (ful2irr(ik) .EQ. rirr2irr(ik_red)) THEN
           !
           nk_ful_red = nk_ful_red + 1
           rful2ful_tmp(nk_ful_red) = ik
           !
        ENDIF
        !
     ENDDO
     !
  ENDDO
  !
  ALLOCATE (rful2ful(nk_ful_red))
  ALLOCATE (ful2rful(nk_ful))
  rful2ful = 0
  ful2rful = 0
  !
  rful2ful(1:nk_ful_red) = rful2ful_tmp(1:nk_ful_red)
  DO ik = 1, nk_ful_red
     ful2rful(rful2ful(ik)) = ik
  ENDDO
  !
  !
  ! red-ful to red-irr
  ALLOCATE(rful2rirr(nk_ful))
  rful2rirr = 0
  DO ik = 1, nk_ful_red
     DO ik0 = 1, nk_irr_red
        IF (ful2irr(rful2ful(ik)) .EQ. rirr2irr(ik0)) rful2rirr(ik) = ik0
     ENDDO
  ENDDO
  !
  !
  ! output file
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     ! binary
     OPEN (99999,FILE='BTE/META/nk_irr_red',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='replace')
        WRITE (99999,REC=1) nk_irr_red
     CLOSE (99999)
     !
     OPEN (99999,FILE='BTE/META/nk_ful_red',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='replace')
        WRITE (99999,REC=1) nk_ful_red
     CLOSE (99999)
     !
     OPEN (99999,FILE='BTE/META/rirr2irr',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='replace')
     DO ik = 1, nk_irr_red
        WRITE (99999,REC=ik) rirr2irr(ik)
     ENDDO
     CLOSE (99999)
     !
     OPEN (99999,FILE='BTE/META/rful2ful',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='replace')
     DO ik = 1, nk_ful_red
        WRITE (99999,REC=ik) rful2ful(ik)
     ENDDO
     CLOSE (99999)
     !
     OPEN (99999,FILE='BTE/META/ful2rful',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='replace')
     DO ik = 1, nk_ful
        WRITE (99999,REC=ik) ful2rful(ik)
     ENDDO
     CLOSE (99999)  
     !
     OPEN (99999,FILE='BTE/META/rful2rirr',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='replace')
     DO ik = 1, nk_ful_red
        WRITE (99999,REC=ik) rful2rirr(ik)
     ENDDO
     CLOSE (99999)
     !
  ENDIF
  !
  CALL CPU_TIME (t1)
  !
  !
  IF (nptype .EQ. 'n' .OR. nptype .EQ. 'p') THEN
     WRITE (stdout,'(/5x,a,f8.2,a)') 'Time for n/p-type selective scheme:', t1-t0, ' s'
  ELSE
     WRITE (stdout,'(/5x,a,f8.2,a)') 'Time for regular selective scheme:', t1-t0, ' s'
  ENDIF
  WRITE (stdout,'(/5x,a,i12,a,i10)') 'k points in ful-BZ: ', nk_ful, ' -> ', nk_ful_red
  WRITE (stdout,'(5x,a,i12,a,i10/)')  '            irr-BZ: ', nk_irr, ' -> ', nk_irr_red
  !
  !
  DEALLOCATE (rirr2irr_tmp)
  DEALLOCATE (rful2ful_tmp)
  DEALLOCATE (rirr2irr)
  DEALLOCATE (rful2ful)
  DEALLOCATE (ful2rful)
  DEALLOCATE (rful2rirr)
  !
  !
END SUBROUTINE reduce_index



!---------------------------------------------------------------------------------
SUBROUTINE reduce_index_ph (wf)
!---------------------------------------------------------------------------------
#INCLUDE "f_defs.h"
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : at, bg
  USE epwcom,    ONLY : phwmax, phkmax
  USE ions_base,     ONLY : nat
  USE phcom,         ONLY : nmodes
  USE constants_epw, ONLY : ryd2thz
  USE bte_var
#ifdef __PARA
  USE mp,        ONLY : mp_barrier, mp_sum, mp_min, mp_bcast
  USE io_global, ONLY : ionode_id, stdout
  USE mp_global, ONLY : my_pool_id, inter_pool_comm
#ENDIF
  !
  IMPLICIT NONE
  !
  REAL(KIND=DP), INTENT(IN) :: wf(nmodes,nq_irr)
  REAL(KIND=DP)             :: xq(3), qleng, t0, t1, qleng_max=-9.9d30, wf_max=-9.9d30
  INTEGER                   :: iq, iq0, iq_red, imode
  LOGICAL                   :: within_range
  INTEGER, ALLOCATABLE      :: rirr2irr_tmp(:), rful2ful_tmp(:)
  !
  !
  CALL CPU_TIME (t0)
  !
  ALLOCATE (rirr2irr_tmp(nq_irr))
  rirr2irr_tmp = 0
  !
  nq_irr_red = 0
  DO iq = 1, nq_irr
     !
     within_range = .FALSE.
     !
     xq(:) = xqf_irr(:,iq)
     CALL cryst_to_cart (1,xq,bg,1)
     qleng = SQRT(DOT_PRODUCT(xq,xq))
     !
     IF (qleng .GT. qleng_max) qleng_max = qleng
     !
     DO imode = 1, nmodes
        IF (wf(imode,iq) .GT. wf_max) wf_max = wf(imode,iq)
        IF (wf(imode,iq) .LE. phwmax .AND. qleng .LE. phkmax) within_range = .TRUE.
     ENDDO
     !
     IF (within_range) THEN
        !
        nq_irr_red = nq_irr_red + 1
        rirr2irr_tmp(nq_irr_red) = iq
        !
     ENDIF
     !
  ENDDO
  !
  ALLOCATE (rirr2irr_q(nq_irr_red))
  rirr2irr_q(1:nq_irr_red) = rirr2irr_tmp(1:nq_irr_red)
  !
  !
  ! mapping to ful-BZ
  ALLOCATE(rful2ful_tmp(nq_ful))
  rful2ful_tmp = 0
  !
  nq_ful_red = 0
  DO iq_red = 1, nq_irr_red
     !
     DO iq = 1, nq_ful
        !
        IF (ful2irr_q(iq) .EQ. rirr2irr_q(iq_red)) THEN
           !
           nq_ful_red = nq_ful_red + 1
           rful2ful_tmp(nq_ful_red) = iq
           !
        ENDIF
        !
     ENDDO
     !
  ENDDO
  !
  ALLOCATE (rful2ful_q(nq_ful_red))
  ALLOCATE (ful2rful_q(nq_ful))
  rful2ful_q = 0
  ful2rful_q = 0
  !
  rful2ful_q(1:nq_ful_red) = rful2ful_tmp(1:nq_ful_red)
  DO iq = 1, nq_ful_red
     ful2rful_q(rful2ful_q(iq)) = iq
  ENDDO
  !
  ! 
  ! red-ful to red-irr
  ALLOCATE(rful2rirr_q(nq_ful))
  rful2rirr_q = 0
  DO iq = 1, nq_ful_red
     DO iq0 = 1, nq_irr_red
        IF (ful2irr_q(rful2ful_q(iq)) .EQ. rirr2irr_q(iq0)) rful2rirr_q(iq) = iq0
     ENDDO
  ENDDO
  !
  !
  ! output file
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     ! binary
     OPEN (99999,FILE='BTE/META/nq_irr_red',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='replace')
        WRITE (99999,REC=1) nq_irr_red
     CLOSE (99999)
     !
     OPEN (99999,FILE='BTE/META/nq_ful_red',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='replace')
        WRITE (99999,REC=1) nq_ful_red
     CLOSE (99999)
     !
     OPEN (99999,FILE='BTE/META/rirr2irr_q',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='replace')
     DO iq = 1, nq_irr_red
        WRITE (99999,REC=iq) rirr2irr_q(iq)
     ENDDO
     CLOSE (99999)
     !
     OPEN (99999,FILE='BTE/META/rful2ful_q',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='replace')
     DO iq = 1, nq_ful_red
        WRITE (99999,REC=iq) rful2ful_q(iq)
     ENDDO
     CLOSE (99999)
     !
     OPEN (99999,FILE='BTE/META/ful2rful_q',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='replace')
     DO iq = 1, nq_ful
        WRITE (99999,REC=iq) ful2rful_q(iq)
     ENDDO
     CLOSE (99999)  
     !
     OPEN (99999,FILE='BTE/META/rful2rirr_q',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='replace')
     DO iq = 1, nq_ful_red
        WRITE (99999,REC=iq) rful2rirr_q(iq)
     ENDDO
     CLOSE (99999)
     !
  ENDIF
  !
  CALL CPU_TIME (t1)
  !
  !
  WRITE (stdout,'(/5x,a,f8.2,a)') 'Time for regular selective scheme:', t1-t0, ' s'
  WRITE (stdout,'(/5x,a,es13.4,a,es13.4,a)') 'Max. phonon frequency  : ', wf_max*ryd2thz,    ' [THz]   ; Cutoff frequency  :', phwmax*ryd2thz, ' [THz]'
  WRITE (stdout,'(5x,a,es13.4,a,es13.4,a)') 'Max. phonon wavevector : ', qleng_max, ' [2pi/a] ; Cutoff wavevector :', phkmax, ' [2pi/a]'
  WRITE (stdout,'(/5x,a/)') 'q point will be selected if frequency < phwmax and wavevector < phkmax'  
  WRITE (stdout,'(5x,a,i12,a,i10)')  'q points in ful-BZ: ', nq_ful, ' -> ', nq_ful_red
  WRITE (stdout,'(5x,a,i12,a,i10/)') '            irr-BZ: ', nq_irr, ' -> ', nq_irr_red
  !
  !
  DEALLOCATE (rirr2irr_tmp)
  DEALLOCATE (rful2ful_tmp)
  ! not to deallocate in this case
  !DEALLOCATE (rirr2irr_q)
  !DEALLOCATE (rful2ful_q)
  !DEALLOCATE (ful2rful_q)
  !DEALLOCATE (rful2rirr_q)
  !
  !
END SUBROUTINE reduce_index_ph

!----------------------------------------------------------------------------
SUBROUTINE cpu_index ()
!----------------------------------------------------------------------------
!
! ik: ID of kpoint in order of serial
! seq2nscat(ik): ID of kpoint in order of nscat
!
#INCLUDE "f_defs.h"
  USE kinds,     ONLY : DP
  USE bte_var
#ifdef __PARA
  USE io_global, ONLY : ionode_id, stdout
  USE mp_global, ONLY : my_pool_id
  USE mp_world,  ONLY : nproc
#ENDIF
  !
  IMPLICIT NONE
  ! 
  INTEGER              :: ik, ik_irr, nk, imx, mxloc(1), ipool, cnt_1
  INTEGER              :: quot, rem
  INTEGER, ALLOCATABLE :: nk_cpu(:), cpu_stop(:)
  INTEGER, ALLOCATABLE :: nscat_ful(:), seq2nscat_tmp(:), nscat_tmp(:)
  REAL(KIND=DP)        :: t0, t1
  ! map
  INTEGER              :: ik_n0, ik_0, nk_n0
  INTEGER, ALLOCATABLE :: nscat_n0(:), n02ful(:)
  !
  !
  CALL CPU_TIME (t0)
  !
  WRITE (stdout,'(/5x,a)') 'Prepare the new sequence between CPU and kpoint'
  !
  ik_n0 = 0
  ik_0  = 0
  nk_n0 = 0
  !
  ! tmp. used
  ALLOCATE (nscat_ful(nk_ful_red))
  ALLOCATE (seq2nscat_tmp(nk_ful_red))
  ALLOCATE (nscat_tmp(nk_ful_red))
  nscat_ful     = 0
  seq2nscat_tmp = 0
  nscat_tmp     = 0
  !
  ALLOCATE (seq2nscat(nk_ful_red))
  ALLOCATE (nscat_new(nk_ful_red))
  ALLOCATE (nk_cpu(nproc))
  ALLOCATE (cpu_stop(nproc))
  seq2nscat = 0
  nscat_new = 0
  nk_cpu    = 0
  cpu_stop  = 0
  !
  ! nk_ful: total number of kpoints
  ! nproc: total number of CPUs
  ! quot: number of kpoints on a CPU
  ! rem: number of leading CPUs which have one more kpoint 
  rem = MOD(nk_ful_red,nproc)
  quot = (nk_ful_red-rem)/nproc
  !
  !
  DO ik = 1, nk_ful_red
     ik_irr = rful2rirr(ik)
     nscat_ful(ik) = nscat_all(ik_irr)
  ENDDO
  !
  ! retain the non-zero elements
  DO ik = 1, nk_ful_red
     IF (nscat_ful(ik) .GT. 0) nk_n0 =  nk_n0 + 1
  ENDDO 
  !
  ALLOCATE (nscat_n0(nk_n0))
  ALLOCATE (n02ful(nk_n0)) ! mapping
  nscat_n0 = 0
  n02ful = 0
  !
  DO ik = 1, nk_ful_red
     IF (nscat_ful(ik) .GT. 0) THEN
        ik_n0 = ik_n0 + 1
        nscat_n0(ik_n0) = nscat_ful(ik)
        n02ful(ik_n0) = ik
     ENDIF
  ENDDO 
  !
  ! ordering the sequence
  DO ik = 1, nk_n0
     !
     mxloc = MAXLOC(nscat_n0)        
     imx = mxloc(1)
     !
     seq2nscat_tmp(ik) = n02ful(imx)
     nscat_tmp(ik) = nscat_n0(imx)
     nscat_n0(imx) = -100
     !
  ENDDO
  !
  ! fill the zero term
  DO ik = 1, nk_ful_red
     IF (nscat_ful(ik) .EQ. 0) THEN
        ik_0 = ik_0 + 1
        seq2nscat_tmp(nk_n0+ik_0) = ik
     ENDIF
  ENDDO 
  !
  ! calculate the number of kpoints on each CPU
  DO ipool = 1, nproc
     !
     IF (ipool .LE. rem) THEN
        nk_cpu(ipool) = quot + 1
     ELSE
        nk_cpu(ipool) = quot
     ENDIF
     !
  ENDDO
  !
  ! 
  ik = 0
  DO ipool = 1, nproc
     !
     DO nk = 1, nk_cpu(ipool)
        !
        ik = ik + 1
        seq2nscat(ik) = seq2nscat_tmp(ipool+(nk-1)*nproc)
        nscat_new(ik) = nscat_tmp(ipool+(nk-1)*nproc)
        !
     ENDDO
     !
  ENDDO
  !
  CALL CPU_TIME (t1)
  WRITE (stdout,'(/13x,a,f8.2,a)') 'Time : ', (t1-t0), ' s'
  !
  !OPEN (99999,FILE='BTE/EPCHECK/seq2nscat.dat')
  !WRITE (99999,'(a)') 'ik_ful_red(seq) | ik_ful_red(nscat) | ik_ful(seq) | ik_ful(nscat) | nscat'
  !DO ik = 1, nk_ful_red
  !   WRITE (99999,'(5i10)') ik, seq2nscat(ik), rful2ful(ik), rful2ful(seq2nscat(ik)), nscat_new(ik)
  !ENDDO 
  !
  ! deallocate
  DEALLOCATE (nscat_ful)
  DEALLOCATE (seq2nscat_tmp)
  DEALLOCATE (nscat_tmp)
  DEALLOCATE (nk_cpu)
  DEALLOCATE (cpu_stop)
  DEALLOCATE (nscat_n0)
  DEALLOCATE (n02ful)
  !
END SUBROUTINE cpu_index



!---------------------------------------------------------------------------------
SUBROUTINE rotate_vel
!---------------------------------------------------------------------------------
#INCLUDE "f_defs.h"
  USE kinds,         ONLY : DP
  USE cell_base,     ONLY : at, bg
  USE epwcom,        ONLY : nbndsub, bte, nkf1, nkf2, nkf3
  USE elph2,         ONLY : ibndmin, ibndmax, nbnd_red, vel_ful, etf_all, vel_all
  USE bte_var
  USE constants_epw, ONLY : ryd2ev, rydcm1, au2m, au2s
#ifdef __PARA
  USE mp,            ONLY : mp_barrier, mp_sum
  USE io_global,     ONLY : ionode_id, stdout
  USE mp_global,     ONLY : my_pool_id, inter_pool_comm
#endif
  !
  IMPLICIT NONE
  !
  INTEGER :: ik_ful_red, ik_irr_red, ibnd, ibnd0, ir, ik_ful, ik_irr, isym
  REAL(KIND=DP) :: t0, t1
  ! para
  INTEGER :: nk_pol, ik_star, ik_stop
  !
  !
  ! Form here, vel_ful belongs to nk_ful_red
  !
  CALL CPU_TIME (t0)
  !
  WRITE (stdout,'(/5x,a)') 'Rotate electron velocity from red-irr-BZ to red-ful-BZ according to crystal symmetry'
  !
  ALLOCATE (vel_ful(3,nbnd_red,nk_ful_red))
  vel_ful = 0.0d0
  !
  CALL para_bounds (ik_star,ik_stop,nk_ful_red)
  !
  ! note: vel_all, etf_all are in dimension of nbndsub, while vel_ful is of nbnd_red
  OPEN (2222,FILE='BTE/META/irr2ful_sym',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='old')  
  DO ik_ful_red = ik_star, ik_stop
     !
     ik_ful = rful2ful(ik_ful_red)
     ik_irr_red = rful2rirr(ik_ful_red)
     !
     READ (2222,REC=ik_ful) isym
     !
     DO ibnd = 1, nbnd_red
        !
        ibnd0 = ibnd+ibndmin-1
        CALL cryst_to_cart (1,vel_all(:,ibnd0,ik_irr_red),at,-1)
        !
        DO ir = 1, 3
           ! now the vel_ful is in crystal coordinate
           vel_ful(ir,ibnd,ik_ful_red) = vel_all(1,ibnd0,ik_irr_red)*symmat_lat(ir,1,isym) + &
                                         vel_all(2,ibnd0,ik_irr_red)*symmat_lat(ir,2,isym) + &
                                         vel_all(3,ibnd0,ik_irr_red)*symmat_lat(ir,3,isym)
           !
        ENDDO
        !
        CALL cryst_to_cart (1,vel_all(:,ibnd0,ik_irr_red),bg,1)
        CALL cryst_to_cart (1,vel_ful(:,ibnd,ik_ful_red),bg,1)
        !
        ! apply symmetrizer to velocity 
        CALL symmetrizer (vel_ful(:,ibnd,ik_ful_red),ik_ful,'k')
        !
     ENDDO
     !
  ENDDO
  CLOSE (2222)
  !  
#ifdef __PARA
  CALL mp_barrier (inter_pool_comm)
  CALL mp_sum (vel_ful,inter_pool_comm)
#ENDIF
  !
  ! output and check
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     OPEN (99999,FILE='BTE/EPCHECK/electron_ful_red.dat')
     DO ik_ful_red = 1, nk_ful_red
        !
        ik_ful = rful2ful(ik_ful_red)
        ik_irr_red = rful2rirr(ik_ful_red)
        !
        DO ibnd = 1, nbnd_red
           ibnd0 = ibnd+ibndmin-1
           WRITE (99999,'(i12,i5,f15.4,3f17.4)') ik_ful, ibnd0, etf_all(ibnd0,ik_irr_red)*ryd2ev, vel_ful(1:3,ibnd,ik_ful_red)*(au2m/au2s)
        ENDDO
     ENDDO
     CLOSE (99999)
     !
  ENDIF
  !
  CALL CPU_TIME (t1)
  WRITE (stdout,'(/13x,a,f8.2,a)') 'Time : ', (t1-t0), ' s'
  !
END SUBROUTINE rotate_vel



!---------------------------------------------------------------------------------
SUBROUTINE rotate_vph
!---------------------------------------------------------------------------------
#INCLUDE "f_defs.h"
  USE kinds,         ONLY : DP
  USE cell_base,     ONLY : at, bg
  USE phcom,         ONLY : nmodes
  USE epwcom,        ONLY : nbndsub, bte, nqf1, nqf2, nqf3
  USE elph2,         ONLY : vph_ful, wf_all, vph_all
  USE bte_var
  USE constants_epw, ONLY : ryd2ev, rydcm1, au2m, au2s, au2ps
#ifdef __PARA
  USE mp,            ONLY : mp_barrier, mp_sum
  USE io_global,     ONLY : ionode_id, stdout
  USE mp_global,     ONLY : my_pool_id, inter_pool_comm
#endif
  !
  IMPLICIT NONE
  !
  INTEGER :: iq_ful_red, iq_irr_red, imode, ir, iq_ful, iq_irr, isym_q
  REAL(KIND=DP) :: t0, t1
  ! para
  INTEGER :: nq_pol, iq_star, iq_stop
  !
  !
  ! From here, vph_ful belongs to nq_ful_red
  !
  CALL CPU_TIME (t0)
  !
  WRITE (stdout,'(/5x,a)') 'Rotate phonon velocity from red-irr-BZ to red-ful-BZ according to crystal symmetry'
  !
  ALLOCATE (vph_ful(3,nmodes,nq_ful_red))
  vph_ful = 0.0d0
  !
  CALL para_bounds (iq_star,iq_stop,nq_ful_red)
  !
  !
  OPEN (2222,FILE='BTE/META/irr2ful_sym_q',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='old')  
  DO iq_ful_red = iq_star, iq_stop
     !
     iq_ful = rful2ful_q(iq_ful_red)
     iq_irr_red = rful2rirr_q(iq_ful_red)
     !
     READ (2222,REC=iq_ful) isym_q
     !
     DO imode = 1, nmodes
        !
        CALL cryst_to_cart (1,vph_all(:,imode,iq_irr_red),at,-1)
        !
        DO ir = 1, 3
           ! now the vph_ful is in crystal coordinate
           vph_ful(ir,imode,iq_ful_red) = vph_all(1,imode,iq_irr_red)*symmat_lat(ir,1,isym_q) + &
                                          vph_all(2,imode,iq_irr_red)*symmat_lat(ir,2,isym_q) + &
                                          vph_all(3,imode,iq_irr_red)*symmat_lat(ir,3,isym_q)
           !
        ENDDO
        !
        CALL cryst_to_cart (1,vph_all(:,imode,iq_irr_red),bg,1)
        CALL cryst_to_cart (1,vph_ful(:,imode,iq_ful_red),bg,1)
        !
        ! apply symmetrizer to velocity 
        CALL symmetrizer (vph_ful(:,imode,iq_ful_red),iq_ful,'q')
        !
     ENDDO
     !
  ENDDO
  CLOSE (2222)
  !  
#ifdef __PARA
  CALL mp_barrier (inter_pool_comm)
  CALL mp_sum (vph_ful,inter_pool_comm)
#ENDIF
  !
  ! output and check
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     OPEN (99999,FILE='BTE/EPCHECK/phonon_ful_red.dat')
     DO iq_ful_red = 1, nq_ful_red
        !
        iq_ful = rful2ful_q(iq_ful_red)
        iq_irr_red = rful2rirr_q(iq_ful_red)
        !
        DO imode = 1, nmodes
           WRITE (99999,'(i12,i5,f15.4,3f17.4)') iq_ful, imode, wf_all(imode,iq_irr_red)*rydcm1, vph_ful(1:3,imode,iq_ful_red)*(au2m/au2s)
        ENDDO
     ENDDO
     CLOSE (99999)
     !
  ENDIF
  !
  !
  CALL CPU_TIME (t1)
  WRITE (stdout,'(/13x,a,f8.2,a)') 'Time : ', (t1-t0), ' s'
  !
END SUBROUTINE rotate_vph




!---------------------------------------------------------------------------------
SUBROUTINE symmetrizer (v, pt_num, flag)
!---------------------------------------------------------------------------------
#INCLUDE "f_defs.h"
  USE kinds,        ONLY : DP
  USE cell_base,    ONLY : at, bg
  USE phcom,        ONLY : nmodes
  USE epwcom,       ONLY : nbndsub
  USE bte_var,        ONLY : ignore_kpt, ignore_qpt
#ifdef __PARA
  USE io_global,    ONLY : ionode_id
  USE mp_global,    ONLY : my_pool_id
#ENDIF
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)          :: pt_num
  CHARACTER(LEN=*), INTENT(IN) :: flag
  REAL(KIND=DP), INTENT(INOUT) :: v(3)
  ! v should be in cartesian coordinate
  REAL(KIND=DP)                :: v_sym(3), symmat_pt(3,3)
  INTEGER                      :: ir, is, ignore
  !
  !
  IF (flag .EQ. 'k') OPEN (9999,FILE='BTE/META/ignore_kpt',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='old')
  IF (flag .EQ. 'q') OPEN (9999,FILE='BTE/META/ignore_qpt',FORM='unformatted',ACCESS='direct',RECL=4,STATUS='old')
  READ (9999,REC=pt_num) ignore
  CLOSE (9999)
  !
  IF (ignore .NE. 0) THEN
     !
     IF (flag .EQ. 'k') OPEN (9999,FILE='BTE/META/symmat_kpt',FORM='unformatted',ACCESS='direct',RECL=9*DP,STATUS='old')
     IF (flag .EQ. 'q') OPEN (9999,FILE='BTE/META/symmat_qpt',FORM='unformatted',ACCESS='direct',RECL=9*DP,STATUS='old')
     READ (9999,REC=pt_num) symmat_pt(1,1:3), symmat_pt(2,1:3), symmat_pt(3,1:3)
     !
     CALL cryst_to_cart (1,v(:),at,-1)
     !
     DO ir = 1, 3
        ! operate in crystal coordinate
        v_sym(ir) = v(1)*symmat_pt(ir,1) + &
                    v(2)*symmat_pt(ir,2) + &
                    v(3)*symmat_pt(ir,3)
        !
     ENDDO
     !
     CALL cryst_to_cart (1,v_sym(:),bg,1)
     !  
     v = v_sym
     !
     CLOSE (9999)
     !
  ENDIF
  !
END SUBROUTINE symmetrizer

