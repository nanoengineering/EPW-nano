SUBROUTINE ep_check (filename_check)
  !
  USE kinds,         ONLY : DP
  USE pwcom,         ONLY : at, bg, nkstot
  USE ions_base,     ONLY : amass, ityp
  USE phcom,         ONLY : nq1, nq2, nq3, nmodes
  USE epwcom,        ONLY : nbndsub, lpolar, eig_read
  USE constants_epw, ONLY : ryd2ev, ryd2mev, rydcm1, one, two, twopi
  USE elph2,         ONLY : nrr_k, nrr_q, irvec, ndegen_k, ndegen_q, chw, chw_ks, epsi, zstar, &
                            ! THL
                            ef_m, delta_egap, vbnd_num
#ifdef __PARA
  USE mp,            ONLY : mp_barrier
  USE io_global,     ONLY : ionode_id, stdout
  USE mp_global,     ONLY : my_pool_id, inter_pool_comm
#ENDIF
  !
  IMPLICIT NONE
  !
  !
  CHARACTER(LEN=256), INTENT(IN) :: filename_check
  !
  REAL(KIND=DP), PARAMETER      :: eps = 1.0d-2/ryd2mev
  INTEGER                       :: nk, nq, ns
  INTEGER, ALLOCATABLE          :: kk(:), ib(:), jb(:)
  REAL(KIND=DP), ALLOCATABLE    :: xk(:,:), xq(:,:)
  REAL(KIND=DP), ALLOCATABLE    :: et(:,:), wf(:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: cufkk(:,:,:), uf(:,:,:), epmatf(:,:,:,:,:), eph_vogl(:,:)
  !
  INTEGER                       :: ik, iq, is, ibnd, jbnd, imode, nu, mu, na
  REAL(KIND=DP)                 :: xkk(3), xxq(3), xkq(3), x0(3), x1(3), et_ks(nbndsub), vq(3,nmodes)
  COMPLEX(KIND=DP)              :: cufkq(nbndsub, nbndsub), epmatwef(nbndsub,nbndsub,nrr_q,nmodes), bmatf(nbndsub,nbndsub), &
                                   uf_cart(nmodes,nmodes)
  REAL(KIND=DP)                 :: delta   
  CHARACTER(LEN=256)            :: coord, filename, k_id, ib_id, jb_id
  !
  !
  ! load file
  OPEN (1111,FILE=filename_check,STATUS='old')
  !
  READ (1111,*) nk, coord
  WRITE (stdout,'(/5x,a,i4,a)') 'Check electron energy on ', nk, ' k point(s)'
  ALLOCATE (xk(3,nk))
  DO ik = 1, nk
     READ (1111,*) xk(1:3,ik)
     IF (coord .EQ. 'cart' .OR. coord .EQ. 'carte' .OR. coord .EQ. 'cartesian' .OR. &
         coord .EQ. 'Cart' .OR. coord .EQ. 'Carte' .OR. coord .EQ. 'Cartesian' .OR. &
         coord .EQ. '1') CALL cryst_to_cart (1, xk(:,ik), at, -1)
  ENDDO 
  !
  READ (1111,*) nq, coord
  WRITE (stdout,'(/5x,a,i4,a)') 'Check phonon frequency on ', nq, ' q point(s)'
  ALLOCATE (xq(3,nq))
  DO iq = 1, nq
     READ (1111,*) xq(1:3,iq)
     IF (coord .EQ. 'cart' .OR. coord .EQ. 'carte' .OR. coord .EQ. 'cartesian' .OR. &
         coord .EQ. 'Cart' .OR. coord .EQ. 'Carte' .OR. coord .EQ. 'Cartesian' .OR. &
         coord .EQ. '1') CALL cryst_to_cart (1, xq(:,iq), at, -1)
  ENDDO 
  !
  READ (1111,*) ns
  WRITE (stdout,'(/5x,a,i3,a/)') 'Check e-ph coupling matrix on ', ns, ' k point(s)'
  ALLOCATE (kk(ns),ib(ns),jb(ns))
  DO is = 1, ns
     READ (1111,*) kk(is), ib(is), jb(is)
     WRITE (stdout,'(13x,a,3f10.6,a,i4,a,i4)') 'k = (', xk(1:3,kk(is)), '), band', ib(is), ' ->', jb(is)
  ENDDO 
  !
  CLOSE (1111)
  !
  ALLOCATE (et(nbndsub,nk))
  ALLOCATE (cufkk(nbndsub,nbndsub,nk))
  ALLOCATE (wf(nmodes,nq))
  ALLOCATE (uf(nmodes,nmodes,nq))
  ALLOCATE (epmatf(nbndsub,nbndsub,nmodes,nk,nq))
  ALLOCATE (eph_vogl(nmodes,nq))
  et       =  0.0d0
  cufkk    = (0.0d0,0.0d0)
  wf       =  0.0d0
  uf       = (0.0d0,0.0d0)
  epmatf   = (0.0d0,0.0d0)
  eph_vogl = (0.0d0,0.0d0)
  !
  !
  ! electron energy
  DO ik = 1, nk
     !
     xkk(:) = xk(:,ik) ! xk is assumed to be in crys coord
     !
     IF (eig_read) THEN
        CALL hamwan2bloch (nbndsub, nrr_k, irvec, ndegen_k, xkk, cufkk(:,:,ik), et_ks, chw_ks)   
     ENDIF
     !
     CALL hamwan2bloch (nbndsub, nrr_k, irvec, ndegen_k, xkk, cufkk(:,:,ik), et(:,ik), chw)      
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
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     OPEN (2222,FILE='BTE/EPCHECK/eband.dat',STATUS='replace')
     DO ik = 1, nk
        !
        IF (ik .EQ. 1) THEN
           delta = 0.0d0
        ELSE
           x0(:) = xk(:,ik-1)
           x1(:) = xk(:,ik)
           CALL cryst_to_cart (1, x0, bg, 1)
           CALL cryst_to_cart (1, x1, bg, 1)
           delta = SQRT((x1(1)-x0(1))**2.0d0+(x1(2)-x0(2))**2.0d0+(x1(3)-x0(3))**2.0d0) + delta
        ENDIF
        !
        WRITE (2222,'(4f12.6,48f10.4)') xk(1:3,ik), delta, et(1:nbndsub,ik)*ryd2ev
        !
     ENDDO
     CLOSE (2222)
     !
  ENDIF
  !
  !
  CALL mp_barrier (inter_pool_comm)
  !
  !
  ! phonon frequency
  DO iq = 1, nq
     !
     xxq = xq(:,iq)
     !
     CALL dynwan2bloch (nmodes, nrr_q, irvec, ndegen_q, xxq, uf(:,:,iq), wf(:,iq), vq)
     !
     DO nu = 1, nmodes
        !
        IF (wf(nu,iq) .GT. 0.0d0) THEN
           wf(nu,iq) =  SQRT(ABS(wf(nu,iq)))
        ELSE
           wf(nu,iq) = -SQRT(ABS(wf(nu,iq)))
        ENDIF
        !
        DO mu = 1, nmodes
           na = (mu-1)/3 + 1
           uf(mu,nu,iq) = uf(mu,nu,iq) / SQRT(amass(ityp(na)))
        ENDDO
        !
     ENDDO   
     !
  ENDDO
  !
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     OPEN (2222,FILE='BTE/EPCHECK/phband.dat',STATUS='replace')
     DO iq = 1, nq
        !
        IF (iq .EQ. 1) THEN
           delta = 0.0d0
        ELSE
           x0(:) = xq(:,iq-1)
           x1(:) = xq(:,iq)
           CALL cryst_to_cart (1, x0, bg, 1)
           CALL cryst_to_cart (1, x1, bg, 1)
           delta = SQRT((x1(1)-x0(1))**2.0d0+(x1(2)-x0(2))**2.0d0+(x1(3)-x0(3))**2.0d0) + delta
        ENDIF
        !
        WRITE (2222,'(4f12.6,48f10.4)') xq(1:3,iq), delta, wf(1:nmodes,iq)*rydcm1
       !
     ENDDO
     CLOSE (2222)
     !
  ENDIF
  !
  !
  CALL mp_barrier (inter_pool_comm)
  !
  !
  ! prepare the eph-vertex
  IF (lpolar) THEN
     !
     uf_cart  = 0.0d0
     DO imode = 1, nmodes
        uf_cart(imode,imode) = 1.0d0
     ENDDO
     !
     DO iq = 1, nq
        !
        xxq = xq(:,iq)
        CALL cryst_to_cart (1, xxq, bg, 1) 
        !
        IF (MAXVAL(ABS(xxq)) .NE. 0.0d0) THEN
           CALL polar_eph (uf_cart, xxq, eph_vogl(:,iq))
        ELSE
           eph_vogl(:,iq) = 0.0d0
        ENDIF
        !
     ENDDO
     !
  ENDIF
  !
  ! eph matrix
  DO is = 1, ns
     !
     ik = kk(is)
     xkk(:) = xk(:,ik) ! xk is assumed to be in crys coord
     !
     !
     ! ======================================================
     ! epmat : Wannier el and Wannier ph -> Bloch el and Wannier ph
     ! ======================================================
     !
     CALL ephwan2bloche (nmodes, xkk, irvec, ndegen_k, nrr_q, cufkk(:,:,ik), epmatwef, nbndsub, nrr_k)
     !
     !
     DO iq = 1, nq
        ! 
        xxq(:) = xq(:,iq) ! xq is assumed to be in crys coord
        xkq = xkk + xxq
        !
        !
        ! ======================================================
        ! hamiltonian : Wannier -> Bloch
        ! ======================================================
        !
        IF (eig_read) THEN    
           CALL hamwan2bloch (nbndsub, nrr_k, irvec, ndegen_k, xkq, cufkq, et_ks, chw_ks)
        ENDIF
        !   
        CALL hamwan2bloch (nbndsub, nrr_k, irvec, ndegen_k, xkq, cufkq, et_ks, chw) ! eigen value is not a matter here
        !
        !
        ! ======================================================
        ! epmat : Bloch el and Wannier ph -> Bloch el and Bloch ph
        ! ======================================================
        !
        CALL ephwan2bloch2 (nbndsub, nrr_q, irvec, ndegen_q, epmatwef, xxq, uf(:,:,iq), cufkk(:,:,ik), cufkq, epmatf(:,:,:,ik,iq), nmodes)
        ! 
        IF (lpolar) THEN
           !
           CALL compute_bmn_para2 (nbndsub, nkstot, cufkk(:,:,ik), cufkq, bmatf)
           !
           IF ((ABS(xxq(1)) .GT. eps) .OR. (ABS(xxq(2)).GT.eps) .OR. (ABS(xxq(3)).GT.eps)) THEN
              CALL cryst_to_cart (1, xxq, bg, 1)
              DO ibnd = 1, nbndsub
                 DO jbnd = 1, nbndsub
                    !CALL rgd_blk_epw2 (nq1, nq2, nq3, xxq, uf(:,:,iq), epmatf(ibnd,jbnd,:,ik,iq), &
                    !                   nmodes, epsi, zstar, bmatf(ibnd,jbnd), +1.d0)
                    CALL rgd_blk_epw3 (uf(:,:,iq), epmatf(ibnd,jbnd,:,ik,iq), eph_vogl(:,iq), &
                                       nmodes, bmatf(ibnd,jbnd), +1.d0)
                 ENDDO
              ENDDO
              CALL cryst_to_cart (1, xxq, at, -1)
           ENDIF
           !
        ENDIF
        !
     ENDDO ! q loop 
     !
  ENDDO ! k loop 
  !
  IF (my_pool_id .EQ. ionode_id) THEN
     !
     DO is = 1, ns
        !
        WRITE (k_id,*) kk(is)
        WRITE (ib_id,*) ib(is)
        WRITE (jb_id,*) jb(is)
        filename = 'BTE/EPCHECK/ephmat1_'//TRIM(ADJUSTL(k_id))//'_'//TRIM(ADJUSTL(ib_id))//'_'//TRIM(ADJUSTL(jb_id))//'.dat'
        OPEN (2222,FILE=filename,STATUS='replace')
        !
        DO iq = 1, nq
           !
           IF (iq .EQ. 1) THEN
              delta = 0.0d0
           ELSE
              CALL cryst_to_cart (1, xq(:,iq-1), bg, 1)
              CALL cryst_to_cart (1, xq(:,iq), bg, 1)
              delta = SQRT((xq(1,iq)-xq(1,iq-1))**2.0d0+(xq(2,iq)-xq(2,iq-1))**2.0d0+(xq(3,iq)-xq(3,iq-1))**2.0d0) + delta
              CALL cryst_to_cart (1, xq(:,iq-1), at, -1)
              CALL cryst_to_cart (1, xq(:,iq), at, -1)
           ENDIF
           !    
           WRITE (2222,'(4f12.6,48es14.4)') xq(1:3,iq), delta, &
                                            ( (ABS(epmatf(ib(is),jb(is),imode,kk(is),iq))**2.0d0)/2.0d0/wf(imode,iq) ,imode=1,nmodes) ! coupling strength ~ Ry^2
           !
        ENDDO 
        !
        CLOSE (2222)
        !
     ENDDO  
     !
  ENDIF  
  !
  !
  DEALLOCATE (xk)
  DEALLOCATE (xq)
  DEALLOCATE (kk,ib,jb)
  DEALLOCATE (et)
  DEALLOCATE (cufkk)
  DEALLOCATE (wf)
  DEALLOCATE (uf)
  DEALLOCATE (epmatf)
  DEALLOCATE (eph_vogl)
  !
END SUBROUTINE ep_check !
