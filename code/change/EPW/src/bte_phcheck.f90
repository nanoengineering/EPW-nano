SUBROUTINE phcheck_scat (filename_check)
  !
  ! output ph-ph el-ph scattering rate along L-Gamma-X, Qian
  ! modified from THL bte_check
  !
  USE kinds,         ONLY : DP
  USE pwcom,         ONLY : at, bg, nkstot
  USE ions_base,     ONLY : amass, ityp
  USE phcom,         ONLY : nq1, nq2, nq3, nmodes
  USE epwcom,        ONLY : nbndsub, lpolar, eig_read, neptemp, nepdope, phdrag, phkmax
  USE constants_epw, ONLY : ryd2ev, ryd2thz, ryd2mev, rydcm1, one, two, twopi, au2ps
  USE elph2,         ONLY : nrr_k, nrr_q, irvec, ndegen_k, ndegen_q, chw, chw_ks, epsi, zstar, &
                            ! Qian
                            gammai_mode_all, wf_irr, wf_all
  USE bte_var
  USE bte_func
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
  INTEGER, ALLOCATABLE          :: iq_gam(:), iq_ful(:)
  REAL(KIND=DP), ALLOCATABLE    :: xk(:,:), xq(:,:)
  REAL(KIND=DP), ALLOCATABLE    :: wf(:,:)
  !
  INTEGER                       :: ik, iq, iqq, is, ibnd, jbnd, imode, itemp, idope, nu, mu, na
  REAL(KIND=DP)                 :: xkk(3), xxq(3), xkq(3), xq_read(3), x0(3), x1(3), et_ks(nbndsub), vq(3,nmodes)
  COMPLEX(KIND=DP)              :: cufkq(nbndsub, nbndsub)
  REAL(KIND=DP)                 :: delta   
  CHARACTER(LEN=256)            :: coord, filename, k_id, ib_id, jb_id
  !
  CHARACTER(LEN=256)  :: phph_L_G_X, elph_L_G_X
  CHARACTER(LEN=12)   :: txnx
  CHARACTER(LEN=3)    :: itemp_num, idope_num
  !
  !
  ! load file
  OPEN (1111,FILE=filename_check,STATUS='old')
  !
  READ (1111,*) nk, coord
  ALLOCATE (xk(3,nk))
  DO ik = 1, nk
     READ (1111,*) xk(1:3,ik)
     IF (coord .EQ. 'cart' .OR. coord .EQ. 'carte' .OR. coord .EQ. 'cartesian' .OR. &
         coord .EQ. 'Cart' .OR. coord .EQ. 'Carte' .OR. coord .EQ. 'Cartesian' .OR. &
         coord .EQ. '1') CALL cryst_to_cart (1, xk(:,ik), at, -1)
  ENDDO 
  !
  READ (1111,*) nq, coord
  ALLOCATE (xq(3,nq))
  DO iq = 1, nq
     READ (1111,*) xq(1:3,iq)
     IF (coord .EQ. 'cart' .OR. coord .EQ. 'carte' .OR. coord .EQ. 'cartesian' .OR. &
         coord .EQ. 'Cart' .OR. coord .EQ. 'Carte' .OR. coord .EQ. 'Cartesian' .OR. &
         coord .EQ. '1') CALL cryst_to_cart (1, xq(:,iq), at, -1)
  ENDDO 
  !
  CLOSE (1111)
  !
  ALLOCATE (wf(nmodes,nq))
  ALLOCATE (iq_ful(nq))
  ALLOCATE (iq_gam(nq))
  wf       =  0.0d0
  iq_gam   =  0
  iq_ful   =  0
  !
  !
  ! electron scatter phonon rate along L-Gamma-X
  !
  ! phonon scatter phonon rate along L-Gamma-X
  !
  !
  DO iq = 1, nq
     !
     xxq = xq(:,iq)
     !
     OPEN (1234,FILE='BTE/META/xqf_ful_cryst',FORM='unformatted',ACCESS='direct',RECL=3*DP,STATUS='old')
     !
     ! match q coordinates with iq number by searching through xqf_ful_cryst
     !
     nqful: DO iqq = 1, nq_ful
       !
       READ (1234,REC=iqq) xq_read(1:3)
       !
       IF ((ABS(xq_read(1)-xxq(1)) .LT. 1.0d-4) .AND. (ABS(xq_read(2)-xxq(2)) .LT. 1.0d-4) .AND. (ABS(xq_read(3)-xxq(3)) .LT. 1.0d-4)) THEN
       iq_ful(iq) = iqq
       !
       iq_gam(iq) = ful2irr_q(iqq)
       !
       EXIT nqful
       !
       ENDIF
       !
     ENDDO nqful
     CLOSE (1234)
     !
     ! if phdrag, some q points along L-G-X might not be in the range
     wf(:,iq) = wf_irr(:,ful2irr_q(iq_ful(iq)))
     !
     !
  ENDDO
  !
  !
  IF (my_pool_id .EQ. ionode_id) THEN
  DO itemp = 1, neptemp
        DO idope = 1, nepdope
           !
           WRITE(itemp_num,'(i3)') itemp
           WRITE(idope_num,'(i3)') idope
           txnx = 'T'//TRIM(ADJUSTL(itemp_num))//'_N'//TRIM(ADJUSTL(idope_num))
           phph_L_G_X = 'BTE/'//TRIM(ADJUSTL(txnx))//'/phph_L_G_X.dat' 
           ! phph should be same across nepdope, put there for easy comparison with el-ph
           elph_L_G_X = 'BTE/'//TRIM(ADJUSTL(txnx))//'/elph_L_G_X.dat'
           !
        !
        OPEN (2222,FILE=phph_L_G_X,STATUS='replace')
        IF (.NOT. phdrag) OPEN (3333,FILE=elph_L_G_X,STATUS='replace')
        !
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
         ! ph_rate_ful only after calling read_shengbte in phdrag_shuffle.f90
        DO imode = 1, nmodes
         WRITE (2222,'(4f12.6,2f16.6)') xq(1:3,iq), delta, wf(imode,iq)*rydcm1, ph_rate_ful(itemp,imode,iq_ful(iq))/au2ps
         !
         IF (.NOT. phdrag) WRITE (3333,'(4f12.6,2f16.6)') xq(1:3,iq), delta, wf(imode,iq)*rydcm1, 2.0d0*gammai_mode_all(itemp,idope,imode,iq_gam(iq))*ryd2thz*twopi
!         ENDIF
        !
        ENDDO
        !
        ENDDO
        CLOSE (2222)
        IF (.NOT. phdrag) CLOSE (3333)
        !
           !
        ENDDO ! idope
  ENDDO  !itemp
  !
  ENDIF ! inode
  !
  !
  !
  CALL mp_barrier (inter_pool_comm)
  !
  !
  !
  DEALLOCATE (xk)
  DEALLOCATE (xq)
  DEALLOCATE (wf)
  DEALLOCATE (iq_ful)
  DEALLOCATE (iq_gam)
  !
END SUBROUTINE phcheck_scat
