SUBROUTINE trilin_interp (vin, nk1, nk2, nk3, nkf1, nkf2, nkf3, vfout)
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN)       :: nk1, nk2, nk3, nkf1, nkf2, nkf3
   REAL(KIND=8), INTENT(IN)  :: vin(nk1*nk2*nk3)
   REAL(KIND=8), INTENT(OUT) :: vfout(nkf1*nkf2*nkf3)
   INTEGER                   :: nk, nkf
   INTEGER                   :: i, j, k, ii, jj, kk, is
   INTEGER                   :: id(8), idf
   REAL(KIND=8)              :: xi, yi, zi, v(8), xb, yb, zb
   !
   !
   nk    = nk1*nk2*nk3
   nkf   = nkf1*nkf2*nkf3
   vfout = 0.0d0
   !
   !
   DO i = 0, nkf1-1
      DO j = 0, nkf2-1      
         DO k = 0, nkf3-1 
            !
            idf = i*nkf2*nkf3 + j*nkf3 + k + 1
            !
            xi = DBLE(i)/DBLE(nkf1)
            yi = DBLE(j)/DBLE(nkf2)
            zi = DBLE(k)/DBLE(nkf3)
            !
            ii = FLOOR(xi*DBLE(nk1))
            jj = FLOOR(yi*DBLE(nk2))
            kk = FLOOR(zi*DBLE(nk3))
            !
            IF (ABS(DBLE(ii)-xi*DBLE(nk1)) .LT. 1.0d-6 .AND. &  
                ABS(DBLE(jj)-yi*DBLE(nk2)) .LT. 1.0d-6 .AND. &  
                ABS(DBLE(kk)-zi*DBLE(nk3)) .LT. 1.0d-6 ) THEN
               !
               id(1) = ii*nk2*nk3 + jj*nk3 + kk + 1
               vfout(idf) = vin(id(1))
               !
            ELSE
               !
               id(1) = MOD(ii  ,nk1)*nk2*nk3 + MOD(jj  ,nk2)*nk3 + MOD(kk  ,nk3) + 1
               id(2) = MOD(ii+1,nk1)*nk2*nk3 + MOD(jj  ,nk2)*nk3 + MOD(kk  ,nk3) + 1             
               id(3) = MOD(ii  ,nk1)*nk2*nk3 + MOD(jj+1,nk2)*nk3 + MOD(kk  ,nk3) + 1              
               id(4) = MOD(ii  ,nk1)*nk2*nk3 + MOD(jj  ,nk2)*nk3 + MOD(kk+1,nk3) + 1
               id(5) = MOD(ii+1,nk1)*nk2*nk3 + MOD(jj  ,nk2)*nk3 + MOD(kk+1,nk3) + 1
               id(6) = MOD(ii  ,nk1)*nk2*nk3 + MOD(jj+1,nk2)*nk3 + MOD(kk+1,nk3) + 1
               id(7) = MOD(ii+1,nk1)*nk2*nk3 + MOD(jj+1,nk2)*nk3 + MOD(kk  ,nk3) + 1
               id(8) = MOD(ii+1,nk1)*nk2*nk3 + MOD(jj+1,nk2)*nk3 + MOD(kk+1,nk3) + 1
               !
               DO is = 1, 8
                  v(is) = vin(id(is))
               ENDDO
               !
               xb = DBLE(ii+1)/DBLE(nk1)
               yb = DBLE(jj+1)/DBLE(nk2)
               zb = DBLE(kk+1)/DBLE(nk3)
               !
               IF (xi .GT. xb .OR. yi .GT. yb .OR. zi .GT. zb) THEN
                  WRITE (*,'(a)') 'Wrong interpolating coefficients'
                  STOP
               ENDIF
               !
               vfout(idf) = v(1)*(1.0d0-xi/xb)*(1.0d0-yi/yb)*(1.0d0-zi/zb) + &
                            v(2)*(      xi/xb)*(1.0d0-yi/yb)*(1.0d0-zi/zb) + &
                            v(3)*(1.0d0-xi/xb)*(      yi/yb)*(1.0d0-zi/zb) + &
                            v(4)*(1.0d0-xi/xb)*(1.0d0-yi/yb)*(      zi/zb) + &
                            v(5)*(      xi/xb)*(1.0d0-yi/yb)*(      zi/zb) + &
                            v(6)*(1.0d0-xi/xb)*(      yi/yb)*(      zi/zb) + &
                            v(7)*(      xi/xb)*(      yi/yb)*(1.0d0-zi/zb) + &
                            v(8)*(      xi/xb)*(      yi/yb)*(      zi/zb)
                !
            ENDIF
            !
         ENDDO ! k
      ENDDO ! j
   ENDDO ! i
   !
   RETURN
   !
   !
END SUBROUTINE trilin_interp



SUBROUTINE lin_interp (vin0_x, vin0_y, nk, nkf, vfin_x, vfout_y)
   !
#ifdef __PARA
  USE mp,        ONLY : mp_barrier, mp_sum, mp_bcast
  USE io_global, ONLY : ionode, stdout
  USE mp_global, ONLY : inter_pool_comm
#ENDIF
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN)       :: nk, nkf
   REAL(KIND=8), INTENT(IN)  :: vin0_x(nk), vin0_y(nk), vfin_x(nkf)
   REAL(KIND=8), INTENT(OUT) :: vfout_y(nkf)
   INTEGER                   :: ik, ikf, ik_star, ik_stop
   INTEGER                   :: imn, mnloc(1)
   REAL(KIND=8)              :: vin_x(nk), vin_y(nk), vin0_x_tmp(nk)
   LOGICAL                   :: already_interp
   !
   !
   vfout_y = 0.0d0
   !
   vin0_x_tmp = vin0_x
   DO ik = 1, nk
      !
      mnloc = MINLOC(vin_x)        
      imn = mnloc(1)
      !
      vin_x(ik) = vin0_x_tmp(imn)
      vin_y(ik) = vin0_y(imn)
      vin0_x_tmp(imn) = +9.9d40
      !
   ENDDO
   !
   !
#ifdef __PARA
   CALL mp_barrier (inter_pool_comm)
#ENDIF
   !
   CALL para_bounds (ik_star, ik_stop, nkf)
   !
   DO ikf = ik_star, ik_stop
      !
      DO ik = 1, nk-1
         !
         already_interp = .FALSE.
         !
         IF (vfin_x(ikf) .GE. vin_x(ik) .AND. vfin_x(ikf) .LT. vin_x(ik+1)) THEN
            !
            vfout_y(ikf) = vin_y(ik) + (vfin_x(ikf)-vin_x(ik)) * (vin_y(ik+1)-vin_y(ik)) / (vin_x(ik+1)-vin_x(ik))
            already_interp = .TRUE.
            !
         ELSEIF (vfin_x(ikf) .LT. vin_x(1)) THEN
            !
            vfout_y(ikf) = vin_y(1) + (vfin_x(ikf)-vin_x(1)) * (vin_y(2)-vin_y(1)) / (vin_x(2)-vin_x(1))
            already_interp = .TRUE.
            !
         ELSEIF (vfin_x(ikf) .GE. vin_x(nk)) THEN
            !
            vfout_y(ikf) = vin_y(nk-1) + (vfin_x(ikf)-vin_x(nk-1)) * (vin_y(nk)-vin_y(nk-1)) / (vin_x(nk)-vin_x(nk-1))
            already_interp = .TRUE.
            !
         ENDIF
         !
         IF (already_interp) EXIT
         !
      ENDDO
      !
   ENDDO
   !
   !
#ifdef __PARA
   CALL mp_barrier (inter_pool_comm)
   CALL mp_sum (vfout_y,inter_pool_comm)
#ENDIF
   !
END SUBROUTINE lin_interp
