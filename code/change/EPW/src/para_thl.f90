!----------------------------------------------------------------------------------------
MODULE para_thl
!----------------------------------------------------------------------------------------
  USE kinds,     ONLY : DP
#ifdef __PARA
  USE mp,        ONLY : mp_bcast, mp_sum
  USE io_global, ONLY : ionode_id
  USE mp,        ONLY : mp_barrier
  USE mp_global, ONLY : my_pool_id, inter_pool_comm
#ENDIF
  !
  IMPLICIT NONE
  !
  PRIVATE
     !
  PUBLIC :: gether_thl
     !
     INTERFACE gether_thl
        !
        MODULE PROCEDURE gether_thl_r2, gether_thl_r3
        !
     END INTERFACE
     !
  CONTAINS
  !
  !----------------------------------------------------------------------------!
  SUBROUTINE gether_thl_r2 (var, var_pol, ik_star, ik_stop, ik_map)
  !----------------------------------------------------------------------------! 
     ! 
     IMPLICIT NONE
     !
     REAL(KIND=DP), INTENT(INOUT) :: var(:,:), var_pol(:,:)
     INTEGER, INTENT(IN)          :: ik_star, ik_stop, ik_map(:)
     INTEGER                      :: ik
     !
     var = 0.0d0
     !
     DO ik = ik_star, ik_stop
        var(:,ik_map(ik)) = var_pol(:,ik-ik_star+1)
     ENDDO
     !
#ifdef __PARA
     CALL mp_barrier (inter_pool_comm)
     CALL mp_sum (var,inter_pool_comm)
#ENDIF
     !
  END SUBROUTINE gether_thl_r2
  !
  !    
  !----------------------------------------------------------------------------!
  SUBROUTINE gether_thl_r3 (var, var_pol, ik_star, ik_stop, ik_map)
  !----------------------------------------------------------------------------! 
     ! 
     IMPLICIT NONE
     !
     REAL(KIND=DP), INTENT(INOUT) :: var(:,:,:), var_pol(:,:,:)
     INTEGER, INTENT(IN)          :: ik_star, ik_stop, ik_map(:)
     INTEGER                      :: ik
     !
     var = 0.0d0
     !
     DO ik = ik_star, ik_stop
        var(:,:,ik_map(ik)) = var_pol(:,:,ik-ik_star+1)
     ENDDO
     !
#ifdef __PARA
     CALL mp_barrier (inter_pool_comm)
     CALL mp_sum (var,inter_pool_comm)
#ENDIF
     !
  END SUBROUTINE gether_thl_r3
  !
  !    
END MODULE para_thl
