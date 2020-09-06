!--------------------------------------------------------------------------
MODULE epw_explore
!-------------------------------------------------------------------------- 
  !
  USE kinds, ONLY : DP
  !
  SAVE
  !
  INTEGER                    :: &
          ignph_num,            &
          ignph_mode(100),      &
          channel_nb,           &
          channel_band(100),    &
          channel_nk
  !
  REAL(KIND=DP)              :: &
          channel_xk(3,100)
  !
  LOGICAL                    :: &
          ignph,                &
          channel
  !
!-------------------------------------------------------------------------- 
END MODULE epw_explore
!-------------------------------------------------------------------------- 


!-----------------------------------------------------------------------
SUBROUTINE explore_readin
!-----------------------------------------------------------------------
  !
  USE epw_explore,   ONLY : ignph, ignph_mode, ignph_num, &
                            channel, channel_band, channel_xk, channel_nb, channel_nk
  USE epwcom,        ONLY : neptemp
  USE io_global,     ONLY : stdout, ionode_id
#ifdef __PARA
  USE mp,            ONLY : mp_bcast
  USE mp_world,      ONLY : world_comm
  USE mp_global,     ONLY : my_pool_id, me_pool
#ENDIF
#ifdef __NAG
  USE F90_UNIX_ENV,  ONLY : IARGC, GETARG
#ENDIF
  !
  IMPLICIT NONE
  !
#ifndef __NAG
  INTEGER :: IARGC
#ENDIF
  INTEGER :: ios
  INTEGER :: i, ik, imode, ibnd, jbnd
  !
  ! initialization
  NAMELIST / input_explore / &
             ignph, ignph_mode, &
             channel, channel_band, channel_xk
  !
#ifdef __PARA
  IF (me_pool .NE. 0 .OR. my_pool_id .NE. 0) GOTO 400 
#ENDIF
  !
  ! default
  ignph           = .FALSE.
  ignph_mode(:)   = 0
  channel         = .FALSE.
  channel_band(:) = 0
  channel_xk(:,:) = -1000.0d0
  !
  !
#ifdef CRAYY
  READ (5, input_explore)
  ios = 0
#ELSE
  READ (5, input_explore, ERR = 200, IOSTAT = ios)
#ENDIF
200 CALL errore ('epw_readin','reading explore_readin namelist', ABS(ios))
  !
  !
  ! ignph
  ignph_num = 0
  IF (ignph .EQ. .TRUE.) THEN
     !
     DO imode = 1, 100
        IF (ignph_mode(imode) .GT. 0) ignph_num = ignph_num + 1
     ENDDO
     !
     DO imode = 1, ignph_num  
        WRITE (stdout,'(5x,a,i4)') 'Ignored phonon mode ', ignph_mode(imode)
     ENDDO
     !
  ENDIF
  !
  !
  ! channel
  IF (neptemp .GT. 1 .AND. channel .EQ. .TRUE.) CALL errore ('epw_readin','multiple temperature can be only used when channel = .FALSE.',1)
  !
  channel_nb = 0
  channel_nk = 0
  IF (channel .EQ. .TRUE.) THEN
     !
     DO ibnd = 1, 100
        IF (channel_band(ibnd) .GT. 0) channel_nb = channel_nb + 1
     ENDDO
     !
     DO ik = 1, 100
        IF (channel_xk(1,ik) .GT. -1000.0d0 .AND. &
            channel_xk(2,ik) .GT. -1000.0d0 .AND. &
            channel_xk(3,ik) .GT. -1000.0d0) channel_nk = channel_nk + 1
     ENDDO
     !
     DO ibnd = 1, channel_nb
        WRITE (stdout,'(5x,a,i4)') 'Scattering channel band ', channel_band(ibnd)
     ENDDO
     !
     DO ik = 1, channel_nk
        WRITE (stdout,'(5x,a,3f13.8)') 'Scattering channel k-point ', channel_xk(1:3,ik)
     ENDDO
     !
  ENDIF
  !
  !
  ! bcast
400 CONTINUE
#ifdef __PARA
  CALL mp_bcast (ignph, ionode_id, world_comm)
  CALL mp_bcast (ignph_mode, ionode_id, world_comm)
  CALL mp_bcast (ignph_num, ionode_id, world_comm)
  !
  CALL mp_bcast (channel, ionode_id, world_comm)
  CALL mp_bcast (channel_band, ionode_id, world_comm)
  CALL mp_bcast (channel_xk, ionode_id, world_comm)
  CALL mp_bcast (channel_nb, ionode_id, world_comm)
  CALL mp_bcast (channel_nk, ionode_id, world_comm)
#ENDIF
  !
END SUBROUTINE explore_readin



