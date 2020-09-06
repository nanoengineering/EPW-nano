  !
  ! THL: added for computing the weighting of q-points
  !
  !--------------------------------------------------------------------------
  SUBROUTINE wsinit(rws,nrwsx,nrws,atw)
  !--------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER i, ii, ir, jr, kr, nrws, nrwsx, nx
  REAL(KIND=DP) rt, eps, rws(0:3,nrwsx), atw(3,3)
  parameter (eps=1.0d-6,nx=2)
  !
  ii = 1
  DO ir=-nx,nx
     DO jr=-nx,nx
        DO kr=-nx,nx
           DO i=1,3
              rws(i,ii) = atw(i,1)*ir + atw(i,2)*jr + atw(i,3)*kr
           END DO
           rws(0,ii)=rws(1,ii)*rws(1,ii)+rws(2,ii)*rws(2,ii)+            &
                               rws(3,ii)*rws(3,ii)
           rws(0,ii)=0.5d0*rws(0,ii)
           IF (rws(0,ii) .GT. eps) ii = ii + 1
        END DO
     END DO
  END DO
  nrws = ii - 1
  RETURN
  !
  END SUBROUTINE wsinit
  !
  !--------------------------------------------------------------------------
  FUNCTION wsweight(r,rws,nrws)
  !--------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER ir, nreq, nrws
  REAL(KIND=DP) r(3), rrt, ck, eps, rws(0:3,nrws), wsweight
  parameter (eps=1.0d-6)
  !
  wsweight = 0.d0
  nreq = 1
  DO ir =1,nrws
     rrt = r(1)*rws(1,ir) + r(2)*rws(2,ir) + r(3)*rws(3,ir)
     ck = rrt-rws(0,ir)
     IF ( ck .GT. eps ) RETURN
     IF ( abs(ck) .LT. eps ) nreq = nreq + 1
  END DO
  wsweight = 1.d0/DBLE(nreq)
  RETURN
  !
  END FUNCTION wsweight
