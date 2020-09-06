  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !                                                                            
  ! 
  !-----------------------------------------------------------------------
  SUBROUTINE evcbloch2wan ( nbnd, nbndsub, nks, nkstot, xk, &
       cu, evck, nrr, irvec, evcw, wslen, evckw)
  !-----------------------------------------------------------------------
  !
  !  From the electron-phonon matrix elements in Bloch representation (coarse 
  !  mesh), find the corresponding matrix elements in Wannier representation
  !
  !
  !-----------------------------------------------------------------------
  !
#include "f_defs.h"
  !
  USE kinds,     ONLY : DP
  USE pwcom,     ONLY : at, bg, celldm
  USE constants_epw, ONLY : bohr2ang, twopi, ci, czero, cone
  USE noncollin_module, ONLY : noncolin, npol
#ifdef __PARA
  USE io_global, ONLY : ionode_id
  USE mp_global, ONLY : inter_pool_comm,my_pool_id
  USE mp       , ONLY : mp_sum 
  USE mp_world,  ONLY : mpime
#endif
  implicit none
  !
  !  input variables
  !
  integer :: nbnd, nbndsub, nks, nrr, irvec (3, nrr), nkstot
  ! number of bands
  ! number of bands in the optimal subspace
  ! number of kpoints
  ! number of kpoint blocks, in the pool
  ! number of kpoint blocks, total
  ! number of WS points and coordinates
  real(kind=DP) :: xk (3, nks), wslen (nrr)
  ! kpoint coordinates (cartesian in units of 2piba)
  complex(kind=DP) :: cu (nbnd, nbndsub, nks), evck (nbnd, npol, nks)
  ! rotation matrix from wannier code
  ! e-p matrix in bloch representation, coarse mesh
  !
  ! output variables 
  !
  complex(kind=DP) :: evcw (nbndsub, npol, nrr), evckw (nks, npol)
  !  e-p matrix  in wannier basis 
  !
  ! work variables 
  !
  complex(kind=DP) :: evcs (nbndsub, npol, nks), evctmp(1, nbnd)
  !  e-p matrix  in smooth Bloch basis, coarse mesh
  !  e-p matrix, temporary
  !
  integer :: ik, ir
  real(kind=DP) :: rdotk, tmp
  complex(kind=DP) :: cfac
  !
  !
  !----------------------------------------------------------
  !  STEP 1: rotation to optimally smooth Bloch states
  !----------------------------------------------------------
  !
  !  g~ = U_k+q^\dagger g U_k
  !
  !  g   is epmatk (ibnd, jbnd, ik)
  !  g~  is epmats (ibnd, jbnd, ik)
  !
  CALL start_clock ( 'ep: step 1' )
  !
  !
  DO ik = 1, nks
     !
     ! the two zgemm calls perform the following ops:
     ! epmats  = [ cu(ikq)^\dagger * epmatk ] * cu(ikk)
     ! [here we have a size-reduction from nbnd*nbnd to nbndsub*nbndsub] 
     !
     evctmp(1,:) = evck(:,1,ik)
     CALL zgemm ('n', 'n', 1, nbndsub, nbnd, cone, evctmp,     &
                  1, cu(:,:,ik), nbnd, czero, evcs(:,1,ik), 1)

     if (noncolin) then
        evctmp(1,:) = evck(:,2,ik)
        CALL zgemm ('n', 'n', 1, nbndsub, nbnd, cone, evctmp,     &
                     1, cu(:,:,ik), nbnd, czero, evcs(:,2,ik), 1)
     endif
     !
  ENDDO
  !
  CALL stop_clock ( 'ep: step 1' )
  !
!  evckw(1:nks) = evcs(1,1:nks)
  evckw(1:nks,:) = 0.0d0
  !
  !----------------------------------------------------------------------
  !  STEP 2: Fourier transform to obtain matrix elements in wannier basis
  !----------------------------------------------------------------------
  !
  !  g (R) = (1/nkc) sum_k e^{-ikR} g~(k)
  !
  !  epmatw (nbndsub,nbndsub,ir) is g(R)
  !
  CALL start_clock ( 'ep: step 2' )
  !
  evcw (:,:,:) = czero
  !
  ! bring xk in crystal coordinates
  !
  CALL cryst_to_cart (nks, xk, at, -1)
  !
  DO ir = 1, nrr
     !
     DO ik = 1, nks
       !
       rdotk = twopi * dot_product( xk ( :, ik), dble(irvec( :, ir) ))
       cfac = exp( -ci*rdotk ) / dble(nkstot)
       evcw (:,1,ir) = evcw (:,1,ir) + cfac * evcs (:,1,ik)
       if (noncolin) &
          evcw (:,2,ir) = evcw (:,2,ir) + cfac * evcs (:,2,ik)
       !
     ENDDO
     !
  ENDDO
  !
#ifdef __PARA
  CALL mp_sum(evcw,inter_pool_comm)  
#endif
  !
  ! bring xk back into cart coord
  !
  CALL cryst_to_cart (nks, xk, bg, 1)
  !
  !
  !  check spatial decay of matrix elements in Wannier basis
  !  the unit in r-space is angstrom, and I am plotting 
  !  the matrix for the first mode only
  !
#ifdef __PARA
    IF (mpime.eq.ionode_id) THEN
#endif
      OPEN (unit=301,file='decay.evcwan')
      WRITE(301, '(/3x,a/)') '#Spatial decay of e-p matrix elements in Wannier basis'
      DO ir = 1, nrr
        ! 
        tmp =  maxval ( abs(evcw(:,1,ir)) ) 
        WRITE(301, *) wslen(ir) * celldm (1) * bohr2ang, tmp
        !
      ENDDO
      CLOSE(301)
#ifdef __PARA
    ENDIF
#endif
  !
  CALL stop_clock ( 'ep: step 2' )
  !
  END SUBROUTINE evcbloch2wan

  
