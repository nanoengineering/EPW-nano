  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2001-2008 Quantum-Espresso group
  ! This file is distributed under the terms of the
  ! GNU General Public License. See the file `License'
  ! in the root directory of the present distribution,
  ! or http://www.gnu.org/copyleft/gpl.txt .
  !
  !-----------------------------------------------------------------------
  SUBROUTINE rgd_blk (nr1,nr2,nr3,nat,dyn,ntemp,ndope,ddyn_l,q,tau,epsil,zeu,bg,omega,sign)
  !-----------------------------------------------------------------------
  ! This is adapted from QE PH/rigid.f90 
  !
  ! compute the rigid-ion (long-range) term for q 
  ! The long-range term used here, to be added to or subtracted from the
  ! dynamical matrices, is exactly the same of the formula introduced
  ! in:
  ! X. Gonze et al, PRB 50. 13035 (1994) . Only the G-space term is 
  ! implemented: the Ewald parameter alpha must be large enough to 
  ! have negligible r-space contribution
  !
  USE kinds, only: dp
  USE cell_base, ONLY : alat
  USE epwcom,    ONLY : vg_ph, L_D, screen_polar
  use pwcom, only : tpiba2, tpiba
  USE constants_epw, only: pi, fpi, e2, ci
  implicit none
  integer ::  nr1, nr2, nr3    !  FFT grid
  integer ::  nat              ! number of atoms 
  integer ::  ntemp, ndope
  complex(DP) :: dyn(3*nat,3*nat,ntemp,ndope) ! dynamical matrix
  real(DP) :: &
       q(3),           &! q-vector
       tau(3,nat),     &! atomic positions
       epsil(3,3),     &! dielectric constant tensor
       zeu(3,3,nat),   &! effective charges tensor
       bg(3,3),        &! reciprocal lattice basis vectors
       omega,          &! unit cell volume
       sign             ! sign=+/-1.0 ==> add/subtract rigid-ion term
  !
  ! THL-
  COMPLEX(KIND=DP) :: ddyn_l(3,3*nat,3*nat)
  INTEGER          :: is
  REAL(KIND=DP)    :: dgeg(3)                   
  ! THL---
  !
  ! local variables
  !
  real(DP) :: geg                    !  <q+G| epsil | q+G>
  real(DP) :: sgeg, dielec
  integer :: na,nb, i,j, m1, m2, m3, itemp, idope
  integer :: nrx1, nrx2, nrx3
  real(DP) :: alph, fac,g1,g2,g3, facgd(ntemp,ndope), arg, gmax
  real(DP) :: zag(3),zbg(3),zcg(3), fnat(3)
  complex(dp) :: facg(ntemp,ndope)
  real(DP) :: eps=1.0d-6
  !
  ! alph is the Ewald parameter, geg is an estimate of G^2
  ! such that the G-space sum is convergent for that alph
  ! very rough estimate: geg/4/alph > gmax = 14 
  ! (exp (-14) = 10^-6)
  !
  gmax= 14.d0
  alph= 1.0d0
  geg = gmax*alph*4.0d0
  !
  ! Estimate of nrx1,nrx2,nrx3 generating all vectors up to G^2 < geg
  ! Only for dimensions where periodicity is present, e.g. if nr1=1 
  ! and nr2=1, then the G-vectors run along nr3 only.
  ! (useful if system is in vacuum, e.g. 1D or 2D)
  !
  IF (nr1 == 1) THEN 
     nrx1=0
  ELSE
     nrx1 = int ( sqrt (geg) / &
                  sqrt (bg (1, 1) **2 + bg (2, 1) **2 + bg (3, 1) **2) ) + 1
  ENDIF
  IF (nr2 == 1) THEN 
     nrx2=0
  ELSE
     nrx2 = int ( sqrt (geg) / &
                  sqrt (bg (1, 2) **2 + bg (2, 2) **2 + bg (3, 2) **2) ) + 1
  ENDIF
  IF (nr3 == 1) THEN 
     nrx3=0
  ELSE
     nrx3 = int ( sqrt (geg) / &
                  sqrt (bg (1, 3) **2 + bg (2, 3) **2 + bg (3, 3) **2) ) + 1
  ENDIF
  !
  IF (abs(abs(sign) - 1.0) > eps) &
       CALL errore ('rgd_blk',' wrong value for sign ',1)
  !
!ERROR: here [dielec] only works for isotropic materials
  dielec = (epsil(1,1) + epsil(2,2) + epsil(3,3)) / 3.d0
  fac = sign*e2*fpi/omega
  !
  DO m1 = -nrx1,nrx1
  DO m2 = -nrx2,nrx2
  DO m3 = -nrx3,nrx3
     !
     g1 = m1*bg(1,1) + m2*bg(1,2) + m3*bg(1,3)
     g2 = m1*bg(2,1) + m2*bg(2,2) + m3*bg(2,3)
     g3 = m1*bg(3,1) + m2*bg(3,2) + m3*bg(3,3)
     !
     geg = (g1*(epsil(1,1)*g1+epsil(1,2)*g2+epsil(1,3)*g3)+      &
            g2*(epsil(2,1)*g1+epsil(2,2)*g2+epsil(2,3)*g3)+      &
            g3*(epsil(3,1)*g1+epsil(3,2)*g2+epsil(3,3)*g3))
     !
     ! -- Consider the carrier screening effect
     !
     facgd = 0.d0
     if ((.not. screen_polar) .or. (sign < 0)) then
        if (geg > 0.d0) &
           facgd = fac*exp(-geg/alph/4.0d0)/geg
     else
!        write(stdout,*) 'L_D_check:'
        do itemp = 1, ntemp
           do idope = 1, ndope
!              write(stdout,*) itemp,ief,' - L_D =', L_D(itemp,ief)*0.529177d0, 'Ang'
!              write(stdout,*) '           - g^2 =',sqrt(g1**2+g2**2+g3**2)*tpiba

              sgeg  = ((g1**2+g2**2+g3**2) * tpiba2 + 1.d0/(L_D(itemp,idope)**2.d0)) * dielec / tpiba2
              ! only consider screening for non-zero G point
              facgd(itemp,idope) = fac*exp(-geg/alph/4.0d0)/sgeg
           enddo
        enddo
     endif

     IF (geg > 0.0_DP .and. geg/alph/4.0_DP < gmax ) THEN
        !
!        facgd = fac*exp(-geg/alph/4.0d0)/geg
        !
        DO na = 1,nat
           zag(:)=g1*zeu(1,:,na)+g2*zeu(2,:,na)+g3*zeu(3,:,na)
           fnat(:) = 0.d0
           DO nb = 1,nat
              arg = 2.d0*pi* (g1 * (tau(1,na)-tau(1,nb))+             &
                              g2 * (tau(2,na)-tau(2,nb))+             &
                              g3 * (tau(3,na)-tau(3,nb)))
              zcg(:) = g1*zeu(1,:,nb) + g2*zeu(2,:,nb) + g3*zeu(3,:,nb)
              fnat(:) = fnat(:) + zcg(:)*cos(arg)
           ENDDO
           DO j=1,3 
              DO i=1,3 
                 dyn( (na-1)*3+i,(na-1)*3+j,:,: ) = dyn((na-1)*3+i,(na-1)*3+j,:,:) - facgd(:,:) * &
                                  zag(i) * fnat(j) 
              ENDDO
           ENDDO
        ENDDO
     ENDIF
     !
     g1 = g1 + q(1)
     g2 = g2 + q(2)
     g3 = g3 + q(3)
     !
     geg = (g1*(epsil(1,1)*g1+epsil(1,2)*g2+epsil(1,3)*g3)+      &
            g2*(epsil(2,1)*g1+epsil(2,2)*g2+epsil(2,3)*g3)+      &
            g3*(epsil(3,1)*g1+epsil(3,2)*g2+epsil(3,3)*g3))
     !
     ! -- Consider the carrier screening effect
     !
     facgd = 0.d0
     if ((.not. screen_polar) .or. (sign < 0)) then
        if (geg > 0.d0) &
           facgd = fac*exp(-geg/alph/4.0d0)/geg
     else
        do itemp = 1, ntemp
           do idope = 1, ndope
              sgeg  = ((g1**2+g2**2+g3**2) * tpiba2 + 1.d0/(L_D(itemp,idope)**2.d0)) * dielec / tpiba2
              ! only consider screening for non-zero G point
              facgd(itemp,idope) = fac*exp(-geg/alph/4.0d0)/sgeg
           enddo
        enddo
     endif
     !
     IF (ABS(sign-1.0d0) .LE. eps) THEN
        !
        dgeg(1) = g1*(epsil(1,1)+epsil(1,1)) + g2*(epsil(1,2)+epsil(2,1)) + g3*(epsil(1,3)+epsil(3,1))
        dgeg(2) = g1*(epsil(2,1)+epsil(1,2)) + g2*(epsil(2,2)+epsil(2,2)) + g3*(epsil(2,3)+epsil(3,2))
        dgeg(3) = g1*(epsil(3,1)+epsil(1,3)) + g2*(epsil(3,2)+epsil(2,3)) + g3*(epsil(3,3)+epsil(3,3))
        !
     ENDIF

     IF (geg > 0.0_DP .and. geg/alph/4.0_DP < gmax ) THEN
        !
!        facgd = fac*exp(-geg/alph/4.0d0)/geg
        !
        DO nb = 1,nat
           zbg(:)=g1*zeu(1,:,nb)+g2*zeu(2,:,nb)+g3*zeu(3,:,nb)
           DO na = 1,nat
              zag(:)=g1*zeu(1,:,na)+g2*zeu(2,:,na)+g3*zeu(3,:,na)
              arg = 2.d0*pi* (g1 * (tau(1,na)-tau(1,nb))+             &
                              g2 * (tau(2,na)-tau(2,nb))+             &
                              g3 * (tau(3,na)-tau(3,nb)))
              !
              facg(:,:) = facgd(:,:) * CMPLX(cos(arg),sin(arg),kind=DP)
              ! NOTE: ddyn_l does not include carrier screening effect
              DO j=1,3 
                 DO i=1,3 
                    dyn( (na-1)*3+i,(nb-1)*3+j,:,: ) = dyn((na-1)*3+i,(nb-1)*3+j,:,:) + facg(:,:) *      &
                                     zag(i) * zbg(j)
                    !
                    ! THL
                    IF (vg_ph .EQ. 'matrix' .AND. ABS(sign-1.0d0) .LE. eps) THEN
                       DO is = 1, 3
                          !
                          ddyn_l(is,(na-1)*3+i,(nb-1)*3+j) = ddyn_l(is,(na-1)*3+i,(nb-1)*3+j) +                  &
                                                             facg(1,1) * (                                            &
                                                             zbg(j)*zeu(is,i,na) + zag(i)*zeu(is,j,nb)        +  &
                                                             zag(i)*zbg(j)*ci*alat*(tau(is,na)-tau(is,nb))    -  &
                                                             zag(i)*zbg(j)*(dgeg(is)/alph/4.0d0+dgeg(is)/geg)    )          
                          !
                       ENDDO
                    ENDIF
                    !
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDIF
  !
  ENDDO
  ENDDO
  ENDDO
  !
END SUBROUTINE rgd_blk
!
!
!-------------------------------------------------------------------------------
SUBROUTINE rgd_blk_epw(q,uq,epmat,nmodes,epsil,zeu,bmat,sign,aaa)
!-------------------------------------------------------------------------------
  !
  !  Compute the long range term for the e-ph vertex
  !  to be added or subtracted from the vertex;
  !  only the major G contribution is implemented
  !
  USE kinds, only: dp
  USE pwcom, ONLY : at, bg, omega, alat
  USE ions_base, ONLY : amass, tau, nat, ntyp => nsp, ityp
  USE constants_epw, only: twopi, fpi, e2, ci, czero, cone, two, ryd2mev
  !
  implicit none
  !
  ! input variables
  !
  integer :: &
       nmodes
  complex(DP) :: &
       uq(nmodes, nmodes),      &! phonon eigenvec
       epmat(nmodes),     &! e-ph matrix
       bmat                ! <u_mk+q|u_nk>
  real(DP) :: &
       q(3),           &! q-vector
       epsil(3,3),     &! dielectric constant tensor
       zeu(3,3,nat),   &! effective charges tensor
       sign             ! sign=+/-1.0 ==> add/subtract long range term
  !
  ! work variables
  !
  real(DP) :: &
       qeq,            &! <q+G| epsil | q+G>
       arg, zaq,       &
       g1, g2, g3, gmax, alph, geg
  integer :: na, nb, ipol, im, m1, m2, m3
  complex(dp) :: fac, facqd, facq, aaa(nmodes)
  !
  IF (abs(sign) /= 1.0) &
       CALL errore ('rgd_blk',' wrong value for sign ',1)
  !
  aaa=czero
  gmax= 14.d0
  alph= 1.0d0
  geg = gmax*alph*4.0d0
  fac = sign*e2*fpi/omega * ci
  !
  call cryst_to_cart (1, q, at, -1)
  m1=-nint(q(1))
  m2=-nint(q(2))
  m3=-nint(q(3)) 
  CALL cryst_to_cart (1, q, bg, 1)
  !
  g1 = m1*bg(1,1) + m2*bg(1,2) + m3*bg(1,3) + q(1)
  g2 = m1*bg(2,1) + m2*bg(2,2) + m3*bg(2,3) + q(2)
  g3 = m1*bg(3,1) + m2*bg(3,2) + m3*bg(3,3) + q(3)
  !
  qeq = (g1*(epsil(1,1)*g1+epsil(1,2)*g2+epsil(1,3)*g3 )+      &
         g2*(epsil(2,1)*g1+epsil(2,2)*g2+epsil(2,3)*g3 )+      &
         g3*(epsil(3,1)*g1+epsil(3,2)*g2+epsil(3,3)*g3 ))*twopi/alat
  !
  facqd = fac*exp(-qeq/alph/4.0d0)/qeq
  !
  DO na = 1,nat
     arg = -twopi* ( g1*tau(1,na)+ g2*tau(2,na)+ g3*tau(3,na) )
     facq = facqd * CMPLX(cos(arg),sin(arg),kind=DP)
     DO ipol=1,3
        zaq=g1*zeu(1,ipol,na)+g2*zeu(2,ipol,na)+g3*zeu(3,ipol,na)
        !
        DO im=1,nmodes 
           epmat(im) = epmat(im) +   &
                 facq * zaq * uq(3*(na-1)+ipol,im)*bmat
           aaa(im) = aaa(im) + facq*zaq*uq(3*(na-1)+ipol,im)*bmat
        ENDDO
        !
     ENDDO !ipol
  ENDDO !nat
  !
  !
END SUBROUTINE rgd_blk_epw
!
!
!-------------------------------------------------------------------------------
SUBROUTINE rgd_blk_epw2(nr1,nr2,nr3,q,uq,epmat,nmodes,epsil,zeu,bmat,sign)
!-------------------------------------------------------------------------------
  !
  !  Compute the long range term for the e-ph vertex
  !  to be added or subtracted from the vertex
  !
  USE kinds, only: dp
  USE pwcom, ONLY : at, bg, omega, alat
  USE ions_base, ONLY : amass, tau, nat, ntyp => nsp, ityp
  USE constants_epw, only: twopi, fpi, e2, ci, czero, cone, two, ryd2mev
  !
  implicit none
  !
  ! input variables
  !
  integer :: &
       nmodes, &
       nr1, nr2, nr3
  complex(DP) :: &
       uq(nmodes, nmodes),      &! phonon eigenvec
       epmat(nmodes),     &! e-ph matrix
       bmat                ! <u_mk+q|u_nk>
  real(DP) :: &
       q(3),           &! q-vector
       epsil(3,3),     &! dielectric constant tensor
       zeu(3,3,nat),   &! effective charges tensor
       sign             ! sign=+/-1.0 ==> add/subtract long range term
  !
  ! work variables
  !
  real(DP) :: &
       qeq,            &! <q+G| epsil | q+G>
       arg, zaq,       &
       g1, g2, g3, gmax, alph, geg, aaa1
  integer :: na, nb, ipol, im, m1,m2,m3, nrx1,nrx2,nrx3
  complex(dp) :: fac, facqd, facq, aaa(nmodes)
  !
  IF (abs(sign) /= 1.0) &
       CALL errore ('rgd_blk',' wrong value for sign ',1)
  !
  aaa=czero
  gmax= 14.d0
  alph= 1.0d0
  geg = gmax*alph*4.0d0
  fac = sign*e2*fpi/omega * ci
  !
  IF (nr1 == 1) THEN 
     nrx1=0
  ELSE
     nrx1 = int ( sqrt (geg) / &
                  sqrt (bg (1, 1) **2 + bg (2, 1) **2 + bg (3, 1) **2) ) + 1
  ENDIF
  IF (nr2 == 1) THEN 
     nrx2=0
  ELSE
     nrx2 = int ( sqrt (geg) / &
                  sqrt (bg (1, 2) **2 + bg (2, 2) **2 + bg (3, 2) **2) ) + 1
  ENDIF
  IF (nr3 == 1) THEN 
     nrx3=0
  ELSE
     nrx3 = int ( sqrt (geg) / &
                  sqrt (bg (1, 3) **2 + bg (2, 3) **2 + bg (3, 3) **2) ) + 1
  ENDIF
  !
  DO m1 = -nrx1,nrx1
  DO m2 = -nrx2,nrx2
  DO m3 = -nrx3,nrx3
  !
  g1 = m1*bg(1,1) + m2*bg(1,2) + m3*bg(1,3) + q(1)
  g2 = m1*bg(2,1) + m2*bg(2,2) + m3*bg(2,3) + q(2)
  g3 = m1*bg(3,1) + m2*bg(3,2) + m3*bg(3,3) + q(3)
  !
  qeq = (g1*(epsil(1,1)*g1+epsil(1,2)*g2+epsil(1,3)*g3 )+      &
         g2*(epsil(2,1)*g1+epsil(2,2)*g2+epsil(2,3)*g3 )+      &
         g3*(epsil(3,1)*g1+epsil(3,2)*g2+epsil(3,3)*g3 )) !*twopi/alat
  !
  IF (qeq > 0.0_DP .and. qeq/alph/4.0_DP < gmax ) THEN
  !
  qeq=qeq*twopi/alat
  facqd = fac*exp(-qeq/alph/4.0d0)/qeq!/(two*wq)
  !
  DO na = 1,nat
     arg = -twopi* ( g1*tau(1,na)+ g2*tau(2,na)+ g3*tau(3,na) )
     facq = facqd * CMPLX(cos(arg),sin(arg),kind=DP)
     DO ipol=1,3
        zaq=g1*zeu(1,ipol,na)+g2*zeu(2,ipol,na)+g3*zeu(3,ipol,na)
        !
        DO im=1,nmodes 
           epmat(im) = epmat(im) +   &
                 facq * zaq * uq(3*(na-1)+ipol,im)*bmat
           aaa(im) = aaa(im) + facq*zaq*uq(3*(na-1)+ipol,im)*bmat
        ENDDO
        !
     ENDDO !ipol
  ENDDO !nat
  ENDIF
  !
  ENDDO
  ENDDO
  ENDDO
  !
  !
END SUBROUTINE rgd_blk_epw2
!
!
!-------------------------------------------------------------------------------
SUBROUTINE rgd_blk_epw3(uq,epmat,ntemp,ndope,ep_polar,nmodes,bmat,sign)
!-------------------------------------------------------------------------------
  USE kinds, only: DP
  USE ions_base, ONLY : nat

  implicit none
  integer :: &
       nmodes, ntemp, ndope
  complex(DP) :: &
       uq(nmodes, nmodes,ntemp,ndope),      &! phonon eigenvec
       epmat(nmodes,ntemp,ndope),     &! e-ph matrix
       ep_polar(nmodes,ntemp,ndope),  &
       bmat                ! <u_mk+q|u_nk>
  real(DP) :: &
       sign             ! sign=+/-1.0 ==> add/subtract long range term

  integer :: na, ipol, im, itemp, idope
  !
  !  Compute the long range term for the e-ph vertex
  !  to be added or subtracted from the vertex
  !
  do itemp = 1, ntemp
     do idope = 1, ndope
        DO na = 1,nat
           DO ipol = 1,3
              DO im = 1,nmodes 
                 epmat(im,itemp,idope) = epmat(im,itemp,idope) + ep_polar(3*(na-1)+ipol,itemp,idope) * uq(3*(na-1)+ipol,im,itemp,idope)*bmat
              ENDDO
           ENDDO !ipol
        ENDDO !nat
     enddo
  enddo
  !
END SUBROUTINE rgd_blk_epw3
!
!
!-----------------------------------------------------------------------------
SUBROUTINE dyndia_epw (nmodes, xq, dyn, uq, wq)
!-----------------------------------------------------------------------------
!
#include "f_defs.h"
  USE kinds, ONLY : DP
  USE pwcom, ONLY : at, bg, celldm, omega
  USE ions_base, ONLY : amass, tau, nat, ntyp => nsp, ityp
  USE constants_epw, ONLY : twopi, ci, czero, rydcm1
  USE io_global, ONLY : stdout
  !
  implicit none
  !
  ! input variables
  integer :: nmodes
  ! number of branches
  real(kind=DP) :: xq(3)
  complex(kind=DP) :: dyn (nmodes, nmodes)
  ! dynamical matrix in bloch representation (Cartesian coordinates)
  !
  ! output variables
  complex(kind=DP) uq(nmodes,nmodes)
  real(kind=DP) wq(nmodes)
  !
  ! work variables 
  complex(kind=DP) dyn1(nmodes, nmodes), dynp( nmodes*(nmodes+1)/2 )
  complex(kind=DP) :: cwork( 2*nmodes ), u1(nmodes,nmodes)
  integer :: neig, info, ifail( nmodes ), iwork( 5*nmodes )
  real(kind=DP) :: w1(nmodes), rwork( 7*nmodes ), massfac
  integer :: imode, jmode, na, nb, i, mu, nu
  !
  ! ------------------------------------------------------------------
  ! diagonalize dyn mat --> uq, wq
  ! ------------------------------------------------------------------
  !
  ! divide by the square root of masses 
  !
  DO na = 1, nat
   DO nb = 1, nat
     massfac = 1.d0 / sqrt ( amass(ityp(na)) * amass(ityp(nb)) )
     dyn1(3*(na-1)+1:3*na, 3*(nb-1)+1:3*nb) = &
     dyn (3*(na-1)+1:3*na, 3*(nb-1)+1:3*nb) * massfac
   END DO
  END DO
  !
  DO jmode = 1, nmodes
   DO imode = 1, jmode
      dynp (imode + (jmode - 1) * jmode/2 ) = &
      ( dyn1 ( imode, jmode) + conjg ( dyn1 ( jmode, imode) ) ) / 2.d0
   ENDDO
  ENDDO
  !
  ! diagonalize
  !
  CALL zhpevx ('V', 'A', 'U', nmodes, dynp , 0.0, 0.0, &
               0, 0,-1.0, neig, w1, u1, nmodes, cwork, &
               rwork, iwork, ifail, info)
  !
  ! frequencies and mass-scaled eigenvectors
  !
  DO nu = 1, nmodes
    IF ( w1 (nu) .gt. 0.d0 ) then
       wq(nu) =  sqrt(abs( w1 (nu) ))
    ELSE
       wq(nu) = -sqrt(abs( w1 (nu) ))
    ENDIF
    DO mu = 1, nmodes
       na = (mu - 1) / 3 + 1
       uq (mu, nu) = u1 (mu, nu) / sqrt(amass(ityp(na)))
    ENDDO
  ENDDO
  !WRITE ( stdout, '(5x,a,3f10.6)') "Frequencies, q = ", (xq(i),i=1,3)
  !WRITE ( stdout, '(6(2x,f10.5))' ) (wq(nu)*rydcm1, nu = 1, nmodes)
  !WRITE(6,'("uq")')
  !DO nu=1,nmodes
   !WRITE(6,'(2f12.8)') (u1(mu,nu), mu=1,nmodes)
   !WRITE(6,'(/,5X)')
  !ENDDO
  !
  !
END SUBROUTINE dyndia_epw
!
!
!-----------------------------------------------------------------------
SUBROUTINE compute_bmn_para2 (nbnd, nkstot, cufkk, cufkq, bmatf)
!-----------------------------------------------------------------------
!
! calculates overlap U_k+q U_k^\dagger
!
#include "f_defs.h"
  USE kinds,     ONLY : DP
  USE constants_epw, ONLY : twopi, ci, czero, cone
  USE mp_global, ONLY : my_pool_id
  implicit none
  !
  !  input variables
  integer :: nbnd, nkstot
  ! number of bands (possibly in the optimal subspace)
  ! number of kpoint blocks, total
  complex(kind=DP) :: cufkk (nbnd, nbnd), cufkq (nbnd, nbnd)
  ! rotation matrix U(k)
  ! rotation matrix U(k+q)
  !
  !  output variables
  complex(kind=DP) :: bmatf (nbnd, nbnd)
  ! overlap wfcs in Bloch representation, fine grid
  !
  !  work variables 
  complex(kind=DP) :: btmp( nbnd, nbnd)
  !
  !  every pool works with its own subset of k points on the fine grid
  !
  bmatf = czero
  !
  !  U_q^ (k') U_k^dagger (k')
  !
  CALL zgemm ('n', 'c', nbnd, nbnd, nbnd, cone, cufkq, &
       nbnd, cufkk, nbnd, czero, bmatf, nbnd)
  !
  !bmatf = bmatf / dble(nkstot)
  !
  !
END SUBROUTINE compute_bmn_para2
!
!-----------------------------------------------------------------------
SUBROUTINE compute_bmn_para3 (nbnd, nbndsub, nks, cuk, cukq, bmat)
!-----------------------------------------------------------------------
!
! calculates <u_mk+q|e^iGr|u_nk> in the approximation q+G->0
!
#include "f_defs.h"
  USE kinds,     ONLY : DP
  USE constants_epw, ONLY : twopi, ci, czero, cone
  USE mp_global, ONLY : my_pool_id
  USE io_global, ONLY : stdout
  implicit none
  !
  !  input variables
  integer :: nbnd, nbndsub, nks
  ! number of bands
  ! number of kpoint blocks, per pool
  complex(kind=DP) :: cuk (nbnd,nbndsub,nks), cukq (nbnd,nbndsub,nks)
  ! rotation matrix U(k)
  ! rotation matrix U(k+q)
  !
  !  output variables
  complex(kind=DP) :: bmat (nbnd, nbnd, nks)
  ! overlap wfcs in Bloch representation, fine grid
  !
  !  work variables 
  integer :: ik,ibnd,jbnd
  !
  !  every pool works with its own subset of k points on the fine grid
  !
  bmat = czero
  !
  !  U_q^ (k') U_k^dagger (k')
  !
  DO ik=1,nks
     CALL zgemm ('n', 'c', nbnd, nbnd, nbndsub, cone, cukq(:,:,ik), &
          nbnd, cuk(:,:,ik), nbnd, czero, bmat(:,:,ik), nbnd)
    !DO ibnd=1,nbnd
    !DO jbnd=1,nbnd
      !write(stdout,'(3i4,2f18.12)') jbnd,ibnd,ik,bmat(jbnd,ibnd,ik)
    !ENDDO
    !ENDDO
  ENDDO
  !
  WRITE(stdout,'(5x,a)') 'BMN calculated'
  WRITE(stdout,*)
  !
  !
END SUBROUTINE compute_bmn_para3
!
!
!-----------------------------------------------------------------
SUBROUTINE polar_eph (uact, xq0, eph_vertex, ntemp, ndope)
!-----------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE ions_base, ONLY : nat, tau
  USE pwcom,     ONLY : tpiba2, tpiba, e2, ngms, tpi, &
                        g, fpi, omega, bg
  use epwcom,    only : L_D, screen_polar
  USE elph2,     ONLY : epsi, zstar
  use io_global, only : stdout

  implicit none

  integer :: ntemp, ndope
  complex(DP) :: uact (3*nat, 3*nat), eph_vertex (3*nat,ntemp,ndope), &
                 eigqts, eigts
  real(DP) :: xq0(3)
  !
  integer :: na, mu, ig, imode, i, j, k, ipol, ibnd, ir, &
             itemp, idope
  real(DP) :: fac, gmg, geg, alpha, g0(3), arg, &
              vpolarq(ntemp,ndope), dielec, sgeg
  complex(DP) :: gtau, fact, u1, u2, u3, g1, g2, g3, uzg
  complex(DP), allocatable :: aux1_c(:), aux2_c(:)
  !
  eph_vertex = 0.d0
  !
  ! ERROR: this is only true for isotropic materials  
  dielec = (epsi(1,1) + epsi(2,2) + epsi(3,3)) / 3.d0

  ! convergence factor
  alpha = 1.d0 * tpiba2
  fac = e2 / tpiba2
  !
  ! substract the long-range perturbed potential due to polar LO phonon
  DO ig = 1, ngms
     g1 = g (1,ig) + xq0 (1)
     g2 = g (2,ig) + xq0 (2)
     g3 = g (3,ig) + xq0 (3)
     gmg = (g1**2 + g2**2 + g3**2)
     geg = (g1*(epsi(1,1)*g1+epsi(1,2)*g2+epsi(1,3)*g3)+ &
            g2*(epsi(2,1)*g1+epsi(2,2)*g2+epsi(2,3)*g3)+ &
            g3*(epsi(3,1)*g1+epsi(3,2)*g2+epsi(3,3)*g3))
     !
     if (.not. screen_polar) then
        vpolarq = - fac * (exp ( - geg * tpiba2 * 0.25d0 / alpha) / geg) * fpi / omega
     else
        do itemp = 1, ntemp
           do idope = 1, ndope
              sgeg  = (gmg * tpiba2 + 1.d0/(L_D(itemp,idope)**2.d0)) * dielec / tpiba2
              vpolarq(itemp,idope) = - fac * (exp ( - geg * tpiba2 * 0.25d0 / alpha) / sgeg) * fpi / omega
           enddo
        enddo
     endif
     !
     DO na = 1, nat
        arg = tpi*(xq0(1)*tau(1,na) + xq0(2)*tau(2,na) + xq0(3)*tau(3,na))
        eigqts = cmplx (cos(arg), -sin(arg))
        fact = tpiba * (0.d0, -1.d0) * eigqts
        mu = 3 * (na - 1)
        !
        arg = tpi*(g(1,ig)*tau(1,na)+g(2,ig)*tau(2,na)+g(3,ig)*tau(3,na))
        gtau = cmplx (cos(arg), -sin(arg))
        !
        DO imode = 1, 3*nat
           !
           IF (abs (uact (mu + 1,imode) ) + abs (uact (mu + 2,imode) ) + abs (uact (mu + &
                3,imode) ) .gt.1.0d-12) THEN
              u1 = uact (mu + 1,imode)
              u2 = uact (mu + 2,imode)
              u3 = uact (mu + 3,imode)
              uzg = (u1*(zstar(1,1,na)*g1+zstar(2,1,na)*g2+zstar(3,1,na)*g3)+ &
                     u2*(zstar(1,2,na)*g1+zstar(2,2,na)*g2+zstar(3,2,na)*g3)+ &
                     u3*(zstar(1,3,na)*g1+zstar(2,3,na)*g2+zstar(3,3,na)*g3))
              !
              eph_vertex(imode,:,:) = eph_vertex(imode,:,:) + vpolarq(:,:) * uzg * fact * gtau
              !
           ENDIF
           !
        enddo
        !
     ENDDO
     !
  ENDDO
  !
END SUBROUTINE
