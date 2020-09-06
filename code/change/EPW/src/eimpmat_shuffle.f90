  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino  
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                        
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE eimpmat_shuffle ( iq_irr, nqc_irr, iq, gmapsym, eigv, isym, invs0, xq0, timerev )
  !-----------------------------------------------------------------------
  !
  ! Electron-impurity calculation from data saved in dv_tot.dat
  !
  !-----------------------------------------------------------------------
  !
#include "f_defs.h"
  !
#ifdef __PARA
  USE mp,        ONLY : mp_barrier, mp_sum
  USE io_global, ONLY : stdout
  USE mp_global, ONLY : my_pool_id, nproc_pool,npool,kunit,&
                          inter_pool_comm
  USE mp_images, ONLY : nproc_image
  USE pwcom,     ONLY : nkstot
#endif
  !
  USE pwcom,     ONLY : nbnd, ngm, doublegrid, nks
  USE kinds,     ONLY : DP
  USE modes,     ONLY : nmodes, nirr, npert, u
  USE elph2,     ONLY : dvimp, eimpmatq, el_imp_mat, bg_max
  USE constants_epw, ONLY : czero, cone
  USE fft_base,  ONLY : dfftp, dffts
  USE noncollin_module,     ONLY : nspin_mag
!  USE noncollin_module,     ONLY : noncolin
  !
  implicit none
  !
  integer :: irr, imode0, ipert, is, iq, iq_irr, nqc_irr
  ! counter on the representations
  ! counter on the modes
  ! the change of Vscf due to perturbations
  ! the current qpoint in the uniform grid
  ! the current ireducible qpoint
  ! the total number of irreducible qpoints in the list
  logical :: timerev
  !  true if we are using time reversal
  !
  integer :: tmp_pool_id, ik0, ik, ibnd, jbnd
  ! 
#ifdef __PARA
  integer :: iks, nkl, nkr
#endif
  !
  integer :: gmapsym ( ngm, 48 ), isym, invs0 (48)
  ! the correspondence G-->S(G)
  ! the symmetry which generates the current q in the star
  ! the index of the inverse operations
  complex(kind=DP) :: eigv (ngm, 48)
  ! e^{ iGv} for 1...nsym ( v the fractional translation)
  real(kind=DP) :: xq0(3)
  ! the first q in the star (cartesian)
  !
  !
  ik0 = 0
  tmp_pool_id = 0
  !
#ifdef __PARA
  !
  npool =  nproc_image / nproc_pool
  IF (npool.gt.1) THEN
    !
    ! number of kpoint blocks, kpoints per pool and reminder
    kunit = 1 
    nkl   = kunit * ( nkstot / npool )
    nkr   = ( nkstot - nkl * npool ) / kunit
    ! the reminder goes to the first nkr pools
    IF ( my_pool_id < nkr ) nkl = nkl + kunit
    !
    iks = nkl * my_pool_id + 1
    IF ( my_pool_id >= nkr ) iks = iks + nkr * kunit
    !
    !  the index of the first k point block in this pool - 1
    !  (I will need the index of ik, not ikk)
    !
    ik0 = ( iks - 1 ) / kunit
    !
  ENDIF
  !
#endif
  !    
  ! Read impurity potential
  !
  write(stdout,*) 'dfftp%nnr, dffts%nnr =',dfftp%nnr, dffts%nnr
  write(stdout,*) 'nspin_mag =', nspin_mag
  !
  ! Calculate electron-impurity scattering matrix
  !
  CALL elimp_shuffle (dvimp, gmapsym, eigv, isym, invs0, xq0, iq, timerev)
  !
  !
#ifdef __PARA
  CALL mp_barrier(inter_pool_comm)
#endif
  !
  !  the output e-p matrix in the pattern representation
  !  must be transformed in the cartesian basis
  !  epmat_{CART} = conjg ( U ) * epmat_{PATTERN}
  !
  !  note it is not U^\dagger ! Have a look to symdyn_munu.f90 
  !  for comparison
  !
  DO ibnd = 1, nbnd
     DO jbnd = 1, nbnd
        DO ik = 1, nks
           ! 
           eimpmatq (ibnd,jbnd,ik,iq) = el_imp_mat (ibnd,jbnd,ik)
           !
        ENDDO
     ENDDO
  ENDDO
  !
  !
  END SUBROUTINE eimpmat_shuffle


  !---------------------------------------------------------------------
  SUBROUTINE elimp_shuffle (dvimp, gmapsym, eigv, isym, invs, xq0, iq, timerev)
  !---------------------------------------------------------------------
  !
  !      Calculation of the electron-impurity(sr) matrix elements
  !      <\psi(k+q)|dV_imp|\psi(k)>
  !
  !---------------------------------------------------------------------
  !
#include "f_defs.h"
  !
#ifdef __PARA
  USE mp_global,     ONLY : my_pool_id, nproc_pool,    &
                            intra_image_comm, intra_pool_comm, &
                            me_pool, root_pool, inter_pool_comm
  USE mp,            ONLY : mp_barrier, mp_bcast, mp_put,mp_sum
#endif
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE wavefunctions_module,  ONLY: evc
  USE io_files,      ONLY : iunigk, diropn, seqopn
  USE wvfct,         ONLY : npwx
  USE pwcom,         ONLY : current_spin, isk, tpiba, g, &
                            lsda, nbnd, npw, xk, ngm, &
                            igk, nks, omega, bg
  USE uspp,          ONLY : vkb
  USE symm_base,     ONLY : s
  USE modes,         ONLY : u
  USE phcom,         ONLY : iuwfc
  USE qpoint,        ONLY : igkq, xq, npwq
  USE eqv,           ONLY : dvpsi, evq
  USE units_ph,      ONLY : lrwfc
  USE phus,          ONLY : alphap
  USE lrus,          ONLY : becp1
  USE becmod,        ONLY : calbec
  USE elph2,         ONLY : shift, gmap, el_imp_mat, umat, umatq, &
                            umat_all, xk_all, et_all, xkq, etq, &
                            nfftmesh, scell_vol, rp, imp_meshmap, &
                            bg_max, dvGr, map_l2g, map_g2l, dvimp_q, &
                            ng_max
  USE fft_base,      ONLY : dffts
  USE constants_epw, ONLY : czero, cone, ci, twopi
  USE control_flags, ONLY : iverbosity
  USE control_lr,    ONLY : lgamma
  USE klist,         ONLY : nkstot
  USE noncollin_module,     ONLY : noncolin, npol, nspin_mag
  ! 
  implicit none
  !
  real(kind=DP) :: dvimp (nfftmesh)
  !  delta scf potential

  logical :: exst!, addnlcc
  !
  ! work variables
  !
  integer :: ik, ipert, mode, ibnd, jbnd, ig, nkq, ipool, &
       ik0, igkq_tmp (npwx), imap, &
       ipooltmp, nkq_abs, ipol, imesh, iq
  complex(kind=DP), ALLOCATABLE :: aux1 (:,:), eimpmat (:,:), eptmp (:,:), &
                                   aux2(:,:), evc_r(:,:,:), evq_r(:,:,:), &
                                   dvimpevc(:,:), evq_rmesh(:,:), &
                                   dvevc(:,:), eemat(:,:,:), eemat_check(:,:,:), &
                                   eemat_r(:,:,:), eimpmat_r(:,:)
!DBSP - NAG complains ...
  COMPLEX(DP),EXTERNAL :: ZDOTC
  real(kind=DP) :: xktmp(3), sxk(3)
  !
  ! variables for folding of k+q grid
  !
  REAL(kind=DP) :: g0vec_all_r(3,125)

  !   G-vectors needed to fold the k+q grid into the k grid, cartesian coord.
  INTEGER :: ng0vec, ngxx
  !   number of inequivalent such translations
  !   bound for the allocation of the array gmap
  !
  ! variables for rotating the wavefunctions (in order to USE q in the irr
  ! wedge)
  !
  INTEGER :: gmapsym ( ngm, 48 ), isym, invs(48)
  ! the map G->S(G)
  ! the symmetry which generates the current q in the star
  ! index of the inverse operation
  complex(kind=DP) :: eigv (ngm, 48), r0_phase
  real(kind=DP) :: xq0(3)
  real(kind=DP) :: zero_vect(3)
  integer :: nkk, nkk_abs
  !  the fractional traslation
  !  work variable
  !  the first q in the star (cartesian)
  logical :: timerev
  !  true if we are using time reversal
  ! 
  ! JZ: wavefunction check
  integer  :: ir, ir_evc, ind_g, igx, igy, igz, &
              index_igk(dffts%nnr), ig_shift, &
              integration_method
  integer, allocatable :: map_kq2k(:,:)
  real(DP) :: evc_checksum, rho_r, xk_check(3), norm_evc, norm_evq, &
              xk_gamma(3), xk_x(3), xk_l(3), norm_
  complex(DP) :: evc_check(npwx*npol,nbnd), &
                 evc_remap(npwx*npol)
  logical :: xk_find, x_find, gamma_find, l_find

  real(kind=DP) :: xq_cart(3), sc_un_ratio
  real(kind=DP) :: t1, t2, t3, t4, t5, tt1, tt2, tt3, tt4, &
                   t4_1, t4_2, t4_3, t4_4

  !
  IF ( .not. ALLOCATED (eimpmat) ) ALLOCATE ( eimpmat( nbnd, nbnd) )
  IF ( .not. ALLOCATED (eptmp) )   ALLOCATE ( eptmp ( nbnd, nbnd) )
  IF ( .not. ALLOCATED (aux1) )    ALLOCATE ( aux1(dffts%nnr, npol) )
  IF ( .not. ALLOCATED (aux2) )    ALLOCATE ( aux2( npwx*npol, nbnd) )
  IF ( .not. ALLOCATED (eemat) ) ALLOCATE ( eemat( nbnd, nbnd, ng_max) )
  IF ( .not. ALLOCATED (eemat_check) ) ALLOCATE ( eemat_check( nbnd, nbnd, ng_max) )
  IF ( .not. ALLOCATED (eemat_r) ) ALLOCATE ( eemat_r( nbnd, nbnd, ng_max) )
  IF ( .not. ALLOCATED (dvevc) )    ALLOCATE ( dvevc(dffts%nnr, npol) )
  IF ( .not. ALLOCATED (eimpmat_r) ) ALLOCATE ( eimpmat_r( nbnd, nbnd) )
  allocate (map_kq2k (npwx*npol,(2*bg_max+1)**3))
  eimpmat(:,:) = (0.d0,0.d0)
  eptmp(:,:) = (0.d0,0.d0)
  aux1(:,:) = (0.d0,0.d0)
  aux2(:,:) = (0.d0,0.d0)
  dvevc(:,:) = (0.d0,0.d0)
  eimpmat_r(:,:) = (0.d0,0.d0)
  zero_vect = 0.d0
  
  ! added
  IF ( .not. ALLOCATED (evc_r) )    ALLOCATE ( evc_r(dffts%nnr, npol, nbnd) )
  IF ( .not. ALLOCATED (evq_r) )    ALLOCATE ( evq_r(dffts%nnr, npol, nbnd) )
  IF ( .not. ALLOCATED (dvimpevc) )    ALLOCATE ( dvimpevc(nfftmesh, nbnd) )
  IF ( .not. ALLOCATED (evq_rmesh) )    ALLOCATE ( evq_rmesh(nfftmesh, nbnd) )

  !
  IF (ALLOCATED(xkq) ) DEALLOCATE (xkq)
  IF (.not. ALLOCATED(xkq) ) ALLOCATE (xkq (3, nkstot) )
  xkq(:,:) = 0.d0
#ifdef __PARA
  IF (nproc_pool>1) call errore &
    ('elphel2_shuffle', 'ONLY one proc per pool in shuffle mode', 1)
#endif
  IF (.not.lgamma) THEN
     !
     ! setup for k+q folding
     !
     CALL kpointdivision ( ik0 )
     CALL readgmap ( nkstot, ngxx, ng0vec, g0vec_all_r )
     !
     IF (iverbosity.eq.1) WRITE(stdout, 5) ngxx
5    FORMAT (5x,'Estimated size of gmap: ngxx =',i5)
     !
  ENDIF
  !
  ! close all sequential files in order to re-open them as direct access
  ! close all .wfc files in order to prepare shuffled read
  !
  CLOSE (unit = iunigk, status = 'keep')
  CLOSE (unit = iuwfc,  status = 'keep')
#ifdef __PARA
  ! never remove this barrier
  CALL mp_barrier(inter_pool_comm)
#endif
  !
  xk_gamma = (/0.d0, 0.d0, 0.d0/)
  xk_x     = (/0.d0, 1.d0, 0.d0/)
  xk_l     = (/-0.5d0,0.5d0,-0.5d0/)
  !
  tt1 = 0.d0
  tt2 = 0.d0
  tt3 = 0.d0
  tt4 = 0.d0
  !
  DO ik = 1, nks
     !
     call cpu_time(t1)

#ifdef __PARA
     ipooltmp= my_pool_id+1
#endif
     !
     CALL ktokpmq ( xk (:, ik), xq, +1, ipool, nkq, nkq_abs )
     !
     !   we define xkq(:,ik) and etq(:,ik) for the current xq
     !
     IF (ALLOCATED(etq)) DEALLOCATE(etq)
     IF (.not. ALLOCATED(etq) ) ALLOCATE (etq (nbnd, nks) )
     etq(:,:) = 0.d0
     xkq(:, ik) = xk_all(:, nkq_abs)
     etq(:, ik) = et_all(:, nkq_abs)
     !
     ! in serial execution ipool is not USEd in the called SUBROUTINEs, 
     ! in parallel mypool is for k and ipool is for k+q
     !
     CALL readwfc (ipooltmp, ik, evc)
     CALL readigk (ipooltmp, ik, npw, igk)
     !
     CALL readwfc (ipool, nkq, evq)
     CALL readigk (ipool, nkq, npwq, igkq)
     !
#ifdef __PARA
     IF (.not.lgamma .and. nks.gt.1 .and. maxval(igkq(1:npwq)).gt.ngxx) &
          CALL errore ('elphel2_shuffle', 'ngxx too small', 1 )
#endif
     !
     ! ----------------------------------------------------------------
     ! Set the gauge for the eigenstates: unitary transform and phases
     ! ----------------------------------------------------------------
     !
     ! With this option, different compilers and different machines
     ! should always give the same wavefunctions.
     !
     CALL ktokpmq ( xk(:,ik),  zero_vect, +1, ipool, nkk, nkk_abs )
     CALL ktokpmq ( xkq(:,ik), zero_vect, +1, ipool, nkk, nkq_abs )
     IF ( .not. ALLOCATED (umat) )  ALLOCATE ( umat(nbnd,nbnd,nks) )
     IF ( .not. ALLOCATED (umatq) ) ALLOCATE ( umatq(nbnd,nbnd,nks) )
     umat(:,:,ik)  = umat_all(:,:, nkk_abs)
     umatq(:,:,ik) = umat_all(:,:, nkq_abs)
     !
     ! the k-vector needed for the KB projectors
     xktmp = xkq (:, ik)
     !
     ! --------------------------------------------------
     !   Fourier translation of the G-sphere igkq
     ! --------------------------------------------------
     !
     !
     !  Translate by G_0 the G-sphere where evq is defined, 
     !  none of the G-points are lost.
     !
     DO ig = 1, npwq
        imap = ng0vec * ( igkq (ig) - 1 ) + shift( ik+ik0 )
        igkq_tmp (ig) = gmap( imap )
        !  the old matrix version... 
        !  igkq_tmp (ig) = gmap( igkq (ig),  shift( ik+ik0 ) )
     ENDDO
     igkq = igkq_tmp
     !
     !  find k+q from k+q+G_0
     !  (this is needed in the calculation of the KB terms
     !  for nonlocal pseudos)
     !
     xktmp = xkq (:, ik) - g0vec_all_r ( :, shift( ik+ik0 ) )
     !
     !
     ! ---------------------------------------------------------------------
     ! phase factor arising from fractional traslations
     ! ---------------------------------------------------------------------
     !
     !  u_{k+q+G_0} carries and additional factor e^{i G_0 v}
     !
     CALL fractrasl ( nbnd, npw,  npwx, ngm, igk,  evc, eigv (:,isym), cone)
     CALL fractrasl ( nbnd, npwq, npwx, ngm, igkq, evq, eigv (:,isym), cone)
     !
     ! ---------------------------------------------------------------------
     ! wave function rotation to generate matrix elements for the star of q
     ! ---------------------------------------------------------------------
     !
     ! ps. don't use npwx instead of npw, npwq since the unused elements
     ! may be large and blow up gmapsym (personal experience)
     !
     igk (1:npw ) = gmapsym ( igk (1:npw ), isym )
     igkq(1:npwq) = gmapsym ( igkq(1:npwq), isym )
     !
     ! In dvqpsi_us_only3 we need becp1 and alphap for the rotated wfs. 
     ! The other quantities (deeq and qq) do not depend on the wfs, in
     ! particular in the KB case (not ultrasoft), the deeq's are the
     ! unscreened coefficients, and the qq's are zero.
     !
     ! For the KB part, remember dV_NL[q_0] ~ |S^-1(k)+q_0> <S^-1(k)|
     ! the total momentum transfer must be q_0 and the rotation 
     ! tranforms k+Sq_0 into S^-1(k)+q_0, k into S^-1(k)
     ! [see Eqs. (A9),(A14) Baroni et al. RMP]
     ! 
     ! Since in QE a normal rotation s is defined as S^-1 we have here
     ! sxk = S(k).  
     !
     CALL rotate_cart( xk(:,ik), s (:, :, isym), sxk)
     !
     ! here we generate vkb on the igk() set and for k ...
     CALL init_us_2 (npw, igk, sxk, vkb)
     !
     ! ... and we recompute the becp terms with the wfs (rotated through igk)
     !
     IF (noncolin) THEN
        CALL calbec (npw, vkb, evc, becp1(ik)%nc(:,:,:) )
     ELSE
        CALL calbec (npw, vkb, evc, becp1(ik)%k(:,:) )
     ENDIF
     !
     ! we also recompute the derivative of the becp terms with the (rotated) wfs
     !
     DO ipol = 1, 3
        aux2=(0.d0,0.d0)
        DO ibnd = 1, nbnd
           DO ig = 1, npw
              aux2(ig,ibnd) = evc(ig,ibnd) * tpiba * ci * ( sxk(ipol) + g(ipol,igk(ig)) )
           END DO
           IF (noncolin) THEN
              DO ig = 1, npw
                 aux2(ig+npwx,ibnd) = evc(ig+npwx,ibnd) * tpiba * ci * &
                           ( sxk(ipol) + g(ipol,igk(ig)) )
              ENDDO
           ENDIF
        ENDDO
        IF (noncolin) THEN
           CALL calbec (npw, vkb, aux2, alphap(ipol,ik)%nc(:,:,:))
        ELSE
           CALL calbec (npw, vkb, aux2, alphap(ipol,ik)%k(:,:) )
        ENDIF
     ENDDO
     !
     ! identify band edge k-point
     !
     if (.true.) then
        !
        gamma_find = .false.
        if ((abs(xk(1,ik)-xk_gamma(1))<1d-5) .and. (abs(xk(2,ik)-xk_gamma(2))<1d-5) .and. &
            (abs(xk(3,ik)-xk_gamma(3))<1d-5)) then
           gamma_find = .true.
           write(*,*) ' find gamma =', xk(:,ik), ' @ pool#', my_pool_id
        endif
        !
        x_find = .false.
        if ((abs(xk(1,ik)-xk_x(1))<1d-5) .and. (abs(xk(2,ik)-xk_x(2))<1d-5) .and. &
            (abs(xk(3,ik)-xk_x(3))<1d-5)) then
           x_find = .true.
           write(*,*) ' find X =', xk(:,ik), ' @ pool#', my_pool_id
        endif 
        !
        l_find = .false.
        if ((abs(xk(1,ik)-xk_l(1))<1d-5) .and. (abs(xk(2,ik)-xk_l(2))<1d-5) .and. &
            (abs(xk(3,ik)-xk_l(3))<1d-5)) then
           l_find = .true.
           write(*,*) ' find L =', xk(:,ik), ' @ pool#', my_pool_id
        endif 
        !
     endif
     !
     ! wave-function check (JZ)
     !
     if (.false.) then
        !
        ! evc( npwx*npol, nbnd )
        evc_check = evc
        !
        ibnd = 9
        CALL cft_wave (evc_check(:, ibnd), aux1, +1)
        !
!        write(stdout,*) ' inside elimp:, dffts%nnr, npol =',dffts%nnr, npol
!        write(stdout,*) ' xk =', xk(:,ik)

        !
!        evc_checksum = 0.d0
!        do ir = 1, dffts%nnr
!           evc_checksum = evc_checksum + real(aux1(ir,1))**2.d0 + &
!                                         dimag(aux1(ir,1))**2.d0
!        enddo
!        write(stdout,*) ' evc_checksum = ', evc_checksum
        !
        ! output
        if (l_find .and. (.true.)) then
           open(unit = 1003, file = 'evc_check_r.dat')
           write(1003, *) xk(:,ik)
           write(1003, *) ibnd
           write(1003, *) dffts%nnr
           do ir = 1, dffts%nnr
              rho_r = real(aux1(ir,1))**2.d0 + dimag(aux1(ir,1))**2.d0
              write(1003,'(3(es15.7,1x))') real(aux1(ir,1)), dimag(aux1(ir,1)), &
                                           rho_r
           enddo
           close(1003)
!           stop
        endif
        !           
     endif
     !
     ! now we generate vkb on the igkq() set because dvpsi is needed on that set
     ! we need S(k)+q_0 in the KB projector: total momentum transfer must be q_0
     !
     xktmp = sxk + xq0
     CALL init_us_2 (npwq, igkq, xktmp, vkb)
     !
     !
     IF (lsda) current_spin = isk (ik)

     call cpu_time(t2)

     ! =================================================================================
     !   Calculation of the matrix element
     ! =================================================================================

     integration_method = 1   ! 1: real space, 2: reciprocal space (faster)

     if (integration_method == 1) then

        ! ===================================================
        ! real-space method - matrix element calculation
        ! ===================================================
        !
        ! IMPORTANT NOTE: this brings xq to first-BZ, but it is only centered
        ! around Gamma, not yet fully symmetric

        ! NOTE: the q-point below should use xq, xq0 is only the first q in the star
        xq_cart(1) = xq(1) !- NINT(xq(1))
        xq_cart(2) = xq(2) !- NINT(xq(2))
        xq_cart(3) = xq(3) !- NINT(xq(3))
        !
!ERROR
        CALL cryst_to_cart ( 1, xq_cart, bg, 1 )   ! xq_cart in [2*pi/alat]

        ! supercell / unit_cell volume multiple

        sc_un_ratio = nint(scell_vol/omega)
        if ((ik == 1) .and. (xq(1)*xq(2)*xq(3) < 1e-5)) &
           write(stdout,*) '  sc_un_ratio =', sc_un_ratio
        !
        !  calculate dvscf_q*psi_k
        !
        DO ibnd = 1, nbnd !, incr
           !
           ! this transforms evc/evq(k) into evc_r/evq_r(r) with mesh given by nffts%nr1..3
           ! NOTE: evc(k) seems to represent only the periodic part u(r), JW
           CALL cft_wave (evc(:, ibnd), aux1, +1)
           evc_r(:,:,ibnd) = aux1(:,:)

           CALL cft_wave_2 (evq(:, ibnd), aux1, +1)
           evq_r(:,:,ibnd) = aux1(:,:)
           !
           ! rotate aux1 based on aux1(r=0) to zero phase
!NOTE,  we do not perform phase rotation here, and thus
!       later on we cannot simply add long-range and short-range scatterings
!       together (?)
!        r0_phase = evc_r(1,1,ibnd) / abs(evc_r(1,1,ibnd))
!        evc_r(:,:,ibnd) = evc_r(:,:,ibnd) / r0_phase
        !
!        r0_phase = evq_r(1,1,ibnd) / abs(evq_r(1,1,ibnd))
!        evq_r(:,:,ibnd) = evq_r(:,:,ibnd) / r0_phase
           !
           ! normalize for later integration on discrete mesh
           ! 1) we find out |phi|^2 integration in a unit cell
           ! 2) we properly normalize phi to require |phi|^2 integration in a unit cell gives one
           !
           norm_evc = 0.d0
           norm_evq = 0.d0
           do imesh = 1, nfftmesh
              norm_evc = norm_evc + abs(evc_r(imp_meshmap(imesh,4),1,ibnd))**2.d0
              norm_evq = norm_evq + abs(evq_r(imp_meshmap(imesh,4),1,ibnd))**2.d0
           enddo

           evc_r(:,:,ibnd) = evc_r(:,:,ibnd) / sqrt(norm_evc/sc_un_ratio)
           evq_r(:,:,ibnd) = evq_r(:,:,ibnd) / sqrt(norm_evq/sc_un_ratio)
           !
        ENDDO
   
        call cpu_time(t3)

        if (x_find .and. .false.) then
           !
           ibnd = 10
           ir_evc = imp_meshmap(14522,4)
           write(*,*) ' check meshmap(14522):', imp_meshmap(14522,1:4)

           write(*,'(a,3(1x,f15.7))') ' check_evc_r at xk =', xk(:,ik)
           write(*,*) evc_r(ir_evc,1,ibnd)

           write(*,'(a,3(1x,f15.7))') ' check_evc_r at (xq, xq0) =', xq, xq0
           write(*,*) evq_r(ir_evc,1,ibnd)
           write(*,*) ' module: ',conjg(evq_r(ir_evc,1,ibnd))*evc_r(ir_evc,1,ibnd)
 
           write(*,*) ' check dvimp(14522):', dvimp(14522)
        endif
        !
        ! calculate eimpmat(j,i)=<psi_{k+q,j}|dV_imp*psi_{k,i}> for this perturbation
        ! should check normalization: if evc/evq_r goes to 1, this reduces to plane
        ! wave form
        !
        eimpmat = (0.d0, 0.d0)
        dvimpevc = (0.d0, 0.d0)
        evq_rmesh = (0.d0, 0.d0)
        !
        do imesh = 1, nfftmesh
           !
           ! NOTE: exp(-i*q*r) term is included
           dvimpevc(imesh,1:nbnd)  = evc_r(imp_meshmap(imesh,4),1,1:nbnd) * dvimp(imesh) * &
                                     exp(-ci*twopi*dot_product(xq_cart,rp(:,imesh)))

           ! TEST
!           dvimpevc(imesh,1:nbnd)  = evc_r(imp_meshmap(imesh,4),1,1:nbnd)
!                                     exp(-ci*twopi*dot_product(xq_cart,rp(:,imesh)))

           evq_rmesh(imesh,1:nbnd) = evq_r(imp_meshmap(imesh,4),1,1:nbnd)
        enddo
        !

        do jbnd = 1, nbnd
           do ibnd = 1, nbnd
               eimpmat (jbnd, ibnd) = eimpmat (jbnd, ibnd) + &
                        ZDOTC (nfftmesh, evq_rmesh(1:nfftmesh, jbnd), 1, dvimpevc(1:nfftmesh, ibnd), 1)
           enddo
        enddo
        !
     else

        ! ===================================================
        ! reciprocal-space method - matrix element calculation
        ! ===================================================

        call cpu_time(t4_1)

        ! normalize wavefunctions

        do ibnd = 1, nbnd
           CALL cft_wave (evc(:, ibnd), aux1, +1)
           evc_r(:,:,ibnd) = aux1(:,:)

           CALL cft_wave_2 (evq(:, ibnd), aux1, +1)
           evq_r(:,:,ibnd) = aux1(:,:)
           norm_evq = 0.d0
           do imesh = 1, dffts%nnr
              norm_evq = norm_evq + abs(evq_r(imesh,1,ibnd))**2.d0
           enddo
           evq_r(:,:,ibnd) = evq_r(:,:,ibnd)/sqrt(norm_evq)
        enddo

        call cpu_time(t4_2)

        ! set up mapping between local and global G vector indices

        index_igk = 0
        do ig = 1, npw
           index_igk(igk(ig)) = ig
        enddo
        !
        ind_g = 0
        map_kq2k = 0
        do igz = -bg_max, bg_max
           do igy = -bg_max, bg_max
              do igx = -bg_max, bg_max
                 ind_g = ind_g + 1
                 ig_shift = igz*dffts%nr2*dffts%nr1 + igy*dffts%nr1 + igx

                 do ig = 1, npwq
                    if (((map_l2g(igkq(ig))-ig_shift) > 0) .and. ((map_l2g(igkq(ig))-ig_shift) < dffts%nnr)) then
                       if (map_g2l(map_l2g(igkq(ig))-ig_shift) /= 0) then
                          map_kq2k(ig,ind_g) = index_igk(map_g2l(map_l2g(igkq(ig))-ig_shift))
                       endif
                    endif
                 enddo
              enddo
           enddo
        enddo
        !
        call cpu_time(t4_3)

        ! calculate products of electron wave functions in reciprocal space

        eemat_r = 0
        do ibnd = 1, nbnd
           !
           ind_g = 0
           do igz = -bg_max, bg_max
              do igy = -bg_max, bg_max
                 do igx = -bg_max, bg_max
                    !
                    ind_g = ind_g + 1
                    ig_shift = igz*dffts%nr2*dffts%nr1 + igy*dffts%nr1 + igx
                    !
                    evc_remap = 0.d0
                    do ig = 1, npwq
                       if (map_kq2k(ig,ind_g) /= 0) &
                          evc_remap(ig) = evc(map_kq2k(ig,ind_g),ibnd)
                    enddo
                    !
                    do jbnd = 1, nbnd
                       eemat_r(jbnd,ibnd,ind_g) = &
                              ZDOTC (npwq, evq(1,jbnd), 1, evc_remap(1), 1)
                    enddo
                       ! test
!                       if ((ibnd==1) .and. (jbnd==1) .and. (ig==1)) then
!                          write(stdout,*) ' --- G test ------'
!                          write(stdout,*) ' G:', igx*bg(:,1)+igy*bg(:,2)+igz*bg(:,3)
!                          write(stdout,*) ' G for k+q:', g(:,igkq(ig))
!                          write(stdout,*) ' G for k:', g(:,map_g2l(map_l2g(igkq(ig))-ig_shift))
!                       endif
                    !
                 enddo
              enddo
           enddo
        enddo

        ! test impurity matrix element calculation
        !
        eimpmat_r = (0.d0, 0.d0)
        do jbnd = 1, nbnd
           do ibnd = 1, nbnd
              ind_g = 0
              do igz = -bg_max, bg_max
                 do igy = -bg_max, bg_max
                    do igx = -bg_max, bg_max
                       ind_g = ind_g + 1
                       eimpmat_r (jbnd,ibnd) = eimpmat_r (jbnd,ibnd) + &
                               dvimp_q(iq,igx+bg_max+1,igy+bg_max+1,igz+bg_max+1) * eemat_r(jbnd,ibnd,ind_g)
                    enddo
                 enddo
              enddo
           enddo
        enddo
        eimpmat = eimpmat_r
        !
        call cpu_time(t4_4)
        !
     endif

     !
     ! the following uses traditional integration over r-space, which has been
     ! checked to give same result as the reciprocal space integration above
     !
     if (.false.) then
        !
        do ibnd = 1, nbnd
   
           norm_evc = 0.d0
           do imesh = 1, dffts%nnr
              norm_evc = norm_evc + abs(evc_r(imesh,1,ibnd))**2.d0
           enddo
   
           ind_g = 0
           do igz = 1, 2*bg_max+1
              do igy = 1, 2*bg_max+1
                 do igx = 1, 2*bg_max+1
                
                    ind_g = ind_g + 1
   !                 write(stdout,*) ind_g
   
                    aux2 = 0.d0
                    aux1(:,:) = evc_r(:,:,ibnd)
                    CALL apply_dpot(dffts%nnr, aux1, dvGr(:,igx,igy,igz), current_spin)
                    dvevc = aux1/sqrt(norm_evc)
                    CALL cft_wave (aux2(:, ibnd), aux1, -1)
   
                    do jbnd = 1, nbnd
                       eemat (jbnd, ibnd, ind_g) = &
                          ZDOTC (npwq, evq(1, jbnd), 1, aux2(1, ibnd), 1)
!                       eemat_check (jbnd, ibnd, ind_g) = &
!                          ZDOTC (dffts%nnr, evq_r(:,1,jbnd), 1, dvevc(:,1), 1)
                    enddo
              
                 enddo
              enddo
           enddo
        enddo
        !
     endif


     call cpu_time(t4)
     !
     if (x_find .and. .false.) then
        !
        ibnd = 10
        write(*,'(a,3(1x,f15.7))') ' check_eimpmat before mp_sum at xq =', xq
        write(*,*) '       ibnd  = ',ibnd
        write(*,*) '     eimpmat = ',eimpmat(ibnd,ibnd)
        !
     endif
     !
#ifdef __PARA
     CALL mp_sum(eimpmat, intra_pool_comm)
     CALL mp_sum(eimpmat_r, intra_pool_comm)
     CALL mp_sum(eemat, intra_pool_comm)
     CALL mp_sum(eemat_check, intra_pool_comm)
     CALL mp_sum(eemat_r, intra_pool_comm)
#endif
     !
     if (x_find .and. .true.) then
        !
        ibnd = 10
        write(*,'(a,3(1x,f15.7))') ' check_eimpmat after mp_sum at xq =', xq
        write(*,*) '       ibnd  = ',ibnd
        write(*,*) '    eimpmat = ',eimpmat(9:10,9:10)
        write(*,*) '    eimpmat_r = ',eimpmat_r(9:10,9:10)

!        write(*,*) '    eemat =', eemat(ibnd,ibnd,12:16)
!        write(*,*) '    eemat_check =', eemat_check(ibnd,ibnd,12:16)
!        write(*,*) '    eemat_r =', eemat_r(ibnd,ibnd,12:16)
        !
!        write(*,*) '     check matrix eimpmat(8..12,8..12)'
!        write(*,*) '     eimpmat(8, 8..12)', eimpmat(8, 8:12)
!        write(*,*) '     eimpmat(9, 8..12)', eimpmat(9, 8:12)
!        write(*,*) '     eimpmat(10,8..12)', eimpmat(10,8:12)
!        write(*,*) '     eimpmat(11,8..12)', eimpmat(11,8:12)
!        write(*,*) '     eimpmat(12,8..12)', eimpmat(12,8:12)
        !
     endif
     !
     !  Rotate eimpmat with the gauge matrices (this should be equivalent 
     !  to calculate eimpmat with the truely rotated eigenstates)
     ! 
     ! the two zgemm call perform the following ops:
     !  eimpmat = umat(k+q)^\dagger * [ eimpmat * umat(k) ]
     !
     CALL zgemm ('n', 'n', nbnd, nbnd, nbnd, cone, eimpmat(:,:), &
          nbnd, umat(:,:,ik), nbnd, czero, eptmp, nbnd)

     ! the second transforming matrix that we really need is umat at k+q, which
     ! is actually just [umatq]
     CALL zgemm ('c', 'n', nbnd, nbnd, nbnd, cone, umatq(:,:,ik), &
          nbnd, eptmp, nbnd, czero, eimpmat(:,:), nbnd)

     ! =================================================================================
     !
     !
     !  save eph matrix elements into el_ph_mat
     !
     DO jbnd = 1, nbnd
        DO ibnd = 1, nbnd
           el_imp_mat (ibnd, jbnd, ik) = eimpmat (ibnd, jbnd)
        ENDDO
     ENDDO
     !
     call cpu_time(t5)
     !
     write(stdout,*) ' t4_1 -> t4_2     :', t4_2-t4_1, 's'
     write(stdout,*) ' t4_2 -> t4_3     :', t4_3-t4_2, 's'
     write(stdout,*) ' t4_3 -> t4_4     :', t4_4-t4_3, 's'
     write(stdout,*) ' t4_4 -> t4     :', t4-t4_4, 's'

     !
     tt1 = tt1 + (t2-t1)
     tt2 = tt2 + (t3-t2)
     tt3 = tt3 + (t4-t3)
     tt4 = tt4 + (t5-t4)
     !
  ENDDO
  !

  write(stdout,*) ' ene     :', tt1, 's'
  write(stdout,*) ' evc     :', tt2, 's'
  write(stdout,*) ' eimpmat :', tt3, 's'
  write(stdout,*) ' rotate  :', tt4, 's'

  !
  !  restore original configuration of files
  !
  CALL seqopn (iunigk, 'igk', 'unformatted', exst)
  CALL diropn (iuwfc, 'wfc', lrwfc, exst)
#ifdef __PARA
  ! never remove this barrier - > insures that wfcs are restored to each pool
  ! before moving on
  CALL mp_barrier(inter_pool_comm)
#endif
  !
  DEALLOCATE (eimpmat, eptmp, aux1, aux2)
  DEALLOCATE (gmap, shift)
  deallocate (evc_r, evq_r)
  deallocate (dvimpevc, evq_rmesh)
  !
  END SUBROUTINE elimp_shuffle






