!Modules and subroutines for implementation of tetrahedra method



module tetrahedron
USE kinds, ONLY : DP
implicit none
!type point
!real(DP) :: w,c
!w: eigenvalues at this k-point
!c: weights for imaginary part (selfenergy)
!integer :: i  !index of the k-point
!end type point

!type tetratype
!type(point) :: p(4)
!end type tetratype

integer, allocatable :: itetra(:,:,:), tetra_dose(:,:)
!(2,30,nkps) store which tetrahedron a given k-point belongs to, and the index
!inside that tetrahedron
integer :: ntetra, ntetra_dose
!total number of tetrahedron
!type(tetratype), allocatable :: tetra(:)
!array storing information of all tetrahedron
real(DP), allocatable :: wkt(:)

integer, allocatable :: tetra_i(:,:)
!k index of each vertex of the tetrahedron
!(4,ntetra)
real(DP), allocatable :: tetra_w(:,:)
!eigenvalues associated with each vertex
!(4,ntetra)
real(DP), allocatable :: tetra_c(:,:)
!weights for imaginary part assocates with vertex
!(4,ntetra)

end module tetrahedron




subroutine make_kp_reg(nx,ny,nz,wkt,ntet, tet_i, itet)
! generate a regular mesh from 0 to g_i, with eventual shift of origin only
! Tetrahedron are generated. How are Tetra near BZ boundaries treated correctly?
USE kinds, ONLY : DP
implicit none
integer :: i,j,k,l,nx,ny,nz,nk,nt,np,pt(8), ntet
integer :: itet(2,30,nx*ny*nz), nktetra(nx*ny*nz)
!itet stores the tetrahedron that the k point belongs to
!nktetra stores how many tetrahedron the k point belongs to
real(DP) :: q(3),ep(3),kpt(3,nx*ny*nz),wkt(nx*ny*nz)
integer :: tet_i(4,ntet)

ep = 0d0  ! shift by a small amount so that only one boundary K remains in the FBZ after folding
! shft = (-0.5d0)*(g1+g2+g3)
itet=0 !0 means the k-point do not belong to a tetrahedra 
nktetra=0

do i=1,nx*ny*nz
        wkt(i)=1d0/(nx*ny*nz)
enddo

np = 0
nt = 0
do i = 1,nx
        do j = 1,ny
                do k = 1,nz                     
                        do l = 1,6
                                pt(1) = (i-1)*ny*nz + (j-1)*nz + (k-1) +1
                                pt(2) = pt(1)+1                                 
                                pt(3) = pt(1)+nz
                                pt(4) = pt(1)+nz+1
                                pt(5) = pt(1)+nz*ny
                                pt(6) = pt(1)+nz*ny+1
                                pt(7) = pt(1)+nz*ny+nz
                                pt(8) = pt(1)+nz*ny+nz+1
                                if(k .eq. nz) then
                                        pt(2)=pt(2)-nz
                                        pt(4)=pt(4)-nz
                                        pt(6)=pt(6)-nz
                                        pt(8)=pt(8)-nz
                                endif
                                if(j .eq. ny) then
                                        pt(3)=pt(3)-ny*nz
                                        pt(4)=pt(4)-ny*nz
                                        pt(7)=pt(7)-ny*nz
                                        pt(8)=pt(8)-ny*nz
                                endif
                                if(i .eq. nx) then
                                        pt(5)=pt(5)-nx*ny*nz
                                        pt(6)=pt(6)-nx*ny*nz
                                        pt(7)=pt(7)-nx*ny*nz
                                        pt(8)=pt(8)-nx*ny*nz
                                endif

                                nt = nt+1 !number of tetrahedron
                                np = np+1 !number of points?

                                if (l.eq.1) then                        ! 1,2,3,6
                                        tet_i(1,nt) = pt(1) !indices
                                        tet_i(2,nt) = pt(2)
                                        tet_i(3,nt) = pt(3)
                                        tet_i(4,nt) = pt(6)
                                        nktetra(pt(1)) = nktetra(pt(1))+1
                                        itet(1,nktetra(pt(1)),pt(1))=nt
                                        itet(2,nktetra(pt(1)),pt(1))=1
                                        nktetra(pt(2)) = nktetra(pt(2))+1
                                        itet(1,nktetra(pt(2)),pt(2))=nt
                                        itet(2,nktetra(pt(2)),pt(2))=2
                                        nktetra(pt(3)) = nktetra(pt(3))+1
                                        itet(1,nktetra(pt(3)),pt(3))=nt
                                        itet(2,nktetra(pt(3)),pt(3))=3
                                        nktetra(pt(6)) = nktetra(pt(6))+1
                                        itet(1,nktetra(pt(6)),pt(6))=nt
                                        itet(2,nktetra(pt(6)),pt(6))=4
                                elseif (l.eq.2) then                    ! 2,3,4,6
                                        tet_i(1,nt) = pt(2) !indices
                                        tet_i(2,nt) = pt(3)
                                        tet_i(3,nt) = pt(4)
                                        tet_i(4,nt) = pt(6)
                                        nktetra(pt(4)) = nktetra(pt(4))+1
                                        itet(1,nktetra(pt(4)),pt(4))=nt
                                        itet(2,nktetra(pt(4)),pt(4))=3
                                        nktetra(pt(2)) = nktetra(pt(2))+1
                                        itet(1,nktetra(pt(2)),pt(2))=nt
                                        itet(2,nktetra(pt(2)),pt(2))=1
                                        nktetra(pt(3)) = nktetra(pt(3))+1
                                        itet(1,nktetra(pt(3)),pt(3))=nt
                                        itet(2,nktetra(pt(3)),pt(3))=2
                                        nktetra(pt(6)) = nktetra(pt(6))+1
                                        itet(1,nktetra(pt(6)),pt(6))=nt
                                        itet(2,nktetra(pt(6)),pt(6))=4
                                elseif (l.eq.3) then                    ! 1,3,5,6
                                        tet_i(1,nt) = pt(1) !indices
                                        tet_i(2,nt) = pt(3)
                                        tet_i(3,nt) = pt(5)
                                        tet_i(4,nt) = pt(6)
                                        nktetra(pt(1)) = nktetra(pt(1))+1
                                        itet(1,nktetra(pt(1)),pt(1))=nt
                                        itet(2,nktetra(pt(1)),pt(1))=1
                                        nktetra(pt(5)) = nktetra(pt(5))+1
                                        itet(1,nktetra(pt(5)),pt(5))=nt
                                        itet(2,nktetra(pt(5)),pt(5))=3
                                        nktetra(pt(3)) = nktetra(pt(3))+1
                                        itet(1,nktetra(pt(3)),pt(3))=nt
                                        itet(2,nktetra(pt(3)),pt(3))=2
                                        nktetra(pt(6)) = nktetra(pt(6))+1
                                        itet(1,nktetra(pt(6)),pt(6))=nt
                                        itet(2,nktetra(pt(6)),pt(6))=4
                                elseif (l.eq.4) then                    ! 3,4,6,8
                                        tet_i(1,nt) = pt(3) !indices
                                        tet_i(2,nt) = pt(4)
                                        tet_i(3,nt) = pt(6)
                                        tet_i(4,nt) = pt(8)
                                        nktetra(pt(4)) = nktetra(pt(4))+1
                                        itet(1,nktetra(pt(4)),pt(4))=nt
                                        itet(2,nktetra(pt(4)),pt(4))=2
                                        nktetra(pt(8)) = nktetra(pt(8))+1
                                        itet(1,nktetra(pt(8)),pt(8))=nt
                                        itet(2,nktetra(pt(8)),pt(8))=4
                                        nktetra(pt(3)) = nktetra(pt(3))+1
                                        itet(1,nktetra(pt(3)),pt(3))=nt
                                        itet(2,nktetra(pt(3)),pt(3))=1
                                        nktetra(pt(6)) = nktetra(pt(6))+1
                                        itet(1,nktetra(pt(6)),pt(6))=nt
                                        itet(2,nktetra(pt(6)),pt(6))=3
                                elseif (l.eq.5) then                    ! 3,5,6,7
                                        tet_i(1,nt) = pt(3) !indices
                                        tet_i(2,nt) = pt(5)
                                        tet_i(3,nt) = pt(6)
                                        tet_i(4,nt) = pt(7)
                                        nktetra(pt(5)) = nktetra(pt(5))+1
                                        itet(1,nktetra(pt(5)),pt(5))=nt
                                        itet(2,nktetra(pt(5)),pt(5))=2
                                        nktetra(pt(7)) = nktetra(pt(7))+1
                                        itet(1,nktetra(pt(7)),pt(7))=nt
                                        itet(2,nktetra(pt(7)),pt(7))=4
                                        nktetra(pt(3)) = nktetra(pt(3))+1
                                        itet(1,nktetra(pt(3)),pt(3))=nt
                                        itet(2,nktetra(pt(3)),pt(3))=1
                                        nktetra(pt(6)) = nktetra(pt(6))+1
                                        itet(1,nktetra(pt(6)),pt(6))=nt
                                        itet(2,nktetra(pt(6)),pt(6))=3
                                else                            ! 3,6,7,8
                                        tet_i(1,nt) = pt(3) !indices
                                        tet_i(2,nt) = pt(6)
                                        tet_i(3,nt) = pt(7)
                                        tet_i(4,nt) = pt(8)
                                        nktetra(pt(7)) = nktetra(pt(7))+1
                                        itet(1,nktetra(pt(7)),pt(7))=nt
                                        itet(2,nktetra(pt(7)),pt(7))=3
                                        nktetra(pt(8)) = nktetra(pt(8))+1
                                        itet(1,nktetra(pt(8)),pt(8))=nt
                                        itet(2,nktetra(pt(8)),pt(8))=4
                                        nktetra(pt(3)) = nktetra(pt(3))+1
                                        itet(1,nktetra(pt(3)),pt(3))=nt
                                        itet(2,nktetra(pt(3)),pt(3))=1
                                        nktetra(pt(6)) = nktetra(pt(6))+1
                                        itet(1,nktetra(pt(6)),pt(6))=nt
                                        itet(2,nktetra(pt(6)),pt(6))=2
                                endif
                        enddo   
                enddo
        enddo
enddo
!DO i=1,nx*ny*nz
!WRITE(6,'(a,i5)') 'nktetra=', nktetra(i)
!ENDDO
! write(ulog,*)'KP_REG: Number of regular kpoints generated is=',nk
! 
! if (mod(nx,2).eq.0 .and. mod(ny,2).eq.0 .and. mod(nz,2).eq.0 ) then
!       do i=1,nk
!               kpt(:,i) = kpt(:,i)+shft(:)
!       enddo
! endif
! 
! 2  format(i7,2x,3(1x,f12.5),5x,f9.5)
! 3  format(3(i3),2x,i6,2x,3(1x,f12.5),5x,f9.5)
! close(126)

end subroutine make_kp_reg


!this subroutine is modified by BL to include only single-band
subroutine eigen_tet(ntet,eival,tet_i,tet_w,nkp)
USE kinds, ONLY : DP
implicit none
integer j,k,ntet, nkp,ww,cnd
real(DP), intent(in) :: eival(nkp)
integer, intent(in) :: tet_i(4,ntet)
real(DP), intent(out) :: tet_w(4,ntet)

do j=1,ntet !! this many tetra
        do k=1,4   !! corners
                        ww = tet_i(k,j)  !! the index of k points in the original kmesh
                        tet_w(k,j)=eival(ww) !! eigenvalues (band energies of band l, kpoint "w"
        enddo
enddo

end subroutine eigen_tet

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine weight_tet(nktot, ntet,e,tet_i,tet_w,tet_c,wkt)
!! input:  number of bands, the energy e that the weights are going to be computed, tet: tetrahedron
!! the first 5 parameters are easy. Just need to get tet right.
USE kinds, ONLY : DP
implicit none
integer i,k,l,m,nktot, ntet,uns(4),zero,kkk,ii
real(DP) :: e1,e2,e3,e4,e21,e31,e41,e32,e42,e43,en(4),en1(4),a1,a2,a3,a4,b11,b12,b21,b22,b31,b32,b41,b42
real(DP) :: a11,a12,a13,a14,a21,a22,a23,a24,a31,a32,a33,a34,a41,a42,a43,a44
real(DP) :: d11,d12,d13,d21,d22,d23,d31,d32,d33,d41,d42,d43
real(DP) :: c11,c12,c21,c22,c31,c32,c41,c42
real(DP) :: s1,s2,s3,s4,s5,s6
real(DP) :: e,wkt(nktot)
integer, intent(in) :: tet_i(4,ntet)
real(DP), intent(in) :: tet_w(4,ntet)
real(DP), intent(out) :: tet_c(4,ntet)

zero=0

do i=1,ntet    ! allocate weights for imaginary terms
       ! do j=1,nb  !number of bands
                do k=1,4
                        en(k)=tet_w(k,i) !eigenvalues
                enddo
                en1=en !pre-sorting
                call ssort(en,uns,4)  !sorting
                if (en(1).ne.en(2) .and. en(2).ne.en(3) .and. en(3).ne.en(4)) then
                        do k=1,4
                                if (en(k).eq.en1(1)) then
                                        uns(k) = 1
                                elseif (en(k).eq.en1(2)) then
                                        uns(k) = 2
                                elseif (en(k).eq.en1(3)) then
                                        uns(k) = 3
                                elseif (en(k).eq.en1(4)) then
                                        uns(k) = 4
                                endif
                        enddo
                else
                        do k=1,4
                                if (en(k).eq.en1(1)) then
                                        uns(k) = 1
                                elseif (en(k).eq.en1(2)) then
                                        uns(k) = 2
                                elseif (en(k).eq.en1(3)) then
                                        uns(k) = 3
                                elseif (en(k).eq.en1(4)) then
                                        uns(k) = 4
                                endif
                                if (k.ne.1) then
                                        do l=1,k-1
                                                if (uns(l).eq.uns(k) .and. uns(k).eq.1) then
                                                        if (en(k).eq.en1(2)) then
                                                                uns(k) = 2
                                                        elseif (en(k).eq.en1(3)) then
                                                                uns(k) = 3
                                                        elseif (en(k).eq.en1(4)) then
                                                                uns(k) = 4
                                                        endif
                                                elseif (uns(l).eq.uns(k) .and. uns(k).eq.2) then
                                                        if (en(k).eq.en1(3)) then
                                                                uns(k) = 3
                                                        elseif (en(k).eq.en1(4)) then
                                                                uns(k) = 4
                                                        endif
                                                elseif (uns(l).eq.uns(k) .and. uns(k).eq.3) then
                                                        uns(k) = 4
                                                endif
                                        enddo
                                endif
                        enddo  !! k = 1,4
                endif
                e1 = en(4)
                e2 = en(3)
                e3 = en(2)
                e4 = en(1)
                e21 = e2-e1
                e31 = e3-e1
                e41 = e4-e1
                e32 = e3-e2
                e42 = e4-e2
                e43 = e4-e3
!               if(e21 .eq. 0. .or. e31 .eq. 0. .or. e41 .eq. 0. .or. e32 .eq. 0. .or. e42 .eq. 0. .or. e43 .eq. 0.) then 
!                       write(*,*) "Equal e tetra corners: ", i
!               endif
                if (e.lt.e1 .or. e.gt.e4) then
                        tet_c(uns(4),i) = 0
                        tet_c(uns(3),i) = 0
                        tet_c(uns(2),i) = 0
                        tet_c(uns(1),i) = 0
                elseif (e1.le.e .and. e.le.e2) then
                        tet_c(uns(4),i) = ((e2-e)/e21+(e3-e)/e31+(e4-e)/e41)*(e-e1)**2/(e41*e31*e21)
                        tet_c(uns(3),i) = (e-e1)**3/(e21**2*e31*e41)
                        tet_c(uns(2),i) = (e-e1)**3/(e21*e31**2*e41)
                        tet_c(uns(1),i) = (e-e1)**3/(e21*e31*e41**2)
                elseif (e2.le.e .and. e.le.e3) then
                        c11 = (e3-e)/e31**2
                        c12 = (e4-e)/e41**2
                        c21 = (e3-e)/e32**2
                        c22 = (e4-e)/e42**2
                        c31 = (e-e2)/e32**2
                        c32 = (e-e1)/e31**2
                        c41 = (e-e2)/e42**2
                        c42 = (e-e1)/e41**2
                        b11 = (e3-e)*(e-e2)/(e42*e32)+(e4-e)*(e-e1)/(e41*e42)+(e3-e)*(e-e1)/(e32*e41)
                        b12 = (e4-e)*(e-e1)/(e42*e31)+(e4-e)*(e-e2)/(e42*e32)+(e3-e)*(e-e1)/(e31*e32)
                        b21 = (e3-e)*(e-e2)/(e42*e31)+(e4-e)*(e-e2)/(e42*e41)+(e3-e)*(e-e1)/(e31*e41)
                        b22 = (e3-e)*(e-e2)/(e32*e31)+(e4-e)*(e-e1)/(e41*e31)+(e4-e)*(e-e2)/(e32*e41)
                        b31 = (e3-e)*(e-e2)/(e42*e31)+(e4-e)*(e-e2)/(e42*e41)+(e3-e)*(e-e1)/(e31*e41)
                        b32 = (e3-e)*(e-e2)/(e42*e32)+(e4-e)*(e-e1)/(e41*e42)+(e3-e)*(e-e1)/(e32*e41)
                        b41 = (e3-e)*(e-e2)/(e32*e31)+(e4-e)*(e-e1)/(e41*e31)+(e4-e)*(e-e2)/(e32*e41)
                        b42 = (e4-e)*(e-e1)/(e42*e31)+(e4-e)*(e-e2)/(e42*e32)+(e3-e)*(e-e1)/(e31*e32)
                        tet_c(uns(4),i) = .5*(c11*b11+c12*b12)
                        tet_c(uns(3),i) = .5*(c21*b21+c22*b22)
                        tet_c(uns(2),i) = .5*(c31*b31+c32*b32)
                        tet_c(uns(1),i) = .5*(c41*b41+c42*b42)
                elseif (e3.le.e .and. e.le.e4) then
                        tet_c(uns(4),i) = (e4-e)**3/(e41**2*e42*e43)
                        tet_c(uns(3),i) = (e4-e)**3/(e41*e42**2*e43)
                        tet_c(uns(2),i) = (e4-e)**3/(e41*e42*e43**2)
                        tet_c(uns(1),i) = ((e-e3)/e43+(e-e2)/e42+(e-e1)/e41)*(e4-e)**2/(e41*e42*e43)
                endif
                tet_c(uns(1),i)=tet_c(uns(1),i)*wkt(tet_i(uns(1),i))*1.0d0/6.0d0 !! the weight is V_t/V_G which is the inverse of # tetra
                tet_c(uns(2),i)=tet_c(uns(2),i)*wkt(tet_i(uns(2),i))*1.0d0/6.0d0 !! spin=2 should not be included here, it's not relevant
                tet_c(uns(3),i)=tet_c(uns(3),i)*wkt(tet_i(uns(3),i))*1.0d0/6.0d0
                tet_c(uns(4),i)=tet_c(uns(4),i)*wkt(tet_i(uns(4),i))*1.0d0/6.0d0
                do ii=1,4
                        if(isnan(tet_c(uns(ii),i) ) ) then
                                tet_c(uns(ii),i)=0.0d0
                        endif
                enddo
       !enddo
enddo

end subroutine weight_tet

SUBROUTINE SSORT (X, IY, N)
USE kinds, ONLY : DP
IMPLICIT NONE

INTEGER N
REAL(DP) :: X(1:N)
INTEGER IY(N)
REAL(DP) :: TEMP
INTEGER I, ISWAP(1), ITEMP, ISWAP1
INTRINSIC MAXLOC
DO 200 I=1,N-1
        ISWAP=MAXLOC(X(I:N))
        ISWAP1=ISWAP(1)+I-1
        IF(ISWAP1.NE.I) THEN
                TEMP=X(I)
                X(I)=X(ISWAP1)
                X(ISWAP1)=TEMP
                ITEMP=IY(I)
                IY(I)=IY(ISWAP1)
                IY(ISWAP1)=ITEMP
        ENDIF
200 CONTINUE
RETURN
END



!broadcast an array of tetratype - BL
!!!!!    SUBROUTINE mp_bcast_tetra(msg,source,gid)
!!!!!        USE tetrahedron, ONLY : point, tetratype
!!!!!        USE kinds,       ONLY : DP
!!!!!        IMPLICIT NONE
!!!!!        !
!!!!!        TYPE(tetratype) :: msg(:)
!!!!!        INTEGER :: source
!!!!!        INTEGER :: gid
!!!!!        INTEGER :: group
!!!!!        INTEGER :: msglen
!!!!!        INTEGER :: im
!!!!!#if defined(__MPI)
!!!!!        msglen = size(msg)
!!!!!        group = gid
!!!!!        DO im = 1,msglen
!!!!!            CALL bcast_point(msg(im),source,group)
!!!!!        ENDDO
!!!!!       ! mp_call_count( 18 ) = mp_call_count( 18 ) + 1
!!!!!       ! mp_call_sizex( 18 ) = MAX( mp_call_sizex( 18 ), msglen )
!!!!!#endif
!!!!!     END SUBROUTINE mp_bcast_tetra
!!!!!
!!!!!!broadcast an array of point - BL
!!!!!      SUBROUTINE bcast_point(msg,source,gid)
!!!!!         USE tetrahedron, ONLY : point
!!!!!         USE kinds,       ONLY : DP
!!!!!         IMPLICIT NONE
!!!!!#if defined __MPI
!!!!!         INCLUDE 'mpif.h'
!!!!!#endif
!!!!!         TYPE(point) :: msg(:)
!!!!!         INTEGER :: source
!!!!!         INTEGER :: gid
!!!!!         INTEGER :: group
!!!!!         INTEGER :: msglen
!!!!!         INTEGER :: ip
!!!!!         
!!!!!#if defined(__MPI)
!!!!!         msglen = size(msg) 
!!!!!         group  = gid
!!!!!         DO ip = 1,msglen
!!!!!             CALL bcast_real(msg(ip)%w,1,source,group)
!!!!!             CALL bcast_real(msg(ip)%c,1,source,group)
!!!!!             CALL bcast_integer(msg(ip)%i,1,source,group)
!!!!!         ENDDO
!!!!!#endif
!!!!!       END SUBROUTINE bcast_point
!!!!!
!!!!!         
!!!!!      
