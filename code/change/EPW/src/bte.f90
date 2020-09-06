!--------------------------------------------------------------------------
MODULE bte_var
!-------------------------------------------------------------------------- 
  !
  USE kinds, ONLY : DP
  !
  SAVE
  !
  !=========== k-mesh ===========
  REAL(KIND=DP), ALLOCATABLE :: &
          ! bz_index
          xkf_ful(:,:)         ,&
          xkf_ful_red(:,:)     ,&
          xkf_irr(:,:)         ,&
          xkf_irr_red(:,:)     ,&
          xkf_fbz(:,:)         ,&
          wkf_irr(:)           ,&
          symmat_kpt(:,:,:)    ,&
          xk_kmesh(:,:)        ,&
          sigmai_ch(:,:,:,:,:)
  !
  INTEGER, ALLOCATABLE       :: &
          ! bz_index
          ful2irr(:)           ,& ! ful-BZ => irr-BZ
          irr2ful(:)           ,& ! irr-BZ => ful-BZ
          irr2ful_sym(:)       ,&
          ignore_kpt(:)        ,&
          ! cpu_index
          seq2nscat(:)         ,&
          nscat_new(:)         ,&
          ! ephwann_shuffle
          rirr2irr(:)          ,& ! red-irr-BZ =>     irr-BZ
          rful2ful(:)          ,& ! red-ful-BZ =>     ful-BZ
          ful2rful(:)          ,& !     ful-BZ => red-ful-BZ
          rful2rirr(:)         ,& ! red-ful-BZ => red-irr-BZ
          nscat_all(:)         ,&
          band_ch(:)  
  !
  INTEGER                    :: &
          nk_ful               ,&
          nk_irr               ,&
          nk_ful_red           ,&
          nk_irr_red           ,&
          nk_kmesh
  !
  !=========== q-mesh ===========
  REAL(KIND=DP), ALLOCATABLE :: &
          ! bz_index_q
          xqf_ful(:,:)         ,&
          xqf_irr(:,:)         ,&
          xqf_fbz(:,:)         ,&
          wqf_irr(:)           ,&
          symmat_qpt(:,:,:)    ,&
          xq_qmesh(:,:)
  !
  INTEGER, ALLOCATABLE       :: &
          ! bz_index_q
          ful2irr_q(:)         ,&
          irr2ful_q(:)         ,&
          irr2ful_sym_q(:)     ,&
          ignore_qpt(:)        ,&
          ! cpu_index_q
          seq2nscat_q(:)       ,&
          nscat_new_q(:)       ,&
          ! ephwann_shuffle
          rirr2irr_q(:)        ,&
          rful2ful_q(:)        ,&
          ful2rful_q(:)        ,&
          rful2rirr_q(:)       ,& ! red-ful-BZ => red-irr-BZ
          nscat_all_q(:)
  !
  INTEGER                    :: &
          nq_ful               ,&
          nq_irr               ,&
          nq_ful_red           ,&
          nq_irr_red           ,&
          nq_qmesh
  !
  !========== all-mesh ==========
  INTEGER                    :: &
          mkq1                 ,&
          mkq2                 ,&
          mkq3                 ,&
          n_rot
  !
  REAL(KIND=DP), ALLOCATABLE :: &
          ! bz_index
          symmat_lat(:,:,:)
  !
  !========== all-mesh ==========
  INTEGER                    :: &
          nqf1_phd             ,&
          nqf2_phd             ,&
          nqf3_phd             ,&
          nq_ful_phd           ,&
          nmodes_phd
  REAL(KIND=DP)              :: &
          eptemp_phd
  !
  REAL(KIND=DP), ALLOCATABLE :: &
           N_q_ful(:,:,:,:)      ,&
       int_N_q_ful(:,:,:,:)      ,&
          dN_q_ful(:,:,:,:)      ,&
          ph_rate_ful(:,:,:)     ,&
          alelrate(:,:)         ,&  ! read in electron (eV) and average alloy-electron scatteting rate (THz)    
          tau0_q_ful(:,:,:)
  !
!-------------------------------------------------------------------------- 
END MODULE bte_var
!-------------------------------------------------------------------------- 



!----------------------------------------------------------------------------
MODULE bte_func
!----------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  CONTAINS
  !
  !----------------------------------------------------------------------------
  FUNCTION equiv_k (x, y)
  !----------------------------------------------------------------------------
    USE kinds,          ONLY : DP
    !
    ! x, y should br in crystal coordinate
    IMPLICIT NONE
    REAL(KIND=DP), INTENT(IN) :: x(3), y(3)
    LOGICAL                   :: equiv_k
    !
    equiv_k = ( ABS(x(1)-y(1) - NINT(x(1)-y(1))) .LT. 1.0d-5 ) .AND. &
             ( ABS(x(2)-y(2) - NINT(x(2)-y(2))) .LT. 1.0d-5 ) .AND. &
             ( ABS(x(3)-y(3) - NINT(x(3)-y(3))) .LT. 1.0d-5 )
    !
  END FUNCTION equiv_k

  !----------------------------------------------------------------------------
  FUNCTION ijk_fbz (ijk_k, ijk_q)
  !----------------------------------------------------------------------------
    USE epwcom, ONLY : nkf1, nkf2, nkf3
    USE bte_var,  ONLY : mkq1, mkq2, mkq3
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: ijk_k(3), ijk_q(3)
    INTEGER, DIMENSION(3) :: ijk_fbz
    !
    ijk_fbz(1) = MOD( ( ijk_k(1) + mkq1*ijk_q(1) ),nkf1 )
    ijk_fbz(2) = MOD( ( ijk_k(2) + mkq2*ijk_q(2) ),nkf2 )
    ijk_fbz(3) = MOD( ( ijk_k(3) + mkq3*ijk_q(3) ),nkf3 )
    !
    IF (ijk_fbz(1) .LT. 0) ijk_fbz(1) = nkf1 + ijk_fbz(1)
    IF (ijk_fbz(2) .LT. 0) ijk_fbz(2) = nkf2 + ijk_fbz(2)
    IF (ijk_fbz(3) .LT. 0) ijk_fbz(3) = nkf3 + ijk_fbz(3)
    !
  END FUNCTION ijk_fbz


  !----------------------------------------------------------------------------
  FUNCTION ijk2id (ijk)
  !----------------------------------------------------------------------------
    USE epwcom, ONLY : nkf1, nkf2, nkf3
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: ijk(3)
    INTEGER             :: ijk2id
    !
    ijk2id = ijk(1)*nkf2*nkf3 + ijk(2)*nkf3 + ijk(3) + 1
    !
  END FUNCTION ijk2id


  !----------------------------------------------------------------------------
  FUNCTION ijk2id_q (ijk)
  !----------------------------------------------------------------------------
    USE epwcom, ONLY : nqf1, nqf2, nqf3
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: ijk(3)
    INTEGER             :: ijk2id_q
    !
    ijk2id_q = ijk(1)*nqf2*nqf3 + ijk(2)*nqf3 + ijk(3) + 1
    !
  END FUNCTION ijk2id_q


  !----------------------------------------------------------------------------
  FUNCTION id2ijk (id)
  !----------------------------------------------------------------------------
    USE epwcom, ONLY : nkf1, nkf2, nkf3
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN)   :: id
    INTEGER, DIMENSION(3) :: id2ijk
    !
    id2ijk(3) = MOD(id-1,nkf3)
    id2ijk(2) = MOD((id-1)/nkf3,nkf2)
    id2ijk(1) = MOD((id-1)/nkf3/nkf2,nkf1)
    !
  END FUNCTION id2ijk


  !----------------------------------------------------------------------------
  FUNCTION id2ijk_q (id)
  !----------------------------------------------------------------------------
    USE epwcom, ONLY : nqf1, nqf2, nqf3
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN)   :: id
    INTEGER, DIMENSION(3) :: id2ijk_q
    !
    id2ijk_q(3) = MOD(id-1,nqf3)
    id2ijk_q(2) = MOD((id-1)/nqf3,nqf2)
    id2ijk_q(1) = MOD((id-1)/nqf3/nqf2,nqf1)
    !
  END FUNCTION id2ijk_q


  !----------------------------------------------------------------------------
  FUNCTION imat (mat)
  !----------------------------------------------------------------------------
    USE kinds,        ONLY : DP
    !
    IMPLICIT NONE
    REAL(KIND=DP), INTENT(IN)     :: mat(3,3)
    REAL(KIND=DP), DIMENSION(3,3) :: imat(3,3)
    REAL(KIND=DP)                 :: det
    !
    imat(1,1) = mat(2,2)*mat(3,3) - mat(3,2)*mat(2,3)
    imat(1,2) = mat(3,2)*mat(1,3) - mat(1,2)*mat(3,3)
    imat(1,3) = mat(1,2)*mat(2,3) - mat(2,2)*mat(1,3)
    imat(2,1) = mat(2,3)*mat(3,1) - mat(3,3)*mat(2,1)
    imat(2,2) = mat(3,3)*mat(1,1) - mat(1,3)*mat(3,1)
    imat(2,3) = mat(1,3)*mat(2,1) - mat(2,3)*mat(1,1)
    imat(3,1) = mat(2,1)*mat(3,2) - mat(3,1)*mat(2,2)
    imat(3,2) = mat(3,1)*mat(1,2) - mat(1,1)*mat(3,2)
    imat(3,3) = mat(1,1)*mat(2,2) - mat(2,1)*mat(1,2)
    !
    det = mat(1,1)*mat(2,2)*mat(3,3) + mat(1,2)*mat(2,3)*mat(3,1) + mat(1,3)*mat(3,2)*mat(2,1) - &
          mat(1,3)*mat(2,2)*mat(3,1) - mat(2,3)*mat(3,2)*mat(1,1) - mat(3,3)*mat(1,2)*mat(2,1)
    !
    IF (det .NE. 0.0d0) THEN
       imat = imat / det
    ELSE
       CALL errore ('inv_mat','matrix is singular',1)
    ENDIF
    !
  END FUNCTION imat

!----------------------------------------------------------------------------
  FUNCTION FindDet (matrix, n)
  !----------------------------------------------------------------------------
    !Function to find the determinant of a square matrix
    !Author : Louisda16th a.k.a Ashwith J. Rego
    !Description: The subroutine is based on two key points:
    !1] A determinant is unaltered when row operations are performed: Hence, using this principle,
    !row operations (column operations would work as well) are used
    !to convert the matrix into upper traingular form
    !2]The determinant of a triangular matrix is obtained by finding the product of the diagonal elements
    !
    USE kinds,        ONLY : DP
    !
    IMPLICIT NONE
    REAL(KIND=DP), DIMENSION(n,n) :: matrix
    INTEGER, INTENT(IN)           :: n
    REAL(KIND=DP)                 :: FindDet, m, temp
    INTEGER                       :: i, j, k, l
    LOGICAL                       :: DetExists = .TRUE.
    !
    l = 1
    !Convert to upper triangular form
    DO k = 1, n-1
        IF (matrix(k,k) == 0.0d0) THEN
            DetExists = .FALSE.
            DO i = k+1, n
                IF (matrix(i,k) /= 0.0d0) THEN
                    DO j = 1, n
                        temp = matrix(i,j)
                        matrix(i,j)= matrix(k,j)
                        matrix(k,j) = temp
                    END DO
                    DetExists = .TRUE.
                    l=-l
                    EXIT
                ENDIF
            END DO
            IF (DetExists .EQV. .FALSE.) THEN
                FindDet = 0.0d0
                RETURN
            END IF
        ENDIF
        DO j = k+1, n
            m = matrix(j,k)/matrix(k,k)
            DO i = k+1, n
                matrix(j,i) = matrix(j,i) - m*matrix(k,i)
            END DO
        END DO
    END DO
    !
    !Calculate determinant by finding product of diagonal elements
    FindDet = l
    DO i = 1, n
        FindDet = FindDet * matrix(i,i)
    END DO
    !
  END FUNCTION FindDet

!----------------------------------------------------------------------------
END MODULE bte_func
!----------------------------------------------------------------------------
