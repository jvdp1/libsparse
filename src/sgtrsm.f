C
C  =========== DOCUMENTATION ===========
C
C Online html documentation available at
C            http://www.netlib.org/lapack/explore-html/
C
C  Definition:
C  ===========
C
C       SUBROUTINE STRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
C
C       .. Scalar Arguments ..
C       REAL ALPHA
C       INTEGER LDA,LDB,M,N
C       CHARACTER DIAG,SIDE,TRANSA,UPLO
C       ..
C       .. Array Arguments ..
C       REAL A(LDA,*),B(LDB,*)
C       ..
C
C
C> \par Purpose:
C  =============
C>
C> \verbatim
C>
C> STRSM  solves one of the matrix equations
C>
C>    op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
C>
C> where alpha is a scalar, X and B are m by n matrices, A is a unit,
Cor
C> non-unit,  upper or lower triangular matrix  and  op( A )  is one
Cof
C>
C>    op( A ) = A   or   op( A ) = A**T.
C>
C> The matrix X is overwritten on B.
C> \endverbatim
C
C  Arguments:
C  ==========
C
C> \param[in] SIDE
C> \verbatim
C>          SIDE is CHARACTER*1
C>           On entry, SIDE specifies whether op( A ) appears on the
Cleft
C>           or right of X as follows:
C>
C>              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
C>
C>              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
C> \endverbatim
C>
C> \param[in] UPLO
C> \verbatim
C>          UPLO is CHARACTER*1
C>           On entry, UPLO specifies whether the matrix A is an upper
Cor
C>           lower triangular matrix as follows:
C>
C>              UPLO = 'U' or 'u'   A is an upper triangular matrix.
C>
C>              UPLO = 'L' or 'l'   A is a lower triangular matrix.
C> \endverbatim
C>
C> \param[in] TRANSA
C> \verbatim
C>          TRANSA is CHARACTER*1
C>           On entry, TRANSA specifies the form of op( A ) to be used
Cin
C>           the matrix multiplication as follows:
C>
C>              TRANSA = 'N' or 'n'   op( A ) = A.
C>
C>              TRANSA = 'T' or 't'   op( A ) = A**T.
C>
C>              TRANSA = 'C' or 'c'   op( A ) = A**T.
C> \endverbatim
C>
C> \param[in] DIAG
C> \verbatim
C>          DIAG is CHARACTER*1
C>           On entry, DIAG specifies whether or not A is unit
Ctriangular
C>           as follows:
C>
C>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
C>
C>              DIAG = 'N' or 'n'   A is not assumed to be unit
C>                                  triangular.
C> \endverbatim
C>
C> \param[in] M
C> \verbatim
C>          M is INTEGER
C>           On entry, M specifies the number of rows of B. M must be
Cat
C>           least zero.
C> \endverbatim
C>
C> \param[in] N
C> \verbatim
C>          N is INTEGER
C>           On entry, N specifies the number of columns of B.  N must
Cbe
C>           at least zero.
C> \endverbatim
C>
C> \param[in] ALPHA
C> \verbatim
C>          ALPHA is REAL
C>           On entry,  ALPHA specifies the scalar  alpha. When  alpha
Cis
C>           zero then  A is not referenced and  B need not be set
Cbefore
C>           entry.
C> \endverbatim
C>
C> \param[in] A
C> \verbatim
C>          A is REAL array, dimension ( LDA, k ),
C>           where k is m when SIDE = 'L' or 'l'
C>             and k is n when SIDE = 'R' or 'r'.
C>           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by
Ck
C>           upper triangular part of the array  A must contain the
Cupper
C>           triangular matrix  and the strictly lower triangular part
Cof
C>           A is not referenced.
C>           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by
Ck
C>           lower triangular part of the array  A must contain the
Clower
C>           triangular matrix  and the strictly upper triangular part
Cof
C>           A is not referenced.
C>           Note that when  DIAG = 'U' or 'u',  the diagonal elements
Cof
C>           A  are not referenced either,  but are assumed to be
Cunity.
C> \endverbatim
C>
C> \param[in] LDA
C> \verbatim
C>          LDA is INTEGER
C>           On entry, LDA specifies the first dimension of A as
Cdeclared
C>           in the calling (sub) program.  When  SIDE = 'L' or 'l'
Cthen
C>           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or
C'r'
C>           then LDA must be at least max( 1, n ).
C> \endverbatim
C>
C> \param[in,out] B
C> \verbatim
C>          B is REAL array, dimension ( LDB, N )
C>           Before entry,  the leading  m by n part of the array  B
Cmust
C>           contain  the  right-hand  side  matrix  B,  and  on exit
Cis
C>           overwritten by the solution matrix  X.
C> \endverbatim
C>
C> \param[in] LDB
C> \verbatim
C>          LDB is INTEGER
C>           On entry, LDB specifies the first dimension of B as
Cdeclared
C>           in  the  calling  (sub)  program.   LDB  must  be  at
Cleast
C>           max( 1, m ).
C> \endverbatim
C
C  Authors:
C  ========
C
C> \author Univ. of Tennessee
C> \author Univ. of California Berkeley
C> \author Univ. of Colorado Denver
C> \author NAG Ltd.
C
C> \date December 2016
C
C> \ingroup single_blas_level3
C
C> \par Further Details:
C  =====================
C>
C> \verbatim
C>
C>  Level 3 Blas routine.
C>
C>
C>  -- Written on 8-February-1989.
C>     Jack Dongarra, Argonne National Laboratory.
C>     Iain Duff, AERE Harwell.
C>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
C>     Sven Hammarling, Numerical Algorithms Group Ltd.
C> \endverbatim
C>
C  =====================================================================
       SUBROUTINE strsm(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
C
C  -- Reference BLAS level3 routine (version 3.7.0) --
C  -- Reference BLAS is a software package provided by Univ. of
C  Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG
C  Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      REAL ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
C     ..
C     .. Array Arguments ..
      REAL A(lda,*),B(ldb,*)
C     ..
C
C  =====================================================================
C
C     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL lsame
C     ..
C     .. External Subroutines ..
      EXTERNAL xerbla
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC max
C     ..
C     .. Local Scalars ..
      REAL atmp
      REAL TEMP
      INTEGER I,INFO,J,K,NROWA
      LOGICAL LSIDE,NOUNIT,UPPER
C     ..
C     .. Parameters ..
      REAL ONE,ZERO
      REAL TOL
      parameter(one=1.0e+0,zero=0.0e+0)
      parameter(tol=1.0e-10)
C     ..
C
C     Test the input parameters.
C
      lside = lsame(side,'L')
      IF (lside) THEN
          nrowa = m
      ELSE
          nrowa = n
      END IF
      nounit = lsame(diag,'N')
      upper = lsame(uplo,'U')
C
      info = 0
      IF ((.NOT.lside) .AND. (.NOT.lsame(side,'R'))) THEN
          info = 1
      ELSE IF ((.NOT.upper) .AND. (.NOT.lsame(uplo,'L'))) THEN
          info = 2
      ELSE IF ((.NOT.lsame(transa,'N')) .AND.
     +         (.NOT.lsame(transa,'T')) .AND.
     +         (.NOT.lsame(transa,'C'))) THEN
          info = 3
      ELSE IF ((.NOT.lsame(diag,'U')) .AND. (.NOT.lsame(diag,'N')))THEN
          info = 4
      ELSE IF (m.LT.0) THEN
          info = 5
      ELSE IF (n.LT.0) THEN
          info = 6
      ELSE IF (lda.LT.max(1,nrowa)) THEN
          info = 9
      ELSE IF (ldb.LT.max(1,m)) THEN
          info = 11
      END IF
      IF (info.NE.0) THEN
          CALL xerbla('STRSM ',info)
          RETURN
      END IF
C
C     Quick return if possible.
C
      IF (m.EQ.0 .OR. n.EQ.0) RETURN
C
C     And when  alpha.eq.zero.
C
      IF (alpha.EQ.zero) THEN
          DO 20 j = 1,n
              DO 10 i = 1,m
                  b(i,j) = zero
   10         CONTINUE
   20     CONTINUE
          RETURN
      END IF
C
C     Start the operations.
C
      IF (lside) THEN
          IF (lsame(transa,'N')) THEN
C
C           Form  B := alpha*inv( A )*B.
C
              IF (upper) THEN
                  DO 60 j = 1,n
                      IF (alpha.NE.one) THEN
                          DO 30 i = 1,m
                              b(i,j) = alpha*b(i,j)
   30                     CONTINUE
                      END IF
                      DO 50 k = m,1,-1
                          IF (b(k,j).NE.zero) THEN
!                              IF (nounit) b(k,j) = b(k,j)/a(k,k)
                              IF (nounit)then
                               atmp=0.0e+0
                               if(a(k,k).gt.tol)atmp=1.0e+0/a(k,k)
                               b(k,j) = b(k,j)*atmp
                              endif
                              DO 40 i = 1,k - 1
                                  b(i,j) = b(i,j) - b(k,j)*a(i,k)
   40                         CONTINUE
                          END IF
   50                 CONTINUE
   60             CONTINUE
              ELSE
                  DO 100 j = 1,n
                      IF (alpha.NE.one) THEN
                          DO 70 i = 1,m
                              b(i,j) = alpha*b(i,j)
   70                     CONTINUE
                      END IF
                      DO 90 k = 1,m
                          IF (b(k,j).NE.zero) THEN
!                              IF (nounit) b(k,j) = b(k,j)/a(k,k)
                              IF (nounit)then
                               atmp=0.0e+0
                               if(a(k,k).gt.tol)atmp=1.0e+0/a(k,k)
                               b(k,j) = b(k,j)*atmp
                              endif
                              DO 80 i = k + 1,m
                                  b(i,j) = b(i,j) - b(k,j)*a(i,k)
   80                         CONTINUE
                          END IF
   90                 CONTINUE
  100             CONTINUE
              END IF
          ELSE
C
C           Form  B := alpha*inv( A**T )*B.
C
              IF (upper) THEN
                  DO 130 j = 1,n
                      DO 120 i = 1,m
                          temp = alpha*b(i,j)
                          DO 110 k = 1,i - 1
                              temp = temp - a(k,i)*b(k,j)
  110                     CONTINUE
!                          IF (nounit) temp = temp/a(i,i)
                          IF (nounit)then
                           atmp=0.0e+0
                           if(a(i,i).gt.tol)atmp=1.0e+0/a(i,i)
                           temp = temp*atmp
                          endif
                          b(i,j) = temp
  120                 CONTINUE
  130             CONTINUE
              ELSE
                  DO 160 j = 1,n
                      DO 150 i = m,1,-1
                          temp = alpha*b(i,j)
                          DO 140 k = i + 1,m
                              temp = temp - a(k,i)*b(k,j)
  140                     CONTINUE
!                          IF (nounit) temp = temp/a(i,i)
                          IF (nounit)then
                           atmp=0.0e+0
                           if(a(i,i).gt.tol)atmp=1.0e+0/a(i,i)
                           temp = temp*atmp
                          endif
                          b(i,j) = temp
  150                 CONTINUE
  160             CONTINUE
              END IF
          END IF
      ELSE
          IF (lsame(transa,'N')) THEN
C
C           Form  B := alpha*B*inv( A ).
C
              IF (upper) THEN
                  DO 210 j = 1,n
                      IF (alpha.NE.one) THEN
                          DO 170 i = 1,m
                              b(i,j) = alpha*b(i,j)
  170                     CONTINUE
                      END IF
                      DO 190 k = 1,j - 1
                          IF (a(k,j).NE.zero) THEN
                              DO 180 i = 1,m
                                  b(i,j) = b(i,j) - a(k,j)*b(i,k)
  180                         CONTINUE
                          END IF
  190                 CONTINUE
                      IF (nounit) THEN
!                          temp = one/a(j,j)
                          atmp=0.0e+0
                          if(a(j,j).gt.tol)atmp=1.0e+0/a(j,j)
                          temp = one*atmp
                          DO 200 i = 1,m
                              b(i,j) = temp*b(i,j)
  200                     CONTINUE
                      END IF
  210             CONTINUE
              ELSE
                  DO 260 j = n,1,-1
                      IF (alpha.NE.one) THEN
                          DO 220 i = 1,m
                              b(i,j) = alpha*b(i,j)
  220                     CONTINUE
                      END IF
                      DO 240 k = j + 1,n
                          IF (a(k,j).NE.zero) THEN
                              DO 230 i = 1,m
                                  b(i,j) = b(i,j) - a(k,j)*b(i,k)
  230                         CONTINUE
                          END IF
  240                 CONTINUE
                      IF (nounit) THEN
!                          temp = one/a(j,j)
                          atmp=0.0e+0
                          if(a(j,j).gt.tol)atmp=1.0e+0/a(j,j)
                          temp = one*atmp
                          DO 250 i = 1,m
                              b(i,j) = temp*b(i,j)
  250                     CONTINUE
                      END IF
  260             CONTINUE
              END IF
          ELSE
C
C           Form  B := alpha*B*inv( A**T ).
C
              IF (upper) THEN
                  DO 310 k = n,1,-1
                      IF (nounit) THEN
!                          temp = one/a(k,k)
                          atmp=0.e+0
                          if(a(k,k).gt.tol)atmp=1.e+0/a(k,k)
                          temp = one*atmp
                          DO 270 i = 1,m
                              b(i,k) = temp*b(i,k)
  270                     CONTINUE
                      END IF
                      DO 290 j = 1,k - 1
                          IF (a(j,k).NE.zero) THEN
                              temp = a(j,k)
                              DO 280 i = 1,m
                                  b(i,j) = b(i,j) - temp*b(i,k)
  280                         CONTINUE
                          END IF
  290                 CONTINUE
                      IF (alpha.NE.one) THEN
                          DO 300 i = 1,m
                              b(i,k) = alpha*b(i,k)
  300                     CONTINUE
                      END IF
  310             CONTINUE
              ELSE
                  DO 360 k = 1,n
                      IF (nounit) THEN
!                          temp = one/a(k,k)
                          atmp=0.0e+0
                          if(a(k,k).gt.tol)atmp=1.0e+0/a(k,k)
                          temp = one*atmp
                          DO 320 i = 1,m
                              b(i,k) = temp*b(i,k)
  320                     CONTINUE
                      END IF
                      DO 340 j = k + 1,n
                          IF (a(j,k).NE.zero) THEN
                              temp = a(j,k)
                              DO 330 i = 1,m
                                  b(i,j) = b(i,j) - temp*b(i,k)
  330                         CONTINUE
                          END IF
  340                 CONTINUE
                      IF (alpha.NE.one) THEN
                          DO 350 i = 1,m
                              b(i,k) = alpha*b(i,k)
  350                     CONTINUE
                      END IF
  360             CONTINUE
              END IF
          END IF
      END IF
C
      RETURN
C
C     End of STRSM .
C
      END
