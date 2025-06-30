c#######################################################################
c#      this is subroutine "GSFCT" as given in : 'Computer Solution    #
c#      of Sparse Linear Systems' by A. George, J. Liu and             #
c#      E. Ng, 1994, pp. 180-183;                                      #
c#######################################################################

C----- SUBROUTINE GSFCT
C***************************************************************
C***************************************************************
C******     GSFCT ..... GENERAL SPARSE SYMMETRIC FACT     ******
C***************************************************************
C***************************************************************
C
C     PURPOSE - THIS SUBROUTINE PERFORMS THE SYMMETRIC
C        FACTORIZATION FOR A GENERAL SPARSE SYSTEM, STORED IN
C        THE COMPRESSED SUBSCRIPT DATA FORMAT.
C
C     INPUT PARAMETERS -
C        NEQNS - NUMBER OF EQUATIONS.
C        XLNZ - INDEX VECTOR FOR LNZ. XLNZ(I) POINTS TO THE
C               START OF NONZEROS IN COLUMN I OF FACTOR L.
C        (XNZSUB, NZSUB) - THE COMPRESSED SUBSCRIPT DATA
C               STRUCTURE FOR FACTOR L.
C
C     UPDATED PARAMETERS -
C        LNZ - ON INPUT, CONTAINS NONZEROS OF A, AND ON
C              RETURN, THE NONZEROS OF L.
C        DIAG - THE DIAGONAL OF L OVERWRITES THAT OF A.
C        IFLAG - THE ERROR FLAG. IT IS SET TO 1 IF A ZERO OR
C               NEGATIVE SQUARE ROOT OCCURS DURING THE
C               FACTORIZATION.
C        OPS - A DOUBLE PRECISION COMMON PARAMETER THAT IS
C              INCREMENTED BY THE NUMBER OF OPERATIONS
C              PERFORMED BY THE SUBROUTINE.
C
C     WORKING PARAMETERS -
C        LINK - AT STEP J, THE LIST IN
C                  LINK(J), LINK(LINK(J)), ...........
C               CONSISTS OF THOSE COLUMNS THAT WILL MODIFY
C               THE COLUMN L(*,J).
C        FIRST - TEMPORARY VECTOR TO POINT TO THE FIRST
C                NONZERO IN EACH COLUMN THAT WILL BE USED
C                NEXT FOR MODIFICATION.
C        TEMP - A TEMPORARY VECTOR TO ACCUMULATE MODIFICATIONS.
C
C***************************************************************
C
      SUBROUTINE GSFCT ( NEQNS, XLNZ, LNZ, XNZSUB, NZSUB, DIAG,
     1                   LINK, FIRST, TEMP, IFLAG )
C
C***************************************************************
C
         DOUBLE PRECISION COUNT, OPS
         COMMON /SPKOPS/ OPS
         REAL DIAG(1), LNZ(1), TEMP(1), DIAGJ, LJK
         INTEGER LINK(1), NZSUB(1)
         INTEGER FIRST(1), XLNZ(1), XNZSUB(1),
     1           I, IFLAG, II, ISTOP, ISTRT, ISUB, J,
     1           K, KFIRST, NEQNS, NEWK
C
C***************************************************************
C
C       ------------------------------
C       INITIALIZE WORKING VECTORS ...
C       ------------------------------
        DO 100 I = 1, NEQNS
           LINK(I) = 0
           TEMP(I) = 0.0E0
  100   CONTINUE
C       --------------------------------------------
C       COMPUTE COLUMN L(*,J) FOR J = 1,...., NEQNS.
C       --------------------------------------------
        DO 600 J = 1, NEQNS
C          -------------------------------------------
C          FOR EACH COLUMN L(*,K) THAT AFFECTS L(*,J).
C          -------------------------------------------
           DIAGJ = 0.0E0
           NEWK = LINK(J)
  200      K    = NEWK
           IF ( K .EQ. 0 ) GO TO 400
              NEWK = LINK(K)
C             ---------------------------------------
C             OUTER PRODUCT MODIFICATION OF L(*,J) BY
C             L(*,K) STARTING AT FIRST(K) OF L(*,K).
C             ---------------------------------------
              KFIRST = FIRST(K)
             LJK     = LNZ(KFIRST)
             DIAGJ = DIAGJ + LJK*LJK
             OPS = OPS + 1.0D0
             ISTRT = KFIRST + 1
             ISTOP = XLNZ(K+1) - 1
             IF ( ISTOP .LT. ISTRT ) GO TO 200
C               ------------------------------------------
C               BEFORE MODIFICATION, UPDATE VECTORS FIRST,
C               AND LINK FOR FUTURE MODIFICATION STEPS.
C               ------------------------------------------
                FIRST(K) = ISTRT
                I = XNZSUB(K) + (KFIRST-XLNZ(K)) + 1
                ISUB = NZSUB(I)
                LINK(K) = LINK(ISUB)
                LINK(ISUB) = K
C               ---------------------------------------
C               THE ACTUAL MOD IS SAVED IN VECTOR TEMP.
C               ---------------------------------------
                DO 300 II = ISTRT, ISTOP
                   ISUB = NZSUB(I)
                   TEMP(ISUB) = TEMP(ISUB) + LNZ(II)*LJK
                   I = I + 1
  300           CONTINUE
                COUNT = ISTOP - ISTRT + 1
                OPS = OPS + COUNT
             GO TO 200
C         ----------------------------------------------
C         APPLY THE MODIFICATIONS ACCUMULATED IN TEMP TO
C         COLUMN L(*,J).
C         ----------------------------------------------
  400     DIAGJ = DIAG(J) - DIAGJ
          IF ( DIAGJ .LE. 0.0E0 ) GO TO 700
          DIAGJ = SQRT(DIAGJ)
          DIAG(J) = DIAGJ
          ISTRT = XLNZ(J)
          ISTOP = XLNZ(J+1) - 1
          IF ( ISTOP .LT. ISTRT ) GO TO 600
             FIRST(J) = ISTRT
             I = XNZSUB(J)
             ISUB = NZSUB(I)
             LINK(J) = LINK(ISUB)
             LINK(ISUB) = J
             DO 500 II = ISTRT, ISTOP
                ISUB = NZSUB(I)
                LNZ(II) = ( LNZ(II)-TEMP(ISUB) ) / DIAGJ
                TEMP(ISUB) = 0.0E0
                I = I + 1
  500        CONTINUE
             COUNT = ISTOP - ISTRT + 1
             OPS = OPS + COUNT
  600   CONTINUE
        RETURN
C       ------------------------------------------------------
C       ERROR - ZERO OR NEGATIVE SQUARE ROOT IN FACTORIZATION.
C       ------------------------------------------------------
  700   IFLAG = 1
        RETURN
        END
