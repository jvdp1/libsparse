c#######################################################################
c#      this is subroutine "SMBFCT" as given in : 'Computer Solutions  #
c#      of Large Sparse Positive Definite Systems' by A. George and    #
c#      J.W.-H. Liu, 1981, Prentice Hall, Inc. Englewood Cliffs,       #
c#      New Jersey 07632, pp. 149-151;                                 #
c#      Modified for my own                                            #
c#######################################################################


C----- SUBROUTINE SMBFCT
C****************************************************************          1.
C****************************************************************          2.
C*********     SMBFCT ..... SYMBOLIC FACTORIZATION       ********          3.
C****************************************************************          4.
C****************************************************************          5.
C                                                                          6.
C     PURPOSE - THIS ROUTINE PERFORMS SYMBOLIC FACTORIZATION               7.
C        ON A PERMUTED LINEAR SYSTEM AND IT ALSO SETS UP THE               8.
C        COMPRESSED DATA STRUCTURE FOR THE SYSTEM.                         9.
C                                                                         10.
C     INPUT PARAMETERS -                                                  11.
C        NEQNS - NUMBER OF EQUATIONS.                                     12.
C        (XADJ, ADJNCY) - THE ADJACENCY STRUCTURE.                        13.
C        (PERM, INVP) - THE PERMUTATION VECTOR AND ITS INVERSE.           14.
C                                                                         15.
C     UPDATED PARAMETERS -                                                16.
C        MAXSUB - SIZE OF THE SUBSCRIPT ARRAY NZSUB.  ON RETURN,          17.
C               IT CONTAINS THE NUMBER OF SUBSCRIPTS USED                 18.
C                                                                         19.
C     OUTPUT PARAMETERS -                                                 20.
C        XLNZ - INDEX INTO THE NONZERO STORAGE VECTOR LNZ.                21.
C        (XNZSUB, NZSUB) - THE COMPRESSED SUBSCRIPT VECTORS.              22.
C        MAXLNZ - THE NUMBER OF NONZEROS FOUND.                           23.
C        FLAG - ERROR FLAG.  POSITIVE VALUE INDICATES THAT.               24.
C               NZSUB ARRAY IS TOO SMALL.                                 25.
C                                                                         26.
C     WORKING PARAMETERS -                                                27.
C        MRGLNK - A VECTOR OF SIZE NEQNS.  AT THE KTH STEP,               28.
C               MRGLNK(K), MRGLNK(MRGLNK(K)) , .........                  29.
C               IS A LIST CONTAINING ALL THOSE COLUMNS L(*,J)             30.
C               WITH J LESS THAN K, SUCH THAT ITS FIRST OFF-              31.
C               DIAGONAL NONZERO IS L(K,J).  THUS, THE                    32.
C               NONZERO STRUCTURE OF COLUMN L(*,K) CAN BE FOUND           33.
C               BY MERGING THAT OF SUCH COLUMNS L(*,J) WITH               34.
C               THE STRUCTURE OF A(*,K).                                  35.
C        RCHLNK - A VECTOR OF SIZE NEQNS.  IT IS USED TO ACCUMULATE       36.
C               THE STRUCTURE OF EACH COLUMN L(*,K).  AT THE              37.
C               END OF THE KTH STEP,                                      38.
C                   RCHLNK(K), RCHLNK(RCHLNK(K)), ........                39.
C               IS THE LIST OF POSITIONS OF NONZEROS IN COLUMN K          40.
C               OF THE FACTOR L.                                          41.
C        MARKER  - AN INTEGER VECTOR OF LENGTH NEQNS. IT IS USED          42.
C               TO TEST IF MASS SYMBOLIC ELIMINATION CAN BE               43.
C               PERFORMED.  THAT IS, IT IS USED TO CHECK WHETHER          44.
C               THE STRUCTURE OF THE CURRENT COLUMN K BEING               45.
C               PROCESSED IS COMPLETELY DETERMINED BY THE SINGLE          46.
C               COLUMN MRGLNK(K).                                         47.
C                                                                         48.
C****************************************************************         49.
C                                                                         50.
      SUBROUTINE  SMBFCT ( NEQNS, XADJ, ADJNCY, PERM, INVP,               51.
     1                     XLNZ, MAXLNZ, XNZSUB, NZSUB, MAXSUB,           52.
     1                     RCHLNK, MRGLNK, MARKER, FLAG )                 53.
C                                                                         54.
C****************************************************************         55.
C                                                                         56.
         INTEGER ADJNCY(*), INVP(*), MRGLNK(*), NZSUB(*),                 57.
     1           PERM(*), RCHLNK(*), MARKER(*)                            58.
         INTEGER XADJ(*), XLNZ(*), XNZSUB(*),                             59.
     1           FLAG, I, INZ, J, JSTOP, JSTRT, K, KNZ,                   60.
     1           KXSUB, MRGK, LMAX, M, MAXLNZ, MAXSUB,                    61.
     1           NABOR, NEQNS, NODE, NP1, NZBEG, NZEND,                   62.
     1           RCHM, MRKFLG                                             63.
C                                                                         64.
C****************************************************************         65.
C                                                                         66.
C       ------------------                                                67.
C       INITIALIZATION ...                                                68.
C       ------------------                                                69.
        NZBEG = 1                                                         70.
        NZEND = 0                                                         71.
        XLNZ(1) = 1                                                       72.
           MRGLNK(1:NEQNS) = 0                                            74.
           MARKER(1:NEQNS) = 0                                            75.
C       --------------------------------------------------                77.
C       FOR EACH COLUMN ......... .  KNZ COUNTS THE NUMBER                78.
C       OF NONZEROS IN COLUMN K ACCUMULATED IN RCHLNK.                    79.
C       --------------------------------------------------                80.
        NP1 = NEQNS + 1                                                   81.
        DO 1500  K = 1, NEQNS                                             82.
           KNZ = 0                                                        83.
           MRGK = MRGLNK(K)                                               84.
           MRKFLG = 0                                                     85.
           MARKER(K) = K                                                  86.
           IF (MRGK .NE. 0 ) MARKER(K) = MARKER(MRGK)                     87.
           XNZSUB(K) = NZEND                                              88.
           NODE = PERM(K)                                                 89.
           JSTRT = XADJ(NODE)                                             90.
           JSTOP = XADJ(NODE+1) - 1                                       91.
           IF (JSTRT.GT.JSTOP)  GO TO 1500                                92.
C          -------------------------------------------                    93.
C          USE RCHLNK TO LINK THROUGH THE STRUCTURE OF                    94.
C          A(*,K) BELOW DIAGONAL                                          95.
C          -------------------------------------------                    96.
           RCHLNK(K) = NP1                                                97.
           DO_300: DO J = JSTRT, JSTOP                                    98.
              NABOR = ADJNCY(J)                                           99.
              NABOR = INVP(NABOR)                                        100.
              IF ( NABOR .LE. K )  cycle DO_300                          101.
                 RCHM = K                                                102.
                 DO_200: do
                 M = RCHM                                                103.
                 RCHM = RCHLNK(M)                                        104.
                 IF ( RCHM .gt. NABOR )  exit DO_200                     105.
                 enddo DO_200
                    KNZ = KNZ+1                                          106.
                    RCHLNK(M) = NABOR                                    107.
                    RCHLNK(NABOR) = RCHM                                 108.
                    IF ( MARKER(NABOR) .NE. MARKER(K) )  MRKFLG = 1      109.
           end do DO_300                                                 110.
C          --------------------------------------                        111.
C          TEST FOR MASS SYMBOLIC ELIMINATION ...                        112.
C          --------------------------------------                        113.
           LMAX = 0                                                      114.
           IF ( MRKFLG .NE. 0 .OR. MRGK .EQ. 0 )then                     115.
           else
           IF ( MRGLNK(MRGK) .NE. 0 )then                                116.
           else
           XNZSUB(K) = XNZSUB(MRGK) + 1                                  117.
           KNZ = XLNZ(MRGK+1) - (XLNZ(MRGK) + 1)                         118.
           GO TO 1400                                                    119.
           endif
           endif
C          -----------------------------------------------               120.
C          LINK THROUGH EACH COLUMN I THAT AFFECTS L(*,K).               121.
C          -----------------------------------------------               122.
           I = K                                                         123.
           DO_400: do
           I = MRGLNK(I)                                                 124.
           IF (I.EQ.0)  exit DO_400 !800                                 125.
              INZ = XLNZ(I+1) - (XLNZ(I)+1)                              126.
              JSTRT = XNZSUB(I) +  1                                     127.
              JSTOP = XNZSUB(I) + INZ                                    128.
              if_500: IF (INZ.LE.LMAX)then                               129.
              else
                 LMAX = INZ                                              130.
                 XNZSUB(K) = JSTRT                                       131.
              endif if_500
C             -----------------------------------------------            132.
C             MERGE STRUCTURE OF L(*,I) IN NZSUB INTO RCHLNK.            133.
C             -----------------------------------------------            134.
              RCHM = K                                                   135.
              DO_700: DO J = JSTRT, JSTOP                                    136.
                 NABOR = NZSUB(J)                                        137.
                 DO_600: do
                 M = RCHM                                                138.
                 RCHM = RCHLNK(M)                                        139.
                 IF (RCHM.ge.NABOR)  exit DO_600                         140.
                 enddo DO_600
                 IF (RCHM.EQ.NABOR)  cycle DO_700                        141.
                    KNZ = KNZ+1                                          142.
                    RCHLNK(M) = NABOR                                    143.
                    RCHLNK(NABOR) = RCHM                                 144.
                    RCHM = NABOR                                         145.
              enddo DO_700                                               146.
              enddo DO_400                                                  147.
C          ------------------------------------------------------        148.
C          CHECK IF SUBSCRIPTS DUPLICATE THOSE OF ANOTHER COLUMN.        149.
C          ------------------------------------------------------        150.
           IF (KNZ.EQ.LMAX)then                                          151.
           else
C             -----------------------------------------------            152.
C             OR IF TAIL OF K-1ST COLUMN MATCHES HEAD OF KTH.            153.
C             -----------------------------------------------            154.
              IF (NZBEG.GT.NZEND)then                                    155.
              else
                 I = RCHLNK(K)                                           156.
                 DO 900 JSTRT=NZBEG,NZEND                                157.
                    IF (NZSUB(JSTRT)-I)  900, 1000, 1200                 158.
  900            CONTINUE                                                159.
                 GO TO 1200                                              160.
 1000            XNZSUB(K) = JSTRT                                       161.
                 DO J=JSTRT,NZEND                                        162.
                    IF (NZSUB(J).NE.I)  GO TO 1200                       163.
                    I = RCHLNK(I)                                        164.
                    IF (I.GT.NEQNS)  GO TO 1400                          165.
                 END DO
                 NZEND = JSTRT - 1                                       167.
               endif
C             ----------------------------------------                   168.
C             COPY THE STRUCTURE OF L(*,K) FROM RCHLNK                   169.
C             TO THE DATA STRUCTURE (XNZSUB, NZSUB).                     170.
C             ----------------------------------------                   171.
 1200         NZBEG = NZEND +  1                                         172.
              NZEND = NZEND + KNZ                                        173.
              IF (NZEND.GT.MAXSUB) then
                 FLAG = 1
                 RETURN
              END IF
              I = K                                                      175.
              DO J=NZBEG,NZEND                                           176.
                 I = RCHLNK(I)                                           177.
                 NZSUB(J) = I                                            178.
                 MARKER(I) = K                                           179.
              END DO
              XNZSUB(K) = NZBEG                                          181.
              MARKER(K) = K                                              182.
           endif
C          --------------------------------------------------------      183.
C          UPDATE THE VECTOR MRGLNK.  NOTE COLUMN L(*,K) JUST FOUND      184.
C          IS REQUIRED TO DETERMINE COLUMN L(*,J), WHERE                 185.
C          L(J,K) IS THE FIRST NONZERO IN L(*,K) BELOW DIAGONAL.         186.
C          --------------------------------------------------------      187.
 1400      IF (KNZ.LE.1)  GO TO 1500                                     188.
              KXSUB = XNZSUB(K)                                          189.
              I = NZSUB(KXSUB)                                           190.
              MRGLNK(K) = MRGLNK(I)                                      191.
              MRGLNK(I) = K                                              192.
 1500      XLNZ(K+1) = XLNZ(K) + KNZ                                     193.
        MAXLNZ = XLNZ(NEQNS) - 1                                         194.
        MAXSUB = XNZSUB(NEQNS)                                           195.
        XNZSUB(NEQNS+1) = XNZSUB(NEQNS)                                  196.
        FLAG = 0                                                         197.
        RETURN                                                           198.
C       ----------------------------------------------------             199.
C       ERROR - INSUFFICIENT STORAGE FOR NONZERO SUBSCRIPTS.             200.
C       ----------------------------------------------------             201.
        FLAG = 1                                                         202.
        RETURN                                                           203.
        END                                                              204.



