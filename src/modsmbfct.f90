module modsmbfct
   implicit none(external, type)
   private
   public :: smbfct

   interface smbfct
      module procedure smbfct_int32
   end interface
contains

!c#######################################################################
!c#      this is subroutine "SMBFCT" as given in : 'Computer Solutions  #
!c#      of Large Sparse Positive Definite Systems' by A. George and    #
!c#      J.W.-H. Liu, 1981, Prentice Hall, Inc. Englewood Cliffs,       #
!c#      New Jersey 07632, pp. 149-151;                                 #
!c#      Translated to Modern Fotran by Jeremie Vandenplas 07/02/2024   #
!c#######################################################################

!C----- SUBROUTINE SMBFCT
!C****************************************************************          1.
!C****************************************************************          2.
!C*********     SMBFCT ..... SYMBOLIC FACTORIZATION       ********          3.
!C****************************************************************          4.
!C****************************************************************          5.
!C                                                                          6.
!C     PURPOSE - THIS ROUTINE PERFORMS SYMBOLIC FACTORIZATION               7.
!C        ON A PERMUTED LINEAR SYSTEM AND IT ALSO SETS UP THE               8.
!C        COMPRESSED DATA STRUCTURE FOR THE SYSTEM.                         9.
!C                                                                         10.
!C     INPUT PARAMETERS -                                                  11.
!C        NEQNS - NUMBER OF EQUATIONS.                                     12.
!C        (XADJ, ADJNCY) - THE ADJACENCY STRUCTURE.                        13.
!C        (PERM, INVP) - THE PERMUTATION VECTOR AND ITS INVERSE.           14.
!C                                                                         15.
!C     UPDATED PARAMETERS -                                                16.
!C        MAXSUB - SIZE OF THE SUBSCRIPT ARRAY NZSUB.  ON RETURN,          17.
!C               IT CONTAINS THE NUMBER OF SUBSCRIPTS USED                 18.
!C                                                                         19.
!C     OUTPUT PARAMETERS -                                                 20.
!C        XLNZ - INDEX INTO THE NONZERO STORAGE VECTOR LNZ.                21.
!C        (XNZSUB, NZSUB) - THE COMPRESSED SUBSCRIPT VECTORS.              22.
!C        MAXLNZ - THE NUMBER OF NONZEROS FOUND.                           23.
!C        FLAG - ERROR FLAG.  POSITIVE VALUE INDICATES THAT.               24.
!C               NZSUB ARRAY IS TOO SMALL.                                 25.
!C                                                                         26.
!C     WORKING PARAMETERS -                                                27.
!C        MRGLNK - A VECTOR OF SIZE NEQNS.  AT THE KTH STEP,               28.
!C               MRGLNK(K), MRGLNK(MRGLNK(K)) , .........                  29.
!C               IS A LIST CONTAINING ALL THOSE COLUMNS L(*,J)             30.
!C               WITH J LESS THAN K, SUCH THAT ITS FIRST OFF-              31.
!C               DIAGONAL NONZERO IS L(K,J).  THUS, THE                    32.
!C               NONZERO STRUCTURE OF COLUMN L(*,K) CAN BE FOUND           33.
!C               BY MERGING THAT OF SUCH COLUMNS L(*,J) WITH               34.
!C               THE STRUCTURE OF A(*,K).                                  35.
!C        RCHLNK - A VECTOR OF SIZE NEQNS.  IT IS USED TO ACCUMULATE       36.
!C               THE STRUCTURE OF EACH COLUMN L(*,K).  AT THE              37.
!C               END OF THE KTH STEP,                                      38.
!C                   RCHLNK(K), RCHLNK(RCHLNK(K)), ........                39.
!C               IS THE LIST OF POSITIONS OF NONZEROS IN COLUMN K          40.
!C               OF THE FACTOR L.                                          41.
!C        MARKER  - AN INTEGER VECTOR OF LENGTH NEQNS. IT IS USED          42.
!C               TO TEST IF MASS SYMBOLIC ELIMINATION CAN BE               43.
!C               PERFORMED.  THAT IS, IT IS USED TO CHECK WHETHER          44.
!C               THE STRUCTURE OF THE CURRENT COLUMN K BEING               45.
!C               PROCESSED IS COMPLETELY DETERMINED BY THE SINGLE          46.
!C               COLUMN MRGLNK(K).                                         47.
!C                                                                         48.
!C****************************************************************         49.
!C                                                                         50.
   subroutine smbfct_int32(neqns, xadj, adjncy, perm, invp, &
                           xlnz, maxlnz, xnzsub, nzsub, maxsub, &
                           rchlnk, mrglnk, marker, flag)
!C                                                                         54.
!C****************************************************************         55.
!C                                                                         56.
      integer(int32), intent(in) :: neqns, xadj(:), adjncy(:), perm(:), invp(:)
      integer(int32), intent(out) :: xlnz(:), maxlnz, xnzsub(:), nzsub(:) &
                                     , rchlnk(:), mrglnk(:), marker(:) &
                                     , flag
      integer(int32), intent(inout) :: maxsub

      integer(int32) :: i, inz, j, jstop, jstrt, k, knz, &
                        kxsub, mrgk, lmax, m, &
                        nabor, node, np1, nzbeg, nzend, &
                        rchm, mrkflg
      logical :: go_to_1400
!C                                                                         64.
!C****************************************************************         65.
!C                                                                         66.
!C       ------------------                                                67.
!C       INITIALIZATION ...                                                68.
!C       ------------------                                                69.
      nzbeg = 1
      nzend = 0
      xlnz(1) = 1
      mrglnk(1:neqns) = 0
      marker(1:neqns) = 0
!C       --------------------------------------------------                77.
!C       FOR EACH COLUMN ......... .  KNZ COUNTS THE NUMBER                78.
!C       OF NONZEROS IN COLUMN K ACCUMULATED IN RCHLNK.                    79.
!C       --------------------------------------------------                80.
      np1 = neqns + 1
      do_1500: do k = 1, neqns
         knz = 0
         mrgk = mrglnk(k)
         mrkflg = 0
         marker(k) = k
         if (mrgk .ne. 0) marker(k) = marker(mrgk)
         xnzsub(k) = nzend
         node = perm(k)
         jstrt = xadj(node)
         jstop = xadj(node + 1) - 1
         if_1500_1: if (jstrt .le. jstop) then
!C          -------------------------------------------                    93.
!C          USE RCHLNK TO LINK THROUGH THE STRUCTURE OF                    94.
!C          A(*,K) BELOW DIAGONAL                                          95.
!C          -------------------------------------------                    96.
            rchlnk(k) = np1
            do_300: do j = jstrt, jstop
               nabor = adjncy(j)
               nabor = invp(nabor)
               if (nabor .le. k) cycle do_300
               rchm = k
               do_200: do
                  m = rchm
                  rchm = rchlnk(m)
                  if (rchm .gt. nabor) exit do_200
               end do do_200
               knz = knz + 1
               rchlnk(m) = nabor
               rchlnk(nabor) = rchm
               if (marker(nabor) .ne. marker(k)) mrkflg = 1
            end do do_300
!C          --------------------------------------                        111.
!C          TEST FOR MASS SYMBOLIC ELIMINATION ...                        112.
!C          --------------------------------------                        113.
            lmax = 0
            go_to_1400 = .false.
            if (mrkflg .ne. 0 .or. mrgk .eq. 0) then
            else
               if (mrglnk(mrgk) .eq. 0) then
                  xnzsub(k) = xnzsub(mrgk) + 1
                  knz = xlnz(mrgk + 1) - (xlnz(mrgk) + 1)
                  go_to_1400 = .true.
               end if
            end if
            if_go_to_1400: if (.not. go_to_1400) then
!C          -----------------------------------------------               120.
!C          LINK THROUGH EACH COLUMN I THAT AFFECTS L(*,K).               121.
!C          -----------------------------------------------               122.
               i = k
               do_400: do
                  i = mrglnk(i)
                  if (i .eq. 0) exit do_400 !800
                  inz = xlnz(i + 1) - (xlnz(i) + 1)
                  jstrt = xnzsub(i) + 1
                  jstop = xnzsub(i) + inz
                  if_500: if (inz .gt. lmax) then
                     lmax = inz
                     xnzsub(k) = jstrt
                  end if if_500
!C             -----------------------------------------------            132.
!C             MERGE STRUCTURE OF L(*,I) IN NZSUB INTO RCHLNK.            133.
!C             -----------------------------------------------            134.
                  rchm = k
                  do_700: do j = jstrt, jstop
                     nabor = nzsub(j)
                     do_600: do
                        m = rchm
                        rchm = rchlnk(m)
                        if (rchm .ge. nabor) exit do_600
                     end do do_600
                     if (rchm .eq. nabor) cycle do_700
                     knz = knz + 1
                     rchlnk(m) = nabor
                     rchlnk(nabor) = rchm
                     rchm = nabor
                  end do do_700
               end do do_400
!C          ------------------------------------------------------        148.
!C          CHECK IF SUBSCRIPTS DUPLICATE THOSE OF ANOTHER COLUMN.        149.
!C          ------------------------------------------------------        150.
               if (knz .ne. lmax) then
!C             -----------------------------------------------            152.
!C             OR IF TAIL OF K-1ST COLUMN MATCHES HEAD OF KTH.            153.
!C             -----------------------------------------------            154.
                  if_pre_1200: if (nzbeg .le. nzend) then
                     i = rchlnk(k)
                     do_900: do jstrt = nzbeg, nzend
                        if (nzsub(jstrt) - i .lt. 0) then
                           cycle do_900
                        elseif (nzsub(jstrt) - i .eq. 0) then
                           xnzsub(k) = jstrt
                           do j = jstrt, nzend
                              if (nzsub(j) .ne. i) exit if_pre_1200
                              i = rchlnk(i)
                              if (i .gt. neqns) exit if_go_to_1400
                           end do
                           nzend = jstrt - 1
                           exit if_pre_1200
                        elseif (nzsub(jstrt) - i .gt. 0) then
                           exit if_pre_1200
                        end if
                     end do do_900
                     exit if_pre_1200
                  end if if_pre_1200
!C             ----------------------------------------                   168.
!C             COPY THE STRUCTURE OF L(*,K) FROM RCHLNK                   169.
!C             TO THE DATA STRUCTURE (XNZSUB, NZSUB).                     170.
!C             ----------------------------------------                   171.
                  nzbeg = nzend + 1
                  nzend = nzend + knz
                  if (nzend .gt. maxsub) then
!C       ----------------------------------------------------             199.
!C       ERROR - INSUFFICIENT STORAGE FOR NONZERO SUBSCRIPTS.             200.
!C       ----------------------------------------------------             201.
                     flag = 1
                     return
                  end if
                  i = k
                  do j = nzbeg, nzend
                     i = rchlnk(i)
                     nzsub(j) = i
                     marker(i) = k
                  end do
                  xnzsub(k) = nzbeg
                  marker(k) = k
               end if
            end if if_go_to_1400
!C          --------------------------------------------------------      183.
!C          UPDATE THE VECTOR MRGLNK.  NOTE COLUMN L(*,K) JUST FOUND      184.
!C          IS REQUIRED TO DETERMINE COLUMN L(*,J), WHERE                 185.
!C          L(J,K) IS THE FIRST NONZERO IN L(*,K) BELOW DIAGONAL.         186.
!C          --------------------------------------------------------      187.
            if (knz .gt. 1) then
               kxsub = xnzsub(k)
               i = nzsub(kxsub)
               mrglnk(k) = mrglnk(i)
               mrglnk(i) = k
            end if
         end if if_1500_1
         xlnz(k + 1) = xlnz(k) + knz
      end do do_1500
      maxlnz = xlnz(neqns) - 1
      maxsub = xnzsub(neqns)
      xnzsub(neqns + 1) = xnzsub(neqns)
      flag = 0
   end subroutine
end module
