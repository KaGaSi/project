PROGRAM traject

!***********************************************************************
!
!     daresbury laboratory ccp5 program for post-processing HISTORY
!     files created from dl_meso_dpd to produce trajectory file
!
!     copyright stfc daresbury laboratory
!     author - m. a. seaton april 2012
!
!     xyz trajectory file modification (M. Lisal, March 2013)
!
!***********************************************************************

      IMPLICIT none
      INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND (15, 307)
      INTEGER, PARAMETER :: nuser=5,nrtout=9,ntraj=10

      CHARACTER(80) :: text
      CHARACTER(8), ALLOCATABLE :: namspe (:), nammol (:)
      CHARACTER(8) :: name1
      CHARACTER(4) :: chan

      INTEGER, ALLOCATABLE :: ltp (:), ltm (:), mole (:), bndtbl (:,:), beads (:), bonds (:), nspe (:)
      INTEGER :: chain, ispe, imol, ibond, imxspe, ioerror, i, numtraj, itime, j, k, first, last, nmoldef
      INTEGER :: maxnum, numspe, numnodes, nbeads, nubeads, numbeads, nbonds, numbonds, global, species, molecule
      INTEGER :: steps, lfrzn

      REAL(KIND=dp), ALLOCATABLE :: xxx (:), yyy (:), zzz (:), amass (:), bbb (:)
      REAL(KIND=dp) :: volm, dimx, dimy, dimz, ndx, ndy, ndz, time, mbeads, mglobal, x, y, z, vx, vy, vz

      LOGICAL :: volchange, eof, bigend

!     get number of nodes

      CALL get_command_argument (1, text)
      numnodes = parseint (text)

      IF (numnodes==0) THEN
        PRINT *,"Number of nodes used in calculations :"
        READ (nuser,*) numnodes
      END IF

      ALLOCATE (beads (numnodes), bonds (numnodes))

!     determine if HISTORY files exist

      IF (numnodes>1) THEN
        INQUIRE (file = 'HISTORY0000', EXIST = eof)
      ELSE
        INQUIRE (file = 'HISTORY', EXIST = eof)
      END IF

      IF (.NOT. eof) THEN
        PRINT *, "ERROR: cannot find HISTORY files"
        STOP
      END IF

!     determine endianness of HISTORY files

      IF (numnodes>1) THEN
!       OPEN (ntraj, file = 'HISTORY0000', access = 'sequential', form = 'unformatted', &
!                   &status = 'unknown', convert = 'big_endian')
        OPEN (ntraj, file = 'HISTORY0000', access = 'sequential', form = 'unformatted', & ! ML
                    &status = 'unknown', convert = 'native')
      ELSE
 !      OPEN (ntraj, file = 'HISTORY', access = 'sequential', form = 'unformatted', &
 !                  &status = 'unknown', convert = 'big_endian')
        OPEN (ntraj, file = 'HISTORY', access = 'sequential', form = 'unformatted', & ! ML
                    &status = 'unknown', convert = 'native')
      END IF

      READ (ntraj, IOSTAT=ioerror) numspe, nmoldef, nubeads, numbeads, nbeads, nbonds
      READ (ntraj, IOSTAT=ioerror) dimx, dimy, dimz, volm

      IF (ioerror/=0) THEN
        bigend = .false.
      ELSE
        bigend = .true.
      END IF

      CLOSE (ntraj)

!     scan HISTORY file for process 0 to determine species/molecule names
!     and properties, then find number of available timesteps

      IF (numnodes>1) THEN
        IF (bigend) THEN
!         OPEN (ntraj, file = 'HISTORY0000', access = 'sequential', form = 'unformatted', &
!                     &status = 'unknown', convert = 'big_endian')
          OPEN (ntraj, file = 'HISTORY0000', access = 'sequential', form = 'unformatted', & ! ML
                      &status = 'unknown', convert = 'native')
        ELSE
!         OPEN (ntraj, file = 'HISTORY0000', access = 'sequential', form = 'unformatted', &
!                     &status = 'unknown', convert = 'little_endian')
          OPEN (ntraj, file = 'HISTORY0000', access = 'sequential', form = 'unformatted', & ! ML
                      &status = 'unknown', convert = 'native')
        END IF
      ELSE
        IF (bigend) THEN
!         OPEN (ntraj, file = 'HISTORY', access = 'sequential', form = 'unformatted', &
!                     &status = 'unknown', convert = 'big_endian')
          OPEN (ntraj, file = 'HISTORY', access = 'sequential', form = 'unformatted', & ! ML
                      &status = 'unknown', convert = 'native')
        ELSE
!         OPEN (ntraj, file = 'HISTORY', access = 'sequential', form = 'unformatted', & ! ML
!                     &status = 'unknown', convert = 'little_endian')
          OPEN (ntraj, file = 'HISTORY', access = 'sequential', form = 'unformatted', &
                      &status = 'unknown', convert = 'native')
        END IF
      END IF

      READ (ntraj) numspe, nmoldef, nubeads, numbeads, nbeads, nbonds
      READ (ntraj) dimx, dimy, dimz, volm

      ALLOCATE (namspe (numspe), nammol (nmoldef))
      ALLOCATE (amass (numspe), bbb (numspe))

      DO i = 1, numspe
        READ (ntraj) namspe (i), amass (i), bbb (i), lfrzn
      END DO

      DO i = 1, nmoldef
        READ (ntraj) nammol (i)
      END DO

      READ (ntraj) text

      DO i = 1, nbeads
        READ (ntraj) global, species, molecule, chain
      END DO

      IF (nbonds>0) THEN
        DO i = 1, nbonds
          READ (ntraj) molecule, chain
        END DO
      END IF

      numtraj = 0

      DO WHILE (.true.)

        READ (ntraj, IOSTAT=ioerror) time, mbeads, ndx, ndy, ndz

        IF (ioerror/=0) THEN
          EXIT
        ELSE
          numtraj = numtraj + 1
          nbeads = NINT (mbeads)
          DO i = 1, nbeads
            READ (ntraj) mglobal, x, y, z, vx, vy, vz
          END DO
        END IF

      END DO

      CLOSE (ntraj)

      IF (numtraj==0) THEN
        PRINT *, 'ERROR: cannot find trajectory data in HISTORY files'
        STOP
      END IF

!     read in and initialize data, including numbers of beads and bonds

      numbonds = 0

      DO j = 1, numnodes

        WRITE (chan, '(i4.4)') j-1
        IF (numnodes>1) THEN
          IF (bigend) THEN
!           OPEN ((ntraj+j-1), file = 'HISTORY'//chan, access = 'sequential', form = 'unformatted', &
!                             &status = 'unknown', convert = 'big_endian')
            OPEN ((ntraj+j-1), file = 'HISTORY'//chan, access = 'sequential', form = 'unformatted', & ! ML
                              &status = 'unknown', convert = 'native')
          ELSE
 !          OPEN ((ntraj+j-1), file = 'HISTORY'//chan, access = 'sequential', form = 'unformatted', &
 !                            &status = 'unknown', convert = 'little_endian')
            OPEN ((ntraj+j-1), file = 'HISTORY'//chan, access = 'sequential', form = 'unformatted', & ! ML
                              &status = 'unknown', convert = 'native')
          END IF
        ELSE
          IF (bigend) THEN
!           OPEN ((ntraj+j-1), file = 'HISTORY', access = 'sequential', form = 'unformatted', &
!                             &status = 'unknown', convert = 'big_endian')
            OPEN ((ntraj+j-1), file = 'HISTORY', access = 'sequential', form = 'unformatted', & ! ML
                              &status = 'unknown', convert = 'native')
          ELSE
!           OPEN ((ntraj+j-1), file = 'HISTORY', access = 'sequential', form = 'unformatted', &
!                             &status = 'unknown', convert = 'little_endian')
            OPEN ((ntraj+j-1), file = 'HISTORY', access = 'sequential', form = 'unformatted', & ! ML
                              &status = 'unknown', convert = 'native')
          END IF
        END IF

        READ (ntraj+j-1) numspe, nmoldef, nubeads, numbeads, nbeads, nbonds
        READ (ntraj+j-1) dimx, dimy, dimz, volm

        beads (j) = nbeads
        bonds (j) = nbonds
        numbonds = numbonds + nbonds

      END DO

      ALLOCATE (bndtbl (numbonds, 2))

      DO j = 1, numnodes

        DO i = 1, numspe
          READ (ntraj+j-1) name1, vx, vy
        END DO

        DO i = 1, nmoldef
          READ (ntraj+j-1) name1
        END DO

        READ (ntraj+j-1) text

      END DO

!     get number of beads to include in trajectory file

      first = 1
      last = numbeads
      steps = numtraj

      ALLOCATE (xxx (first:last), yyy (first:last), zzz (first:last))
      ALLOCATE (ltp (1:numbeads), ltm (1:numbeads), mole (1:numbeads), nspe (numspe))

!     read in data for beads and bonds

      ibond = 0
      nspe = 0

      DO j = 1, numnodes

        nbeads = beads (j)
        nbonds = bonds (j)

        DO i = 1, nbeads
          READ (ntraj+j-1) global, species, molecule, chain
          IF (global>=first .AND. global<=last) THEN
            ltp (global) = species
            ltm (global) = molecule
            mole (global) = chain
            nspe (species) = nspe (species) + 1
          END IF
        END DO

        IF (nbonds>0) THEN
          DO i = 1, nbonds
            ibond = ibond + 1
            READ (ntraj+j-1) bndtbl (ibond, 1), bndtbl (ibond, 2)
          END DO
        END IF

      END DO

!     write output file (.VTF for VMD)

      OPEN (nrtout, file = 'dl_meso.vsf')

      maxnum = 0
      DO i = 1, numspe
        IF (nspe (i)>maxnum) THEN
          maxnum = nspe (i)
          imxspe = i
        END IF
      END DO

!     define beads

      IF ((last-first+1)==maxnum) THEN
        WRITE (nrtout, '("atom 0:",I8.8,"    radius ",F10.6," name ",A8)') last-first, bbb (imxspe), namspe (imxspe)
      ELSE
        WRITE (nrtout, '("atom default    radius ",F10.6," name ",A8)') bbb (imxspe), namspe (imxspe)
        DO i = first, last
          IF ((ltp (i)/=imxspe .OR. i==last) .AND. ltm(i)==0) THEN
            WRITE (nrtout, '("atom ",I8,"    radius ",F10.6," name ",A8)') i-first, bbb (ltp (i)), namspe (ltp (i))
          ELSE IF (ltm (i)/=0) THEN
            WRITE (nrtout, '("atom ",I8,"    radius ",F10.6," name ",A8," resid ",I8)') i-first, bbb (ltp (i)), &
                  &namspe (ltp (i)), mole (i)
          END IF
        END DO
      END IF

!     define bonds

      IF (numbonds>0 .AND. last>nubeads) THEN
        WRITE (nrtout, '()')
        DO i = 1, numbonds
          IF ((bndtbl (i, 1)-first)>=0) WRITE (nrtout, '("bond ",I8,":",I8)') (bndtbl (i, 1)-first), (bndtbl (i, 2)-first)
        END DO
      END IF

      CLOSE (nrtout)

!     obtain positions and velocities for all beads at each time step

      volchange = .true.
      eof = .false.

      OPEN (nrtout, file = 'All.vcf')

      WRITE (nrtout, '("pbc ", 3F10.6)') dimx, dimy, dimz

      DO k = 1, steps

        DO j = 1, numnodes

          READ (ntraj+j-1, IOSTAT=ioerror) time, mbeads, ndx, ndy, ndz

          IF (ioerror/=0) THEN
            PRINT *, 'ERROR: End of file reached prematurely - ', k-1, ' timesteps written to VTF file'
            eof = .true.
            EXIT
          END IF

          nbeads = NINT (mbeads)
          volchange = (ndx/=dimx .OR. ndy/=dimy .OR. ndz/=dimz)
          dimx = ndx
          dimy = ndy
          dimz = ndz

          DO i = 1, nbeads

            READ (ntraj+j-1) mglobal, x, y, z, vx, vy, vz
            global = NINT (mglobal)

            IF (global>=first .AND. global <=last) THEN
              xxx (global) = x
              yyy (global) = y
              zzz (global) = z
            END IF

          END DO

          IF (eof) EXIT

        END DO

        IF (eof) EXIT

        IF (volchange .OR. k==1) THEN
          WRITE (nrtout, '(/,"#",(I6),/,"timestep")') k
        ELSE
          WRITE (nrtout, '(/,"#",(I6),/,"timestep")') k
        END IF

        DO i = first, last
          WRITE (nrtout, '(3F8.3)') xxx (i), yyy (i), zzz (i)
        END DO

        PRINT *, "Timestep ",k,"/",numtraj," : time = ",time

      END DO

      CLOSE (nrtout)

      DO j = 1, numnodes
        CLOSE (ntraj+j-1)
      END DO

CONTAINS

      INTEGER FUNCTION parseint (word)

!*********************************************************************
!
!     copyright - daresbury laboratory
!     author    - w.smith 2001
!
!*********************************************************************

      IMPLICIT none

      CHARACTER(1) :: u
      CHARACTER(16) :: word
      INTEGER :: i, m, k, s

      k = 0
      s = 1
      m = LEN (word)

      DO i = 1, m

        u = word (i:i)

        IF (u=="-") THEN

          s = -1

        ELSE IF (u=="0") THEN

          k = 10 * k

        ELSE IF (u=="1") THEN

          k = 10 * k + 1

        ELSE IF (u=="2") THEN

          k = 10 * k + 2

        ELSE IF (u=="3") THEN

          k = 10 * k + 3

        ELSE IF (u=="4") THEN

          k = 10 * k + 4

        ELSE IF (u=="5") THEN

          k = 10 * k + 5

        ELSE IF (u=="6") THEN

          k = 10 * k + 6

        ELSE IF (u=="7") THEN

          k = 10 * k + 7

        ELSE IF (u=="8") THEN

          k = 10 * k + 8

        ELSE IF (u=="9") THEN

          k = 10 * k + 9

        END IF

      END DO

      parseint = s * k

      RETURN
      END FUNCTION parseint

END PROGRAM traject
