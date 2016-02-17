PROGRAM traject

!***********************************************************************
!
!     daresbury laboratory ccp5 program for post-processing HISTORY
!     files created from dl_meso_dpd to produce trajectory file
!
!     copyright stfc daresbury laboratory
!     author - m. a. seaton october 2015
!
!***********************************************************************

      IMPLICIT none
      INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND (15, 307)
      INTEGER, PARAMETER :: nuser=5,nrtout=9,ntraj=11

      CHARACTER(80) :: text
      CHARACTER(8), ALLOCATABLE :: namspe (:)
      CHARACTER(8) :: name1, nammol
      CHARACTER(6) :: chan

      INTEGER, ALLOCATABLE :: ltp (:), ltm (:), mole (:), bndtbl (:,:), beads (:), bonds (:), nspe (:)
      INTEGER :: chain, ibond, imxspe, ioerror, i, j, k, nmoldef
      INTEGER :: maxnum, numspe, numnodes, nbeads, nubeads, numbeads, nbonds, numbonds, global, species, molecule
      INTEGER :: lfrzn, keytrj, srfx, srfy, srfz

      REAL(KIND=dp), ALLOCATABLE :: bbb (:)
      REAL(KIND=dp) :: volm, dimx, dimy, dimz, halfx, halfy, halfz, shrdx, shrdy, shrdz, amass
      REAL(KIND=dp) :: ndx, ndy, ndz, time, mbeads, mglobal, x, y, z, vx, vy, vz, fx, fy, fz

      LOGICAL :: volchange, eof

!     get number of nodes

      CALL get_command_argument (1, text)
      numnodes = parseint (text)

      IF (numnodes==0) THEN
        PRINT *,"Number of nodes used in calculations :"
        READ (nuser,*) numnodes
      END IF

      ALLOCATE (beads (numnodes), bonds (numnodes))

!     CALL get_command_argument (2, text)

!     determine if HISTORY files exist

      IF (numnodes>1) THEN
        INQUIRE (file = 'HISTORY000000', EXIST = eof)
      ELSE
        INQUIRE (file = 'HISTORY', EXIST = eof)
      END IF

      IF (.NOT. eof) THEN
        PRINT *, "ERROR: cannot find HISTORY files"
        STOP
      END IF

!     scan HISTORY file for process 0 to determine species/molecule names and properties

      IF (numnodes>1) THEN
        OPEN (ntraj, file = 'HISTORY000000', access = 'sequential', form = 'unformatted', status = 'unknown')
      ELSE
        OPEN (ntraj, file = 'HISTORY', access = 'sequential', form = 'unformatted', status = 'unknown')
      END IF

      READ (ntraj) numspe, nmoldef, nubeads, numbeads, nbeads, nbonds
      READ (ntraj) dimx, dimy, dimz, volm
      READ (ntraj) keytrj, srfx, srfy, srfz

      ALLOCATE (namspe (numspe))
      ALLOCATE (bbb (numspe))

      DO i = 1, numspe
        READ (ntraj) namspe (i), amass, bbb (i), lfrzn
      END DO

      DO i = 1, nmoldef
        READ (ntraj) nammol
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

      CLOSE (ntraj)

!     read in and initialize data, including numbers of beads and bonds

      numbonds = 0

      DO j = 1, numnodes

        WRITE (chan, '(i6.6)') j-1
        IF (numnodes>1) THEN
          OPEN ((ntraj+j-1), file = 'HISTORY'//chan, access = 'sequential', form = 'unformatted', status = 'unknown')
        ELSE
          OPEN (ntraj, file = 'HISTORY', access = 'sequential', form = 'unformatted', status = 'unknown')
        END IF

        READ (ntraj+j-1) numspe, nmoldef, nubeads, numbeads, nbeads, nbonds
        READ (ntraj+j-1) dimx, dimy, dimz, volm
        READ (ntraj+j-1) keytrj, srfx, srfy, srfz

        beads (j) = nbeads
        bonds (j) = nbonds
        numbonds = numbonds + nbonds

      END DO

      ALLOCATE (bndtbl (numbonds, 2))

      DO j = 1, numnodes

        DO i = 1, numspe
          READ (ntraj+j-1) name1, vx, vy, lfrzn
        END DO

        DO i = 1, nmoldef
          READ (ntraj+j-1) name1
        END DO

        READ (ntraj+j-1) text

      END DO

!     get number of beads to include in trajectory file

      ALLOCATE (ltp (1:numbeads), ltm (1:numbeads), mole (1:numbeads), nspe (numspe))

!     read in data for beads and bonds

      ibond = 0
      nspe = 0

      DO j = 1, numnodes

        nbeads = beads (j)
        nbonds = bonds (j)

        DO i = 1, nbeads
          READ (ntraj+j-1) global, species, molecule, chain
          ltp (global) = species
          ltm (global) = molecule
          mole (global) = chain
          nspe (species) = nspe (species) + 1
        END DO

        IF (nbonds>0) THEN
          DO i = 1, nbonds
            ibond = ibond + 1
            READ (ntraj+j-1) bndtbl (ibond, 1), bndtbl (ibond, 2)
          END DO
        END IF

      END DO

      DEALLOCATE (beads, bonds)

!     write output file (.VTF for VMD)

!     OPEN (nrtout, file = 'traject.vtf')
      OPEN (nrtout, file = 'dl_meso.vsf')

      maxnum = 0
      DO i = 1, numspe
        IF (nspe (i)>maxnum) THEN
          maxnum = nspe (i)
          imxspe = i
        END IF
      END DO

!     define beads

      IF (numbeads==maxnum .AND. numbonds==0) THEN
        WRITE (nrtout, '("atom 0:",I8.8,"    radius ",F10.6," name ",A8)') numbeads-1, bbb (imxspe), namspe (imxspe)
      ELSE
        WRITE (nrtout, '("atom default    radius ",F10.6," name ",A8)') bbb (imxspe), namspe (imxspe)
        DO i = 1, numbeads
          IF ((ltp (i)/=imxspe .OR. i==numbeads) .AND. ltm(i)==0) THEN
            WRITE (nrtout, '("atom ",I8,"    radius ",F10.6," name ",A8)') i-1, bbb (ltp (i)), namspe (ltp (i))
          ELSE IF (ltm (i)/=0) THEN
            WRITE (nrtout, '("atom ",I8,"    radius ",F10.6," name ",A8," resid ",I8)') i-1, bbb (ltp (i)), &
                  &namspe (ltp (i)), mole (i)
          END IF
        END DO
      END IF

!     define bonds

      IF (numbonds>0 .AND. numbeads>nubeads) THEN
        WRITE (nrtout, '()')
        DO i = 1, numbonds
          IF ((bndtbl (i, 1)-1)>=0) WRITE (nrtout, '("bond ",I8,":",I8)') (bndtbl (i, 1)-1), (bndtbl (i, 2)-1)
        END DO
      END IF

      CLOSE (nrtout)

      DEALLOCATE (ltp, ltm, mole, bndtbl, nspe, namspe, bbb)

      OPEN (nrtout, file = 'All.vcf')

!     obtain positions and velocities for all beads at each time step

      volchange = .true.
      eof = .false.
      k = 0

      DO WHILE (.true.)

        READ (ntraj, IOSTAT=ioerror) time, mbeads, ndx, ndy, ndz, shrdx, shrdy, shrdz

        IF (ioerror/=0) THEN
          eof = .true.
          IF (k==0) THEN
            PRINT *, 'ERROR: cannot find trajectory data in HISTORY files'
            STOP
          END IF
          EXIT
        END IF

        k = k + 1
        volchange = (ndx/=dimx .OR. ndy/=dimy .OR. ndz/=dimz)
        dimx = ndx
        dimy = ndy
        dimz = ndz
!        halfx = 0.5_dp * dimx
!        halfy = 0.5_dp * dimy
!        halfz = 0.5_dp * dimz
        halfx = 0.0_dp
        halfy = 0.0_dp
        halfz = 0.0_dp

        IF (k==1) THEN
          WRITE (nrtout, '("pbc ", 3F10.6,//,"#"(I10),/,"timestep")') dimx, dimy, dimz, k
        ELSE
          WRITE (nrtout, '(/,"#",(I6),/,"timestep")') k
        END IF

        DO j = 1, numnodes

          IF (j>1) THEN
            READ (ntraj+j-1, IOSTAT=ioerror) time, mbeads, ndx, ndy, ndz, shrdx, shrdy, shrdz
            IF (ioerror/=0) THEN
              PRINT *, 'ERROR: End of file reached prematurely - ', k-1, ' timesteps written to VTF file'
              eof = .true.
              EXIT
            END IF
          END IF

          IF (eof) EXIT

          nbeads = NINT (mbeads)

          SELECT CASE (keytrj)
          CASE (0)
            DO i = 1, nbeads
              READ (ntraj+j-1) mglobal, x, y, z
              global = NINT (mglobal)
              WRITE (nrtout, '(3F7.3)')  x-halfx, y-halfy, z-halfz
            END DO
          CASE (1)
            DO i = 1, nbeads
              READ (ntraj+j-1) mglobal, x, y, z, vx, vy, vz
              global = NINT (mglobal)
              WRITE (nrtout, '(3F7.3)')  x-halfx, y-halfy, z-halfz
            END DO
          CASE (2)
            DO i = 1, nbeads
              READ (ntraj+j-1) mglobal, x, y, z, vx, vy, vz, fx, fy, fz
              global = NINT (mglobal)
              WRITE (nrtout, '(3F7.3)')  x-halfx, y-halfy, z-halfz
            END DO
          END SELECT

        END DO

        IF (eof) EXIT

        PRINT *, "Timestep ",k," : time = ",time

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
