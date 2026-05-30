! MT2D package for OCCAM 3.0
! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
!
! Revision History:
!     Oct 2006 David Myer - move to F90, allocated memory, cmdline params,
!           modules, and general code reorg.  From mathematical perspective
!           this version includes options for model value limits and
!           diagonal penalties.  PW2D derivatives now returned with respect
!           to linear conductivity.                     (Version 3.0)
!     Dec-13-2005, Kerry Key
!           Modified to use different regularization options :
!     July 25, 2001 -  (Catherine deGroot-Hedlin)
!           alteration to fordiv, PW's new adjoint method,
!           which returns derviatives with respect to log(resistivity)
!     Feb 23, 1994  revision to allow for impedance data - cdh -
!     January 1993  Steven Constable's re-write of Catherine deGroot-Hedlin's
!           thesis code under Geotools sponsorship     (Version 1.01) 
!
! REFERENCES: CONSTABLE, PARKER & CONSTABLE, 1987: GEOPHYSICS 52, 289-300.
!             DEGROOT-HEDLIN & CONSTABLE, 1990: GEOPHYSICS 55, 1613-1624.
!             CONSTABLE, 1991: GEOPHYS. J. INT. 106, 387-388.
!             CONSTABLE, 1992: OCCAM DISTRIBUTION NOTES.
!             Myer et al., 2007: OCCAM 3.0 release notes.
!
! See http://marineemlab.ucsd.edu/Projects/Occam/
!
!-----------------------------------------------------------------------
SUBROUTINE OBJMAT(IRUF, model, DEL, nPTerms, bHaveLink, linkpr)
!-----------------------------------------------------------------------

! OCCAM MT2D 3.0 Package
! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
! Subroutine Revision 3.0, November 2006
!  Dec-13-2005, Kerry Key
!       Modified to use different regularization options 3-2XX:
! Subroutine Revision 1.01, 13 Jan 1993
!
! CONSTRUCTS PENALTY MATRIX (DEL) FOR 2DMT CASE.
! WE USE A SIMPLE PENALTY MEASURE FOR THE TIME BEING WHICH ASSUMES BLOCK
! HEIGHT INCREASES EXPONENTIALLY WITH DEPTH.  LATER WE CAN USE ACTUAL 
! DIMENSIONS OF THE COMMON EDGES OF THE BRICKS TO WORK OUT THE PENALTY 
! MEASURE.

! ON INPUT:

!  IRUF valid values:
!           1    1ST DERIVATIVE ROUGHNESS
!           2    2nd Deriv is NOT supported in the 2D occam
!           3    1st Deriv + PARAMETER ASPECT RATIO WEIGHTING
!           4    1st Deriv + DEPTH WEIGHTED SMOOTHING (LOG10(DEPTH))
!           1XX  HORIZONTALLY BIASED, XX IS FACTOR WEIGHT APPLIED 
!                TO HORIZ ROUGHNESS, VERT REMAINS AS FIRST DIFFERENCES 
!           2XX  VERTICALLY BIASED, XX IS FACTOR WEIGHT APPLIED 
!                TO VERT ROUGHNESS, HORIZ REMAINS AS FIRST DIFFERENCES
!   nPTerms = (INTEGER) NUMBER OF NON-ZERO ROWS IN PENALTY MATRIX
! ON OUTPUT:
!   DEL() = (REAL) PENALTY MATRIX
!   linkpr = brick #s of pairs of bricks linked in the penalty matrix.
!           Used in deltdel to calculate (DEL)Transpose * DEL using
!           essentially sparse matrix methods.
!-----------------------------------------------------------------------
    USE Occam_Mod
    USE Fwd_mod

    IMPLICIT NONE
    
    ! ARGUMENTS:
    INTEGER, intent(in)     :: IRUF
    real, intent(in)        :: model(nParams)
    REAL, intent(out)       :: DEL(nPTerms,nParams)
    INTEGER, intent(in)     :: nPTerms    ! # of penalty terms (for dim of DEL & linkpr)
    logical, intent(in)     :: bHaveLink  ! true if linkpr was passed.  NB: intrinsic fctn present() doesn't appear to work.
    INTEGER, OPTIONAL, intent(out) :: linkpr(nPTerms,2)        ! DGM Oct 2006 - not always needed.  Quite large!
    
    ! LOCAL VARIABLES
    INTEGER I, J,  K, nPen
    logical :: bChkPair, bExceptUsed(nExep)
    REAL PENALT
    !   bExceptUsed() = A TABLE TO INDICATE WHETHER PENALTY 'EXCEPTIONS' ARE MOD-
    !      IFICATIONS OF THE DEFAULT PENALTIES OR ADDITIONS TO THE PENALTY
    !      MATRIX
    !   PENALT = WEIGHTING OF PENALTY

!------------------------------------------------
    ! Init working variables
    DEL = 0.0     ! Array math
    nPen = 0
    PENALT = 0    ! Oct 2006 - was not initialized - if IRUF 0 or 2, wld be random.
    bExceptUsed = .false.
    
    ! Count the horizontal & vertical penalty blocks
    do i = 1,NBRICK
        ! Look for blocks connected to block i
        ! NB: ICORN's 2nd index: 1=top, 2=left, 3=bottom, 4=right
        do j = 1,NBRICK
            bChkPair = .false.
            
            ! Horizontal or Diagonal
            if (ICORN(i,4) == ICORN(j,2)) then
                if (ICORN(i,1) == ICORN(j,1)) then          ! Horizontal
                    nPen = nPen + 1
                    bChkPair = .true.
                    
                    IF (IRUF .EQ. 1) then                           ! FIRST DIFFERENCE ONLY
                       PENALT = 1.0
                    ELSEIF (IRUF .EQ. 3) then                       ! ASPECT RATIO WEIGHTING 
                       PENALT = 2.*THICK(I)/(WIDTH(I)+WIDTH(J))  
                    ELSEIF (IRUF .EQ. 4) then                       ! DEPTH WEIGHTED
                       PENALT = log10( sum( DLZ(1:ICORN(I,1)) ) )
                    elseif (IRUF >= 100 .and. IRUF <= 199) then     ! 1XX HORIZONTALLY BIASED
                       PENALT = MOD(IRUF,100)
                    elseif (IRUF >= 200 .and. IRUF <= 299) then     ! 2XX VERTICALLY BIASED
                       PENALT = 1.0
                    ENDIF
                    
                elseif (gbAddDiags) then      ! if we care about diagonals
                    if ( ICORN(i,1) == ICORN(j,3) &         ! diag up
                    .or. ICORN(i,3) == ICORN(j,1) ) then    ! diag down
                        bChkPair = .true.
                        nPen = nPen + 1

                        IF (IRUF .EQ. 1) then                           ! first difference only:
                           PENALT = 1.0
                        ELSEIF (IRUF .EQ. 3) then                       ! Aspect ratio weighting
                           PENALT = 2.*THICK(I)/(WIDTH(I)+WIDTH(J))  
                        ELSEIF (IRUF .EQ. 4) then                       ! DEPTH WEIGHTED
                           PENALT = log10( sum( DLZ(1:ICORN(I,1)) ) )
                        elseif (IRUF >= 100 .and. IRUF <= 199) then     ! 1XX HORIZONTALLY BIASED
                           PENALT = 1.0
                        elseif (IRUF >= 200 .and. IRUF <= 299) then     ! 2XX VERTICALLY BIASED
                           PENALT = 1.0
                        ENDIF
                    endif
                endif
            elseif ( (ICORN(i,3) == ICORN(j,1)) &           ! Vertical?
              .AND. &
            ((ICORN(i,2) == ICORN(j,2)).OR.(ICORN(i,4) == ICORN(j,4))) &
                ) then
                bChkPair = .true.
                nPen = nPen + 1

                IF (IRUF .EQ. 1) then                           ! first difference only:
                   PENALT = 1.0
                ELSEIF (IRUF .EQ. 3) then                       ! Aspect ratio weighting
                   PENALT = 2.*THICK(I)/(WIDTH(I)+WIDTH(J))  
                ELSEIF (IRUF .EQ. 4) then                       ! DEPTH WEIGHTED
                   PENALT = log10( sum( DLZ(1:ICORN(I,1)) ) )
                elseif (IRUF >= 100 .and. IRUF <= 199) then     ! 1XX HORIZONTALLY BIASED
                   PENALT = 1.0
                elseif (IRUF >= 200 .and. IRUF <= 299) then     ! 2XX VERTICALLY BIASED
                   PENALT = MOD(IRUF,200)
                ENDIF
            endif
            
            ! If we found a linkage between blocks i and j...
            if (bChkPair) then
                ! Is an exception applied to this pair?  If so, modify the penalty
                do k = 1, NEXEP
                    if ( ((IJEXEP(k,1) == i) .AND. (IJEXEP(k,2) == j)) &
                    .OR. &
                         ((IJEXEP(k,1) == j) .AND. (IJEXEP(k,2) == i)) &
                        ) then
                        bExceptUsed(k) = .true.
                        PENALT = PENALT*ABS(EXPEN(K))
                        exit
                    endif
                enddo
                
                ! Record the penalty in the DEL matrix & the relationship in linkpr
                DEL(nPen, I) = -1. * PENALT
                DEL(nPen, J) =  1. * PENALT
                if (bHaveLink) then     ! didn't work: (present(linkpr)) then
                    linkpr(nPen,1)=i
                    linkpr(nPen,2)=j
                endif
            endif
            
        enddo
    enddo
    
    ! Any exceptions in the array that were not used to modify bricks
    ! that would normally get penalties are used to link bricks that
    ! are not normally linked.  Add the linkages to the DEL array now.
    DO I = 1, NEXEP
        IF (.not. bExceptUsed(I)) THEN
            nPen = nPen + 1
            PENALT = ABS(EXPEN(I))
            DEL(nPen, IJEXEP(I,1)) = -1.*PENALT
            DEL(nPen, IJEXEP(I,2)) = 1.*PENALT
            if (bHaveLink) then     ! didn't work: (present(linkpr)) then
                linkpr(nPen,1) = ijexep(i,1)
                linkpr(nPen,2) = ijexep(i,2)
            endif
        END IF
    enddo

    ! SET THE PREJUDICE WEIGHTS (i.e. block on itself)
    ! Change the weighting so that as the blocks get taller but
    ! no wider, structures are not abnormally elongated.
    PREWTS(1:NBRICK) = PREWT(1:NBRICK) !kwk * WIDTH(1:NBRICK) / THICK(1:NBRICK)
      
    ! Copy over the weights for the free static shifts too.
    PREWTS(NBRICK+1:nParams) = PREWT(NBRICK+1:nParams)

    return
end subroutine objmat

!-----------------------------------------------------------------------
INTEGER FUNCTION CountPenaltyTerms( iRuf ) result(nTerms)
!-----------------------------------------------------------------------
! OCCAM MT2D 3.0 Package
! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
! Subroutine Revision 3.0, November 2006
! DGM Oct 2006 - # of penalty terms may change depending on whether the
!   roughness type requested does just horiz & vert or also diagonal.
    use Occam_mod
    use Fwd_mod
    implicit none
    
    ! Parameters
    integer, intent(in)     :: iRuf
    
    ! Local vars
    integer :: i, j, k
    logical :: bChkPair, bExceptUsed(nExep)
    
    ! Init various working variables
    nTerms = 0
    bExceptUsed = .false.
    
    ! Count the horizontal & vertical penalty blocks
    do i = 1,NBRICK
        ! Look for blocks connected to block i
        ! NB: ICORN's 2nd index: 1=top, 2=left, 3=bottom, 4=right
        do j = 1,NBRICK
            bChkPair = .false.
            
            ! Horizontal or Diagonal
            if (ICORN(i,4) == ICORN(j,2)) then
                if (ICORN(i,1) == ICORN(j,1)) then          ! Horizontal
                    nTerms = nTerms + 1
                    bChkPair = .true.
                elseif (gbAddDiags) then      ! if we care about diagonals
                    if ( ICORN(i,1) == ICORN(j,3) &         ! diag up
                    .or. ICORN(i,3) == ICORN(j,1) ) then    ! diag down
                        nTerms = nTerms + 1
                        bChkPair = .true.
                    endif
                endif
            elseif ( (ICORN(i,3) == ICORN(j,1)) &           ! Vertical?
              .AND. &
            ((ICORN(i,2) == ICORN(j,2)).OR.(ICORN(i,4) == ICORN(j,4))) &
                ) then
                nTerms = nTerms + 1
                bChkPair = .true.
            endif
                
            ! Is an exception applied to this pair?
            if (bChkPair) then
                do k = 1, NEXEP
                    if ( ((IJEXEP(k,1) == i) .AND. (IJEXEP(k,2) == j)) &
                    .OR. &
                         ((IJEXEP(k,1) == j) .AND. (IJEXEP(k,2) == i)) &
                        ) then
                        bExceptUsed(k) = .true.
                        exit
                    endif
                enddo
            endif
            
        enddo
    enddo
    
    ! Count exceptions that were not applied above.  This allows you
    ! to use the exceptions to create penalties between blocks that
    ! otherwise would not have penalties (i.e. disconnected blocks)
    do i = 1, NEXEP
        if (.not. bExceptUsed(i)) then
            nTerms = nTerms + 1
        endif
    enddo

    return
end FUNCTION CountPenaltyTerms     



!-----------------------------------------------------------------------
      subroutine inputd(datfil)
!-----------------------------------------------------------------------
!
! OCCAM MT2D 3.0 Package
! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
! Subroutine Revision 3.0, November 2006
! revision to allow for Zrot value of input data - cdh -March 29, 1999
! revision to allow for impedance data - cdh -February 23, 1994
! Subroutine Revision 1.02, 13 Jan 1993
!
! inputd opens the data file datfil and reads in the data array
!
! includes:
      USE Occam_Mod
      USE Fwd_mod
      IMPLICIT NONE

! arguments:
      character(*) datfil
      
! local variables
      logical iffte, ifftm
!   flags to indicate which modes are required
      integer lerr, i, j, nAllocErr
!   error flag, loop index and handy integer
      character(80)  string
!   string buffer
      integer idcode
!   idcode= function that performs the equivalent of a free-format internal
!            read for an integer
!
! parameters:
      integer iof 
      parameter (iof = 15)
!   iof = io unit number for file operations
!----------------------------------------------------------
    ! open data file 
      open (iof, file=datfil, status='old', iostat=lerr)
      if (lerr .ne. 0) then
        write(*,*) ' Error opening data file'
        stop
      end if
      
    ! read stuff
      read(iof,'(18x,a)', end=198, err=199) string
      if (string(1:16) .ne. 'OCCAM2MTDATA_1.0') call filerr( &
    &   ' Data file not supported',iof)
      read(iof,'(18x,a)', end=198, err=199) string
      
      ! # of receivers
      read(iof,'(18x,a)', end=198, err=199) string
      nrc = idcode(string, iof)
      write(*,*) '# Receivers:  ', nrc

      ! Allocate things dependent on # rcvrs
      allocate( nrl(nrc)            &
              , sitloc(nrc)         &
              , isx1(2*nrc)         &   ! Static shift stuff
              , isx3(2*nrc)         &
              , isfx(2*nrc)         &
              , shift(2*nrc)        &
              , cSiteName(nrc)      & 
              , stat=nAllocErr )
      if (nAllocErr .ne. 0) then
        write(*,*) 'Out of memory.  Try reducing the model size.'
        stop
      endif

      ! Read the receiver names.
      do i = 1, nrc
         read(iof,'(a)', end=198, err=199) cSiteName(i)
      enddo
      read(iof,'(a)', end=198, err=199) string
      
      ! Rcvr locations
      read(iof,*, end=198, err=199) (sitloc(i), i=1,nrc)
      
      ! # of Frequencies
      read(iof,'(18x,a)', end=198, err=199) string
      nfre = idcode(string, iof)

      ! Allocate freq related vars
      allocate( frek(nfre)          &
              , stat=nAllocErr )
      if (nAllocErr .ne. 0) then
        write(*,*) 'Out of memory.  Try reducing the model size'
        stop
      endif
      
      ! Read frequency list
      read(iof,*, end=198, err=199) (frek(i), i=1,nfre)
      
      ! # of data
      read(iof,'(18x,a)', end=198, err=199) string
      ndata = idcode(string, iof)

      ! Allocate vars related to # of data
      allocate( d(ndata)                &
              , sd(ndata)               &
              , dp(ndata,cnDataParams)  &
              , stat=nAllocErr )
      if (nAllocErr .ne. 0) then
        write(*,*) 'Out of memory.  Try reducing the model size'
        stop 
      endif
      
      ! Read the data
      read(iof,'(a)', end=198, err=199) string  ! header line
      read(iof,*, end=198, err=199) &
        (dp(i,1), dp(i,2), dp(i,3), d(i), sd(i), i=1,ndata)
      dp(:,4) = 0   ! init the extra column for rotation (not used at the moment)
      
      ! set idx, the mode definition for PW2DI
      iffte = .false.
      ifftm = .false.
      do 20 i = 1, ndata
        j = int(dp(i,3)+0.01)
        select case (j)
        case (1, 2, 3, 4, 9, 13, 14)      ! See Fwd_mod for list of valid values
            iffte = .true.
        case (5, 6, 7, 8, 10, 15, 16) 
            ifftm = .true.
        case (11, 12, 17, 18)
            iffte = .true.
            ifftm = .true.
        case default
            write(*,*) ' Invalid data type: ', j
            stop
        end select
20      continue

      if ((iffte .and. ifftm)) then
        idx = 1
      else if (iffte) then
        idx = 2
      else if (ifftm) then
        idx = 3
      end if

      if (idebug .ge. 1) then
        write(*,*) '# Data:       ', ndata
      end if

      ! initially set nd for occam as number of MT data
      !   May be modified later to account for static shift parameters
      nd = ndata

      return

198      call filerr (' Data file ended prematurely', iof)
199      call filerr (' Error reading data file', iof)
!
      end
!     
!-----------------------------------------------------------------------
      subroutine inputm(modfil)
!----------------------------------------------------------------------
!
! OCCAM MT2D 3.0 Package
! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
! Subroutine Revision 3.0, November 2006
! Subroutine Revision 1.04, 14 Jan 1993
!
! inputm() opens and reads the model file modfil, opens and reads the 
! mesh and static shift files pointed at by modfil, and sets up the 
! environments needed by OCCAM and PW2DI.
!
! Note: inputm() assumes that inputd() has already been called
!
! on input:
!   modfil = (character) model file name 
! on output:
!   incredible amounts of stuff in common blocks (see include files)
! calls:
!   filerr
!
      USE Occam_Mod
      USE Fwd_mod
      IMPLICIT NONE

!
! arguments:
      character(*) modfil
!      
! local variables:
      INTEGER, ALLOCATABLE:: IRZ(:), IRY(:,:), IRINDX(:,:)
      INTEGER, ALLOCATABLE:: LOCEXP(:,:), ICONB(:)
      
      integer lerr, i, l, j
      integer nsites, lidx, lfre, nexec, nAllocErr
      integer ii, jj, inz1, inz2, iny1, iny2
      integer nprej, nColMax
!   irz() = no of vert. FEs in each row
!   iry() = no of horiz FEs in each col for each row
!   irindx() = layer/column to parameter pointers.  0 if fixed block
!   lerr = iostat from file operations
!   i,l,j = loop counters
!   nsites = no of sites read from file.  Ignored after input by OCCAM
!   lidx = type of modelling; should be 0 for OCCAM
!   lfre = number of frequencies read from file. Ignored by OCCAM
!   nexec = type of action required - not used by inversion process.
!   ii, jj = index counters for regularization bricks
!   inz1, inz2 = top and  bottom edges of bricks defined by FE mesh
!   iny1, iny2 = left and right edges     ..    ditto          ..
!   locexp() = brick layer/col indices of exeption pairs 
!   nprej = number of entries in model prejudice 
!   iconb() = flag to indicate nodes on conductivity boundaries
      real dummy, tesum, tmsum, sserr
      real temp1,temp2
!   dummy = dummy variable to read frequencies into
!   tesum = sum of logs of TE free static shifts 
!   tmsum = sum of logs of TM free static shifts 
!   sserr = error term to use in conjuction with static shift constraint

      character(1), allocatable :: apt(:,:,:)
!   apt() = character array to store 4 resistivity codes for each element

      character(80) mshfil, mshtyp, string, ssfile, prfile
!   mshfil = mesh file name
!   mshtyp = mesh type
!   string = string buffer for free format reads
!   ssfile = static shift file name
!   prfile = file containing model prejudice and weights
!
      logical iffree
!   iffree = flag to indicate free element in regularization brick
      integer idcode
!   idcode = function to perform the equivalent of a free-format internal
!            read for an integer
      real rdcode
!   rdcode = function to perform the equivalent of a free-format internal
!            read for a real
!
! parameters:
      integer iof 
      parameter (iof = 15)
!   iof = io unit number for file operations
!
!-----------------------------------------------------------------------      

!
! Allocate arrays needed:
!
    allocate( thick(nParams)          &
            , width(nParams)          &
            , y(nParams+cnMaxFixedRes) &   ! + for all possible fixed resistivities (0-9,A-Z)
            , icorn(nParams,4)        &   ! 4 sides: 1=top, 2=left, 3=bottom, 4=right 
            , laprm(nParams)          &
            , prewt(nParams)          &      
            , stat=nAllocErr )
            
    if (nAllocErr .ne. 0) then
        write(*,*) 'Out of memory.  Too many free parameters (', nParams, ')'
        stop 
    endif           



!ccccccccccccccccccc model file reading cccccccccccccccccccccccccccccccc
!
      open (iof, file=modfil, status='old', iostat=lerr)
      if (lerr .ne. 0) then
        write(*,*) ' Error opening model file'
        stop 
      end if
      
    ! read character descriptors      
      read (iof,'(18x,a)', end=198, err=199) string
      if (string(1:15) .ne. 'OCCAM2MTMOD_1.0') call filerr( &
     &  ' Model file not supported',iof)
      read (iof,'(18x,a)', end=198, err=199) string
      read (iof,'(18x,a)', end=198, err=199) string
      read (iof,'(18x,a)', end=198, err=199) mshfil
      read (iof,'(18x,a)', end=198, err=199) mshtyp
      read (iof,'(18x,a)', end=198, err=199) ssfile
      read (iof,'(18x,a)', end=198, err=199) prfile
      
    ! read in the binding offset and number of layers
      read (iof,'(18x,a)', end=198, err=199) string
      boffs = rdcode(string, iof)
      read (iof,'(18x,a)', end=198, err=199) string
      nlay = idcode(string, iof)
      
      ! Allocate vars related to the # of layers
      Allocate( irz(nlay)           &
              , nrcol(nlay)         &
              , stat=nAllocErr )
      if (nAllocErr .ne. 0) then
        write(*,*) 'Out of memory.  Try reducing the model size.'
        stop
      endif
      
      ! Don't know how many cols there are.  So pre-read through the layers
      ! and keep the max #.  Then allocate local vars to hold the data, and
      ! re-read through the layer data.
      do i = 1,nlay
        read (iof,*,end=198,err=199) irz(i), nrcol(i)
        read (iof,*,end=198,err=199) (ii, j=1,nrcol(i)) ! discard these for now
      enddo
      close(iof)
      open (iof, file=modfil, status='old', iostat=lerr)
      if (lerr .ne. 0) then
        write(*,*) ' Error opening model file (2nd pass)'
        stop 
      end if
      do i = 1,9
          read (iof,*,end=198,err=199) string
      enddo
      
      nColMax = maxval(nrcol)
      Allocate( iry(nlay,nColMax), irindx(nlay,nColMax) &
              , stat=nAllocErr )
      if (nAllocErr .ne. 0) then
        write(*,*) 'Out of memory.  Try reducing the model size.'
        stop
      endif
      
    ! read in regularization mesh description
      do i = 1,nlay
          read (iof,*, end=198, err=199) irz(i),nrcol(i)
          read (iof,*, end=198, err=199) (iry(i,j), j=1,nrcol(i))
      enddo
      
    ! read in exceptions to penalty matrix
      read (iof,'(18x,a)', end=198, err=199) string
      nexep = idcode(string, iof)

    ! Allocate vars related to the # of exceptions
      Allocate( locexp(abs(nexep),4)     &   ! (local) exception bricks: row1,col1, row2,col2
              , ijexep(abs(nexep),2)     &   ! (kept) exceptions: brick#1, brick#2
              , expen(abs(nexep))        &   ! (kept) exception expense value
              , stat=nAllocErr )
      if (nAllocErr .ne. 0) then
        write(*,*) 'Out of memory.  Try reducing the model size.'
        stop 
      endif
     
    ! If # exceptions is positive, then use a scheme where the regularization 
    !   row & column numbers are used to identify the blocks.
    ! An expen(i) value of zero means no penalty for those two blocks to be 
    !   different.
      if (nexep .gt. 0) then
        do i = 1,nexep
          read (iof,*, end=198, err=199) (locexp(i,j),j=1,4),expen(i)
        end do
    ! If exception # is negative, use actual brick numbers.
      else if (nexep .lt. 0) then
        do i = 1,-nexep
          read (iof,*, end=198, err=199) (ijexep(i,j),j=1,2),expen(i)
        end do
      end if
      
    ! we are done with the model file
      close (iof)

!cccccccccccccccccccc mesh file reading cccccccccccccccccccccccccccccccc
!
      if (mshtyp(1:4) .ne. 'PW2D') call filerr( &
     &  ' Mesh type not supported', iof)
      open (iof, file=mshfil, status='old', iostat=lerr)
      if (lerr .ne. 0) then
        write(*,*) ' Error opening mesh file'
        stop 
      end if

    ! read in mesh sizes, etc...
      read (iof,'(a)', end=198, err=199) string
      read (iof,*, end=188, err=189) lidx,nodey,nodez,nrfix,lfre,nexec
    ! errors:
      if (lidx .ne. 0) call filerr(' Mesh file not OCCAM file', iof)
      if (nrfix .ge. cnMaxFixedRes) call filerr('Too many fixed resistivities',iof)
      
      ! Allocate things related to node values
      Allocate( iconb(nodeY)                &
              , yy(nodeY), dly(nodeY)       &
              , zz(nodeZ + cnAirLayers)     &
              , dlz(nodeZ + cnAirLayers)    &
              , npt(nodeZ + cnAirLayers,4,nodeY)    &
              , apt(nodeZ,4,nodeY)          &   ! local var doesn't need air layers
              , stat=nAllocErr )
      if (nAllocErr .ne. 0) then
        write(*,*) 'Out of memory.  Try reducing the model size.'
        stop
      endif
      
    ! input values for fixed resistivities, if any
      if (nrfix .ge. 1) then 
        read (iof,*, end=188, err=189) (y(i), i = 2,nrfix+1)
      end if
    ! input frequencies, if any (these are ignored)
      if (lfre .ge. 1) then
        read (iof,*, end=188, err=189) (dummy, i=1,lfre)
      end if
    ! input the element widths and heights
      read (iof,*, end=188, err=189) (dly(i), i = 1,nodey-1)
      read (iof,*, end=188, err=189) (dlz(i), i = 1,nodez-1)
    ! input site information, if any (these are ignored)      
      read (iof,*, end=188, err=189) nsites, (dummy, i = 1,nsites)
    ! input the resistivity codes
      do i = 1, nodez-1
        do l = 1,4
          read (iof,'(1000a1)',end=188,err=189) (apt(i,l,j), j=1,nodey-1)
        enddo
      enddo
      
! We are done with the mesh file
      close (iof)

!ccccccccccccccccccc statics file reading cccccccccccccccccccccccccccccc
!
      if (ssfile(1:4) .eq. 'none' .or. ssfile(1:4) .eq. 'NONE') then
        iscon = 0
        nsb = 0
      else
        open (iof, file=ssfile, status='old', iostat=lerr)
        if (lerr .ne. 0) then
            write(*,*) ' Error opening static shift file'
            stop 
        end if
! Read the stuff in it.
        read (iof,'(18x,a)', end=208, err=209) string
        if (string(1:17) .ne. 'OCCAM2MTSHIFT_1.0') call filerr( &
     &    ' Static shift file not supported', iof)
        read (iof,'(18x,a)', end=208, err=209) string
        read (iof,'(18x,a)', end=208, err=209) string
        read (iof,'(18x,a)', end=208, err=209) string
        iscon = idcode(string, iof)
        read (iof,'(18x,a)', end=208, err=209) string
        sserr = rdcode(string, iof)
        read (iof,'(18x,a)', end=208, err=209) string
        nsb = idcode(string, iof)
        read (iof,'(a)', end=208, err=209) string
        do i = 1,nsb
          read (iof,*,end=208,err=209) isx1(i),isx3(i),shift(i),isfx(i)
        enddo
        close (iof)
      end if
      
!ccccccccccccccccccc prejudice file reading cccccccccccccccccccccccccccc
!
! no prejudice: zero out the model pejudice and associated weights
      premod = 1.0
      prewt = 0
      if (prfile(1:4) .eq. 'none' .or. prfile(1:4) .eq. 'NONE') then
        nprej = 0
        
      else
! read the prejudice file
        open (iof, file=prfile, status='old', iostat=lerr)
        
        if (lerr .ne. 0) then
            write(*,*) ' Error opening model prejudice file'
            stop 
        else
           write(*,*) 'Reading in prejudice file: ',prfile
        end if
! Read the stuff in it.
        read (iof,'(18x,a)', end=218, err=219) string
!KWK 01-24-2006:
! old style prejudice file, needs full prejudice model
        if  (string(1:16) .eq. 'OCCAM2MTPREJ_1.0') then
           read (iof,'(18x,a)', end=218, err=219) string
           read (iof,'(18x,a)', end=218, err=219) string
           nprej = idcode(string, iof)
           if (nprej .gt. nParams) then
             call filerr(' Too much prejudice', iof)
           endif
           read (iof,*,end=218,err=219) (premod(i), i=1,nprej)
           read (iof,*,end=218,err=219) (prewt(i), i=1,nprej)
           close (iof)
        elseif  (string(1:16) .eq. 'OCCAM2MTPREJ_2.0')  then
!KWK new version just lists parameter number, premod and weight 
! for any parameters to be prejudiced
           write(*,*) 'Prejudice file format: OCCAM2MTPREJ_2.0'
           read (iof,'(18x,a)', end=218, err=219) string
           nprej = idcode(string, iof)
           if (nprej .gt. nParams) then
             call filerr(' Too much prejudice', iof)
           endif
! read in prejudiced parameters line by line   
!KWK: still need to clean up the error messages
           write(*,*) 'Prejudiced parameters: ', nprej
           do i = 1,nprej
              read (iof,*,end=218,err=219) j,temp1,temp2
              if (j.gt.nParams) then
                 write(*,*) 'prejudice for parameter:??',j
                 call filerr(' Too much prejudice', iof)
              endif
              premod(j) = temp1
              prewt(j)  = temp2
!              write(*,*) '#,premod,prewt: ',j,premod(j),prewt(j)
           enddo
           close (iof)
        else
           call filerr(' Model prejudice file not supported', iof)
        endif
      
      end if
      

!
!ccccccccccccccc set up mesh and penalty matrix info ccccccccccccccccccc
!


! npt() runs from 0 (0) to 35 (Z), '?' is flagged as -99, to be overwritten 
! from the parameter array later.
      do j = 1,nodey-1
        do l = 1,4
          do i = 1,nodez-1
            npt(i,l,j) = ichar(apt(i,l,j)) - 48
            if (npt(i,l,j) .eq. 15) npt(i,l,j) = -99  ! sends ? to 99
            if (npt(i,l,j) .eq. 42) npt(i,l,j) = -1   !sends Z to seawater
            if (npt(i,l,j) .gt.  9) npt(i,l,j) = npt(i,l,j) - 7 !sends A-Z to 10-35
          enddo
        enddo
      enddo
!
! Set up the free parameter part of npt() and also laprm().
! The following block of code finds the outline of every parameter brick,
! and if any part of it is free, adds it to npt() and laprm().
! This block of code also computes the widths and thickness of the bricks,
! which will be needed later during the computation of the penalty matrix.
!
    ! initialize the number of parameters and the top edge of bricks
    nprm = 0
    inz1 = 1
    ! for every layer in the reg. mesh
    do i = 1,nlay
        ! find the bottom edge of bricks and initialize the right edge of first brick
        inz2 = inz1 + irz(i) - 1
        iny1 = 1
        ! for every brick in that layer 
        do j = 1,nrcol(i)
            ! set irindx() assuming that this brick is fixed (will be tested later)
            irindx(i,j) = 0
            ! find the right edge of brick
            iny2 = iny1 + abs(iry(i,j)) - 1
            ! test to see if the brick is free 
            iffree = .false.
            do ii = inz1,inz2
              do jj = iny1, iny2
                do l = 1,4
                  if (apt(ii,l,jj) .eq. '?') iffree = .true.
                enddo
              enddo
            enddo
            if (iffree) then
                ! it is a free param: incr count and set a new element of laprm
                nprm = nprm + 1
                laprm(nprm) = nprm + nrfix 
                
                ! set irindx() to point to the new parameter
                irindx(i,j) = nprm
                
                ! find the thickness of the brick
                thick(nprm) = 0.0
                do ii = inz1, inz2
                    thick(nprm) = thick(nprm) + dlz(ii)
                enddo
                
                ! find the width of the brick
                width(nprm) = 0.0
                do jj = iny1, iny2
                    width(nprm) = width(nprm) + dly(jj)
                enddo
                
                ! record the FE nodes for the four corners of the block (top, left, bottom, right)
                icorn(nprm,1) = inz1
                icorn(nprm,2) = iny1
                icorn(nprm,3) = inz2+1
                icorn(nprm,4) = iny2+1
                
                ! over the entire brick, including sub-triangles,
                ! test to see if triangular element is free and set npt() if it is
                do ii = inz1,inz2
                  do jj = iny1, iny2
                    do l = 1,4
                      if (apt(ii,l,jj) .eq. '?') then
                        npt(ii,l,jj) = nprm + nrfix 
                      end if
                    enddo
                  enddo
                enddo
            end if
            ! increment the left edge of brick
            iny1 = iny2 + 1
        enddo
        ! increment the top edge of bricks
        inz1 = inz2 + 1
    enddo
!
! copy nprm into nbrick to give routines the number of parameters
! which are actually resistivity bricks (c.f. static shifts)
      nbrick = nprm
!
! Now, the bricks on the edges of the mesh and the bottom of the mesh will have
! widths and thicknesses that are too large to use in regularization penalty.  Set
! them equal to their inside neighbours.  The indexing in this code is the same as 
! above.
      do i = 1,nlay
        do j = 1, nrcol(i)
          ii = irindx(i,j)
          if (ii .gt. 0) then
            ! DGM Oct 2006 - if the row is one block wide, don't do the width adjustments!
            if (nrcol(i) > 1) then
                if (j .eq. 1)        width(ii) = width(ii+1)
                if (j .eq. nrcol(i)) width(ii) = width(ii-1)
            endif
            
            ! DGM Oct 2006 - the line below assumed the last two rows had the same
            !   number of columns ... which was probably BOGUS.  
            !   So now I enforce the assumption.
            if (i .eq. nlay) then
                if(nrcol(i) == nrcol(i-1)) thick(ii) = thick(ii-nrcol(i))
            endif
          end if
        enddo
      enddo
!
! search the penalty exceptions list to find the block numbers
!
      if (nexep .gt. 0) then
!        print*, 'Exceptions between blocks'
        do i = 1,nexep
          ijexep(i,1) = irindx(locexp(i,1),locexp(i,2))
          ijexep(i,2) = irindx(locexp(i,3),locexp(i,4))
          if ((ijexep(i,1).le.0) .and. (ijexep(i,2).le.0)) call filerr( &
     &  ' Penalty requested for two fixed blocks; use prejudice file',0)
!          print*, ijexep(i,1), ijexep(i,2)
        end do
      else
        nexep = -nexep
      end if
      
      write(*,*) '# Exceptions: ', nexep

! search the penalty exceptions list to find the block numbers
!      do 97 i = 1,nexep
!        ijexep(i,1) = irindx(locexp(i,1),locexp(i,2)) 
!        ijexep(i,2) = irindx(locexp(i,3),locexp(i,4)) 
!        if ((ijexep(i,1).le.0) .and. (ijexep(i,2).le.0)) call filerr( 
!     * ' Penalty requested for two fixed blocks; use prejudice file',0)
!97    continue
!c 
!ccccccccccccccccccc set up receiver locations ccccccccccccccccccccccccc
!
! first find the position of the right edge of the left terminating layer: 
      dummy = 0.0
      do i = 1, abs(iry(1,1))
         dummy = dummy + dly(i)
      enddo
! now find cumulative positions of all nodes (borrowing an array 
! used for similar thing in PW2DI), zero out iconb() 
      yy(1) = boffs - dummy
      iconb(1) = 0
      do i = 1,nodey-1
        iconb(i+1) = 0
        yy(i+1) = yy(i) + dly(i)
      enddo
! find the node locations for surface conductivity boundaries
!kk      l = 1
!kk      do 103 i = 1, nrcol(1)
!kk        l = l + iry(1,i)
!kk        iconb(l) = 1
!kk103      continue
! now find the nodes closest to the station locations, avoiding conductivity
! boundaries
      do 110 i = 1, nrc
        do 105 j = 2,nodey
          if (yy(j) .ge. sitloc(i)) then
            if (yy(j)-sitloc(i) .le. sitloc(i)-yy(j-1)) then
              nrl(i) = j
              if (iconb(j) .eq. 1) nrl(i) = j-1
            else
              nrl(i) = j-1
              if (iconb(j-1) .eq. 1) nrl(i) = j
            end if
            goto 110
          end if
105     continue
110   continue
      if (idebug .ge. 1) then
        write(*,*) ' '
        write(*,*) 'station position    node number     actual position'
        do i = 1, nrc
            write(*,'(4X,f10.0,10X,i5,11X,f10.0)') sitloc(i), nrl(i), yy(nrl(i))
        enddo
        write(*,*) ' '
      end if
! now check to see if two stations occupy the same location
!kk who cares!
!kk      do 115 i = 1, nrc
!kk       do 114 j = 1, nrc
!kk          if ((nrl(i) .eq. nrl(j)) .and. (i.ne.j)) then
!kk            call filerr(' Two stations on same node', 0)
!kk         end if
!kk114        continue
!kk115      continue

!ccccccccccccccccccc set up static shift arrays cccccccccccccccccccccccc
!
! Examine the statics blocks, allocate new parameters and set up the
! data constriants
      nfss = 0
      tesum = 0.0
      tmsum = 0.0
      if (nsb .eq. 0) iscon = 0
      do 140 i = 1,nsb
        if (isfx(i) .ne. 0) then
! this static shift is a free parameter: add one more to the list
          nfss = nfss + 1
          if (nbrick + nfss .gt. nParams) call filerr( &
     &     ' No room for free statics in parameter array',0)
          isfx(i) = nbrick + nfss
! add a term to the constraint sums
          if (isx3(i) .eq. 1) then
            tesum = tesum + shift(i)
          else if (isx3(i) .eq. 5) then
            tmsum = tmsum + shift(i)
          else 
            call filerr(' Can only apply statics to resistivities',0)
          end if
        end if
140   continue
!
! Set up the statics constraint
      if (iscon.eq.0) then
! no constraint: do nothing
        continue
      else if (iscon .eq. 1) then
! constraint is to sum both TE and TM shifts together to whatever was the
! sum in the statics file
        nd = ndata + 1
        call Realloc_dp_d_sd
        dp(ndata+1,3) = -1
        d(ndata+1) = tesum + tmsum
        sd(ndata+1) = sserr
      else if (iscon .eq. 2) then
! constraint is to sum both TE and TM shifts separately to whatever was the
! sum in the statics file
        nd = ndata + 2
        call Realloc_dp_d_sd
        dp(ndata+1,3) = -2
        dp(ndata+2,3) = -3
        d(ndata+1) = tesum
        d(ndata+2) = tmsum
        sd(ndata+1) = sserr
        sd(ndata+2) = sserr
      else
        call filerr ( &
     &   ' Statics constraints type greater than 2 not supported', 0)
      end if
!
      !nParams = nbrick + nfss
      if (nParams .ne. nbrick + nfss) then
        if (iscon == 0) then
            call filerr( 'Model file error: # params <> # free bricks', 0 )
        else
            call filerr( 'Model or Static Shift error: # params <> # free bricks + # static shift params', 0 )
        endif
      endif


!cccccccccccccccccccccc clean up and leave cccccccccccccccccccccccccccc
!
    ! Allocate things related to the number of data (and potentially
    !   extra elements for statics).
    ! Use nd instead of ndata because it will reflect the extra elements.
    allocate( dm(nd)               &
            , dwk1(nd)             &
            , dwk3(nd)             &
            , dwk4(nd)             &
            , stat=nAllocErr )
    if (nAllocErr .ne. 0) then
        write(*,*) 'Out of memory.  Try reducing the model size'
        stop
    endif

! We should have read the iteration file already and set nParams.  Test for
! consistency with nbrick and nfss
    if (nParams .ne. nbrick+nfss) &
     &  call filerr (' Model file and nParams do not match', 0)
     
     
    if (abs(idebug) .ge. 1) then
        write(*,*) nodey, ' x ', nodez, ' FE mesh read'
        write(*,*) nParams, ' model parameters defined'
    end if

    return

188      call filerr (' Mesh file ended prematurely', iof)
189      call filerr (' Error reading mesh file', iof)
198      call filerr (' Model file ended prematurely', iof)
199      call filerr (' Error reading model file', iof)
208      call filerr (' Statics file ended prematurely', iof)
209      call filerr (' Error reading statics file', iof)
218      call filerr (' Prejudice file ended prematurely', iof)
219      call filerr (' Error reading prejudice file', iof)
!
      end
      
!-----------------------------------------------------------------------
subroutine Realloc_dp_d_sd
!-----------------------------------------------------------------------
! OCCAM MT2D 3.0 Package
! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
! Subroutine Revision 3.0, November 2006
! When statics params are included in the inversion, they 
! *may* add one or two new rows to 'nd' so that nd = ndata + 1 or 2.
! Various data arrays have to be reallocated to accommodate this.
    use Occam_mod
    use Fwd_mod
    implicit none

    ! Local vars to hold existing values during the realloc.
    real, dimension(:), allocatable     :: d_old, sd_old
    real, dimension(:,:), allocatable   :: dp_old
    integer :: nDim, nAllocErr
    
    nDim = size(d,1)
    
    ! Allocate buffers to hold the old & copy the values
    allocate( d_old(nDim), sd_old(nDim), dp(nDim,cnDataParams) &
            , stat=nAllocErr )
    if (nAllocErr .ne. 0) then
        write(*,*) 'Out of memory.  Try reducing the model size.'
        stop 
    endif
    
    d_old  = d
    sd_old = sd
    dp_old = dp
    
    ! Reallocate the targets & copy the data back over
    deallocate( d, sd, dp )
    allocate( d(nd), sd(nd), dp(nd,cnDataParams) &
            , stat=nAllocErr )
    if (nAllocErr .ne. 0) then
        write(*,*) 'Out of memory.  Try reducing the model size.'
        stop 
    endif
    
    d(:nDim)    = d_old
    sd(:nDim)   = sd_old
    dp(:nDim,:) = dp_old
    
    ! release our local buffers
    deallocate( d_old, sd_old, dp_old )
    
    return
end subroutine Realloc_dp_d_sd


!-----------------------------------------------------------------------
subroutine FwdCalc( bDoPartials, np, nd, pm, dp, dm, d, wj )
!-----------------------------------------------------------------------
! OCCAM MT2D 3.0 Package
! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
! Subroutine Revision 3.0, November 2006
! DGM Oct 2006: Routine combined from old formod & fordiv routines which
!   were identical except for addition of 'wj' handling (jacobian) in
!   fordiv.  Parameter 'wj' is now optional: reqd only if (bDoPartials).
! Notes from previous function headers:
    ! revision to allow for impedance data - cdh -Feb 23, 1994
    ! Impedances are broken up into real and imaginary
    !   parts so that all data can be treated as real
    !   allows for multiplication by decomposition matrix
    ! KWK Feb 23,2006 input d(nd) to compare for  360 degree phase wraps

    ! includes:
    Use Occam_mod, only: idebug
    USE Fwd_mod
    
    IMPLICIT NONE

    ! arguments:
    logical, intent(in) :: bDoPartials
    integer             :: np, nd       ! # parameters, # data
    real, intent(in)    :: pm(np)       ! input model parameters
    real, intent(in)    :: dp(nd,4)     ! 4 cols: see defn in module
    real, intent(out)   :: dm(nd)       ! output data: apparent resistivity, phase, etc....
    real, intent(in)    :: d(nd)        ! original input data
    real, optional      :: wj(nd,np)    ! output jacobians
    
    ! local variables
    real :: tesum, tmsum, p0, pm360, pp360, e10
    !   tesum = sum of logs of TE free static shifts 
    !   tmsum = sum of logs of TM free static shifts 
    integer :: i, j, nres, itype, isite, ifreq
    !   i,j = loop counters
    !   nres = number of resistivities in y() array
    !   itype, isite, ifreq = integer versions of dp() entries
    real :: nRMin, nRMax
    !   nRMin/ax - used for tracking input to PW2DI.  If the range is too large,
    !       the forward model will choke big time.

    !------------------------------------------------
    

    ! total number of resistivities = fixed resistivities + free res'ties 
    ! (but not air)
    nres  = nrfix + nbrick
    tesum = 0.0
    tmsum = 0.0
    e10   = alog(10.)   ! 1/log10(e)

    ! Make sure that the model params are converted into resistivities.
    ! Input resistivities should be placed into array 'y'
    y(1) = 1/3.25   ! init the "water" entry.  Fixed resistivities set in inputd()
    forall( i=1:nbrick ) y(nrfix+1+i) = 10.**pm(i)
    
    ! copy across any static shifts in the free parameters
    if (nfss .gt. 0) then
        do j = 1,nsb
          if (isfx(j) .ne. 0) shift(j) = pm(isfx(j))
        enddo
    end if
    
    ! As an info step, look at the range of the min & max resistivities that
    !   are being passed to the forward model.  If this range is too large
    !   the forward model will fail or produce unstable results.
    if (abs(idebug) >= 1) then
        nRMax = maxval( y(nrfix+2:nrfix+1+nbrick) )
        nRMin = minval( y(nrfix+2:nrfix+1+nbrick) )
        write(*,*) ' '
        write(*,*) '  Fwd model resistivity range:', nRMin, ' ', nRMax
        if (nRMax > 10e7) then
            write(*,*) '!!WARNING!! RESISTIVITIES TOO HIGH.  FORWARD MODEL MIGHT FAIL!'
            write(*,*) 'Moving resistivity back down to finite range!'
            where (y>10e7) y = 10e7
            nRMax = maxval( y(nrfix+2:nrfix+1+nbrick) )
        	nRMin = minval( y(nrfix+2:nrfix+1+nbrick) )
			write(*,*) '  Fwd model resistivity range:', nRMin, ' ', nRMax            
        endif
        if (nRMin < 10e-7) then
            write(*,*) '!!WARNING!! RESISTIVITIES TOO LOW.  FORWARD MODEL MIGHT FAIL!'
            write(*,*) 'Moving resistivity back up to finite range!'
            where (y<10e-7) y = 10e-7
            nRMax = maxval( y(nrfix+2:nrfix+1+nbrick) )
        	nRMin = minval( y(nrfix+2:nrfix+1+nbrick) )
			write(*,*) '  Fwd model resistivity range:', nRMin, ' ', nRMax            
        endif
        
    endif
    
    ! compute the MT response
    if (bDoPartials) then
        call PW2DI(nres, 2) ! 2 = Compute resp + compute jacobian
    else
        call PW2DI(nres, 1) ! 1 = Compute response
    endif
    
    ! For every datum:
    ! Note: As of Oct 2006, derivatives from PW2DI are NOT wrt 
    !   log10(resistivity).  They need to be converted to log space.
    do i = 1, ndata

        ! strip off data from dp(:,4)
        !  dp(:,1) = site number
        !  dp(:,2) = frequency number
        !  dp(:,3) = data type (see defn of dp in module)
        !  dp(:,4) = input strike angle
        isite = int(dp(i,1)+0.1)
        ifreq = int(dp(i,2)+0.1)
        itype = int(dp(i,3)+0.1)

        select case (iType)
        case (1)                                ! TE Resistivity (log10)
            dm(i) = alog10(rhte(isite,ifreq))
            if (bDoPartials) then
                ! from drho_a/dsigma_m to dlog10rho_a/dlog10rho_m
                forall( j=1:nbrick ) wj(i,j) = -rted(j,iSite,iFreq) / rhte(iSite,iFreq) / y(j+1)
            endif
            
        case (9)                               ! TE Resistivity (linear)
            dm(i) = rhte(isite,ifreq)
            if (bDoPartials) then
                wj(i,1:nbrick) = rted(1:nbrick,iSite,iFreq)
            endif
            
        case (2)                                ! TE Phase
            dm(i) = phte(isite,ifreq)
            
            ! KWK Feb 23, 2006, phase shifting:
            p0    = abs(dm(i) - d(i))
            pm360 = abs(dm(i) - 360. -d(i))
            pp360 = abs(dm(i) + 360. -d(i))
            
            if (pm360.lt.p0) then  ! +- 360 degrees if misfit is better
                dm(i) = dm(i) - 360.
            elseif (pp360.lt.p0) then
                dm(i) = dm(i) + 360.
            endif


            if (bDoPartials) then
                wj(i,1:nbrick) = pted(1:nbrick,iSite,iFreq)
            endif

        case (3)                                ! real(tipper)
            dm(i) = real(kzte(isite,ifreq))
            if (bDoPartials) then
                wj(i,1:nbrick) = real(kted(1:nbrick,iSite,iFreq))
            endif
            
        case (4)                                ! imag(tipper)
            dm(i) = aimag(kzte(isite,ifreq))
            if (bDoPartials) then
                wj(i,1:nbrick) = aimag(kted(1:nbrick,iSite,iFreq))
            endif
            
        case (5)                                ! TM resistivity (log10)
            dm(i) = alog10(rhtm(isite,ifreq))
            if (bDoPartials) then
                ! from drho_a/dsigma_m to dlog10rho_a/dlog10rho_m
                forall( j=1:nbrick ) wj(i,j) = -rtmd(j,iSite,iFreq) / rhtm(iSite,iFreq) / y(j+1)
            endif
            
        case (10)                               ! TM resistivity (linear)
            dm(i) = rhtm(isite,ifreq)
            if (bDoPartials) then
                wj(i,1:nbrick) = rtmd(1:nbrick,iSite,iFreq)
            endif
            
        case (6)                                ! TM phase
            dm(i) = phtm(isite,ifreq)
            
            ! KWK Feb 23, 2006, phase shifting:                    
            p0    = abs(dm(i) - d(i))
            pm360 = abs(dm(i) - 360. -d(i))
            pp360 = abs(dm(i) + 360. -d(i))
            if (pm360.lt.p0) then  ! +- 360 degrees if misfit is better
                dm(i) = dm(i) - 360.
            elseif (pp360.lt.p0) then
                dm(i) = dm(i) + 360.
            endif
            
            if (bDoPartials) then
                wj(i,1:nbrick) = ptmd(1:nbrick,iSite,iFreq)
            endif

            
        case (11, 12, 17, 18)                     ! real(Zxx), imag(Zxx), ... Zyy
            dm(i) = 0.
            if (bDoPartials) &
                wj(i,1:nbrick) = 0

        case (13)                                ! real(Zxy)
            dm(i) = real(zimpe(isite,ifreq))
            if (bDoPartials) then
                wj(i,1:nbrick) = real(zimped(1:nbrick,iSite,iFreq))
            endif
            
        case (14)                               ! imag(Zxy)
            dm(i) = aimag(zimpe(isite,ifreq))
            if (bDoPartials) then
                wj(i,1:nbrick) = aimag(zimped(1:nbrick,iSite,iFreq))
            endif
            
        case (15)                               ! real(Zyx)
            dm(i) = real(zimpm(isite,ifreq))
            if (bDoPartials) then
                wj(i,1:nbrick) = real(zimpmd(1:nbrick,iSite,iFreq))
            endif
            
        case (16)                               ! imag(Zyx)
            dm(i) = aimag(zimpm(isite,ifreq))
            if (bDoPartials) then
                wj(i,1:nbrick) = aimag(zimpmd(1:nbrick,iSite,iFreq))
            endif
            
        case default
            write(*,*) 'Unknown Data type:', iType
            write(*,*) 'Need to implement this type in FwdCalc().  Stopping.'
            stop
        end select
        
        ! Convert from d.../dsigma_m to d.../dlog10rho_m.  The conversion is the same
        !   for almost all derivatives because only the denominator of the derivative
        !   (which doesn't depend on the data type) is being converted:
        !   dG/dlog10rho_m = dG/dsigma_m * dsigma_m/drho_m * drho_m/dlog10rho_m
        !       <want>         <have>            <conversion factors>
        !   Similar with bounded constraints where denom is not dlog10rho_m but some
        !   other function of rho_m.
        if (bDoPartials) then
            ! TE & TM log10 resistivity are done in the case above because 
            !   their transformations are slightly different.
            if (iType .ne. 1 .and. iType .ne. 5) then
                forall( j=1:nbrick ) 
                	wj(i,j) = -1*wj(i,j)*e10 / y(j+1)
                end forall
            endif
        endif

        ! Handle the static shifts
        if (bDoPartials) then
            ! Zero out any additional columns of the Jacobian matrix associated
            ! with free static shifts (there might not be a free shift for this
            ! datum)
            do j = 1, nfss
                wj(i,nbrick+j) = 0.0
            enddo
        endif
        
        ! Search the static shifts to see if one applies to this datum:
        do j = 1, nsb
            if ((isx1(j).eq.isite) .and. (isx3(j).eq.itype)) then
                ! apply the static shift
                dm(i) = dm(i) + shift(j)
                
                ! If this shift is a free parameter set the augmentation of WJ to 1.0
                if (bDoPartials .and. isfx(j) .ne. 0) wj(i,isfx(j)) = 1.0
            end if
        enddo

    enddo   ! all data
    
    ! Now compute the sums used in the static shift constraint
    do j = 1, nsb
        if (isfx(j) .ne. 0) then
            ! compute the sum of shifts for the shift constrains
            if (isx3(j) .eq. 1) then
                tesum = tesum + shift(j)
            else if (isx3(j) .eq. 5) then
                tmsum = tmsum + shift(j)
            end if
        end if
    enddo

    ! Now add rows to the Jacobian matrix for static shift 
    ! constraints and create the synthetic constraint data
    select case (iscon)
    case (0)        ! no constraint: do nothing
        continue
    
    case (1)        ! constraint is to sum both TE and TM shifts together
        dm(ndata+1) = tesum + tmsum
        
        if (bDoPartials) then
            ! zero out those cols of wj associated with model resistivities
            wj(nData+1,1:nBrick) = 0.0
            
            ! set those cols of wj associated with free static shifts 
            ! to all ones for TE + TM summing
            do j = 1, nfss
                wj(ndata+1,nbrick+j) = 1.0
            enddo
        endif
        
    case (2)        ! constraint is to sum both TE and TM shifts separately 
        dm(ndata+1) = tesum
        dm(ndata+2) = tmsum
        
        if (bDoPartials) then
            ! zero out those cols of wj associated with model resistivities
            wj(ndata+1:ndata+2,1:nBrick) = 0.0
            
            ! set those cols of wj associated with free static shifts 
            ! (need to test for either TE or TM sum)
            do j = 1, nsb
                if (isfx(j) .ne. 0) then
                    if (isx3(j) .eq. 1) then        ! TE constraint
                        wj(ndata+1,isfx(j)) = 1.0
                        wj(ndata+2,isfx(j)) = 0.0
                    else if (isx3(j) .eq. 5) then   ! TM constraint
                        wj(ndata+1,isfx(j)) = 0.0
                        wj(ndata+2,isfx(j)) = 1.0
                    end if
                end if
            enddo
        endif
    end select
    
    return
end subroutine FwdCalc


!-----------------------------------------------------------------------
subroutine OutFwd( cFile )
!-----------------------------------------------------------------------
! OCCAM MT2D 3.0 Package
! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
! Subroutine Revision 3.0, November 2006
    use Occam_Mod
    use Fwd_mod
    implicit none
    
    ! Arguments
    character(*), intent(in)   :: cFile
    
    ! Local vars
    integer :: i, j
    
    open( unit=42, file=cFile )
    write(42,'(A)')     'FORMAT:           OCCAM2MTDATA_1.0 '
    write(42,'(A)')     'TITLE:            Forward model response from Occam2D'
    write(42,'(A,I4)')  'SITES:            ', nRc
    write(42,'(A)') (trim(cSiteName(i)), i=1,nRc)
    write(42,'(A)')     'OFFSETS (M):'
    write(42,'(G14.6)') (SitLoc(i), i=1, nRc)
    write(42,'(A,I4)')  'FREQUENCIES:      ', nFre
    write(42,'(G14.6)') (frek(i), i=1, nFre)
    write(42,'(A,I4)')  'DATA BLOCKS:      ', nd
    write(42,'(A)')     'SITE   FREQ   TYPE          DATUM          ERROR'
    do i=1,nd
        write(42,'(3I6,2(2X,G14.6))') (INT(DP(i,j)), j=1,3),DM(i),SD(i)
    enddo
    
    close(42)
    
    return
end subroutine OutFwd


! **********************************************************************
      SUBROUTINE PW2DI(NRES,NEXEC)
! **********************************************************************
! OCCAM MT2D 3.0 Package
! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
! Subroutine Revision 3.0, November 2006
! Original code written by Phil Wannamaker.
! July 25, 2001 - alteration to FwdCalc, PW's new adjoint method
!       used, which returns derviatives with respect to log(resistivity)
! Additions & edits by Catherine DeGroot-Hedlin, Kerry Key, & David Myer.
!
!========================== INCLUDE FILES ==============================
!
      use Occam_Mod
      USE Fwd_mod


      IMPLICIT NONE

!
!============================ INTERFACE ================================
!
! WRITTEN BY:  P. WANNAMAKER, UURI.
!
! DESCRIPTION: 
!
!     FOR THE FORWARD PROBLEM AND JACOBIANS,
!     FINITE ELEMENTS ARE USED WITH LINEAR BASIS FUNCTIONS.  TOPO-
!     GRAPHIC OR BATHYMETRIC VARIATIONS ACROSS THE PROFILE CAN BE
!     INCLUDED AS WELL (WANNAMAKER ET AL, 1986, GEOPHYSICS).  MKS
!     UNITS AND EXP(+JWT) TIME DEPENDENCE ARE ASSUMED.  THE + COORD-
!     INATE DIRECTIONS FOR THE PROBLEM ARE X- NORTH, Y-EAST, Z-DOWN.
!     STRIKE DIRECTION IS ASSUMED N-S.  ORIGIN IS AT A USER-SPECIFIED
!     NODAL COLUMN AT THE HIGHEST TOPOGRAPHIC POINT.  A SECONDARY
!     FIELD FORMULATION OF THE FINITE ELEMENT EQUATIONS HAS BEEN
!     UTILIZED HERE TO IMPROVE ACCURACY (WANNAMAKER ET AL, 1987, GEOP.
!     J. ROY. ASTR. SOC.).  THE LAYERED HOST FOR THE INHOMOGENEITY BY
!     CONVENTION IS TAKEN TO BE THAT WITH THE HIGHER RESISTIVITY IN
!     THE UPPERMOST LAYER.  IN THIS VERSION, THE PARAMETER JACOBIANS
!     ARE COMPUTED BY TAKING ADVANTAGE OF RECIPROCITY.  THERE IS A
!     REPORT DOCUMENTATION BY WANNAMAKER (1990) DESCRIBING THE METHOD.
!
!     ROUTINES NEEDED:  PW2DI, PARFLD, COAMP, ZLFLD, GAUSS, AUXFLD,
!                       BNORM, DERIV, APPARMT, addair, RECIJAC, PARJAC, 
!                       GAUSSD, AUXADJ, APPRMTJ, DRVADJ
!
! ON INPUTS:
!         NRES - NO. OF UNIQUE RESISTIVITY MEDIA IN THE MESH
!                (INCLUDES ANY SEAWATER BUT NOT THE AIR)
!         NFRE - NO. OF FREQUENCIES TO BE MODELED
!         NEXEC - IF .EQ. 0, JUST PRINT INPUT PARAMETERS AND MESH
!                 TO CHECK CORRECTNESS
!               - IF .EQ. 1, COMPUTE RESPONSE
!               - IF .EQ. 2, INVERT DATA
!         Y - ARRAY OF RESISTIVITIES OF MEDIA MAKING UP THE MESH;
!             LAST VALUE IS THAT OF ANY SEAWATER.
!         FREK - ARRAY OF FREQUENCIES
!         NRL - SITE LOCATIONS BY FE MESH NODE NUMBER
!
!         /BLKI/ IDX - PARAMETER DEFINING MODE OF PROBLEM TO BE SOLVED
!             - IF .EQ. 1, SOLVE BOTH TE AND TM
!             - IF .EQ. 2, SOLVE TE ONLY
!             - IF .EQ. 3, SOLVE TM ONLY
!         /BLKI/ NODEY - NO. OF NODES IN THE Y-DIRECTION OF THE FINITE
!                 ELEMENT MESH
!         /BLKI/ NODEZ - NO. OF NODES IN THE Z-DIRECTION OF THE FINITE
!                 ELEMENT MESH
!         /BLKI/ NPRM - NO. OF POLYGONAL MEDIA OF INVERSION MODEL
!         /BLKI/ LAPRM - INTEGER LABELS OF EACH INVERSION MEDIUM
!         /BLKD/ DLY - ARRAY OF ELEMENT COLUMN WIDTHS
!         /BLKD/ DLZ - ARRAY OF ELEMENT ROW HEIGHTS
!         /BLKE/ NPT - INTEGER LABELS FOR EACH TRIANGULAR ELEMENT OF EACH 
!               RECTANGULAR ELEMENT OF THE MESH.  THESE LABELS ARE 
!               ASSIGNED THE VALUES IN THE RESISTIVITY ARRAY Y.  
!               LABEL 0 IS RESERVED FOR AIR AND LABEL -1 IS RESERVED
!               FOR SEAWATER.  TRIANGULAR ELEMENT LABELING ORDER IS 
!               TOP, LEFT, BOTTOM AND RIGHT.
!
!============== PARAMETER AND LOCAL VARIABLE DECLARATIONS ==============
!
!     ---- PASS PARAMETERS
      INTEGER NRES,NEXEC
!
!     ---- LOCAL VARIABLES
      INTEGER lidx,NODEY1,MAIR,NODEZ1
      INTEGER I,J,K,L,M,III,NYRB
      real conds(0:NRES+1)
!
!
!=======================================================================
      
      ! Make sure the things we need are allocated.  If already allocated
      ! this call will do no harm.  Make sure that there are enough extra
      ! Z nodes to account for the air layers added (them removed) by the
      ! calls to AddAir().
      call Fwd_Allocate( NRC, NRES, nParams, nfre, nlay, nodey, nodez + cnAirLayers)
      
      NODEY1 = NODEY - 1
!
!     ---- DEFINE LOCATION OF Y=0 AT CENTER OF MESH
      YY(1) = 0.0
      DO 20 I = 2,NODEY
         YY(I) = YY(I-1) + DLY(I-1)
20    CONTINUE
!
      DO 30 I = 1,NODEY
        YY(I) = YY(I) - YY(NODEY)/2.
30    CONTINUE
!
! add the air layer if necessary
      lidx = idx
      call addair(lidx, mair)   ! if TE involved, will extra air layers and modify nodez
      
      NODEZ1 = NODEZ - 1
      NTR = MAIR + 1
!
!     ---- CHANGE FROM RESISTIVITIES INPUT TO CONDUCTIVITIES AND
!          ASSOCIATE CONDUCTIVITIES WITH THE MESH INPUT CODE
      conds(1) = 1.0E-18  ! air conductivity
      conds(0) = 3.25     ! seawater conductivity
      DO 90 I = 2,NRES+1
        conds(i) = 1./Y(I)
90    CONTINUE
!
!      ISEA=NRES+1
      ISEA = -1   !Changed by Kerry, Aug 17, 2001
!      conds(ISEA+1) = 1./Y(NRES+1)
      DO 100 J = 1,NODEZ1
        DO 101 L = 1,NODEY1
          DO 102 M = 1,4
            K = NPT(J,M,L)
            SIG(J,M,L) = conds(K+1) ! K=0 is air, K=-1 is seawater
102       CONTINUE
101     CONTINUE
100   CONTINUE
!
!     ---- POSITION THE DATUM (HIGHEST TOPOGRAPHIC POINT)
      L = NTR
      ZZ(L) = 0.
      DO 110 I = NTR,NODEZ1
        ZZ(L+1) = ZZ(L) + DLZ(I)
        L = L + 1
110   CONTINUE
!
      IF(MAIR .NE. 0)THEN
        L = NTR
        DO 130 III = 1,MAIR
          I = MAIR + 1 - III
          ZZ(L-1) = ZZ(L) - DLZ(I)
          L = L - 1
130     CONTINUE
      END IF
!
!     ---- SPECIFY NODAL ROW DEPRESSIONS OF TOPO/BATH SURFACE.
      DO 140 I = 2,NODEY1
        NTAO(I) = 0
        DO 141 J = NTR,NODEZ1
          IF(NPT(J,4,I-1) .EQ. 0  .AND. NPT(J,2,I) .EQ. 0 )GO TO 142
          IF(NPT(J,4,I-1) .EQ. ISEA .AND. NPT(J,2,I) .EQ. ISEA)GO TO 142
          GO TO 140
142       CONTINUE
          NTAO(I) = J - MAIR
141     CONTINUE
140   CONTINUE
      NTAO(1) = NTAO(2)
      NTAO(NODEY) = NTAO(NODEY1)
!      write(*,*),'ntao',(NTAO(i),i=1,nodey)
!
!     ---- DISCERN BOUNDARY LAYERING PARAMETERS FROM INPUT MESH
      NYRB = 1
      DO 150 I = 1,2
        IF(I .EQ. 2)NYRB = NODEY1
        NOLA(I) = 1
        NZZB(I) = NTR
        DO 160 J = NTR,NODEZ1-1
!
!         ---- PATCH TO PREVENT ARRAY OVERRUN
          IF (NOLA(I) .GE. nLay)GOTO 160
          IF(NPT(J,1,NYRB) .EQ. 0)THEN
            NZZB(I) = J + 1
          ELSE
            RLR(NOLA(I),I) = 1./SIG(J,2*I,NYRB)
            IF(NPT(J,2*I,NYRB) .NE. NPT(J+1,2*I,NYRB))THEN
              DLR(NOLA(I),I) = ZZ(J+1) - ZZ(NZZB(I))
              NOLA(I) = NOLA(I) + 1
            END IF
          END IF
160     CONTINUE
        RLR(NOLA(I),I) = 1./SIG(NODEZ1,2*I,NYRB)
150   CONTINUE
! calling new version of PARFLD, etc
!
      NDX = IDX
      NMODE = 1
      IF(NDX .EQ. 1) NMODE = 2
      if (abs(idebug) >= 1) write(*,'("   Freq #: ",$)') ! suppress CR
      DO 180 LFRQ = 1,NFRE
        if (abs(idebug) >= 1) write(*,'(I3,$)') LFRQ       ! suppress CR
!
        DO 200 IMD = 1,NMODE
!
!         ---- FORM AND SOLVE GLOBAL SYSTEM OF EQUATIONS
          CALL PARFLD
!
!         ---- CALCULATE AUXILIARY FIELDS, REAL AND IMAGINARY, 
!              AND APPARENT RESISTIVITY AND PHASE AND TIPPER
!       
          CALL APPARMT
!
!  CALCULATE DATA-PARAMETER JACOBIAN IF INVERSION IS SELECTED.
!
          IF(NEXEC .EQ. 2)THEN
             CALL RECIJAC(NRES)
          END IF
!
200     CONTINUE
180   CONTINUE
      if (abs(idebug) >= 1) write (*,*) ' '   ! force CR at end of line
      IDX=NDX
!
! remove the air layer if necessary
      call addair(0, mair)      ! if TE involved, will remove extra air layers added above

      RETURN
END SUBROUTINE PW2DI

!
!-----------------------------------------------------------------------
      subroutine addair(lidx, mair)
!-----------------------------------------------------------------------
!
! OCCAM MT2D 3.0 Package
! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
! Subroutine Revision 3.0, November 2006
! Subroutine Revision 1.01, 13 Jan 1993
! 
! addair adds or removes an air layer of 10 nodes 
! from the FE mesh (nodez, dlz(), npt() in common blocks)
! if a TE mode calculation is being done.  Much of this code has been
! abstracted from Phil Wannamaker's PW2DI
!
! on input:
!   lidx = flag for direction and action (local version of one in PW2DI): 
!           = 1 or 2 causes addition of air layer
!           = 3 sets mair = 0 and does nothing else
!           = otherwise removes air layer  
!   mair = unchanged from any previous call to addair
!
! on ouput:
!   mair = number of air layers, 0 for TM, 10 for TE
!   /blki/nodez, /blkd/dlz() and /blke/npt() are modified for TE calculation
!
! calls:
!   no other routines
!
!-----------------------------------------------------------------------
! includes:
      USE Fwd_mod

      IMPLICIT NONE

!
! arguments 
      integer lidx, mair
! local variables
      REAL DZF,XA,XB,FXA,FXM,XM
      integer i,j,l,m

!-----------------------------------------------------------------------
      if (lidx .eq. 3) then
        MAIR = 0
        go to 9999
      else if (lidx.eq.1  .or. lidx.eq.2) then
!
        MAIR    = cnAirLayers
        NODEZ   = NODEZ + MAIR
!
! make room in dlz() and npt() for air layers (push things down mair elements)
        do  j = nodez-mair-1, 1, -1
          dlz(j+mair) = dlz(j)
          do  l = 1, nodey-1
            do  m = 1,4
               npt(j+mair,m,l) = npt(j,m,l)
            enddo
          enddo
        enddo
!
! compute the thickness of the elements in the air layer
! TOTAL AIR THICKNESS IS ONE-HALF MESH WIDTH, I.E., 
! YY(NODEY) ABOVE.  ELEMENT ROW JUST ABOVE SURFACE 
! EQUALS THAT JUST BELOW AND ROW THICKNESSES INCREASE
! EXPONENTIALLY UPWARD.  IN ALL, THERE ARE MAIR ELEMENT 
! ROWS OF AIR.
        DZF = YY(NODEY)/DLZ(MAIR + 1)
        XA = 1.
        XB = DZF**0.25
        DO 40 I = 1,50
          XM = (XA + XB)/2.
          IF((XM-XA)/XM .LT. 1.0E-04)GO TO 41
          FXA = 1.
          FXM = 1.
          DO 42 J = 1,MAIR-1
            FXA = 1. + XA*FXA
            FXM = 1. + XM*FXM
 42       CONTINUE
          FXA = FXA - DZF
          FXM = FXM - DZF
!
          IF(SIGN(1.,FXA)*SIGN(1.,FXM) .LE. 0)THEN
            XB = XM
          ELSE
            XA = XM
          END IF
 41       CONTINUE
 40     CONTINUE
!
        DLZ(MAIR) = DLZ(MAIR + 1)
        DO 43 I = 1,MAIR-1
          J = MAIR - I
          DLZ(J) = XM*DLZ(J+1)
43      CONTINUE
!
! fill the air layer lookup table with 0 (air)
        DO  j = 1,MAIR
          DO  l = 1,NODEY-1
            DO  m = 1,4
               NPT(j,m,l) = 0
            enddo
          enddo
        enddo
!
      else
! remove air layers in dlz() and npt() (push things up mair elements)
        if (mair .ne. 0) then
          do  j = 1, nodez-mair-1
            dlz(j) = dlz(j+mair)
            do  l = 1, nodey-1
              do  m = 1,4
                npt(j,m,l) = npt(j+mair,m,l)
              enddo
            enddo
          enddo
          nodez = nodez - mair
        end if
      end if
      
9999  continue
      return
end subroutine addair
      
!-----------------------------------------------------------------------
      SUBROUTINE PARFLD
!-----------------------------------------------------------------------
! OCCAM MT2D 3.0 Package
! Subroutine Revision 3.0, November 2006 (notional only to track future fixes)
! Original code written by Phil Wannamaker.
!
!     THIS ROUTINE COMPUTES THE GLOBAL FINITE ELEMENT SYSTEM TO GIVE
!     THE FIELD PARALLEL TO STRIKE.  TASKS PERFORMED INCLUDE THE
!     PRIMARY FIELDS FOR THE SOURCE VECTOR AND AT EACH NODAL COLUMN,
!     LOADING THE MATRIX ELEMENTS INTO THE GLOBAL SYSTEM, INCORPORA-
!     TING ANY DIRICHLET BOUNDARY CONDITIONS, AND SOLVING THE SYSTEM.
!
      USE Fwd_mod

      IMPLICIT NONE
      
      COMPLEX CP(4),CP12,CP23,CP34,CP14,CP1234,A11,A12,A13,A15,A22, &
     &A24,A25,A33,A34,A35,A44,A45,A55,E11,G11,C11,H11,E22,H22,C22, &
     &E33,G33,E44,A55I,ZLFLD,FSE,FSH 
      COMPLEX CK(4),CK12,CK13,CK14,CK23,CK24,CK34
      COMPLEX FP1(nodeZ),FP2(nodeZ),FP3(nodeZ),FP4(nodeZ),DJ1,DJ2,DJ3,DJ4, &
     &  DSS1,DSS2,DSS3,DSS4,FP5(nodeZ),FP6(nodeZ),FP7(nodeZ),FP8(nodeZ),SM1, &
     &  SM2,SM3,SM4,SM5,RU1,RU2,RU3,RU4,RU5,RC1,RC2,RC3,RC4
! miscellaneous stuff that needs to be defined
      integer i,nodey1,nodez1,nband,nban1,nlyr,nlyr1,lzp,j,jm1,k,l1,l, &
     &   m,l2,l3,l4,nji,iyn,it,ib
      real pi,rad,w,wu,we,ra,zap,rhal,phal,rhar,phar,zt,z1,z2,z3,z5,z6,&
     &   z7,z8,zl,dlzy,dlyz,a,b,c,c4,dy4

      IF(NDX .EQ. 1)IDX = 2
      IF(IMD .EQ. 2)IDX = 3
      IF(IDX .EQ. 2)THEN
       NTQ = 1
       NDZ = NODEZ
       NDZ1 = NDZ - 1
       DO 10 I = 1,NODEY
        NTA(I) = NTAO(I) + NTR + NDZ*(I-1)
!   WRITE(*,*),'NTA(I),TE MODE: ', NTA(I)
10     CONTINUE
      ELSE IF(IDX .EQ. 3)THEN
       NTQ = NTR
       NDZ = NODEZ - NTR + 1
       NDZ1 = NDZ - 1
       DO 15 I = 1,NODEY
        NTA(I) = NTAO(I) + 1 + NDZ*(I-1)
!   WRITE(*,*),'NTA(I),TM MODE: ', NTA(I)
15     CONTINUE
!      write(*,*),'nta: ',(nta(i),i=1,nodey)
      END IF
      NODEY1 = NODEY - 1
      NODEZ1 = NODEZ - 1
      NBAND = NDZ + 2
      NNODE = NODEY*NDZ
      NBAN1 = NBAND - 1
!
!     COMPUTE APPARENT RESISTIVITIES AND IMPEDANCE PHASES FOR LEFT
!     AND RIGHT SIDE LAYER SEQUENCES.
!
      PI = 3.1415927
      RAD = 180./PI
      W = 2.*PI*FREK(LFRQ)
      WU = W*PI*4.E-07
      WE = W*8.85433E-12
      RA = 1.0E+18
      NLYR = NOLA(1)
      NLYR1 = NLYR - 1
      IF (NLYR .EQ. 1) GO TO 20
       DO 30 I = 1,NLYR1
        DDP(I) = DLR(I,1)
        RL(I) = RLR(I,1)
30     CONTINUE
20    CONTINUE
      RL(NLYR) = RLR(NLYR,1)
      CALL COAMP(WU,WE,NLYR)
      IF(IMD .EQ. 1) THEN
       ZAP = 0.
       IF(NPT(NZZB(1),1,1) .EQ. ISEA)ZAP = DDP(1)
       FSE = ZLFLD(WU,ZAP,1,NLYR)
       FSH = ZLFLD(WU,ZAP,2,NLYR)
       RHAL = CABS(FSE/FSH)**2/WU
       PHAL = RAD*ATAN2(AIMAG(FSE/FSH),REAL(FSE/FSH))
      END IF
!
      NLYR = NOLA(2)
      NLYR1 = NLYR - 1
      IF(NLYR .EQ. 1)GO TO 40
       DO 50 I = 1,NLYR1
        DDP(I) = DLR(I,2)
        RL(I) = RLR(I,2)
50     CONTINUE
40    CONTINUE
      RL(NLYR) = RLR(NLYR,2)
      CALL COAMP(WU,WE,NLYR)
      IF (IMD .EQ. 1) THEN
       ZAP = 0.
       IF(NPT(NZZB(2),1,NODEY1) .EQ. ISEA)ZAP = DDP(1)
       FSE = ZLFLD(WU,ZAP,1,NLYR)
       FSH = ZLFLD(WU,ZAP,2,NLYR)
       RHAR = CABS(FSE/FSH)**2/WU
       PHAR = RAD*ATAN2(AIMAG(FSE/FSH),REAL(FSE/FSH))
      END IF
!
!     THIS SECTION CALCULATES THE PRIMARY E- OR H-FIELDS FOR THE TE
!     AND TM PROBLEMS.  FIELDS IN THE LAYERED HOST ARE CALCULATED FOR
!     AN INCIDENT E-FIELD AMPLITUDE OF (1.,0.) AT THE EARTH'S SURFACE
!     (NOT AT THE TOP OF THE MESH) POLARIZED IN THE +X DIRECTION FOR
!     THE TE MODE AND IN THE +Y DIRECTION FOR THE TM MODE.
!     ZLFLD(Z,1) IS FOR E AND ZLFLD(Z,2) IS FOR H.  
!
      LSID = 1
      IF(RLR(1,1) .LT. RLR(1,2))LSID = 2         !WHICH SIDE IS HOST
!
      NLYR = NOLA(LSID)
      NLYR1 = NLYR - 1
      IF(NLYR .EQ. 1)GO TO 60
       DO 70 I = 1,NLYR1
        DDP(I) = DLR(I,LSID)
        RL(I) = RLR(I,LSID)
70     CONTINUE
60    CONTINUE
      RL(NLYR) = RLR(NLYR,LSID)
      CALL COAMP(WU,WE,NLYR)
      DO 80 I = NTQ,NODEZ1
       ZT = ZZ(I) - ZZ(NZZB(LSID))
       Z1 = ZT + DLZ(I)/6.
       Z2 = ZT + DLZ(I)/2.
       Z3 = ZT + DLZ(I)*(5./6.)
       LZP = 0
       IF(Z2 .GE. 0.)LZP = 1
       IF(NLYR .EQ. 1)GO TO 100
        DO 110 J = 2,NLYR
         JM1 = J - 1
         IF(Z2 .GT. DDP(JM1))LZP = J
110     CONTINUE
100    CONTINUE
       IF(LZP .EQ. 0)GO TO 120
        SIGH(I) = CMPLX(1./RL(LZP),WE)
!        write(*,*) 'RL(LZP): ', RL(LZP)
        GO TO 130
120    CONTINUE
        SIGH(I) = CMPLX(1./RA,WE)
!        write(*,*) 'RA: ', RA
130    CONTINUE
       IF(IDX .EQ. 3)GO TO 140
        FP1(I) = ZLFLD(WU,Z1,1,NLYR)
        FP2(I) = ZLFLD(WU,Z2,1,NLYR)
        FP3(I) = ZLFLD(WU,Z3,1,NLYR)
        FP4(I) = FP2(I)
        GO TO 90
140    CONTINUE
        FP1(I) = -ZLFLD(WU,Z1,2,NLYR)
        FP2(I) = -ZLFLD(WU,Z2,2,NLYR)
        FP3(I) = -ZLFLD(WU,Z3,2,NLYR)
        FP4(I) = FP2(I)
        Z5 = ZT
        Z6 = ZT + DLZ(I)/4.
        Z7 = ZT + DLZ(I)*(3./4.)
        Z8 = ZT + DLZ(I)
        FP5(I) = ZLFLD(WU,Z5,1,NLYR)
        FP6(I) = ZLFLD(WU,Z6,1,NLYR)
        FP7(I) = ZLFLD(WU,Z7,1,NLYR)
        FP8(I) = ZLFLD(WU,Z8,1,NLYR)
90     CONTINUE
80    CONTINUE
      DO 150 I = NTQ,NODEZ
       ZL = ZZ(I) - ZZ(NZZB(LSID))
       IF(IDX .EQ. 2)THEN
        FPA(I) = ZLFLD(WU,ZL,1,NLYR)
        FPB(I) = ZLFLD(WU,ZL,2,NLYR)
       ELSE
        FPA(I) = -ZLFLD(WU,ZL,2,NLYR)
        FPB(I) = ZLFLD(WU,ZL,1,NLYR)
       END IF
150   CONTINUE
!
!     THIS SECTION CALCULATES A 4 X 4 ELEMENT MATRIX FOR RECTANGULAR
!     ELEMENT MADE UP OF 4 TRIANGULAR ELEMENTS. THE 5 X 5 MATRIX IS
!     OBTAINED FROM PROPER COMBINATION OF THE FOUR 3 X 3 TRIANGULAR
!     ELEMENT MATRICIES AND REDUCED TO A 4 X 4 MATRIX THROUGH STATIC
!     CONDENSATION. THIS ELIMINATES THE DEGREE OF FREEDOM ASSOCIATED
!     WITH THE INTERNAL NODE OF THE RECTANGLE, THUS REDUCING THE
!     SIZE OF THE GLOBAL SYSTEM MATRIX.
!
      DO 160 I = 1,NNODE
       RM(I) = (0.,0.)
160   CONTINUE
      DO 170 K = 1,NBAND
       DO 180 I = 1,NNODE
        S(K,I) = (0.,0.)
180    CONTINUE
170   CONTINUE
!
      L1 = 0
      DO 190 L = 1,NODEY1
       L1 = L1 + 1
       DO 200 M = NTQ,NODEZ1
        DLZY = DLZ(M)/(2.*DLY(L))
        DLYZ = DLY(L)/(2.*DLZ(M))
        A = -(DLZY + DLYZ)/2.
        B = (DLZY - DLYZ)/2.
        C = DLY(L)*DLZ(M)/48.
        C4 = C*4.
        DY4 = DLY(L)/4.
!
!     CALCULATE THE APPROPRIATE COEFFICIENTS OF THE HELMHOLTZ
!     EQUATION FOR INSERTION INTO THE ELEMENT MATRICIES
!     ACCORDING TO WHETHER THE TE OR TM MODE IS BEING SOLVED.
!
        IF(IDX .EQ. 2)THEN
         DO 210 I = 1,4
          CK(I) = 1./CMPLX(0.,WU)
          CP(I) = -CMPLX(SIG(M,I,L),WE)
210      CONTINUE
        ELSE
         DO 220 I = 1,4
          CK(I) = 1./CMPLX(SIG(M,I,L),WE)
          CP(I) = -CMPLX(0.,WU)
220      CONTINUE
        END IF
        CK12 = CK(1) + CK(2)
        CK13 = CK(1) + CK(3)
        CK14 = CK(1) + CK(4)
        CK23 = CK(2) + CK(3)
        CK24 = CK(2) + CK(4)
        CK34 = CK(3) + CK(4)
        CP12 = CP(1) + CP(2)
        CP23 = CP(2) + CP(3)
        CP34 = CP(3) + CP(4)
        CP14 = CP(1) + CP(4)
        CP1234 = CP12 + CP34
!
!       THESE ARE THE ELEMENTS OF THE 5 X 5 MATRIX
!
        A11 = A*CK12 + 2.*C*CP12
        A12 = -B*CK(2) + C*CP(2)
        A13 = B*CK(1) + C*CP(1)
        A15 = DLYZ*CK(1) + DLZY*CK(2) + C*CP12
        A22 = A*CK23 + 2.*C*CP23
        A24 = B*CK(3) + C*CP(3)
        A25 = DLYZ*CK(3) + DLZY*CK(2) + C*CP23
        A33 = A*CK14 + 2.*C*CP14
        A34 = -B*CK(4) + C*CP(4)
        A35 = DLYZ*CK(1) + DLZY*CK(4) + C*CP14
        A44 = A*CK34 + 2.*C*CP34
        A45 = DLYZ*CK(3) + DLZY*CK(4) + C*CP34
        A55 = 2.*(C*CP1234 - DLYZ*CK13 - DLZY*CK24)
!
!       THESE ARE THE ELEMENTS OF THE 5 X 1 SOURCE VECTOR
!
        DJ1 = (CMPLX(SIG(M,1,L),WE) - SIGH(M))*FP1(M)
        DJ2 = (CMPLX(SIG(M,2,L),WE) - SIGH(M))*FP2(M)
        DJ3 = (CMPLX(SIG(M,3,L),WE) - SIGH(M))*FP3(M)
        DJ4 = (CMPLX(SIG(M,4,L),WE) - SIGH(M))*FP4(M)
        SM1 = (0.,0.)
        SM2 = (0.,0.)
        SM3 = (0.,0.)
        SM4 = (0.,0.)
        SM5 = (0.,0.)
        IF(IDX .EQ. 2)GO TO 230
         DSS1 = 1. - SIGH(M)/CMPLX(SIG(M,1,L),WE)
         DSS2 = 1. - SIGH(M)/CMPLX(SIG(M,2,L),WE)
         DSS3 = 1. - SIGH(M)/CMPLX(SIG(M,3,L),WE)
         DSS4 = 1. - SIGH(M)/CMPLX(SIG(M,4,L),WE)
         DJ1 = CMPLX(0.,WU)*DSS1*FP1(M)
         DJ2 = CMPLX(0.,WU)*DSS2*FP2(M)
         DJ3 = CMPLX(0.,WU)*DSS3*FP3(M)
         DJ4 = CMPLX(0.,WU)*DSS4*FP4(M)
         SM1 = DSS1*(2.*FP5(M) - FP6(M)) + DSS2*FP6(M)
         SM2 = -DSS2*FP7(M) - DSS3*(2.*FP8(M) - FP7(M))
         SM3 = DSS1*(2.*FP5(M) - FP6(M)) + DSS4*FP6(M)
         SM4 = -DSS3*(2.*FP8(M) - FP7(M)) - DSS4*FP7(M)
         SM5 = -FP6(M)*(2.*DSS1 - DSS2 - DSS4) &
     &    - FP7(M)*(DSS2 - 2.*DSS3 + DSS4)
230     CONTINUE
        RU1 = C4*(DJ1 + DJ2) + DY4*SM1
        RU2 = C4*(DJ2 + DJ3) + DY4*SM2
        RU3 = C4*(DJ1 + DJ4) + DY4*SM3
        RU4 = C4*(DJ3 + DJ4) + DY4*SM4
        RU5 = C4*(DJ1 + DJ2 + DJ3 + DJ4) + DY4*SM5
!
!       ELEMENTS OF THE 4 X 4 CONDENSED MATRIX AND SOURCE
!
        A55I = 1./A55
        AA55(1,M,L) = A15*A55I
        AA55(2,M,L) = A25*A55I
        AA55(3,M,L) = A35*A55I
        AA55(4,M,L) = A45*A55I
        E11 = A11 - A15*AA55(1,M,L)
        G11 = A12 - A15*AA55(2,M,L)
        C11 = A13 - A15*AA55(3,M,L)
        H11 = -A15*AA55(4,M,L)
        E22 = A22 - A25*AA55(2,M,L)
        H22 = -A25*AA55(3,M,L)
        C22 = A24 - A25*AA55(4,M,L)
        E33 = A33 - A35*AA55(3,M,L)
        G33 = A34 - A35*AA55(4,M,L)
        E44 = A44 - A45*AA55(4,M,L)
        RC1 = RU1 - A15*(RU5/A55)
        RC2 = RU2 - A25*(RU5/A55)
        RC3 = RU3 - A35*(RU5/A55)
        RC4 = RU4 - A45*(RU5/A55)
!
!       THIS SECTION LOADS ELEMENT MATRICIES INTO GLOBAL MATRIX
!
        S(1,L1) = S(1,L1) + E11
        S(2,L1) = S(2,L1) + G11
        S(NBAN1,L1) = S(NBAN1,L1) + C11
        S(NBAND,L1) = S(NBAND,L1) + H11
        L2 = L1 + 1
        S(1,L2) = S(1,L2) + E22
        S(NDZ,L2) = S(NDZ,L2) + H22
        S(NBAN1,L2) = S(NBAN1,L2) + C22
        L3 = L1 + NDZ
        S(1,L3) = S(1,L3) + E33
        S(2,L3) = S(2,L3) + G33
        L4 = L3 + 1
        S(1,L4) = S(1,L4) + E44
        RM(L1) = RM(L1) + RC1
        RM(L2) = RM(L2) + RC2
        RM(L3) = RM(L3) + RC3
        RM(L4) = RM(L4) + RC4
        L1 = L1 + 1
200    CONTINUE
190   CONTINUE
!
!     INCORPORATE BOUNDARY CONDITIONS INTO GLOBAL EQUATION SYSTEM.
!     METHOD 1 OF HUEBNER (1983, P. 52) IS USED, REDUCING UNDERFLOWS.
!     ZERO CONDITIONS ARE USED AT THE MESH TOP FOR THE TM MODE AND AT
!     THE MESH BOTTOM FOR THE TE MODE.  NEUMANN CONDITIONS ELSEWHERE.
!
      NJI = NNODE - NDZ + 1
      IYN = 1
      DO 240 IT = 1,NJI,NDZ
       IF(IDX .EQ. 3)THEN
        S(1,IT) = (1.,0.)
        IF(IYN .EQ. 1)GO TO 250
         S(2,IT-1) = (0.,0.)
         S(NDZ,IT-NDZ+1) = (0.,0.)
         S(NBAN1,IT-NDZ) = (0.,0.)
         IF(IYN .EQ. 2)GO TO 250
         S(NBAND,IT-NBAN1) = (0.,0.)
250     CONTINUE
        S(2,IT) = (0.,0.)
        S(NDZ,IT) = (0.,0.)
        S(NBAN1,IT) = (0.,0.)
        S(NBAND,IT) = (0.,0.)
        RM(IT) = (0.,0.)
       END IF
       IF(IDX .EQ. 2)THEN
        IB = IT + NDZ1
        S(1,IB) = (1.,0.)
        IF(IYN .EQ. NODEY)GO TO 260
         S(2,IB) = (0.,0.)
         S(NDZ,IB) = (0.,0.)
         S(NBAN1,IB) = (0.,0.)
         IF(IYN .EQ. NODEY1)GO TO 260
         S(NBAND,IB) = (0.,0.)
260     CONTINUE
        S(2,IB-1) = (0.,0.)
        S(NDZ,IB-NDZ+1) = (0.,0.)
        IF(IYN .EQ. 1)GO TO 270
         S(NBAN1,IB-NDZ) = (0.,0.)
         S(NBAND,IB-NBAN1) = (0.,0.)
270     CONTINUE
        RM(IB) = (0.,0.)
       END IF
       IYN = IYN + 1
240   CONTINUE
!
!     SOLVE THE GLOBAL SYSTEM USING GAUSSIAN ELIMINATION
!
      CALL GAUSS(NNODE,NBAND)
      
      RETURN
END SUBROUTINE PARFLD
      
!-----------------------------------------------------------------------
      SUBROUTINE COAMP(WM,WP,NLYR)
!-----------------------------------------------------------------------
! OCCAM MT2D 3.0 Package
! Subroutine Revision 3.0, November 2006 (notional only to track future fixes)
! Original code written by Phil Wannamaker.
!     SUBROUTINE TO COMPUTE AMPLITUDE AND REFLECTION COEFFICIENTS FOR
!     ELECTRIC FIELDS PROPAGATING IN LAYERED EARTH.  COEFFICIENTS ARE
!     PASSED TO FUNCTION ZLFLD IN COMMON BLOCK BLKP.
!

      USE Fwd_mod

      IMPLICIT NONE

      COMPLEX AEXI,HTN(nLay-1),A2EX(nLay-1),CNUM,DEN,TL1
      integer nlyr,nlyr1,i,j,k,kk1
      real wm,wp,ra

      RA = 1.0E+30
      KL10 = CSQRT(-CMPLX(0.,WM)*CMPLX(1./RA,WP))
      DO 10 I = 1,NLYR
       KL1(I) = CSQRT(-CMPLX(0.,WM)*CMPLX(1./RL(I),WP))
10    CONTINUE
      IF(NLYR .GT. 1)THEN
       NLYR1 = NLYR - 1
       HL(1) = DDP(1)
       IF(NLYR1 .EQ. 1)GO TO 30
        DO 40 I = 2,NLYR1
         HL(I) = DDP(I) - DDP(I-1)
40      CONTINUE
30     CONTINUE
       DO 50 J = 1,NLYR1
        ARG(J) = (0.,1.)*KL1(J)*HL(J)
        IF(REAL(ARG(J)) .GT. 100.)GO TO 60
         AEX(J) = CEXP(-ARG(J))
         GO TO 70
60      CONTINUE
         AEX(J) = (0.,0.)
70      CONTINUE
        IF(REAL(ARG(J)) .GT. 40.)GO TO 80
         AEXI = 1./AEX(J)
         HTN(J) = (AEXI - AEX(J))/(AEXI + AEX(J))
         GO TO 81
80      CONTINUE
         HTN(J) = (1.,0.)
81      CONTINUE
        A2EX(J) = AEX(J)*AEX(J)
50     CONTINUE
       ZLA = WM/KL1(NLYR)
       DO 90 I = 1,NLYR1
        J = NLYR - I + 1
        ZL1 = WM/KL1(J-1)
        RFC(J) = (ZLA - ZL1)/(ZLA + ZL1)
        ZLA = ZL1*(ZLA + ZL1*HTN(J-1))/(ZL1 + ZLA*HTN(J-1))
90     CONTINUE
       ZL1 = WM/KL10
       COA(1) = (2.*ZLA/(ZLA + ZL1))/(1. + RFC(2)*A2EX(1))
       TL1 = COA(1)*(1. + RFC(2)*A2EX(1))
       RFC(1) = TL1 - 1.
       IF(NLYR1 .EQ. 1)GO TO 100
        DO 110 K = 2,NLYR1
         KK1 = K - 1
         DEN = 1. + RFC(K+1)*A2EX(K)
         CNUM = (1. + RFC(K))*AEX(KK1)
         COA(K) = COA(KK1)*CNUM/DEN
110     CONTINUE
100    CONTINUE
       COA(NLYR) = COA(NLYR1)*AEX(NLYR1)*(1. + RFC(NLYR))
       GO TO 120
      ELSE
       COA(1) = 2.*KL10/(KL1(1) + KL10)
       TL1 = COA(1)
       RFC(1) = (KL10 - KL1(1))/(KL10 + KL1(1))
      END IF
120   CONTINUE

      RETURN
END SUBROUTINE COAMP

!-----------------------------------------------------------------------
      COMPLEX FUNCTION ZLFLD(WM,Z,IFL,NLYR)
!-----------------------------------------------------------------------
! OCCAM MT2D 3.0 Package
! Subroutine Revision 3.0, November 2006 (notional only to track future fixes)
! Original code written by Phil Wannamaker.
!     FUNCTION COMPUTES ELECTRIC (IFL .EQ. 1)  OR MAGNETIC (IFL .EQ.
!     2) FIELD AT A DEPTH Z DUE TO UNIT-AMPLITUDE, DOWNGOING ELECTRIC
!     VECTOR AT SURFACE.  REFLECTION AND TRANSMISSION COEFFICIENTS
!     FOR THE LAYERING ARE PASSED FROM COAMP IN COMMON BLOCK BLKP.
!
      USE Fwd_mod

      IMPLICIT NONE

      COMPLEX K1WC,AR,AE,AEI,AR1,AR2
      integer ifl,nlyr,lzp,lzp1,i,im1
      real wm,z,zd1,zd2
      
      IF(Z .LT. 0.)THEN
       AR = (0.,1.)*KL10*Z
       IF(IFL .EQ. 1)ZLFLD = CEXP(-AR)*(1.+RFC(1)*CEXP(2.*AR))
!      IF(IFL .EQ. 1)ZLFLD = CMPLX(CDEXP(DCMPLX(-AR))*
!    &  (1.+RFC(1)*CDEXP(DCMPLX(2.*AR))))
       IF(IFL .EQ. 2)ZLFLD = (KL10/WM)*CEXP(-AR)* &
     &  (1.-RFC(1)*CEXP(2.*AR))
!      IF(IFL .EQ. 2)ZLFLD = CMPLX((KL10/WM)*CDEXP(DCMPLX(-AR))*
!    &  (1.-RFC(1)*CDEXP(DCMPLX(2.*AR))))
       GO TO 10
      ELSE
       IF(NLYR .GT. 1)THEN
        LZP = 1
        DO 30 I = 2,NLYR
         IM1 = I - 1
         IF(Z .GT. DDP(IM1))LZP = I
30      CONTINUE
        K1WC = (KL1(LZP)/WM)*COA(LZP)
        IF(LZP .EQ. NLYR)GO TO 40
         IF(REAL(ARG(LZP)) .GT. 40.)GO TO 50
          AR = (0.,1.)*KL1(LZP)*(Z-DDP(LZP))
          AE = CEXP(-AR)
          IF(IFL .EQ. 1)ZLFLD = COA(LZP)*(AE + RFC(LZP+1)/AE)*AEX(LZP)
          IF(IFL .EQ. 2)ZLFLD = K1WC*(AE - RFC(LZP+1)/AE)*AEX(LZP)
         GO TO 10
50       CONTINUE
         IF(LZP .EQ. 1)GO TO 60
          LZP1 = LZP - 1
          ZD1 = Z - DDP(LZP1)
          ZD2 = Z - 2.*DDP(LZP) + DDP(LZP1)
         GO TO 70
60       CONTINUE
          ZD1 = Z
          ZD2 = Z - 2.*DDP(LZP)
70       CONTINUE
         AR1 = (0.,1.)*KL1(LZP)*ZD1
         AE = (0.,0.)
         AEI = (0.,0.)
         IF(REAL(AR1) .GT. 100.)GO TO 80
          AE = CEXP(-AR1)
          AR2 = (0.,1.)*KL1(LZP)*ZD2
         IF(REAL(AR2) .LT. -100.)GO TO 80
          AEI = CEXP(AR2)
80       CONTINUE
          IF(IFL .EQ. 1)ZLFLD = COA(LZP)*(AE + RFC(LZP+1)*AEI)
          IF(IFL .EQ. 2)ZLFLD = K1WC*(AE - RFC(LZP+1)*AEI)
         GO TO 10
40      CONTINUE
        AR = (0.,1.)*KL1(LZP)*(Z-DDP(LZP-1))
        IF(REAL(AR) .GT. 100.)GO TO 90
         IF(IFL .EQ. 1)ZLFLD = COA(LZP)*CEXP(-AR)
         IF(IFL .EQ. 2)ZLFLD = K1WC*CEXP(-AR)
        GO TO 10
90      CONTINUE
         ZLFLD = (0.,0.)
        GO TO 10
       ELSE
        AR = (0.,1.)*KL1(1)*Z
        IF(REAL(AR) .GT. 100.)GO TO 100
         IF(IFL .EQ. 1)ZLFLD = COA(1)*CEXP(-AR)
         IF(IFL .EQ. 2)ZLFLD = (KL1(1)/WM)*COA(1)*CEXP(-AR)
        GO TO 10
100     CONTINUE
         ZLFLD = (0.,0.)
        GO TO 10
       END IF
      END IF
10    CONTINUE

      RETURN
END FUNCTION ZLFLD

!-----------------------------------------------------------------------
      SUBROUTINE GAUSS(NNODE,NBAND)
!-----------------------------------------------------------------------
! OCCAM MT2D 3.0 Package
! Subroutine Revision 3.0, November 2006 (notional only to track future fixes)
! Original code written by Phil Wannamaker.
!     ROUTINE TO REDUCE FINITE ELEMENT MATRIX.  GAUSSIAN ELIMINATION
!     WITHOUT PIVOTING BUT SPECIALIZED TO BANDED SYMMETRIC MATRICES
!     IS USED.  WRITTEN BY L. RIJO, 1977.  SUBSCRIPTING REVERSED IN
!     MATRIX S TO IMPROVE RUN-TIME, BY P. LUGAO, 1994.  ARITHMETIC IF
!     STATEMENTS REPLACED BY P. WANNAMAKER, 1994.
!
      USE Fwd_mod, NNODE_dont => NNODE
      IMPLICIT NONE
      
      COMPLEX D
      COMPLEX(8) C
      INTEGER NNODE,NBAND,N,NBLIM,L,I,J,K,M

      DO 10 N = 1,NNODE
       NBLIM = NBAND
       IF(N .GT. NNODE-NBAND+1)NBLIM = NNODE - N + 1
       IF(NBLIM .LT. 2)GO TO 10
       DO 20  L = 2,NBLIM
        I = N + L - 1
         C = S(L,N)/S(1,N)
         J = 0
         DO 40 K = L,NBAND
          J = J + 1
          S(J,I) = S(J,I) - C*S(K,N)
40       CONTINUE
         S(L,N) = C
20     CONTINUE
10    CONTINUE
!
!     FORWARD REDUCTION OF SOURCES
!
      DO 50 N = 1,NNODE
       NBLIM = NBAND
       IF(N .GT. NNODE-NBAND+1)NBLIM = NNODE - N + 1
       IF(NBLIM .LT. 2)GO TO 70
       DO 60 L = 2,NBLIM
        I = N + L - 1
        C = S(L,N)
        RM(I) = RM(I) - RM(N)*C
60     CONTINUE
70     CONTINUE
       D = 1./S(1,N)
       RM(N) = RM(N)*D
50    CONTINUE
!
!     BACK SUBSTITUTION OF SOURCES
!
      DO 80 M = 2,NNODE
       N = NNODE + 1 - M
       NBLIM = NBAND
       IF(N .GT. NNODE-NBAND+1)NBLIM = NNODE - N + 1
       IF(NBLIM .LT. 2)GO TO 80
       DO 90 L = 2,NBLIM
        K = N + L - 1
        C = S(L,N)
        RM(N) = RM(N) - RM(K)*C
90     CONTINUE
80    CONTINUE

      RETURN
END SUBROUTINE GAUSS

!-----------------------------------------------------------------------
      SUBROUTINE APPARMT
!-----------------------------------------------------------------------
! OCCAM MT2D 3.0 Package
! Subroutine Revision 3.0, November 2006 (notional only to track future fixes)
! Original code written by Phil Wannamaker.
!     CALCULATE AUXILIARY FIELDS, REAL AND IMAGINARY, AND APPARENT
!     RESISTIVITY, IMPEDANCE PHASE AND TIPPER.  AUXILIARY FIELDS
!     RETURNED FROM CALL TO AUXFLD.
!
      USE Fwd_mod

      IMPLICIT NONE

      COMPLEX FT,PF,HF,VF,TP,ZIMP
! miscellaneous stuff
      integer nzzn,in,i,nto,iii
      real pi,rad,w,wu

      PI = 3.1415927
      RAD = 180./PI
      W = 2.*PI*FREK(LFRQ)
      WU = W*PI*4.E-07
      NZZN = NTAO(1) + NTR
!ccKERRY: WRITE OUT X FIELD FOR ENTIRE MODEL:
      IF (loutfields) THEN
       write(*,*)'loutfields is true...'
       IF(IDX .EQ. 2)THEN
          WRITE(*,*) 'WRITING OUT TE ELECTRIC FIELD FOR ENTIRE MODEL...'
          WRITE(UNITTE,*) (nodey*nodez)
          WRITE(UNITTE,'(E11.4,1x,E11.4)') (RM(iii),iii=1,nodey*nodez)
          WRITE(UNITTE,*) (nodez)
          WRITE(UNITTE,'(E11.4,1x,E11.4)') (FPA(iii),iii=1,nodez)
       ELSE
          WRITE(*,*)'WRITING OUT TM MAGNETIC FIELD FOR ENTIRE MODEL...'
          WRITE(UNITTM,*)(nodey*nodez)
          WRITE(UNITTM,'(E11.4,1x,E11.4)')(RM(iii),iii=1,nodey*nodez)
          WRITE(UNITTM,*)(nodez)
          WRITE(UNITTM,'(E11.4,1x,E11.4)')(FPA(iii),iii=1,nodez)
       END IF
      END IF
!ccKERRY END            
      DO 10 IN = 1,NRC
       AFN(IN) = (0.,0.)
       AFT(IN) = (0.,0.)
       I = NRL(IN)
       CALL AUXFLD(I,IN)
       NTO = NTAO(I) + NTR
       FT = RM(NTA(I)) + FPA(NTO)
       PF = FT/FPA(NZZN)
       HF = AFT(IN)/FPB(NZZN)
       VF = AFN(IN)/FPB(NZZN)
       TP = VF/HF
       IF(IDX .EQ. 2)THEN
        ZIMP = FT/AFT(IN)
        zimpe(in,lfrq)=zimp
        RHTE(IN,LFRQ) = CABS(ZIMP)**2/WU
!   write(*,*)'RHTE(IN,LFRQ): ',RHTE(IN,LFRQ) 
        PHTE(IN,LFRQ) = RAD*ATAN2(AIMAG(ZIMP),REAL(ZIMP))
        KZTE(IN,LFRQ) = TP
       ELSE
        ZIMP = AFT(IN)/FT
        zimpm(in,lfrq)=zimp
        RHTM(IN,LFRQ) = CABS(ZIMP)**2/WU
        PHTM(IN,LFRQ) = RAD*ATAN2(-AIMAG(ZIMP),-REAL(ZIMP))
        TZTM(IN,LFRQ) = TP
       END IF
!       write(*,*) 'FT: ', FT
10    CONTINUE

      RETURN
END SUBROUTINE APPARMT

!-----------------------------------------------------------------------
      SUBROUTINE AUXFLD(I,IN)
!-----------------------------------------------------------------------
! OCCAM MT2D 3.0 Package
! Subroutine Revision 3.0, November 2006 (notional only to track future fixes)
! Original code written by Phil Wannamaker.
!     CALCULATE THE TOTAL FIELD NORMAL TO THE SLOPE AND THE TOTAL
!     FIELD PARALLEL TO THE SLOPE AND TRANSVERSE TO STRIKE.  THIS IS
!     DONE THROUGH THREE-POINT DIFFERENCE APPROXIMATIONS TO MAXWELL'S
!     EQUATIONS IN SUBROUTINE DERIV.  AUXILIARY FIELDS ARE RETURNED
!     TO MAIN IN ARRAYS AFN (NORMAL FIELD) AND AFT (TRANSVERSE 
!     FIELD).  FOR THE TE MODE , THE TRANSVERSE H-FIELD IS COMPUTED
!     IN A HORIZONTAL DIRECTION.  SIMILAR PROCEDURES ARE CARRIED OUT
!     FOR THE AUXILIARY JACOBIAN QUANTITIES IN AFND AND AFTD.
!
      USE Fwd_mod

      IMPLICIT NONE

      COMPLEX DERIV,DV,DH,UC
      DIMENSION IR(8)
      INTEGER I,IN,IR,NT,NCL,NCL1,NRW,NRW1,IH,IV,NPL,NPR,NPL1,NPR1,IH1,&
     &   IH2,IH3,IV1,IV2,IV3,IY1,IY2,IZ1,IZ2,ISB,NCSH,NHS,MSH,NTZ
      REAL PI,W,WU,WE,AN,SIGU

      PI = 3.1415927
      W = 2.*PI*FREK(LFRQ)
      WU = W*PI*4.E-07
      WE = W*8.85433E-12
!
!     DEFINE NODAL ROW AND COLUMN LOCATIONS OF FIELDPOINT
!
      NT = NTA(I)
      NCL = I
      NCL1 = NCL - 1
      NRW = NTAO(I) + NTR
      NRW1 = NRW - 1
!
!     COMPUTE SHIFT IN INDICES OF NODAL FIELDS ALONG STRIKE USED FOR
!     DETERMINING AUXILIARY FIELDS.  NODES USED ARE KEPT WITHIN THE
!     EARTH.  FOR BATHYMETRY, NODES ARE SHIFTED BY TWO INTO SEAWATER.
!
      IH = 0
      IV = -1
      AN = 0.
      NPL = 0
      NPR = 0
      IF(I .NE. 1)NPL = NPT(NRW,1,NCL1)
      IF(I .NE. NODEY)NPR = NPT(NRW,1,NCL)
!      write(*,*),'NPL,NPR: ',NPL,NPR
      NPARF(IN) = NPT(NRW,2,NCL)      !RESIS. LABEL UNDER RECEIVER
      IF(NPR .EQ. 0 .OR. NPR .EQ. ISEA)NPARF(IN) = NPT(NRW,4,NCL1)
      NPL1 = 0
      NPR1 = 0
      IF(NRW .NE. 1)THEN
       NPL1 = NPT(NRW1,3,NCL1)       
       NPR1 = NPT(NRW1,3,NCL)
!       write(*,*),'NPL1,NPR1: ',NPL1,NPR1
      END IF
      IF(I .EQ. 1 .OR. I .EQ. NODEY)GO TO 30
       CALL BNORM(NCL,NRW,NCL1,NRW1,AN,IR)       !SLOPE AT RECEIVER
       IF(NPL .EQ. 0)IH = 1                      !HORIZONTAL SHIFT TO
       IF(NPR .EQ. 0)IH = -1                     !AVOID AIR OR SEA
       IF(NPL .EQ. 0 .AND. NPR .EQ. 0)IH = 0
       IF(NRW .NE. 1)THEN                        !SEAFLOOR SHIFTS
        IF(NPL1 .EQ. ISEA)IH = -1
        IF(NPR1 .EQ. ISEA)IH = 1
        IF(NPL1 .EQ. ISEA .AND. NPR1 .EQ. ISEA)IH = 0
        IF(NPL1 .EQ. ISEA .OR. NPR1 .EQ. ISEA)IV = 1
       END IF
       GO TO 40
30    CONTINUE
      IF(I .EQ. 1 .OR. I .EQ. 2)THEN
       IH = 1
       IF(NRW .NE. 1 .AND. NPR1 .EQ. ISEA)IV = 1
      END IF
      IF(I .EQ. NODEY .OR. I .EQ. NODEY-1)THEN
       IH = -1
       IF(NRW .NE. 1 .AND. NPL1 .EQ. ISEA)IV = 1
      END IF
      AN = 0.
40    CONTINUE
!      IF(IV .EQ. 1)AN = 0.            !ASSUME SEAFLOOR LOCALLY FLAT
!kwk commented out for seafloor use slope  
!kwk  IF(IDX .EQ. 2)AN = 0.                      !TE HORIZ. & VERT.
      IF(IV .EQ. 1)NPARF(IN) = -1           !NO E-PRIM. TERM IN SEA
!
!     INDICES FOR FIELDPOINT LOCATIONS AND ELEMENT DIMENSIONS
!
      IH1 = NT + NDZ*(IH - 1)
      IH2 = NT + NDZ*IH
      IH3 = NT + NDZ*(IH + 1)
      IV1 = NT - IV - 1
      IV2 = NT - IV
      IV3 = NT - IV + 1
      IY1 = NCL1 + IH
      IY2 = NCL + IH
      IZ1 = NRW1 - IV
      IZ2 = NRW - IV
      ISB = IZ1 + (1 + IV)/2
!      write(*,*)'ih1,ih2,ih3,iv1,iv2,iv3: ',ih1,ih2,ih3,iv1,iv2,iv3 
!
!     COMPUTE VERTICAL AND HORIZONTAL FIELD DERIVATIVES (DV AND DH).
!     AVOID JUMPS IN SIG OR SIGH ONE NODE AWAY WITH 1ST-ORDER DIFF'S.
!
! DGM Jan 20, 2006 - added warning messages below.
!
! IV = -1 for land, +1 for sea
! ISB = the row that the site is at the bottom of (for sea)
! IZ1 = the row above the row containing the site.
! I,NCL = the column to the right of the site (the site sits above its left side)
! NPT is the parameter number in the regularization grid.  This number
!   will be 0=air, -1=water, or >0 for a real parameter.

      NCSH = 1
      IF(LSID .EQ. 2)NCSH = NODEY-1
      DV = -IV*(RM(IV2) - RM(IV1))/DLZ(ISB)       !DEFAULT 1ST-ORDER
!      write(*,*) 'IZ1: ', IZ1
      IF(IZ1 .LT. 1) THEN
        write(*,*) 'WARNING: Site ', IN, ' using linear derivative ' &
     & , 'because it is at the top of the model.'
        GO TO 50
      END IF
!       write(*,*) 'I, NODEY', I, NODEY
      IF(I .LT. NODEY) THEN
!        write(*,*) '(I<NodeY) ISB,NCL,IV,NPT...', ISB, NCL, IV
!     & , NPT(ISB,2,NCL), NPT(ISB-IV,2,NCL)
        IF(NPT(ISB,2,NCL) .NE. NPT(ISB-IV,2,NCL)) THEN
          IF (IV .gt. 0) THEN       ! Land site
            write(*,*) 'WARNING: Site ', IN, ' using linear derivative ' &
     &  , 'because triangle #2 on two blocks below are in different ' &
     &  , 'model blocks.  Chg model.'
          ELSE                      ! Sea site
            write(*,*) 'WARNING: Site ', IN, ' using linear derivative ' &
     &  , 'because triangle #2 on two blocks above are in different ' &
     &  , 'model blocks.  Chg model.'
          END IF
          GO TO 50
        END IF
      END IF
      IF(I .GT. 1) THEN
!        write(*,*) '(I>1) ISB,NCL,IV,NPT etc', ISB, NCL, IV
!     & , NPT(ISB,4,NCL-1), NPT(ISB-IV,4,NCL-1)
        IF(NPT(ISB,4,NCL-1) .NE. NPT(ISB-IV,4,NCL-1)) THEN
          IF (IV .gt. 0) THEN       ! Land site
            write(*,*) 'WARNING: Site ', IN, ' using linear derivative ' &
     &  , 'because triangle #4 on two blocks below are in different ' &
     &  , 'model blocks.  Chg model.'
          ELSE                      ! Sea site
            write(*,*) 'WARNING: Site ', IN, ' using linear derivative ' &
     &  , 'because triangle #4 on two blocks above are in different ' &
     &  , 'model blocks.  Chg model.'
          END IF
          GO TO 50
        END IF
      END IF
!       write(*,*) 'ISB,NCSH,IV,NPT...', ISB,NCSH,IV
!     & , NPT(ISB,3,NCSH), NPT(ISB-IV,1,NCSH)

! DGM Jan 20, 2006
! The logic below is looking at the parameters on the very left of
!   the model because that is where the 1D response is calculated ... maybe.
!   Or maybe that's just left-overs from some old process...
!   The purpose of the last condition is to avoid a change in the 
!   background conductivity model.
! This IF doesn't work for models with bathymetry where the lefthand
!   side of the model is higher than some of the sites.  In the case
!   where a site is deeper than the lefthand side AND the two blocks 
!   above the site are in different regularization blocks, then this
!   IF is triggered and the parabolic derivative is skipped.  This 
!   causes the occam inversion to diverge.
!      IF(NPT(ISB,3,NCSH) .NE. NPT(ISB-IV,1,NCSH)) THEN
!       write(*,*) 'WARNING: Site ', IN, ' using linear derivative '
!       GO TO 50 !AVOID SIGH
!      END IF
      
      DV = DERIV(RM(IV1),RM(IV2),RM(IV3),DLZ(IZ1),DLZ(IZ2),IV)
!       write(*,*) 'RM(IV1),RM(IV2),RM(IV3),DLZ(IZ1),DLZ(IZ2),IV: '
!     & ,RM(IV1),RM(IV2),RM(IV3),DLZ(IZ1),DLZ(IZ2),IV
50    CONTINUE
      DH = (0.,0.)
      IF(I .EQ. 1 .OR. I .EQ. NODEY)GO TO 52
       IF(IH .NE. 0)THEN
        NHS = 1
        MSH = 1
        IF(IH .EQ. 1)NHS = 0
        IF(IH .EQ. 1)MSH = -1
        IF(NPT(ISB,2+IV,NCL-NHS) .NE. NPT(ISB,2+IV,NCL+1-3*NHS))THEN
         NTZ = NT + NDZ*(1-2*NHS)
         DH = MSH*(RM(NT) - RM(NTZ))/DLY(NCL-NHS) !DEFAULT 1ST-ORDER
         GO TO 52
        END IF
       END IF
       DH = DERIV(RM(IH1),RM(IH2),RM(IH3),DLY(IY1),DLY(IY2),IH)
52    CONTINUE
!
!     TE MODE AUXILIARY FIELDS
!
      IF(IDX .EQ. 2)THEN
!kwk commented out for using seafloor slope       
!kwk   AN = 0.
       UC = -CMPLX(0.,WU)
       AFN(IN) = -(FPB(NRW) + DV/UC)*SIN(AN) - (DH/UC)*COS(AN)
       AFT(IN) = (FPB(NRW) + DV/UC)*COS(AN) - (DH/UC)*SIN(AN)
       CSIGU(IN) = UC                  !SAVED FOR PARJAC AND AUXADJ
!
!     TM MODE AUXILIARY FIELDS
!
      ELSE
       SIGU = SIG(ISB,2,NCL)
       IF(I .EQ. NODEY) SIGU = SIG(ISB,4,NCL-1)
       UC = CMPLX(SIGU,WE)
!       AN = 0.
       AFN(IN) = (-(SIGH(ISB)*FPB(NRW) +DV)*SIN(AN) - DH*COS(AN))/UC
       AFT(IN) = ((SIGH(ISB)*FPB(NRW) + DV)*COS(AN) - DH*SIN(AN))/UC
!       write(*,*) 'IN, AN, AFN, AFT: ',IN, AN, AFN(IN), AFT(IN)
!       write(*,*) 'SIGH(ISB), FPB(NRW), DV, UC: ', SIGH(ISB), FPB(NRW), 
!     &  DV, UC
       CSIGU(IN) = UC                  !SAVED FOR PARJAC AND AUXADJ
      END IF
      
      RETURN
END SUBROUTINE AUXFLD

!-----------------------------------------------------------------------
      SUBROUTINE BNORM(NCL,NRW,NCL1,NRW1,AN,IR)
!-----------------------------------------------------------------------
! OCCAM MT2D 3.0 Package
! Subroutine Revision 3.0, November 2006 (notional only to track future fixes)
! Original code written by Phil Wannamaker.
!     CALCULATE DOWNWARD BISECTOR FROM FIELD NODE AT BREAK IN TOPO-
!     GRAPHY.  AR, AL CLOCKWISE FROM UP, AN CLOCKWISE FROM DOWN.
!
      USE Fwd_mod

      IMPLICIT NONE

      COMPLEX TAL,TAR
      INTEGER NCL,NRW,NCL1,NRW1,IR,I
      REAL AN,SM,AL,AR
      DIMENSION IR(8)
      
!     IDENTIFY RESISTIVITY CODES OF TRIANGULAR ELEMENTS SURROUNDING
!     FIELD NODE.  UPPERMOST LEFT OCTANT IS IR(1), RIGHT IS IR(8).
!
      IR(3) = NPT(NRW,1,NCL1)
      IR(4) = NPT(NRW,4,NCL1)
      IR(5) = NPT(NRW,2,NCL)
      IR(6) = NPT(NRW,1,NCL)
      IF(NRW .EQ. 1)GO TO 10
       IR(1) = NPT(NRW1,4,NCL1)
       IR(2) = NPT(NRW1,3,NCL1)
       IR(7) = NPT(NRW1,3,NCL)
       IR(8) = NPT(NRW1,2,NCL)
      GO TO 20
10    CONTINUE
       IR(1) = 0
       IR(2) = 0
       IR(7) = 0
       IR(8) = 0
20    CONTINUE
      DO 40 I = 1,8
       IF (IR(I) .EQ. ISEA) IR(I) = 0
40    CONTINUE
!
!     DETERMINE SLOPE TO LEFT OF FIELD NODE
!
      SM = 1.0E-05
      TAL = 100.*CMPLX(1.,-SM)
      IF (NRW .EQ. 1) GO TO 50
      TAL = DLZ(NRW1)*CMPLX(1.,-SM)
      IF (IR(1).EQ.0 .AND. IR(2).NE.0)TAL = CMPLX(DLZ(NRW1),-DLY(NCL1))
50    CONTINUE
      IF (IR(2).EQ.0 .AND. IR(3).NE.0)TAL = DLY(NCL1)*CMPLX(SM,-1.)
      IF (IR(3).EQ.0 .AND. IR(4).NE.0)TAL = CMPLX(-DLZ(NRW),-DLY(NCL1))
      IF (IR(4).EQ.0 .AND. IR(5).NE.0)TAL = DLZ(NRW)*CMPLX(-1.,-SM)
!
!     DETERMINE SLOPE TO RIGHT OF FIELD NODE
!
      TAR = 100.*CMPLX(1.,SM)
      IF (NRW .EQ. 1) GO TO 60
      TAR = DLZ(NRW1)*CMPLX(1.,SM)
      IF (IR(8).EQ.0 .AND. IR(7).NE.0)TAR = CMPLX(DLZ(NRW1),DLY(NCL))
60    CONTINUE
      IF (IR(7).EQ.0 .AND. IR(6).NE.0)TAR = DLY(NCL)*CMPLX(SM,1.)
      IF (IR(6).EQ.0 .AND. IR(5).NE.0)TAR = CMPLX(-DLZ(NRW),DLY(NCL))
      IF (IR(5).EQ.0 .AND. IR(4).NE.0)TAR = DLZ(NRW)*CMPLX(-1.,SM)
      AL = ATAN2(AIMAG(TAL),REAL(TAL))
      AR = ATAN2(AIMAG(TAR),REAL(TAR))
      AN = (AL + AR)/2.
!      write(*,*),'BNORM SLOPE:  ',AN*180/3.16

      RETURN
END SUBROUTINE BNORM

!-----------------------------------------------------------------------
      COMPLEX FUNCTION DERIV(R1,R2,R3,D1,D2,IF)
!-----------------------------------------------------------------------
! OCCAM MT2D 3.0 Package
! Subroutine Revision 3.0, November 2006 (notional only to track future fixes)
! Original code written by Phil Wannamaker.
!     FUNCTION ESTIMATES DERIVATIVE AT LEFT END, MIDDLE OR RIGHT END
!     OF A COMPLEX THREE-POINT SEQUENCE BY FITTING THE PARABOLA
!     R = A + B*D + C*D**2.  PARABOLA COEFFICIENTS ARE CALCULATED
!     FROM AMPLITUDES AND SEPARATIONS OF THE THREE POINTS.
!
      IMPLICIT NONE
      COMPLEX R1,R2,R3,B,C
      INTEGER IF
      REAL D1,D2,DEN

      DEN = D1*D2*D2 + D1*D1*D2
!     A = R2
      B = ((R3-R2)*D1*D1 - (R1-R2)*D2*D2)/DEN
      C = ((R3-R2)*D1 + (R1-R2)*D2)/DEN
      IF(IF)10,20,30
!     DERIVATIVE AT LEFT POINT
10    CONTINUE
       DERIV = B - 2.*C*D1
      GO TO 40
!     DERIVATIVE AT MIDDLE POINT
20    CONTINUE
       DERIV = B
      GO TO 40
!     DERIVATIVE AT RIGHT POINT
30    CONTINUE
       DERIV = B + 2.*C*D2
40    CONTINUE

      RETURN
END FUNCTION DERIV

!-----------------------------------------------------------------------
      SUBROUTINE RECIJAC(NRES)
!-----------------------------------------------------------------------
! OCCAM MT2D 3.0 Package
! Subroutine Revision 3.0, November 2006 (notional only to track future fixes)
! Original code written by Phil Wannamaker.
!     SUBROUTINE TO CALCULATE DATA-PARAMETER JACOBIANS USING
!     RECIPROCITY.  PARALLEL FIELD QUANTITIES COMPUTED FOR TE MODE
!     ONLY BY LOADING UNIT SOURCES AT MT RECEIVER POINTS.  TM MODE
!     JACOBIAN ALONG STRIKE IS ZERO.  TRANSVERSE AND VERTICAL FIELD
!     QUANTITIES COMPUTED BY DEFINING A DISTRIBUTED SOURCE AT MT
!     FIELD POINT AND NEIGHBOURING POINTS FOR FIRST OR SECOND-ORDER
!     DIFFERENCING.
!
      USE Fwd_mod

      DO 10 KP = 1,NRES
       DO 11 I = 1,NRC
        FLDJ(I,KP) = (0.,0.)
        FLDJH(I,KP) = (0.,0.)
        FLDJV(I,KP) = (0.,0.)
11     CONTINUE 
10    CONTINUE 
!      
!     ALONG-STRIKE COMPONENT COMPUTED BY LOADING ARRAY OF RECIPROCAL
!     SOURCES AT MT RECEIVERS AND SOLVING IN GAUSSD.
!
      DO 15 JY = 1,NNODE
       DO 16 IY = 1,NODEY
        RJ(IY,JY) = (0.,0.)
16     CONTINUE  
15    CONTINUE  
      IX = 0
      DO 17 ILD = 1,NRC
       I = NTA(NRL(ILD))
       IF(IDX .EQ. 2 .OR. NPARF(ILD) .EQ. -1)RJ(ILD,I) = (1.,0.)
17    CONTINUE
      NBAND=NDZ+2
      CALL GAUSSD(NNODE,NBAND,NRC)
      CALL PARJAC(IX)
!
!     CALCULATE AUXILIARY FIELD JACOBIANS, USING DISTRIBUTED SOURCES.
!     FOR TE MODE, CALCULATE BOTH HORIZONTAL AND VERTICAL AUXILIARY
!     JACOBIANS.  FOR TM MODE, CALCULATE ONLY JACOBIAN OF AUXILIARY
!     ELECTRIC FIELD ALONG SLOPE.
!
      DO 20 IX = 1,2
       DO 21 JY = 1,NNODE
        DO 22 IY = 1,NODEY
         RJ(IY,JY) = (0.,0.)
22      CONTINUE  
21     CONTINUE  
       IF(IX .EQ. 2 .AND. IDX .EQ. 3)GO TO 20              !NO EZ
       DO 23 ILD = 1,NRC
        I = NRL(ILD)
        CALL AUXADJ(I,ILD,IX)     !DEFINE DISTRIBUTED SOURCE VECTORS
23     CONTINUE
       CALL GAUSSD(NNODE,NBAND,NRC)
       CALL PARJAC(IX)
20    CONTINUE
!
!     CONVERT E AND H-FIELD JACOBIANS TO THOSE OF LOG APP RHO, PHASE
!     AND TIPPER.  JACOBIANS ARE WITH RESPECT TO LOG RHO.
!
      CALL APPRMTJ

      RETURN
END SUBROUTINE RECIJAC

!-----------------------------------------------------------------------
      SUBROUTINE PARJAC(IX)
!-----------------------------------------------------------------------
! OCCAM MT2D 3.0 Package
! Subroutine Revision 3.0, November 2006 (notional only to track future fixes)
! Original code written by Phil Wannamaker.
!     THIS ROUTINE COMPUTES THE DERIVATIVE OF THE FIELD PARALLEL TO
!     STRIKE WITH RESPECT TO THE CONDUCTIVITY OF A POLYGONAL REGION
!     IN THE MESH USING RECIPROCITY.  THIS IS ACCOMPLISHED FOLLOWING
!     SOLUTION OF UNIT OR DISTRIBUTED RECIPROCAL SOURCES BY WEIGHTING
!     WITH ELEMENT PROPERTIES AND TOTAL FIELDS.
!
      USE Fwd_mod


      COMPLEX CPD(4),CPD12,CPD23,CPD34,CPD14,CPD1T4,AD11,AD12,AD13,&
     &AD15,AD22,AD24,AD25,AD33,AD34,AD35,AD44,AD45,AD55,ED11,GD11,&
     &CD11,HD11,ED22,HD22,CD22,ED33,GD33,ED44,&
     &FT1,FT2,FT3,FT4,RJC1,RJC2,RJC3,RJC4
      COMPLEX DSIG,CKD(4),CKD12,CKD13,CKD14,CKD23,CKD24,CKD34

      PI = 3.1415927
      W = 2.*PI*FREK(LFRQ)
      WE = W*8.85433E-12
!
!     LOOP OVER ELEMENTS, CONSTRUCTING DERIVATIVE SOURCE TERMS
!
      NODEY1 = NODEY - 1
      NODEZ1 = NODEZ - 1
      L1 = 0
      DO 140 L = 1,NODEY1
       L1 = L1 + 1
       DO 150 M = NTQ,NODEZ1
        L2 = L1 + 1
        L3 = L1 + NDZ 
        L4 = L3 + 1
!
!       CALCULATE THE APPROPRIATE COEFFICIENTS FOR INSERTION INTO
!       THE SOURCE. IF TRIANGLES WITHIN RECTANGULAR ELEMENT ARE NOT
!       ALL THE SAME, SUM THEIR CONTRIBUTIONS SEPARATELY (NPE = 4). 
!
        NPT1 = NPT(M,1,L)
        NPE = 1
        DO 120 IC = 2,4
         IF(NPT(M,IC,L) .NE. NPT1)NPE = 4
120     CONTINUE
        IF(NPE .EQ. 1 .AND. NPT1 .EQ. 0)GO TO 131
        DLZY = DLZ(M)/(2.*DLY(L))
        DLYZ = DLY(L)/(2.*DLZ(M))
        A = -(DLZY + DLYZ)/2.
        B = (DLZY - DLYZ)/2.
        C = DLY(L)*DLZ(M)/48.
!
!       CALCULATE THE APPROPRIATE COEFFICIENTS FOR INSERTION INTO
!       SOURCE VECTOR ACCORDING TO WHETHER TE OR TM MODE SOLVED.
!
        DO 130 IE = 1,NPE
         NPTE = NPT(M,IE,L)
         IF(NPE .EQ. 1)THEN                 !HOMOGENEOUS RECTANGLE
          IF(IDX .EQ. 2)THEN                !TE MODE
           DO 132 I = 1,4
            CKD(I) = CMPLX(0.,0.)
            CPD(I) = CMPLX(1.,0.)
132        CONTINUE
          ELSE                              !TM MODE
           DO 133 I = 1,4
            DSIG = 1./CMPLX(SIG(M,I,L),WE)
            CKD(I)= DSIG*DSIG
            CPD(I) = CMPLX(0.,0.)
133        CONTINUE
          END IF
         ELSE                               !INHOMOGENEOUS RECTANGLE
          DO 135 I = 1,4
           CKD(I) = 0.
           CPD(I) = 0.
135       CONTINUE
!
!         CALCULATE THE APPROPRIATE COEFFICIENTS FOR INSERTION INTO
!         SOURCE VECTOR ACCORDING TO WHETHER TE OR TM MODE SOLVED.

          IF(IDX .EQ. 2)THEN
           CKD(IE) = CMPLX(0.,0.)
           CPD(IE) = CMPLX(1.,0.)
          ELSE
           DSIG = 1./CMPLX(SIG(M,IE,L),WE)
           CKD(IE)= DSIG*DSIG
           CPD(IE) = CMPLX(0.,0.)
          END IF
         END IF
         CKD12 = CKD(1) + CKD(2)
         CKD13 = CKD(1) + CKD(3)
         CKD14 = CKD(1) + CKD(4)
         CKD23 = CKD(2) + CKD(3)
         CKD24 = CKD(2) + CKD(4)
         CKD34 = CKD(3) + CKD(4)
         CPD12 = CPD(1) + CPD(2)
         CPD23 = CPD(2) + CPD(3)
         CPD34 = CPD(3) + CPD(4)
         CPD14 = CPD(1) + CPD(4)
         CPD1T4 = CPD12 + CPD34
!
!        ELEMENTS OF THE 5 X 5 DIFFERENTIATED MATRIX
!
         AD11 = A*CKD12 + 2.*C*CPD12
         AD12 = -B*CKD(2) + C*CPD(2)
         AD13 = B*CKD(1) + C*CPD(1)
         AD15 = DLYZ*CKD(1) + DLZY*CKD(2) + C*CPD12
         AD22 = A*CKD23 + 2.*C*CPD23
         AD24 = B*CKD(3) + C*CPD(3)
         AD25 = DLYZ*CKD(3) + DLZY*CKD(2) + C*CPD23
         AD33 = A*CKD14 + 2.*C*CPD14
         AD34 = -B*CKD(4) + C*CPD(4)
         AD35 = DLYZ*CKD(1) + DLZY*CKD(4) + C*CPD14
         AD44 = A*CKD34 + 2.*C*CPD34
         AD45 = DLYZ*CKD(3) + DLZY*CKD(4) + C*CPD34
         AD55 = 2.*(C*CPD1T4 - DLYZ*CKD13 - DLZY*CKD24)

!        ELEMENTS OF THE CONDENSED DIFFERENTIATED MATRIX

         DO 170 IA = 1,4
          AAD55(IA) = AD55*AA55(IA,M,L)
170      CONTINUE
         ED11 = AD11 - 2.*AD15*AA55(1,M,L) + AAD55(1)*AA55(1,M,L)
         GD11 = AD12 - AD15*AA55(2,M,L) - AD25*AA55(1,M,L) &
     &               + AAD55(1)*AA55(2,M,L)
         CD11 = AD13 - AD15*AA55(3,M,L) - AD35*AA55(1,M,L) &
     &               + AAD55(1)*AA55(3,M,L)
         HD11 =  - AD15*AA55(4,M,L) - AD45*AA55(1,M,L) &
     &               + AAD55(1)*AA55(4,M,L)
         ED22 = AD22 - 2.*AD25*AA55(2,M,L) + AAD55(2)*AA55(2,M,L)
         HD22 =  - AD25*AA55(3,M,L) - AD35*AA55(2,M,L) &
     &               + AAD55(2)*AA55(3,M,L)
         CD22 = AD24 - AD25*AA55(4,M,L) - AD45*AA55(2,M,L) &
     &               + AAD55(2)*AA55(4,M,L)
         ED33 = AD33 - 2.*AD35*AA55(3,M,L) + AAD55(3)*AA55(3,M,L)
         GD33 = AD34 - AD35*AA55(4,M,L) - AD45*AA55(3,M,L) &
     &               + AAD55(3)*AA55(4,M,L)
         ED44 = AD44 - 2.*AD45*AA55(4,M,L) + AAD55(4)*AA55(4,M,L)
!
!        MULTIPLYING THE CONDENSED DIFFERENTIATED MATRIX BY THE TOTAL
!        FIELDS ALONG STRIKE AT THE CORNERS OF THE RECTANGULAR
!        ELEMENT YIELDS THE CONDENSED SOURCE. EVALUATION OF SOURCE TERMS
!        UTILIZES FIELDS FROM FORWARD PROBLEM AND RECIPROCITY.
! 
! DGM Oct 9, 2006 - NPTE may also equal -1 (seawater)
         IF (NPTE .NE. 0 .and. NPTE .NE. ISEA) THEN
          FT1 = RM(L1) + FPA(M)  
          FT2 = RM(L2) + FPA(M+1)
          FT3 = RM(L3) + FPA(M)
          FT4 = RM(L4) + FPA(M+1)
          RJC1 = ED11*FT1 + GD11*FT2 + CD11*FT3 + HD11*FT4
          RJC2 = GD11*FT1 + ED22*FT2 + HD22*FT3 + CD22*FT4
          RJC3 = CD11*FT1 + HD22*FT2 + ED33*FT3 + GD33*FT4
          RJC4 = HD11*FT1 + CD22*FT2 + GD33*FT3 + ED44*FT4
!
          DO 107 KR = 1,NRC
!
           IF(IX .EQ. 0)THEN           !TE PARALLEL JACOBIAN
            FLDJ(KR,NPTE) = FLDJ(KR,NPTE) + RJC1*RJ(KR,L1) &
     &       + RJC2*RJ(KR,L2) + RJC3*RJ(KR,L3) + RJC4*RJ(KR,L4)
           ELSE IF(IX .EQ. 1)THEN      !HORIZONTAL AUX. JACOBIAN
            FLDJH(KR,NPTE) = FLDJH(KR,NPTE) + RJC1*RJ(KR,L1) &
     &       + RJC2*RJ(KR,L2) + RJC3*RJ(KR,L3) + RJC4*RJ(KR,L4)
           ELSE IF(IX .EQ. 2)THEN      !VERTICAL AUX. JACOBIAN
            FLDJV(KR,NPTE) = FLDJV(KR,NPTE) + RJC1*RJ(KR,L1) &
     &       + RJC2*RJ(KR,L2) + RJC3*RJ(KR,L3) + RJC4*RJ(KR,L4)
           END IF
!
107       CONTINUE
         END IF
130     CONTINUE
131     CONTINUE
        L1 = L1 + 1
150    CONTINUE
140   CONTINUE
!
!     FOR TM MODE, ADD E-FIELD TO THE ALONG-SLOPE ELECTRIC JACOBIAN
!     WHEN RECEIVER POINT IS IN THE PARAMETER.
!
      IF(IDX .EQ. 3 .AND. IX .EQ. 1) THEN
       DO 177 KS = 1,NRC
        IF(NPARF(KS) .GT. 0)THEN
         FLDJH(KS,NPARF(KS)) = FLDJH(KS,NPARF(KS)) &
     &                         - AFT(KS)/CSIGU(KS)
        END IF
177    CONTINUE
      END IF
      
      RETURN
END SUBROUTINE PARJAC

!-----------------------------------------------------------------------
      SUBROUTINE GAUSSD(NNODE,NBAND,NSRC)
!-----------------------------------------------------------------------
! OCCAM MT2D 3.0 Package
! Subroutine Revision 3.0, November 2006 (notional only to track future fixes)
! Original code written by Phil Wannamaker.
      use Fwd_mod, NNODE_dont => NNODE

      COMPLEX C,D
      
      DO 50 N = 1,NNODE
       NBLIM = NBAND
       IF(N .GT. NNODE-NBAND+1)NBLIM = NNODE - N + 1
       IF(NBLIM .LT. 2)GO TO 70
        DO 60 L = 2,NBLIM
         C = S(L,N)
         I = N + L - 1
         DO 80 KP = 1,NSRC
          RJ(KP,I) = RJ(KP,I) - C*RJ(KP,N)
80       CONTINUE
60      CONTINUE
70     CONTINUE
       D = 1./S(1,N)
       DO 90 KP = 1,NSRC
        RJ(KP,N) = RJ(KP,N)*D
90     CONTINUE
50    CONTINUE
!
      DO 100 M = 2,NNODE
       N = NNODE + 1 - M
       NBLIM = NBAND
       IF(N .GT. NNODE-NBAND+1)NBLIM = NNODE - N + 1
       IF(NBLIM .LT. 2)GO TO 100
       DO 110 L = 2,NBLIM
        C = S(L,N)
        K = N + L - 1
        DO 120 KP = 1,NSRC
         RJ(KP,N) = RJ(KP,N) - C*RJ(KP,K)
120     CONTINUE
110    CONTINUE
100   CONTINUE

      RETURN
END SUBROUTINE GAUSSD

!-----------------------------------------------------------------------
      SUBROUTINE AUXADJ(I,ILD,IX)
!-----------------------------------------------------------------------
! OCCAM MT2D 3.0 Package
! Subroutine Revision 3.0, November 2006 (notional only to track future fixes)
! Original code written by Phil Wannamaker.
!     LOADS THE AUXILIARY FIELD JACOBIAN SOURCE VECTORS ACCORDING TO
!     RECIPROCITY.  FOR EACH RECEIVER, SUBROUTINE DRVADJ COMPUTES
!     PARABOLA COEFFICIENTS A, B, AND C WHICH WILL BE LOADED INTO
!     SOURCE VECTOR RJ(,) TO BE SOLVED IN SUBROUTINE GAUSSD.  COEFF-
!     ICIENT VALUES DEPEND UPON TOPOGRAPHIC SLOPE AND ELEMENT PROPER-
!     TIES.  FOR TM MODE, PROJECTED VERTICAL AND HORIZONTAL DERIVA-
!     TIVES ARE ADDED TO GIVE E-FIELD ALONG SLOPE.
!
      USE Fwd_mod

      COMPLEX UC
      DIMENSION IR(8)

!     DEFINE NODAL ROW AND COLUMN LOCATIONS OF FIELDPOINT
!
      NT = NTA(I)
      NCL = I
      NCL1 = NCL - 1
      NRW = NTAO(I) + NTR
      NRW1 = NRW - 1
      UC = CSIGU(ILD)
!
!     COMPUTE SHIFT IN INDICES OF NODAL FIELDS ALONG STRIKE USED FOR
!     DETERMINING AUXILIARY FIELDS.  NODES USED ARE KEPT WITHIN THE
!     EARTH.  FOR BATHYMETRY, NODES ARE SHIFTED BY TWO INTO SEAWATER.
!
      IH = 0
      IV = -1
      AN = 0.
      NPL = 0
      NPR = 0
      IF(I .NE. 1)NPL = NPT(NRW,1,NCL1)
      IF(I .NE. NODEY)NPR = NPT(NRW,1,NCL)
      NPL1 = 0
      NPR1 = 0
      IF(NRW .NE. 1)THEN
       NPL1 = NPT(NRW1,3,NCL1)
       NPR1 = NPT(NRW1,3,NCL)
      END IF
      IF(I .EQ. 1 .OR. I .EQ. NODEY)GO TO 30
       CALL BNORM(NCL,NRW,NCL1,NRW1,AN,IR)       !SLOPE AT RECEIVER
       IF(NPL .EQ. 0)IH = 1                      !HORIZONTAL SHIFT TO
       IF(NPR .EQ. 0)IH = -1                     !AVOID AIR OR SEA
       IF(NPL .EQ. 0 .AND. NPR .EQ. 0)IH = 0
       IF(NRW .NE. 1)THEN                        !SEAFLOOR SHIFTS
        IF(NPL1 .EQ. ISEA)IH = -1
        IF(NPR1 .EQ. ISEA)IH = 1
        IF(NPL1 .EQ. ISEA .AND. NPR1 .EQ. ISEA)IH = 0
        IF(NPL1 .EQ. ISEA .OR. NPR1 .EQ. ISEA)IV = 1
       END IF
       GO TO 40
30    CONTINUE
      IF(I .EQ. 1 .OR. I .EQ. 2)THEN
       IH = 1
       IF(NRW .NE. 1 .AND. NPR1 .EQ. ISEA)IV = 1
      END IF
      IF(I .EQ. NODEY .OR. I .EQ. NODEY-1)THEN
       IH = -1
       IF(NRW .NE. 1 .AND. NPL1 .EQ. ISEA)IV = 1
      END IF
      AN = 0.
40    CONTINUE
!      IF(IV .EQ. 1)AN = 0.             !ASSUME SEAFLOOR LOCALLY FLAT
!kwk commented out for using seafloor slope:
!kwk      IF(IDX .EQ. 2)AN = 0.                      !TE HORIZ. & VERT.
!
!     INDICES FOR FIELDPOINT LOCATIONS AND ELEMENT DIMENSIONS
!
      IH1 = NT + NDZ*(IH - 1)
      IH2 = NT + NDZ*IH
      IH3 = NT + NDZ*(IH + 1)
      IV1 = NT - IV - 1
      IV2 = NT - IV
      IV3 = NT - IV + 1
      IY1 = NCL1 + IH
      IY2 = NCL + IH
      IZ1 = NRW1 - IV
      IZ2 = NRW - IV
      ISB = IZ1 + (1 + IV)/2
!
      IF(IX .EQ. 1)THEN                !CROSS-STRIKE FIELD JACOBIAN
!
!      VERTICAL DERIVATIVE COEFFICIENTS A, B AND C FOR HORIZONTAL
!      AUXILIARY JACOBIAN.
!
       A = -1./DLZ(IZ1)                 !AVOID JUMPS IN SIG ONE NODE
       B = 1./DLZ(IZ1)                  !AWAY WITH 1ST-ORDER DIFF'S.
       C = 0.
       IF(IZ1 .LT. 1)GO TO 50
        IF(I .LT. NODEY) THEN
         IF(NPT(ISB,2,NCL) .NE. NPT(ISB-IV,2,NCL))GO TO 50
        END IF
        IF(I .GT. 1) THEN
         IF(NPT(ISB,4,NCL-1) .NE. NPT(ISB-IV,4,NCL-1))GO TO 50
        END IF
!
        CALL DRVADJ(DLZ(IZ1),DLZ(IZ2),IV,A,B,C)
!
50     CONTINUE
!
!      LOAD COEFFICIENTS A, B AND C INTO RHS VECTOR RJ, WEIGHTED BY
!      PROPERTY FACTOR FOR EITHER TE OR TM.  FOR TM MODE, COEFFICIENT
!      W.R.T. RECEIVER IS ALWAYS ZERO DUE TO NO SECONDARY H-FIELD AT
!      SURFACE.
!
       IF(IDX .EQ. 3 .AND. IV .EQ. -1)A = 0.
       RJ(ILD,IV1) = COS(AN)*CMPLX(A,0.)/UC
       RJ(ILD,IV2) = COS(AN)*CMPLX(B,0.)/UC
       RJ(ILD,IV3) = COS(AN)*CMPLX(C,0.)/UC
      END IF
!
!     HORIZONTAL DERIVATIVES FOR VERTICAL FIELD COMPONENTS.  THESE
!     ARE DONE FOR THE TE MODE GENERALLY AND THE TM FOR TOPOGRAPHY.
!
      NVGO = 0
      IF(IX .EQ. 1 .AND. IDX .EQ. 3)NVGO = 1
      IF(IX .EQ. 2 .AND. IDX .EQ. 2)NVGO = 1
      IF(NVGO .EQ. 1)THEN
       IF(IH .EQ. 1)THEN                    !FIRST-ORDER DEFAULTS
        A = -1./DLY(IY1)
        B = 1./DLY(IY1)
        C = 0.
       ELSE IF(IH .EQ. -1)THEN
        A = 0.
        B = -1./DLY(IY2)
        C = 1./DLY(IY2)
       END IF
       IF(I .EQ. 1 .OR. I .EQ. NODEY)GO TO 51
       IF(IH .NE. 0)THEN
        NHS = 1
        IF(IH .EQ. 1)NHS = 0
        IF(NPT(ISB,2+IV,NCL-NHS) .NE. &
     &     NPT(ISB,2+IV,NCL+1-3*NHS))GO TO 51
       END IF
!
       CALL DRVADJ(DLY(IY1),DLY(IY2),IH,A,B,C)
!
51     CONTINUE         
!
       ANM = AN
       IF(IDX .EQ. 2)ANM = AN + 1.570796         !ADD PI/2 FOR TE
       IF(IV .EQ. -1 .AND. IDX .EQ. 3)THEN
        IF(IH .EQ. 1)A = 0.
        IF(IH .EQ. 0)B = 0.
        IF(IH .EQ. -1)C = 0.
       END IF
!
       IF(IH .EQ. 0)THEN
        RJ(ILD,IH1) = -SIN(ANM)*CMPLX(A,0.)/UC
        RJ(ILD,IH2) = RJ(ILD,IH2) - SIN(ANM)*CMPLX(B,0.)/UC
        RJ(ILD,IH3) = -SIN(ANM)*CMPLX(C,0.)/UC
       ELSE IF(IH .EQ. 1)THEN
        RJ(ILD,IH1) = RJ(ILD,IH1) - SIN(ANM)*CMPLX(A,0.)/UC
        RJ(ILD,IH2) = -SIN(ANM)*CMPLX(B,0.)/UC
        RJ(ILD,IH3) = -SIN(ANM)*CMPLX(C,0.)/UC
       ELSE IF(IH .EQ. -1)THEN
        RJ(ILD,IH1) = -SIN(ANM)*CMPLX(A,0.)/UC
        RJ(ILD,IH2) = -SIN(ANM)*CMPLX(B,0.)/UC
        RJ(ILD,IH3) = RJ(ILD,IH3) - SIN(ANM)*CMPLX(C,0.)/UC
       END IF
      END IF
      
      RETURN
END SUBROUTINE AUXADJ

!-----------------------------------------------------------------------
      SUBROUTINE DRVADJ(D1,D2,IF,A,B,C)
!-----------------------------------------------------------------------
! OCCAM MT2D 3.0 Package
! Subroutine Revision 3.0, November 2006 (notional only to track future fixes)
! Original code written by Phil Wannamaker.
!     CALCULATE PARABOLA FITTING COEFFICIENTS A,B AND C WHICH
!     WILL GO INTO RHS VECTOR OF RECIPROCITY SOLUTION FOR THE
!     AUXILIARY FIELD JACOBIANS.  PARABOLA COEFFICIENTS ARE
!     CALCULATED FROM SEPARATIONS OF THE THREE POINTS AT LEFT, 
!     MIDDLE OR RIGHT END.  FORM IS  R = A*R1 + B*R2 + C*R3.  
!
      REAL A,B,C,D1,D2,DEN

      DEN = D1*D2*D2 + D1*D1*D2
      IF(IF)10,20,30
!
!     DERIVATIVE AT LEFT POINT (IF=-1)
!
10    CONTINUE
       A = (-D2*D2 - 2.*D1*D2)/DEN
       B = (D2*D2 + D1*D1 + 2.*D1*D2)/DEN
       C = -D1*D1/DEN
      GO TO 40
!
!     DERIVATIVE AT MIDDLE POINT (IF=0)
!
20    CONTINUE
       A = -D2*D2/DEN
       B = (D2*D2 - D1*D1)/DEN
       C = D1*D1/DEN
      GO TO 40
!
!     DERIVATIVE AT RIGHT POINT (IF=1)
!
30    CONTINUE
       A = D2*D2/DEN
       B = -(D2*D2 + D1*D1 + 2.*D1*D2)/DEN
       C = (D1*D1 + 2.*D1*D2)/DEN
40    CONTINUE
      
      RETURN
END SUBROUTINE DRVADJ

!-----------------------------------------------------------------------
      SUBROUTINE APPRMTJ
!-----------------------------------------------------------------------
! OCCAM MT2D 3.0 Package
! Subroutine Revision 3.0, November 2006 (notional only to track future fixes)
! Original code written by Phil Wannamaker.
!     CALCULATE JACOBIANS OF LOG APP. RES., IMPEDANCE PHASE AND TIPPER
!     FROM E AND H-FIELD JACOBIANS.  FINAL QUANTITIES ARE WITH RESPECT 
!     TO LOG RHO.
!
      USE Fwd_mod

      COMPLEX FT,ZIMP,FTD,TPD,ZIMPD

      PI = 3.1415927
      RAD = 180./PI
      W = 2.*PI*FREK(LFRQ)
      WU = W*PI*4.E-07
      ELOG = 2.3025851
!
!     LOOP OVER RECEIVER NODES
!
      DO 10 IN = 1,NRC
       I = NRL(IN)
       NTO = NTAO(I) + NTR
       FT = RM(NTA(I)) + FPA(NTO)           ! TOTAL PARALLEL FIELD
       IF(IDX .EQ. 2)ZIMP = FT/AFT(IN)      ! TE -> E(par)/H(aux)
       IF(IDX .EQ. 3)ZIMP = AFT(IN)/FT      ! TM -> E(aux)/H(par)
!
!      LOOP OVER PARAMETERS
!
       DO 20 KP = 1,NPRM
        PFC = -1./y(LAPRM(KP)+1)
        FTD = FLDJ(IN,LAPRM(KP))            !PARALLEL FIELD JACOBIAN
        AFTD(KP) = FLDJH(IN,LAPRM(KP))      !HORIZONTAL AUX. JACOBIAN
        AFND(KP) = FLDJV(IN,LAPRM(KP))      !VERTICAL AUX. JACOBIAN
        
        ! DGM Oct 2006 - modified to return the jacobians in non-log
        !       space.  The input is not in log, the output should
        !       also be NOT in log.  Let the caller figure it out.
        !       These derivatives are wrt sigma_model (not sigma_apparent),
        !       e.g. d-rho_apparent / d-sigma_model.
        !   Also note: Phase is returned in degrees.
        IF(IDX .EQ. 2)THEN                            !TE MODE
         ZIMPD = (FTD - ZIMP*AFTD(KP))/AFT(IN)
         
         ! Zxy
!         zimped(kp,in,lfrq) = ELOG*PFC*ZIMPD
         zimped(kp,in,lfrq) = ZIMPD
         
         ! TE Apparent Resistivity
!         RTED(KP,IN,LFRQ) = PFC*(2./WU)*(REAL(ZIMP)*REAL(ZIMPD) &
!     &   + AIMAG(ZIMP)*AIMAG(ZIMPD))/RHTE(IN,LFRQ)
         RTED(KP,IN,LFRQ) = (2./WU) &
            * (REAL(ZIMP)*REAL(ZIMPD) + AIMAG(ZIMP)*AIMAG(ZIMPD))
         
         ! TE Phase
         AIOR = AIMAG(ZIMP)/REAL(ZIMP)
!         PTED(KP,IN,LFRQ) = ELOG*PFC*RAD*(1./(1. + AIOR*AIOR)) &
!     &   *(AIMAG(ZIMPD) - REAL(ZIMPD)*AIOR)/REAL(ZIMP)
         PTED(KP,IN,LFRQ) = RAD*(1./(1. + AIOR*AIOR)) &
                 *(AIMAG(ZIMPD) - REAL(ZIMPD)*AIOR)/REAL(ZIMP)
         
         ! Tipper
!         TPD = ELOG*PFC*(AFND(KP) - KZTE(IN,LFRQ)*AFTD(KP))/AFT(IN)
         TPD = (AFND(KP) - KZTE(IN,LFRQ)*AFTD(KP))/AFT(IN)
         KTED(KP,IN,LFRQ) = TPD
         
        ELSE                                          !TM MODE
         ZIMPD = (AFTD(KP) - ZIMP*FTD)/FT
         
         ! Zyx
!         zimpmd(kp,in,lfrq) = ELOG*PFC*ZIMPD
         zimpmd(kp,in,lfrq) = ZIMPD
         
         ! TM Apparent Resistivity
!         RTMD(KP,IN,LFRQ) = PFC*(2./WU)*(REAL(ZIMP)*REAL(ZIMPD) &
!     &   + AIMAG(ZIMP)*AIMAG(ZIMPD))/RHTM(IN,LFRQ)
         RTMD(KP,IN,LFRQ) = (2./WU) &
            * (REAL(ZIMP)*REAL(ZIMPD) + AIMAG(ZIMP)*AIMAG(ZIMPD))
         
         ! TM Phase
         AIOR = AIMAG(ZIMP)/REAL(ZIMP)
!         PTMD(KP,IN,LFRQ) = ELOG*PFC*RAD*(1./(1. + AIOR*AIOR)) &
!     &   *(AIMAG(ZIMPD) - REAL(ZIMPD)*AIOR)/REAL(ZIMP)
         PTMD(KP,IN,LFRQ) = RAD*(1./(1. + AIOR*AIOR)) &
                *(AIMAG(ZIMPD) - REAL(ZIMPD)*AIOR)/REAL(ZIMP)
         
         ! Tipper
!         TPD = ELOG*PFC*(AFND(KP) - TZTM(IN,LFRQ)*AFTD(KP))/AFT(IN)
         TPD = (AFND(KP) - TZTM(IN,LFRQ)*AFTD(KP))/AFT(IN)
         TTMD(KP,IN,LFRQ) = TPD
         
        END IF
20     CONTINUE
10    CONTINUE

      RETURN
END SUBROUTINE APPRMTJ
