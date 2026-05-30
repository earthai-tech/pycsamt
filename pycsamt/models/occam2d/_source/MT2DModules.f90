
module Fwd_mod
!
!  This is used by Occam for declaring the number of data parameters, not
!  including the data value or standard error
!
	 integer, parameter :: cnDataParams  = 4 ! 4 is dummy value for now
    !   cnDataParams    = number of parameters for each data point 
    !                   (e.g. for MT: site, freq, data type, dummy)
    !  dp() = array of parameters associated with each datum (freq., etc)
    !  dp(:,1) = site number
    !  dp(:,2) = frequency number
    !  dp(:,3) = data type (see below)
    !  dp(:,4) = unused 
    !
    ! Data Types:
    !   1 = TE resistivity (log10)
    !   2 = TE phase
    !   3 = real(tipper)
    !   4 = imag(tipper)
    !   5 = TM resistivity (log10)
    !   6 = TM phase
    !   7 = reserved for TM Ez/Ey tipper
    !   8 = reserved for TM Ez/Ey tipper
    !   9 = TE resistivity (linear)                 DGM Oct 2006
    !  10 = TM resistivity (linear)                 DGM Oct 2006
    !
    !  11 = real(Zxx)  ! not used yet until we incorporate site rotation angle        
    !  12 = imag(Zxx)  ! not used yet until we incorporate site rotation angle  
    !  13 = real(Zxy)          
    !  14 = imag(Zxy)    
    !  15 = real(Zyx)        
    !  16 = imag(Zyx)  
    !  17 = real(Zyy)  ! not used yet until we incorporate site rotation angle        
    !  18 = imag(Zyy)  ! not used yet until we incorporate site rotation angle
    
   
    
    
    !-- internal use -- output only:
    !  -1 = static shift: sum TE & TM together
    !  -2 = static shift: TE sum
    !  -3 = static shift: TM sum
    
    !-----------------------------------------------------------------------
    ! information about the regularization bricks in terms
    ! of the FE mesh, as well as info about the penalty matrix
    integer, allocatable            :: icorn(:,:), nrcol(:), ijexep(:,:)
    real, dimension(:), allocatable :: thick, width, y, prewt, expen
    integer :: nlay, nrfix, nbrick, nexep
    !   nlay = number of layers in regularization grid
    !   nrcol() = number of regularization bricks in each row
    !   thick() = vertical thickness of each reg. brick 
    !   width() = width of each reg. brick in metres
    !   nrfix = number of fixed resistivities in mesh (an offset for the y() vector)
    !   y() = lookup table of resistivities, passed into pw2di
    !   icorn() = top, left, bottom, right coords of every free brick in terms of FE nodes
    !   nbrick = number of free resistivity bricks (my copy of nprm)
    !   nexep = number of exceptions in penalty matrix
    !   ijexep() = brick parameter indices of exception pairs
    !   expen() = penalties for exceptions
    !   prewt() = penalties for model prejudice

    !-----------------------------------------------------------------------
    ! station information and other stuff pertinent to the data array
    character(50), dimension(:), allocatable :: cSiteName
    real, dimension(:), allocatable     :: sitloc, frek
    integer, dimension(:), allocatable  :: nrl
    integer :: nrc, nfre, ndata
    real    :: boffs
    !   nrc = number of receiver sites
    !   sitloc() = site locations by distance in meters from origin
    !   nrl() = site locations by FE mesh node number
    !   boffs = binding offset (location of right edge of the left terminating layers of
    !           regularization mesh in the same coords as sitloc())
    !   nfre = number of frequencies
    !   frek() = frequencies (Hz)
    !   ndata = number of actual MT data (c.f. static shift constraints)

    !-----------------------------------------------------------------------
    ! info on the static shifts, both fixed and free.
    !   Size of the arrays should not be any larger than (2 modes) times
    !   maximum number of sites
    integer, dimension(:), allocatable  :: isx1, isx3, isfx
    real, dimension(:), allocatable     :: shift
    integer :: iscon, nfss, nsb
    !   nsb = number of blocks of static shift data (one per site/mode)
    !   iscon = type of static shift data constraint:
    !     0 = no constraint
    !     1 = TE and TM modes sum to sum of log(initial values) together
    !     2 = TE and TM modes sum to sum of log(initial values) separately
    !   nfss = number of free (invertible) static shifts
    !   isx1() = index of data stations for each shift block
    !   isx3() = intex of data type for each shift block (should be only
    !     1 = TE resisitivity or 5 = TM resistivity)
    !   isfx() = lookup table for invertible shifts. 0 = fixed shift,
    !     n = invertible shift is nth parameter (PM(n))
    !   shift() = log10(static shift corrections) to be applied to the data


    COMPLEX, dimension(:,:), allocatable    :: FLDJ, FLDJH, FLDJV

    !------------------------------------------------
    ! From /BLKC/ embedded in the code
    COMPLEX, allocatable, dimension(:,:,:)  :: AA55
    COMPLEX :: AAD55(4)
      
    !------------------------------------------------
    ! From BLKD.inc
    REAL, allocatable, dimension(:) :: YY, ZZ, DLY, DLZ
    !     NOTE: DLY AND DLZ WERE DL AND DZ IN AUXFLD AND BNORM
    !     YY     - OFFSET OF EACH COLUMN IN METERS
    !     ZZ     - DEPTH OF EACH ROW IN METERS
    !     DLY    - WIDTH OF EACH COLUMN IN METERS
    !     DLZ    - THICKNESS OF EACH ROW IN METERS
    
    !------------------------------------------------
    ! From BLKE.inc
    integer, allocatable, dimension(:,:,:)  :: npt
    real, allocatable, dimension(:,:,:)     :: sig
    complex, allocatable, dimension(:)      :: sigh

    !------------------------------------------------
    ! From BLKF.inc
    COMPLEX, allocatable, dimension(:)  :: FPA, FPB, AFN, AFT, CSIGU
    
    !------------------------------------------------
    ! From BLKI.inc
    INTEGER, allocatable, dimension(:)  :: NTAO, NTA, LAPRM, NPARF
    INTEGER :: IDX,NDX,IMD,NMODE,NODEY,NODEZ,NTR,NTQ,NDZ,NDZ1
    INTEGER :: NZZB(2),NOLA(2),LSID,NPRM,LFRQ,NNODE,ISEA

    !------------------------------------------------
    ! From BLKJ.inc
    ! These are outputs from PW2DI().  Derivative outputs are only
    !   valid when the option to create the Jacobian is specified.
    COMPLEX, allocatable, dimension(:,:)    :: KZTE, TZTM
        ! KZTE = tipper
    COMPLEX, allocatable, dimension(:,:,:)  :: KTED, TTMD
        ! Jacobians of the above
    COMPLEX, allocatable, dimension(:)      :: AFND, AFTD
    REAL, allocatable, dimension(:,:)       :: RHTE, PHTE, RHTM, PHTM
        ! RHo-TE, PHase-TE, etc...
    REAL, allocatable, dimension(:,:,:)     :: RTED, PTED, RTMD, PTMD
        ! Rho-TE-Derivative, Phase-TE-Deriv, etc...
    
    !------------------------------------------------
    ! From BLKJC.inc
    ! These are outputs from PW2DI().  Derivative outputs are only
    !   valid when the option to create the Jacobian is specified.
    COMPLEX, allocatable, dimension(:,:)    :: zimpe, zimpm
        ! XY and YZ transfer functions (XX & YY assumed to be zero).
    COMPLEX, allocatable, dimension(:,:,:)  :: zimped, zimpmd
        ! Jacobians for the above.

    !------------------------------------------------
    ! From BLKL.inc
    REAL, allocatable, dimension(:,:)   :: DLR, RLR
    REAL, allocatable, dimension(:)     :: DDP, RL

    !------------------------------------------------
    ! From /BLKP/ embedded in the code
    COMPLEX, allocatable, dimension(:)  :: KL1, RFC, COA, ARG, AEX
    REAL, allocatable, dimension(:)     :: HL
    COMPLEX                             :: KL10, ZL1, ZLA

    !------------------------------------------------
    ! From /FORSOL/ embedded in the code
    COMPLEX, allocatable :: S(:,:), RM(:), RJ(:,:)
    
    !
    ! Parameters for PW2D:
    !       
        integer, parameter :: cnAirLayers   = 10
        integer, parameter :: cnMaxFixedRes = 36
    !   cnAirLayers     = # of Z layers of air added for TE or TE+TM inversion
    !   cnMaxFixedRes   = max # of fixed resistivities that can be provided
    !       by the user in the mesh file.  They are referred to in the
    !       weird mesh thingee by numbers 0-9 and letters A-Z.  Note that
    !       Z is automatically assigned to seawater.

    
!---------------------------------------------------------------------------
! This is used only for 2D forward modeling calculations (not during inversions)
! contains TE/TM strike fields for entire model
        integer UNITTE, UNITTM
        logical loutfields
        common /outfields/ loutfields,UNITTE, UNITTM
    
    
contains
    ! Routine to allocate the stuff required by the Fwd calcs the first
    !   time through.  Once the vars are allocated, the call does nothing.
    subroutine Fwd_Allocate( nRcvrs, nResists, nParams, nFreqs &
                           , nModLayers, nNodeY, nNodeZ )
        implicit none
        ! Parameters
        integer, intent(in)     :: nRcvrs, nResists, nParams, nFreqs
        integer, intent(in)     :: nModLayers, nNodeY, nNodeZ
        ! Locals
        integer :: nAllocErr
        
        ! Have we already been in here?
        if (allocated(FLDJ)) return
        
        ! Allocate.  Stop suddenly if can't.
        allocate( FLDJ( nRcvrs, nResists )      &
                , FLDJH( nRcvrs, nResists )     &
                , FLDJV( nRcvrs, nResists )     &
                , KZTE(nRcvrs,nFreqs)           &
                , TZTM(nRcvrs,nFreqs)           &
                , AFND(nParams)                 &
                , AFTD(nParams)                 &
                , KTED(nParams,nRcvrs,nFreqs)   &
                , TTMD(nParams,nRcvrs,nFreqs)   &
                , RHTE(nRcvrs,nFreqs)           &
                , PHTE(nRcvrs,nFreqs)           & 
                , RHTM(nRcvrs,nFreqs)           &
                , PHTM(nRcvrs,nFreqs)           &
                , RTED(nParams,nRcvrs,nFreqs)   &
                , PTED(nParams,nRcvrs,nFreqs)   &
                , RTMD(nParams,nRcvrs,nFreqs)   &
                , PTMD(nParams,nRcvrs,nFreqs)   &
                , DLR(nModLayers-1,2)           &
                , stat=nAllocErr)
 !  Maximum number of continuation lines is 19 in Fortran 95,
 ! Fortran 90, FORTRAN 77, and SAA fixed form and 39 in Fortran 95 and Fortran 90 free form.
        if (nAllocErr == 0) then
            allocate( RLR(nModLayers,2)         &
                , DDP(nModLayers-1)             &
                , RL(nModLayers)                &
                , ARG(nModLayers-1)             &
                , AEX(nModLayers-1)             &
                , HL(nModLayers-1)              &
                , KL1(nModLayers)               &
                , RFC(nModLayers)               &
                , COA(nModLayers)               &
                , sig(nNodeZ,4,nNodeY)          &
                , sigh(nNodeZ)                  &
                , FPA(nNodeZ), FPB(nNodeZ)      &
                , AFN(nNodeY), AFT(nNodeY)      &
                , CSIGU(nNodeY), NPARF(nNodeY)  &
                , NTAO(nNodeY), NTA(nNodeY)     &
                , zimpe(nRcvrs,nFreqs)          &
                , zimpm(nRcvrs,nFreqs)          &
                , stat=nAllocErr)
        endif
        if (nAllocErr == 0) then
            allocate( zimped(nParams,nRcvrs,nFreqs) &
                , zimpmd(nParams,nRcvrs,nFreqs) &
                , S(nNodeZ+2,nNodeY*nNodeZ)     &
                , RM(nNodeY*nNodeZ)             &
                , RJ(nNodeY,nNodeY*nNodeZ)      &
                , AA55(4,nNodeZ,nNodeY)         &
                , stat=nAllocErr)
        endif
        
        ! If we didn't survive allocation, explain nicely & stop the program
        if (nAllocErr .ne. 0) then
            write(*,*) 'Out of memory.  Try reducing the number of model blocks.'
            stop 'Out of memory.  Try reducing the number of model blocks.'
        endif
        
        return
    end subroutine Fwd_Allocate
    
    subroutine Fwd_Dealloc
        implicit none
        integer n
        
        deallocate( FLDJ, FLDJH, FLDJV, KZTE, TZTM, AFND, AFTD, KTED &
                  , TTMD, RHTE, PHTE, RHTM, PHTM, RTED, PTED, RTMD, PTMD &
                  , DLR, RLR, DDP, RL, ARG, AEX, HL, KL1, RFC, COA, sig &
                  , sigh, FPA, FPB, AFN, AFT, CSIGU, NPARF, NTAO, NTA &
                  , zimpe, zimpm, zimped, zimpmd, S, RM, RJ, AA55 &
                  , stat=n )  ! stat=  doesn't allow process to blow chunks on dealloc error
        ! Module vars not allocated in Fwd_Allocate but belonging to the module
        !   so deallocated here. (They are only allocated elsewhere because they 
        !   are needed early - usually in an OccamIO.f90 routine.)
        deallocate( LAPRM, YY, DLY, ZZ, DLZ, NPT &
                  , stat=n )  ! stat=  doesn't allow process to blow chunks on dealloc error
 
        deallocate( thick, width, y, icorn, prewt,nrl  &
                  , sitloc, nrcol,isfx, shift, frek, isx1, isx3 &
                  , stat=n )
        
        ! These things are conditionally allocated
        if (nexep > 0) deallocate( ijexep, expen, stat=n )
        
    end subroutine Fwd_Dealloc
    
end module Fwd_mod
