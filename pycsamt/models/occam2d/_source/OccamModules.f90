! Occam_Mod.f90 - module of shared memory objects.
!   Accumulated from a variety of old common blocks and made
!   so that the memory can be allocated at runtime instead of
!   at compile time.
!
! David Myer
! Oct 2006
!----------------------------------------------------------------------

module Occam_mod 
    ! constraint data for the inversion:
    logical                             :: gbAddDiags       ! DGM Oct 2006
    integer                             :: nd
    real, dimension(:), allocatable     :: d, sd
    real, dimension(:,:), allocatable   :: dp

    ! gbAddDiags = true if diagonal penalties should be added to horiz & vert.
    ! nd = number of data.  + 1 or 2 if have statics
    ! d() = data 
    ! sd() = data errors (one standard error) - usually read from data file
    ! dp() = array of parameters associated with each datum (freq., etc)
    !
    !-----------------------------------------------------------------------
    ! model parameters for the current iteration and its response
    integer                         :: nParams
    real, dimension(:), allocatable :: pm, dm
    !   nParams = number of free model parameters
    !   pm() = array of model parameters
    !   dm() = array of model response (data domain)

    !-----------------------------------------------------------------------
    ! stuff to be carried into critical subroutines 
    real, dimension(:,:), allocatable   :: ptp, wjtwj
    real, dimension(:), allocatable     :: wjtwd, pwk1, pwk2, pwk3, prewts, premod, dwk1, dwk3
    real                                :: frac
    integer                             :: iout, idebug, nfor, nFracTimes
    !   ptp() = product of penalty matrix with its transpose
    !   wjtwj() = product of weighted jacobian matrix with its transpose
    !   wjtwd() = product of weighted jacobian matrix with weighted translated data
    !   pwk1,2,3() = vectors to carry around working copies of models
    !   frac = current reduction in step size
    !   nfor = tally for total number of forward calculations
    !   iout = unit number of output (usually LOGFIL)
    !   idebug = current debug level
    !   prewts() = weights for model prejudice
    !   premod() = model prejudice
    !   dwk1,3 = working copies of the model responses - 'dm' is used as dwk2
    !   nFracTimes = # of times to cut stepsize before giving up (DGM Nov 2006 - was hardcoded)

    !-----------------------------------------------------------------------
    ! Working variables used by TOFMU - a function that is called many times
    !   per occam iteration and which forms 90% of the program run time.
    !   These vars were allocated automatically at every call, but it is more
    !   efficient to do it one time and let the memory stand.
    real, dimension(:,:), allocatable   :: aMat
    real, dimension(:), allocatable     :: pwk4, dwk4
end module Occam_mod


!----------------------------------------------------------------------

module Occam_ModelLimits

! DGM Nov 2006 - Values read from the startup iteration file controlling
!   the span of values the model can take as well as the number of steps
!   per decade.  Values are read in ItrIn() and used in TOFMU().

    logical ::  gbModelLimited
    real    ::  gnModelMin, gnModelMax
    logical ::  gbModelStepped
    real    ::  gnModelSteps  
contains

    subroutine ModelLimits_Init
        implicit none
        gbModelLimited  = .false.
        gbModelStepped  = .false.
        gnModelMin      = 0
        gnModelMax      = 0
        gnModelSteps    = 0
    end subroutine ModelLimits_Init
end module Occam_ModelLimits

