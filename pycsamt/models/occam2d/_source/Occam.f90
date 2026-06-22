! OCCAM 3.0 Package
! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
!
! REVISION HISTORY:
!     Oct/Nov 2006 David Myer with help from Kerry Key    (Version 3.0)
! Big re-write! Moved to F90, allocated memory, cmdline params,
! modules, and general code reorg.  
!     July 2001 deltdel for fast roughning matrix multiply Catherine deGroot-Hedlin 
!     incredibly stable for 9 years!
!     APRIL 1992 (MORE THANKS TO CATHERINE AND ALSO GEOTOOLS)  (VERSION 2.0)
! Version 2.0 is a re-write to make OCCAM independent of the dimensionality
! of the forward problem and make maximum use of dynamic memory allocation.
! The calling program is now an integral part of the package since all the
! model-dependent stuff has been removed.
!     SEPTEMBER 1989,                                     (VERSION 1.5.2)
!     JANUARY 1989 (WITH THANKS TO C.DEGROOT-HEDLIN),     (VERSION 1.5)
!     AUGUST 1988,                                        (VERSION 1.4)
!     OCTOBER 1987,                                       (VERSION 1.3)
!     MARCH 1986,                                         (VERSION 1.2)
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
program runocc
!-----------------------------------------------------------------------
!
! OCCAM 3.0 Package
! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
! Program Revision 3.0, November 2006
! Program Revision 2.01, 13 Jan 1993
!
! runocc is the calling program for OCCAM v3.0.
!

    ! includes:
    USE Occam_Mod
    use Fwd_mod, only: Fwd_Dealloc
    IMPLICIT NONE
      
    ! Functions
    real ANORM

! local variables:
    integer lerr, maxitr, iruf, nit, itmax, konv, ifftol, n, i, iargc
!   lerr = error flag for file operations
!   maxitr = maximum number of subsequent iterations by this run
!   iruf = roughness measure -> see OBJMAT for details
!   nit = current iteration number (0=start)
!   i = loop index
!   itmax = maximum iteration number
!   konv = convergence flag from OCCAM
!   ifftol = convergence info from OCCAM (has the desired misfit been achieved?)
    real tolreq, tobt, pmu, stepsz, rlast
!   tolreq = desired misfit, in RMS standard errors
!   tobt = current misfit
!   pmu = current Lagrange multiplier
!   stepsz = current RMS step though model space
!   rlast = last RMS measure of roughness
    character(50) itform, descr, modfil, datfil, cStartup, cRoot
    character(2) :: cTemp
!   itform = format number of iteration files
!   descr = description line in iteration files
!   modfil = model file name
!   datfil = data file name
!   cStartup = name of starting file.  (See notes below)
!   cRoot = root filename for output files.  (See notes below)
    logical   :: bFwdOnly
!   bFwdOnly = set true if -F on cmd line ... saying to produce a fwd model & exit.
    logical :: lrunocc
!   lrunocc = .true. to run the occam iterations.  Some compilers do not
!     like looping over a static .true.
! parameters
    integer iof1, iof2
    parameter (iof1 = 13, iof2 = 15)
!   iof1 = unit number for the LOGFIL (passed to /result/iout)
!   iof2 = unit number for reading and writing ITER files
!-----------------------------------------------------------------------

! DGM Oct 2006 - look for a command line parameter telling the name of
! the startup file.  If not given, then use 'startup'.
! If two params are given, then first is the startup file and the 2nd
! is the root filename for output of iteration, response, & log files.
    write(*,*) '------------------------------------------------------------'
    write(*,*) 'OCCAM Inversion v3.0                                Dec 2006'
    write(*,*) '------------------------------------------------------------'

    cStartup  = 'startup'     ! Establish defaults
    cRoot     = ''
    bFwdOnly  = .false.
    nFracTimes= 8           ! DGM Nov 2006 old way was equivalent to 17.  Found in
                            ! testing that if it happened >= 5 times, it went up to
                            ! 17 then aborted or there was no change in the RMS
                            ! misfit.  Basically, when the inversion gets to the
                            ! point where it needs lots of cuts in step size then
                            ! it has found its convergence point.

    n = iargc()
    if (n >= 1) then
        i = 1
        call getarg( i, cStartup )
        if (cStartup(1:1) == '?') then     ! help requested
            n = 9999
        elseif ( (cStartup(1:1) == '-' .or. &
                  cStartup(1:1) == '/' .or. & 
                  cStartup(1:1) == char(92))  & !KWK:char(92) is backslash, safer than using '\' since some compilers treat this as an escape character
           .and. (cStartup(2:2) == 'f' .or.   &
                  cStartup(2:2) == 'F') )  then
            bFwdOnly = .true.
            i = 2
            if (n >= i) then
                call getarg( i, cStartup )
            else
                cStartup = 'startup'
            endif
        endif
        i = i + 1
        
        if (n >= i) then
            call getarg( i, cRoot )
            
        elseif (n > i) then
            write(*,*) ' '
            write(*,*) 'Occam2D accepts one flag:'
            write(*,*) '      -F    When given, this flag means produce a forward model only.'
            write(*,*) '            The forward model will have the name <output root name>.fwd.'
            write(*,*) '            If no output root is given, the name is Forward.fwd'
            write(*,*) ' '
            write(*,*) 'Occam2D has 2 optional parameters:'
            write(*,*) '      <startup file> -- name of the file with the starting model'
            write(*,*) '                    and runtime parameters.  Traditionally this'
            write(*,*) '                    file has been called ''startup''.  To continue'
            write(*,*) '                    a previous inversion, specify an iteration file.'
            write(*,*) '                    If not given, ''startup'' is assumed.'
            write(*,*) '      <output root name> -- This parameter can only be given if '
            write(*,*) '                    a startup filename is given.  It is the root'
            write(*,*) '                    name to use for the output iteration, response,'
            write(*,*) '                    and log files.  If not given, these files will'
            write(*,*) '                    be named ITERxx.iter, RESPxx.resp, LogFile.logfile.'
            write(*,*) '                    If given, the root name will be used in place of'
            write(*,*) '                    ''ITER'', ''RESP'', and ''LogFile'' in the above names.'
            write(*,*) '                    Do NOT include a path or a ''.''.'
            stop
        end if
    end if


    ! Open startup file and read stuff.
    ! (Startup and iteration files should be machine written, and so we do
    ! not necessarily test the validity of the input data unless arrays are 
    ! involved.)
    write(*,*) 'Reading startup model from ', trim(cStartup), '...'
    open (iof2, file=trim(cStartup), status='old', iostat=lerr)
    if (lerr .ne. 0) then
        write(*,*) ' Error opening startup file'
        stop 
    end if
    call itrin (itform,descr,modfil,datfil,maxitr,tolreq,iruf,    &
          nit,pmu,rlast,tobt,ifftol,iof2)
    close(iof2)

    ! input the data
    write(*,*) ' '
    write(*,*) 'Reading data from ', trim(datfil), '...'
    call inputd(datfil)

    ! input the fixed model characteristics
    write(*,*) ' '
    write(*,*) 'Reading model from ', trim(modfil), '...'
    call inputm(modfil)
    
    ! If we are to produce a forward model only, do so & exit.
    if (bFwdOnly) then
        write(*,*) 'Outputting FORWARD MODEL ONLY...'
        
        call FwdCalc( .false., nParams, nd, pm, dp, dm, d )
        
        ! Show the misfit
        DWK4 = (D - DM)/SD        ! Array math - stack alloc of only nd size
        rlast = SQRT( ANORM( ND, DWK4 ) / ND )
        write(*,*) '  Misfit:', rlast
        
        ! Output the forward response
        if (cRoot == '') then
            call OutFwd( 'Forward.fwd' )
        else
            call OutFwd( trim(cRoot)//'.fwd' )
        endif
        
    else        ! Invert!
    
        ! if this is zeroth iteration copy the startup file into ITER00
        if (nit .eq. 0) then
            ! initialize the inversion parameters
            pmu = 5.0
            rlast = 1.0e+10
            tobt = 1000.0
            ifftol = 0
            call itrout (itform,descr,modfil,datfil,maxitr,tolreq, &
              iruf,nit,pmu,rlast,tobt,ifftol,iof2,cRoot)
        end if

        itmax = maxitr + nit

        iout = iof1
        if (cRoot == '') then
            open (iout, file='LogFile.logfile', iostat=lerr)
        else
            open (iout, file=trim(cRoot)//'.logfile', iostat=lerr)
        end if
        if (lerr .ne. 0) then
            write(*,*) ' Error opening log file'
            stop 
        end if
        write(*,*) 'Major I/O accomplished:'
        write(*,*) '  sit back and relax.......'

        lrunocc= .true.
        
        ! Run the main loop
        do while (lrunocc)
            write (cTemp,'(I2)') nit+1
            call occam(iruf,tolreq,tobt,ifftol,nit,pmu,rlast,stepsz,konv)

            ! write an iteration file
            call itrout (itform,descr,modfil,datfil,maxitr,tolreq,iruf, &
         &      nit,pmu,rlast,tobt,ifftol,iof2,cRoot)

            ! let the console know roughly what is going on
            write(*,*) ' ------------------------------'
            write(*,*) ' ITERATION  =', nit
            write(*,*) ' RMS MISFIT = ', tobt
            if (ifftol == 1) then
                write(*,*) '     Target = ', tolreq, ' reached!'
            else
                write(*,*) '     Target = ', tolreq
            endif
            write(*,*) ' ROUGHNESS  = ', rlast
            write(*,*) ' STEPSIZE   = ', stepsz
            write(*,*) ' ------------------------------'
            write(*,*) ' '

            ! test for convergence: continue inversion if 
            !   1) iteration number does not exceed maximum, and
            !   2) stepsize is still significant, and
            !   3) rms misfit larger than required, and 
            !   4) there are no irrecoverable convergence problems
            if (konv .ne. 0) then
                select case (konv)
                case (1)
                    write(iof1,*)   'Perfectly smooth model found.'
                    write(*,*)      'Perfectly smooth model found.'
                case (2)
                    write(iof1,*)   'Convergence problems from RMS Misfit'
                    write(*,*)      'Convergence problems from RMS Misfit'
                    write(*,*)      'May have reached minimum possible RMS Misfit level!'
                case (3)
                    write(iof1,*)   'Convergence problems from Smoothness'
                    write(*,*)      'Convergence problems from Smoothness'
                    write(*,*)      'May have reached minimum possible RMS Misfit level!'
                case default
                    write(iof1,*)   ' Convergence problems of type ', konv
                    write(*,*)      ' Convergence problems of type ', konv
                end select
                close(iof1)
                exit
                
            else if (nit .ge. itmax) then
            ! Inversion has reached max # of iterations specified in startup file.
            !   Stop now even if not converged.
            ! NOTE: This check MUST be before the "continue" check below.
                write(iof1,*) ' Stop on maximum iterations'
                close(iof1)
                write(*,*) ' Reached maximum iterations.  Done.'
                exit
                
            else if ((tobt .lt. 1.01*tolreq) .and. (stepsz .lt. 0.001)) then
                write(iof1,*) ' Stop on normal convergence'
                close(iof1)
                write(*,*) ' Solution converged!  Done.'
                exit
            !    
            ! else... continue...    
            
            end if
        enddo
    endif
    
    ! Deallocate memory & return control to the user
    call Fwd_Dealloc
    deallocate( pm, ptp, wjtwj, wjtwd, pwk1, pwk2, pwk3, prewts, premod &
              , aMat, pwk4, d, sd, dp, dm, dwk1, dwk3, dwk4 &
              , stat=n )        ! stat=  doesn't allow process to blow chunks on dealloc error
    
end program runocc


!-----------------------------------------------------------------------
SUBROUTINE OCCAM(IRUF,TOLREQ,TOBT,IFFTOL,NIT,PMU,RLAST,STEPSZ,KONV)
!-----------------------------------------------------------------------
! OCCAM 3.0 Package
! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
! Subroutine Revision 3.0, November 2006
! Subroutine Revision 2.01, 20 Jan 1993

! OCCAM EXECUTES ONE ITERATION OF A SMOOTH MODEL FINDER

!  ARGUMENTS:
!     IRUF = see OBJMAT function
!     TOLREQ = REQUIRED RMS MISFIT
!     TOBT = RMS MISFIT FROM LAST ITERATION (LARGE INITIALLY)
!     IFFTOL = RECORD OF WHETHER FEASIBLE MODEL HAS BEEN ATTAINED (0=NO,
!         1=YES).  INITIALLY SET 0 AND THEN LEFT ALONE.
!     NIT = THE NUMBER OF PREVIOUS ITERATIONS IN THIS INVERSION. SET TO 0 ON FIRST 
!         CALL, IT WILL BE UPDATED BY OCCAM. 
!     PMU = LOG10(LAGRANGE) MULTIPLIER TO START SEARCHES WITH.  INITIALLY 
!         SET LARGE AND THEN LEFT FROM LAST CALL FOR EFFICIENCY
!     RLAST = ROUGHNESS OF LAST MODEL.  INITIALLY SET LARGE AND THEN LEFT.        

! ON OUTPUT:

!  Occam_mod CONTAINS
!     PM(nParams) = VECTOR OF UPDATED MODEL PARAMETERS
!     DM(nData+) = VECTOR CONTAINING *APPROXIMATE* RESPONSE OF NEW MODEL
!  ARGUMENTS:
!     TOBT = RMS MISFIT OF NEW MODEL RESPONSE.
!     NIT = NUMBER OF PREVIOUS CALLS TO OCCAM1 DURING THIS INVERSION
!     STEPSZ = RMS MEASURE OF THE CHANGES IN THE MODEL PARAMETERS 
!     KONV = STATUS FLAG:
!        0 = NORMAL EXIT FOR A SUCCESSFUL ITERATION,
!        1 = A PERFECTLY SMOOTH MODEL HAS BEEN FOUND FOR THE REQUIRED TOLERANCE
!        2 = CONVERGENCE PROBLEMS.  UNABLE TO FIND A MODEL WITH AN R.M.S. MISFIT
!          LESS THAN OR EQUAL TO THAT OF THE LAST ITERATION'S
!        3 = CONVERGENCE PROBLEMS.  UNABLE TO FIND A MODEL WITH AN R.M.S. MISFIT
!          OF TOLREQ AND SMOOTHER THAN THE LAST ITERATION'S

! SUBROUTINES REQUIRED AND SUPPLIED:

!   MAKPTP, A SUBROUTINE WHICH CREATES THE ROUGHENING MATRIX TIMES ITS TRANS.
!   TOFMU(AMU), A FUNCTION WHICH RETURNS THE RMS MISFIT OF THE RESPONSE
!      OF THE MODEL CONSTRUCTED USING THE LAGRANGE MULTIPLIER AMU.
!      (CALLS CHOLIN,CHOLSL,FwdCalc,ANORM)
!   TRMULT, MULT, ATAMULT AND ANORM, SUBROUTINES FOR MATRIX OPERATIONS
!   fminocc, A SUBROUTINE WHICH MINIMISES A UNIVARIATE FUNCTION
!   FROOT, A SUBROUTINE WHICH FINDS THE ROOT OF A UNIVARIATE FUNCTION
!   MINBRK, A SUBROUTINE WHICH BRACKETS A UNIVARIATE MINIMUM
!   SCANMU, A SUBROUTINE WHICH GENERATES THE MISFIT TABLE WHEN IBUG .GT. 0

! SUBROUTINES WHICH MUST BE SUPPLIED BY THE USER:

!  FwdCalc: COMPUTES THE FORWARD FUNCTION FOR MODEL PM() AT
!    THE DATA PARAMETERS DP() AND RETURNS IT IN DM().  Optionally
!    RETURNS THE MATRIX OF PARTIALS, WJ(ND,nParams).
!  OBJMAT(), ASSEMBLES THE PENALTY MATRIX (THIS WILL BE DEPENDENT ON MODEL 
!    TYPE AND DIMENSION)


! INCLUDES:
      USE Occam_Mod
      IMPLICIT NONE

! ARGUMENTS:
      INTEGER IRUF,IFFTOL,NIT,KONV
      REAL TOLREQ,TOBT,PMU,RLAST,STEPSZ
! LOCAL VARIABLES:
      REAL TOL0, TMIN, TINT
      REAL AMU1, AMU2, AMU3, TAMU1, TAMU2, TAMU3
      REAL RUF, PMU2
! FUNCTIONS
      REAL TOFMU, fminocc, FROOT, ANORM, FNDRUF
      EXTERNAL TOFMU

! LOCAL PARAMETERS
      REAL TOLM, TOLI
! TOLM AND TOLI ARE THE TOLERANCES FOR fminocc AND FROOT RESPECTIVELY
! DECREASING THEM WILL IMPROVE ACCURACY AT THE COST OF MORE FORWARD 
! CALCULATIONS
      PARAMETER (TOLM= 0.1, TOLI = 0.001)

!-----------------------------------------------------------------------
      IF (abs(idebug) .EQ. 1) WRITE(*,*) ' Entering Occam...'
      
! SET THINGS UP FOR THIS ITERATION:
      KONV = 0
! FRAC CONTROLS THE STEP SIZE; NORMALLY WILL REMAIN AT 1.0 UNLESS WE HAVE
! CONVERGENCE PROBLEMS
      FRAC = 1.0
! CONSTRUCT WJTWJ, WJTWD, AND FIND INITIAL MISFIT
      IF (abs(idebug) .EQ. 1) WRITE(*,*) &
     &  ' Constructing Derivative Dependent Matrices...'
      CALL MAKJTJ(TOL0)

      WRITE(IOUT,*) 'STARTING R.M.S. = ',TOL0
      IF (NIT .EQ. 0) THEN
          if ((TOL0 .le. TOLREQ) .and. (IFFTOL .eq. 0)) then
            IF (abs(idebug) .ge. 1) &
                write(*,*) ' Tolerance met.  This iteration begins smoothing.'
            IFFTOL = 1
          end if
      END IF
      IF (abs(idebug) .EQ. 1) then
        WRITE(*,*) ' Starting Misfit = ', TOL0
        WRITE(*,*) ' Starting Mu     = ', 10**PMU
      endif
      WRITE(IOUT,*) '** ITERATION ',NIT+1,' **'

! CONSTRUCT THE PENALTY MATRIX MULTIPLIED INTO ITS TRANSPOSE
      IF (abs(idebug) .EQ. 1) WRITE(*,*)' Constructing Penalty Matrix...'
      CALL MAKPTP(IRUF, pm)


! THE NEXT BLOCK OF CODE CONTROLS THE SELECTION OF THE LAGRANGE MULTIPLIER,
! PMU:
      IF (abs(idebug) .EQ. 1) &
     & WRITE(*,*) ' Searching Lagrange Multiplier:'
     
120   NFOR = 0
! PRODUCE THE MISFIT FUNCTION IF REQUIRED (THIS WILL BE USED TO BRACKET MIN)
      IF (abs(idebug) .GE. 1) WRITE(*,*) ' ...Bracketing Minimum...'
      IF (abs(idebug) .GE. 2) THEN
        CALL SCANMU(AMU1,AMU2,AMU3,TAMU2,TOLREQ)
      ELSE
! BRACKET THE MINIMUM USING MINBRK AND TWO GUESSES
        AMU1 = PMU - 1.0
        AMU2 = PMU
        
        ! In: MINBRK( Mu guess #1, #2, ..., functional to calc rms)
        ! Out: AMU1,3 - brackets around AMU2 - the minimum MU
        !      TAMUx = TOFMU( AMUx ) -- the value of the rms misfit
        CALL MINBRK(AMU1,AMU2,AMU3,TAMU1,TAMU2,TAMU3,TOFMU)
      END IF
      
! FIND THE MINIMUM
      IF (TAMU2 .LT. TOLREQ) THEN
! WE'VE BEEN LUCKY AND FOUND AN ACCEPTABLE MINIMUM USING MINBRK
        TMIN = TAMU2
        PMU = AMU2
        WRITE(IOUT,*) 'MINIMUM TOL FROM MINBRK IS AT MU =',PMU
      ELSE
        IF (abs(idebug) .GE. 1) WRITE(*,*) ' ...Finding Minimum...'
        
        TMIN = fminocc(AMU1,AMU2,AMU3,TAMU2,TOFMU,TOLM,PMU)
        
        WRITE(IOUT,*) 'MINIMUM TOL FROM fminocc IS AT MU =',PMU
      END IF
      WRITE(IOUT,*) 'AND IS =',TMIN
      WRITE(IOUT,*) 'USING ',NFOR,' EVALUATIONS OF FUNCTION'
! IF THE NEW MINIMUM TOLERANCE IS GREATER THAN THE TOLERANCE FROM THE
! PREVIOUS MODEL, WE ARE HAVING CONVERGENCE PROBLEMS. CUT THE STEP SIZE
! TOL0 is tol of the previous model.
! TOLREQ is the target tol.
      IF( (  ((TMIN .GT. TOL0)   .AND. (IFFTOL .EQ. 0)) &
        .OR. ((TMIN .GT. TOLREQ) .AND. (IFFTOL .EQ. 1)) &
          ) .and. idebug .ge. 0 ) THEN
        WRITE(IOUT,*) 'DIVERGENCE PROBLEMS, CUTTING STEP SIZE'
        IF (abs(idebug) .GE. 1) &
     &    WRITE(*,*) ' ..DIVERGENCE PROBLEMS, CUTTING STEPSIZE..'
        FRAC = FRAC*2.0
        
        ! WE HAVE CUT STEP SIZE A LOT TO NO EFFECT: GIVE UP
        IF (FRAC >= 2**nFracTimes) then   ! .GT. 1.0E+05) THEN
          KONV = 2
          NIT = NIT + 1     ! DGM Nov 2006 - must increment or previous iter file will be overwritten
          RETURN
        ELSE
          GOTO 120
        END IF
      END IF

      IF (TMIN .LT. TOLREQ) THEN
        IF (abs(idebug) .GE. 1) WRITE(*,*) ' ...Finding Intercept...'
! TOLERANCE IS BELOW THAT REQUIRED; FIND INTERCEPT.
! THE LOWER VALUE OF MU BRACKETING THE INTERCEPT IS EASILY FOUND: IT IS
! JUST THE MINIMUM
        AMU1 = PMU
        TAMU1 = TMIN
! THE UPPER BOUND IS FOUND BY TESTING SUCCESSIVELY GREATER MU
        AMU2 = PMU
        TAMU2 = TMIN
        NFOR = 0
        do while (TAMU2 .LT. TOLREQ)
! SUCCESSIVELY DOUBLE MU:
          AMU2 = AMU2 + 0.30103 ! log10(2x) = log10(x) + (0.30103)
          TAMU2 = TOFMU(AMU2)
        ENDDO
        
        PMU2 = FROOT(TOFMU,AMU1,AMU2,TAMU1,TAMU2,TOLREQ,TOLI)
        
        TINT = TAMU2 + TOLREQ
        WRITE(IOUT,*) 'INTERCEPT IS AT MU = ',PMU2
        WRITE(IOUT,*) 'AND IS = ',TINT
        WRITE(IOUT,*) 'USING ',NFOR,' FUNCTION EVALUATIONS'
        
        IF (abs(idebug) .GE. 1) write(*,*) ' Lagrange MU= ', 10**PMU2
      ELSE
! TOLERANCE IS NOT YET SMALL ENOUGH.  WE WILL KEEP THE MINIMUM.
! SINCE fminocc RETURNS WITH PWK1() INSTEAD OF PWK2() WE NEED THIS COPY
        PWK2 = PWK1     ! array math
        dm = dwk1       ! array math
        IF (abs(idebug) .GE. 1) write(*,*) ' Lagrange MU= ', 10**PMU
      END IF
!***END LAGRANGE MULTIPLIER SELECTION

! COMPUTE ROUGHNESS.  WE DO THIS BY A FUNCTION CALL TO AVOID HAVING THE PENALTY
! MATRIX HANG AROUND TAKING UP MEMORY.
      IF (abs(idebug) .GE. 1) WRITE(*,*) 'Computing Roughness...'
      RUF = FNDRUF(IRUF)

! IF WE ATTAINED THE INTERCEPT LAST ITERATION BUT THE MODEL IS GETTING 
! ROUGHER WE HAVE PROBLEMS WITH CONERGENCE.
      IF ((RUF.GT.1.01*RLAST) .AND. (IFFTOL.EQ.1) .and. (idebug.ge.0)) &
     & THEN
        WRITE(IOUT,*) 'ROUGHNESS PROBLEMS, CUTTING STEP SIZE'
        FRAC = FRAC*2.0
! CHECK TO SEE IF ALL IS HOPELESS
        IF (FRAC >= 2**nFracTimes) then   ! .GT. 1.0E+05) THEN
            if (abs(idebug) .ge. 1) &
     &      write(*,*) 'ROUGHNESS PROBLEM NOT RESOLVED BY SMALL STEPS'
            KONV = 3
            NIT = NIT + 1     ! DGM Nov 2006 - must increment or previous iter file will be overwritten
            RETURN
        END IF
! OTHERWISE PLOW ON
        if (abs(idebug) .ge. 1) then
            write(*,*) 'ROUGHNESS GROWING! Cut step size & recalc MU'
            write(*,*) '     Roughness =', RUF
        endif
        GOTO 120
      END IF

! SAVE NEW MODEL AND COMPUTE STEP SIZE
      IF (abs(idebug) .GE. 1) WRITE(*,*) 'Computing Model Stepsize...'
      PWK3 = PWK2 - PM      ! Array math -- may have memory allocation drawback!
      PM = PWK2             ! Array math
! THE STEPSIZE IS THE ACTUAL CHANGE IN THE MODEL, NORMALIZED BY nParams:
      STEPSZ = SQRT(ANORM(nParams,PWK3)/nParams)
      WRITE(IOUT,*) 'STEPSIZE IS = ',STEPSZ
      WRITE(IOUT,*) 'ROUGHNESS IS = ',RUF
      WRITE(IOUT,*) ' '
! TIDY UP
      NIT = NIT + 1
      RLAST = RUF
      IF (TMIN .LT. TOLREQ) THEN
        IF (abs(idebug) .ge. 1 .and. IFFTOL == 0) &
     &      write(*,*) ' Tolerance met.  Next iteration begins smoothing.'
        IFFTOL = 1
        TOBT = TINT
      ELSE
        TOBT = TMIN
      END IF
! SEE IF WE HAVE A PERFECTLY SMOOTH MODEL THAT FITS DATA:
      IF (RUF .LT. 1.0E-5 .AND. TOBT .LE. 1.01*TOLREQ) KONV = 1
      IF (abs(idebug) .GE. 1) WRITE(*,*) 'Leaving Occam'
      RETURN
END SUBROUTINE OCCAM


!-----------------------------------------------------------------------
SUBROUTINE SCANMU(AM1,AM2,AM3,T2,T0)
!-----------------------------------------------------------------------

! OCCAM 3.0 Package
! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
! Subroutine Revision 3.0, November 2006
! Subroutine Revision 2.00, 13 May 1992

! SCANMU PRODUCES A REPORT OF MISFIT VERSUS LAGRANGE MULTIPLIER.
! IT CALLS THE FORWARD ROUTINE A LOT SO IS REALLY ONLY USED FOR DEBUGGING
! AND PARTICULARLY NASTY CASES
! IT LOOKS FOR THE GLOBAL MINIMUM AS FOLLOWS:
!   STARTING AT THE LARGEST MU USED FOR THE DEBUG REPORT, WE LOOK FOR THE 
!   FIRST MINIMUM.  WE ALSO LOOK FOR THE GLOLBAL MIN.  WE WANT THE GLOBAL 
!   MIN UNLESS THE LOCAL MIN AT LARGER MU IS LOWER THAN THE REQUIRED 
!   TOLERANCE. 

! ON INPUT:
!   T0 IS THE REQUIRED TOLERANCE
!   /RESULT/ IDEBUG IS THE PRINT LEVEL FROM OCCAM1
!   /RESULT/ IOUT IS A UNIT NUMBER FOR OUTPUT
! ON OUTPUT:
!   AM1, AM2, AM3 ARE THREE LAGRANGE MULTIPLIERS WHICH BRACKET A MINIMUM
!   T2 IS THE TOLERANCE AT AM2.
!   /RESULT/ ..PMW1().. HAS THE MODEL ASSOCIATED WITH AM2,T2.
! SUBROUTINES CALLED: 
!   TOFMU

! INCLUDES:
      USE Occam_Mod
      IMPLICIT NONE

! ARGUMENTS:
      REAL AM1,AM2,AM3,T2,T0
! LOCAL VARIABLES
      LOGICAL SEARCH
      INTEGER J,K
      REAL T1,TG2,T3,TL2,AML3,AML2,AML1,AMG3,AMG2,AMG1,AMU,T
! FUNCTIONS
      REAL TOFMU

!-----------------------------------------------------------------------
      WRITE(IOUT,*) 'MISFIT AS A FUNCTION OF MU:'
      T1 = 1.0E+09
      T2 = 1.1E+09
      TG2 = 0.9E+09
      AM2 = 16.
      AM1 = 16.
      SEARCH = .TRUE.

! START AT THE TOP (SMOOTHEST AND PRESUMABLY WORST FIT), SWEEPING FROM
! LOG10(MU) = 16 TO LOG10(MU) = -4 AT TWO VALUES PER DECADE
      DO 100 K = 16,-4,-1
! KEEP THE LAST TWO VALUES FOR MINIMUM SEARCH
        T3 = T2
        T2 = T1
        AM3 = AM2
        AM2 = AM1
! KEEP THE MODEL ASSOCIATED WITH THE MIDDLE VALUE IF WE ARE SEARCHING
! FOR THE MINIMUM
        IF (SEARCH) THEN
            PWK1 = PWK2     ! Array Math
            dwk1 = dm
        END IF
! KEEP THE MODEL ASSOCIATED WITH THE MIDDLE VALUE ANYWAY IN CASE
! GLOBAL AND LOCAL MINIMA ARE FOUND
        PWK3 = PWK2     ! Array math
        dwk3 = dm
! GET THE NEXT MU AND ASSOCIATED TOLERANCE
        AM1 = FLOAT(K)/2.
        T1 = TOFMU(AM1)

        IF ((T2 .LE. T1) .AND. (T2 .LE. T3) .AND. SEARCH) THEN
! WE HAVE FOUND THE FIRST MINIMUM. KEEP IT AND TURN OFF SEARCH
! (MODEL IS ALREADY STASHED IN PM1())
          SEARCH = .FALSE.
          TL2 = T2
          AML3 = AM3
          AML2 = AM2
          AML1 = AM1
        END IF
! KEEP AN EYE OUT FOR THE GLOBAL MIN IN CASE IT IS DIFFERENT
        IF (T2 .LE. TG2) THEN
          TG2 = T2
          AMG3 = AM3
          AMG2 = AM2
          AMG1 = AM1
! COPY ACROSS MODEL ASSOCIATED WITH GLOBAL MIN.
          PWK1 = PWK3       ! Array math
          dwk1 = dwk3
        END IF
! RECORD IN THE LOG FILE
        WRITE(IOUT,*) AM1,T1

! OUTPUT THE MODEL AS WELL IF REQUIRED
        IF (abs(idebug) .GE. 3) THEN
          DO J = 1,nParams
            WRITE(IOUT,*) '   ',PWK2(J)
          enddo
        END IF
100   CONTINUE
! WRITE OUT THE GAUSS STEP AT HIGHEST DEBUG LEVEL
      IF (abs(idebug) .GE. 3) THEN
        AMU = -99.0
        T = TOFMU(AMU)
        WRITE(IOUT,*) AM1,T1
        DO J = 1,nParams
          WRITE(IOUT,*) '   ',PWK2(J)
        enddo
      END IF

! NOW CHECK OUT THE MINIMA
!  IF THE FIRST MINIMUM ACHIEVES THE TOLERANCE REQUIRED WE WANT IT
      IF (TL2 .LE. T0) THEN
        T2 = TL2
        AM1 = AML1
        AM2 = AML2
        AM3 = AML3
      ELSE 
! OTHERWISE THE GLOBAL MIN IS JUST FINE
        T2 = TG2
        AM1 = AMG1
        AM2 = AMG2
        AM3 = AMG3
      END IF
      RETURN
END SUBROUTINE SCANMU


!-----------------------------------------------------------------------
REAL FUNCTION TOFMU(ALOGMU)
!-----------------------------------------------------------------------
! OCCAM 3.0 Package
! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
! Subroutine Revision 3.0, November 2006
! Subroutine Revision 2.01, 13 Jan 1993

! FUNCTION TOFMU RETURNS THE RMS MISFIT OF THE RESPONSE OF A MODEL 
! CONSTRUCTED USING THE GIVEN VALUE OF LAGRANGE MULTIPLIER.

! ON INPUT:
!   ALOGMU = LOG10 OF THE LAGRANGE MULTIPLIER 
!            (-99. IF A ZERO LAGRANGE MULTIPLIER NEEDED)
!   Stuff in Occam_mod

! ON OUTPUT
!   TOFMU = R.M.S. MISFIT OR 1.0E+10 IF CHOLESKY DECOMPOSITION FAILED
!   /RESULT/ PWK2() CONTAINS THE MODEL WHICH PRODUCES THE MISFIT TOFMU
!   /MODEL/ DM() CONTIANS THE RESPONSE OF THE MODEL 

! SUBROUTINES CALLED:
!   ANORM, FwdCalc ARE EXPLAINED IN OCCAM

!  AT VERSION 1.4 LOG10(MU) USED TO IMPROVE PERFORMANCE OF 1D OPTIMISATION

    ! INCLUDES:
    USE Occam_Mod
    use Occam_ModelLimits
    IMPLICIT NONE

    ! ARGUMENTS:
    REAL ALOGMU
    
    ! LOCAL VARIABLES
    REAL AMU, PART
    INTEGER ISTAT
!LAPACK        integer, dimension(:), allocatable :: idummy
    
    ! FUNCTIONS:
    REAL ANORM
      
!-----------------------------------------------------------------------
    IF (ALOGMU .EQ. -99.) THEN
        ! THE GAUSS STEP IS REQUIRED
        AMU = 0.0
    ELSE
        AMU = 10.**(ALOGMU)
    END IF
    NFOR = NFOR + 1
      
    ! ADD AMU.TRANS(PARTIAL).PARTIAL TO TRANS(W.J).W.J
    ! ALSO ADD MODEL PREJUDICE TO WJTWD()
    ! NB: memory allocaton drawback to the below array math.  Temp
    !     stack space is allocated for the rhs of real(nParams) size,
    !     possibly two or three times.
    PWK4  = WJTWD + AMU*PREWTS*PREMOD  ! Array math
    aMat  = 0.0
    aMat(:,1:nParams) = AMU * PTP + WJTWJ

    ! SOLVE THE LINEAR SYSTEM
    ! Eqn: [mu.dTd + (WJ)T.WJ] m = (WJ)T.Wd + mu*prejudice-wts
    !    aka:     aMat   *     m = PWK4
    !    NOTE: prejudice also included in diagonal of dTd (follow variable PTP)
    !         Prejudice almost always ZERO!
    ! Soln for m returned in PWK4 -- original value is overwritten.
    CALL CHOLIN(nParams,nParams,aMat,1,ISTAT)       ! Comment out if using LAPACK
! MAC Users:
!   To use the built-in Lapack libraries of OS/X (which are significantly
!   faster than B.Parker's Cholesky solver), do the following:
!
! - Change allocation of aMat in itrin() - remove the "+1" from dimension
! - Comment out the CHOLIN and CHOLSL calls (above & a little further down)
! - Uncomment the "!LAPACK" lines here and in the variable section above.
! - add "-framework vecLib" to the linking step of your make file.
!
!LAPACK    allocate(idummy(nParams))
!LAPACK    idummy = 0
!LAPACK    call sgesv( size(aMat,1), 1, aMat, size(aMat,1), idummy &
!LAPACK              , pwk4, size(aMat,1), ISTAT )
!LAPACK    deallocate(idummy)
!LAPACK    pwk2 = pwk4    ! Rest of code expects soln to be in pwk2.  So accommodate it.
    

    IF (ISTAT .NE. 0) THEN
        ! IF DECOMP FAILED WE GIVE WARNING BUT MIGHT STILL GET BY IF TOFMU
        ! IS SET TO A LARGE VALUE
        if (AMU > 1e20) then
            write(*,*) 'Error: Lagrange param went to infinity.  Stopping.'
            stop 'Error: Lagrange param went to infinity.  Stopping.'
        endif
        WRITE(*,*) 'Warning:- Cholesky decomp. failed at mu = ',AMU
        TOFMU = 1.0E+10

    else
        CALL CHOLSL(nParams,nParams,aMat,PWK2,PWK4,0)   ! Comment out if using LAPACK
        ! aMat is mtx of cholesky factors - upper right triangle + diagonal
        ! PWK4 is rhs of system of equations
        ! PWK2 is output: model
        ! 0=do full solution
        
        ! CUT STEP SIZE IF NECESSARY
        IF (FRAC .GT. 1.001) THEN
            PART = 1./FRAC
            PWK2 = (1.-PART)*PM + PART*PWK2     ! Array math - possible memory allocation problem
        END IF
        
        ! DGM Nov 2005
        ! If the system is configured to only allow discrete values
        !   in the model space, limit them now.  Do this *before* 
        !   applying mins & maxes as the rounding may affect this.
        if (gbModelStepped) then
            PWK2 = anint( PWK2 / gnModelSteps ) * gnModelSteps
        endif
        
        ! DGM Nov 2005
        ! If the system is configured for model limits, apply them.
        if (gbModelLimited) then
            where( PWK2 < gnModelMin ) PWK2 = gnModelMin
            where( PWK2 > gnModelMax ) PWK2 = gnModelMax
        endif
          
        ! CALCULATE THE MODEL RESPONSE AND MISFIT
        CALL FwdCalc(.false., nParams, ND, PWK2, DP, DM, D )
        DWK4 = (D - DM)/SD        ! Array math - stack alloc of only nd size
        TOFMU = SQRT(ANORM(ND,DWK4)/ND)
         
        IF (abs(idebug) .GE. 1) then
            write(*,*) '  Misfit:', tofmu, '  Mu:', amu
            WRITE(IOUT,*) 'TOFMU: MISFIT =', TOFMU,' AMU =', ALOGMU
        endif
          
        ! If the misfit is NaN then exit early.  Inversion has failed.
        !if (isnan(tofmu)) then
        if (huge(tofmu).eq.tofmu) then
            write(*,*) 'ERROR: NaN misfit.  Process aborted.'
            write(iout,*) 'ERROR: NaN misfit.  Process aborted.'
            stop
        endif
    endif
      
    RETURN
END FUNCTION TOFMU


!-----------------------------------------------------------------------
SUBROUTINE MAKJTJ(TOL0)
!-----------------------------------------------------------------------
! OCCAM 3.0 Package
! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
! Subroutine Revision 3.0, November 2006
! Subroutine Revision 2.00, 13 May 1992

! CONSTRUCTS MATRICES THAT DEPEND ON THE JACOBIAN (WJTWJ AND WJTWD).
! THE JACOBIAN ONLY NEEDS TO BE AROUND LONG ENOUGH TO DO
! THIS, SO BY STUFFING THIS INTO A SUBROUTINE DYNAMIC MEMORY ALLOCATION
! WILL MAKE WJ DISAPPEAR AFTER USE.

! ON INPUT:
!   CURRENT MODEL AND DATA STORED IN COMMON /MODEL/ AND /DATA/ 
! ON OUTPUT:
!   /RESULT/ WJTWJ = W.J.TRANS(W.J)
!   /RESULT/ WJTWD = WEIGHTED, TRANSLATED DATA PREMULTIPLIED BY TRANS(W.J)
!   TOL0 = RMS MISFIT ASSOCIATED WITH CURRENT MODEL
! CALLS:
!   ATAMUL, MULT, TRMULT, FwdCalc, ANORM

! INCLUDES:
      USE Occam_Mod
      IMPLICIT NONE

! ARGUMENTS:
      REAL TOL0
! LOCAL VARIABLES
      REAL, Allocatable :: WJ(:,:), DHAT(:)
      INTEGER I !, J
! FUNCTIONS 
      REAL ANORM

      Allocate( WJ(nd,nParams), DHAT(nd), STAT=I )
      if (I .ne. 0) then
        write(*,*) 'Out of memory.  Too many free parameters (', nParams, ')'
        stop 
      endif
!----------------------------------------------------------
! CALCULATE THE MATRIX OF PARTIALS AND MODEL RESPONSE
      call FwdCalc( .true., nParams, ND, PM, DP, DM, D, WJ )
      ! CALC MISFIT VECTOR AND MISFIT
      DWK1 = (D - DM) / SD      ! Array math - stack alloc of only nd size
      TOL0 = SQRT(ANORM(ND,DWK1)/ND)

      ! WEIGHT THE JACOBIAN MATRIX
      DO I = 1,ND
        WJ(I,:) = WJ(I,:) / SD(I)   ! Array math - possible alloc problem
      enddo

      ! FORM W.J.TRANS(W.J)
      CALL ATAMUL(ND,ND,nParams,nParams,WJ,WJTWJ)  
      
      ! FORM THE WEIGHTED, TRANSLATED DATA AND PREMULTIPLY BY TRANS(W.J)

     ! CALL MULT(ND,ND,nParams,nParams,1,WJ,PM,DHAT)
     
       DHAT = matmul(WJ,PM)
     
     
      DHAT = DHAT + DWK1

      CALL TRMULT(ND,ND,nParams,nParams,1,WJ,DHAT,WJTWD)
      
      Deallocate( WJ, DHAT )
      RETURN
END SUBROUTINE MAKJTJ


!-----------------------------------------------------------------------
SUBROUTINE MAKPTP(IRUF, model)
!-----------------------------------------------------------------------
! OCCAM 3.0 Package
! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
! Subroutine Revision 3.0, November 2006
! Subroutine Revision 2.00, 13 May 1992

! CONSTRUCTS PENALTY MATRIX (DEL) MULTIPLIED BY ITS TRANSPOSE.
! THE PENALTY MATRIX ONLY NEEDS TO BE AROUND LONG ENOUGH TO DO
! THIS, SO BY STUFFING THIS INTO A SUBROUTINE DYNAMIC MEMORY ALLOCATION
! CAN GET RID OF IT AFTERWARDS

! ON INPUT:
!   IRUF = see OBJMAT function
! ON OUTPUT:
!   /RESULT/ PTP() = DEL(TRANS)*DEL

    ! INCLUDES:
    USE Occam_Mod
    IMPLICIT NONE

    ! ARGUMENTS:
    INTEGER, intent(in)     :: IRUF
    real, intent(in)        :: model(nParams)
    
    ! LOCAL VARIABLES
    INTEGER                 :: nPTerms, i
    INTEGER, Allocatable    :: linkpr(:,:)
    REAL, Allocatable       :: DEL(:,:)
    
    ! Functions
    integer CountPenaltyTerms
    
!------------------------------------------------
    
    ! How many penalty terms are there?
    nPTerms = CountPenaltyTerms( IRUF )
    IF (abs(idebug) .GE. 1) then
        write(*,"('  # Penalty Terms = ',I8,' (Possible =',I8,')')") &
              nPTerms, nParams**2
    endif
    
    Allocate( DEL(nPTerms,nParams), linkpr(nPTerms,2), STAT=I )
    if (I .ne. 0) then
        write(*,*) 'Out of memory.  Too many free parameters (', nParams, ')'
        stop 
    endif

    ! MAKE PARTIALS MATRIX:
    CALL OBJMAT(IRUF, model, DEL, nPTerms, .true., linkpr)
    
    ! FORM TRANS(DEL).DEL
    call deltdel(nPTerms,nParams,del,ptp,linkpr)
    
    ! ADD THE WEIGHTS FOR THE MODEL PREJUDICE to the diagonal
    forall (I = 1:nParams, PREWTS(I) .ne. 0)        ! PREWTS is mostly zero!
        PTP(I,I) = PTP(I,I) + PREWTS(I)**2
    end forall

    Deallocate( DEL, linkpr )
    
    return
end SUBROUTINE MAKPTP


!-----------------------------------------------------------------------
subroutine deltdel(ma,na,a,c,lpair)
!-----------------------------------------------------------------------
! OCCAM 3.0 Package
! Catherine deGroot-Hedlin IGPP/SIO La Jolla CA 92093-0225
! Subroutine Revision 3.0, November 2006

! multiplies the transpose of del by del 
! takes advantage of sparseness of the del matrix 
! (ie multiply out only connected bricks)
      IMPLICIT NONE

      integer ma,na,i,j,kk,k,lpair(ma,2)
      real a(ma,na),c(na,na), cij

      c = 0.0
! connected pairs
      do kk=1,ma
         i=lpair(kk,1)
         j=lpair(kk,2)
         cij = 0.0
         do k=1,ma
            cij = a(k,i)*a(k,j) + cij
         enddo
         c(i,j) = cij
         c(j,i) = cij
      enddo
! diagonal component
      do i=1,na
         cij = 0.0
         do k=1,ma
            cij = a(k,i)*a(k,i) + cij
         enddo
         c(i,i) = cij
      enddo
      return
end subroutine deltdel


!-----------------------------------------------------------------------
REAL FUNCTION FNDRUF(IRUF)      
!-----------------------------------------------------------------------
! OCCAM 3.0 Package
! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
! Subroutine Revision 3.0, November 2006
! Subroutine Revision 2.00, 13 May 1992

! COMPUTES THE ROUGHNESS OF THE MODEL. FOR MEMORY ALLOCATION CONSIDERATIONS,
! THE PENALTY MATRIX HAS NOT BEEN KEPT AROUND, SO WE NEED TO REFORM IT.

! ON INPUT:
!   IRUF = see OBJMAT function
!   /RESULT/ PWK2() = LOCATION OF THE MODEL WE WANT TO MEASURE
! ON OUTPUT:
!   FNDRUF = ROUGHNESS OF MODEL
! CALLS:
!   OBJMAT, TRMULT

! INCLUDES:
    USE Occam_Mod
    IMPLICIT NONE

    ! ARGUMENTS:
    INTEGER IRUF
    
    ! Functions
    integer CountPenaltyTerms
    REAL ANORM
    
    ! LOCAL VARIABLES
    INTEGER I, nPTerms
    REAL, Allocatable :: DEL(:,:), VECT(:)
    
!------------------------------------------------
    ! How many penalty terms are there?
    nPTerms = CountPenaltyTerms( IRUF )
    Allocate( DEL(nPTerms,nParams), VECT(nPTerms), STAT=I )
    if (I .ne. 0) then
        write(*,*) 'Out of memory.  Too many free parameters (', nParams, ')'
        stop 
    endif

    ! MAKE PARTIALS MATRIX:
    ! Must use model1 (pm) not model2 (pwk2) to calculate the objective matrix
    ! or we get "roughness growing" troubles. This only affects modified
    ! roughness schemes.  Nothing else uses the 2nd param to ObjMat().
!    CALL OBJMAT(IRUF, pwk2, DEL, nPTerms, .false. )
    CALL OBJMAT(IRUF, pm, DEL, nPTerms, .false. )
    
    ! Calculate the new roughness (del * model2)
    DO I = 1,nPTerms
        VECT(I) = SUM( DEL(I,1:nParams)*PWK2(1:nParams) ) ! Array math
    ENDDO
    FNDRUF = ANORM( nPTerms, VECT )
      
    Deallocate( DEL, VECT )
      
    RETURN
END FUNCTION FNDRUF


!-----------------------------------------------------------------------
      SUBROUTINE MINBRK(AX,BX,CX,FA,FB,FC,FUNC)
!-----------------------------------------------------------------------
! OCCAM 3.0 Package
! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
! Subroutine Revision 3.0, November 2006
! Subroutine Revision 2.00, 13 May 1992
!
! MINBRK BRACKETS A UNIVARIATE MINIMUM OF A FUNCTION. TO BE USED PRIOR TO
! A UNIVARIATE MINIMISATION ROUTINE.
! BASED ON A ROUTINE IN NUMERICAL RECIPES BY PRESS ET AL
! MODIFIED SO THAT THE MODEL ASSOCIATED WITH THE MISFITS IS CARRIED AROUND
! FOR USE IN THE MINIMISATION ROUTINES AND POSSIBLY ULTIMATELY KEPT AS THE
! RESULT OF THIS ITERATION;
!
! ON INPUT:
!   AX, BX = TWO DISTINCT ESTIMATES OF THE MINIMUM'S WHEREABOUTS
!   FUNC = THE FUNCTION IN QUESTION
! ON OUTPUT:
!   AX,BX,CX = THREE NEW POINTS WHICH BRACKET THE MINIMUM
!   FA,FB,FC = FUNC(AX), FUNC(BX) ETC.
!       Apparently the minimum value is in Fb, such that Fb <= Fa and Fb < Fc.
!       No guaranteed relationship between Fa and Fc.
!   /RESULT/ PWK1() = THE MODEL ASSOCIATED WITH BX AND FB
!
! INCLUDES:
      USE Occam_mod
      IMPLICIT NONE
      EXTERNAL FUNC
!
! ARGUMENTS:
      REAL AX,BX,CX,FA,FB,FC
! LOCAL VARIABLES:
!      INTEGER I
      REAL DUM,R,Q,ULIM,U,FU
! LOCAL PARAMETERS
      REAL GOLD,GLIMIT,TINY
      PARAMETER (GOLD=1.618034, GLIMIT=100., TINY=1.E-32)
! FUNCTIONS:
      REAL FUNC
!
!-----------------------------------------------------------------------
      
      FB = FUNC(BX)
! MAKE PWK1 ASSOCIATED WITH B
      PWK1 = PWK2   ! Array math
      dwk1 = dm     ! Array math

! PM2 WILL BE ASSOCIATED WITH A
      FA = FUNC(AX)
      IF (FB.GT.FA) THEN
! SWITCH A AND B
        DUM = AX
        AX = BX
        BX = DUM
        DUM = FB
        FB = FA
        FA = DUM
! KEEP PWK1 ASSOCIATED WITH B (A WILL BE LOST SO DOESN'T MATTER)
        PWK1 = PWK2     ! Array math
        dwk1 = dm       ! Array math
      ENDIF
      ! At this point, Fb <= Fa
      
      CX = BX + GOLD*(BX - AX)
      FC = FUNC(CX)
! MAIN TEST. IF IT FAILS WE EXIT WITH B
! I.e. if Fb <= Fa AND < Fc then exit.
! Note: nothing inferred about Fc wrt Fa.
10    IF (FB.GE.FC) THEN
!      but in here: Fc <= Fb <= Fa
! KEEP PWK1 ASSOCIATED WITH C 
        PWK1 = PWK2     ! Array math
        dwk1 = dm       ! Array math
        ! Guess a new value (U) near Bx for Mu to try.
        R = (BX - AX)*(FB - FC)
        Q = (BX - CX)*(FB - FA)
        U = BX-((BX-CX)*Q-(BX-AX)*R)/(2.*SIGN(MAX(ABS(Q-R),TINY),Q-R))
        ULIM = BX + GLIMIT*(CX - BX)
        ! If "delta" (stuff subtracted from Bx above) is positive
        ! then this is IF U > Cx.
        ! If delta is negative then this is IF Cx > U
        ! I "think" this is trickery to account for AX,BX,CX being either directly
        ! or inversely related to the output of its functional FA,FB,FC
        IF ((BX - U)*(U - CX).GT.0.) THEN
          FU = FUNC(U)
          IF (FU.LT.FC) THEN
          ! Have: Fu < Fc <= Fb <= Fa
          ! Make: A=B & B=U (want Fb < Fc & Fb <= Fa)
            AX = BX
            FA = FB
            BX = U
            FB = FU
! KEEP PWK1 ASSOCIATED WITH B 
            PWK1 = PWK2     ! Array math
            dwk1 = dm       ! Array math
            GO TO 10
          ELSE IF (FU.GT.FB) THEN
          ! Have: Fc <= Fb < Fu and <= Fa (Fu may be >,=,< Fa)
          ! Make: Fc=U so Fb < Fc and Fb <= Fa
            CX = U
            FC = FU
            GO TO 10
          ENDIF
          ! Have: Fc <= Fu, Fb <= Fu
          ! So, take the low value (Fc) and try again
          U = CX + GOLD*(CX - BX)
          FU = FUNC(U)
          
        ELSE IF ((CX - U)*(U - ULIM).GT.0.) THEN
          FU = FUNC(U)
          IF (FU.LT.FC) THEN
            BX = CX
            CX = U
            U = CX + GOLD*(CX - BX)
            FB = FC
            FC = FU
! KEEP PWK1 ASSOCIATED WITH B 
            PWK1 = PWK2     ! Array Math
            dwk1 = dm       ! Array Math
            FU = FUNC(U)
          ENDIF
        ELSE IF ((U - ULIM)*(ULIM - CX).GE.0.) THEN
          U = ULIM
          FU = FUNC(U)
        ELSE
          U = CX + GOLD*(CX - BX)
          FU = FUNC(U)
        ENDIF
        
        ! It appears that a somewhat linear relationship is being assumed
        ! between successive values of out = func(in).  Above bits that 
        ! have fallen through have used Fc to generate U and then called
        ! FU = func(U).  And here we presume this generates a lower value
        ! than FC and "pop" A off the top.
        AX = BX
        BX = CX
        CX = U
        FA = FB
        FB = FC
        FC = FU
        GO TO 10
      ENDIF
      RETURN
END


!-----------------------------------------------------------------------
      REAL FUNCTION fminocc(AX,BX,CX,FBX,F,TOL,XMIN)
!-----------------------------------------------------------------------
! OCCAM 3.0 Package
! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
! Subroutine Revision 3.0, November 2006
! Subroutine Revision 2.00, 13 May 1992
!
! fminocc RETURNS THE MINIMUM VALUE OF A FUNCTION
! BASED ON A ROUTINE IN NUMERICAL RECIPES BY PRESS ET AL 
! MODIFIED SO THAT THE MODEL ASSOCIATED WITH THE MISFITS IS CARRIED AROUND
! FOR USE POSSIBLE ULTIMATE USE AS THE
! RESULT OF THIS ITERATION;
!
! ON INPUT:
!  AX,BX,CX = INDEPENDENT VARIABLE WHICH BRACKET THE MINIMUM
!  FBX = F(BX) (USUALLY AVAILABLE FROM THE BRACKETING PROCEDURE)
!  F = FUNCTION IN QUESTION
!  TOL = FRACTIONAL TOLERANCE REQUIRED IN THE INDEPENDENT VARIABLE
!  /RESULT/ PWK1 = MODEL ASSOCIATED WITH BX AND FBX
! ON OUTPUT:
!  XMIN = ABSCISSA OF MINIMUM
!  fminocc = F(XMIN)
!  /RESULT/ PWK1 = MODEL ASSOCIATED WITH XMIN AND fminocc
!
! INCLUDES:
      USE Occam_Mod
      IMPLICIT NONE
      EXTERNAL F
!
! ARGUMENTS:
      REAL AX,BX,CX,FBX,TOL,XMIN
! LOCAL VARIABLES
      REAL AA,B,V,W,X,E,FX,FV,FW,TOL1,Q,R,P,ETEMP,DD,FU,U,XM,TOL2
      INTEGER ITER !, I
! LOCAL PARAMETERS
      INTEGER ITMAX
      REAL CGOLD,ZEPS
      PARAMETER (ITMAX=100,CGOLD=.3819660,ZEPS=1.0E-10)
! FUNCTION:
      REAL F
!
!-----------------------------------------------------------------------
      AA = MIN(AX,CX)
      B = MAX(AX,CX)
      V = BX
      W = V
      X = V
      E = 0.
      FX = FBX
      FV = FX
      FW = FX
      DO 11 ITER = 1,ITMAX
        ! Calc a new value (U) to run the function on.
        XM = 0.5*(AA + B)
        TOL1 = TOL*ABS(X) + ZEPS
        TOL2 = 2.*TOL1
        IF (ABS(X - XM).LE.(TOL2 - .5*(B - AA))) then
            GOTO 30     ! exit!
        end if
        IF (ABS(E).GT.TOL1) THEN
          R = (X - W)*(FX - FV)
          Q = (X - V)*(FX - FW)
          P = (X - V)*Q - (X - W)*R
          Q = 2.*(Q - R)
          IF (Q.GT.0.) P =  - P
          Q = ABS(Q)
          ETEMP = E
          E = DD
          IF (ABS(P).GE.ABS(.5*Q*ETEMP).OR.P.LE.Q*(AA - X) .OR. &
     &        P.GE.Q*(B - X)) GOTO 1
          DD = P/Q
          U = X + DD
          IF (U - AA.LT.TOL2 .OR. B - U.LT.TOL2) DD = SIGN(TOL1,XM - X)
          GOTO 2
        ENDIF
1       IF (X.GE.XM) THEN
          E = AA - X
        ELSE
          E = B - X
        ENDIF
        DD = CGOLD*E
2       IF (ABS(DD).GE.TOL1) THEN
          U = X + DD
        ELSE
          U = X + SIGN(TOL1,DD)
        ENDIF
        
        ! Run the function
        FU = F(U)
        
        ! Evaluate the result - looking for a minimum
        IF (FU.LE.FX) THEN
          IF (U.GE.X) THEN
            AA = X
          ELSE
            B = X
          ENDIF
          V = W
          FV = FW
          W = X
          FW = FX
          X = U
          FX = FU
          PWK1 = PWK2   ! Array math
          dwk1 = dm     ! Array math
        ELSE
          IF(U.LT.X) THEN
            AA = U
          ELSE
            B = U
          ENDIF
          IF (FU.LE.FW .OR. W.EQ.X) THEN
            V = W
            FV = FW
            W = U
            FW = FU
          ELSE IF (FU.LE.FV .OR. V.EQ.X .OR. V.EQ.W) THEN
            V = U
            FV = FU
          ENDIF
        ENDIF
11    CONTINUE
      WRITE(*,*) 'MAXIMUM ITERATIONS EXCEEDED IN fminocc'
      
30    XMIN = X
      fminocc = FX
      RETURN
END


!-----------------------------------------------------------------------
      REAL FUNCTION FROOT(FUNC,X1,X2,FA,FB,OFF,TOL)
!-----------------------------------------------------------------------
! OCCAM 3.0 Package
! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
! Subroutine Revision 3.0, November 2006
! Subroutine Revision 2.00, 13 May 1992
!
! FROOT FINDS THE POINT AT WHICH A UNIVARIATE FUNCTION ATTAINS A GIVEN
!  VALUE
! BASED ON A ROUTINE IN NUMERICAL RECIPES BY PRESS ET AL 
! MODIFIED SO THAT THE MODEL ASSOCIATED WITH THE MISFITS IS CARRIED AROUND
! FOR USE POSSIBLE ULTIMATE USE AS THE
! RESULT OF THIS ITERATION;
!
! ON INPUT:
!  FUNC = THE FUNCTION IN QUESTION
!  X1,X2 = INDEPENDENT VARIABLES BRACKETING THE ROOT
!  FA,FB = FUNC(X1),FUNC(X2) (USUALLY AVAILABLE FROM THE BRACKETING)
!  OFF = VALUE OF THE FUNCTIONAL AT THE REQUIRED ABSCISSA
!  TOL = FRACTIONAL TOLERANCE REQUIRED FOR ABSCISSA
!  /RESULT/ PWK2 = MODEL ASSOCIATED WITH X2 AND FB
!  /RESULT/ PWK1 = MODEL ASSOCIATED WITH X1 AND FA
!
! ON OUTPUT:
!  FB = FUNC(FROOT) - OFF
!  FROOT = ABSCISSA REQUIRED
!  /RESULT/ PWK2 = MODEL ASSOCIATED WITH FB, FROOT
!
! INCLUDES:
      USE Occam_Mod
      IMPLICIT NONE
      EXTERNAL FUNC
!
! ARGUMENTS:
      REAL X1,X2,FA,FB,OFF,TOL
! LOCAL VARIABLES
      REAL AA,B,FC,C,E,DD,TOL1,P,Q,R,S,XM
      INTEGER ITER !, I
! LOCAL PARAMETERS
      INTEGER ITMAX
      REAL EPS
      PARAMETER (ITMAX = 100,EPS = 3.E-8)
! FUNCTIONS:
      REAL FUNC
!
!-----------------------------------------------------------------------
      AA = X1
      B = X2
      FA = FA - OFF
      FB = FB - OFF
      IF (FB*FA.GT.0.) WRITE(*,*) 'ROOT NOT BRACKETED IN FROOT'
      FC = FB
! MODELS WILL BE CARRIED AROUND AS:
!   AA: PWK1()
!   B: PWK2()
!   C: PWK3()
! SO COPY PM2 (WITH B) INTO WK (WITH C)
      PWK3 = PWK2   ! Array math
      dwk3 = dm     ! Array math
      DO 11 ITER = 1,ITMAX
! IF B AND C ARE NOT ON DIFFERENT SIDES OF ROOT THEN MAKE THEM  SO
        IF (FB*FC.GT.0.) THEN
          C = AA
          FC = FA
          DD = B - AA
          E = DD
! COPY PMF (WITH AA) INTO WK (WITH C)
          PWK3 = PWK1   ! Array math
          dwk3 = dwk1   ! Array math
        ENDIF
        IF (ABS(FC).LT.ABS(FB)) THEN
! DGM Oct 2006 - found this code was rotating values
! using the old A=B, B=C, C=A bug.  However, SC confirmed that
! this is the way it is in Numerical Recipes.  Inversion test
! showed that it made no difference.  So leaving it the confusing
! way just in case it does make a difference in some odd case.
! IF B IS NOT CLOSER THAN C MAKE IT SO
          AA = B
          B = C
          C = AA
          FA = FB
          FB = FC
          FC = FA
! SWITCH PM2 (WITH B) AND WK (WITH C) (AA WILL BE LOST IF WE DON'T EXIT)
          PWK1 = PWK2   ! Array math
          PWK2 = PWK3
          PWK3 = PWK1
          dwk1 = dm
          dm   = dwk3
          dwk3 = dwk1
        ENDIF
        TOL1 = 2.*EPS*ABS(B) + 0.5*TOL
        XM = .5*(C - B)
        IF (ABS(XM).LE.TOL1 .OR. FB.EQ.0.)  THEN
          FROOT = B
          RETURN
        ENDIF
        IF (ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
          S = FB/FA
          IF (AA.EQ.C) THEN
            P = 2.*XM*S
            Q = 1. - S
          ELSE
            Q = FA/FC
            R = FB/FC
            P = S*(2.*XM*Q*(Q - R) - (B - AA)*(R - 1.))
            Q = (Q - 1.)*(R - 1.)*(S - 1.)
          ENDIF
          IF (P.GT.0.) Q =  - Q
          P = ABS(P)
          IF (2.*P .LT. MIN(3.*XM*Q - ABS(TOL1*Q),ABS(E*Q))) THEN
            E = DD
            DD = P/Q
          ELSE
            DD = XM
            E = DD
          ENDIF
        ELSE
          DD = XM
          E = DD
        ENDIF
        AA = B
        FA = FB
! COPY PM2 (WITH B) INTO PMF (WITH AA)
        PWK1 = PWK2     ! Array math
        dwk1 = dm
        IF (ABS(DD) .GT. TOL1) THEN
          B = B + DD
        ELSE
          B = B + SIGN(TOL1,XM)
        ENDIF
! THE NEW FUNCION EVALUATION AUTOMATICALLY USES PM2 FOR B
        FB = FUNC(B) - OFF
11    CONTINUE
      WRITE(*,*) 'MAXIMUM ITERATIONS EXCEEDED IN FROOT'
      FROOT = B
      RETURN
END


!-----------------------------------------------------------------------
      SUBROUTINE CHOLSL(LA, N, A, X, B, ISTAGE)
!-----------------------------------------------------------------------
! OCCAM 3.0 Package
! R.L. Parker IGPP/SIO La Jolla CA 92093-0225
! Subroutine Revision 3.0, November 2006
!
!$$$$$ CALLS NO OTHER ROUTINES
!  ROUTINE FOR BACK-SUBSTITUTION SOLUTION OF A LINEAR SYSTEM THAT HAS
!  ALREADY BEEN FACTORED INTO CHOLESKY FACTORS BY 'CHOLIN' (Q.V.)
!  LA      LEADING DIMENSION OF A
!  N       THE ORDER OF MATRIX
!  A       ARRAY CONTAINING ON AND ABOVE ITS DIAGONAL THE TRANSPOSE OF
!          THE CHOLESKY FACTOR.  THE NORMAL WAY OF FINDING THE FACTORS
!          IS TO CALL  CHOLIN  AS FOLLOWS
!     CALL CHOLIN(N, A, 1, IFBAD)
!     CALL CHOLSL(N, A, X, B, 0)
!  X       THE UNKNOWN VECTOR.  SOLUTION TO EQUATIONS IS RETURNED HERE
!  B       THE RIGHT SIDE OF SYSTEM
!  ISTAGE  THE MODE OF THE SOLUTION.  IF=1 DO NOT COMPLETE BACK-
!     SUBSTITUTION.  INSTEAD RETURN INVERSE(L-TRANS)*B .  OTHER VALUES
!     FORM THE COMPLETE SOLUTION.
! R.L. PARKER
!
! INCLUDES:
      IMPLICIT NONE
! ARGUMENTS
      INTEGER LA,N,ISTAGE
      REAL A(LA,N+1),X(N),B(N)
! LOCAL VARIABLES
      INTEGER I,II !,K,KK
      REAL Z
!
!-----------------------------------------------------------------------
    DO I=1,N
        Z=B(I)
        IF (I .ne. 1) then
            Z = Z - sum( A(1:I-1,I+1) * X(1:I-1) )
!            DO KK=2,I
!                K=I-KK+1
!                Z=Z - A(K,I+1)*X(K)
!            enddo
        endif
        X(I)=Z/A(I,I+1)
    enddo
    IF (ISTAGE.EQ.1) RETURN
    
!  COMPLETE BACK-SUBSTITUTION
    DO II=1,N
        I=N-II+1
        Z=X(I)
        if (I .ne. N) then
            Z = Z - sum( A(I,I+2:N+1) * X(I+1:N) )
!            DO K=I+1,N
!                Z= Z - A(I,K+1)*X(K)
!            enddo
        end if
        X(I)=Z/A(I,I+1)
    enddo
    RETURN
END


!-----------------------------------------------------------------------
      SUBROUTINE CHOLIN(LA,N, A, JSTAGE, IFBAD)
!-----------------------------------------------------------------------
! OCCAM 3.0 Package
! R.L. Parker IGPP/SIO La Jolla CA 92093-0225
! Subroutine Revision 3.0, November 2006
!
!$$$$$ CALLS NO OTHER ROUTINES
!  INVERSION AND FACTORIZATION OF POSITIVE DEFINITE SYMMETRIC MATRICES
!  FROM ALGOL ROUTINE OF WILKINSON & REINSCH 'LINEAR ALGEBRA' PP16-17.
!  EXPLANATION
!  CHOLIN  CALCULATES THE FACTOR  L  IN THE CHOLESKY FACTORIZATION
!  A = L L-TRANS, OR THE INVERSE OF L , OR THE INVERSE OF  A.
!  TRANSPOSED STORAGE IS FOR THE CONVENIENCE OF FORTRAN USERS.  HOWEVER,
!  IT MEANS THAT L-TRANS AND INVERSE OF L-TRANS ARE COMPUTED.
!  ARGUMENTS
!    A     AN ARRAY DIMENSION N,N+1.  ON ENTRY THE INPUT MATRIX SHOULD
!          SHOULD OCCUPY THE 1ST  N COLS, OR IF DESIRED, JUST THE LOWER
!          LOWER LEFT TRIANGLE.  ON EXIT THE UPPER TRIANGLE OF THE ARRAY
!          DEFINED BY COLS 2,3...N+1 CONTAINS THE DESIRED RESULT.
!          THIS CAN BE CONVENIENTLY ACCESSED AS A(1,2), OR A(N+1) IF
!          SINGLY DIMENSIONED IN CALLING PROGRAM.
!    N     THE ORDER OF THE MATRIX
!    LA    LEADING DIMENSION OF A
!  JSTAGE  INDICATES DESIRED CALCULATION - IF=1 UPPER TRIANGLE CONTAINS
!          L-TRANS.  IF=2 UPPER TRIANGLE CONTAINS INVERSE OF L-TRANS.
!          IF=0 OR 3 CONTAINS THE UPPER HALF OF THE INVERSE OF A.  NOTE
!          LOWER TRIANGLE IS NOT DESTROYED. IF=-2,-3, THE CHOLESKY
!          FACTORS ARE ASSUMED TO BE IN PLACE ON ENTRY.
!  IFBAD   ERROR FLAG, SET TO 1 IF MATRIX NOT POSITIVE DEFINITE. SET TO
!          0 IF OK.
! R.L. PARKER
!
    IMPLICIT NONE
    
    ! PASS PARAMETERS
    INTEGER LA,N,JSTAGE,IFBAD
    REAL A(LA,N+1)
    
    ! LOCAL VARIABLES
    INTEGER I,J,K,KBACK,ISTAGE,I1,J1,J2,N1
    REAL X,Y,Z
!
!-----------------------------------------------------------------------
!  STAGE 1 - FORM FACTORS (HERE L-TRANS IS STORED ABOVE DIAGONAL)
    IFBAD  = 0
    ISTAGE = IABS(JSTAGE)
    IF (JSTAGE.GE.0) THEN
        DO I=1,N
            I1=I + 1
            DO J=I,N
                J1=J + 1
                X=A(J,I)
                IF (I.NE.1) THEN
                    ! DGM: Array math here is found to be 25% faster in tests!
                    ! Drawback: stack space allocated equal to 3 times A(1:I-1,1)
                    !   - once for each sub-array, then for the sum.
                    X = X - SUM( A(1:I-1,J1)*A(1:I-1,I1) )
!                    DO KBACK=2,I
!                        K=I-KBACK+1
!                        X=X - A(K,J1)*A(K,I1)
!                    END DO
                ENDIF
                IF (I.NE.J) THEN
                    A(I,J1)=X*Y
                ELSE
                    IF (X.LE.0.0) GO TO 4000
                    Y=1.0/SQRT(X)
                    A(I,I1)=Y
                ENDIF
            ENDDO
        ENDDO
        IF (ISTAGE.NE.1) GO TO 2000
    ENDIF
    
    DO I=1,N
        A(I,I+1)=1.0/A(I,I+1)
    ENDDO
    IF (ISTAGE.EQ.1) RETURN
    
    ! NB: As of Oct 2006, rest of routine never used by Occam2D.
    ! Always called with JSTAGE==1
    
!  STAGE 2 - INVERSION OF L-TRANS
 2000 DO 2600 I=1,N
      I1=I+1
      IF (I1.GT.N) GO TO 2600
      DO 2500 J=I1,N
      Z=0.0
      J1=J + 1
      J2=J - 1
      DO 2200 KBACK=I,J2
      K=J2-KBACK+I
 2200 Z=Z - A(K,J1)*A(I,K+1)
      A(I,J1)=Z*A(J,J1)
 2500 CONTINUE
 2600 CONTINUE
      IF (ISTAGE.EQ.2) RETURN
!  STAGE 3  - CONSTRUCTION OF INVERSE OF  A  ABOVE DIAGONAL.
      DO 3500 I=1,N
      DO 3500 J=I,N
      Z=0.0
      N1=N + 1
      J1=J + 1
      DO 3200 K=J1,N1
 3200 Z=Z + A(J,K)*A(I,K)
      A(I,J+1)=Z
 3500 CONTINUE
      RETURN
!  ERROR EXIT
 4000 IFBAD=1
! ERROR MESSAGE SUPRESSED - S.CONSTABLE
!      WRITE(*,400)
! 400  FORMAT(44H0CHOLIN FAILS - MATRIX NOT POSITIVE DEFINITE           )
      RETURN
END  

!-----------------------------------------------------------------------
      SUBROUTINE TRMULT(MAD,MA,NAD,NA,NB,A,B,C)
!-----------------------------------------------------------------------
! OCCAM 3.0 Package
! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
! Subroutine Revision 3.0, November 2006
! MULTIPLIES THE TRANSPOSE OF A 2D MATRIX BY ANOTHER MATRIX
      IMPLICIT NONE
!
      INTEGER MAD, MA, NAD, NA, NB, I, J !, K
      REAL A(MAD,NA),B(MAD,NB),C(NAD,NB) !, CIJ
      DO I = 1,NA
        DO J = 1,NB
           C(I,J) = SUM( A(1:MA,I) * B(1:MA,J) )
        enddo
      enddo
      RETURN
END


!-----------------------------------------------------------------------
      REAL FUNCTION ANORM(N,D)
!-----------------------------------------------------------------------
! OCCAM 3.0 Package
! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
! Subroutine Revision 3.0, November 2006
! RETURNS THE SQUARE OF THE EUCLIDEAN NORM OF A VECTOR
      IMPLICIT NONE
!
      INTEGER N 
      REAL D(N)
!
      ANORM = SUM( D * D )
      RETURN
END


!-----------------------------------------------------------------------
      SUBROUTINE ATAMUL(MAD,MA,NAD,NA,A,C)
!-----------------------------------------------------------------------
      IMPLICIT NONE
! OCCAM 3.0 Package
! C. deGroot-Hedlin IGPP/SIO La Jolla CA 92093-0225
! Subroutine Revision 3.0, November 2006
!
! DESCRIPTION: 
!
! MULTIPLIES THE TRANSPOSE OF A 2D MATRIX BY ITSELF
! TAKES ADVANTAGE OF SYMMETRY ABOUT DIAGONAL
!
!
! ARGUMENTS
      INTEGER MAD,MA,NAD,NA
      REAL A(MAD,NA),C(NAD,NA)
!
! LOCAL VARIABLES
      INTEGER I,J !,K
!      REAL CIJ
!-----------------------------------------------------------------------
! NOTE: Cannot simply do the below.  It allocs A-size on the stack 
!       twice: once as transpose, then as output of MATMUL.
!    C = MATMUL( TRANSPOSE(A), A )
!    RETURN
! This only allocs one A-size, but this is still too much when 
!   A is very large - which it usually is.
!    C = TRANSPOSE(A)
!    C = MATMUL( C, A )
!    RETURN
    
      DO I = 1,NA
        DO J = I,NA
          ! NOTE assumption of symmetry!!
          C(I,J) = SUM( A(1:MA,I)*A(1:MA,J) )
        ENDDO
      ENDDO
      ! Copy into the symmetric part of the matrix
      DO I=2,NA
        C(I,1:I-1) = C(1:I-1,I)
      ENDDO
      RETURN
END


!-----------------------------------------------------------------------
      subroutine itrin (itform,descr,modfil,datfil,maxitr,tolreq, &
     &    iruf,nit,pmu,rlast,tobt,ifftol,iof2)      
!-----------------------------------------------------------------------
! OCCAM 3.0 Package
! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
! Subroutine Revision 3.0, November 2006
! Subroutine Revision 2.01, 13 Jan 1993
! DGM Nov 2006 - new startup file type (OCCAMITER_FLEX) allows
!   flexible placement & specification of header items.
!
! Reads a STARTUP or ITER file for OCCAM.
!
!
    ! includes:
    USE Occam_Mod
    use Occam_ModelLimits

    IMPLICIT NONE

      
    ! arguments:
    integer maxitr, iruf, nit, ifftol,iof2
    real tolreq, pmu, rlast, tobt
    character(50) modfil, datfil, itform, descr
    ! local variables
    integer i, idcode, nAllocErr
!   i = loop index
!   idcode = function to perform the equivalent of a free-format internal
!            read for an integer
    real rdcode
!   rdcode = function to perform the equivalent of a free-format internal
!            read for a real
    character(180)  sLine, sCode, sValue
    logical         bComment
!
!-----------------------------------------------------------------------
    ! Set default values to the file params so that any not in the file
    !   which are optional don't cause a crash.
    maxitr      = 10
    tolreq      = 1.0
    iruf        = 1
    gbAddDiags  = .false.
    idebug      = 1     ! print status msgs by default.  Everyone likes them.
    nit         = 0
    pmu         = 5.0
    rlast       = 1.0e+10
    tobt        = 1000
    ifftol      = 0
    nParams     = 0
    
    call ModelLimits_Init
    
    ! Read the file header.
    do while (.true.) !KWK nov 16,2006 eof not supported generally (.not. EOF(iof2))
        ! Get the next code/value pair
        ! ParseCode forces the code portion to be all lowercase with no
        !   padding and ending colon stripped off.  User comments are 
        !   stripped from the value portion.
        read( iof2, '(A)', end=198, err=199 ) sLine
        call ParseCode( len(sLine), sLine, sCode, sValue, bComment )
        if( bComment ) cycle
        
        ! What do we have?
        select case (trim(sCode))
        case ('format')
            itform = sValue(1:len(itform))
            
        case ('description')
            descr = sValue(1:len(descr))
            
        case ('model file')
            modfil = sValue(1:len(modfil))
            
        case ('data file')
            datfil = sValue(1:len(datfil))
            
        case ('date/time')
            ! ignored!
            
        case ('max iter', 'iterations to run')
            maxitr = idcode( sValue, iof2 )
            
        case ('req tol', 'target misfit')
            tolreq = rdcode( sValue, iof2 )
            
        case ('iruf', 'roughness type')
            iruf = abs( idcode( sValue, iof2 ) )
            
        case ('diagonal penalties')
            if( idcode( sValue, iof2 ) == 0 ) then
                gbAddDiags  = .false.
            else
                gbAddDiags  = .true.
            endif
            
        case ('debug level')
            idebug = idcode( sValue, iof2 )
            
        case ('iteration')
            ! Actually the # of the previous iteration - so that at startup an iter00.iter
            !   file will be written reflecting the startup values.
            nit = idcode( sValue, iof2 )
            
        case ('pmu', 'lagrange value')
            pmu = rdcode( sValue, iof2 )
            
        case ('rlast', 'roughness value')
            rlast = rdcode( sValue, iof2 )
            
        case ('tlast', 'misfit value')
            tobt = rdcode( sValue, iof2 )
            
        case ('ifftol', 'misfit reached')
            ifftol = idcode( sValue, iof2 )
            
        case ('stepsize cut count')
            nFracTimes = idcode( sValue, iof2 )
            
        case ('model limits')   ! DGM Nov 2006 see module Occam_ModelLimits
            ! Limit value string must have two values: min, max
            ! Comma is required.
            i = index( sValue, ',' )
            if (i == 0) goto 196
            gnModelMin = rdcode( sValue(:i-1), iof2 )
            
            ! Skip space padding
            do i=i+1,len_trim(sValue)
                if (sValue(i:i) .ne. ' ') exit
            enddo
            gnModelMax = rdcode( sValue(i:len_trim(sValue)), iof2 )
            
            ! Min & Max must be reasonable
            if (gnModelMax < gnModelMin .or. gnModelMax == gnModelMin) then
                goto 196
            endif
            
            gbModelLimited = .true.
            write(*,*) '  Limiting model values from', gnModelMin, 'to', gnModelMax
            
        case ('model value steps')   ! DGM Nov 2006 see module Occam_ModelLimits
            gnModelSteps = rdcode( sValue, iof2 )
            if (gnModelSteps > 0) then
                gbModelStepped  = .true.
                write(*,*) '  Limiting model values to steps of', gnModelSteps
            endif
            
        case ('modify roughness')
            ! Left-over from previous investigation of ways to make 
            ! inversion sharper.  Any startup file with this tag in
            ! it will no longer function properly.  DGM Dec 2006.
            write(*,*) '  Statement: Modify Roughness :no longer supported'
            
            
        case ('no. parms', 'param count')
            nParams = idcode( sValue, iof2 )
            exit    ! The next thing in the file is the data param list!
            
        case default
            write(*,*) 'Error reading startup/iteration file!'
            write(*,*) 'Unknown or unsupported code:', sCode
            stop
        end select
    enddo
    
    ! Did we end early without a params statement??
    if (nParams <= 0) then
        write(*,*) 'Cannot have zero parameters.'
        write(*,*) 'Is header of startup file formatted properly?'
        stop
    endif
    
    ! Write out some values
    write(*,*) '  Target Misfit:', tolreq
    write(*,*) '  Converged Flag:', ifftol
    
    if (gbAddDiags) then
        write(*,*) '  Adding diagonal penalties to horizontal & vertical.'
    endif
    
    ! DGM Oct 2006 - allocate what we need for the rest of this routine.
    ! Most vars allocated in the main program.
    ! If we don't survive allocation, explain nicely & stop the program.
    allocate( pm(nParams)             &
            , ptp(nParams,nParams)    &
            , wjtwj(nParams,nParams)  &
            , wjtwd(nParams)          &
            , pwk1(nParams)           &
            , pwk2(nParams)           &
            , pwk3(nParams)           &
            , prewts(nParams)         &
            , premod(nParams)         &
            , aMat(nParams,nParams+1) &   ! extra dim for use in cholin(); remove "+1" if using LAPACK
            , pwk4(nParams)           &   ! work var or TOFMU
            , stat=nAllocErr )   

    if (nAllocErr .ne. 0) then
        write(*,*) 'Out of memory.  Too many free parameters (', nParams, ')'
        stop 
    endif
      
    ! Read the data
    read (unit=iof2,fmt= *, end=198, err=199) (pm(i), i=1,nParams)
    
      
    return
! Various error conditions
196 continue
    write(*,*) 'Invalid "Model Limits" entry in startup file.'
    write(*,*) 'The format for this line is'
    write(*,*) '   Model Limits: <min>,<max>'
    write(*,*) 'Do not include angle brackets.  Comma required.'
    write(*,*) 'Min must be less than max.'
    stop
    
198 call filerr (' Startup file ended prematurely', iof2)
199 call filerr (' Error reading startup file', iof2)

end subroutine itrin


!-----------------------------------------------------------------------
subroutine itrout (itform,descr,modfil,datfil,maxitr,tolreq, &
     &    iruf,nit,pmu,rlast,tobt,ifftol,iof2,cRoot)      
!-----------------------------------------------------------------------
! OCCAM 3.0 Package
! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
! Subroutine Revision 3.0, November 2006
! Subroutine Revision 2.00, 13 May 1992
! Nov 2006 - DGM updated to new FLEX format & to add roughness setup.
!
! writes an iteration file for OCCAM
!
    ! includes:
    USE Occam_Mod
    use Occam_ModelLimits
    use Fwd_mod
    IMPLICIT NONE
    
    ! arguments:
    integer maxitr, iruf, nit, ifftol,iof2
    real tolreq, pmu, rlast, tobt
    character(50) modfil, datfil, itform, descr, datetm, cRoot
    ! local variables
    integer lerr, i, j
    !   lerr = error flag on file operation
    !   i = loop index
    character(2) cNum
    
!-----------------------------------------------------------------------
    ! create file name suffix
    write (cNum, '(i2)') nit
    if (nit .le. 9) cNum(1:1) = '0'
    
    ! open iteration file
    if (cRoot == '') then
        open (unit=iof2, file='ITER'//cNum//'.iter', iostat=lerr)
    else
        open (unit=iof2, file=trim(cRoot)//cNum//'.iter', iostat=lerr)
    end if
    if (lerr .ne. 0) then
        ! we could be in the middle of an iteration, so don't die on error, but
        ! send a message and forge on regardless
        write(*,*) ' Error opening iteration file'
        return
    end if
    
    ! Get the current date & time
    call datime(datetm)

    ! Write the header
    write(iof2,'(a)')       'Format:             OCCAMITER_FLEX'
    write(iof2,'(a,a)')     'Description:        ', trim(descr)
    write(iof2,'(a,a)')     'Model File:         ', trim(modfil)
    write(iof2,'(a,a)')     'Data File:          ', trim(datfil)
    write(iof2,'(a,a)')     'Date/Time:          ', trim(datetm)
    write(iof2,'(a,i3)')    'Iterations to run:  ', maxitr
    write(iof2,'(a,f6.3)')  'Target Misfit:      ', tolreq
    write(iof2,'(a,i3)')    'Roughness Type:     ', abs(iruf)
    if (gbAddDiags) then
        write(iof2,'(a,i3)')'Diagonal Penalties: 1'
    else
        write(iof2,'(a,i3)')'Diagonal Penalties: 0'
    endif
    write(iof2,'(a,i3)')    'Stepsize Cut Count: ', nFracTimes
    if (gbModelLimited) then
        write(iof2,'(a,f8.6,",",f8.6)') &
                            'Model Limits:       ', gnModelMin, gnModelMax
    else
        write(iof2,'(a)')   '!Model Limits:        min,max'
    endif
    if (gbModelStepped) then
        write(iof2,'(a,f8.6)') &
                            'Model Value Steps:  ', gnModelSteps
    else
        write(iof2,'(a)')   '!Model Value Steps:   stepsize (e.g. 0.2 or 0.1)'
    endif
    write(iof2,'(a,i3)')    'Debug Level:        ', idebug
    write(iof2,'(a,i3)')    'Iteration:          ', nit
    write(iof2,'(a,g15.7)') 'Lagrange Value:     ', pmu
    write(iof2,'(a,g15.7)') 'Roughness Value:    ', rlast
    write(iof2,'(a,g15.7)') 'Misfit Value:       ', tobt
    write(iof2,'(a,i3)')    'Misfit Reached:     ', ifftol
    
    ! Write the end of the header and all of the param data.
    write(unit=iof2,fmt='(a,i4)')    'Param Count:        ', nParams
    write(unit=iof2,fmt='(4(2x,g15.7))') (pm(i), i=1,nParams)
    close(unit=iof2)
       
    ! open response file
    if (cRoot == '') then
        open (unit=iof2, file='RESP'//cNum//'.resp', iostat=lerr)
    else
        open (unit=iof2, file=trim(cRoot)//cNum//'.resp', iostat=lerr)
    end if
    if (lerr .ne. 0) then
        ! we could be in the middle of an iteration, so don't die on error, but
        ! send a message and forge on regardless
        write(*,*) ' Error opening response file'
        return
    end if
    do i = 1, nd
        write(iof2,'(7(2x,g13.5))') &
            (dp(i,j),j=1,cnDataParams), d(i), dm(i), (d(i)-dm(i))/sd(i)
    end do
    close(iof2)
    return
    
end subroutine itrout


!-----------------------------------------------------------------------
subroutine ParseCode( nLen, sLine, sCode, sValue, bComment )
!-----------------------------------------------------------------------
! OCCAM 3.0 Package
! David Myer IGPP/SIO La Jolla CA 92093-0225
! Subroutine Revision 3.0, November 2006
! DGM Nov 2006 - parse a line read from a file into a code & value.
! Force the code to be all lowercase with no ending colon.  Terminate
! the line at a '%' or '!' sign (these allow for user comments!)
!
    implicit none
    
    ! Args
    integer, intent(in)         :: nLen
    character(nLen)             :: sLine
    character(nLen), intent(out) :: sCode, sValue
    logical, intent(out)        :: bComment
    
    ! Local vars
    integer :: iFrom, iTo
    
    ! Init returns
    bComment = .false.
    sCode = ' '
    sValue = ' '
    
    ! Convert all tab characters to spaces
    forall( iTo = 1:nLen, ichar(sLine(iTo:iTo)) == 9 ) sLine(iTo:iTo) = ' '
    
    ! Skip any beginning blanks
    do iFrom = 1,nLen
        if (sLine(iFrom:iFrom) .ne. ' ') exit
    enddo
    
    ! If the first char is a comment char, then the whole line is a comment.
    if( sLine(iFrom:iFrom) == '%' .or. sLine(iFrom:iFrom) == '!' ) then
        bComment = .true.
        return
    endif
    
    ! Pull off the code value.  Cvt to lowercase as we go.
    iTo = index(sLine,':') - 1
    if (iTo < iFrom) then
        write(*,*) 'Parsing Error: missing colon in line below:'
        write(*,*) sLine
        return
    endif
    sCode = sLine(iFrom:iTo)
    call Lower(sCode)
    
    ! Skip spaces after the colon
    do iFrom = iTo+2,nLen
        if (sLine(iFrom:iFrom) .ne. ' ') exit
    enddo
    
    ! Get the rest, up to any comment
    sValue = sLine(iFrom:)
    iTo = len_trim(sValue)
    iFrom = min(index(sValue,'%'), index(sValue,'!'))
    if (iFrom > 0 .and. iFrom < iTo) then
        sValue(iFrom:iTo) = ' '
    endif
    !call Lower(sValue)     ! No: Some values are filenames which are case-sensitive on UNIX!
    
end subroutine ParseCode


!-----------------------------------------------------------------------
subroutine Lower( s )
!-----------------------------------------------------------------------
! OCCAM 3.0 Package
! David Myer IGPP/SIO La Jolla CA 92093-0225
! Subroutine Revision 3.0, November 2006
! DGM Nov 2006 - convert string to lower case
    implicit none
    character(*), intent(out)  :: s
    integer i
!    forall( i=1:len_trim(s), s(i:i) >= 'A' .and. s(i:i) <= 'Z' ) 
!        s(i:i) = char(ichar(s(i:i)) + 32)
!    end forall
    do  i=1,len_trim(s)
      if  ( s(i:i) >= 'A' .and. s(i:i) <= 'Z' ) then
        s(i:i) = char(ichar(s(i:i)) + 32)
        endif
    enddo
    
end subroutine Lower


!-----------------------------------------------------------------------
      real function rdcode(string, iof2)
!-----------------------------------------------------------------------
! OCCAM 3.0 Package
! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
! Subroutine Revision 3.0, November 2006
! Subroutine Revision 2.00, 13 Jan 1993
!
! rdcode provides a patch to get around the lack of an internal
! free-format (list-directed) read on some compilers, notably the
! cray.  This version for reals.  Butchered from decode by Bob Parker.
! 
! on input:
!    string = character string containing number to be read
!    iof2 = unit number of open file, to be closed on error 
! on output:
!    rdcode contains a real number from string
! calls:
!    filerr
!
! includes:
      IMPLICIT NONE
! input arguments:
      character(*) string
      integer iof2
! local variables
      character(40) local
      integer k, k1, kn, l
!
!  Terminate line at % sign
      kn=index(string//'%', '%') - 1
      k1=1
      do 1100 k=k1, kn
        if (string(k:k) .ne. ' ') goto 1200
 1100 continue
      call filerr (' No number found by rdcode', iof2)
 1200 do 1300 l=k, kn
        if (string(l:l).eq. ',' .or. string(l:l) .eq. ' ') goto 1500
 1300 continue
 1500 local=string(k:l-1)

      ! DGM - some Fortran compilers have trouble with CR/LF
      ! characters.  Remove that trouble now.
      k1 = index(local,char(13))  ! carriage return position
      if (k1.gt.0) local(k1:k1) = char(32)
      k1 = index(local,char(10))  ! line feed position
      if (k1.gt.0) local(k1:k1) = char(32)
 
      read (local, '(f40.0)', err=1900) rdcode
    return
 1900 call filerr (' Error reading string in rdcode', iof2)
end
      

!-----------------------------------------------------------------------
      integer function idcode(string, iof2)
!-----------------------------------------------------------------------
! OCCAM 3.0 Package
! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
! Subroutine Revision 3.0, November 2006
! Subroutine Revision 2.00, 13 Jan 1993
!
! idcode provides a patch to get around the lack of an internal
! free-format (list-directed) read on some compilers, notably the
! cray.  This version for integers.  Butchered from decode by Bob Parker.
! 
! on input:
!    string = character string containing number to be read
!    iof2 = unit number of open file, to be closed on error 
! on output:
!    idcode contains an integer number from string
! calls:
!    filerr
!
! includes:
      IMPLICIT NONE
! input arguments:
      character(*) string
      integer iof2
! DGM - code is same as rdcode, except output is made into real
      real rdcode
      
! local variables
!      character(40) local
!      integer k, k1, kn, l
!      real realno
!
!!  Terminate line at % sign
!      kn=index(string//'%', '%') - 1
!      k1=1
!      do 1100 k=k1, kn
!        if (string(k:k) .ne. ' ') goto 1200
! 1100 continue
!      call filerr (' No number found by idcode', iof2)
! 1200 do 1300 l=k, kn
!        if (string(l:l).eq. ',' .or. string(l:l) .eq. ' ') goto 1500
! 1300 continue
! 1500 local=string(k:l-1)
!
!      ! DGM - some Fortran compilers have trouble with CR/LF
!      ! characters.  Remove that trouble now.
!      k1 = index(local,char(13))  ! carriage return position
!      if (k1.gt.0) local(k1:k1) = char(32)
!      k1 = index(local,char(10))  ! line feed position
!      if (k1.gt.0) local(k1:k1) = char(32)
! 
!      read (local, '(f40.0)', err=1900) realno
!      idcode = int(realno)
        idcode = int( rdcode( string, iof2 ) )
    return
! 1900 call filerr (' Error reading string in idcode', iof2)
end


!-----------------------------------------------------------------------
SUBROUTINE indexx( n, arr, indx) 
!-----------------------------------------------------------------------
! OCCAM MT2D 3.0 Package
! Subroutine Revision 3.0, November 2006 (notional only to track future fixes)
! Original FORTRAN program from Numerical Recipes.
!
! Indexes an array arr(1:n), i.e., outputs the array indx(1:n) such that 
! arr(indx(j)) is in ascending order for j=1,2,...,N. 
! The input quantities n and arr are not changed.
      implicit none
      INTEGER n,indx(n),M,NSTACK
      REAL(4) arr(n)
      PARAMETER (M=7,NSTACK=1000)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL(4) a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,l,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=l-1
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(l+1)))then
          itemp=indx(l)
          indx(l)=indx(l+1)
          indx(l+1)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l+1)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l+1)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK)stop 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
END SUBROUTINE indexx

      
!-----------------------------------------------------------------------
      subroutine filerr(mssg, io1)
!-----------------------------------------------------------------------
! OCCAM 3.0 Package
! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
! Subroutine Revision 3.0, November 2006
! Subroutine Revision 2.00, 13 May 1992
!
! filerr prints an error message, closes file io1, and stops
!
! on input:
!    mssg = character string containing error message
!    i01 = unit number of file to close (0 = no file open)
! on output:
!    outputs nothing
! calls:
!    no other routines
!
!
! includes:
      IMPLICIT NONE
! input arguments:
      character(*) mssg
      integer io1
!
!-----------------------------------------------------------------------
      write(*,*) mssg
      if (io1 .gt. 0) close (io1)
      stop
end


!-----------------------------------------------------------------------
      subroutine datime(datetm)
!-----------------------------------------------------------------------
! OCCAM 3.0 Package
! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
! Subroutine Revision 3.0, November 2006
! Subroutine Revision 1.00, 18 May 1992
!
! Date and time utility for OCCAM 2.0.  
!
! on input:
!   nothing
! on output:
!   datetm = 80 character string containing message saying that date and time
!           are not available
!
      implicit none
      character(80) datetm
      integer icr
!      
!-----------------------------------------------------------------------
!
! Default is no time given (this makes compilation easier for novices)
!
      datetm = 'No date and time available on this system.'      
!
! To instead get the time, uncomment the command for your platform below: 
!
! Portland Group PGF77
!      call time(datetm)
! IBM XLF (Mac OS/X)
!      call fdate_(datetm)
! Intel Fortran 9.1
!      call fdate(datetm)
      
      
      ! DGM some fdate()s put CR and/or LF in the string.  Replace
      ! these with spaces for easy trimming.
      icr = index(datetm,char(13))  ! carriage return position
      if (icr.gt.0) datetm(icr:icr) = char(32)
      icr = index(datetm,char(10))  ! line feed position
      if (icr.gt.0) datetm(icr:icr) = char(32)
      return
end
