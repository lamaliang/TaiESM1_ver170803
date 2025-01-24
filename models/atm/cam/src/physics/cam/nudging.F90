module nudging
!=====================================================================
!
! Purpose: Implement Nudging of the model state of U,V,T,Q, and/or PS
!          toward specified values from analyses.
!
! Author: Patrick Callaghan
!
! Description:
!    This module assumes that the user has {U,V,T,Q,PS} analyses which 
!    have been preprocessed onto the current model grid and are stored 
!    in individual files which are indexed with respect to year, month, 
!    day, and second of the day. When the model is inbetween the given 
!    begining and ending times, forcing is added to nudge the model toward
!    the appropriate analyses values. After the model passes the ending 
!    analyses time, the forcing discontinues.
!
! Revisions:
!    01/14/13 - Modified to manage 'GAPS' in analyses data. For now the
!               approach is to coast through the gaps...  If a given
!               analyses file is missing, nudging is turned off for 
!               that interval of time. Once an analyses file is found, 
!               the Nudging is switched back on.
!    02/22/13 - Modified to add functionality for FV and EUL dynamical 
!               cores.
!    03/03/13 - For ne120 runs, the automatic arrays used for reading in
!               U,V,T,Q,PS values were putting too much of a burden on the
!               stack memory. Until Parallel I/O is implemented, the impact
!               on the stack was reduced by using only one automatic array
!               to read in and scatter the data.
!    04/01/13 - Added Heaviside window function for localized nudging
!    04/10/13 - Modified call to physics_ptend_init() to accomodate the
!               new interface (in CESM1_2_BETA05).
!    05/06/13 - 'WRAP_NF' was modified from a generic interface so that 
!               now it can only read in 1D arrays from netCDF files. 
!               To eliminate errors from future meddling of this sort, all 
!               refenences to the 'wrap_nf' module were removed and replaced 
!               with direct nf90 calls.
!    08/19/13 - Add optional forms for Nudging force.
!    10/16/13 - Add option for Nudging Diagnostic outputs.
!               Move application Nudging tendency from the end of tphysbc()
!               to the end of tphysac() [DONE IN PHYSPKG.F90]
!    11/11/13 - Remove the FV kludge to use staggered winds (US,VS) 
!               instead of (U,V) - for FV the input datasets are assumed 
!               to contain both (US,VS) and (U,V).
!    11/12/13 - Add diurnal filter forcing options.
!               ** Forcing options 1 and 2 have swapped from what they were
!               ** before this date.
!    11/27/13 - Add routine to calc Dry Static Energy and modify the
!               tendency values from temperature only to DSE values.
!               Added 'Nudge_TSmode' (internal only) to toggle between 
!               nudging DSE or temperature only.
!     4/01/14 - Fixed Radian->Degree error with windows (Pedro DiNezio)
!
! Input/Output Values:
!    Forcing contributions are available for history file output by 
!    the names:    {'Nudge_U','Nudge_V','Nudge_T',and 'Nudge_Q'}
!
!    The nudging of the model toward the analyses data is controlled by 
!    the 'nudging_nl' namelist in 'user_nl_cam'; whose variables control the
!    time interval over which nudging is applied, the strength of the nudging
!    tendencies, and its spatial distribution. The strength of the nudging is
!    specified as a fractional coeffcient between [0,1]. The spatial distribution 
!    is specified with a profile index:
!
!        (U,V,T,Q) Profiles:      0 == OFF      (No Nudging of this variable)
!        -------------------      1 == CONSTANT (Spatially Uniform Nudging)
!                                 2 == HEAVISIDE WINDOW FUNCTION
!
!        (PS) Profiles:           0 == OFF (Not Implemented)
!        -------------------      1 == N/A (Not Implemented)
!                  
!    The Heaviside window function is the product of separate horizonal and vertical 
!    windows that are controled via 14 parameters:
!        Nudge_Hwin_lat0:     Provide the horizontal center of the window in degrees. 
!        Nudge_Hwin_lon0:     The longitude must be in the range [0,360] and the 
!                             latitude should be [-90,+90].
!
!        Nudge_Hwin_latWidth: Specify the lat and lon widths of the window as positive 
!        Nudge_Hwin_lonWidth: values in degrees.Setting a width to a large value (e.g. 999) 
!                             renders the window a constant in that direction.
!                                                 
!        Nudge_Hwin_latDelta: Controls the sharpness of the window transition with a 
!        Nudge_Hwin_lonDelta: length in degrees. Small non-zero values yeild a step 
!                             function while a large value leads to a smoother transition.
!
!        Nudge_Vwin_Lindex:   In the vertical, the window is specified in terms of model 
!        Nudge_Vwin_Ldelta:   level indcies. The High and Low transition levels should 
!        Nudge_Vwin_Hindex:   range from [0,(NCOL+1)]. The transition lengths are also 
!        Nudge_Vwin_Hdelta:   specified in terms of model indices. For a window function 
!                             constant in the vertical, the Low index should be set to 0,
!                             the High index should be set to (NCOL+1), and the transition 
!                             lengths should be set to 0.1 
!
!        Nudge_Hwin_lo:       For a given set of spatial parameters, the raw window 
!        Nudge_Hwin_hi:       function may not span the range [0,1], so those values are 
!        Nudge_Vwin_lo:       mapped to the range of values specified in by the user. 
!        Nudge_Vwin_hi:       The 'hi' values are mapped to the maximum of the raw window 
!                             function and 'lo' values are mapped to its minimum. 
!                             Typically the 'hi' values will be set equal to 1, and the 
!                             'lo' values set equal 0 or the desired window minimum. 
!                             Specifying the 'lo' value as 1 and the 'hi' value as 0 acts 
!                             to invert the window function. For a properly specified
!                             window its maximum should be equal to 1: MAX('lo','hi')==1
!
!        EXAMPLE: For a channel window function centered at the equator and independent 
!                 of the vertical (30 levels):
!                        Nudge_Hwin_lo = 0.               Nudge_Vwin_lo = 0.
!                        Nudge_Hwin_hi = 1.               Nudge_Vwin_hi = 1.
!                        Nudge_Hwin_lat0     = 0.         Nudge_Vwin_Lindex = 0.
!                        Nudge_Hwin_latWidth = 30.        Nudge_Vwin_Ldelta = 0.1
!                        Nudge_Hwin_latDelta = 5.0        Nudge_Vwin_Hindex = 31.
!                        Nudge_Hwin_lon0     = 180.       Nudge_Vwin_Hdelta = 0.1
!                        Nudge_Hwin_lonWidth = 999.
!                        Nudge_Hwin_lonDelta = 1.0
!
!                 If on the other hand one desired to apply nudging at the poles and
!                 not at the equator, the settings would be similar but with:
!                        Nudge_Hwin_lo = 1.
!                        Nudge_Hwin_hi = 0.
!
!    &nudging_nl
!      Nudge_Model         - LOGICAL toggle to activate nudging.
!      Nudge_Path          - CHAR path to the analyses files.
!      Nudge_File_Template - CHAR Analyses filename with year, month, day, and second
!                                 values replaced by %y, %m, %d, and %s respectively.
!      Nudge_Force_Opt     - INT Index to select the nudging force with the form:
!
!                                       F_nudge = Alpha*((Target-Model(t_curr))/TimeScale
!           
!            (default)     0 -->  Target=Anal(t'_next)            , TimeScale=Tdlt_Anal
!                          1 -->  Target=Anal(t'_next)            , TimeScale=(t'_next - t_curr )
!                          2 -->  Target=Anal(t'_curr)            , TimeScale=Tdlt_Anal
!                          3 -->  Target=Anal(t'_curr)            , TimeScale=(t_curr  - t'_curr)
!                          4 -->  Target=Diurnal_Estimate(t'_next), TimeScale=Tdlt_Anal
!                          5 -->  Target=Diurnal_Estimate(t'_next), TimeScale=(t'_next - t_curr )
!                          6 -->  Target=*STABLE*_Diurnal(t'_next), TimeScale=Tdlt_Anal
!                          7 -->  Target=*STABLE*_Diurnal(t'_next), TimeScale=(t'_next - t_curr )
!
!                                where (t'==Analysis times ; t==Model Times) and Diurnal estimates 
!                                are calcualted using 1 cycle of previous values[Nudge_Times_Per_Day]. 
!
!      Nudge_Diag_Opt      - INT Index to select diagnostic output.
!            (default)     0 -->  No diagnostic outputs.
!                          1 -->  10 [U,V,T,Q] outputs in tphysbc().
!                          2 -->  10 [U,V,T,Q] outputs mostly in tphysac().
!                          3 -->  What do you want??
!
!      Nudge_Times_Per_Day - INT Number of analyses files available per day.
!      Model_Times_Per_Day - INT Number of times to update the model state (used for nudging) 
!                                each day. The value is restricted to be longer than the 
!                                current model timestep and shorter than the analyses 
!                                timestep. As this number is increased, the nudging
!                                force has the form of newtonian cooling.
!      Nudge_Uprof         - INT index of profile structure to use for U.  [0,1,2]
!      Nudge_Vprof         - INT index of profile structure to use for V.  [0,1,2]
!      Nudge_Tprof         - INT index of profile structure to use for T.  [0,1,2]
!      Nudge_Qprof         - INT index of profile structure to use for Q.  [0,1,2]
!      Nudge_PSprof        - INT index of profile structure to use for PS. [0,N/A]
!      Nudge_Ucoef         - REAL fractional nudging coeffcient for U. 
!                                    Utau=(Nudge_Ucoef/analyses_timestep)
!      Nudge_Vcoef         - REAL fractional nudging coeffcient for V. 
!                                    Vtau=(Nudge_Vcoef/analyses_timestep)
!      Nudge_Tcoef         - REAL fractional nudging coeffcient for T. 
!                                    Stau=(Nudge_Tcoef/analyses_timestep)
!      Nudge_Qcoef         - REAL fractional nudging coeffcient for Q. 
!                                    Qtau=(Nudge_Qcoef/analyses_timestep)
!      Nudge_PScoef        - REAL fractional nudging coeffcient for PS. 
!                                    PStau=(Nudge_PScoef/analyses_timestep)
!      Nudge_Beg_Year      - INT nudging begining year.
!      Nudge_Beg_Month     - INT nudging begining month.
!      Nudge_Beg_Day       - INT nudging begining day.
!      Nudge_End_Year      - INT nudging ending year.
!      Nudge_End_Month     - INT nudging ending month.
!      Nudge_End_Day       - INT nudging ending day.
!      Nudge_Hwin_lo       - REAL value mapped to RAW horizontal window minimum. [0]
!      Nudge_Hwin_hi       - REAL value mapped to RAW horizontal window maximum. [1]
!      Nudge_Vwin_lo       - REAL value mapped to RAW vertical window minimum.   [0]
!      Nudge_Vwin_hi       - REAL value mapped to RAW vertical window maximum.   [1]
!      Nudge_Hwin_lat0     - REAL latitudinal center of window in degrees.
!      Nudge_Hwin_lon0     - REAL longitudinal center of window in degrees.
!      Nudge_Hwin_latWidth - REAL latitudinal width of window in degrees.
!      Nudge_Hwin_lonWidth - REAL longitudinal width of window in degrees.
!      Nudge_Hwin_latDelta - REAL latitudinal transition length of window in degrees.
!      Nudge_Hwin_lonDelta - REAL longitudinal transition length of window in degrees.
!      Nudge_Vwin_Lindex   - REAL LO model index of transition
!      Nudge_Vwin_Hindex   - REAL HI model index of transition
!      Nudge_Vwin_Ldelta   - REAL LO transition length 
!      Nudge_Vwin_Hdelta   - REAL HI transition length 
!    /
!
!================
!
! TO DO:
! -----------
!    ** Currently the surface pressure is read in, but there is no forcing
!       meachnism implemented.
!    ** Analyses data is read in and then distributed to processing elements 
!       via 'scatted_field_to_chunk' calls. The SE's want this to be changed
!       to parallel I/O calls.
!    ** Possibly implement time variation to nudging coeffcients, so that 
!       rather than just bashing the model with a sledge hammer, the user has the
!       option to ramp up the nudging coefs over a startup time frame via a 
!       heavyside step function.
!          
!=====================================================================
  ! Useful modules
  !------------------
  use shr_kind_mod,only:r8=>SHR_KIND_R8,cs=>SHR_KIND_CS,cl=>SHR_KIND_CL
  use time_manager,only:timemgr_time_ge,timemgr_time_inc,get_curr_date,dtime
  use phys_grid   ,only:scatter_field_to_chunk
  use abortutils  ,only:endrun
  use spmd_utils  ,only:masterproc
  use cam_logfile ,only:iulog
#ifdef SPMD
  use mpishorthand
#endif

  ! Set all Global values and routines to private by default 
  ! and then explicitly set their exposure.
  !----------------------------------------------------------
  implicit none
  private

  public:: Nudge_Model,Nudge_ON,Nudge_Diag_Opt
  public:: nudging_readnl
  public:: nudging_init
  public:: nudging_timestep_init
  public:: nudging_timestep_tend
  public:: nudging_diag_init
  public:: nudging_diag
  private::nudging_update_analyses_se
  private::nudging_update_analyses_eul
  private::nudging_update_analyses_fv
  private::nudging_set_PSprofile
  private::nudging_set_profile
  private::calc_DryStaticEnergy

  ! Nudging Parameters
  !--------------------
  logical::         Nudge_Model       =.false.
  logical::         Nudge_ON          =.false.
  logical::         Nudge_Initialized =.false.
  character(len=cl) Nudge_Path
  character(len=cs) Nudge_File,Nudge_File_Template
  integer           Nudge_Force_Opt
  integer           Nudge_Diag_Opt
  integer           Nudge_TSmode
  integer           Nudge_Times_Per_Day
  integer           Model_Times_Per_Day
  real(r8)          Nudge_Ucoef,Nudge_Vcoef
  integer           Nudge_Uprof,Nudge_Vprof
  real(r8)          Nudge_Qcoef,Nudge_Tcoef
  integer           Nudge_Qprof,Nudge_Tprof
  real(r8)          Nudge_PScoef
  integer           Nudge_PSprof
  integer           Nudge_Beg_Year ,Nudge_Beg_Month
  integer           Nudge_Beg_Day  ,Nudge_Beg_Sec
  integer           Nudge_End_Year ,Nudge_End_Month
  integer           Nudge_End_Day  ,Nudge_End_Sec
  integer           Nudge_Curr_Year,Nudge_Curr_Month
  integer           Nudge_Curr_Day ,Nudge_Curr_Sec
  integer           Nudge_Next_Year,Nudge_Next_Month
  integer           Nudge_Next_Day ,Nudge_Next_Sec
  integer           Nudge_Step
  integer           Model_Curr_Year,Model_Curr_Month
  integer           Model_Curr_Day ,Model_Curr_Sec
  integer           Model_Next_Year,Model_Next_Month
  integer           Model_Next_Day ,Model_Next_Sec
  integer           Model_Step
  real(r8)          Nudge_Hwin_lo
  real(r8)          Nudge_Hwin_hi
  real(r8)          Nudge_Hwin_lat0
  real(r8)          Nudge_Hwin_latWidth
  real(r8)          Nudge_Hwin_latDelta
  real(r8)          Nudge_Hwin_lon0
  real(r8)          Nudge_Hwin_lonWidth
  real(r8)          Nudge_Hwin_lonDelta
  real(r8)          Nudge_Vwin_lo
  real(r8)          Nudge_Vwin_hi
  real(r8)          Nudge_Vwin_Hindex
  real(r8)          Nudge_Vwin_Hdelta
  real(r8)          Nudge_Vwin_Lindex
  real(r8)          Nudge_Vwin_Ldelta
  real(r8)          Nudge_Hwin_latWidthH
  real(r8)          Nudge_Hwin_lonWidthH
  real(r8)          Nudge_Hwin_max
  real(r8)          Nudge_Hwin_min

  ! Nudging State Arrays
  !-----------------------
  integer Nudge_nlon,Nudge_nlat,Nudge_ncol,Nudge_nlev
  real(r8),allocatable::Target_U(:,:,:)     !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_V(:,:,:)     !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_T(:,:,:)     !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_S(:,:,:)     !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_Q(:,:,:)     !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_PS(:,:)      !(pcols,begchunk:endchunk)
  real(r8),allocatable::Model_U(:,:,:)      !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Model_V(:,:,:)      !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Model_T(:,:,:)      !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Model_S(:,:,:)      !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Model_Q(:,:,:)      !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Model_PS(:,:)       !(pcols,begchunk:endchunk)
  real(r8),allocatable::Nudge_Utau(:,:,:)   !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Nudge_Vtau(:,:,:)   !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Nudge_Stau(:,:,:)   !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Nudge_Qtau(:,:,:)   !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Nudge_PStau(:,:)    !(pcols,begchunk:endchunk)
  real(r8),allocatable::Nudge_Ustep(:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Nudge_Vstep(:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Nudge_Sstep(:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Nudge_Qstep(:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Nudge_PSstep(:,:)   !(pcols,begchunk:endchunk)

  ! Nudging Observation Arrays
  !-----------------------------
  integer               Nudge_NumObs
  integer,allocatable:: Nudge_ObsInd(:)
  logical ,allocatable::Nudge_File_Present(:)
  real(r8)              Nudge_Acoef
  real(r8),allocatable::Nudge_Bcoef(:)
  real(r8),allocatable::Nudge_Ccoef(:)
  real(r8),allocatable::Nobs_U(:,:,:,:) !(pcols,pver,begchunk:endchunk,Nudge_NumObs)
  real(r8),allocatable::Nobs_V(:,:,:,:) !(pcols,pver,begchunk:endchunk,Nudge_NumObs)
  real(r8),allocatable::Nobs_T(:,:,:,:) !(pcols,pver,begchunk:endchunk,Nudge_NumObs)
  real(r8),allocatable::Nobs_Q(:,:,:,:) !(pcols,pver,begchunk:endchunk,Nudge_NumObs)
  real(r8),allocatable::Nobs_PS(:,:,:)  !(pcols,begchunk:endchunk,Nudge_NumObs)
  real(r8),allocatable::Mobs_U(:,:,:,:) !(pcols,pver,begchunk:endchunk,Nudge_NumObs)
  real(r8),allocatable::Mobs_V(:,:,:,:) !(pcols,pver,begchunk:endchunk,Nudge_NumObs)
  real(r8),allocatable::Mobs_T(:,:,:,:) !(pcols,pver,begchunk:endchunk,Nudge_NumObs)
  real(r8),allocatable::Mobs_Q(:,:,:,:) !(pcols,pver,begchunk:endchunk,Nudge_NumObs)
  real(r8),allocatable::Mobs_PS(:,:,:)  !(pcols,begchunk:endchunk,Nudge_NumObs)

contains
  !================================================================
  subroutine nudging_readnl(nlfile)
   ! 
   ! NUDGING_READNL: Initialize default values controlling the Nudging 
   !                 process. Then read namelist values to override 
   !                 them.
   !===============================================================
   use ppgrid        ,only: pver
   use namelist_utils,only:find_group_name
   use units         ,only:getunit,freeunit
   !
   ! Arguments
   !-------------
   character(len=*),intent(in)::nlfile
   !
   ! Local Values
   !---------------
   integer ierr,unitn

   namelist /nudging_nl/ Nudge_Model,Nudge_Path,                       &
                         Nudge_File_Template,Nudge_Force_Opt,          &
                         Nudge_Diag_Opt,                               &
                         Nudge_Times_Per_Day,Model_Times_Per_Day,      &
                         Nudge_Ucoef ,Nudge_Uprof,                     &
                         Nudge_Vcoef ,Nudge_Vprof,                     &
                         Nudge_Qcoef ,Nudge_Qprof,                     &
                         Nudge_Tcoef ,Nudge_Tprof,                     &
                         Nudge_PScoef,Nudge_PSprof,                    &
                         Nudge_Beg_Year,Nudge_Beg_Month,Nudge_Beg_Day, &
                         Nudge_End_Year,Nudge_End_Month,Nudge_End_Day, &
                         Nudge_Hwin_lo,Nudge_Hwin_hi,                  &
                         Nudge_Vwin_lo,Nudge_Vwin_hi,                  &
                         Nudge_Hwin_lat0,Nudge_Hwin_lon0,              &
                         Nudge_Hwin_latWidth,Nudge_Hwin_lonWidth,      &
                         Nudge_Hwin_latDelta,Nudge_Hwin_lonDelta,      &
                         Nudge_Vwin_Lindex,Nudge_Vwin_Hindex,          &
                         Nudge_Vwin_Ldelta,Nudge_Vwin_Hdelta           

   ! Nudging is NOT initialized yet, For now
   ! Nudging will always begin/end at midnight.
   !--------------------------------------------
   Nudge_Initialized =.false.
   Nudge_ON          =.false.
   Nudge_Beg_Sec=0
   Nudge_End_Sec=0

   ! Set Default Namelist values
   !-----------------------------
   Nudge_Model         =.false.
   Nudge_Path          ='./Data/YOTC_ne30np4_001/'
   Nudge_File_Template ='YOTC_ne30np4_L30.cam2.i.%y-%m-%d-%s.nc'
   Nudge_Force_Opt=0
   Nudge_Diag_Opt =0
   Nudge_TSmode   =0
   Nudge_Times_Per_Day=4
   Model_Times_Per_Day=4
   Nudge_Ucoef  =0._r8
   Nudge_Vcoef  =0._r8
   Nudge_Qcoef  =0._r8
   Nudge_Tcoef  =0._r8
   Nudge_PScoef =0._r8
   Nudge_Uprof  =0
   Nudge_Vprof  =0
   Nudge_Qprof  =0
   Nudge_Tprof  =0
   Nudge_PSprof =0
   Nudge_Beg_Year =2008
   Nudge_Beg_Month=5
   Nudge_Beg_Day  =1
   Nudge_End_Year =2008
   Nudge_End_Month=9
   Nudge_End_Day  =1
   Nudge_Hwin_lo      =0.0_r8
   Nudge_Hwin_hi      =1.0_r8
   Nudge_Hwin_lat0    =0._r8
   Nudge_Hwin_latWidth=9999._r8
   Nudge_Hwin_latDelta=1.0_r8
   Nudge_Hwin_lon0    =180._r8
   Nudge_Hwin_lonWidth=9999._r8
   Nudge_Hwin_lonDelta=1.0_r8
   Nudge_Vwin_lo      =0.0_r8
   Nudge_Vwin_hi      =1.0_r8
   Nudge_Vwin_Hindex  =float(pver+1)
   Nudge_Vwin_Hdelta  =0.1_r8
   Nudge_Vwin_Lindex  =0.0_r8
   Nudge_Vwin_Ldelta  =0.1_r8

   ! Read in namelist values
   !------------------------
   if(masterproc) then
     unitn = getunit()
     open(unitn,file=trim(nlfile),status='old')
     call find_group_name(unitn,'nudging_nl',status=ierr)
     if(ierr.eq.0) then
       read(unitn,nudging_nl,iostat=ierr)
       if(ierr.ne.0) then
         call endrun('nudging_readnl:: ERROR reading namelist')
       endif
     endif
     close(unitn)
     call freeunit(unitn)
   endif

   ! Check for valid namelist values 
   !----------------------------------
   if((max(Nudge_Hwin_lo,Nudge_Hwin_hi).ne.1.0).or. &
      (max(Nudge_Vwin_lo,Nudge_Vwin_hi).ne.1.0)   ) then
     write(iulog,*) 'NUDGING: The window function must have a maximum value of 1'
     write(iulog,*) 'NUDGING:  Nudge_Hwin_lo=',Nudge_Hwin_lo
     write(iulog,*) 'NUDGING:  Nudge_Hwin_hi=',Nudge_Hwin_hi
     write(iulog,*) 'NUDGING:  Nudge_Vwin_lo=',Nudge_Vwin_lo
     write(iulog,*) 'NUDGING:  Nudge_Vwin_hi=',Nudge_Vwin_hi
     call endrun('nudging_readnl:: ERROR in namelist')
   endif

   if((Nudge_Hwin_lat0.lt.-90.).or.(Nudge_Hwin_lat0.gt.+90.)) then
     write(iulog,*) 'NUDGING: Window lat0 must be in [-90,+90]'
     write(iulog,*) 'NUDGING:  Nudge_Hwin_lat0=',Nudge_Hwin_lat0
     call endrun('nudging_readnl:: ERROR in namelist')
   endif

   if((Nudge_Hwin_lon0.lt.0.).or.(Nudge_Hwin_lon0.ge.360.)) then
     write(iulog,*) 'NUDGING: Window lon0 must be in [0,+360)'
     write(iulog,*) 'NUDGING:  Nudge_Hwin_lon0=',Nudge_Hwin_lon0
     call endrun('nudging_readnl:: ERROR in namelist')
   endif

   if((Nudge_Vwin_Lindex.gt.Nudge_Vwin_Hindex)                         .or. &
      (Nudge_Vwin_Hindex.gt.float(pver+1)).or.(Nudge_Vwin_Hindex.lt.0.).or. &
      (Nudge_Vwin_Lindex.gt.float(pver+1)).or.(Nudge_Vwin_Lindex.lt.0.)   ) then
     write(iulog,*) 'NUDGING: Window Lindex must be in [0,pver+1]'
     write(iulog,*) 'NUDGING: Window Hindex must be in [0,pver+1]'
     write(iulog,*) 'NUDGING: Lindex must be LE than Hindex'
     write(iulog,*) 'NUDGING:  Nudge_Vwin_Lindex=',Nudge_Vwin_Lindex
     write(iulog,*) 'NUDGING:  Nudge_Vwin_Hindex=',Nudge_Vwin_Hindex
     call endrun('nudging_readnl:: ERROR in namelist')
   endif

   if((Nudge_Hwin_latDelta.le.0.).or.(Nudge_Hwin_lonDelta.le.0.).or. &
      (Nudge_Vwin_Hdelta  .le.0.).or.(Nudge_Vwin_Ldelta  .le.0.)    ) then
     write(iulog,*) 'NUDGING: Window Deltas must be positive'
     write(iulog,*) 'NUDGING:  Nudge_Hwin_latDelta=',Nudge_Hwin_latDelta
     write(iulog,*) 'NUDGING:  Nudge_Hwin_lonDelta=',Nudge_Hwin_lonDelta
     write(iulog,*) 'NUDGING:  Nudge_Vwin_Hdelta=',Nudge_Vwin_Hdelta
     write(iulog,*) 'NUDGING:  Nudge_Vwin_Ldelta=',Nudge_Vwin_Ldelta
     call endrun('nudging_readnl:: ERROR in namelist')

   endif

   if((Nudge_Hwin_latWidth.le.0.).or.(Nudge_Hwin_lonWidth.le.0.)) then
     write(iulog,*) 'NUDGING: Window widths must be positive'
     write(iulog,*) 'NUDGING:  Nudge_Hwin_latWidth=',Nudge_Hwin_latWidth
     write(iulog,*) 'NUDGING:  Nudge_Hwin_lonWidth=',Nudge_Hwin_lonWidth
     call endrun('nudging_readnl:: ERROR in namelist')
   endif

   ! Broadcast namelist variables
   !------------------------------
#ifdef SPMD
   call mpibcast(Nudge_Path         ,len(Nudge_Path)         ,mpichar,0,mpicom)
   call mpibcast(Nudge_File_Template,len(Nudge_File_Template),mpichar,0,mpicom)
   call mpibcast(Nudge_Model        , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_Initialized  , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_ON           , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_Force_Opt    , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Diag_Opt     , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_TSmode       , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Times_Per_Day, 1, mpiint, 0, mpicom)
   call mpibcast(Model_Times_Per_Day, 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Ucoef        , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Vcoef        , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Tcoef        , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Qcoef        , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_PScoef       , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Uprof        , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Vprof        , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Tprof        , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Qprof        , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_PSprof       , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Beg_Year     , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Beg_Month    , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Beg_Day      , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Beg_Sec      , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_End_Year     , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_End_Month    , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_End_Day      , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_End_Sec      , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Hwin_lo      , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_hi      , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_lat0    , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_latWidth, 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_latDelta, 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_lon0    , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_lonWidth, 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_lonDelta, 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Vwin_lo      , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Vwin_hi      , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Vwin_Hindex  , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Vwin_Hdelta  , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Vwin_Lindex  , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Vwin_Ldelta  , 1, mpir8 , 0, mpicom)
#endif

   ! End Routine
   !------------
   return
  end subroutine ! nudging_readnl
  !================================================================


  !================================================================
  subroutine nudging_diag_init
   ! 
   ! NUDGING_DIAG_INIT: Register diagnostic outputs fo U,V,T,Q
   !                    so values can be incrementally sampled
   !                    as physics paramterizations are added
   !===============================================================
   use ppgrid        ,only: pver
   use cam_history   ,only: addfld,phys_decomp

   ! Register Diagnostic output fields at 10 points
   !-----------------------------------------------------
   call addfld('UAP0','m/s'  ,pver,'A','U AfterPhysics Diag 0',phys_decomp)
   call addfld('VAP0','m/s'  ,pver,'A','V AfterPhysics Diag 0',phys_decomp)
   call addfld('TAP0','K'    ,pver,'A','T AfterPhysics Diag 0',phys_decomp)
   call addfld('SAP0','J'    ,pver,'A','S AfterPhysics Diag 0',phys_decomp)
   call addfld('QAP0','kg/kg',pver,'A','Q AfterPhysics Diag 0',phys_decomp)

   call addfld('UAP1','m/s'  ,pver,'A','U AfterPhysics Diag 1',phys_decomp)
   call addfld('VAP1','m/s'  ,pver,'A','V AfterPhysics Diag 1',phys_decomp)
   call addfld('TAP1','K'    ,pver,'A','T AfterPhysics Diag 1',phys_decomp)
   call addfld('SAP1','J'    ,pver,'A','S AfterPhysics Diag 1',phys_decomp)
   call addfld('QAP1','kg/kg',pver,'A','Q AfterPhysics Diag 1',phys_decomp)

   call addfld('UAP2','m/s'  ,pver,'A','U AfterPhysics Diag 2',phys_decomp)
   call addfld('VAP2','m/s'  ,pver,'A','V AfterPhysics Diag 2',phys_decomp)
   call addfld('TAP2','K'    ,pver,'A','T AfterPhysics Diag 2',phys_decomp)
   call addfld('SAP2','J'    ,pver,'A','S AfterPhysics Diag 2',phys_decomp)
   call addfld('QAP2','kg/kg',pver,'A','Q AfterPhysics Diag 2',phys_decomp)

   call addfld('UAP3','m/s'  ,pver,'A','U AfterPhysics Diag 3',phys_decomp)
   call addfld('VAP3','m/s'  ,pver,'A','V AfterPhysics Diag 3',phys_decomp)
   call addfld('TAP3','K'    ,pver,'A','T AfterPhysics Diag 3',phys_decomp)
   call addfld('SAP3','J'    ,pver,'A','S AfterPhysics Diag 3',phys_decomp)
   call addfld('QAP3','kg/kg',pver,'A','Q AfterPhysics Diag 3',phys_decomp)

   call addfld('UAP4','m/s'  ,pver,'A','U AfterPhysics Diag 4',phys_decomp)
   call addfld('VAP4','m/s'  ,pver,'A','V AfterPhysics Diag 4',phys_decomp)
   call addfld('TAP4','K'    ,pver,'A','T AfterPhysics Diag 4',phys_decomp)
   call addfld('SAP4','J'    ,pver,'A','S AfterPhysics Diag 4',phys_decomp)
   call addfld('QAP4','kg/kg',pver,'A','Q AfterPhysics Diag 4',phys_decomp)

   call addfld('UAP5','m/s'  ,pver,'A','U AfterPhysics Diag 5',phys_decomp)
   call addfld('VAP5','m/s'  ,pver,'A','V AfterPhysics Diag 5',phys_decomp)
   call addfld('TAP5','K'    ,pver,'A','T AfterPhysics Diag 5',phys_decomp)
   call addfld('SAP5','J'    ,pver,'A','S AfterPhysics Diag 5',phys_decomp)
   call addfld('QAP5','kg/kg',pver,'A','Q AfterPhysics Diag 5',phys_decomp)
   call addfld('UAP6','m/s'  ,pver,'A','U AfterPhysics Diag 6',phys_decomp)
   call addfld('VAP6','m/s'  ,pver,'A','V AfterPhysics Diag 6',phys_decomp)
   call addfld('TAP6','K'    ,pver,'A','T AfterPhysics Diag 6',phys_decomp)
   call addfld('SAP6','J'    ,pver,'A','S AfterPhysics Diag 6',phys_decomp)
   call addfld('QAP6','kg/kg',pver,'A','Q AfterPhysics Diag 6',phys_decomp)

   call addfld('UAP7','m/s'  ,pver,'A','U AfterPhysics Diag 7',phys_decomp)
   call addfld('VAP7','m/s'  ,pver,'A','V AfterPhysics Diag 7',phys_decomp)
   call addfld('TAP7','K'    ,pver,'A','T AfterPhysics Diag 7',phys_decomp)
   call addfld('SAP7','J'    ,pver,'A','S AfterPhysics Diag 7',phys_decomp)
   call addfld('QAP7','kg/kg',pver,'A','Q AfterPhysics Diag 7',phys_decomp)

   call addfld('UAP8','m/s'  ,pver,'A','U AfterPhysics Diag 8',phys_decomp)
   call addfld('VAP8','m/s'  ,pver,'A','V AfterPhysics Diag 8',phys_decomp)
   call addfld('TAP8','K'    ,pver,'A','T AfterPhysics Diag 8',phys_decomp)
   call addfld('SAP8','J'    ,pver,'A','S AfterPhysics Diag 8',phys_decomp)
   call addfld('QAP8','kg/kg',pver,'A','Q AfterPhysics Diag 8',phys_decomp)

   call addfld('UAP9','m/s'  ,pver,'A','U AfterPhysics Diag 9',phys_decomp)
   call addfld('VAP9','m/s'  ,pver,'A','V AfterPhysics Diag 9',phys_decomp)
   call addfld('TAP9','K'    ,pver,'A','T AfterPhysics Diag 9',phys_decomp)
   call addfld('SAP9','J'    ,pver,'A','S AfterPhysics Diag 9',phys_decomp)
   call addfld('QAP9','kg/kg',pver,'A','Q AfterPhysics Diag 9',phys_decomp)

   ! End Routine
   !------------
   return
  end subroutine ! nudging_diag_init
  !================================================================


  !================================================================
  subroutine nudging_diag(state,indx)
   ! NUDGING_DIAG: Write out the current state of U,V,T,Q to the
   !               specified increment index.
   !===============================================================
   use physics_types, only: physics_state
   use ppgrid,        only: pcols
   use cam_history   ,only: outfld

   ! Arguments
   !-----------
   type(physics_state),intent(in):: state
   integer            ,intent(in):: indx

   ! Local values
   !----------------
   integer  lchnk

   ! Write out the current state to the variables of the given index
   !------------------------------------------------------------------
   lchnk = state%lchnk

   if(indx.eq.0) then
     call outfld('UAP0',state%u       ,pcols,lchnk)
     call outfld('VAP0',state%v       ,pcols,lchnk)
     call outfld('TAP0',state%t       ,pcols,lchnk)
     call outfld('SAP0',state%s       ,pcols,lchnk)
     call outfld('QAP0',state%q(1,1,1),pcols,lchnk)
   elseif(indx.eq.1) then
     call outfld('UAP1',state%u       ,pcols,lchnk)
     call outfld('VAP1',state%v       ,pcols,lchnk)
     call outfld('TAP1',state%t       ,pcols,lchnk)
     call outfld('SAP1',state%s       ,pcols,lchnk)
     call outfld('QAP1',state%q(1,1,1),pcols,lchnk)
   elseif(indx.eq.2) then
     call outfld('UAP2',state%u       ,pcols,lchnk)
     call outfld('VAP2',state%v       ,pcols,lchnk)
     call outfld('TAP2',state%t       ,pcols,lchnk)
     call outfld('SAP2',state%s       ,pcols,lchnk)
     call outfld('QAP2',state%q(1,1,1),pcols,lchnk)
   elseif(indx.eq.3) then
     call outfld('UAP3',state%u       ,pcols,lchnk)
     call outfld('VAP3',state%v       ,pcols,lchnk)
     call outfld('TAP3',state%t       ,pcols,lchnk)
     call outfld('SAP3',state%s       ,pcols,lchnk)
     call outfld('QAP3',state%q(1,1,1),pcols,lchnk)
   elseif(indx.eq.4) then
     call outfld('UAP4',state%u       ,pcols,lchnk)
     call outfld('VAP4',state%v       ,pcols,lchnk)
     call outfld('TAP4',state%t       ,pcols,lchnk)
     call outfld('SAP4',state%s       ,pcols,lchnk)
     call outfld('QAP4',state%q(1,1,1),pcols,lchnk)
   elseif(indx.eq.5) then
     call outfld('UAP5',state%u       ,pcols,lchnk)
     call outfld('VAP5',state%v       ,pcols,lchnk)
     call outfld('TAP5',state%t       ,pcols,lchnk)
     call outfld('SAP5',state%s       ,pcols,lchnk)
     call outfld('QAP5',state%q(1,1,1),pcols,lchnk)
   elseif(indx.eq.6) then
     call outfld('UAP6',state%u       ,pcols,lchnk)
     call outfld('VAP6',state%v       ,pcols,lchnk)
     call outfld('TAP6',state%t       ,pcols,lchnk)
     call outfld('SAP6',state%s       ,pcols,lchnk)
     call outfld('QAP6',state%q(1,1,1),pcols,lchnk)
   elseif(indx.eq.7) then
     call outfld('UAP7',state%u       ,pcols,lchnk)
     call outfld('VAP7',state%v       ,pcols,lchnk)
     call outfld('TAP7',state%t       ,pcols,lchnk)
     call outfld('SAP7',state%s       ,pcols,lchnk)
     call outfld('QAP7',state%q(1,1,1),pcols,lchnk)
   elseif(indx.eq.8) then
     call outfld('UAP8',state%u       ,pcols,lchnk)
     call outfld('VAP8',state%v       ,pcols,lchnk)
     call outfld('TAP8',state%t       ,pcols,lchnk)
     call outfld('SAP8',state%s       ,pcols,lchnk)
     call outfld('QAP8',state%q(1,1,1),pcols,lchnk)
   elseif(indx.eq.9) then
     call outfld('UAP9',state%u       ,pcols,lchnk)
     call outfld('VAP9',state%v       ,pcols,lchnk)
     call outfld('TAP9',state%t       ,pcols,lchnk)
     call outfld('SAP9',state%s       ,pcols,lchnk)
     call outfld('QAP9',state%q(1,1,1),pcols,lchnk)
   else
     write(iulog,*) 'ERROR: nudging_diag(): indx=',indx
     call endrun('NUDGING: Unknown index for nudging_diag()')
   endif

   ! End Routine
   !------------
   return
  end subroutine ! nudging_diag
  !================================================================


  !================================================================
  subroutine nudging_init
   ! 
   ! NUDGING_INIT: Allocate space and initialize Nudging values
   !===============================================================
   use ppgrid        ,only: pver,pcols,begchunk,endchunk
   use error_messages,only: alloc_err
   use dycore        ,only: dycore_is
   use dyn_grid      ,only: get_horiz_grid_dim_d
   use phys_grid     ,only: get_rlat_p,get_rlon_p,get_ncols_p
   use cam_history   ,only: addfld,phys_decomp
   use shr_const_mod ,only: SHR_CONST_PI

   ! Local values
   !----------------
   integer  Year,Month,Day,Sec
   integer  YMD1,YMD
   logical  After_Beg,Before_End
   integer  istat,lchnk,ncol,icol,ilev
   integer  hdim1_d,hdim2_d
   real(r8) rlat,rlon
   real(r8) Wprof(pver)
   real(r8) lonp,lon0,lonn,latp,lat0,latn
   real(r8) Val1_p,Val2_p,Val3_p,Val4_p
   real(r8) Val1_0,Val2_0,Val3_0,Val4_0
   real(r8) Val1_n,Val2_n,Val3_n,Val4_n
   integer               nn
   real(r8)              NumObs,Freq
   real(r8),allocatable::CosVal(:)
   real(r8),allocatable::SinVal(:)

   ! Allocate Space for Nudging data arrays
   !-----------------------------------------
   allocate(Target_U(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Target_U',pcols*pver*((endchunk-begchunk)+1))
   allocate(Target_V(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Target_V',pcols*pver*((endchunk-begchunk)+1))
   allocate(Target_T(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Target_T',pcols*pver*((endchunk-begchunk)+1))
   allocate(Target_S(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Target_S',pcols*pver*((endchunk-begchunk)+1))
   allocate(Target_Q(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Target_Q',pcols*pver*((endchunk-begchunk)+1))
   allocate(Target_PS(pcols,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Target_PS',pcols*((endchunk-begchunk)+1))

   allocate(Model_U(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Model_U',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_V(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Model_V',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_T(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Model_T',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_S(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Model_S',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_Q(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Model_Q',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_PS(pcols,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Model_PS',pcols*((endchunk-begchunk)+1))

   ! Allocate Space for spatial dependence of 
   ! Nudging Coefs and Nudging Forcing.
   !-------------------------------------------
   allocate(Nudge_Utau(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Nudge_Utau',pcols*pver*((endchunk-begchunk)+1))
   allocate(Nudge_Vtau(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Nudge_Vtau',pcols*pver*((endchunk-begchunk)+1))
   allocate(Nudge_Stau(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Nudge_Stau',pcols*pver*((endchunk-begchunk)+1))
   allocate(Nudge_Qtau(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Nudge_Qtau',pcols*pver*((endchunk-begchunk)+1))
   allocate(Nudge_PStau(pcols,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Nudge_PStau',pcols*((endchunk-begchunk)+1))

   allocate(Nudge_Ustep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Nudge_Ustep',pcols*pver*((endchunk-begchunk)+1))
   allocate(Nudge_Vstep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Nudge_Vstep',pcols*pver*((endchunk-begchunk)+1))
   allocate(Nudge_Sstep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Nudge_Sstep',pcols*pver*((endchunk-begchunk)+1))
   allocate(Nudge_Qstep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Nudge_Qstep',pcols*pver*((endchunk-begchunk)+1))
   allocate(Nudge_PSstep(pcols,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Nudge_PSstep',pcols*((endchunk-begchunk)+1))

   ! Register output fields with the cam history module
   !-----------------------------------------------------
   call addfld('Nudge_U','m/s/s'  ,pver,'A','U Nudging Tendency',phys_decomp)
   call addfld('Nudge_V','m/s/s'  ,pver,'A','V Nudging Tendency',phys_decomp)
   call addfld('Nudge_T','K/s'    ,pver,'A','T Nudging Tendency',phys_decomp)
   call addfld('Nudge_Q','kg/kg/s',pver,'A','Q Nudging Tendency',phys_decomp)

   ! Add diagnistic output fileds
   !-------------------------------
   if(Nudge_Diag_Opt.ne.0) then
     call nudging_diag_init
   endif
   
   !-----------------------------------------
   ! Values initialized only by masterproc
   !-----------------------------------------
   if(masterproc) then

     ! Set the Stepping intervals for Model and Nudging values
     ! Ensure that the Model_Step is not smaller then one timestep
     !  and not larger then the Nudge_Step.
     !--------------------------------------------------------
     Model_Step=86400/Model_Times_Per_Day
     Nudge_Step=86400/Nudge_Times_Per_Day
     if(Model_Step.lt.dtime) then
       write(iulog,*) ' '
       write(iulog,*) 'NUDGING: Model_Step cannot be less than a model timestep'
       write(iulog,*) 'NUDGING:  Setting Model_Step=dtime , dtime=',dtime
       write(iulog,*) ' '
       Model_Step=dtime
     endif
     if(Model_Step.gt.Nudge_Step) then
       write(iulog,*) ' '
       write(iulog,*) 'NUDGING: Model_Step cannot be more than Nudge_Step'
       write(iulog,*) 'NUDGING:  Setting Model_Step=Nudge_Step, Nudge_Step=',Nudge_Step
       write(iulog,*) ' '
       Model_Step=Nudge_Step
     endif

     ! Initialize column and level dimensions
     !--------------------------------------------------------
     call get_horiz_grid_dim_d(hdim1_d,hdim2_d)
     Nudge_nlon=hdim1_d
     Nudge_nlat=hdim2_d
     Nudge_ncol=hdim1_d*hdim2_d
     Nudge_nlev=pver

     ! Check the time relative to the nudging window
     !------------------------------------------------
     call get_curr_date(Year,Month,Day,Sec)
     YMD=(Year*10000) + (Month*100) + Day
     YMD1=(Nudge_Beg_Year*10000) + (Nudge_Beg_Month*100) + Nudge_Beg_Day
     call timemgr_time_ge(YMD1,Nudge_Beg_Sec,         &
                          YMD ,Sec          ,After_Beg)
     YMD1=(Nudge_End_Year*10000) + (Nudge_End_Month*100) + Nudge_End_Day
     call timemgr_time_ge(YMD ,Sec          ,          &
                          YMD1,Nudge_End_Sec,Before_End)
  
     if((After_Beg).and.(Before_End)) then
       ! Set Time indicies so that the next call to 
       ! timestep_init will initialize the data arrays.
       !--------------------------------------------
       Model_Next_Year =Year
       Model_Next_Month=Month
       Model_Next_Day  =Day
       Model_Next_Sec  =(Sec/Model_Step)*Model_Step
       Nudge_Next_Year =Year
       Nudge_Next_Month=Month
       Nudge_Next_Day  =Day
       Nudge_Next_Sec  =(Sec/Nudge_Step)*Nudge_Step
     elseif(.not.After_Beg) then
       ! Set Time indicies to Nudging start,
       ! timestep_init will initialize the data arrays.
       !--------------------------------------------
       Model_Next_Year =Nudge_Beg_Year
       Model_Next_Month=Nudge_Beg_Month
       Model_Next_Day  =Nudge_Beg_Day
       Model_Next_Sec  =Nudge_Beg_Sec
       Nudge_Next_Year =Nudge_Beg_Year
       Nudge_Next_Month=Nudge_Beg_Month
       Nudge_Next_Day  =Nudge_Beg_Day
       Nudge_Next_Sec  =Nudge_Beg_Sec
     elseif(.not.Before_End) then
       ! Nudging will never occur, so switch it off
       !--------------------------------------------
       Nudge_Model=.false.
       Nudge_ON   =.false.
       write(iulog,*) ' '
       write(iulog,*) 'NUDGING: WARNING - Nudging has been requested by it will'
       write(iulog,*) 'NUDGING:           never occur for the given time values'
       write(iulog,*) ' '
     endif

     ! Initialize values for window function  
     !----------------------------------------
     lonp= 180.
     lon0=   0.
     lonn=-180.
     latp=  90.-Nudge_Hwin_lat0
     lat0=   0.
     latn= -90.-Nudge_Hwin_lat0
    
     Nudge_Hwin_lonWidthH=Nudge_Hwin_lonWidth/2.
     Nudge_Hwin_latWidthH=Nudge_Hwin_latWidth/2.

     Val1_p=(1.+tanh((Nudge_Hwin_lonWidthH+lonp)/Nudge_Hwin_lonDelta))/2.
     Val2_p=(1.+tanh((Nudge_Hwin_lonWidthH-lonp)/Nudge_Hwin_lonDelta))/2.
     Val3_p=(1.+tanh((Nudge_Hwin_latWidthH+latp)/Nudge_Hwin_latDelta))/2.
     Val4_p=(1.+tanh((Nudge_Hwin_latWidthH-latp)/Nudge_Hwin_latDelta))/2.

     Val1_0=(1.+tanh((Nudge_Hwin_lonWidthH+lon0)/Nudge_Hwin_lonDelta))/2.
     Val2_0=(1.+tanh((Nudge_Hwin_lonWidthH-lon0)/Nudge_Hwin_lonDelta))/2.
     Val3_0=(1.+tanh((Nudge_Hwin_latWidthH+lat0)/Nudge_Hwin_latDelta))/2.
     Val4_0=(1.+tanh((Nudge_Hwin_latWidthH-lat0)/Nudge_Hwin_latDelta))/2.

     Val1_n=(1.+tanh((Nudge_Hwin_lonWidthH+lonn)/Nudge_Hwin_lonDelta))/2.
     Val2_n=(1.+tanh((Nudge_Hwin_lonWidthH-lonn)/Nudge_Hwin_lonDelta))/2.
     Val3_n=(1.+tanh((Nudge_Hwin_latWidthH+latn)/Nudge_Hwin_latDelta))/2.
     Val4_n=(1.+tanh((Nudge_Hwin_latWidthH-latn)/Nudge_Hwin_latDelta))/2.

     Nudge_Hwin_max=     Val1_0*Val2_0*Val3_0*Val4_0
     Nudge_Hwin_min=min((Val1_p*Val2_p*Val3_n*Val4_n), &
                        (Val1_p*Val2_p*Val3_p*Val4_p), &
                        (Val1_n*Val2_n*Val3_n*Val4_n), &
                        (Val1_n*Val2_n*Val3_p*Val4_p))

     ! Initialize number of nudging observation values to keep track of.
     ! Allocate and initialize observation indices 
     !-----------------------------------------------------------------
     if((Nudge_Force_Opt.ge.0).and.(Nudge_Force_Opt.le.3)) then
       Nudge_NumObs=2
     else
       Nudge_NumObs=Nudge_Times_Per_Day
       if(Nudge_NumObs.lt.4) then
         write(iulog,*) 'NUDGING: Nudge_NumObs=',Nudge_NumObs
         write(iulog,*) 'NUDGING: The Diurnal Filter was forumlated for a minimum'
         write(iulog,*) 'NUDGING: of 4 observations per day. What your doing may'
         write(iulog,*) 'NUDGING: work, but you better check it first before you'
         write(iulog,*) 'NUDGING: remove this stop.'
         call endrun('NUDGING DIURNAL FILTER ONLY CONFIGURED FOR >4 OBS/DAY')
       endif
     endif
     allocate(Nudge_ObsInd(Nudge_NumObs),stat=istat)
     call alloc_err(istat,'nudging_init','Nudge_ObsInd',Nudge_NumObs)
     allocate(Nudge_File_Present(Nudge_NumObs),stat=istat)
     call alloc_err(istat,'nudging_init','Nudge_File_Present',Nudge_NumObs)
     do nn=1,Nudge_NumObs
       Nudge_ObsInd(nn) = Nudge_NumObs+1-nn
     end do
     Nudge_File_Present(:)=.false.

     ! Allocate/Initialize values for Diurnal Filter Coefs
     !------------------------------------------------------
     allocate(Nudge_Bcoef(Nudge_NumObs),stat=istat)
     call alloc_err(istat,'nudging_init','Nudge_Bcoef',Nudge_NumObs)
     allocate(Nudge_Ccoef(Nudge_NumObs),stat=istat)
     call alloc_err(istat,'nudging_init','Nudge_Ccoef',Nudge_NumObs)

     if((Nudge_Force_Opt.ge.0).and.(Nudge_Force_Opt.le.3)) then
       ! These coefs are never used.
       !------------------------------
       Nudge_Acoef   =1._r8
       Nudge_Bcoef(:)=0._r8
       Nudge_Ccoef(:)=0._r8
     elseif((Nudge_Force_Opt.ge.4).and.(Nudge_Force_Opt.le.7)) then
       ! Load Sin/Cos for 1 Diurnal Cycle
       !----------------------------------
       NumObs=float(Nudge_NumObs)
       Freq = (2._r8*SHR_CONST_PI)/NumObs
       allocate(CosVal(Nudge_NumObs),stat=istat)
       call alloc_err(istat,'nudging_init','CosVal',Nudge_NumObs)
       allocate(SinVal(Nudge_NumObs),stat=istat)
       call alloc_err(istat,'nudging_init','SinVal',Nudge_NumObs)
       do nn=1,Nudge_NumObs
         CosVal(nn)=cos(Freq*(1.5_r8-float(nn)))
         SinVal(nn)=sin(Freq*(1.5_r8-float(nn)))
       end do

       ! Load Diurnal Filter Coefs
       !----------------------------------
       Nudge_Acoef   = (NumObs-2._r8*(1._r8 - (CosVal(1)*CosVal(2)         &
                                              +SinVal(1)*SinVal(2))))/NumObs
       Nudge_Bcoef(:)=0._r8
       Nudge_Ccoef(:)=0._r8
       do nn=2,Nudge_NumObs
         Nudge_Bcoef(nn)=(Nudge_Acoef*(1._r8 + 4._r8*(CosVal(nn)*CosVal(1)        &
                                                     +SinVal(nn)*SinVal(1)))/5._r8)
         Nudge_Ccoef(nn)=2._r8*(CosVal(nn)*(CosVal(2)-CosVal(1))       &
                               +SinVal(nn)*(SinVal(2)-SinVal(1)))/NumObs
       end do

       ! For forcing options 6,7 force the square peg into the round hole...
       !  -the diurnal filter is *made* stable by using abs(Bcoef).
       !-------------------------------------------------------------------
       if((Nudge_Force_Opt.eq.6).or.(Nudge_Force_Opt.eq.7)) then
         Nudge_Bcoef(:)=abs(Nudge_Bcoef(:))
       endif

       ! For coding simplicity later on...
       !-----------------------------------
       Nudge_Ccoef(:)=Nudge_Ccoef(:)-Nudge_Bcoef(:)
       
       deallocate(CosVal)
       deallocate(SinVal)
     endif ! ((Nudge_Force_Opt.ge.4).and.(Nudge_Force_Opt.le.7)) then

     ! Initialization is done, 
     !--------------------------
     Nudge_Initialized=.true.

     ! Check that this is a valid DYCORE model
     !------------------------------------------
     if((.not.dycore_is('UNSTRUCTURED')).and. &
        (.not.dycore_is('EUL')         ).and. &
        (.not.dycore_is('LR')          )      ) then
       call endrun('NUDGING IS CURRENTLY ONLY CONFIGURED FOR CAM-SE, FV, or EUL')
     endif

     ! Informational Output
     !---------------------------
     write(iulog,*) ' '
     write(iulog,*) '---------------------------------------------------------'
     write(iulog,*) '  MODEL NUDGING INITIALIZED WITH THE FOLLOWING SETTINGS: '
     write(iulog,*) '---------------------------------------------------------'
     write(iulog,*) 'NUDGING: Nudge_Model=',Nudge_Model
     write(iulog,*) 'NUDGING: Nudge_Path=',Nudge_Path
     write(iulog,*) 'NUDGING: Nudge_File_Template =',Nudge_File_Template
     write(iulog,*) 'NUDGING: Nudge_Force_Opt=',Nudge_Force_Opt    
     write(iulog,*) 'NUDGING: Nudge_Diag_Opt =',Nudge_Diag_Opt    
     write(iulog,*) 'NUDGING: Nudge_TSmode=',Nudge_TSmode
     write(iulog,*) 'NUDGING: Nudge_Times_Per_Day=',Nudge_Times_Per_Day
     write(iulog,*) 'NUDGING: Model_Times_Per_Day=',Model_Times_Per_Day
     write(iulog,*) 'NUDGING: Nudge_Step=',Nudge_Step
     write(iulog,*) 'NUDGING: Model_Step=',Model_Step
     write(iulog,*) 'NUDGING: Nudge_Ucoef  =',Nudge_Ucoef
     write(iulog,*) 'NUDGING: Nudge_Vcoef  =',Nudge_Vcoef
     write(iulog,*) 'NUDGING: Nudge_Qcoef  =',Nudge_Qcoef
     write(iulog,*) 'NUDGING: Nudge_Tcoef  =',Nudge_Tcoef
     write(iulog,*) 'NUDGING: Nudge_PScoef =',Nudge_PScoef
     write(iulog,*) 'NUDGING: Nudge_Uprof  =',Nudge_Uprof
     write(iulog,*) 'NUDGING: Nudge_Vprof  =',Nudge_Vprof
     write(iulog,*) 'NUDGING: Nudge_Qprof  =',Nudge_Qprof
     write(iulog,*) 'NUDGING: Nudge_Tprof  =',Nudge_Tprof
     write(iulog,*) 'NUDGING: Nudge_PSprof =',Nudge_PSprof
     write(iulog,*) 'NUDGING: Nudge_Beg_Year =',Nudge_Beg_Year
     write(iulog,*) 'NUDGING: Nudge_Beg_Month=',Nudge_Beg_Month
     write(iulog,*) 'NUDGING: Nudge_Beg_Day  =',Nudge_Beg_Day
     write(iulog,*) 'NUDGING: Nudge_End_Year =',Nudge_End_Year
     write(iulog,*) 'NUDGING: Nudge_End_Month=',Nudge_End_Month
     write(iulog,*) 'NUDGING: Nudge_End_Day  =',Nudge_End_Day
     write(iulog,*) 'NUDGING: Nudge_Hwin_lo       =',Nudge_Hwin_lo
     write(iulog,*) 'NUDGING: Nudge_Hwin_hi       =',Nudge_Hwin_hi
     write(iulog,*) 'NUDGING: Nudge_Hwin_lat0     =',Nudge_Hwin_lat0
     write(iulog,*) 'NUDGING: Nudge_Hwin_latWidth =',Nudge_Hwin_latWidth
     write(iulog,*) 'NUDGING: Nudge_Hwin_latDelta =',Nudge_Hwin_latDelta
     write(iulog,*) 'NUDGING: Nudge_Hwin_lon0     =',Nudge_Hwin_lon0
     write(iulog,*) 'NUDGING: Nudge_Hwin_lonWidth =',Nudge_Hwin_lonWidth
     write(iulog,*) 'NUDGING: Nudge_Hwin_lonDelta =',Nudge_Hwin_lonDelta
     write(iulog,*) 'NUDGING: Nudge_Vwin_lo       =',Nudge_Vwin_lo
     write(iulog,*) 'NUDGING: Nudge_Vwin_hi       =',Nudge_Vwin_hi
     write(iulog,*) 'NUDGING: Nudge_Vwin_Hindex   =',Nudge_Vwin_Hindex
     write(iulog,*) 'NUDGING: Nudge_Vwin_Hdelta   =',Nudge_Vwin_Hdelta
     write(iulog,*) 'NUDGING: Nudge_Vwin_Lindex   =',Nudge_Vwin_Lindex
     write(iulog,*) 'NUDGING: Nudge_Vwin_Ldelta   =',Nudge_Vwin_Ldelta
     write(iulog,*) 'NUDGING: Nudge_Hwin_latWidthH=',Nudge_Hwin_latWidthH
     write(iulog,*) 'NUDGING: Nudge_Hwin_lonWidthH=',Nudge_Hwin_lonWidthH
     write(iulog,*) 'NUDGING: Nudge_Hwin_max      =',Nudge_Hwin_max
     write(iulog,*) 'NUDGING: Nudge_Hwin_min      =',Nudge_Hwin_min
     write(iulog,*) 'NUDGING: Nudge_Initialized   =',Nudge_Initialized
     write(iulog,*) ' '
     write(iulog,*) 'NUDGING: Nudge_NumObs=',Nudge_NumObs
     write(iulog,*) 'NUDGING: Nudge_Acoef =',Nudge_Acoef
     write(iulog,*) 'NUDGING: Nudge_Bcoef =',Nudge_Bcoef
     write(iulog,*) 'NUDGING: Nudge_Ccoef =',Nudge_Ccoef
     write(iulog,*) ' '

   endif ! (masterproc) then

   ! Broadcast other variables that have changed
   !---------------------------------------------
#ifdef SPMD
   call mpibcast(Model_Step          ,            1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Step          ,            1, mpir8 , 0, mpicom)
   call mpibcast(Model_Next_Year     ,            1, mpiint, 0, mpicom)
   call mpibcast(Model_Next_Month    ,            1, mpiint, 0, mpicom)
   call mpibcast(Model_Next_Day      ,            1, mpiint, 0, mpicom)
   call mpibcast(Model_Next_Sec      ,            1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Next_Year     ,            1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Next_Month    ,            1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Next_Day      ,            1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Next_Sec      ,            1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Model         ,            1, mpilog, 0, mpicom)
   call mpibcast(Nudge_ON            ,            1, mpilog, 0, mpicom)
   call mpibcast(Nudge_Initialized   ,            1, mpilog, 0, mpicom)
   call mpibcast(Nudge_ncol          ,            1, mpiint, 0, mpicom)
   call mpibcast(Nudge_nlev          ,            1, mpiint, 0, mpicom)
   call mpibcast(Nudge_nlon          ,            1, mpiint, 0, mpicom)
   call mpibcast(Nudge_nlat          ,            1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Hwin_max      ,            1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_min      ,            1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_lonWidthH,            1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_latWidthH,            1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_NumObs        ,            1, mpiint, 0, mpicom)
#endif

   ! All non-masterproc processes also need to allocate space
   ! before the broadcast of Nudge_NumObs dependent data.
   !------------------------------------------------------------
   if(.not.masterproc) then
     allocate(Nudge_ObsInd(Nudge_NumObs),stat=istat)
     call alloc_err(istat,'nudging_init','Nudge_ObsInd',Nudge_NumObs)
     allocate(Nudge_File_Present(Nudge_NumObs),stat=istat)
     call alloc_err(istat,'nudging_init','Nudge_File_Present',Nudge_NumObs)
     allocate(Nudge_Bcoef(Nudge_NumObs),stat=istat)
     call alloc_err(istat,'nudging_init','Nudge_Bcoef',Nudge_NumObs)
     allocate(Nudge_Ccoef(Nudge_NumObs),stat=istat)
     call alloc_err(istat,'nudging_init','Nudge_Ccoef',Nudge_NumObs)
   endif
#ifdef SPMD
   call mpibcast(Nudge_ObsInd        , Nudge_NumObs, mpiint, 0, mpicom)
   call mpibcast(Nudge_File_Present  , Nudge_NumObs, mpilog, 0, mpicom)
   call mpibcast(Nudge_Acoef         ,            1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Bcoef         , Nudge_NumObs, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Ccoef         , Nudge_NumObs, mpir8 , 0, mpicom)
#endif

   ! Allocate Space for Nudging observation arrays, initialize with 0's
   !---------------------------------------------------------------------
   allocate(Nobs_U(pcols,pver,begchunk:endchunk,Nudge_NumObs),stat=istat)
   call alloc_err(istat,'nudging_init','Nobs_U',pcols*pver*((endchunk-begchunk)+1)*Nudge_NumObs)
   allocate(Nobs_V(pcols,pver,begchunk:endchunk,Nudge_NumObs),stat=istat)
   call alloc_err(istat,'nudging_init','Nobs_V',pcols*pver*((endchunk-begchunk)+1)*Nudge_NumObs)
   allocate(Nobs_T(pcols,pver,begchunk:endchunk,Nudge_NumObs),stat=istat)
   call alloc_err(istat,'nudging_init','Nobs_T',pcols*pver*((endchunk-begchunk)+1)*Nudge_NumObs)
   allocate(Nobs_Q(pcols,pver,begchunk:endchunk,Nudge_NumObs),stat=istat)
   call alloc_err(istat,'nudging_init','Nobs_Q',pcols*pver*((endchunk-begchunk)+1)*Nudge_NumObs)
   allocate(Nobs_PS(pcols,begchunk:endchunk,Nudge_NumObs),stat=istat)
   call alloc_err(istat,'nudging_init','Nobs_PS',pcols*((endchunk-begchunk)+1)*Nudge_NumObs)
   allocate(Mobs_U(pcols,pver,begchunk:endchunk,Nudge_NumObs),stat=istat)

   call alloc_err(istat,'nudging_init','Mobs_U',pcols*pver*((endchunk-begchunk)+1)*Nudge_NumObs)
   allocate(Mobs_V(pcols,pver,begchunk:endchunk,Nudge_NumObs),stat=istat)
   call alloc_err(istat,'nudging_init','Mobs_V',pcols*pver*((endchunk-begchunk)+1)*Nudge_NumObs)
   allocate(Mobs_T(pcols,pver,begchunk:endchunk,Nudge_NumObs),stat=istat)
   call alloc_err(istat,'nudging_init','Mobs_T',pcols*pver*((endchunk-begchunk)+1)*Nudge_NumObs)
   allocate(Mobs_Q(pcols,pver,begchunk:endchunk,Nudge_NumObs),stat=istat)
   call alloc_err(istat,'nudging_init','Mobs_Q',pcols*pver*((endchunk-begchunk)+1)*Nudge_NumObs)
   allocate(Mobs_PS(pcols,begchunk:endchunk,Nudge_NumObs),stat=istat)
   call alloc_err(istat,'nudging_init','Mobs_PS',pcols*((endchunk-begchunk)+1)*Nudge_NumObs)

   Nobs_U(:pcols,:pver,begchunk:endchunk,:Nudge_NumObs)=0._r8
   Nobs_V(:pcols,:pver,begchunk:endchunk,:Nudge_NumObs)=0._r8
   Nobs_T(:pcols,:pver,begchunk:endchunk,:Nudge_NumObs)=0._r8
   Nobs_Q(:pcols,:pver,begchunk:endchunk,:Nudge_NumObs)=0._r8
   Nobs_PS(:pcols     ,begchunk:endchunk,:Nudge_NumObs)=0._r8
   Mobs_U(:pcols,:pver,begchunk:endchunk,:Nudge_NumObs)=0._r8
   Mobs_V(:pcols,:pver,begchunk:endchunk,:Nudge_NumObs)=0._r8
   Mobs_T(:pcols,:pver,begchunk:endchunk,:Nudge_NumObs)=0._r8
   Mobs_Q(:pcols,:pver,begchunk:endchunk,:Nudge_NumObs)=0._r8
   Mobs_PS(:pcols     ,begchunk:endchunk,:Nudge_NumObs)=0._r8

!!DIAG
   if(masterproc) then
     write(iulog,*) 'NUDGING: nudging_init() OBS arrays allocated and initialized'
     write(iulog,*) 'NUDGING: nudging_init() SIZE#',(9*pcols*pver*((endchunk-begchunk)+1)*Nudge_NumObs)
     write(iulog,*) 'NUDGING: nudging_init() MB:',float(8*9*pcols*pver*((endchunk-begchunk)+1)*Nudge_NumObs)/(1024.*1024.)
     write(iulog,*) 'NUDGING: nudging_init() pcols=',pcols,' pver=',pver
     write(iulog,*) 'NUDGING: nudging_init() begchunk:',begchunk,' endchunk=',endchunk
     write(iulog,*) 'NUDGING: nudging_init() chunk:',(endchunk-begchunk+1),' Nudge_NumObs=',Nudge_NumObs
     write(iulog,*) 'NUDGING: nudging_init() Nudge_ObsInd=',Nudge_ObsInd
     write(iulog,*) 'NUDGING: nudging_init() Nudge_File_Present=',Nudge_File_Present
     write(iulog,*) 'NUDGING: nudging_init() Nudge_Acoef=',Nudge_Acoef
     write(iulog,*) 'NUDGING: nudging_init() Nudge_Bcoef=',Nudge_Bcoef
     write(iulog,*) 'NUDGING: nudging_init() Nudge_Ccoef=',Nudge_Ccoef
   endif
!!DIAG

   ! Initialize Nudging Coeffcient profiles in local arrays
   ! Load zeros into nudging arrays
   !------------------------------------------------------
   do lchnk=begchunk,endchunk
     ncol=get_ncols_p(lchnk)
     do icol=1,ncol
       rlat=get_rlat_p(lchnk,icol)*180._r8/SHR_CONST_PI
       rlon=get_rlon_p(lchnk,icol)*180._r8/SHR_CONST_PI

       call nudging_set_profile(rlat,rlon,Nudge_Uprof,Wprof,pver)
       Nudge_Utau(icol,:,lchnk)=Wprof(:)
       call nudging_set_profile(rlat,rlon,Nudge_Vprof,Wprof,pver)
       Nudge_Vtau(icol,:,lchnk)=Wprof(:)
       call nudging_set_profile(rlat,rlon,Nudge_Tprof,Wprof,pver)
       Nudge_Stau(icol,:,lchnk)=Wprof(:)
       call nudging_set_profile(rlat,rlon,Nudge_Qprof,Wprof,pver)
       Nudge_Qtau(icol,:,lchnk)=Wprof(:)

       Nudge_PStau(icol,lchnk)=nudging_set_PSprofile(rlat,rlon,Nudge_PSprof)
     end do
     Nudge_Utau(:ncol,:pver,lchnk) =                             &
     Nudge_Utau(:ncol,:pver,lchnk) * Nudge_Ucoef/float(Nudge_Step)
     Nudge_Vtau(:ncol,:pver,lchnk) =                             &
     Nudge_Vtau(:ncol,:pver,lchnk) * Nudge_Vcoef/float(Nudge_Step)
     Nudge_Stau(:ncol,:pver,lchnk) =                             &
     Nudge_Stau(:ncol,:pver,lchnk) * Nudge_Tcoef/float(Nudge_Step)
     Nudge_Qtau(:ncol,:pver,lchnk) =                             &
     Nudge_Qtau(:ncol,:pver,lchnk) * Nudge_Qcoef/float(Nudge_Step)
     Nudge_PStau(:ncol,lchnk)=                             &
     Nudge_PStau(:ncol,lchnk)* Nudge_PScoef/float(Nudge_Step)

     Nudge_Ustep(:pcols,:pver,lchnk)=0._r8
     Nudge_Vstep(:pcols,:pver,lchnk)=0._r8
     Nudge_Sstep(:pcols,:pver,lchnk)=0._r8
     Nudge_Qstep(:pcols,:pver,lchnk)=0._r8
     Nudge_PSstep(:pcols,lchnk)=0._r8
     Target_U(:pcols,:pver,lchnk)=0._r8
     Target_V(:pcols,:pver,lchnk)=0._r8
     Target_T(:pcols,:pver,lchnk)=0._r8
     Target_S(:pcols,:pver,lchnk)=0._r8
     Target_Q(:pcols,:pver,lchnk)=0._r8
     Target_PS(:pcols,lchnk)=0._r8
   end do

!DIAG
!  if(masterproc) then
!    write(iulog,*) 'NUDGING: exiting nudging_init()'
!  endif
!DIAG

   ! End Routine
   !------------
   return
  end subroutine ! nudging_init
  !================================================================


  !================================================================
  subroutine nudging_timestep_init(phys_state)
   ! 
   ! NUDGING_TIMESTEP_INIT: 
   !                 Check the current time and update Model/Nudging 
   !                 arrays when necessary. Toggle the Nudging flag
   !                 when the time is withing the nudging window.
   !===============================================================
   use physconst    ,only: cpair
   use physics_types,only: physics_state
   use constituents ,only: cnst_get_ind
   use dycore       ,only: dycore_is
   use ppgrid       ,only: pver,pcols,begchunk,endchunk
   use filenames    ,only: interpret_filename_spec
   use ESMF

   ! Arguments
   !-----------
   type(physics_state),intent(in):: phys_state(begchunk:endchunk)

   ! Local values
   !----------------
   integer Year,Month,Day,Sec
   integer YMD1,YMD2,YMD
   logical Update_Model,Update_Nudge,Sync_Error
   logical After_Beg   ,Before_End
   integer lchnk,ncol,indw

   type(ESMF_Time)         Date1,Date2
   type(ESMF_TimeInterval) DateDiff
   integer                 DeltaT
   real(r8)                Tscale
   integer                 rc
   integer                 nn
   integer                 kk
   real(r8)                Sbar,Qbar,Wsum

!DIAG
!  if(masterproc) then
!    write(iulog,*) 'NUDGING: entering nudging_timestep_init()'
!  endif
!DIAG
   ! Check if Nudging is initialized
   !---------------------------------
   if(.not.Nudge_Initialized) then
     call endrun('nudging_timestep_init:: Nudging NOT Initialized')
   endif

   ! Get Current time
   !--------------------
   call get_curr_date(Year,Month,Day,Sec)
   YMD=(Year*10000) + (Month*100) + Day

   !-------------------------------------------------------
   ! Determine if the current time is AFTER the begining time
   ! and if it is BEFORE the ending time.
   !-------------------------------------------------------
   YMD1=(Nudge_Beg_Year*10000) + (Nudge_Beg_Month*100) + Nudge_Beg_Day
   call timemgr_time_ge(YMD1,Nudge_Beg_Sec,         &
                        YMD ,Sec          ,After_Beg)

   YMD1=(Nudge_End_Year*10000) + (Nudge_End_Month*100) + Nudge_End_Day
   call timemgr_time_ge(YMD ,Sec,                    &
                        YMD1,Nudge_End_Sec,Before_End)

   !--------------------------------------------------------------
   ! When past the NEXT time, Update Model Arrays and time indices
   !--------------------------------------------------------------
   YMD1=(Model_Next_Year*10000) + (Model_Next_Month*100) + Model_Next_Day
   call timemgr_time_ge(YMD1,Model_Next_Sec,            &
                        YMD ,Sec           ,Update_Model)

   if((Before_End).and.(Update_Model)) then
     ! Increment the Model times by the current interval
     !---------------------------------------------------
     Model_Curr_Year =Model_Next_Year
     Model_Curr_Month=Model_Next_Month
     Model_Curr_Day  =Model_Next_Day
     Model_Curr_Sec  =Model_Next_Sec
     YMD1=(Model_Curr_Year*10000) + (Model_Curr_Month*100) + Model_Curr_Day
     call timemgr_time_inc(YMD1,Model_Curr_Sec,              &
                           YMD2,Model_Next_Sec,Model_Step,0,0)

     ! Check for Sync Error where NEXT model time after the update
     ! is before the current time. If so, reset the next model 
     ! time to a Model_Step after the current time.
     !--------------------------------------------------------------
     call timemgr_time_ge(YMD2,Model_Next_Sec,            &
                          YMD ,Sec           ,Sync_Error)
     if(Sync_Error) then
       Model_Curr_Year =Year
       Model_Curr_Month=Month
       Model_Curr_Day  =Day
       Model_Curr_Sec  =Sec
       call timemgr_time_inc(YMD ,Model_Curr_Sec,              &
                             YMD2,Model_Next_Sec,Model_Step,0,0)
       write(iulog,*) 'NUDGING: WARNING - Model_Time Sync ERROR... CORRECTED'
     endif
     Model_Next_Year =(YMD2/10000)
     YMD2            = YMD2-(Model_Next_Year*10000)
     Model_Next_Month=(YMD2/100)
     Model_Next_Day  = YMD2-(Model_Next_Month*100)

     ! Load values at Current into the Model arrays
     !-----------------------------------------------
     call cnst_get_ind('Q',indw)
     do lchnk=begchunk,endchunk
       ncol=phys_state(lchnk)%ncol
       Model_U(:ncol,:pver,lchnk)=phys_state(lchnk)%u(:ncol,:pver)
       Model_V(:ncol,:pver,lchnk)=phys_state(lchnk)%v(:ncol,:pver)
       Model_T(:ncol,:pver,lchnk)=phys_state(lchnk)%t(:ncol,:pver)
       Model_Q(:ncol,:pver,lchnk)=phys_state(lchnk)%q(:ncol,:pver,indw)
       Model_PS(:ncol,lchnk)=phys_state(lchnk)%ps(:ncol)
     end do

     ! Load Dry Static Energy values for Model
     !-----------------------------------------
     if(Nudge_TSmode.eq.0) then
       ! DSE tendencies from Temperature only
       !---------------------------------------
       do lchnk=begchunk,endchunk
         ncol=phys_state(lchnk)%ncol
         Model_S(:ncol,:pver,lchnk)=cpair*Model_T(:ncol,:pver,lchnk)
       end do
     elseif(Nudge_TSmode.eq.1) then
       ! Caluculate DSE tendencies from Temperature, Water Vapor, and Surface Pressure
       !------------------------------------------------------------------------------
       do lchnk=begchunk,endchunk
         ncol=phys_state(lchnk)%ncol
         call calc_DryStaticEnergy(Model_T(:,:,lchnk)  , Model_Q(:,:,lchnk), &
                                 phys_state(lchnk)%phis,  Model_PS(:,lchnk), &
                                                  Model_S(:,:,lchnk), ncol)
       end do
     endif 
!PFCDIAG
!    if(.FALSE.) then
!      ! OPTIONALLY remove vertical mean from Model_S
!      !-----------------------------------------------
!      do lchnk=begchunk,endchunk
!        ncol=phys_state(lchnk)%ncol
!        do nn=1,ncol
!          Sbar=0._r8
!          do kk=1,pver
!            Sbar=Sbar+Model_S(nn,kk,lchnk)
!          end do
!          Sbar=Sbar/float(pver)
!          Model_S(nn,:,lchnk)=Model_S(nn,:,lchnk)-Sbar
!        end do
!      end do
!    elseif(.TRUE.) then
!      ! OPTIONALLY remove weighted vertical mean from Model_S
!      !-------------------------------------------------------
!      do lchnk=begchunk,endchunk
!        ncol=phys_state(lchnk)%ncol
!        do nn=1,ncol
!          Sbar=0._r8
!          Qbar=0._r8
!          Wsum=0._r8
!          do kk=1,pver
!            Sbar=Sbar+Model_S(nn,kk,lchnk)*(phys_state(lchnk)%pdel(nn,kk)/phys_state(lchnk)%pmid(nn,kk))
!            Qbar=Qbar+Model_Q(nn,kk,lchnk)*(phys_state(lchnk)%pdel(nn,kk)/phys_state(lchnk)%pmid(nn,kk))
!            Wsum=Wsum+(phys_state(lchnk)%pdel(nn,kk)/phys_state(lchnk)%pmid(nn,kk))
!          end do
!          Sbar=Sbar/Wsum
!          Qbar=Qbar/Wsum
!          Model_S(nn,:,lchnk)=Model_S(nn,:,lchnk)-Sbar
!          Model_Q(nn,:,lchnk)=Model_Q(nn,:,lchnk)-Qbar
!        end do
!      end do
!    endif
!PFCDIAG

   endif ! ((Before_End).and.(Update_Model)) then

   !----------------------------------------------------------------
   ! When past the NEXT time, Update Nudging Arrays and time indices
   !----------------------------------------------------------------
   YMD1=(Nudge_Next_Year*10000) + (Nudge_Next_Month*100) + Nudge_Next_Day
   call timemgr_time_ge(YMD1,Nudge_Next_Sec,            &
                        YMD ,Sec           ,Update_Nudge)

   if((Before_End).and.(Update_Nudge)) then
     ! Increment the Nudge times by the current interval
     !---------------------------------------------------
     Nudge_Curr_Year =Nudge_Next_Year
     Nudge_Curr_Month=Nudge_Next_Month
     Nudge_Curr_Day  =Nudge_Next_Day
     Nudge_Curr_Sec  =Nudge_Next_Sec
     YMD1=(Nudge_Curr_Year*10000) + (Nudge_Curr_Month*100) + Nudge_Curr_Day
     call timemgr_time_inc(YMD1,Nudge_Curr_Sec,              &
                           YMD2,Nudge_Next_Sec,Nudge_Step,0,0)
     Nudge_Next_Year =(YMD2/10000)
     YMD2            = YMD2-(Nudge_Next_Year*10000)
     Nudge_Next_Month=(YMD2/100)
     Nudge_Next_Day  = YMD2-(Nudge_Next_Month*100)

     ! Set the analysis filename at the NEXT time.
     !---------------------------------------------------------------
     Nudge_File=interpret_filename_spec(Nudge_File_Template      , &
                                         yr_spec=Nudge_Next_Year , &
                                        mon_spec=Nudge_Next_Month, &
                                        day_spec=Nudge_Next_Day  , &
                                        sec_spec=Nudge_Next_Sec    )
     if(masterproc) then
      write(iulog,*) 'NUDGING: Reading analyses:',trim(Nudge_Path)//trim(Nudge_File)
     endif

     ! Rotate Nudge_ObsInd() indices for new data, then update 
     ! the Nudge observation arrays with analysis data at the 
     ! NEXT==Nudge_ObsInd(1) time.
     !----------------------------------------------------------
     if(dycore_is('UNSTRUCTURED')) then
       call nudging_update_analyses_se (trim(Nudge_Path)//trim(Nudge_File))
     elseif(dycore_is('EUL')) then
       call nudging_update_analyses_eul(trim(Nudge_Path)//trim(Nudge_File))
     else !if(dycore_is('LR')) then
       call nudging_update_analyses_fv (trim(Nudge_Path)//trim(Nudge_File))
     endif

     ! Update the Model observation arrays with model data at 
     ! the CURR==Nudge_ObsInd(2) time.
     !---------------------------------------------------------------
     call cnst_get_ind('Q',indw)
     do lchnk=begchunk,endchunk
       ncol=phys_state(lchnk)%ncol
       Mobs_U(:ncol,:pver,lchnk,Nudge_ObsInd(2))=phys_state(lchnk)%u(:ncol,:pver)
       Mobs_V(:ncol,:pver,lchnk,Nudge_ObsInd(2))=phys_state(lchnk)%v(:ncol,:pver)
       Mobs_T(:ncol,:pver,lchnk,Nudge_ObsInd(2))=phys_state(lchnk)%t(:ncol,:pver)
       Mobs_Q(:ncol,:pver,lchnk,Nudge_ObsInd(2))=phys_state(lchnk)%q(:ncol,:pver,indw)
       Mobs_PS(:ncol     ,lchnk,Nudge_ObsInd(2))=phys_state(lchnk)%ps(:ncol)
     end do

     ! Now Load the Target values for nudging tendencies
     !---------------------------------------------------
     if((Nudge_Force_Opt.eq.0).or.(Nudge_Force_Opt.eq.1)) then
       ! Target is OBS data at NEXT time
       !----------------------------------
       do lchnk=begchunk,endchunk
         ncol=phys_state(lchnk)%ncol
         Target_U(:ncol,:pver,lchnk)=Nobs_U(:ncol,:pver,lchnk,Nudge_ObsInd(1))
         Target_V(:ncol,:pver,lchnk)=Nobs_V(:ncol,:pver,lchnk,Nudge_ObsInd(1))
         Target_T(:ncol,:pver,lchnk)=Nobs_T(:ncol,:pver,lchnk,Nudge_ObsInd(1))
         Target_Q(:ncol,:pver,lchnk)=Nobs_Q(:ncol,:pver,lchnk,Nudge_ObsInd(1))
         Target_PS(:ncol     ,lchnk)=Nobs_PS(:ncol     ,lchnk,Nudge_ObsInd(1))
       end do
     elseif((Nudge_Force_Opt.eq.2).or.(Nudge_Force_Opt.eq.3)) then
       ! Target is OBS data at CURR time
       !----------------------------------
       do lchnk=begchunk,endchunk
         ncol=phys_state(lchnk)%ncol
         Target_U(:ncol,:pver,lchnk)=Nobs_U(:ncol,:pver,lchnk,Nudge_ObsInd(2))
         Target_V(:ncol,:pver,lchnk)=Nobs_V(:ncol,:pver,lchnk,Nudge_ObsInd(2))
         Target_T(:ncol,:pver,lchnk)=Nobs_T(:ncol,:pver,lchnk,Nudge_ObsInd(2))
         Target_Q(:ncol,:pver,lchnk)=Nobs_Q(:ncol,:pver,lchnk,Nudge_ObsInd(2))
         Target_PS(:ncol     ,lchnk)=Nobs_PS(:ncol     ,lchnk,Nudge_ObsInd(2))
       end do
     elseif((Nudge_Force_Opt.ge.4).or.(Nudge_Force_Opt.le.7)) then
       ! Target is Diurnal Estimate at NEXT time
       !-----------------------------------------
       do lchnk=begchunk,endchunk
         ncol=phys_state(lchnk)%ncol
         Target_U(:ncol,:pver,lchnk)=Nudge_Acoef*Nobs_U(:ncol,:pver,lchnk,Nudge_ObsInd(1))
         Target_V(:ncol,:pver,lchnk)=Nudge_Acoef*Nobs_V(:ncol,:pver,lchnk,Nudge_ObsInd(1))
         Target_T(:ncol,:pver,lchnk)=Nudge_Acoef*Nobs_T(:ncol,:pver,lchnk,Nudge_ObsInd(1))
         Target_Q(:ncol,:pver,lchnk)=Nudge_Acoef*Nobs_Q(:ncol,:pver,lchnk,Nudge_ObsInd(1))
         Target_PS(:ncol     ,lchnk)=Nudge_Acoef*Nobs_PS(:ncol     ,lchnk,Nudge_ObsInd(1))
         do nn=2,Nudge_NumObs
           Target_U(:ncol,:pver,lchnk) = Target_U(:ncol,:pver,lchnk)                   &
                         +(Nudge_Bcoef(nn)*Nobs_U(:ncol,:pver,lchnk,Nudge_ObsInd(nn))) &
                         +(Nudge_Ccoef(nn)*Mobs_U(:ncol,:pver,lchnk,Nudge_ObsInd(nn)))
           Target_V(:ncol,:pver,lchnk) = Target_V(:ncol,:pver,lchnk)                   &
                         +(Nudge_Bcoef(nn)*Nobs_V(:ncol,:pver,lchnk,Nudge_ObsInd(nn))) &
                         +(Nudge_Ccoef(nn)*Mobs_V(:ncol,:pver,lchnk,Nudge_ObsInd(nn)))
           Target_T(:ncol,:pver,lchnk) = Target_T(:ncol,:pver,lchnk)                   &
                         +(Nudge_Bcoef(nn)*Nobs_T(:ncol,:pver,lchnk,Nudge_ObsInd(nn))) &
                         +(Nudge_Ccoef(nn)*Mobs_T(:ncol,:pver,lchnk,Nudge_ObsInd(nn)))
           Target_Q(:ncol,:pver,lchnk) = Target_Q(:ncol,:pver,lchnk)                   &
                         +(Nudge_Bcoef(nn)*Nobs_Q(:ncol,:pver,lchnk,Nudge_ObsInd(nn))) &
                         +(Nudge_Ccoef(nn)*Mobs_Q(:ncol,:pver,lchnk,Nudge_ObsInd(nn)))
           Target_PS(:ncol     ,lchnk) =Target_PS(:ncol      ,lchnk)                   &
                        +(Nudge_Bcoef(nn)*Nobs_PS(:ncol      ,lchnk,Nudge_ObsInd(nn))) &
                        +(Nudge_Ccoef(nn)*Mobs_PS(:ncol      ,lchnk,Nudge_ObsInd(nn)))
         end do
       end do ! lchnk=begchunk,endchunk
     else
       write(iulog,*) 'NUDGING: Unknown Nudge_Force_Opt=',Nudge_Force_Opt
       call endrun('nudging_timestep_init:: ERROR unknown Nudging_Force_Opt')
     endif

     ! Now load Dry Static Energy values for Target
     !---------------------------------------------
     if(Nudge_TSmode.eq.0) then
       ! DSE tendencies from Temperature only
       !---------------------------------------
       do lchnk=begchunk,endchunk
         ncol=phys_state(lchnk)%ncol
         Target_S(:ncol,:pver,lchnk)=cpair*Target_T(:ncol,:pver,lchnk)
       end do
     elseif(Nudge_TSmode.eq.1) then
       ! Caluculate DSE tendencies from Temperature, Water Vapor, and Surface Pressure
       !------------------------------------------------------------------------------
       do lchnk=begchunk,endchunk
         ncol=phys_state(lchnk)%ncol
         call calc_DryStaticEnergy(Target_T(:,:,lchnk), Target_Q(:,:,lchnk), &
                                 phys_state(lchnk)%phis, Target_PS(:,lchnk), &
                                                  Target_S(:,:,lchnk), ncol)
       end do
     endif
!PFCDIAG
!    if(.FALSE.) then
!      ! OPTIONALLY remove vertical mean from Target_S
!      !-----------------------------------------------
!      do lchnk=begchunk,endchunk
!        ncol=phys_state(lchnk)%ncol
!        do nn=1,ncol
!          Sbar=0._r8
!          do kk=1,pver
!            Sbar=Sbar+Target_S(nn,kk,lchnk)
!          end do
!          Sbar=Sbar/float(pver)
!          Target_S(nn,:,lchnk)=Target_S(nn,:,lchnk)-Sbar
!        end do
!      end do
!    elseif(.TRUE.) then
!      ! OPTIONALLY remove weighted vertical mean from Target_S
!      !-------------------------------------------------------
!      do lchnk=begchunk,endchunk
!        ncol=phys_state(lchnk)%ncol
!        do nn=1,ncol
!          Sbar=0._r8
!          Qbar=0._r8
!          Wsum=0._r8
!          do kk=1,pver
!            Sbar=Sbar+Target_S(nn,kk,lchnk)*(phys_state(lchnk)%pdel(nn,kk)/phys_state(lchnk)%pmid(nn,kk))
!            Qbar=Qbar+Target_Q(nn,kk,lchnk)*(phys_state(lchnk)%pdel(nn,kk)/phys_state(lchnk)%pmid(nn,kk))
!            Wsum=Wsum+(phys_state(lchnk)%pdel(nn,kk)/phys_state(lchnk)%pmid(nn,kk))
!          end do
!          Sbar=Sbar/Wsum
!          Qbar=Qbar/Wsum
!          Target_S(nn,:,lchnk)=Target_S(nn,:,lchnk)-Sbar
!          Target_Q(nn,:,lchnk)=Target_Q(nn,:,lchnk)-Qbar
!        end do
!      end do
!    endif
!PFCDIAG

   endif ! ((Before_End).and.(Update_Nudge)) then

   !----------------------------------------------------------------
   ! Toggle Nudging flag when the time interval is between 
   ! beginning and ending times, and all of the analyses files exist.
   !----------------------------------------------------------------
   if((After_Beg).and.(Before_End)) then
     if(    (Nudge_Force_Opt.eq.0).or.(Nudge_Force_Opt.eq.1)) then
       ! Verify that the NEXT analyses are available
       !---------------------------------------------
       Nudge_ON=Nudge_File_Present(Nudge_ObsInd(1))
     elseif((Nudge_Force_Opt.eq.2).or.(Nudge_Force_Opt.eq.3)) then
       ! Verify that the CURR analyses are available
       !---------------------------------------------
       Nudge_ON=Nudge_File_Present(Nudge_ObsInd(2))
     else
       ! Verify that the ALL analyses are available
       !---------------------------------------------
       Nudge_ON=.true.
       do nn=1,Nudge_NumObs
         if(.not.Nudge_File_Present(nn)) Nudge_ON=.false.
       end do
     endif
     if(.not.Nudge_ON) then
       if(masterproc) then
         write(iulog,*) 'NUDGING: WARNING - analyses file NOT FOUND. Switching '
         write(iulog,*) 'NUDGING:           nudging OFF to coast thru the gap. '
       endif
     endif
   else
     Nudge_ON=.false.
   endif

   !-------------------------------------------------------
   ! HERE Implement time dependence of Nudging Coefs HERE
   !-------------------------------------------------------


   !---------------------------------------------------
   ! If Data arrays have changed update stepping arrays
   !---------------------------------------------------
   if((Before_End).and.((Update_Nudge).or.(Update_Model))) then

     ! Set Tscale for the specified Forcing Option 
     !-----------------------------------------------
     if((Nudge_Force_Opt.eq.0).or.(Nudge_Force_Opt.eq.2).or. &
        (Nudge_Force_Opt.eq.4).or.(Nudge_Force_Opt.eq.6)     ) then
       Tscale=1._r8
     elseif(Nudge_Force_Opt.eq.3) then
       call ESMF_TimeSet(Date1,YY=Year,MM=Month,DD=Day,S=Sec)
       call ESMF_TimeSet(Date2,YY=Nudge_Curr_Year,MM=Nudge_Curr_Month, &
                               DD=Nudge_Curr_Day , S=Nudge_Curr_Sec    )
       DateDiff =Date1-Date2
       call ESMF_TimeIntervalGet(DateDiff,S=DeltaT,rc=rc)
       DeltaT=DeltaT+dtime
       Tscale=float(Nudge_Step)/float(DeltaT)
     elseif((Nudge_Force_Opt.eq.1).or.(Nudge_Force_Opt.eq.5).or. &
            (Nudge_Force_Opt.eq.7)                               ) then
       call ESMF_TimeSet(Date1,YY=Year,MM=Month,DD=Day,S=Sec)
       call ESMF_TimeSet(Date2,YY=Nudge_Next_Year,MM=Nudge_Next_Month, &
                               DD=Nudge_Next_Day , S=Nudge_Next_Sec    )
       DateDiff =Date2-Date1
       call ESMF_TimeIntervalGet(DateDiff,S=DeltaT,rc=rc)
       Tscale=float(Nudge_Step)/float(DeltaT)
     endif

     ! Update the nudging tendencies
     !--------------------------------
     do lchnk=begchunk,endchunk
       ncol=phys_state(lchnk)%ncol
       Nudge_Ustep(:ncol,:pver,lchnk)=(  Target_U(:ncol,:pver,lchnk)      &
                                         -Model_U(:ncol,:pver,lchnk))     &
                                      *Tscale*Nudge_Utau(:ncol,:pver,lchnk)
       Nudge_Vstep(:ncol,:pver,lchnk)=(  Target_V(:ncol,:pver,lchnk)      &
                                         -Model_V(:ncol,:pver,lchnk))     &
                                      *Tscale*Nudge_Vtau(:ncol,:pver,lchnk)
       Nudge_Sstep(:ncol,:pver,lchnk)=(  Target_S(:ncol,:pver,lchnk)      &
                                         -Model_S(:ncol,:pver,lchnk))     &
                                      *Tscale*Nudge_Stau(:ncol,:pver,lchnk)
       Nudge_Qstep(:ncol,:pver,lchnk)=(  Target_Q(:ncol,:pver,lchnk)      &
                                         -Model_Q(:ncol,:pver,lchnk))     &
                                      *Tscale*Nudge_Qtau(:ncol,:pver,lchnk)
       Nudge_PSstep(:ncol,     lchnk)=(  Target_PS(:ncol,lchnk)      &
                                         -Model_PS(:ncol,lchnk))     &
                                      *Tscale*Nudge_PStau(:ncol,lchnk)
     end do

     !******************
     ! DIAG
     !******************
!    if(masterproc) then
!      write(iulog,*) 'PFC: Target_T(1,:pver,begchunk)=',Target_T(1,:pver,begchunk)  
!      write(iulog,*) 'PFC:  Model_T(1,:pver,begchunk)=',Model_T(1,:pver,begchunk)
!      write(iulog,*) 'PFC: Target_S(1,:pver,begchunk)=',Target_S(1,:pver,begchunk)  
!      write(iulog,*) 'PFC:  Model_S(1,:pver,begchunk)=',Model_S(1,:pver,begchunk)
!      write(iulog,*) 'PFC:      Target_PS(1,begchunk)=',Target_PS(1,begchunk)  
!      write(iulog,*) 'PFC:       Model_PS(1,begchunk)=',Model_PS(1,begchunk)
!      write(iulog,*) 'PFC: Nudge_Sstep(1,:pver,begchunk)=',Nudge_Sstep(1,:pver,begchunk)
!      write(iulog,*) 'PFC: Nudge_Xstep arrays updated:'
!    endif
   endif ! ((Before_End).and.((Update_Nudge).or.(Update_Model))) then

!DIAG
!  if(masterproc) then
!    write(iulog,*) 'NUDGING: exiting nudging_timestep_init()'
!  endif
!DIAG

   ! End Routine
   !------------
   return
  end subroutine ! nudging_timestep_init
  !================================================================


  !================================================================
  subroutine nudging_timestep_tend(phys_state,phys_tend)
   ! 
   ! NUDGING_TIMESTEP_TEND: 
   !                If Nudging is ON, return the Nudging contributions 
   !                to forcing using the current contents of the Nudge 
   !                arrays. Send output to the cam history module as well.
   !===============================================================
   use physconst    ,only: cpair
   use physics_types,only: physics_state,physics_ptend,physics_ptend_init
   use constituents ,only: cnst_get_ind,pcnst
   use ppgrid       ,only: pver,pcols,begchunk,endchunk
   use cam_history  ,only: outfld

   ! Arguments
   !-------------
   type(physics_state), intent(in) :: phys_state
   type(physics_ptend), intent(out):: phys_tend

   ! Local values
   !--------------------
   integer indw,ncol,lchnk
   logical lq(pcnst)

!DIAG
!  if(masterproc) then
!    write(iulog,*) 'NUDGING: entering nudging_timestep_tend()'
!  endif
!DIAG
   call cnst_get_ind('Q',indw)
   lq(:)   =.false.
   lq(indw)=.true.
   call physics_ptend_init(phys_tend,phys_state%psetcols,'nudging',lu=.true.,lv=.true.,ls=.true.,lq=lq)

   if(Nudge_ON) then
     lchnk=phys_state%lchnk
     ncol =phys_state%ncol
     phys_tend%u(:ncol,:pver)     =Nudge_Ustep(:ncol,:pver,lchnk)
     phys_tend%v(:ncol,:pver)     =Nudge_Vstep(:ncol,:pver,lchnk)
     phys_tend%s(:ncol,:pver)     =Nudge_Sstep(:ncol,:pver,lchnk)
     phys_tend%q(:ncol,:pver,indw)=Nudge_Qstep(:ncol,:pver,lchnk)

     call outfld('Nudge_U',phys_tend%u          ,pcols,lchnk)
     call outfld('Nudge_V',phys_tend%v          ,pcols,lchnk)
     call outfld('Nudge_T',phys_tend%s          ,pcols,lchnk)
     call outfld('Nudge_Q',phys_tend%q(1,1,indw),pcols,lchnk)
   endif

!DIAG
!  if(masterproc) then
!    write(iulog,*) 'NUDGING: exiting nudging_timestep_tend()'
!  endif
!DIAG
   ! End Routine
   !------------
   return
  end subroutine ! nudging_timestep_tend
  !================================================================


  !================================================================
  subroutine nudging_update_analyses_se(anal_file)
   ! 
   ! NUDGING_UPDATE_ANALYSES_SE: 
   !                 Open the given analyses data file, read in 
   !                 U,V,T,Q, and PS values and then distribute
   !                 the values to all of the chunks.
   !===============================================================
   use ppgrid ,only: pver,begchunk
   use netcdf

   ! Arguments
   !-------------
   character(len=*),intent(in):: anal_file

   ! Local values
   !-------------
   integer lev
   integer ncol,plev,istat
   integer ncid,varid
   real(r8) Xanal(Nudge_ncol,Nudge_nlev)
   real(r8) PSanal(Nudge_ncol)
   real(r8) Lat_anal(Nudge_ncol)
   real(r8) Lon_anal(Nudge_ncol)
   integer  nn,Nindex

!DIAG
!  if(masterproc) then
!    write(iulog,*) 'NUDGING: entering nudging_update_analyses_se()'
!  endif
!DIAG
   ! Rotate Nudge_ObsInd() indices, then check the existence of the analyses 
   ! file; broadcast the updated indices and file status to all the other MPI nodes. 
   ! If the file is not there, then just return.
   !------------------------------------------------------------------------
   if(masterproc) then
     Nindex=Nudge_ObsInd(Nudge_NumObs)
     do nn=Nudge_NumObs,2,-1
       Nudge_ObsInd(nn)=Nudge_ObsInd(nn-1)
     end do
     Nudge_ObsInd(1)=Nindex
     inquire(FILE=trim(anal_file),EXIST=Nudge_File_Present(Nudge_ObsInd(1)))
     write(iulog,*)'NUDGING: Nudge_ObsInd=',Nudge_ObsInd
     write(iulog,*)'NUDGING: Nudge_File_Present=',Nudge_File_Present
   endif
#ifdef SPMD
   call mpibcast(Nudge_File_Present, Nudge_NumObs, mpilog, 0, mpicom)
   call mpibcast(Nudge_ObsInd      , Nudge_NumObs, mpiint, 0, mpicom)
#endif
   if(.not.Nudge_File_Present(Nudge_ObsInd(1))) return

   ! masterporc does all of the work here
   !-----------------------------------------
   if(masterproc) then
   
     ! Open the given file
     !-----------------------
     istat=nf90_open(trim(anal_file),NF90_NOWRITE,ncid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*)'NF90_OPEN: failed for file ',trim(anal_file)
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif

     ! Read in Dimensions
     !--------------------
     istat=nf90_inq_dimid(ncid,'ncol',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=ncol)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif

     istat=nf90_inq_dimid(ncid,'lev',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=plev)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif

     istat=nf90_inq_varid(ncid,'lon',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_get_var(ncid,varid,Lon_anal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif

     istat=nf90_inq_varid(ncid,'lat',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_get_var(ncid,varid,Lat_anal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif

     if((Nudge_ncol.ne.ncol).or.(plev.ne.pver)) then
      write(iulog,*) 'ERROR: nudging_update_analyses_se: ncol=',ncol,' Nudge_ncol=',Nudge_ncol
      write(iulog,*) 'ERROR: nudging_update_analyses_se: plev=',plev,' pver=',pver
      call endrun('nudging_update_analyses_se: analyses dimension mismatch')
     endif

     ! Read in and scatter data arrays
     !----------------------------------
     istat=nf90_inq_varid(ncid,'U',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_ncol,Xanal,    &
                               Nobs_U(1,1,begchunk,Nudge_ObsInd(1)))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'V',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_ncol,Xanal,    &
                               Nobs_V(1,1,begchunk,Nudge_ObsInd(1)))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'T',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_ncol,Xanal,    &
                               Nobs_T(1,1,begchunk,Nudge_ObsInd(1)))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'Q',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_ncol,Xanal,    &
                               Nobs_Q(1,1,begchunk,Nudge_ObsInd(1)))

   if(masterproc) then
    istat=nf90_inq_varid(ncid,'PS',varid)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_SE')
    endif
    istat=nf90_get_var(ncid,varid,PSanal)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_SE')
    endif

     ! Close the analyses file
     !-----------------------
     istat=nf90_close(ncid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,1,1,Nudge_ncol,PSanal,           &
                               Nobs_PS(1,begchunk,Nudge_ObsInd(1)))

!DIAG
!  if(masterproc) then
!    write(iulog,*) 'NUDGING: exiting nudging_update_analyses_se()'
!  endif
!DIAG
   ! End Routine
   !------------
   return
  end subroutine ! nudging_update_analyses_se
  !================================================================


  !================================================================
  subroutine nudging_update_analyses_eul(anal_file)
   ! 
   ! NUDGING_UPDATE_ANALYSES_EUL: 
   !                 Open the given analyses data file, read in 
   !                 U,V,T,Q, and PS values and then distribute
   !                 the values to all of the chunks.
   !===============================================================
   use ppgrid ,only: pver,begchunk
   use netcdf

   ! Arguments
   !-------------
   character(len=*),intent(in):: anal_file

   ! Local values
   !-------------
   integer lev
   integer nlon,nlat,plev,istat
   integer ncid,varid
   integer ilat,ilon,ilev
   real(r8) Xanal(Nudge_nlon,Nudge_nlat,Nudge_nlev)
   real(r8) PSanal(Nudge_nlon,Nudge_nlat)
   real(r8) Lat_anal(Nudge_nlat)
   real(r8) Lon_anal(Nudge_nlon)
   real(r8) Xtrans(Nudge_nlon,Nudge_nlev,Nudge_nlat)
   integer  nn,Nindex

   ! Rotate Nudge_ObsInd() indices, then check the existence of the analyses 
   ! file; broadcast the updated indices and file status to all the other MPI nodes. 
   ! If the file is not there, then just return.
   !------------------------------------------------------------------------
   if(masterproc) then
     Nindex=Nudge_ObsInd(Nudge_NumObs)
     do nn=Nudge_NumObs,2,-1
       Nudge_ObsInd(nn)=Nudge_ObsInd(nn-1)
     end do
     Nudge_ObsInd(1)=Nindex
     inquire(FILE=trim(anal_file),EXIST=Nudge_File_Present(Nudge_ObsInd(1)))
   endif
#ifdef SPMD
   call mpibcast(Nudge_File_Present, Nudge_NumObs, mpilog, 0, mpicom)
   call mpibcast(Nudge_ObsInd      , Nudge_NumObs, mpiint, 0, mpicom)
#endif
   if(.not.Nudge_File_Present(Nudge_ObsInd(1))) return

   ! masterporc does all of the work here
   !-----------------------------------------
   if(masterproc) then
   
     ! Open the given file
     !-----------------------
     istat=nf90_open(trim(anal_file),NF90_NOWRITE,ncid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*)'NF90_OPEN: failed for file ',trim(anal_file)
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif

     ! Read in Dimensions
     !--------------------
     istat=nf90_inq_dimid(ncid,'lon',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=nlon)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif

     istat=nf90_inq_dimid(ncid,'lat',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=nlat)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif

     istat=nf90_inq_dimid(ncid,'lev',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=plev)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif

     istat=nf90_inq_varid(ncid,'lon',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     istat=nf90_get_var(ncid,varid,Lon_anal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif

     istat=nf90_inq_varid(ncid,'lat',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     istat=nf90_get_var(ncid,varid,Lat_anal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif

     if((Nudge_nlon.ne.nlon).or.(Nudge_nlat.ne.nlat).or.(plev.ne.pver)) then
      write(iulog,*) 'ERROR: nudging_update_analyses_eul: nlon=',nlon,' Nudge_nlon=',Nudge_nlon
      write(iulog,*) 'ERROR: nudging_update_analyses_eul: nlat=',nlat,' Nudge_nlat=',Nudge_nlat
      write(iulog,*) 'ERROR: nudging_update_analyses_eul: plev=',plev,' pver=',pver
      call endrun('nudging_update_analyses_eul: analyses dimension mismatch')
     endif

     ! Read in, transpose lat/lev indices, 
     ! and scatter data arrays
     !----------------------------------
     istat=nf90_inq_varid(ncid,'U',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     do ilat=1,nlat
     do ilev=1,plev
     do ilon=1,nlon
       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
     end do
     end do
     end do
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_nlon,Xtrans,   &
                               Nobs_U(1,1,begchunk,Nudge_ObsInd(1)))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'V',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     do ilat=1,nlat
     do ilev=1,plev
     do ilon=1,nlon
       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
     end do
     end do
     end do
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_nlon,Xtrans,   &
                               Nobs_V(1,1,begchunk,Nudge_ObsInd(1)))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'T',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     do ilat=1,nlat
     do ilev=1,plev
     do ilon=1,nlon
       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
     end do
     end do
     end do
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_nlon,Xtrans,   &
                               Nobs_T(1,1,begchunk,Nudge_ObsInd(1)))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'Q',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     do ilat=1,nlat
     do ilev=1,plev
     do ilon=1,nlon
       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
     end do
     end do
     end do
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_nlon,Xtrans,   &
                               Nobs_Q(1,1,begchunk,Nudge_ObsInd(1)))

   if(masterproc) then
    istat=nf90_inq_varid(ncid,'PS',varid)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_SE')
    endif
    istat=nf90_get_var(ncid,varid,PSanal)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_SE')
    endif

     ! Close the analyses file
     !-----------------------
     istat=nf90_close(ncid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,1,1,Nudge_nlon,PSanal,           &
                               Nobs_PS(1,begchunk,Nudge_ObsInd(1)))

   ! End Routine
   !------------
   return
  end subroutine ! nudging_update_analyses_eul
  !================================================================


  !================================================================
  subroutine nudging_update_analyses_fv(anal_file)
   ! 
   ! NUDGING_UPDATE_ANALYSES_FV: 
   !                 Open the given analyses data file, read in 
   !                 U,V,T,Q, and PS values and then distribute
   !                 the values to all of the chunks.
   !===============================================================
   use ppgrid ,only: pver,begchunk
   use netcdf

   ! Arguments
   !-------------
   character(len=*),intent(in):: anal_file

   ! Local values
   !-------------
   integer lev
   integer nlon,nlat,plev,istat
   integer ncid,varid
   integer ilat,ilon,ilev
   real(r8) Xanal(Nudge_nlon,Nudge_nlat,Nudge_nlev)
   real(r8) PSanal(Nudge_nlon,Nudge_nlat)
   real(r8) Lat_anal(Nudge_nlat)
   real(r8) Lon_anal(Nudge_nlon)
   real(r8) Xtrans(Nudge_nlon,Nudge_nlev,Nudge_nlat)
   integer  nn,Nindex

   ! Rotate Nudge_ObsInd() indices, then check the existence of the analyses 
   ! file; broadcast the updated indices and file status to all the other MPI nodes. 
   ! If the file is not there, then just return.
   !------------------------------------------------------------------------
   if(masterproc) then
     Nindex=Nudge_ObsInd(Nudge_NumObs)
     do nn=Nudge_NumObs,2,-1
       Nudge_ObsInd(nn)=Nudge_ObsInd(nn-1)
     end do
     Nudge_ObsInd(1)=Nindex
     inquire(FILE=trim(anal_file),EXIST=Nudge_File_Present(Nudge_ObsInd(1)))
   endif
#ifdef SPMD
   call mpibcast(Nudge_File_Present, Nudge_NumObs, mpilog, 0, mpicom)
   call mpibcast(Nudge_ObsInd      , Nudge_NumObs, mpiint, 0, mpicom)
#endif
   if(.not.Nudge_File_Present(Nudge_ObsInd(1))) return

   ! masterporc does all of the work here
   !-----------------------------------------
   if(masterproc) then
   
     ! Open the given file
     !-----------------------
     istat=nf90_open(trim(anal_file),NF90_NOWRITE,ncid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*)'NF90_OPEN: failed for file ',trim(anal_file)
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

     ! Read in Dimensions
     !--------------------
     istat=nf90_inq_dimid(ncid,'lon',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=nlon)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

     istat=nf90_inq_dimid(ncid,'lat',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=nlat)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

     istat=nf90_inq_dimid(ncid,'lev',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=plev)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

     istat=nf90_inq_varid(ncid,'lon',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid,Lon_anal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

     istat=nf90_inq_varid(ncid,'lat',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid,Lat_anal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

     if((Nudge_nlon.ne.nlon).or.(Nudge_nlat.ne.nlat).or.(plev.ne.pver)) then
      write(iulog,*) 'ERROR: nudging_update_analyses_fv: nlon=',nlon,' Nudge_nlon=',Nudge_nlon
      write(iulog,*) 'ERROR: nudging_update_analyses_fv: nlat=',nlat,' Nudge_nlat=',Nudge_nlat
      write(iulog,*) 'ERROR: nudging_update_analyses_fv: plev=',plev,' pver=',pver
      call endrun('nudging_update_analyses_fv: analyses dimension mismatch')
     endif

     ! Read in, transpose lat/lev indices, 
     ! and scatter data arrays
     !----------------------------------
     istat=nf90_inq_varid(ncid,'U',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     do ilat=1,nlat
     do ilev=1,plev
     do ilon=1,nlon
       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
     end do
     end do
     end do
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_nlon,Xtrans,   &
                               Nobs_U(1,1,begchunk,Nudge_ObsInd(1)))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'V',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     do ilat=1,nlat
     do ilev=1,plev
     do ilon=1,nlon
       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
     end do
     end do
     end do
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_nlon,Xtrans,   &
                               Nobs_V(1,1,begchunk,Nudge_ObsInd(1)))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'T',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     do ilat=1,nlat
     do ilev=1,plev
     do ilon=1,nlon
       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
     end do
     end do
     end do
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_nlon,Xtrans,   &
                               Nobs_T(1,1,begchunk,Nudge_ObsInd(1)))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'Q',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     do ilat=1,nlat
     do ilev=1,plev
     do ilon=1,nlon
       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
     end do
     end do
     end do
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_nlon,Xtrans,   &
                               Nobs_Q(1,1,begchunk,Nudge_ObsInd(1)))

   if(masterproc) then
    istat=nf90_inq_varid(ncid,'PS',varid)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_SE')
    endif
    istat=nf90_get_var(ncid,varid,PSanal)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_SE')
    endif

     ! Close the analyses file
     !-----------------------
     istat=nf90_close(ncid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,1,1,Nudge_nlon,PSanal,           &
                               Nobs_PS(1,begchunk,Nudge_ObsInd(1)))

   ! End Routine
   !------------
   return
  end subroutine ! nudging_update_analyses_fv
  !================================================================


  !================================================================
  subroutine nudging_set_profile(rlat,rlon,Nudge_prof,Wprof,nlev)
   ! 
   ! NUDGING_SET_PROFILE: for the given lat,lon, and Nudging_prof, set
   !                      the verical profile of window coeffcients.
   !                      Values range from 0. to 1. to affect spatial
   !                      variations on nudging strength.
   !===============================================================

   ! Arguments
   !--------------
   integer  nlev,Nudge_prof
   real(r8) rlat,rlon
   real(r8) Wprof(nlev)

   ! Local values
   !----------------
   integer  ilev
   real(r8) Hcoef,latx,lonx,Vmax,Vmin
   real(r8) lon_lo,lon_hi,lat_lo,lat_hi,lev_lo,lev_hi

   !---------------
   ! set coeffcient
   !---------------
   if(Nudge_prof.eq.0) then
     ! No Nudging
     !-------------
     Wprof(:)=0.0
   elseif(Nudge_prof.eq.1) then
     ! Uniform Nudging
     !-----------------
     Wprof(:)=1.0
   elseif(Nudge_prof.eq.2) then
     ! Localized Nudging with specified Heaviside window function
     !------------------------------------------------------------
     if(Nudge_Hwin_max.le.Nudge_Hwin_min) then
       ! For a constant Horizontal window function, 
       ! just set Hcoef to the maximum of Hlo/Hhi.
       !--------------------------------------------
       Hcoef=max(Nudge_Hwin_lo,Nudge_Hwin_hi)
     else
       ! get lat/lon relative to window center
       !------------------------------------------
       latx=rlat-Nudge_Hwin_lat0
       lonx=rlon-Nudge_Hwin_lon0
       if(lonx.gt. 180.) lonx=lonx-360.
       if(lonx.le.-180.) lonx=lonx+360.

       ! Calcualte RAW window value
       !-------------------------------
       lon_lo=(Nudge_Hwin_lonWidthH+lonx)/Nudge_Hwin_lonDelta
       lon_hi=(Nudge_Hwin_lonWidthH-lonx)/Nudge_Hwin_lonDelta
       lat_lo=(Nudge_Hwin_latWidthH+latx)/Nudge_Hwin_latDelta
       lat_hi=(Nudge_Hwin_latWidthH-latx)/Nudge_Hwin_latDelta
       Hcoef=((1.+tanh(lon_lo))/2.)*((1.+tanh(lon_hi))/2.) &
            *((1.+tanh(lat_lo))/2.)*((1.+tanh(lat_hi))/2.)

       ! Scale the horizontal window coef for specfied range of values.
       !--------------------------------------------------------
       Hcoef=(Hcoef-Nudge_Hwin_min)/(Nudge_Hwin_max-Nudge_Hwin_min)
       Hcoef=(1.-Hcoef)*Nudge_Hwin_lo + Hcoef*Nudge_Hwin_hi
     endif

     ! Load the RAW vertical window
     !------------------------------
     do ilev=1,nlev
       lev_lo=(float(ilev)-Nudge_Vwin_Lindex)/Nudge_Vwin_Ldelta
       lev_hi=(Nudge_Vwin_Hindex-float(ilev))/Nudge_Vwin_Hdelta
       Wprof(ilev)=((1.+tanh(lev_lo))/2.)*((1.+tanh(lev_hi))/2.)
     end do 

     ! Scale the Window function to span the values between Vlo and Vhi:
     !-----------------------------------------------------------------
     Vmax=maxval(Wprof)
     Vmin=minval(Wprof)
     if(Vmax.le.Vmin) then
       ! For a constant Vertical window function, 
       ! load maximum of Vlo/Vhi into Wprof()
       !--------------------------------------------
       Vmax=max(Nudge_Vwin_lo,Nudge_Vwin_hi)
       Wprof(:)=Vmax
     else
       ! Scale the RAW vertical window for specfied range of values.
       !--------------------------------------------------------
       Wprof(:)=(Wprof(:)-Vmin)/(Vmax-Vmin)
       Wprof(:)=Nudge_Vwin_lo + Wprof(:)*(Nudge_Vwin_hi-Nudge_Vwin_lo)
     endif

     ! The desired result is the product of the vertical profile 
     ! and the horizontal window coeffcient.
     !----------------------------------------------------
     Wprof(:)=Hcoef*Wprof(:)
   else
     call endrun('nudging_set_profile:: Unknown Nudge_prof value')
   endif

   ! End Routine
   !------------
   return
  end subroutine ! nudging_set_profile
  !================================================================


  !================================================================
  real(r8) function nudging_set_PSprofile(rlat,rlon,Nudge_PSprof)
   ! 
   ! NUDGING_SET_PSPROFILE: for the given lat and lon set the surface
   !                      pressure profile value for the specified index.
   !                      Values range from 0. to 1. to affect spatial
   !                      variations on nudging strength.
   !===============================================================

   ! Arguments
   !--------------
   real(r8) rlat,rlon
   integer  Nudge_PSprof

   ! Local values
   !----------------

   !---------------
   ! set coeffcient
   !---------------
   if(Nudge_PSprof.eq.0) then
     ! No Nudging
     !-------------
     nudging_set_PSprofile=0.0
   elseif(Nudge_PSprof.eq.1) then
     ! Uniform Nudging
     !-----------------
     nudging_set_PSprofile=1.0
   else
     call endrun('nudging_set_PSprofile:: Unknown Nudge_prof value')
   endif

   ! End Routine
   !------------
   return
  end function ! nudging_set_PSprofile
  !================================================================


  !================================================================
  subroutine calc_DryStaticEnergy(t, q, phis, ps, dse, ncol)
   ! 
   ! calc_DryStaticEnergy: Given the temperature, specific humidity, surface pressure,
   !                       and surface geopotential for a chunk containing 'ncol' columns, 
   !                       calculate and return the corresponding dry static energy values.
   !--------------------------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid,       only: pver, pverp
   use dycore,       only: dycore_is
   use hycoef,       only: hyai, hybi, ps0, hyam, hybm
   use physconst,    only: zvir, gravit, cpair, rair
   !
   ! Input/Output arguments
   !-----------------------
   integer , intent(in) :: ncol      ! Number of columns in chunk
   real(r8), intent(in) :: t(:,:)    ! (pcols,pver) - temperature
   real(r8), intent(in) :: q(:,:)    ! (pcols,pver) - specific humidity
   real(r8), intent(in) :: ps(:)     ! (pcols)      - surface pressure 
   real(r8), intent(in) :: phis(:)   ! (pcols)      - surface geopotential
   real(r8), intent(out):: dse(:,:)  ! (pcols,pver)  - dry static energy
   !
   ! Local variables
   !------------------
   logical  :: fvdyn                 ! finite volume dynamics
   integer  :: ii,kk                 ! Lon, level, level indices
   real(r8) :: tvfac                 ! Virtual temperature factor
   real(r8) :: hkk(ncol)             ! diagonal element of hydrostatic matrix
   real(r8) :: hkl(ncol)             ! off-diagonal element
   real(r8) :: pint(ncol,pverp)      ! Interface pressures
   real(r8) :: pmid(ncol,pver )      ! Midpoint pressures
   real(r8) :: zi(ncol,pverp)        ! Height above surface at interfaces
   real(r8) :: zm(ncol,pver )        ! Geopotential height at mid level

   ! Set dynamics flag
   !-------------------
   fvdyn = dycore_is ('LR')

   ! Load Pressure values and midpoint pressures 
   !----------------------------------------------
   do kk=1,pverp
     do ii=1,ncol
       pint(ii,kk)=(hyai(kk)*ps0)+(hybi(kk)*ps(ii))
     end do
   end do
   do kk=1,pver
     do ii=1,ncol
       pmid(ii,kk)=(hyam(kk)*ps0)+(hybm(kk)*ps(ii))
     end do
   end do

   ! The surface height is zero by definition.
   !-------------------------------------------
   do ii = 1,ncol
     zi(ii,pverp) = 0.0_r8
   end do

   ! Compute the dry static energy, zi, zm from bottom up
   ! Note, zi(i,k) is the interface above zm(i,k)
   !---------------------------------------------------------
   do kk=pver,1,-1

     ! First set hydrostatic elements consistent with dynamics
     !--------------------------------------------------------
     if(fvdyn) then
       do ii=1,ncol
         hkl(ii)=log(pint(ii,kk+1))-log(pint(ii,kk))
         hkk(ii)=1._r8-(hkl(ii)*pint(ii,kk)/(pint(ii,kk+1)-pint(ii,kk)))
       end do
     else
       do ii=1,ncol
         hkl(ii)=(pint(ii,kk+1)-pint(ii,kk))/pmid(ii,kk)
         hkk(ii)=0.5_r8*hkl(ii)
       end do
     endif

     ! Now compute zm, zi, and dse  (WACCM-X vars rairv/zairv/cpairv not used!)
     !------------------------------------------------------------------------
     do ii=1,ncol
       tvfac=t(ii,kk)*rair*(1._r8+(zvir*q(ii,kk)))/gravit
       zm (ii,kk)=zi(ii,kk+1) + (tvfac*hkk(ii))
       zi (ii,kk)=zi(ii,kk+1) + (tvfac*hkl(ii))
       dse(ii,kk)=(t(ii,kk)*cpair) + phis(ii) + (gravit*zm(ii,kk))
     end do

   end do ! kk=pver,1,-1

   ! End Routine
   !-----------
   return
  end subroutine calc_DryStaticEnergy
  !================================================================

end module nudging

module runtime_opts

!----------------------------------------------------------------------- 
! 
! Purpose: This module is responsible for reading CAM namelist cam_inparm 
!          and broadcasting namelist values if needed.  
! 
! Author:
!   Original routines:  CMS
!   Module:             T. Henderson, September 2003
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!- use statements ------------------------------------------------------
!-----------------------------------------------------------------------
use shr_kind_mod,    only: r8 => shr_kind_r8, SHR_KIND_CL
use spmd_utils,      only: masterproc
use namelist_utils,  only: find_group_name
use pmgrid,          only: plat, plev, plon
use cam_instance,    only: inst_suffix
use cam_history
use cam_control_mod
use cam_diagnostics, only: inithist_all
use cam_logfile,     only: iulog
use pspect
use units
use constituents,    only: pcnst, readtrace
use tracers,         only: tracers_flag
use time_manager,    only: dtime
use filenames,       only: ncdata, bnd_topo, &
                           absems_data, &
                           caseid, &
                           brnch_retain_casename
use dycore,          only: dycore_is
use abortutils,      only: endrun
use rayleigh_friction, only: rayk0, raykrange, raytau0

!-----------------------------------------------------------------------
!- module boilerplate --------------------------------------------------
!-----------------------------------------------------------------------
implicit none
private
save


!-----------------------------------------------------------------------
! Public interfaces ----------------------------------------------------
!-----------------------------------------------------------------------
public read_namelist        ! Set and/or get all runtime options

!-----------------------------------------------------------------------
! Private data ---------------------------------------------------------
!-----------------------------------------------------------------------

character(len=SHR_KIND_CL), private :: nlfilename = 'atm_in' ! Namelist filename

!-----------------------------------------------------------------------
!
! SOMEWHAT ALPHABETICAL listing of variables in the cam_inparm namelist:
!
! variable                description
! --------             -----------------
!
! bnd_topo             Path and filename of topography dataset
! 
! absems_data          Dataset with absorption and emissivity factors.
!
! dtime = nnnn,        Model time step in seconds. Default is dycore dependent.
! 
! fincl1 = 'field1', 'field2',...
!                      List of fields to add to the primary history file.
! fincl1lonlat = 'longitude by latitude','longitude by latitude',...
!                      List of columns ('longitude_latitude') or contiguous 
!                      columns ('longitude:longitude_latitude:latitude') at 
!                      which the fincl1 fields will be output. Individual 
!                      columns are specified as a string using a longitude
!                      degree (greater or equal to 0.) followed by a single 
!                      character (e)ast/(w)est identifer, an
!                      underscore '_' , and a latitude degree followed by a 
!                      single character (n)orth/(s)outh identifier.
!                      example '10e_20n' would pick the model column closest
!                      to 10 degrees east longitude by 20 degrees north 
!                      latitude.  A group of contiguous columns can be 
!                      specified by using lon lat ranges with their single
!                      character east/west or north/south identifiers
!                      example '10e:20e_15n:20n'.  Would outfield all 
!                      fincl1 fields at the model columns which fall
!                      with in the longitude range from 10 east to 20 east
!                      and the latitude range from 15 north to 20 north
!
! fincl[2..6] = 'field1', 'field2',...
!                      List of fields to add to the auxiliary history file.
!
! fincl2..6]lonlat = 'longitude by latitude','longitude by latitude',...
!                      List of columns ('longitude_latitude') or contiguous 
!                      columns ('longitude:longitude_latitude:latitude') at 
!                      which the fincl[2..6] fields will be output. Individual 
!                      columns are specified as a string using a longitude
!                      degree (greater or equal to 0.) followed by a single 
!                      character (e)ast/(w)est identifer, an
!                      underscore '_' , and a latitude degree followed by a 
!                      singel character (n)orth/(s)outh identifier.
!                      example '10e_20n' would pick the model column closest
!                      to 10 degrees east longitude by 20 degrees north 
!                      latitude.  A group of contiguous columns can be 
!                      specified by using lon lat ranges with their single
!                      character east/west or north/south identifiers
!                      example '10e:20e_15n:20n'.  Would outfield all 
!                      fincl[2..6] fields at the model columns which fall
!                      with in the longitude range from 10 east to 20 east
!                      and the latitude range from 15 north to 20 north
!
! fexcl1 = 'field1','field2',... 
!                      List of field names to exclude from default
!                      primary history file (default fields on the 
!                      Master Field List).
! 
! fexcl[2..6] = 'field1','field2',... 
!                      List of field names to exclude from
!                      auxiliary history files.
! 
! lcltod_start = nn,nn,nn,...
!                      Array containing the starting time of day for local time history
!                      averaging. Used in conjuction with lcltod_stop. If lcltod_stop
!                      is less than lcltod_start, then the time range wraps around
!                      24 hours. The start time is included in the interval. Time is
!                      in seconds and defaults to 39600 (11:00 AM).
!                      The first value applies to the primary hist. file,
!                      the second to the first aux. hist. file, etc.
! 
! lcltod_stop = nn,nn,nn,...
!                      Array containing the stopping time of day for local time history
!                      averaging. Used in conjuction with lcltod_start. If lcltod_stop
!                      is less than lcltod_start, then the time range wraps around
!                      24 hours. The stop time is not included in the interval. Time is
!                      in seconds and defaults to 0 (midnight).
!                      The first value applies to the primary hist. file,
!                      the second to the first aux. hist. file, etc.
! 
! lcltod_start = nn,nn,nn,...
!                      Array containing the starting time of day for local time history
!                      averaging. Used in conjuction with lcltod_stop. If lcltod_stop
!                      is less than lcltod_start, then the time range wraps around
!                      24 hours. The start time is included in the interval. Time is
!                      in seconds and defaults to 39600 (11:00 AM).
!                      The first value applies to the primary hist. file,
!                      the second to the first aux. hist. file, etc.
! 
! lcltod_stop = nn,nn,nn,...
!                      Array containing the stopping time of day for local time history
!                      averaging. Used in conjuction with lcltod_start. If lcltod_stop
!                      is less than lcltod_start, then the time range wraps around
!                      24 hours. The stop time is not included in the interval. Time is
!                      in seconds and defaults to 0 (midnight).
!                      The first value applies to the primary hist. file,
!                      the second to the first aux. hist. file, etc.
! 
! mfilt = nn,nn,nn     Array containing the maximum number of time 
!                      samples per disk history file. Defaults to 5.
!                      The first value applies to the primary hist. file,
!                      the second to the first aux. hist. file, etc.
! 
! ncdata               Path and filename of initial condition dataset.
! 
! nhtfrq = nn,nn,nn,.. Output history frequency for each tape
!
!                      If = 0 : monthly average
!                      If > 0 : output every nhtfrq time steps.
!                      If < 0 : output every abs(nhtfrq) hours.
! 
! nlvdry = nn,         Number of layers over which to do dry
!                      adjustment. Defaults to 3.
! 
! cam_branch_file      Filepath of restart file to branch from (nsrest=3)
!                      Full pathname required.
character(len=256) :: cam_branch_file = ' '
!
! use_64bit_nc         True if new 64-bit netCDF formit, false otherwise (default false)
! 

!------------------------------------------------------------------
! The following 3 are specific to Rayleigh friction
! integer rayk0         vertical level at which rayleigh friction term is centered
! real(r8) raykrange    range of rayleigh friction profile; if 0, range is set automatically
! real(r8) raytau0      approximate value of decay time at model top (days);
!                       if 0., no rayleigh friction is applied
!------------------------------------------------------------------
!
!
! hfilename_spec       Flexible filename specifier for history files
!
! 
! pertlim = n.n        Max size of perturbation to apply to initial
!                      temperature field.
!
! phys_alltoall        Dynamics/physics transpose option. See phys_grid module.
!
integer :: phys_alltoall
! 
! phys_loadbalance     Load balance option for performance tuning of 
!                      physics chunks.  See phys_grid module.  
integer :: phys_loadbalance
! 
! phys_twin_algorithm  Load balance option for performance tuning of 
!                      physics chunks.  See phys_grid module.  
integer :: phys_twin_algorithm
! 
! phys_chnk_per_thd    Performance tuning option for physics chunks.  See 
!                      phys_grid module.  
integer :: phys_chnk_per_thd
! 
! tracers_flag = .F.    If true, implement tracer test code. Number of tracers determined
!                      in tracers_suite.F90 must agree with PCNST
!
! readtrace = .T.      If true, tracer initial conditions obtained from 
!                      initial file. 
!
! inithist             Generate initial dataset as auxillary history file
!                      can be set to '6-HOURLY', 'DAILY', 'MONTHLY', 'YEARLY' or 'NONE'. 
!                      default: 'YEARLY'
!
! empty_htapes         true => no fields by default on history tapes
!
! print_step_cost      true => print per timestep cost info
!
! avgflag_pertape      A, I, X, or M means avg, instantaneous, max or min for all fields on
!                      that tape
!
!
!   logical indirect     
!                    ! true => include indirect radiative effects of
!                    ! sulfate aerosols.  Default is false.
!
! inithist_all         .false.:  include only REQUIRED fields on IC file
!                      .true. :  include required AND optional fields on IC file
!                      default:  .false.
!
! met_data_file        name of file that contains the offline meteorology data
! met_data_path        name of directory that contains the offline meteorology data
!
! met_filenames_list   name of file that contains names of the offline 
!                      meteorology data files
!
! met_remove_file      true => the offline meteorology file will be removed
!
! met_cell_wall_winds  true => the offline meteorology winds are defined on the model
!                      grid cell walls
! Physics buffer
logical :: pbuf_global_allocate       ! allocate all buffers as global (default: .true.)


! Conservation checks

logical            :: print_energy_errors ! switch for diagnostic output from check_energy module

! Radiative heating rate calculation options

integer :: iradsw        ! freq. of shortwave radiation calc in time steps (positive)
                         ! or hours (negative).  Default: -1
integer :: iradlw        ! frequency of longwave rad. calc. in time steps (positive)
                         ! or hours (negative).  Default: -1
integer :: iradae        ! frequency of absorp/emis calc in time steps (positive)
                         ! or hours (negative).  Default: -12
integer :: irad_always   ! Specifies length of time in timesteps (positive)
                         ! or hours (negative) SW/LW radiation will be run continuously
                         ! from the start of an initial run.  Default: 0
logical :: spectralflux  ! calculate fluxes (up and down) per band. Default: FALSE

#if (defined WACCM_PHYS)
! iondrag / efield
character(len=256) :: efield_lflux_file
character(len=256) :: efield_hflux_file
character(len=256) :: efield_wei96_file
! waccm qbo data variables
character(len=256) :: qbo_forcing_file
logical            :: qbo_use_forcing
logical            :: qbo_cyclic
#endif

! Upper atmosphere radiative processes (waccm phys)
logical :: nlte_use_mo              ! Determines which constituents are used from NLTE calculations
                                    !  = .true. uses MOZART constituents
                                    !  = .false. uses constituents from bnd dataset cftgcm

! SCM Options
logical  :: single_column
real(r8) :: scmlat,scmlon
integer, parameter :: max_chars = 128
character(len=max_chars) iopfile
character(len=200) :: scm_clubb_iop_name
logical  :: scm_iop_srf_prop
logical  :: scm_relaxation
logical  :: scm_diurnal_avg
logical  :: scm_crm_mode

contains

!=======================================================================

  subroutine read_namelist(single_column_in, scmlon_in, scmlat_in, nlfilename_in )

   !----------------------------------------------------------------------- 
   ! 
   ! Purpose: 
   ! Read data from namelist cam_inparm to define the run. Process some of the
   ! namelist variables to determine history and restart/branch file path 
   ! names.  Check input namelist variables for validity and print them
   ! to standard output. 
   ! 
   ! Method: 
   ! Important Note for running on SUN systems: "implicit automatic (a-z)"
   ! will not work because namelist data must be static.
   !
   ! Author: 
   ! Original version:  CCM1
   ! Standardized:      L. Bath, June 1992
   !                    T. Acker, March 1996
   !     
   !-----------------------------------------------------------------------

   ! Note that the following interfaces are prototypes proposed by Henderson 
   ! and Eaton.  They minimize coupling with other modules.  Design of these 
   ! interfaces should be refined via review by other CAM developers.  
   ! Interface *_defaultopts() gets default values from the responsible 
   ! module (Expert) prior to namelist read.  
   ! Interface *_setopts() sends values to the responsible module (Expert) 
   ! after namelist read.  Erroneous values are handled by Experts.  
   ! TBH  9/8/03 
   !
   use phys_grid,        only: phys_grid_defaultopts, phys_grid_setopts
   
#if (defined WACCM_PHYS)
   use iondrag,          only: iondrag_defaultopts, iondrag_setopts
   use qbo,              only: qbo_defaultopts, qbo_setopts
   use waccm_forcing,    only: waccm_forcing_readnl
#endif

   use chem_surfvals,    only: chem_surfvals_readnl
   use check_energy,     only: check_energy_defaultopts, check_energy_setopts
   use radiation,        only: radiation_defaultopts, radiation_setopts, radiation_printopts
   use cam_restart,      only: restart_defaultopts, restart_setopts, restart_printopts
   use radheat,          only: radheat_defaultopts, radheat_setopts
   use carma_flags_mod,  only: carma_readnl
   use co2_cycle,        only: co2_cycle_readnl
   use shr_string_mod,   only: shr_string_toUpper
   use scamMod,          only: scam_setopts,scam_default_opts

   ! Some modules read their own namelist input.
   use spmd_utils,          only: spmd_utils_readnl
   use physconst,           only: physconst_readnl
   use phys_control,        only: phys_ctl_readnl
   use wv_saturation,       only: wv_sat_readnl
   use ref_pres,            only: ref_pres_readnl
   use cam3_aero_data,      only: cam3_aero_data_readnl
   use cam3_ozone_data,     only: cam3_ozone_data_readnl
   use macrop_driver,       only: macrop_driver_readnl
   use microp_driver,       only: microp_driver_readnl
   use microp_aero,         only: microp_aero_readnl
   use cloud_fraction,      only: cldfrc_readnl
   use cldwat,              only: cldwat_readnl
   use zm_conv,             only: zmconv_readnl
   use hk_conv,             only: hkconv_readnl
   use uwshcu,              only: uwshcu_readnl
   use pkg_cld_sediment,    only: cld_sediment_readnl
   use gw_drag,             only: gw_drag_readnl
   use phys_debug_util,     only: phys_debug_readnl
   use rad_constituents,    only: rad_cnst_readnl
   use radiation_data,      only: rad_data_readnl
   use modal_aer_opt,       only: modal_aer_opt_readnl
   use chemistry,           only: chem_readnl
   use prescribed_volcaero, only: prescribed_volcaero_readnl
   use aerodep_flx,         only: aerodep_flx_readnl
   use solar_data,          only: solar_data_readnl
   use tropopause,          only: tropopause_readnl
   use aoa_tracers,         only: aoa_tracers_readnl
   use prescribed_ozone,    only: prescribed_ozone_readnl
   use prescribed_aero,     only: prescribed_aero_readnl
   use prescribed_ghg,      only: prescribed_ghg_readnl
   use aircraft_emit,       only: aircraft_emit_readnl
   use cospsimulator_intr,  only: cospsimulator_intr_readnl
   use sat_hist,            only: sat_hist_readnl
   use vertical_diffusion,  only: vd_readnl
   use cam_history_support, only: fieldname_len, fieldname_lenp2
   use cam_diagnostics,     only: diag_readnl
   use nudging,             only: nudging_readnl
#if ( defined OFFLINE_DYN )
   use metdata,             only: metdata_readnl
#endif

!---------------------------Arguments-----------------------------------

   logical , intent(in), optional :: single_column_in 
   real(r8), intent(in), optional :: scmlon_in
   real(r8), intent(in), optional :: scmlat_in
   character(len=*)    , optional :: nlfilename_in
!-----------------------------------------------------------------------

   include 'netcdf.inc'

!---------------------------Local variables-----------------------------
   character(len=*), parameter ::  subname = "read_namelist"
! 
   character ctemp*8      ! Temporary character strings
   integer ntspdy         ! number of timesteps per day
   integer t              ! history tape index
   integer lastchar       ! index to last char of a char variable
   integer ierr           ! error code
   integer unitn          ! namelist unit number

   integer f, i
   integer, parameter :: max_chars = 128

   character(len=fieldname_lenp2) fincl1(pflds)
   character(len=fieldname_lenp2) fincl2(pflds)
   character(len=fieldname_lenp2) fincl3(pflds)
   character(len=fieldname_lenp2) fincl4(pflds)
   character(len=fieldname_lenp2) fincl5(pflds)
   character(len=fieldname_lenp2) fincl6(pflds)

   character(len=max_chars) fincl1lonlat(pflds)
   character(len=max_chars) fincl2lonlat(pflds)
   character(len=max_chars) fincl3lonlat(pflds)
   character(len=max_chars) fincl4lonlat(pflds)
   character(len=max_chars) fincl5lonlat(pflds)
   character(len=max_chars) fincl6lonlat(pflds)

   character(len=fieldname_len) fexcl1(pflds)
   character(len=fieldname_len) fexcl2(pflds)
   character(len=fieldname_len) fexcl3(pflds)
   character(len=fieldname_len) fexcl4(pflds)
   character(len=fieldname_len) fexcl5(pflds)
   character(len=fieldname_len) fexcl6(pflds)


   character(len=fieldname_lenp2) fwrtpr1(pflds)
   character(len=fieldname_lenp2) fwrtpr2(pflds)
   character(len=fieldname_lenp2) fwrtpr3(pflds)
   character(len=fieldname_lenp2) fwrtpr4(pflds)
   character(len=fieldname_lenp2) fwrtpr5(pflds)
   character(len=fieldname_lenp2) fwrtpr6(pflds)

!
! Define the cam_inparm namelist
! ***NOTE*** If a namelist option is not described in the CAM Users Guide,
!            it is not supported.  

  namelist /cam_inparm/ ncdata, bnd_topo, &
                    cam_branch_file  ,ndens   ,nhtfrq  , &
                    mfilt   ,absems_data, &
                    lcltod_start, lcltod_stop, &
                    fincl1  ,fincl2  ,fincl3  ,fincl4  ,fincl5  , &
                    fincl1lonlat,fincl2lonlat,fincl3lonlat, &
                    fincl4lonlat  ,fincl5lonlat  , fincl6lonlat , &
                    collect_column_output, &
                    fincl6  ,fexcl1  ,fexcl2  ,fexcl3  ,fexcl4  , &
                    fexcl5  ,fexcl6  ,hfilename_spec, &
                    fwrtpr1 ,fwrtpr2 ,fwrtpr3, fwrtpr4 ,fwrtpr5 ,fwrtpr6 , &
                    dtime, &
                    nlvdry,  &
                    pertlim ,&
                    readtrace, rayk0, raykrange, raytau0, &
                    tracers_flag, &
                    inithist, indirect, &
                    empty_htapes, use_64bit_nc, &
                    print_step_cost, avgflag_pertape, &
                    phys_alltoall, phys_loadbalance, phys_twin_algorithm, &
                    phys_chnk_per_thd, &
                    inithist_all

  ! physics buffer
  namelist /cam_inparm/ pbuf_global_allocate

  ! conservation checks
  namelist /cam_inparm/ print_energy_errors

  ! radiative heating calculation options
  namelist /cam_inparm/ iradsw, iradlw, iradae, irad_always, spectralflux

#if (defined WACCM_PHYS)
  ! iondrag / efield options
  namelist /cam_inparm/ efield_lflux_file, efield_hflux_file, efield_wei96_file
  ! waccm qbo namelist variables
  namelist /cam_inparm/ qbo_use_forcing, qbo_forcing_file, qbo_cyclic
#endif

  ! upper atmosphere radiative processes
  namelist /cam_inparm/ nlte_use_mo

  ! scam
  namelist /cam_inparm/ iopfile,scm_iop_srf_prop,scm_relaxation, &
                        scm_diurnal_avg,scm_crm_mode, scm_clubb_iop_name

! 
!-----------------------------------------------------------------------
  if (present(nlfilename_in)) then
     nlfilename = nlfilename_in
  end if
!
! Determine preset values (this is currently being phased out)
!
   call preset ()
!
! Preset sulfate aerosol related variables

   indirect  = .false.

   ! restart write interval
   call restart_defaultopts( &
      cam_branch_file_out          =cam_branch_file            )

   ! Get default values of runtime options for physics chunking.
   call phys_grid_defaultopts(                      &
      phys_loadbalance_out    =phys_loadbalance,    &
      phys_twin_algorithm_out =phys_twin_algorithm, &
      phys_alltoall_out       =phys_alltoall,       &
      phys_chnk_per_thd_out   =phys_chnk_per_thd    )

   ! conservation
   call check_energy_defaultopts( &
      print_energy_errors_out = print_energy_errors )

   ! radiative heating calcs
   call radiation_defaultopts( &
      iradsw_out      = iradsw,     &
      iradlw_out      = iradlw,     &
      iradae_out      = iradae,     &
      irad_always_out = irad_always, &
      spectralflux_out = spectralflux )

#if (defined WACCM_PHYS)
   ! iondrag / efield
   call iondrag_defaultopts( &
      efield_lflux_file_out =efield_lflux_file, &
      efield_hflux_file_out =efield_hflux_file, &
      efield_wei96_file_out =efield_wei96_file )
   ! qbo forcing
   call qbo_defaultopts( &
      qbo_use_forcing_out  = qbo_use_forcing, &
      qbo_forcing_file_out = qbo_forcing_file,&
      qbo_cyclic_out       = qbo_cyclic       )
#endif

   ! Upper atmosphere radiative processes
   call radheat_defaultopts( nlte_use_mo_out =nlte_use_mo )

   if (present(single_column_in)) then
      call scam_default_opts(scmlat_out=scmlat,scmlon_out=scmlon, &
        single_column_out=single_column, &
        scm_iop_srf_prop_out=scm_iop_srf_prop,&
        scm_relaxation_out=scm_relaxation, &
        scm_diurnal_avg_out=scm_diurnal_avg, &
        scm_crm_mode_out=scm_crm_mode, &
        scm_clubb_iop_name_out=scm_clubb_iop_name)
   end if

   do f = 1, pflds
      fincl1(f) = ' '         
      fincl2(f) = ' '         
      fincl3(f) = ' '         
      fincl4(f) = ' '         
      fincl5(f) = ' '         
      fincl6(f) = ' '         
      fincl1lonlat(f) = ' '
      fincl2lonlat(f) = ' '
      fincl3lonlat(f) = ' '
      fincl4lonlat(f) = ' '
      fincl5lonlat(f) = ' '
      fincl6lonlat(f) = ' '
      fexcl1(f) = ' '
      fexcl2(f) = ' '
      fexcl3(f) = ' '
      fexcl4(f) = ' '
      fexcl5(f) = ' '
      fexcl6(f) = ' '
      fwrtpr1(f) = ' '
      fwrtpr2(f) = ' '
      fwrtpr3(f) = ' '
      fwrtpr4(f) = ' '
      fwrtpr5(f) = ' '
      fwrtpr6(f) = ' '
   enddo

   ! Read in the cam_inparm namelist from input filename

   if (masterproc) then
      write(iulog,*) 'Read in cam_inparm namelist from: ', trim(nlfilename)
      unitn = getunit()
      open( unitn, file=trim(nlfilename), status='old' )

      ! Look for cam_inparm group name in the input file.  If found, leave the
      ! file positioned at that namelist group.
      call find_group_name(unitn, 'cam_inparm', status=ierr)
      if (ierr == 0) then  ! found cam_inparm
         read(unitn, cam_inparm, iostat=ierr)  ! read the cam_inparm namelist group
         if (ierr /= 0) then
            call endrun( subname//':: namelist read returns an'// &
                          ' error condition for cam_inparm' )
         end if
      else
         call endrun(subname // ':: can''t find cam_inparm in file ' // trim(nlfilename))
      end if
      close( unitn )
      call freeunit( unitn )
      !
      ! Check CASE namelist variable
      !
      if (caseid==' ') then
         call endrun ('READ_NAMELIST: Namelist variable CASEID must be set')
      end if

      lastchar = len(caseid)
      if (caseid(lastchar:lastchar) /= ' ') then
         write(iulog,*)'READ_NAMELIST: CASEID must not exceed ', len(caseid)-1, ' characters'
         call endrun
      end if

      do f=1, pflds
         fincl(f, 1) = fincl1(f)
         fincl(f, 2) = fincl2(f)
         fincl(f, 3) = fincl3(f)
         fincl(f, 4) = fincl4(f)
         fincl(f, 5) = fincl5(f)
         fincl(f, 6) = fincl6(f)
         
         fincllonlat(f, 1) = fincl1lonlat(f)
         fincllonlat(f, 2) = fincl2lonlat(f)
         fincllonlat(f, 3) = fincl3lonlat(f)
         fincllonlat(f, 4) = fincl4lonlat(f)
         fincllonlat(f, 5) = fincl5lonlat(f)
         fincllonlat(f, 6) = fincl6lonlat(f)
         if(dycore_is('UNSTRUCTURED') ) then
            do i=1,6
               if (fincllonlat(f,i) .ne. ' ') then
                  call endrun('READ_NAMELIST: Column output is not supported in Unstructered Grids')
               end if
            end do
         end if


         fexcl(f, 1) = fexcl1(f)
         fexcl(f, 2) = fexcl2(f)
         fexcl(f, 3) = fexcl3(f)
         fexcl(f, 4) = fexcl4(f)
         fexcl(f, 5) = fexcl5(f)
         fexcl(f, 6) = fexcl6(f)

         fwrtpr(f, 1) = fwrtpr1(f)
         fwrtpr(f, 2) = fwrtpr2(f)
         fwrtpr(f, 3) = fwrtpr3(f)
         fwrtpr(f, 4) = fwrtpr4(f)
         fwrtpr(f, 5) = fwrtpr5(f)
         fwrtpr(f, 6) = fwrtpr6(f)
      enddo
   end if
!
! Scatter namelist data to all processes
#if ( defined SPMD )
   call distnl ( )
#endif
!
! Auxiliary history files:
! Store input auxf values in array aux (from common block /comhst/).
!
! If generate an initial conditions history file as an auxillary tape:
!
   ctemp = shr_string_toUpper(inithist) 
   inithist = trim(ctemp)
   if (inithist /= '6-HOURLY' .and. inithist /= 'DAILY' .and. &
       inithist /= 'MONTHLY'  .and. inithist /= 'YEARLY' .and. &
       inithist /= 'CAMIOP'   .and. inithist /= 'ENDOFRUN') then
      inithist = 'NONE'
   endif
! 
! History file write up times
! Convert write freq. of hist files from hours to timesteps if necessary.
! 
   do t=1,ptapes
      if (nhtfrq(t) < 0) then
         nhtfrq(t) = nint((-nhtfrq(t)*3600._r8)/dtime)
      end if
   end do
!
! Initialize the filename specifier if not already set
! This is the format for the history filenames:
! %c= caseid, %t=tape no., %y=year, %m=month, %d=day, %s=second, %%=%
! See the filenames module for more information
!
   do t = 1, ptapes
      if ( len_trim(hfilename_spec(t)) == 0 )then
         if ( nhtfrq(t) == 0 )then
            hfilename_spec(t) = '%c.cam' // trim(inst_suffix) // '.h%t.%y-%m.nc'        ! Monthly files
         else
            hfilename_spec(t) = '%c.cam' // trim(inst_suffix) // '.h%t.%y-%m-%d-%s.nc'
         end if
      end if
!
! Only one time sample allowed per monthly average file
! 
      if (nhtfrq(t) == 0) mfilt(t) = 1
   end do

   ! Print per-tape averaging flags
   if (masterproc) then
      do t=1,ptapes
         if (avgflag_pertape(t) /= ' ') then
            write(iulog,*)'Unless overridden by namelist input on a per-field basis (FINCL),'
            write(iulog,*)'All fields on history file ',t,' will have averaging flag ',avgflag_pertape(t)
         end if
      end do
   end if

   ! restart write interval
   call restart_setopts( nsrest,            &
      cam_branch_file_in          =cam_branch_file            )


   ! Set runtime options for physics chunking.
   call phys_grid_setopts(                          &
       phys_loadbalance_in    =phys_loadbalance,    &
       phys_twin_algorithm_in =phys_twin_algorithm, &
       phys_alltoall_in       =phys_alltoall,       &
       phys_chnk_per_thd_in   =phys_chnk_per_thd    )

   ! conservation
   call check_energy_setopts( &
      print_energy_errors_in = print_energy_errors )

   call radiation_setopts( dtime, nhtfrq(1), &
      iradsw_in      = iradsw,     &
      iradlw_in      = iradlw,     &
      iradae_in      = iradae,     &
      irad_always_in = irad_always, &
      spectralflux_in = spectralflux )

#if (defined WACCM_PHYS)
   ! iondrag / efield
   call iondrag_setopts( &
        efield_lflux_file_in =efield_lflux_file, &
        efield_hflux_file_in =efield_hflux_file, &
        efield_wei96_file_in =efield_wei96_file)
   ! qbo forcing
   call qbo_setopts( &
        qbo_use_forcing_in  = qbo_use_forcing, &
        qbo_forcing_file_in = qbo_forcing_file,&
        qbo_cyclic_in       = qbo_cyclic       )
#endif

   ! Upper atmosphere radiative processes
   call radheat_setopts( nlte_use_mo_in =nlte_use_mo )
! 
! Set runtime options for single column mode
!
   if (present(single_column_in) .and. present(scmlon_in) .and. present(scmlat_in)) then 
      if (single_column_in) then
         single_column = single_column_in
         scmlon = scmlon_in
         scmlat = scmlat_in
         call scam_setopts( scmlat_in=scmlat,scmlon_in=scmlon, &
                            iopfile_in=iopfile,single_column_in=single_column,&
                            scm_iop_srf_prop_in=scm_iop_srf_prop,&
                            scm_relaxation_in=scm_relaxation, &
                            scm_diurnal_avg_in=scm_diurnal_avg, &
                            scm_crm_mode_in=scm_crm_mode, &
                            scm_clubb_iop_name_in=scm_clubb_iop_name)
      end if
   endif

   ! Call subroutines for modules to read their own namelist.
   ! In some cases namelist default values may depend on settings from
   ! other modules, so there may be an order dependence in the following
   ! calls.
   ! ***N.B.*** In particular, physconst_readnl should be called before
   !            the other readnl methods in case that method is used to set
   !            physical constants, some of which are set at runtime
   !            by the physconst_readnl method.
   ! Modules that read their own namelist are responsible for making sure
   ! all processes receive the values.

   call spmd_utils_readnl(nlfilename)
   call physconst_readnl(nlfilename)
   call chem_surfvals_readnl(nlfilename)
   call phys_ctl_readnl(nlfilename)
   call wv_sat_readnl(nlfilename)
   call ref_pres_readnl(nlfilename)
   call cam3_aero_data_readnl(nlfilename)
   call cam3_ozone_data_readnl(nlfilename)
   call macrop_driver_readnl(nlfilename)
   call microp_driver_readnl(nlfilename)
   call microp_aero_readnl(nlfilename)
   call cldfrc_readnl(nlfilename)
   call zmconv_readnl(nlfilename)
   call cldwat_readnl(nlfilename)
   call hkconv_readnl(nlfilename)
   call uwshcu_readnl(nlfilename)
   call cld_sediment_readnl(nlfilename)
   call gw_drag_readnl(nlfilename)
   call phys_debug_readnl(nlfilename)
   call rad_cnst_readnl(nlfilename)
   call rad_data_readnl(nlfilename)
   call modal_aer_opt_readnl(nlfilename)
   call chem_readnl(nlfilename)
   call prescribed_volcaero_readnl(nlfilename)
   call solar_data_readnl(nlfilename)
   call carma_readnl(nlfilename)
   call tropopause_readnl(nlfilename)
   call aoa_tracers_readnl(nlfilename)
   call aerodep_flx_readnl(nlfilename)
   call prescribed_ozone_readnl(nlfilename)
   call prescribed_aero_readnl(nlfilename)
   call prescribed_ghg_readnl(nlfilename)
   call co2_cycle_readnl(nlfilename)
   call aircraft_emit_readnl(nlfilename)
   call cospsimulator_intr_readnl(nlfilename)
   call sat_hist_readnl(nlfilename, hfilename_spec, mfilt, fincl, nhtfrq, avgflag_pertape)
   call diag_readnl(nlfilename)
   call nudging_readnl(nlfilename)
#if (defined WACCM_PHYS)
   call waccm_forcing_readnl(nlfilename)
#endif
   call vd_readnl(nlfilename)
#if ( defined OFFLINE_DYN )
   call metdata_readnl(nlfilename)
#endif

! 
! Print cam_inparm input variables to standard output
! 
   if (masterproc) then
      write(iulog,*)' ------------------------------------------'
      write(iulog,*)'     *** INPUT VARIABLES (CAM_INPARM) ***'
      write(iulog,*)' ------------------------------------------'
      if (nsrest/=0) then
         write(iulog,*) '  Continuation of an earlier run'
      else
         write(iulog,*) '         Initial run'
      end if
      write(iulog,*) ' ********** CASE = ',trim(caseid),' **********'
      write(iulog,'(1x,a)') ctitle
      if (len_trim(ncdata) > 0) then
         write(iulog,*) 'Initial dataset is: ',trim(ncdata)
      end if
      write(iulog,*)'Topography dataset is: ', trim(bnd_topo)
      write(iulog,*)'Time-invariant (absorption/emissivity) factor dataset is: ', trim(absems_data)

      ! Type of run
      write(iulog,*)'Run type flag (NSREST) 0=initial, 1=restart, 3=branch ',nsrest

      call restart_printopts()

   end if
!
! History file info 
!
   if (masterproc) then
      if (inithist == '6-HOURLY' ) then
         write(iulog,*)'Initial conditions history files will be written 6-hourly.'
      else if (inithist == 'DAILY' ) then
         write(iulog,*)'Initial conditions history files will be written daily.'
      else if (inithist == 'MONTHLY' ) then
         write(iulog,*)'Initial conditions history files will be written monthly.'
      else if (inithist == 'YEARLY' ) then
         write(iulog,*)'Initial conditions history files will be written yearly.'
      else if (inithist == 'CAMIOP' ) then
         write(iulog,*)'Initial conditions history files will be written for IOP.'
      else if (inithist == 'ENDOFRUN' ) then
         write(iulog,*)'Initial conditions history files will be written at end of run.'
      else
         write(iulog,*)'Initial conditions history files will not be created'
      end if

!
! Write physics variables from namelist cam_inparm to std. output
!
      write(iulog,9108) nlvdry
9108 format('Lowest level for dry adiabatic adjust (NLVDRY)',i10)


      call radiation_printopts()

      if ( (adiabatic .and. ideal_phys) .or. (adiabatic .and. aqua_planet) .or. &
           (ideal_phys .and. aqua_planet) ) then
         call endrun ('READ_NAMELIST: Only one of ADIABATIC, IDEAL_PHYS, or AQUA_PLANET can be .true.')
      end if

#ifdef COUP_SOM
      if (adiabatic .or. ideal_phys .or. aqua_planet )then
         call endrun ('READ_NAMELIST: adiabatic, ideal_phys or aqua_planet can not be used with SOM')
      end if
#else
      if (adiabatic)   write(iulog,*) 'Model will run ADIABATICALLY (i.e. no physics)'
      if (ideal_phys)  write(iulog,*) 'Run ONLY the "idealized" dynamical core of the ', &
                                  'model  (dynamics + Held&Suarez-specified physics)'
      if (aqua_planet) write(iulog,*) 'Run model in "AQUA_PLANET" mode'
#endif
   end if

   ! set public data in cam_control_mod
   moist_physics = (.not. adiabatic) .and. (.not. ideal_phys)

#ifdef PERGRO
   if (masterproc) then
      write(iulog,*)'pergro for cloud water is true'
   end if
#endif

   ntspdy = nint(86400._r8/dtime) ! no. timesteps per day


end subroutine read_namelist


!=======================================================================

#ifdef SPMD
subroutine distnl
!-----------------------------------------------------------------------
!     
! Purpose:     
! Distribute namelist data to all processors.
!
! The cpp SPMD definition provides for the funnelling of all program i/o
! through the master processor. Processor 0 either reads restart/history
! data from the disk and distributes it to all processors, or collects
! data from all processors and writes it to disk.
!     
!---------------------------Code history-------------------------------
!
! Original version:  CCM2
! Standardized:      J. Rosinski, Oct 1995
!                    J. Truesdale, Feb. 1996
!
!-----------------------------------------------------------------------
   use mpishorthand
!-----------------------------------------------------------------------

!
!-----------------------------------------------------------------------
! 
   call mpibcast (dtime,       1,mpiint,0,mpicom)
   call mpibcast (ndens   ,ptapes,mpiint,0,mpicom)
   call mpibcast (nhtfrq  ,ptapes,mpiint,0,mpicom)
   call mpibcast (mfilt   ,ptapes,mpiint,0,mpicom)
   call mpibcast (lcltod_start ,ptapes,mpiint,0,mpicom)
   call mpibcast (lcltod_stop  ,ptapes,mpiint,0,mpicom)
   call mpibcast (nsrest  ,1,mpiint,0,mpicom)
   call mpibcast (nlvdry  ,1,mpiint,0,mpicom)

   call mpibcast (rayk0    ,1,mpiint,0,mpicom)
   call mpibcast (raykrange,1,mpir8,0,mpicom)
   call mpibcast (raytau0  ,1,mpir8,0,mpicom)

   call mpibcast (collect_column_output,ptapes,mpilog,0,mpicom)

   call mpibcast (tracers_flag,1,mpilog,0,mpicom)
   call mpibcast (readtrace   ,1,mpilog,0,mpicom)
   call mpibcast (adiabatic   ,1,mpilog,0,mpicom)
   call mpibcast (ideal_phys  ,1,mpilog,0,mpicom)
   call mpibcast (aqua_planet ,1,mpilog,0,mpicom)

   call mpibcast (empty_htapes,1,mpilog,0,mpicom)
   call mpibcast (use_64bit_nc,1,mpilog,0,mpicom)
   call mpibcast (print_step_cost,1,mpilog,0,mpicom)
   call mpibcast (inithist_all   ,1,mpilog,0,mpicom)
   call mpibcast (pertlim     ,1, mpir8,  0, mpicom )

   call mpibcast (caseid  ,len(caseid) ,mpichar,0,mpicom)
   call mpibcast (avgflag_pertape, ptapes, mpichar,0,mpicom)
   call mpibcast (ctitle  ,len(ctitle),mpichar,0,mpicom)
   call mpibcast (ncdata  ,len(ncdata) ,mpichar,0,mpicom)
   call mpibcast (bnd_topo  ,len(bnd_topo) ,mpichar,0,mpicom)
   call mpibcast (absems_data,len(absems_data),mpichar,0,mpicom)
   call mpibcast (cam_branch_file  ,len(cam_branch_file) ,mpichar,0,mpicom)
   call mpibcast (inithist,len(inithist)  ,mpichar,0,mpicom)
   call mpibcast (hfilename_spec, len(hfilename_spec(1))*ptapes, mpichar, 0, mpicom)
   call mpibcast (fincl   ,len(fincl (1,1))*pflds*ptapes,mpichar,0,mpicom)
   call mpibcast (fexcl   ,len(fexcl (1,1))*pflds*ptapes,mpichar,0,mpicom)

   call mpibcast (fincllonlat   ,len(fincllonlat (1,1))*pflds*ptapes,mpichar,0,mpicom)

   call mpibcast (fwrtpr  ,len(fwrtpr(1,1))*pflds*ptapes,mpichar,0,mpicom)

   call mpibcast (indirect     , 1 ,mpilog, 0,mpicom)

   ! Physics chunk tuning
   call mpibcast (phys_loadbalance   ,1,mpiint,0,mpicom)
   call mpibcast (phys_twin_algorithm,1,mpiint,0,mpicom)
   call mpibcast (phys_alltoall      ,1,mpiint,0,mpicom)
   call mpibcast (phys_chnk_per_thd  ,1,mpiint,0,mpicom)

   ! Physics buffer
   call mpibcast (pbuf_global_allocate, 1, mpilog, 0, mpicom)

   ! Conservation
   call mpibcast (print_energy_errors, 1, mpilog, 0, mpicom)

   ! Radiative heating calculation
   call mpibcast (iradsw,     1, mpiint, 0, mpicom)
   call mpibcast (iradlw,     1, mpiint, 0, mpicom)
   call mpibcast (iradae,     1, mpiint, 0, mpicom)
   call mpibcast (irad_always,1, mpiint, 0, mpicom)
   call mpibcast (spectralflux,1, mpilog, 0, mpicom)

#if (defined WACCM_PHYS)
   ! iondrag / efield options
   call mpibcast (efield_lflux_file, len(efield_lflux_file), mpichar, 0, mpicom)
   call mpibcast (efield_hflux_file, len(efield_hflux_file), mpichar, 0, mpicom)
   call mpibcast (efield_wei96_file, len(efield_wei96_file), mpichar, 0, mpicom)
   ! qbo variables
   call mpibcast (qbo_forcing_file,  len(qbo_forcing_file ), mpichar, 0, mpicom)
   call mpibcast (qbo_use_forcing,   1,                      mpilog,  0, mpicom)
   call mpibcast (qbo_cyclic,        1,                      mpilog,  0, mpicom)
#endif

   call mpibcast (nlte_use_mo,            1,  mpilog, 0, mpicom)

end subroutine distnl
#endif



subroutine preset
!----------------------------------------------------------------------- 
! 
! Purpose: Preset namelist CAM_INPARM input variables and initialize some other variables
! 
! Method: Hardwire the values
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
   use cam_history,  only: fincl, fexcl, fwrtpr, fincllonlat, collect_column_output
   use rgrid
!-----------------------------------------------------------------------
   include 'netcdf.inc'
!-----------------------------------------------------------------------
!
! Preset character history variables here because module initialization of character arrays
! does not work on all machines
! $$$ TBH:  is this still true?  12/14/03
!
   fincl(:,:)  = ' '
   fincllonlat(:,:)  = ' '
   fexcl(:,:)  = ' '
   fwrtpr(:,:) = ' '
!
! Flags
!
   print_step_cost = .false.   ! print per timestep cost info
   collect_column_output = .false.
!
! rgrid: set default to full grid
!
   nlon(:) = plon
!!
!! Unit numbers: set to invalid
!!
!   ncid_ini = -1
!   ncid_sst = -1
!   ncid_trc = -1
!
   return
end subroutine preset

end module runtime_opts

module physpkg
  !-----------------------------------------------------------------------
  ! Purpose:
  !
  ! Provides the interface to CAM physics package
  !
  ! Revision history:
  ! Aug  2005,  E. B. Kluzek,  Creation of module from physpkg subroutine
  ! 2005-10-17  B. Eaton       Add contents of inti.F90 to phys_init().  Add
  !                            initialization of grid info in phys_state.
  ! Nov 2010    A. Gettelman   Put micro/macro physics into separate routines
  ! Oct 2014    Yi-Chi Wang    revise this scheme to put a flag (ideep) as
  !                            communication for deep and shallow scheme.
  !                            After CAM5.2, this file combines previous tphysbc.
  ! Jun 2015    Yi-Chi Wang    Modify the physpkg with nudging code v1.0 from Patrick of NCAR.
  !-----------------------------------------------------------------------

  use shr_kind_mod,     only: r8 => shr_kind_r8
  use spmd_utils,       only: masterproc
  use physconst,        only: latvap, latice, rh2o
  use physics_types,    only: physics_state, physics_tend, physics_state_set_grid, &
       physics_ptend, physics_tend_init, physics_update,    &
       physics_type_alloc, physics_ptend_dealloc
  use phys_grid,        only: get_ncols_p
  use phys_gmean,       only: gmean_mass
  use ppgrid,           only: begchunk, endchunk, pcols, pver, pverp
  use constituents,     only: pcnst, cnst_name, cnst_get_ind
  use camsrfexch,       only: cam_out_t, cam_in_t

  use cam_control_mod,  only: ideal_phys, adiabatic
  use phys_control,     only: phys_do_flux_avg, waccmx_is
  use scamMod,          only: single_column, scm_crm_mode
  use flux_avg,         only: flux_avg_init
  use infnan,           only: posinf, assignment(=)
#ifdef SPMD
  use mpishorthand
#endif
  use perf_mod
  use cam_logfile,     only: iulog
  use camsrfexch,      only: cam_export
  use phys_control,    only: do_waccm_phys

  implicit none
  private

  !  Physics buffer index
  integer ::  teout_idx          = 0  

  integer ::  tini_idx           = 0 
  integer ::  qini_idx           = 0 
  integer ::  cldliqini_idx      = 0 
  integer ::  cldiceini_idx      = 0 

  integer ::  prec_str_idx       = 0
  integer ::  snow_str_idx       = 0
  integer ::  prec_sed_idx       = 0
  integer ::  snow_sed_idx       = 0
  integer ::  prec_pcw_idx       = 0
  integer ::  snow_pcw_idx       = 0
  integer ::  prec_dp_idx        = 0
  integer ::  snow_dp_idx        = 0
  integer ::  prec_sh_idx        = 0
  integer ::  snow_sh_idx        = 0


  save

  ! Public methods
  public phys_register ! was initindx  - register physics methods
  public phys_init   ! Public initialization method
  public phys_run1   ! First phase of the public run method
  public phys_run2   ! Second phase of the public run method
  public phys_final  ! Public finalization method
  !
  ! Private module data
  !

  !======================================================================= 
contains


subroutine phys_register
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Register constituents and physics buffer fields.
    ! 
    ! Author:    CSM Contact: M. Vertenstein, Aug. 1997
    !            B.A. Boville, Oct 2001
    !            A. Gettelman, Nov 2010 - put micro/macro physics into separate routines
    ! 
    !-----------------------------------------------------------------------
    use physics_buffer, only: pbuf_init_time
    use shr_kind_mod,       only: r8 => shr_kind_r8
    use spmd_utils,         only: masterproc
    use constituents,       only: pcnst, cnst_add, cnst_chk_dim, cnst_name

    use cam_control_mod,    only: moist_physics
    use phys_control,       only: phys_do_flux_avg, phys_getopts, waccmx_is
    use chemistry,          only: chem_register
    use cloud_fraction,     only: cldfrc_register
    use stratiform,         only: stratiform_register
    use microp_driver,      only: microp_driver_register
    use microp_aero,        only: microp_aero_register
    use macrop_driver,      only: macrop_driver_register
    use clubb_intr,         only: clubb_register_cam
    use conv_water,         only: conv_water_register
    use physconst,          only: mwdry, cpair, mwh2o, cpwv
    use tracers,            only: tracers_register
    use check_energy,       only: check_energy_register
    use aerosol_intr,       only: aerosol_register_cnst
    use carma_intr,         only: carma_register
    use cam3_aero_data,     only: cam3_aero_data_on, cam3_aero_data_register
    use cam3_ozone_data,    only: cam3_ozone_data_on, cam3_ozone_data_register
    use ghg_data,           only: ghg_data_register
    use vertical_diffusion, only: vd_register
    use convect_deep,       only: convect_deep_register
    use convect_shallow,    only: convect_shallow_register
    use radiation,          only: radiation_register
    use co2_cycle,          only: co2_register
    use flux_avg,           only: flux_avg_register
    use exbdrift,           only: exbdrift_register
    use gw_drag,            only: gw_drag_register
    use iondrag,            only: iondrag_register
    use ionosphere,         only: ionos_register
    use string_utils,       only: to_lower
    use prescribed_ozone,   only: prescribed_ozone_register
    use prescribed_volcaero,only: prescribed_volcaero_register
    use prescribed_aero,    only: prescribed_aero_register
    use prescribed_ghg,     only: prescribed_ghg_register
    use sslt_rebin,         only: sslt_rebin_register
    use aoa_tracers,        only: aoa_tracers_register
    use aircraft_emit,      only: aircraft_emit_register
    use cam_diagnostics,    only: diag_register
    use cloud_diagnostics,  only: cloud_diagnostics_register
    use physics_buffer,     only: pbuf_add_field, dtype_r8

    implicit none
    !---------------------------Local variables-----------------------------
    !
    integer  :: m        ! loop index
    integer  :: mm       ! constituent index 
    !-----------------------------------------------------------------------

    character(len=16) :: microp_scheme
    logical           :: do_clubb_sgs

    call phys_getopts( microp_scheme_out = microp_scheme )
    call phys_getopts( do_clubb_sgs_out  = do_clubb_sgs )

    ! Initialize pbuf_times
    call pbuf_init_time()

    ! Register water vapor.
    ! ***** N.B. ***** This must be the first call to cnst_add so that
    !                  water vapor is constituent 1.
    if (moist_physics) then
       call cnst_add('Q', mwh2o, cpwv, 1.E-12_r8, mm, &
            longname='Specific humidity', readiv=.true., is_convtran1=.true.)
    else
       call cnst_add('Q', mwh2o, cpwv, 0.0_r8, mm, &
            longname='Specific humidity', readiv=.false., is_convtran1=.true.)
    end if

    ! Fields for physics package diagnostics
    call pbuf_add_field('TINI',      'physpkg', dtype_r8, (/pcols,pver/), tini_idx)
    call pbuf_add_field('QINI',      'physpkg', dtype_r8, (/pcols,pver/), qini_idx)
    call pbuf_add_field('CLDLIQINI', 'physpkg', dtype_r8, (/pcols,pver/), cldliqini_idx)
    call pbuf_add_field('CLDICEINI', 'physpkg', dtype_r8, (/pcols,pver/), cldiceini_idx)

    ! check energy package
    call check_energy_register

    ! If using an ideal/adiabatic physics option, the CAM physics parameterizations 
    ! aren't called.
    if (moist_physics) then

       ! register fluxes for saving across time
       if (phys_do_flux_avg()) call flux_avg_register()

       call cldfrc_register()

       ! cloud water
       if( microp_scheme == 'RK' ) then
          call stratiform_register()
       elseif( microp_scheme == 'MG' ) then
          if (.not. do_clubb_sgs) call macrop_driver_register()
          call microp_aero_register()
          call microp_driver_register()
       end if
       
       ! Register CLUBB_SGS here
       if (do_clubb_sgs) call clubb_register_cam()
       

       call pbuf_add_field('PREC_STR',  'physpkg',dtype_r8,(/pcols/),prec_str_idx)
       call pbuf_add_field('SNOW_STR',  'physpkg',dtype_r8,(/pcols/),snow_str_idx)
       call pbuf_add_field('PREC_PCW',  'physpkg',dtype_r8,(/pcols/),prec_pcw_idx)
       call pbuf_add_field('SNOW_PCW',  'physpkg',dtype_r8,(/pcols/),snow_pcw_idx)
       call pbuf_add_field('PREC_SED',  'physpkg',dtype_r8,(/pcols/),prec_sed_idx)
       call pbuf_add_field('SNOW_SED',  'physpkg',dtype_r8,(/pcols/),snow_sed_idx)

       call conv_water_register()

       ! chemical constituents
       call chem_register()

       ! co2 constituents
       call co2_register()

       ! register data model ozone with pbuf
       if (cam3_ozone_data_on) then
          call cam3_ozone_data_register()
       end if
       call prescribed_volcaero_register()
       call prescribed_ozone_register()
       call prescribed_aero_register()
       call prescribed_ghg_register()
       call sslt_rebin_register

       ! CAM3 prescribed aerosols
       if (cam3_aero_data_on) then
          call cam3_aero_data_register()
       end if

       ! register various data model gasses with pbuf
       call ghg_data_register()

       ! Initialize e and b fields
       if (do_waccm_phys()) call exbdrift_register()

       ! waccm gravity wave drag
       call gw_drag_register()

       ! carma microphysics
       ! 
       ! NOTE: Needs to come before aerosol_register_cnst, so that the CARMA
       ! flags are defined by then.
       call carma_register()

       if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then
          ! Register iondrag variables with pbuf
          call iondrag_register()
          ! Register ionosphere variables with pbuf if mode set to ionosphere
          if( waccmx_is('ionosphere') ) then
             call ionos_register()
          endif
       endif

       ! aerosols
       call aerosol_register_cnst()

       call aircraft_emit_register()

       ! deep convection
       call convect_deep_register

       !  shallow convection
       call convect_shallow_register

       ! radiation
       call radiation_register
       call cloud_diagnostics_register

       ! vertical diffusion
       if (.not. do_clubb_sgs) call vd_register()
    end if

    ! Register diagnostics PBUF
    call diag_register()

    ! Register age of air tracers
    call aoa_tracers_register()

    ! Register test tracers
    ! ***** N.B. ***** This is the last call to register constituents because
    !                  the test tracers fill the remaining available slots up
    !                  to constituent number PCNST -- regardless of what PCNST is set to.
    call tracers_register()

    ! All tracers registered, check that the dimensions are correct
    call cnst_chk_dim()

    ! ***NOTE*** No registering constituents after the call to cnst_chk_dim.

end subroutine phys_register



  !======================================================================= 

subroutine phys_inidat( cam_out, pbuf2d )
    use abortutils, only : endrun

    use physics_buffer, only : pbuf_get_index, pbuf_get_field, physics_buffer_desc, pbuf_set_field, pbuf_times


    use cam_initfiles,       only: initial_file_get_id, topo_file_get_id
    use pio,                 only: file_desc_t
    use ncdio_atm,           only: infld
    use dycore,              only: dycore_is
    use polar_avg,           only: polar_average
    use short_lived_species, only: initialize_short_lived_species
    use comsrf,              only: landm, sgh, sgh30
    use cam_control_mod,     only: aqua_planet

    type(cam_out_t),     intent(inout) :: cam_out(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    integer :: lchnk, m, n, i, k, ncol
    type(file_desc_t), pointer :: fh_ini, fh_topo
    character(len=8) :: fieldname
    real(r8), pointer :: cldptr(:,:,:,:), convptr_3d(:,:,:,:)
    real(r8), pointer :: tptr(:,:), tptr3d(:,:,:), tptr3d_2(:,:,:)
    real(r8), pointer :: qpert(:,:)

    character*11 :: subname='phys_inidat' ! subroutine name
    integer :: tpert_idx, qpert_idx, pblh_idx

    logical :: found=.false., found2=.false.
    integer :: ierr
    character(len=4) :: dim1name
    integer :: ixcldice, ixcldliq
    nullify(tptr,tptr3d,tptr3d_2,cldptr,convptr_3d)

    fh_ini=>initial_file_get_id()

    !   dynamics variables are handled in dyn_init - here we read variables needed for physics 
    !   but not dynamics

    if(dycore_is('UNSTRUCTURED')) then  
       dim1name='ncol'
    else
       dim1name='lon'
    end if
    if(aqua_planet) then
       sgh = 0._r8
       sgh30 = 0._r8
       landm = 0._r8
    else
       fh_topo=>topo_file_get_id()
       call infld('SGH', fh_topo, dim1name, 'lat', 1, pcols, begchunk, endchunk, &
            sgh, found, grid_map='PHYS')
       if(.not. found) call endrun('ERROR: SGH not found on topo file')

       call infld('SGH30', fh_topo, dim1name, 'lat', 1, pcols, begchunk, endchunk, &
            sgh30, found, grid_map='PHYS')
       if(.not. found) then
          if (masterproc) write(iulog,*) 'Warning: Error reading SGH30 from topo file.'
          if (masterproc) write(iulog,*) 'The field SGH30 will be filled using data from SGH.'
          sgh30 = sgh
       end if

       call infld('LANDM_COSLAT', fh_topo, dim1name, 'lat', 1, pcols, begchunk, endchunk, &
            landm, found, grid_map='PHYS')

       if(.not.found) call endrun(' ERROR: LANDM_COSLAT not found on topo dataset.')
    end if

    allocate(tptr(1:pcols,begchunk:endchunk))

    call infld('PBLH', fh_ini, dim1name, 'lat', 1, pcols, begchunk, endchunk, &
         tptr(:,:), found, grid_map='PHYS')
    if(.not. found) then
       tptr(:,:) = 0._r8
       if (masterproc) write(iulog,*) 'PBLH initialized to 0.'
    end if
    pblh_idx = pbuf_get_index('pblh')

    call pbuf_set_field(pbuf2d, pblh_idx, tptr)

    call infld('TPERT', fh_ini, dim1name, 'lat', 1, pcols, begchunk, endchunk, &
         tptr(:,:), found, grid_map='PHYS')
    if(.not. found) then
       tptr(:,:) = 0._r8
       if (masterproc) write(iulog,*) 'TPERT initialized to 0.'
    end if
    tpert_idx = pbuf_get_index( 'tpert')
    call pbuf_set_field(pbuf2d, tpert_idx, tptr)

    fieldname='QPERT'  
    qpert_idx = pbuf_get_index( 'qpert',ierr)
    if (qpert_idx > 0) then
       call infld(fieldname, fh_ini, dim1name, 'lat', 1, pcols, begchunk, endchunk, &
            tptr, found, grid_map='PHYS')
       if(.not. found) then
          tptr=0_r8
          if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.'
       end if

       allocate(tptr3d_2(pcols,pcnst,begchunk:endchunk))
       tptr3d_2 = 0_r8
       tptr3d_2(:,1,:) = tptr(:,:)

       call pbuf_set_field(pbuf2d, qpert_idx, tptr3d_2)
       deallocate(tptr3d_2)
    end if

    fieldname='CUSH'
    m = pbuf_get_index('cush')
    call infld(fieldname, fh_ini, dim1name, 'lat', 1, pcols, begchunk, endchunk, &
         tptr, found, grid_map='PHYS')
    if(.not.found) then
       if(masterproc) write(iulog,*) trim(fieldname), ' initialized to 1000.'
       tptr=1000._r8
    end if
    do n=1,pbuf_times
       call pbuf_set_field(pbuf2d, m, tptr, start=(/1,n/), kount=(/pcols,1/))
    end do
    deallocate(tptr)

    do lchnk=begchunk,endchunk
       cam_out(lchnk)%tbot(:) = posinf
    end do

    !
    ! 3-D fields
    !

    allocate(tptr3d(pcols,pver,begchunk:endchunk))

    fieldname='CLOUD'
    m = pbuf_get_index('CLD')
    call infld(fieldname, fh_ini, dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
         tptr3d, found, grid_map='PHYS')
    if(found) then
       do n = 1, pbuf_times
          call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
       end do
    else
       call pbuf_set_field(pbuf2d, m, 0._r8)
       if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.'
    end if

    fieldname='QCWAT'
    m = pbuf_get_index(fieldname,ierr)
    if (m > 0) then
       call infld(fieldname, fh_ini, dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
            tptr3d, found, grid_map='PHYS')
       if(.not. found) then
          call infld('Q',fh_ini,dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
               tptr3d, found, grid_map='PHYS')
          if (found) then
             if (masterproc) write(iulog,*) trim(fieldname), ' initialized with Q'
             if(dycore_is('LR')) call polar_average(pver, tptr3d) 	
          else
             call endrun('  '//trim(subname)//' Error:  Q must be on Initial File')
          end if
       end if
       do n = 1, pbuf_times
          call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
       end do
    end if

    fieldname = 'ICCWAT'
    m = pbuf_get_index(fieldname, ierr)
    if (m > 0) then
       call infld(fieldname, fh_ini, dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
          tptr3d, found, grid_map='phys')
       if(found) then
          do n = 1, pbuf_times
             call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
          end do
       else
          call cnst_get_ind('CLDICE', ixcldice)
          call infld('CLDICE',fh_ini,dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
             tptr3d, found, grid_map='PHYS')
          if(found) then
             do n = 1, pbuf_times
                call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
             end do
          else
             call pbuf_set_field(pbuf2d, m, 0._r8)
          end if
          if (masterproc) then
             if (found) then
                write(iulog,*) trim(fieldname), ' initialized with CLDICE'
             else
                write(iulog,*) trim(fieldname), ' initialized to 0.0'
             end if
          end if
       end if
    end if

    fieldname = 'LCWAT'
    m = pbuf_get_index(fieldname,ierr)
    if (m > 0) then
       call infld(fieldname, fh_ini, dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
            tptr3d, found, grid_map='phys')
       if(found) then
          do n = 1, pbuf_times
             call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
          end do
       else
          allocate(tptr3d_2(pcols,pver,begchunk:endchunk))     
          call cnst_get_ind('CLDICE', ixcldice)
          call cnst_get_ind('CLDLIQ', ixcldliq)
          call infld('CLDICE',fh_ini,dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
               tptr3d, found, grid_map='PHYS')
          call infld('CLDLIQ',fh_ini,dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
               tptr3d_2, found2, grid_map='PHYS')
          if(found .and. found2) then
             tptr3d(:,:,:)=tptr3d(:,:,:)+tptr3d_2(:,:,:)
             if (masterproc) write(iulog,*) trim(fieldname), ' initialized with CLDICE + CLDLIQ'
          else if (found) then ! Data already loaded in tptr3d
             if (masterproc) write(iulog,*) trim(fieldname), ' initialized with CLDICE only'
          else if (found2) then
             tptr3d(:,:,:)=tptr3d_2(:,:,:)
             if (masterproc) write(iulog,*) trim(fieldname), ' initialized with CLDLIQ only'
          end if

          if (found .or. found2) then
             do n = 1, pbuf_times
                call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
             end do
             if(dycore_is('LR')) call polar_average(pver, tptr3d) 	
          else
             call pbuf_set_field(pbuf2d, m, 0._r8)
             if (masterproc)  write(iulog,*) trim(fieldname), ' initialized to 0.0'
          end if
          deallocate(tptr3d_2)
       end if
    end if

    deallocate(tptr3d)
    allocate(tptr3d(pcols,pver,begchunk:endchunk))

    fieldname = 'TCWAT'
    m = pbuf_get_index(fieldname,ierr)
    if (m > 0) then
       call infld(fieldname, fh_ini, dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
            tptr3d, found, grid_map='phys')
       if(.not.found) then
          call infld('T', fh_ini, dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
               tptr3d, found, grid_map='phys')
          if(dycore_is('LR')) call polar_average(pver, tptr3d) 	
          if (masterproc) write(iulog,*) trim(fieldname), ' initialized with T'
       end if
       do n = 1, pbuf_times
          call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
       end do
    end if

    deallocate(tptr3d)
    allocate(tptr3d(pcols,pverp,begchunk:endchunk))

    fieldname = 'TKE'
    m = pbuf_get_index( 'tke')
    call infld(fieldname, fh_ini, dim1name, 'ilev', 'lat', 1, pcols, 1, pverp, begchunk, endchunk, &
         tptr3d, found, grid_map='phys')
    if (found) then
       call pbuf_set_field(pbuf2d, m, tptr3d)
    else
       call pbuf_set_field(pbuf2d, m, 0.01_r8)
       if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.01'
    end if


    fieldname = 'KVM'
    m = pbuf_get_index('kvm')
    call infld(fieldname, fh_ini, dim1name, 'ilev', 'lat', 1, pcols, 1, pverp, begchunk, endchunk, &
         tptr3d, found, grid_map='phys')
    if (found) then
       call pbuf_set_field(pbuf2d, m, tptr3d)
    else
       call pbuf_set_field(pbuf2d, m, 0._r8)
       if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.'
    end if


    fieldname = 'KVH'
    m = pbuf_get_index('kvh')
    call infld(fieldname, fh_ini, dim1name, 'ilev', 'lat', 1, pcols, 1, pverp, begchunk, endchunk, &
         tptr3d, found, grid_map='phys')
    if (found) then
       call pbuf_set_field(pbuf2d, m, tptr3d)
    else
       call pbuf_set_field(pbuf2d, m, 0._r8)
       if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.'
    end if

    deallocate(tptr3d)
    allocate(tptr3d(pcols,pver,begchunk:endchunk))

    fieldname = 'CONCLD'
    m = pbuf_get_index('CONCLD')
    call infld(fieldname, fh_ini, dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
         tptr3d, found, grid_map='phys')
    if(found) then
       do n = 1, pbuf_times
          call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
       end do
    else
       call pbuf_set_field(pbuf2d, m, 0._r8)
       if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.'
    end if

    deallocate (tptr3d)

    call initialize_short_lived_species(fh_ini, pbuf2d)
end subroutine phys_inidat


subroutine phys_init( phys_state, phys_tend, pbuf2d, cam_out )

    !----------------------------------------------------------------------- 
    ! 
    ! Initialization of physics package.
    ! 
    !-----------------------------------------------------------------------

    use physics_buffer,     only: physics_buffer_desc, pbuf_initialize, pbuf_get_index
    use physconst,          only: rair, cpair, cpwv, gravit, stebol, tmelt, &
                                  latvap, latice, rh2o, rhoh2o, pstd, zvir,         &
                                  karman, rhodair, physconst_init 
    use ref_pres,           only: pref_edge, pref_mid

    use aerosol_intr,       only: aerosol_init
    use carma_intr,         only: carma_init
    use cloud_rad_props,    only: cloud_rad_props_init
    use cam_control_mod,    only: nsrest  ! restart flag
    use check_energy,       only: check_energy_init
    use chemistry,          only: chem_init
    use prescribed_ozone,   only: prescribed_ozone_init
    use prescribed_ghg,     only: prescribed_ghg_init
    use prescribed_aero,    only: prescribed_aero_init
    use aerodep_flx,        only: aerodep_flx_init
    use aircraft_emit,      only: aircraft_emit_init
    use prescribed_volcaero,only: prescribed_volcaero_init
    use cloud_fraction,     only: cldfrc_init
    use co2_cycle,          only: co2_init, co2_transport
    use convect_deep,       only: convect_deep_init
    use convect_shallow,    only: convect_shallow_init
    use cam_diagnostics,    only: diag_init
    use gw_drag,            only: gw_inti
    use cam3_aero_data,     only: cam3_aero_data_on, cam3_aero_data_init
    use cam3_ozone_data,    only: cam3_ozone_data_on, cam3_ozone_data_init
    use radheat,            only: radheat_init
    use radiation,          only: radiation_init
    use cloud_diagnostics,  only: cloud_diagnostics_init
    use stratiform,         only: stratiform_init
    use phys_control,       only: phys_getopts, waccmx_is
    use wv_saturation,      only: wv_sat_init
    use microp_driver,      only: microp_driver_init
    use microp_aero,        only: microp_aero_init
    use macrop_driver,      only: macrop_driver_init
    use conv_water,         only: conv_water_init
    use tracers,            only: tracers_init
    use aoa_tracers,        only: aoa_tracers_init
    use rayleigh_friction,  only: rayleigh_friction_init
    use pbl_utils,          only: pbl_utils_init
    use vertical_diffusion, only: vertical_diffusion_init
    use dycore,             only: dycore_is
    use phys_debug_util,    only: phys_debug_init
    use rad_constituents,   only: rad_cnst_init
    use aer_rad_props,      only: aer_rad_props_init
#if ( defined WACCM_PHYS )
    use qbo,                only: qbo_init
    use iondrag,            only: iondrag_init
#endif
#if ( defined OFFLINE_DYN )
    use metdata,            only: metdata_phys_init
#endif
    use ionosphere,	   only: ionos_init  ! Initialization of ionosphere module (WACCM-X)
    use majorsp_diffusion,  only: mspd_init   ! Initialization of major species diffusion module (WACCM-X)
    use clubb_intr,         only: clubb_ini_cam
    use sslt_rebin,         only: sslt_rebin_init
    use tropopause,         only: tropopause_init
    use solar_data,         only: solar_data_init
    use rad_solar_var,      only: rad_solar_var_init
    ! +++ Nudging coe +++ !
    use nudging,            only: Nudge_Model,nudging_init
    ! ------------------- !

    ! Input/output arguments
    type(physics_state), pointer       :: phys_state(:)
    type(physics_tend ), pointer       :: phys_tend(:)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    type(cam_out_t),intent(inout)      :: cam_out(begchunk:endchunk)

    ! local variables
    integer :: lchnk

    character(len=16) :: microp_scheme 
    logical           :: do_clubb_sgs

    !-----------------------------------------------------------------------

    ! Get microphysics option
    call phys_getopts(microp_scheme_out = microp_scheme)
    call phys_getopts(do_clubb_sgs_out  = do_clubb_sgs )

    call physics_type_alloc(phys_state, phys_tend, begchunk, endchunk, pcols)

    do lchnk = begchunk, endchunk
       call physics_state_set_grid(lchnk, phys_state(lchnk))
    end do

    !-------------------------------------------------------------------------------------------
    ! Initialize any variables in physconst which are not temporally and/or spatially constant
    !------------------------------------------------------------------------------------------- 
    call physconst_init()

    ! Initialize debugging a physics column
    call phys_debug_init()

    call pbuf_initialize(pbuf2d)

    ! diag_init makes addfld calls for dynamics fields that are output from
    ! the physics decomposition
    call diag_init()

    call check_energy_init()

    call tracers_init()

    ! age of air tracers
    call aoa_tracers_init()

    teout_idx = pbuf_get_index( 'TEOUT')

    ! For adiabatic or ideal physics don't need to initialize any of the
    ! parameterizations below:
    if (adiabatic .or. ideal_phys) return

    if (nsrest .eq. 0) then
       call phys_inidat(cam_out, pbuf2d) 
    end if
    
    ! wv_saturation is relatively independent of everything else and
    ! low level, so init it early. Must at least do this before radiation.
    call wv_sat_init

    ! CAM3 prescribed aerosols
    if (cam3_aero_data_on) call cam3_aero_data_init(phys_state)

    ! Initialize rad constituents and their properties
    call rad_cnst_init()
    call aer_rad_props_init()
    call cloud_rad_props_init()

    ! Initialize some aerosol code
    call aerosol_init(pbuf2d)

    ! initialize carma
    call carma_init()

    ! solar irradiance data modules
    call solar_data_init()


    ! Prognostic chemistry.
    call chem_init(phys_state,pbuf2d)

    ! Prescribed tracers
    call prescribed_ozone_init()
    call prescribed_ghg_init()
    call prescribed_aero_init()
    call aerodep_flx_init()
    call aircraft_emit_init()
    call prescribed_volcaero_init()

    ! co2 cycle            
    if (co2_transport()) then
       call co2_init()
    end if

    ! CAM3 prescribed ozone
    if (cam3_ozone_data_on) call cam3_ozone_data_init(phys_state)

    call gw_inti(cpair, cpwv, gravit, rair, pref_edge)

    call rayleigh_friction_init()

    call pbl_utils_init(gravit, karman, cpair, rair, zvir)
    if (.not. do_clubb_sgs) call vertical_diffusion_init(pbuf2d)

    if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then
       call mspd_init ()
       ! Initialization of ionosphere module if mode set to ionosphere
       if( waccmx_is('ionosphere') ) then
          call ionos_init()
       endif
    endif

    call tsinti(tmelt, latvap, rair, stebol, latice)

    call radiation_init

    call rad_solar_var_init()

    call cloud_diagnostics_init

    call radheat_init(pref_mid)

    call convect_shallow_init(pref_edge)

    call cldfrc_init

    call convect_deep_init(pref_edge)

    if( microp_scheme == 'RK' ) then
       call stratiform_init()
    elseif( microp_scheme == 'MG' ) then 
       if (.not. do_clubb_sgs) call macrop_driver_init()
       call microp_aero_init()
       call microp_driver_init(pbuf2d)
       call conv_water_init
    end if

    ! initiate CLUBB within CAM
    if (do_clubb_sgs) call clubb_ini_cam(pbuf2d)

#if ( defined WACCM_PHYS )
    call iondrag_init(pref_mid)
    call qbo_init
#endif

#if ( defined OFFLINE_DYN )
    call metdata_phys_init()
#endif
    call sslt_rebin_init()
    call tropopause_init()

    prec_dp_idx  = pbuf_get_index('PREC_DP')
    snow_dp_idx  = pbuf_get_index('SNOW_DP')
    prec_sh_idx  = pbuf_get_index('PREC_SH')
    snow_sh_idx  = pbuf_get_index('SNOW_SH')

    ! +++ Nudging code +++ !
    ! Initialize Nudging Parameters
    !--------------------------------
    if(Nudge_Model) call nudging_init
    ! -------------------- !

end subroutine phys_init

  !
  !-----------------------------------------------------------------------
  !

subroutine phys_run1(phys_state, ztodt, phys_tend, pbuf2d,  cam_in, cam_out)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! First part of atmospheric physics package before updating of surface models
    ! 
    !-----------------------------------------------------------------------
    use time_manager,   only: get_nstep
    use cam_diagnostics,only: diag_allocate, diag_physvar_ic
    use check_energy,   only: check_energy_gmean

    use physics_buffer,         only: physics_buffer_desc, pbuf_get_chunk, pbuf_allocate
#if (defined BFB_CAM_SCAM_IOP )
    use cam_history,    only: outfld
#endif
    use comsrf,         only: fsns, fsnt, flns, sgh30, flnt, landm, fsds
    use abortutils,     only: endrun
#if ( defined OFFLINE_DYN )
     use metdata,       only: get_met_srf1
#endif
    !
    ! Input arguments
    !
    real(r8), intent(in) :: ztodt            ! physics time step unless nstep=0
    !
    ! Input/Output arguments
    !
    type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
    type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend

    type(physics_buffer_desc), pointer, dimension(:,:) :: pbuf2d
    type(cam_in_t),                     dimension(begchunk:endchunk) :: cam_in
    type(cam_out_t),                    dimension(begchunk:endchunk) :: cam_out
    !-----------------------------------------------------------------------
    !
    !---------------------------Local workspace-----------------------------
    !
    integer :: c                                 ! indices
    integer :: ncol                              ! number of columns
    integer :: nstep                             ! current timestep number
#if (! defined SPMD)
    integer  :: mpicom = 0
#endif
    type(physics_buffer_desc), pointer :: phys_buffer_chunk(:)

    call t_startf ('physpkg_st1')
    nstep = get_nstep()

#if ( defined OFFLINE_DYN )
    !
    ! if offline mode set SNOWH and TS for micro-phys
    !
    call get_met_srf1( cam_in )
#endif

    ! The following initialization depends on the import state (cam_in)
    ! being initialized.  This isn't true when cam_init is called, so need
    ! to postpone this initialization to here.
    if (nstep == 0 .and. phys_do_flux_avg()) call flux_avg_init(cam_in,  pbuf2d)

    ! Compute total energy of input state and previous output state
    call t_startf ('chk_en_gmean')
    call check_energy_gmean(phys_state, pbuf2d, ztodt, nstep)
    call t_stopf ('chk_en_gmean')

    call t_stopf ('physpkg_st1')

    if ( adiabatic .or. ideal_phys )then
       call t_startf ('bc_physics')
       call phys_run1_adiabatic_or_ideal(ztodt, phys_state, phys_tend,  pbuf2d)
       call t_stopf ('bc_physics')
    else
       call t_startf ('physpkg_st1')

       call pbuf_allocate(pbuf2d, 'physpkg')
       call diag_allocate()

       !-----------------------------------------------------------------------
       ! Advance time information
       !-----------------------------------------------------------------------

       call phys_timestep_init( phys_state, cam_out, pbuf2d)

       call t_stopf ('physpkg_st1')

#ifdef TRACER_CHECK
       call gmean_mass ('before tphysbc DRY', phys_state)
#endif


       !-----------------------------------------------------------------------
       ! Tendency physics before flux coupler invocation
       !-----------------------------------------------------------------------
       !

#if (defined BFB_CAM_SCAM_IOP )
       do c=begchunk, endchunk
          call outfld('Tg',cam_in(c)%ts,pcols   ,c     )
       end do
#endif

       call t_barrierf('sync_bc_physics', mpicom)
       call t_startf ('bc_physics')
       call t_adj_detailf(+1)

!$OMP PARALLEL DO PRIVATE (C, phys_buffer_chunk)
       do c=begchunk, endchunk
          !
          ! Output physics terms to IC file
          !
          phys_buffer_chunk => pbuf_get_chunk(pbuf2d, c)

          call t_startf ('diag_physvar_ic')
          call diag_physvar_ic ( c,  phys_buffer_chunk, cam_out(c), cam_in(c) )
          call t_stopf ('diag_physvar_ic')

          call tphysbc (ztodt, fsns(1,c), fsnt(1,c), flns(1,c), flnt(1,c), phys_state(c),        &
                       phys_tend(c), phys_buffer_chunk,  fsds(1,c), landm(1,c),          &
                       sgh30(1,c), cam_out(c), cam_in(c) )

       end do

       call t_adj_detailf(-1)
       call t_stopf ('bc_physics')

       ! Don't call the rest in CRM mode
       if(single_column.and.scm_crm_mode) return

#ifdef TRACER_CHECK
       call gmean_mass ('between DRY', phys_state)
#endif
    end if

end subroutine phys_run1

  !
  !-----------------------------------------------------------------------
  !

subroutine phys_run1_adiabatic_or_ideal(ztodt, phys_state, phys_tend,  pbuf2d)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Physics for adiabatic or idealized physics case.
    ! 
    !-----------------------------------------------------------------------
    use physics_buffer, only : physics_buffer_desc, pbuf_set_field, pbuf_get_chunk, pbuf_old_tim_idx
    use time_manager,     only: get_nstep
    use cam_diagnostics,  only: diag_phys_writeout
    use check_energy,     only: check_energy_fix, check_energy_chng
    use dycore,           only: dycore_is

    !
    ! Input arguments
    !
    real(r8), intent(in) :: ztodt            ! physics time step unless nstep=0
    !
    ! Input/Output arguments
    !
    type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
    type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend

    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    !-----------------------------------------------------------------------
    !---------------------------Local workspace-----------------------------
    !
    integer             :: c               ! indices
    integer             :: nstep           ! current timestep number
    type(physics_ptend) :: ptend(begchunk:endchunk) ! indivdual parameterization tendencies
    real(r8)            :: flx_heat(pcols) ! effective sensible heat flux
    real(r8)            :: zero(pcols)     ! array of zeros

    ! physics buffer field for total energy
    integer itim
    real(r8), pointer, dimension(:) :: teout
    logical, SAVE :: first_exec_of_phys_run1_adiabatic_or_ideal  = .TRUE.
    !-----------------------------------------------------------------------

    nstep = get_nstep()
    zero  = 0._r8

    ! Associate pointers with physics buffer fields
    itim = pbuf_old_tim_idx()
    if (first_exec_of_phys_run1_adiabatic_or_ideal) then
       first_exec_of_phys_run1_adiabatic_or_ideal  = .FALSE.
    endif

!$OMP PARALLEL DO PRIVATE (C, FLX_HEAT)
    do c=begchunk, endchunk

       ! Initialize the physics tendencies to zero.
       call physics_tend_init(phys_tend(c))

       ! Dump dynamics variables to history buffers
       call diag_phys_writeout(phys_state(c))

       if (dycore_is('LR')) then
          call check_energy_fix(phys_state(c), ptend(c), nstep, flx_heat)
          call physics_update(phys_state(c), ptend(c), ztodt, phys_tend(c))
          call check_energy_chng(phys_state(c), phys_tend(c), "chkengyfix", nstep, ztodt, &
               zero, zero, zero, flx_heat)
          call physics_ptend_dealloc(ptend(c))
       end if

       if ( ideal_phys )then
          call t_startf('tphysidl')
          call tphysidl(ztodt, phys_state(c), phys_tend(c))
          call t_stopf('tphysidl')
       end if

       ! Save total enery after physics for energy conservation checks
       call pbuf_set_field(pbuf_get_chunk(pbuf2d, c), teout_idx, phys_state(c)%te_cur)


    end do

end subroutine phys_run1_adiabatic_or_ideal

  !
  !-----------------------------------------------------------------------
  !

subroutine phys_run2(phys_state, ztodt, phys_tend, pbuf2d,  cam_out, &
       cam_in )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Second part of atmospheric physics package after updating of surface models
    ! 
    !-----------------------------------------------------------------------
    use physics_buffer,         only: physics_buffer_desc, pbuf_get_chunk, pbuf_deallocate, pbuf_update_tim_idx
    use mo_lightning,   only: lightning_no_prod


    use cam_diagnostics,only: diag_deallocate, diag_surf
    use comsrf,         only: trefmxav, trefmnav, sgh, sgh30, fsds 
    use physconst,      only: stebol, latvap
    use carma_intr,     only: carma_accumulate_stats
#if ( defined OFFLINE_DYN )
    use metdata,        only: get_met_srf2
#endif
    !
    ! Input arguments
    !
    real(r8), intent(in) :: ztodt                       ! physics time step unless nstep=0
    !
    ! Input/Output arguments
    !
    type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
    type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend
    type(physics_buffer_desc),pointer, dimension(:,:)     :: pbuf2d

    type(cam_out_t),     intent(inout), dimension(begchunk:endchunk) :: cam_out
    type(cam_in_t),      intent(inout), dimension(begchunk:endchunk) :: cam_in
    !
    !-----------------------------------------------------------------------
    !---------------------------Local workspace-----------------------------
    !
    integer :: c                                 ! chunk index
    integer :: ncol                              ! number of columns
#if (! defined SPMD)
    integer  :: mpicom = 0
#endif
    type(physics_buffer_desc),pointer, dimension(:)     :: phys_buffer_chunk
    !
    ! If exit condition just return
    !

    if(single_column.and.scm_crm_mode) return

    if ( adiabatic .or. ideal_phys ) return
    !-----------------------------------------------------------------------
    ! Tendency physics after coupler 
    ! Not necessary at terminal timestep.
    !-----------------------------------------------------------------------
    !
#if ( defined OFFLINE_DYN )
    !
    ! if offline mode set SHFLX QFLX TAUX TAUY for vert diffusion
    !
    call get_met_srf2( cam_in )
#endif
    ! Set lightning production of NO
    call t_startf ('lightning_no_prod')
    call lightning_no_prod( phys_state, pbuf2d,  cam_in )
    call t_stopf ('lightning_no_prod')

    call t_barrierf('sync_ac_physics', mpicom)
    call t_startf ('ac_physics')
    call t_adj_detailf(+1)

!$OMP PARALLEL DO PRIVATE (C, NCOL, phys_buffer_chunk)

    do c=begchunk,endchunk
       ncol = get_ncols_p(c)
       phys_buffer_chunk => pbuf_get_chunk(pbuf2d, c)
       !
       ! surface diagnostics for history files
       !
       call t_startf('diag_surf')
       call diag_surf(cam_in(c), cam_out(c), phys_state(c)%ps,trefmxav(1,c), trefmnav(1,c))
       call t_stopf('diag_surf')

       call tphysac(ztodt, cam_in(c),  &
            sgh(1,c), sgh30(1,c), cam_out(c),                              &
            phys_state(c), phys_tend(c), phys_buffer_chunk,&
            fsds(1,c))
    end do                    ! Chunk loop

    call t_adj_detailf(-1)
    call t_stopf('ac_physics')

#ifdef TRACER_CHECK
    call gmean_mass ('after tphysac FV:WET)', phys_state)
#endif

    call t_startf ('carma_accumulate_stats')
    call carma_accumulate_stats()
    call t_stopf ('carma_accumulate_stats')

    call t_startf ('physpkg_st2')
    call pbuf_deallocate(pbuf2d, 'physpkg')

    call pbuf_update_tim_idx()
    call diag_deallocate()
    call t_stopf ('physpkg_st2')

end subroutine phys_run2

  !
  !----------------------------------------------------------------------- 
  !

subroutine phys_final( phys_state, phys_tend, pbuf2d )
    use physics_buffer, only : physics_buffer_desc, pbuf_deallocate
    use chemistry, only : chem_final
    use carma_intr, only : carma_final
    use wv_saturation, only : wv_sat_final
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Finalization of physics package
    ! 
    !-----------------------------------------------------------------------
    ! Input/output arguments
    type(physics_state), pointer :: phys_state(:)
    type(physics_tend ), pointer :: phys_tend(:)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    if(associated(pbuf2d)) then
       call pbuf_deallocate(pbuf2d,'global')
       deallocate(pbuf2d)
    end if
    deallocate(phys_state)
    deallocate(phys_tend)
    call chem_final
    call carma_final
    call wv_sat_final

end subroutine phys_final


subroutine tphysac (ztodt,   cam_in,  &
       sgh,     sgh30,                                     &
       cam_out,  state,   tend,    pbuf,            &
       fsds    )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Tendency physics after coupling to land, sea, and ice models.
    ! Computes the following:
    !   o Radon surface flux and decay (optional)
    !   o Vertical diffusion and planetary boundary layer
    !   o Multiple gravity wave drag
    ! 
    ! Method: 
    ! <Describe the algorithm(s) used in the routine.> 
    ! <Also include any applicable external references.> 
    ! 
    ! Author: CCM1, CMS Contact: J. Truesdale
    ! 
    !-----------------------------------------------------------------------
    use physics_buffer, only: physics_buffer_desc, pbuf_set_field, pbuf_get_index, pbuf_get_field, pbuf_old_tim_idx
    use shr_kind_mod,       only: r8 => shr_kind_r8
    use chemistry,          only: chem_is_active, chem_timestep_tend
    use cam_diagnostics,    only: diag_phys_tend_writeout
    use gw_drag,            only: gw_intr
    use vertical_diffusion, only: vertical_diffusion_tend
    use rayleigh_friction,  only: rayleigh_friction_tend
    use constituents,       only: cnst_get_ind
    use physics_types,      only: physics_state, physics_tend, physics_ptend, physics_update,    &
         physics_dme_adjust, set_dry_to_wet, physics_state_check
    use majorsp_diffusion,  only: mspd_intr  ! WACCM-X major diffusion
    use ionosphere,         only: ionos_intr ! WACCM-X ionosphere
    use phys_control,       only: phys_getopts
    use tracers,            only: tracers_timestep_tend
    use aoa_tracers,        only: aoa_tracers_timestep_tend
    use physconst,          only: rhoh2o, latvap,latice
    use aerosol_intr,       only: aerosol_emis_intr, aerosol_drydep_intr
    use carma_intr,         only: carma_emission_tend, carma_timestep_tend
    use carma_flags_mod,    only: carma_do_aerosol, carma_do_emission
    use check_energy,       only: check_energy_chng
    use check_energy,       only: check_tracers_data, check_tracers_init, check_tracers_chng
    use time_manager,       only: get_nstep
    use abortutils,         only: endrun
    use dycore,             only: dycore_is
    use cam_control_mod,    only: aqua_planet 
    use mo_gas_phase_chemdr,only: map2chm
    use clybry_fam,         only: clybry_fam_set
    use charge_neutrality,  only: charge_fix
#if ( defined WACCM_PHYS )
    use iondrag,            only: iondrag_calc, do_waccm_ions
    use qbo,                only: qbo_relax
#endif
    use clubb_intr,         only: clubb_surface
    use perf_mod
    use phys_control,       only: phys_do_flux_avg, waccmx_is
    use flux_avg,           only: flux_avg_run
    ! +++ Nudging code +++ !
    use nudging,            only: Nudge_Model,Nudge_ON,nudging_timestep_tend,nudging_diag,Nudge_Diag_Opt
    ! -------------------- !

    implicit none

    !
    ! Arguments
    !
    real(r8), intent(in) :: ztodt                  ! Two times model timestep (2 delta-t)
    real(r8), intent(in) :: fsds(pcols)            ! down solar flux
    real(r8), intent(in) :: sgh(pcols)             ! Std. deviation of orography for gwd
    real(r8), intent(in) :: sgh30(pcols)           ! Std. deviation of 30s orography for tms

    type(cam_in_t),      intent(inout) :: cam_in
    type(cam_out_t),     intent(inout) :: cam_out
    type(physics_state), intent(inout) :: state
    type(physics_tend ), intent(inout) :: tend
    type(physics_buffer_desc), pointer :: pbuf(:)


    type(check_tracers_data):: tracerint             ! tracer mass integrals and cummulative boundary fluxes

    !
    !---------------------------Local workspace-----------------------------
    !
    type(physics_ptend)     :: ptend               ! indivdual parameterization tendencies

    integer  :: nstep                              ! current timestep number
    real(r8) :: zero(pcols)                        ! array of zeros

    integer :: lchnk                                ! chunk identifier
    integer :: ncol                                 ! number of atmospheric columns
    integer i,k,m                 ! Longitude, level indices
    integer :: yr, mon, day, tod       ! components of a date
    integer :: ixcldice, ixcldliq      ! constituent indices for cloud liquid and ice water.

    logical :: labort                            ! abort flag

    real(r8) tvm(pcols,pver)           ! virtual temperature
    real(r8) prect(pcols)              ! total precipitation
    real(r8) surfric(pcols)            ! surface friction velocity
    real(r8) obklen(pcols)             ! Obukhov length
    real(r8) :: fh2o(pcols)            ! h2o flux to balance source from methane chemistry
    real(r8) :: tmp_q     (pcols,pver) ! tmp space
    real(r8) :: tmp_cldliq(pcols,pver) ! tmp space
    real(r8) :: tmp_cldice(pcols,pver) ! tmp space
    real(r8) :: tmp_t     (pcols,pver) ! tmp space

    ! physics buffer fields for total energy and mass adjustment
    integer itim, ifld

    real(r8), pointer, dimension(:,:) :: tini
    real(r8), pointer, dimension(:,:) :: cld
    real(r8), pointer, dimension(:,:) :: qini
    real(r8), pointer, dimension(:,:) :: cldliqini
    real(r8), pointer, dimension(:,:) :: cldiceini
    real(r8), pointer, dimension(:,:) :: dtcore
    real(r8), pointer, dimension(:,:) :: ast     ! relative humidity cloud fraction 

    logical :: do_clubb_sgs 

    ! Debug physics_state.
    logical :: state_debug_checks
    !
    !-----------------------------------------------------------------------
    !
    lchnk = state%lchnk
    ncol  = state%ncol

    nstep = get_nstep()
    
    call phys_getopts( do_clubb_sgs_out       = do_clubb_sgs, &
                       state_debug_checks_out = state_debug_checks)

    ! Adjust the surface fluxes to reduce instabilities in near sfc layer
    if (phys_do_flux_avg()) then 
       call flux_avg_run(state, cam_in,  pbuf, nstep, ztodt)
    endif

    ! Validate the physics state.
    if (state_debug_checks) &
         call physics_state_check(state, name="before tphysac")

    call t_startf('tphysac_init')
    ! Associate pointers with physics buffer fields
    itim = pbuf_old_tim_idx()


    ifld = pbuf_get_index('DTCORE')
    call pbuf_get_field(pbuf, ifld, dtcore, start=(/1,1,itim/), kount=(/pcols,pver,1/) )

    call pbuf_get_field(pbuf, tini_idx, tini)
    call pbuf_get_field(pbuf, qini_idx, qini)
    call pbuf_get_field(pbuf, cldliqini_idx, cldliqini)
    call pbuf_get_field(pbuf, cldiceini_idx, cldiceini)

    ifld = pbuf_get_index('CLD')
    call pbuf_get_field(pbuf, ifld, cld, start=(/1,1,itim/),kount=(/pcols,pver,1/))

    ifld = pbuf_get_index('AST')
    call pbuf_get_field(pbuf, ifld, ast, start=(/1,1,itim/), kount=(/pcols,pver,1/) )

    !
    ! accumulate fluxes into net flux array for spectral dycores
    ! jrm Include latent heat of fusion for snow
    !
    do i=1,ncol
       tend%flx_net(i) = tend%flx_net(i) + cam_in%shf(i) + (cam_out%precc(i) &
            + cam_out%precl(i))*latvap*rhoh2o &
            + (cam_out%precsc(i) + cam_out%precsl(i))*latice*rhoh2o
    end do

    ! emission of aerosols at surface
    call aerosol_emis_intr (state, cam_in)

    if (carma_do_emission) then
       ! carma emissions
       call carma_emission_tend (state, ptend, cam_in, ztodt)
       call physics_update(state, ptend, ztodt, tend)
    end if

    ! get nstep and zero array for energy checker
    zero = 0._r8
    nstep = get_nstep()
    call check_tracers_init(state, tracerint)

    ! Check if latent heat flux exceeds the total moisture content of the
    ! lowest model layer, thereby creating negative moisture.

    call qneg4('TPHYSAC '       ,lchnk               ,ncol  ,ztodt ,               &
         state%q(1,pver,1),state%rpdel(1,pver) ,cam_in%shf ,         &
         cam_in%lhf , cam_in%cflx )

    call t_stopf('tphysac_init')
    !===================================================
    ! Source/sink terms for advected tracers.
    !===================================================
    call t_startf('adv_tracer_src_snk')
    ! Test tracers

    call tracers_timestep_tend(state, ptend, cam_in%cflx, cam_in%landfrac, ztodt)      
    call physics_update(state, ptend, ztodt, tend)
    call check_tracers_chng(state, tracerint, "tracers_timestep_tend", nstep, ztodt,   &
         cam_in%cflx)

    call aoa_tracers_timestep_tend(state, ptend, cam_in%cflx, cam_in%landfrac, ztodt)      
    call physics_update(state, ptend, ztodt, tend)
    call check_tracers_chng(state, tracerint, "aoa_tracers_timestep_tend", nstep, ztodt,   &
         cam_in%cflx)

    ! Chemistry calculation
    if (chem_is_active()) then
       call chem_timestep_tend(state, ptend, cam_in, cam_out, ztodt, &
            pbuf,  fh2o, fsds)

       call physics_update(state, ptend, ztodt, tend)
       call check_energy_chng(state, tend, "chem", nstep, ztodt, fh2o, zero, zero, zero)
       call check_tracers_chng(state, tracerint, "chem_timestep_tend", nstep, ztodt, &
            cam_in%cflx)
    end if
    call t_stopf('adv_tracer_src_snk')

    !===================================================
    ! Vertical diffusion/pbl calculation
    ! Call vertical diffusion code (pbl, free atmosphere and molecular)
    !===================================================

    ! If CLUBB is called, do not call vertical diffusion, but obukov length and
    !   surface friction velocity still need to be computed.  In addition, 
    !   surface fluxes need to be updated here for constituents 
    if (do_clubb_sgs) then

       call clubb_surface ( state, ptend, ztodt, cam_in, surfric, obklen)
       
       ! Update surface flux constituents 
       call physics_update(state, ptend, ztodt, tend)

    else

       call t_startf('vertical_diffusion_tend')
       call vertical_diffusion_tend (ztodt ,state ,cam_in%wsx, cam_in%wsy,   &
            cam_in%shf     ,cam_in%cflx     ,surfric  ,obklen   ,ptend    ,ast    ,&
            cam_in%ocnfrac  , cam_in%landfrac ,        &
            sgh30    ,pbuf )

    !------------------------------------------
    ! Call major diffusion for extended model
    !------------------------------------------
    if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then
       call mspd_intr (ztodt    ,state    ,ptend)
    endif

       call physics_update(state, ptend, ztodt, tend)
       call t_stopf ('vertical_diffusion_tend')
    
    endif


    !===================================================
    ! Rayleigh friction calculation
    !===================================================
    call t_startf('rayleigh_friction')
    call rayleigh_friction_tend( ztodt, state, ptend)
    call physics_update(state, ptend, ztodt, tend)
    call t_stopf('rayleigh_friction')

    if (do_clubb_sgs) then
      call check_energy_chng(state, tend, "vdiff", nstep, ztodt, zero, zero, zero, zero)
    else
      call check_energy_chng(state, tend, "vdiff", nstep, ztodt, cam_in%cflx(:,1), zero, &
           zero, cam_in%shf)
    endif
    
    call check_tracers_chng(state, tracerint, "vdiff", nstep, ztodt, cam_in%cflx)

    !  aerosol dry deposition processes
    call t_startf('aero_drydep')
    call aerosol_drydep_intr (state, ptend, cam_in, cam_out, ztodt,  &
         fsds, obklen, surfric, prect, pbuf)
    call physics_update(state, ptend, ztodt, tend)
    call t_stopf('aero_drydep')

   ! CARMA microphysics
   !
   ! NOTE: This does both the timestep_tend for CARMA aerosols as well as doing the dry
   ! deposition for CARMA aerosols. It needs to follow vertical_diffusion_tend, so that
   ! obklen and surfric have been calculated. It needs to follow aerosol_drydep_intr, so
   ! that cam_out%xxxdryxxx fields have already been set for CAM aerosols and cam_out
   ! can be added to for CARMA aerosols.
   if (carma_do_aerosol) then
     call t_startf('carma_timestep_tend')
     call carma_timestep_tend(state, cam_in, cam_out, ptend, ztodt, pbuf, obklen=obklen, ustar=surfric)
     call physics_update(state, ptend, ztodt, tend)
   
     call check_energy_chng(state, tend, "carma_tend", nstep, ztodt, zero, zero, zero, zero)
     call t_stopf('carma_timestep_tend')
   end if


    !---------------------------------------------------------------------------------
    !	... enforce charge neutrality
    !---------------------------------------------------------------------------------
    if (do_waccm_phys()) call charge_fix( ncol, state%q(:,:,:) )

    !===================================================
    ! Gravity wave drag
    !===================================================
    call t_startf('gw_intr')

    call gw_intr(state, sgh, pbuf, ztodt, ptend, cam_in%landfrac)

    call physics_update(state, ptend, ztodt, tend)
    ! Check energy integrals
    call check_energy_chng(state, tend, "gwdrag", nstep, ztodt, zero, zero, zero, zero)
    call t_stopf('gw_intr')

#if ( defined WACCM_PHYS )

    ! QBO relaxation
    call qbo_relax(state, ptend)
    call physics_update(state, ptend, ztodt, tend)
    ! Check energy integrals
    call check_energy_chng(state, tend, "qborelax", nstep, ztodt, zero, zero, zero, zero)
    ! Ion drag calculation
    call t_startf ( 'iondrag' )

    if ( do_waccm_ions ) then
       call iondrag_calc( lchnk, ncol, state, ptend, pbuf,  ztodt )
    else
       call iondrag_calc( lchnk, ncol, state, ptend)
    endif
    !----------------------------------------------------------------------------
    ! Call ionosphere routines for extended model if mode is set to ionosphere
    !----------------------------------------------------------------------------
    if( waccmx_is('ionosphere') ) then
       call ionos_intr(state, ptend, pbuf, ztodt)
    endif

    call physics_update(state, ptend, ztodt, tend)
    ! Check energy integrals
    call check_energy_chng(state, tend, "iondrag", nstep, ztodt, zero, zero, zero, zero)
    call t_stopf  ( 'iondrag' )

#endif


    !-------------- Energy budget checks ----------------------------------- 

    call pbuf_set_field(pbuf, teout_idx, state%te_cur, (/1,itim/),(/pcols,1/))       

    !*** BAB's FV heating kludge *** apply the heating as temperature tendency.
    !*** BAB's FV heating kludge *** modify the temperature in the state structure
    tmp_t(:ncol,:pver) = state%t(:ncol,:pver)
    state%t(:ncol,:pver) = tini(:ncol,:pver) + ztodt*tend%dtdt(:ncol,:pver)

    ! store dse after tphysac in buffer
    do k = 1,pver
       dtcore(:ncol,k) = state%t(:ncol,k)
    end do


    !
    ! FV: convert dry-type mixing ratios to moist here because physics_dme_adjust
    !     assumes moist. This is done in p_d_coupling for other dynamics. Bundy, Feb 2004.


    if ( dycore_is('LR') .or. dycore_is('SE')) call set_dry_to_wet(state)    ! Physics had dry, dynamics wants moist


    ! Scale dry mass and energy (does nothing if dycore is EUL or SLD)
    call cnst_get_ind('CLDLIQ', ixcldliq)
    call cnst_get_ind('CLDICE', ixcldice)
    tmp_q     (:ncol,:pver) = state%q(:ncol,:pver,1)
    tmp_cldliq(:ncol,:pver) = state%q(:ncol,:pver,ixcldliq)
    tmp_cldice(:ncol,:pver) = state%q(:ncol,:pver,ixcldice)
    call physics_dme_adjust(state, tend, qini, ztodt)
!!!   REMOVE THIS CALL, SINCE ONLY Q IS BEING ADJUSTED. WON'T BALANCE ENERGY. TE IS SAVED BEFORE THIS
!!!   call check_energy_chng(state, tend, "drymass", nstep, ztodt, zero, zero, zero, zero)

    !-------------- Energy budget checks ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    if (aqua_planet) then
       labort = .false.
       do i=1,ncol
          if (cam_in%ocnfrac(i) /= 1._r8) labort = .true.
       end do
       if (labort) then
          call endrun ('TPHYSAC error:  grid contains non-ocean point')
       endif
    endif
    ! +++ Nudging code +++ !
    !===================================================
    ! Update Nudging values, if needed
    !===================================================
    if((Nudge_Model).and.(Nudge_ON)) then
      call nudging_timestep_tend(state,ptend)
      call physics_update(state,ptend,ztodt,tend)
    endif
    ! ---------------------!

    call diag_phys_tend_writeout (state, pbuf,  tend, ztodt, tmp_q, tmp_cldliq, tmp_cldice, &
         tmp_t, qini, cldliqini, cldiceini)

    call clybry_fam_set( ncol, lchnk, map2chm, state%q, pbuf )

end subroutine tphysac

subroutine tphysbc (ztodt,               &
       fsns,    fsnt,    flns,    flnt,    state,   &
       tend,    pbuf,     fsds,    landm,            &
       sgh30, cam_out, cam_in )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Evaluate and apply physical processes that are calculated BEFORE 
    ! coupling to land, sea, and ice models.  
    !
    ! Processes currently included are: 
    ! dry adjustment, moist convection, stratiform, wet deposition, radiation
    !
    ! Pass surface fields for separate surface flux calculations
    ! Dump appropriate fields to history file.
    ! 
    ! Method: 
    !
    ! Each parameterization should be implemented with this sequence of calls:
    !  1)  Call physics interface
    !  2)  Check energy
    !  3)  Call physics_update
    ! See Interface to Column Physics and Chemistry Packages 
    !   http://www.ccsm.ucar.edu/models/atm-cam/docs/phys-interface/index.html
    ! 
    ! Author: CCM1, CMS Contact: J. Truesdale
    !         modified by A. Gettelman and C. Craig Nov 2010 to separate micro/macro physics
    ! 
    !-----------------------------------------------------------------------

    use physics_buffer,          only : physics_buffer_desc, pbuf_get_field
    use physics_buffer,          only : pbuf_get_index, pbuf_old_tim_idx, pbuf_times
    use shr_kind_mod,    only: r8 => shr_kind_r8

    use stratiform,      only: stratiform_tend
    use phys_control,    only: phys_getopts
    use microp_driver,   only: microp_driver_tend
    use microp_aero,     only: microp_aero_run
    use macrop_driver,   only: macrop_driver_tend
    use physics_types,   only: physics_state, physics_tend, physics_ptend, physics_update, &
         physics_ptend_init, physics_ptend_sum, physics_state_check
    use cam_diagnostics, only: diag_conv_tend_ini, diag_phys_writeout, diag_conv, diag_export, diag_state_b4_phys_write
    use cam_history,     only: outfld
    use physconst,       only: cpair, latvap
    use constituents,    only: pcnst, qmin, cnst_get_ind
    use convect_deep,    only: convect_deep_tend, convect_deep_tend_2, deep_scheme_does_scav_trans
    use time_manager,    only: is_first_step, get_nstep
    use convect_shallow, only: convect_shallow_tend
    use check_energy,    only: check_energy_chng, check_energy_fix
    use check_energy,    only: check_tracers_data, check_tracers_init, check_tracers_chng
    use dycore,          only: dycore_is
    use aerosol_intr,    only: aerosol_wet_intr
    use carma_intr,      only: carma_wetdep_tend, carma_timestep_tend
    use carma_flags_mod, only: carma_do_detrain, carma_do_cldice, carma_do_cldliq,  carma_do_wetdep
    use radiation,       only: radiation_tend
    use cloud_diagnostics, only: cloud_diagnostics_calc
    use perf_mod
#ifdef MODAL_AERO
    use modal_aero_data, only: qneg3_worst_thresh_amode
#endif
    use mo_gas_phase_chemdr,only: map2chm
    use clybry_fam,         only: clybry_fam_adj
    use clubb_intr,      only: clubb_tend_cam
    use sslt_rebin,      only: sslt_rebin_adv
    use tropopause,      only: tropopause_output
    use abortutils,      only: endrun

    implicit none

    !
    ! Arguments
    !
    real(r8), intent(in) :: ztodt                          ! 2 delta t (model time increment)
    real(r8), intent(inout) :: fsns(pcols)                   ! Surface solar absorbed flux
    real(r8), intent(inout) :: fsnt(pcols)                   ! Net column abs solar flux at model top
    real(r8), intent(inout) :: flns(pcols)                   ! Srf longwave cooling (up-down) flux
    real(r8), intent(inout) :: flnt(pcols)                   ! Net outgoing lw flux at model top
    real(r8), intent(inout) :: fsds(pcols)                   ! Surface solar down flux
    real(r8), intent(in) :: landm(pcols)                   ! land fraction ramp
    real(r8), intent(in) :: sgh30(pcols)                   ! Std. deviation of 30 s orography for tms

    type(physics_state), intent(inout) :: state
    type(physics_tend ), intent(inout) :: tend
    type(physics_buffer_desc), pointer :: pbuf(:)

    type(cam_out_t),     intent(inout) :: cam_out
    type(cam_in_t),      intent(in)    :: cam_in


    !
    !---------------------------Local workspace-----------------------------
    !

    type(physics_ptend)   :: ptend            ! indivdual parameterization tendencies
    type(physics_state)   :: state_sc         ! state for sub-columns
    type(physics_ptend)   :: ptend_sc         ! ptend for sub-columns
    type(physics_ptend)   :: ptend_aero       ! ptend for microp_aero
    type(physics_tend)    :: tend_sc          ! tend for sub-columns

    integer :: nstep                          ! current timestep number

    real(r8) :: net_flx(pcols)

    real(r8) :: zdu(pcols,pver)               ! detraining mass flux from deep convection
    real(r8) :: cmfmc(pcols,pverp)            ! Convective mass flux--m sub c

    real(r8) cmfcme(pcols,pver)                ! cmf condensation - evaporation
    real(r8) cmfmc2(pcols,pverp)               ! Moist convection cloud mass flux
    real(r8) dlf(pcols,pver)                   ! Detraining cld H20 from shallow + deep convections
    real(r8) dlf2(pcols,pver)                  ! Detraining cld H20 from shallow convections
    real(r8) pflx(pcols,pverp)                 ! Conv rain flux thru out btm of lev
    real(r8) rtdt                              ! 1./ztodt

    integer lchnk                              ! chunk identifier
    integer ncol                               ! number of atmospheric columns
    integer ierr

    integer  i,k,m                             ! Longitude, level, constituent indices
    integer :: ixcldice, ixcldliq              ! constituent indices for cloud liquid and ice water.

    ! physics buffer fields to compute tendencies for stratiform package
    integer itim, ifld
    real(r8), pointer, dimension(:,:) :: cld        ! cloud fraction


    ! physics buffer fields for total energy and mass adjustment
    real(r8), pointer, dimension(:  ) :: teout
    real(r8), pointer, dimension(:,:) :: tini
    real(r8), pointer, dimension(:,:) :: qini
    real(r8), pointer, dimension(:,:) :: cldliqini
    real(r8), pointer, dimension(:,:) :: cldiceini
    real(r8), pointer, dimension(:,:) :: dtcore

    real(r8), pointer, dimension(:,:,:) :: fracis  ! fraction of transported species that are insoluble

    ! convective precipitation variables
    real(r8),pointer :: prec_dp(:)                ! total precipitation from ZM convection
    real(r8),pointer :: snow_dp(:)                ! snow from ZM convection
    real(r8),pointer :: prec_sh(:)                ! total precipitation from Hack convection
    real(r8),pointer :: snow_sh(:)                ! snow from Hack convection

    ! carma precipitation variables
    real(r8) :: prec_sed_carma(pcols)          ! total precip from cloud sedimentation (CARMA)
    real(r8) :: snow_sed_carma(pcols)          ! snow from cloud ice sedimentation (CARMA)

    ! stratiform precipitation variables
    real(r8),pointer :: prec_str(:)    ! sfc flux of precip from stratiform (m/s)
    real(r8),pointer :: snow_str(:)     ! sfc flux of snow from stratiform   (m/s)
    real(r8),pointer :: prec_pcw(:)     ! total precip from prognostic cloud scheme
    real(r8),pointer :: snow_pcw(:)     ! snow from prognostic cloud scheme
    real(r8),pointer :: prec_sed(:)     ! total precip from cloud sedimentation
    real(r8),pointer :: snow_sed(:)     ! snow from cloud ice sedimentation

    ! energy checking variables
    real(r8) :: zero(pcols)                    ! array of zeros
    real(r8) :: rliq(pcols)                    ! vertical integral of liquid not yet in q(ixcldliq)
    real(r8) :: rliq2(pcols)                   ! vertical integral of liquid from shallow scheme
    real(r8) :: det_s  (pcols)                 ! vertical integral of detrained static energy from ice
    real(r8) :: det_ice(pcols)                 ! vertical integral of detrained ice
    real(r8) :: flx_cnd(pcols)
    real(r8) :: flx_heat(pcols)
    type(check_tracers_data):: tracerint             ! energy integrals and cummulative boundary fluxes
    real(r8) :: zero_tracers(pcols,pcnst)

    logical   :: lq(pcnst)

    ! +++ ycw
    integer :: ideep(pcols)
    ! --- ycw

    !  pass macro to micro
    character(len=16) :: microp_scheme 
    character(len=16) :: macrop_scheme

    ! Debug physics_state.
    logical :: state_debug_checks

    call phys_getopts( microp_scheme_out      = microp_scheme, &
                       macrop_scheme_out      = macrop_scheme, &
                       state_debug_checks_out = state_debug_checks)
    
    !-----------------------------------------------------------------------
    call t_startf('bc_init')

    zero = 0._r8
    zero_tracers(:,:) = 0._r8

    lchnk = state%lchnk
    ncol  = state%ncol

    rtdt = 1._r8/ztodt

    nstep = get_nstep()


    ! Associate pointers with physics buffer fields
    itim = pbuf_old_tim_idx()
    ifld = pbuf_get_index('CLD')
    call pbuf_get_field(pbuf, ifld, cld, (/1,1,itim/),(/pcols,pver,1/))

    call pbuf_get_field(pbuf, teout_idx, teout, (/1,itim/), (/pcols,1/))

    call pbuf_get_field(pbuf, tini_idx, tini)
    call pbuf_get_field(pbuf, qini_idx, qini)
    call pbuf_get_field(pbuf, cldliqini_idx, cldliqini)
    call pbuf_get_field(pbuf, cldiceini_idx, cldiceini)

    ifld   =  pbuf_get_index('DTCORE')
    call pbuf_get_field(pbuf, ifld, dtcore, start=(/1,1,itim/), kount=(/pcols,pver,1/) )

    ifld    = pbuf_get_index('FRACIS')
    call pbuf_get_field(pbuf, ifld, fracis, start=(/1,1,1/), kount=(/pcols, pver, pcnst/)  )

    ! Set physics tendencies to 0
    tend %dTdt(:ncol,:pver)  = 0._r8
    tend %dudt(:ncol,:pver)  = 0._r8
    tend %dvdt(:ncol,:pver)  = 0._r8

    !
    ! Make sure that input tracers are all positive (otherwise,
    ! clybry_fam_adj will crash every few years).
    !

#ifdef MODAL_AERO
    call qneg3_modalx1( &
         'TPHYSBCb',lchnk  ,ncol    ,pcols   ,pver    , &
         1, pcnst, qmin  ,state%q, qneg3_worst_thresh_amode )
#else
    call qneg3('TPHYSBCb',lchnk  ,ncol    ,pcols   ,pver    , &
         1, pcnst, qmin  ,state%q )
#endif

    ! Validate state coming from the dynamics.
    if (state_debug_checks) &
         call physics_state_check(state, name="before tphysbc (dycore?)")

    call clybry_fam_adj( ncol, lchnk, map2chm, state%q, pbuf )

    ! Since clybry_fam_adj operates directly on the tracers, and has no
    ! physics_update call, re-run qneg3.

#ifdef MODAL_AERO
    call qneg3_modalx1( &
         'TPHYSBCc',lchnk  ,ncol    ,pcols   ,pver    , &
         1, pcnst, qmin  ,state%q, qneg3_worst_thresh_amode )
#else
    call qneg3('TPHYSBCc',lchnk  ,ncol    ,pcols   ,pver    , &
         1, pcnst, qmin  ,state%q )
#endif

    ! Validate output of clybry_fam_adj.
    if (state_debug_checks) &
         call physics_state_check(state, name="clybry_fam_adj")

    fracis (:ncol,:,1:pcnst) = 1._r8
    !
    ! Dump out "before physics" state
    !
    call diag_state_b4_phys_write (state)

    ! compute mass integrals of input tracers state
    call check_tracers_init(state, tracerint)

    call t_stopf('bc_init')

    !===================================================
    ! Global mean total energy fixer
    !===================================================
    call t_startf('energy_fixer')

    !*** BAB's FV heating kludge *** save the initial temperature
    tini(:ncol,:pver) = state%t(:ncol,:pver)
    if (dycore_is('LR')) then
       call check_energy_fix(state, ptend, nstep, flx_heat)
       call physics_update(state, ptend, ztodt, tend)
       call check_energy_chng(state, tend, "chkengyfix", nstep, ztodt, zero, zero, zero, flx_heat)
    end if
    ! Save state for convective tendency calculations.
    call diag_conv_tend_ini(state, pbuf)

    call cnst_get_ind('CLDLIQ', ixcldliq)
    call cnst_get_ind('CLDICE', ixcldice)
    qini     (:ncol,:pver) = state%q(:ncol,:pver,       1)
    cldliqini(:ncol,:pver) = state%q(:ncol,:pver,ixcldliq)
    cldiceini(:ncol,:pver) = state%q(:ncol,:pver,ixcldice)


    call outfld('TEOUT', teout       , pcols, lchnk   )
    call outfld('TEINP', state%te_ini, pcols, lchnk   )
    call outfld('TEFIX', state%te_cur, pcols, lchnk   )

    ! set and output the dse change due to dynpkg
    if( nstep > pbuf_times-1 ) then
       do k = 1,pver
          dtcore(:ncol,k) = (tini(:ncol,k) - dtcore(:ncol,k))/(ztodt) + tend%dTdt(:ncol,k)
       end do
       call outfld( 'DTCORE', dtcore, pcols, lchnk )
    end if

    call t_stopf('energy_fixer')
    !
    !===================================================
    ! Dry adjustment
    ! This code block is not a good example of interfacing a parameterization
    !===================================================
    call t_startf('dry_adjustment')

    ! Copy state info for input to dadadj
    ! This is a kludge, so that dadadj does not have to be correctly reformulated in dry static energy

    lq(:) = .FALSE.
    lq(1) = .TRUE.
    call physics_ptend_init(ptend, state%psetcols, 'dadadj', ls=.true., lq=lq)
    ptend%s(:ncol,:pver)   = state%t(:ncol,:pver)
    ptend%q(:ncol,:pver,1) = state%q(:ncol,:pver,1)

    call dadadj (lchnk, ncol, state%pmid,  state%pint,  state%pdel,  &
         ptend%s, ptend%q(1,1,1))
    ptend%s(:ncol,:)   = (ptend%s(:ncol,:)   - state%t(:ncol,:)  )/ztodt * cpair
    ptend%q(:ncol,:,1) = (ptend%q(:ncol,:,1) - state%q(:ncol,:,1))/ztodt
    call physics_update(state, ptend, ztodt, tend)

    call t_stopf('dry_adjustment')
    !
    !===================================================
    ! Moist convection
    !===================================================
    call t_startf('moist_convection')
    !
    ! Since the PBL doesn't pass constituent perturbations, they
    ! are zeroed here for input to the moist convection routine
    !
    call t_startf ('convect_deep_tend')
    call convect_deep_tend(  &
         cmfmc,      cmfcme,             &
         dlf,        pflx,    zdu,       &
         rliq,    &
         ztodt,   &
    ! +++ ycw
         state,   ptend, cam_in%landfrac, pbuf, ideep)
    !     state,   ptend, cam_in%landfrac, pbuf)  ! default
    ! --- ycw
    call t_stopf('convect_deep_tend')

    call physics_update(state, ptend, ztodt, tend)

    call pbuf_get_field(pbuf, prec_dp_idx, prec_dp )
    call pbuf_get_field(pbuf, snow_dp_idx, snow_dp )
    call pbuf_get_field(pbuf, prec_sh_idx, prec_sh )
    call pbuf_get_field(pbuf, snow_sh_idx, snow_sh )
    call pbuf_get_field(pbuf, prec_str_idx, prec_str )
    call pbuf_get_field(pbuf, snow_str_idx, snow_str )
    call pbuf_get_field(pbuf, prec_sed_idx, prec_sed )
    call pbuf_get_field(pbuf, snow_sed_idx, snow_sed )

    ! Check energy integrals, including "reserved liquid"
    flx_cnd(:ncol) = prec_dp(:ncol) + rliq(:ncol)
    call check_energy_chng(state, tend, "convect_deep", nstep, ztodt, zero, flx_cnd, snow_dp, zero)

    !
    ! Call Hack (1994) convection scheme to deal with shallow/mid-level convection
    !
    call t_startf ('convect_shallow_tend')

    ! +++ Yi-Chi : Oct 2014 ++++
    !call convect_shallow_tend (ztodt   , cmfmc,  cmfmc2  ,&
    !     dlf        , dlf2   ,  rliq   , rliq2, & 
    !     state      , ptend  ,  pbuf)
    call convect_shallow_tend (ztodt   , cmfmc,  cmfmc2  ,&
         dlf        , dlf2   ,  rliq   , rliq2, &
         state      , ptend  ,  pbuf, & !)
! cjshiu add three additional parameters needed by HP2011
         cam_in%landfrac, cam_in%shf, cam_in%cflx, ideep)
!cjshiu addEND
    ! ++++++++++++++++++++++++
    call t_stopf ('convect_shallow_tend')

    call physics_update(state, ptend, ztodt, tend)

    flx_cnd(:ncol) = prec_sh(:ncol) + rliq2(:ncol)
    call check_energy_chng(state, tend, "convect_shallow", nstep, ztodt, zero, flx_cnd, snow_sh, zero)

    call check_tracers_chng(state, tracerint, "convect_shallow", nstep, ztodt, zero_tracers)

    call t_stopf('moist_convection')

    ! Rebin the 4-bin version of sea salt into bins for coarse and accumulation
    ! modes that correspond to the available optics data.  This is only necessary
    ! for CAM-RT.  But it's done here so that the microphysics code which is called
    ! from the stratiform interface has access to the same aerosols as the radiation
    ! code.
    call sslt_rebin_adv(pbuf,  state)
    
    !===================================================
    ! Calculate tendencies from CARMA bin microphysics.
    !===================================================
    !
    ! If CARMA is doing detrainment, then on output, rliq no longer represents water reserved
    ! for detrainment, but instead represents potential snow fall. The mass and number of the
    ! snow are stored in the physics buffer and will be incorporated by the MG microphysics.
    !
    ! Currently CARMA cloud microphysics is only supported with the MG microphysics.
    call t_startf('carma_timestep_tend')

    if (carma_do_cldice .or. carma_do_cldliq) then
       call carma_timestep_tend(state, cam_in, cam_out, ptend, ztodt, pbuf, dlf=dlf, rliq=rliq, &
            prec_str=prec_str, snow_str=snow_str, prec_sed=prec_sed_carma, snow_sed=snow_sed_carma)
       call physics_update(state, ptend, ztodt, tend)

       ! Before the detrainment, the reserved condensate is all liquid, but if CARMA is doing
       ! detrainment, then the reserved condensate is snow.
       if (carma_do_detrain) then
          call check_energy_chng(state, tend, "carma_tend", nstep, ztodt, zero, prec_str+rliq, snow_str+rliq, zero)
       else
          call check_energy_chng(state, tend, "carma_tend", nstep, ztodt, zero, prec_str, snow_str, zero)
       end if
    end if

    call t_stopf('carma_timestep_tend')

    if( microp_scheme == 'RK' ) then

       !===================================================
       ! Calculate stratiform tendency (sedimentation, detrain, cloud fraction and microphysics )
       !===================================================
       call t_startf('stratiform_tend')

       call stratiform_tend(state, ptend, pbuf, ztodt, &
            cam_in%icefrac, cam_in%landfrac, cam_in%ocnfrac, &
            landm, cam_in%snowhland, & ! sediment
            dlf, dlf2, & ! detrain
            rliq  , & ! check energy after detrain
            cmfmc,   cmfmc2, &
            cam_in%ts,      cam_in%sst,        zdu)

       call physics_update(state, ptend, ztodt, tend)
       call check_energy_chng(state, tend, "cldwat_tend", nstep, ztodt, zero, prec_str, snow_str, zero)

       call t_stopf('stratiform_tend')

    elseif( microp_scheme == 'MG' ) then

       !===================================================
       ! Calculate macrophysical tendency (sedimentation, detrain, cloud fraction)
       !===================================================

       call t_startf('macrop_tend')

       ! don't call Park macrophysics if CLUBB is called
       if (macrop_scheme .ne. 'CLUBB_SGS') then

          call macrop_driver_tend(state, ptend, ztodt, &
               cam_in%landfrac, cam_in%ocnfrac, &
               cam_in%snowhland, & ! sediment
               dlf, dlf2, & ! detrain
               cmfmc,   cmfmc2, &
               cam_in%ts,      cam_in%sst, zdu,  pbuf, &
               det_s, det_ice)

          !  Since we "added" the reserved liquid back in this routine, we need 
	  !    to account for it in the energy checker
          flx_cnd(:ncol) = -1._r8*rliq(:ncol) 
	  flx_heat(:ncol) = det_s(:ncol)
          
          call physics_update(state, ptend, ztodt, tend)
          call check_energy_chng(state, tend, "macrop_tend", nstep, ztodt, zero, flx_cnd, det_ice, flx_heat)
       
       else ! Calculate CLUBB macrophysics

          ! =====================================================
          !    CLUBB call (PBL, shallow convection, macrophysics)
          ! =====================================================  
   
          call clubb_tend_cam(state,ptend,pbuf,1.0_r8*ztodt,&
             cmfmc, cmfmc2, cam_in, sgh30, dlf, det_s, det_ice)

          !  Since we "added" the reserved liquid back in this routine, we need 
	  !    to account for it in the energy checker
          flx_cnd(:ncol) = -1._r8*rliq(:ncol) 
	  flx_heat(:ncol) = cam_in%shf(:ncol) + det_s(:ncol)

          !    Update physics tendencies and copy state to state_eq, because that is 
          !      input for microphysics              
          call physics_update(state, ptend, ztodt, tend)
          call check_energy_chng(state, tend, "clubb_tend", nstep, ztodt, cam_in%lhf/latvap, flx_cnd, det_ice, flx_heat)
 
       endif 

       call t_stopf('macrop_tend') 

       !===================================================
       ! Calculate cloud microphysics 
       !===================================================

       call t_startf('microp_aero_run')
       call microp_aero_run(state, ptend_aero, ztodt, pbuf)
       call t_stopf('microp_aero_run')

       call t_startf('microp_tend')

       call microp_driver_tend( &
            state, ptend, ztodt, pbuf)

       ! combine aero and micro tendencies
       call physics_ptend_sum(ptend_aero, ptend, ncol)

       call physics_update(state, ptend, ztodt, tend)
       call check_energy_chng(state, tend, "microp_tend", nstep, ztodt, zero, prec_str, snow_str, zero)

       call physics_ptend_dealloc(ptend_aero)
       call t_stopf('microp_tend')

    endif

    ! Add the precipitation from CARMA to the precipitation from stratiform.
    if (carma_do_cldice .or. carma_do_cldliq) then
       prec_sed(:ncol) = prec_sed(:ncol) + prec_sed_carma(:ncol)
       snow_sed(:ncol) = snow_sed(:ncol) + snow_sed_carma(:ncol)
    end if

    if ( .not. deep_scheme_does_scav_trans() ) then

       !===================================================
       !  Aerosol wet chemistry determines scavenging fractions, and transformations
       !
       !
       !  Then do convective transport of all trace species except water vapor and
       !     cloud liquid and ice (we needed to do the scavenging first
       !     to determine the interstitial fraction) 
       !===================================================

       call t_startf('bc_aerosols')
       call aerosol_wet_intr (state, ptend, ztodt, pbuf,  cam_out, dlf)
       call physics_update(state, ptend, ztodt, tend)

       if (carma_do_wetdep) then
          ! CARMA wet deposition
          !
          ! NOTE: It needs to follow aerosol_drydep_intr, so that cam_out%xxxwetxxx
          ! fields have already been set for CAM aerosols and cam_out can be added
          ! to for CARMA aerosols.
          call t_startf ('carma_wetdep_tend')
          call carma_wetdep_tend(state, ptend, ztodt, pbuf, dlf, cam_out)
          call physics_update(state, ptend, ztodt, tend)
          call t_stopf ('carma_wetdep_tend')
       end if

       call t_startf ('convect_deep_tend2')
       call convect_deep_tend_2( state,   ptend,  ztodt,  pbuf ) 
       call t_stopf ('convect_deep_tend2')

       call physics_update(state, ptend, ztodt, tend)

       ! check tracer integrals
       call check_tracers_chng(state, tracerint, "cmfmca", nstep, ztodt,  zero_tracers)

       call t_stopf('bc_aerosols')

   endif

    !===================================================
    ! Moist physical parameteriztions complete: 
    ! send dynamical variables, and derived variables to history file
    !===================================================

    call t_startf('bc_history_write')
    call diag_phys_writeout(state, cam_out%psl)
    call diag_conv(state, ztodt, pbuf)

    call t_stopf('bc_history_write')

    !===================================================
    ! Write cloud diagnostics on history file
    !===================================================

    call t_startf('bc_cld_diag_history_write')

    call cloud_diagnostics_calc(state, pbuf)

    call t_stopf('bc_cld_diag_history_write')

    !===================================================
    ! Radiation computations
    !===================================================
    call t_startf('radiation')

    call radiation_tend(state,ptend, pbuf, &
         cam_out, cam_in, &
         cam_in%landfrac,landm,cam_in%icefrac, cam_in%snowhland, &
         fsns,    fsnt, flns,    flnt,  &
         fsds, net_flx)

    ! Set net flux used by spectral dycores
    do i=1,ncol
       tend%flx_net(i) = net_flx(i)
    end do
    call physics_update(state, ptend, ztodt, tend)
    call check_energy_chng(state, tend, "radheat", nstep, ztodt, zero, zero, zero, net_flx)

    call t_stopf('radiation')

    ! Diagnose the location of the tropopause and its location to the history file(s).
    call t_startf('tropopause')
    call tropopause_output(state)
    call t_stopf('tropopause')

    ! Save atmospheric fields to force surface models
    call t_startf('cam_export')
    call cam_export (state,cam_out,pbuf)
    call t_stopf('cam_export')

    ! Write export state to history file
    call t_startf('diag_export')
    call diag_export(cam_out)
    call t_stopf('diag_export')

end subroutine tphysbc

subroutine phys_timestep_init(phys_state, cam_out, pbuf2d)
!-----------------------------------------------------------------------------------
!
! Purpose: The place for parameterizations to call per timestep initializations.
!          Generally this is used to update time interpolated fields from boundary
!          datasets.
!
!-----------------------------------------------------------------------------------
  use shr_kind_mod,        only: r8 => shr_kind_r8
  use chemistry,           only: chem_timestep_init
  use chem_surfvals,       only: chem_surfvals_set
  use physics_types,       only: physics_state
  use physics_buffer,      only: physics_buffer_desc
  use carma_intr,          only: carma_timestep_init
  use ghg_data,            only: ghg_data_timestep_init
  use cam3_aero_data,      only: cam3_aero_data_on, cam3_aero_data_timestep_init
  use cam3_ozone_data,     only: cam3_ozone_data_on, cam3_ozone_data_timestep_init
  use radiation,           only: radiation_do
  use tracers,             only: tracers_timestep_init
  use aoa_tracers,         only: aoa_tracers_timestep_init
  use vertical_diffusion,  only: vertical_diffusion_ts_init
  use radheat,             only: radheat_timestep_init
  use solar_data,          only: solar_data_advance
  use efield,              only: get_efield
#if ( defined WACCM_PHYS )
  use iondrag,             only: do_waccm_ions
  use qbo,                 only: qbo_timestep_init
#endif
  use perf_mod

  use prescribed_ozone,    only: prescribed_ozone_adv
  use prescribed_ghg,      only: prescribed_ghg_adv
  use prescribed_aero,     only: prescribed_aero_adv
  use aerodep_flx,         only: aerodep_flx_adv
  use aircraft_emit,       only: aircraft_emit_adv
  use prescribed_volcaero, only: prescribed_volcaero_adv
  ! +++ Nudging +++ !
  use nudging,             only: Nudge_Model,nudging_timestep_init
  ! --------------- !

  implicit none

  type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
  type(cam_out_t),     intent(inout), dimension(begchunk:endchunk) :: cam_out
  
  type(physics_buffer_desc), pointer                 :: pbuf2d(:,:)

  !-----------------------------------------------------------------------------

  ! Chemistry surface values
  call chem_surfvals_set(phys_state)

  ! Solar irradiance
  call solar_data_advance()

  ! Time interpolate for chemistry.
  call chem_timestep_init(phys_state, pbuf2d)

  ! Prescribed tracers
  call prescribed_ozone_adv(phys_state, pbuf2d)
  call prescribed_ghg_adv(phys_state, pbuf2d)
  call prescribed_aero_adv(phys_state, pbuf2d)
  call aircraft_emit_adv(phys_state, pbuf2d)
  call prescribed_volcaero_adv(phys_state, pbuf2d)

  ! prescribed aerosol deposition fluxes
  call aerodep_flx_adv(phys_state, pbuf2d, cam_out)

  ! CAM3 prescribed aerosol masses
  if (cam3_aero_data_on) call cam3_aero_data_timestep_init(pbuf2d,  phys_state)

  ! CAM3 prescribed ozone data
  if (cam3_ozone_data_on) call cam3_ozone_data_timestep_init(pbuf2d,  phys_state)

  ! Time interpolate data models of gasses in pbuf2d
  call ghg_data_timestep_init(pbuf2d,  phys_state)

  ! Upper atmosphere radiative processes
  call radheat_timestep_init(phys_state, pbuf2d)
 
  ! Time interpolate for vertical diffusion upper boundary condition
  call vertical_diffusion_ts_init(pbuf2d, phys_state)

#if ( defined WACCM_PHYS )
  if (do_waccm_ions) then
     ! Compute the electric field
     call t_startf ('efield')
     call get_efield
     call t_stopf ('efield')
  endif
  !----------------------------------------------------------------------
  ! update QBO data for this time step
  !----------------------------------------------------------------------
  call qbo_timestep_init
#endif

  call carma_timestep_init()

  ! Time interpolate for tracers, if appropriate
  call tracers_timestep_init(phys_state)

  ! age of air tracers
  call aoa_tracers_timestep_init(phys_state)

  ! +++ Nudging code +++ !
  ! Update Nudging values, if needed
  !----------------------------------
  if(Nudge_Model) call nudging_timestep_init(phys_state)
  ! -------------------- !

end subroutine phys_timestep_init

end module physpkg

end module nudging
