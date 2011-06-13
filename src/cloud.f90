!###############################################################################
!!
!!  NAME
!!  Program "CLOUD"
!!
!!
!!  VERSION 
!!  $Revision: 1.44 $
!!
!!
!!  PURPOSE 
!!  Calculation of radiative characteristics for optically thick layers.
!!  Absorbing case.
!!
!!
!!  PROCEDURE 
!!
!!
!!  ADDITIONAL REQUIREMENTS 
!!  Configuration file "cloud.cfg" for definition of initial parameters.
!!
!!
!!  ATTENTION
!!  If used within SLALOM, several code clocks within the main routine
!!  (program CLOUD) have to be commented out. Look for "SLALOM" to see what to
!!  do. In addition, "use CLOUD_Defs" has to be replaced by "use SLALOM_Defs".
!!
!!
!!  CALLING SEQUENCE 
!!  CLOUD 
!!
!!
!!  COMMENT 
!!  - none 
!!
!!
!!  LANGUAGE 
!!  FORTRAN 90/95
!!
!!
!!  COMPILER 
!!  Intel Visual FORTRAN MS Windows v8.0
!!  Intel FORTRAN Linux v8.0
!!
!!
!!  CVS INFORMATION 
!!  $Id: cloud.f90,v 1.44 2008/05/28 06:41:28 tnauss Exp $ 
!!
!!
!!  CONTACT 
!!  Please send any comments, suggestions, criticism, or (for our sake) bug
!!  reports to nauss@lcrs.de
!!
!! 
!!  Development
!!  2003-2008 
!!  Thomas Nauss, Alexander A. Kokhanovsky
!!
!! 
!!  This program is free software without any warranty and without even the
!!  implied warranty of merchantability or fitness for a particular purpose. 
!!  If you use the software you must acknowledge the software and its authors. 
!!
!!
!!  HISTORY 
!!  $Log: cloud.f90,v $
!!  Revision 1.44  2008/05/28 06:41:28  tnauss
!!  Update
!!  Include pointers for LUTKSSAGrid and LUTRdInfSSAGrid.
!!
!!  Revision 1.43  2008/05/20 10:35:55  tnauss
!!  Update
!!  Delete some variables no longer necessary.
!!
!!  Revision 1.42  2008/05/19 15:12:18  tnauss
!!  Update
!!  Include documentation and adjust some lines.
!!
!!  Revision 1.41  2008/04/15 18:29:17  tnauss
!!  Update
!!  Minor modifications related to the incusion of pointers for the use
!!  within SLALOM with no consequences for the accuracy of the results.
!!
!!  Revision 1.40  2008/04/15 10:21:13  tnauss
!!  Update
!!  Include pointers for variables needed for different wavelengths withing
!!  slalom.
!!
!!  Revision 1.39  2008/04/14 10:50:43  tnauss
!!  Update:
!!  Include 856 micrometer LUTs.
!!
!!  Revision 1.38  2007/07/06 16:30:12  tnauss
!!  Major update and bug fix:
!!  Include handling of new LUTs for RInf and corresponding values of g.
!!  Actual asymmetry parameter will no longer be used for the conversion of
!!  LUTs for RInf.
!!
!!  Revision 1.37  2006/12/11 13:34:35  tnauss
!!  Update and Bugfix:
!!  Fix error and change LUT interpolation routines.
!!
!!  Revision 1.36  2006/05/31 13:47:51  tnauss
!!  Update
!!  Change LUT_Rinf grid to ssa steps of 0.001 between ssa 0.90 and 1.00.
!!
!!  Revision 1.35  2006/04/23 09:38:27  tnauss
!!  Update
!!  Include option for user defined assymetry parameter.
!!
!!  Revision 1.34  2006/04/14 11:57:35  tnauss
!!  Update
!!  Rename asymp to CLOUD.
!!
!!  Revision 1.33  2006/04/13 14:46:49  tnauss
!!  Update
!!  Document some more stuff.
!!
!!  Revision 1.32  2006/03/27 17:38:59  tnauss
!!  Revision
!!  Implement LUT-interpolation with respect to similarity parameter instead
!!  of actual ssa in order to correct for different values of the asymmetry
!!  parameter since LUTs have been computed for g=0.8500.
!!
!!  Revision 1.31  2006/03/27 14:40:15  tnauss
!!  Update
!!  Document some stuff.
!!
!!  Revision 1.30  2006/03/26 10:11:23  tnauss
!!  Update
!!  Include computation for semi-infinite media.
!!
!!  Revision 1.29  2006/03/24 21:45:34  tnauss
!!  Update
!!  Include LUTs for RInf for zenith angles between 0 and 89 degree.
!!
!!  Revision 1.28  2006/03/23 16:41:13  tnauss
!!  Update
!!  Include ice cloud LUT for Rinf.
!!
!!  Revision 1.27  2006/03/20 12:59:09  tnauss
!!  Update
!!  Update computation of reflection function for 0.999<=ssa<1.000.
!!
!!  Revision 1.26  2006/03/20 12:08:49  tnauss
!!  Update
!!  Change data types in subroutine CLOUD_1DInterpolation.
!!
!!  Revision 1.25  2006/03/20 09:41:34  tnauss
!!  Update
!!  Include new routine for reflection function at ssa=0.9990
!!
!!  Revision 1.24  2006/03/19 11:02:56  tnauss
!!  Bugfix
!!  Fix bug at computation for 0.999 < ssa < 1.000.
!!
!!  Revision 1.23  2006/03/18 15:23:50  tnauss
!!  Bugfix and update
!!  Fix computation for reflection function at ssa=1.
!!  Include new routine for reflectoin function for 0.999<=ssa<1.000.
!!
!!  Revision 1.22  2006/03/17 12:36:34  tnauss
!!  Update
!!  Include 650 micrometer LUT.
!!
!!  Revision 1.21  2006/03/16 09:43:17  tnauss
!!  Revision
!!  Change code for compilation with gnu 95.
!!
!!  Revision 1.20  2006/03/15 08:22:09  tnauss
!!  Update
!!  Disable command line configuration.
!!
!!  Revision 1.19  2006/03/13 22:01:27  tnauss
!!  Update
!!  Include initial SACURA equations for ssa larger 0.990.
!!
!!  Revision 1.18  2006/03/13 19:39:30  tnauss
!!  Update
!!  Store everything in one file.
!!
!!  Revision 1.17  2006/03/13 19:20:01  tnauss
!!  Update
!!  Include command line argument configuration.
!!
!!  Revision 1.16  2006/01/14 16:29:52  tnauss
!!  Update
!!  Interpolation of R_infinite is now computed with respect to sqrt(1-ssa).
!!
!!  Revision 1.15  2006/01/14 16:12:32  tnauss
!!  Update
!!  Include print/pause optinons in cfg file.
!!
!!  Revision 1.14  2006/01/13 17:05:23  tnauss
!!  Major update
!!  Split CLOUD.f90 into several files - one for each subroutine/module.
!!  Include 4D interpolation for LUT values of R_infinite.
!!
!!  Revision 1.13  2006/01/04 15:37:26  TNauss
!!  Bug fix
!!  Fix bug related to 3D interpolation of R_infinite.
!!
!!  Revision 1.12  2005/12/31 16:08:19  TNauss
!!  Update
!!  Include LUT interpolation for plane albedo of semi-infinite
!!  layer.
!!
!!  Revision 1.11  2005/12/23 16:30:00  TNauss
!!  Update
!!  Define all variables explicit
!!
!!  Revision 1.10  2005/12/23 10:44:42  TNauss
!!  Bug fix
!!  Delete second copy of subroutine CLOUD_ReadSettings
!!
!!  Revision 1.9  2005/12/22 12:57:15  tnauss
!!  Bug fix.
!!  Change Mu to Mu0 in some equations.
!!
!!  Revision 1.8  2005/12/22 12:07:23  tnauss
!!  Modification
!!  Change of output format.
!!
!!  Revision 1.7  2005/12/22 11:55:50  tnauss
!!  Update
!!  Multi-angle computation is now possible.
!!
!!  Revision 1.6  2005/12/21 16:30:43  tnauss
!!  Update
!!  Include multi-angle computation.
!!
!!  Revision 1.5  2005/12/17 15:21:06  tnauss
!!  Bugfix
!!  2D interpolation for mu had no valid lut array for boundary value
!!  computation.
!!
!!  Revision 1.4  2005/12/08 23:24:46  tnauss
!!  Update
!!  Include LUT for R_infinite(szen,pzen,relazm,ssa).
!!  Include 3D interpolation of LUT values of R_infinite(szen,pzen,relazm).
!!  Include input parameter configuration via configuration file cloud.cfg.
!!
!!  Revision 1.3  2005/12/08 18:35:20  tnauss
!!  Bug fix
!!  Fix bug for the case that given values for the 2D combination equal gridded
!!  values.
!!
!!  Revision 1.2  2005/12/08 14:40:54  tnauss
!!  Bug fix
!!  Check for NaN and division by 0.0
!!
!!  Revision 1.1.1.1  2005/12/08 08:58:07  tnauss
!!  Calculation of radiative characteristics for optically thick layers.
!!  Absorbing case.
!! 
!! 
!###############################################################################

!###############################################################################
!! 
!!  NAME 
!!  Module "CLOUD_Defs" 
!! 
!! 
!!  PURPOSE 
!!  Define global variables. 
!! 
!! 
!!  CONTACT 
!!  Please send any comments, suggestions, criticism, or (for our sake) bug
!!  reports to nauss@lcrs.de
!!
!! 
!!  Development
!!  2003-2008 
!!  Thomas Nauss, Alexander A. Kokhanovsky
!!
!! 
!!  This program is free software without any warranty and without even the
!!  implied warranty of merchantability or fitness for a particular purpose. 
!!  If you use the software you must acknowledge the software and its authors. 
!! 
!###############################################################################
 
    module  CLOUD_Defs
 
!*******************************************************************************
! 
!   Declaration of variables 
! 
!*******************************************************************************

    logical(4) :: bPrint        !! Print intermediate results
    logical(4) :: bPause        !! Hold after printing intermediate results
    logical(4) :: bFinite       !! Compute finite (1) or semi-infinite (0) media

    character(20) :: chVersion    !! Version of the program
    character(300) :: chControlFile !! Name of the control file
    character(300) :: chLUTPath !! Path to LUT directory

    character(50), dimension(10) :: rgchOutFile !! Names of output files

    integer(2) :: iDummy1       !! Dummy
    integer(2) :: iDummy2       !! Dummy
    integer(2) :: iDummy3       !! Dummy
    integer(2) :: iMinPosMu0    !! Temporary variable for 2D interpolation
    integer(2) :: iMinPosMu     !! Temporary variable for 2D interpolation
    integer(2) :: iMinPosSSA    !! Temporary variable for 2D interpolation
    integer(2) :: iSmalerPosMu0 !! Position of next smaler value in mu0 grid
    integer(2) :: iLargerPosMu0 !! Position of next larger value in mu grid
    integer(2) :: iSmalerPosMu  !! Position of next smaler value in mu grid
    integer(2) :: iLargerPosMu  !! Position of next larger value in mu grid
    integer(2) :: iSmalerPosSSA !! Position of next smaler value in ssa grid
    integer(2) :: iLargerPosSSA !! Position of next larger value in ssa grid
    integer(2) :: iStatus       !! Return status of subroutines

    integer(4) :: liArgument     !! Return position of the command line argument
    integer(4) :: liTauMin       !! Same as fTauMin*100 (for do routines)
    integer(4) :: liTauMax       !! Same as fTauMax*100 (for do routines)
    integer(4) :: liTauStep      !! Same as fTauStep*100 (for do routines)
    integer(4) :: liTauCounter   !! Counter for tau
    integer(4) :: liSSAMin       !! Same as fSSAMin*100000 (for do routines)
    integer(4) :: liSSAMax       !! Same as fSSAMax*100000 (for do routines)
    integer(4) :: liSSAStep      !! Same as fSSAStep*100000 (for do routines)
    integer(4) :: liSSACounter   !! Counter for SSA
    integer(4) :: liSZenMin      !! Same as fSZenMin*100 (for do routines)
    integer(4) :: liSZenMax      !! Same as fSZenMax*100 (for do routines)
    integer(4) :: liSZenStep     !! Same as fSZenStep*100 (for do routines)
    integer(4) :: liSZenCounter  !! Counter for sun zenith angles
    integer(4) :: liPZenMin      !! Same as fPZenMin*100 (for do routines)
    integer(4) :: liPZenMax      !! Same as fPZenMax*100 (for do routines)
    integer(4) :: liPZenStep     !! Same as fPZenStep*100 (for do routines)
    integer(4) :: liPZenCounter  !! Counter for pixel zenith angles
    integer(4) :: liRelAzmMin    !! Same as fRelAzmMin*100 (for do routines)
    integer(4) :: liRelAzmMax    !! Same as fRelAzmMax*100 (for do routines)
    integer(4) :: liRelAzmStep   !! Same as fRelAzmStep*100 (for do routines)
    integer(4) :: liRelAzmCounter !! Counter for relative azimuth angles

    real(4) :: fMu0          !! Actual value of cos(sun_zenith)=mu0
    real(4) :: fMu           !! Actual value of cos(pixel_zenith)=mu
    real(4) :: fSSA          !! Actual value of single-scattering albedo (ssa)
    real(4) :: fSSALook      !! Actual value of ssa to be interpolated from LUT
    real(4) :: fSmalerMu0    !! Next smaler grid value to actual mu0
    real(4) :: fLargerMu0    !! Next larger grid value to actual mu0
    real(4) :: fSmalerMu     !! Next smaler grid value to actual mu
    real(4) :: fLargerMu     !! Next larger grid value to actual mu
    real(4) :: fSmalerSSA    !! Next smaler grid value to actual ssa
    real(4) :: fLargerSSA    !! Next larger grid value to actual ssa
    real(4) :: ftestrun      !! Temporary variable for 2D interpolation
    real(4) :: fEscape0      !! Escape function e(mu0)=K/2
    real(4) :: fEscape0_1000 !! Escape function e(mu0)=K/2 for ssa = 1.000
    real(4) :: fEscape       !! Escape function e(mu)=K/2
    real(4) :: fEscape_1000  !! Escape function e(mu)=K/2 for ssa = 1.000
    real(4) :: fSZen         !! Sun zenith angle
    real(4) :: fPZen         !! Pixel zenith angle
    real(4) :: fRelAzm       !! Relative azimuth angle
    real(4) :: fSZenMin      !! Minimum sun zenith angle for computation
    real(4) :: fSZenMax      !! Maximum sun zenith angle for computation
    real(4) :: fSZenStep     !! Computation step for sun zenith angle
    real(4) :: fPZenMin      !! Minimum pixel zenith angle for computation
    real(4) :: fPZenMax      !! Maximum pixel zenith angle for computation
    real(4) :: fPZenStep     !! Computation step for pixel zenith angle
    real(4) :: fRelAzmMin    !! Minimum relative azimuth angle for computation
    real(4) :: fRelAzmMax    !! Maximum relative azimuth angle for computation
    real(4) :: fRelAzmStep   !! Computation step for relative azimuth angle
    real(4) :: fTDiff        !! Diffuse transmittance
    real(4) :: ft            !! Global reflectance
    real(4) :: fRefl         !! Reflectance
    real(4) :: fADiff        !! Absorptance 
    real(4) :: fRDiff        !! Plane albedo
    real(4) :: fBeta         !! Probability of absorption
    real(4) :: fTauMin       !! Minimum tau for computation
    real(4) :: fTauMax       !! Maximum tau for computation
    real(4) :: fTauStep      !! Computation step for tau
    real(4) :: fSSAMin       !! Minimum SSA for computation
    real(4) :: fSSAMax       !! Maximum SSA for computation
    real(4) :: fSSAStep      !! Computation step for SSA
    real(4) :: frs           !! Spherical albedo
    real(4) :: frsinf        !! Spherical albedo of a semi-infinite layer
    real(4) :: fTrans        !! Transmission
    real(4) :: fTau          !! Cloud optical thickness
    real(4) :: fg            !! Asymmetry parameter
    real(4) :: fs            !! Similarity parameter
    real(4) :: fk            !! Parameter (see publication)
    real(4) :: fl            !! Parameter (see publication)
    real(4) :: fm            !! Parameter (see publication)
    real(4) :: fn            !! Parameter (see publication)
    real(4) :: fh            !! Parameter (see publication)
    real(4) :: fa1           !! Parameter (see publication)
    real(4) :: fa2           !! Parameter (see publication)
    real(4) :: fa3           !! Parameter (see publication)
    real(4) :: fa4           !! Parameter (see publication)
    real(4) :: fanew         !! Parameter (see publication)
    real(4) :: fastra        !! Parameter (see publication)
    real(4) :: fP            !! Parameter (see publication)
    real(4) :: fx            !! Parameter (see publication)
    real(4) :: fz            !! Parameter (see publication)
    real(4) :: fy            !! Parameter (see publication)

    real(4),dimension(1) :: rgfLUTKMu0Grid(30) !! Grid (mu0) of lut_k
    real(4),dimension(2) :: rgfLUTKSSA(30,8)   !! Grid (mu0,ssa) from lut_k
    real(4),dimension(1) :: rgfLUTKMuGrid(30)  !! Grid (mu) of lut_k
    real(4),dimension(1) :: rgfLUTRdInfMu0Grid(100) !! Grid (mu0) of lut_RdInf
    real(4),dimension(2) :: rgfLUTRdInfSSA(100,9) !! Grid (mu0,ssa) lut_RdInf

    real(8) :: dRelativeMu0 !! Relative distance actual and gridded mu0 values
    real(8) :: dRelativeMu  !! Relative distance actual and gridded mu values
    real(8) :: dRelativeSSA !! Relative distance actual and gridded ssa values
    real(8) :: dx00         !! Temporary variable for 2D interpolation (c/r)
    real(8) :: dx10         !! Temporary variable for 2D interpolation (c/r)
    real(8) :: dx01         !! Temporary variable for 2D interpolation (c/r)
    real(8) :: dx11         !! Temporary variable for 2D interpolation (c/r)
    real(8) :: dK           !! Temporary variable for 2D interpolation (c/r)
    real(8) :: dRinf        !! Actual R_infinite
    real(8) :: dRinf_1000   !! Actual R_infinite for ssa = 1.000
    real(8) :: dRDInf       !! Plane albedo of semi-infinite layer
    real(8) :: dRDInf_1000  !! Plane albedo of semi-infinite layer for ssa=1.0

!   Targets and pointer (this backdoor is necessary for SLALOM)
    character(300),target ::  chLUT_RInf   !! Name of the LUT for Rinf01
    character(300),pointer :: pochLUT_RInf !! Pointer to chLUT_Rinf
    integer(4),target  :: rgiLUTRInf(90,90,181,30) !! LUT for R_inf (mu0,mu,fi)
    integer(4),pointer :: porgiLUTRInf(:,:,:,:)    !! Pointer to rgiLUTRInf
    real(4),target ::  fgLUT   !! Asymmetry parameter used computing LUT RInf
    real(4),pointer :: pofgLUT !! Pointer to fgLUT
    real(4),target ::  rgfLUTRInfSSAGrid(30)  !! Grid (ssa) of lut_RInf
    real(4),pointer :: porgfLUTRInfSSAGrid(:) !! Pointer to rgfLUTRInfSSAGrid
    real(4),target :: rgfLUTKSSAGrid(8)  !! Grid (ssa) of lut_k
    real(4),pointer :: porgfLUTKSSAGrid(:)  !! Grid (ssa) of lut_k
    real(4),target :: rgfLUTRdInfSSAGrid(9) !! Grid (ssa) lut_RdInf
    real(4),pointer :: porgfLUTRdInfSSAGrid(:) !! Grid (ssa) lut_RdInf

!   Grid of single-scattering albedo array in LUT K
    data rgfLUTRInfSSAGrid/0.800,0.810,0.820,0.830,0.840,0.850,0.860,0.870, &
                           0.880,0.890,0.900,0.910,0.920,0.930,0.940,0.950, &
                           0.960,0.970,0.980,0.990,0.991,0.992,0.993,0.994, &
                           0.995,0.996,0.997,0.998,0.999,1.000/
    data rgfLUTKSSAGrid/0.800,0.900,0.950,0.980,0.990,0.995,0.999,1.000/ 
    data rgfLUTRdInfSSAGrid/0.500,0.600,0.700,0.800,0.900,0.950,0.990,0.999, &
                            0.9999/ 
    end module CLOUD_Defs


!###############################################################################
!! 
!!  NAME 
!!  Program "CLOUD" 
!! 
!! 
!!  VERSION 
!!  $Revision: 1.44 $ 
!! 
!! 
!!  PURPOSE 
!!  Calculation of radiative characteristics for optically thick layers.
!!  Absorbing case.
!! 
!! 
!!  PROCEDURE 
!! 
!! 
!!  CALLING SEQUENCE 
!!  CLOUD 
!! 
!! 
!!  CONTACT 
!!  Please send any comments, suggestions, criticism, or (for our sake) bug
!!  reports to nauss@lcrs.de
!!
!! 
!!  Development
!!  2003-2008 
!!  Thomas Nauss, Alexander A. Kokhanovsky
!!
!! 
!!  This program is free software without any warranty and without even the
!!  implied warranty of merchantability or fitness for a particular purpose. 
!!  If you use the software you must acknowledge the software and its authors. 
!! 
!###############################################################################

    program CLOUD

!############# UNCOMMENT THE NEXT LINE IF USED WITHIN SLALOM ###################
!    subroutine CLOUD

!*******************************************************************************
!
!   Declaration of variables 
!
!*******************************************************************************

    use CLOUD_Defs

    implicit none

!############# COMMENT THE FOLLOWING SECTION IF USED WITHIN SLALOM #############

!*******************************************************************************
!
!   Format
!
!*******************************************************************************

950 format(a8,f15.7)    !! Print final results
960 format(' szen    pzen  relazm  ssa     g       ' &
           'tau    Refl   Trans  RDiff  TDiff  rs     t       ' &
           'dRinf  dRDInf fEsc0  fEsc')
961 format(3f7.2,f11.7,f8.5,f7.2,6f7.4,4f7.4)
970 format(' szen    pzen  relazm  ssa     g       ' &
           'tau    Refl   RDiff  rs     Trans  TDiff  t')
971 format(3f7.2,f11.7,f8.5,f7.2,6f7.4)

!*******************************************************************************
!
!   Say hello...
!
!*******************************************************************************

    chVersion = '$Revision: 1.44 $'
    print*, ' '
    print*, ' '
    print*, 'Program CLOUD'
    print*, trim(chVersion)

!*******************************************************************************
!
!   Read settings
!
!*******************************************************************************

    chControlFile    = 'cloud.cfg'
    call CLOUD_ReadSettings

!*******************************************************************************
!
!   Set pointers
!
!*******************************************************************************

    pochLUT_RInf => chLUT_RInf
    porgiLUTRInf => rgiLUTRInf
    pofgLUT => fgLUT
    porgfLUTRInfSSAGrid => rgfLUTRInfSSAGrid
    porgfLUTKSSAGrid => rgfLUTKSSAGrid
    porgfLUTRdInfSSAGrid => rgfLUTRdInfSSAGrid

!*******************************************************************************
!
!   Define parameters
!
!*******************************************************************************

    chControlFile    = 'cloud.cfg'
    call CLOUD_ReadSettings

!   Set asymmetry parameter used for computation of actual LUT for RInf
   if(pochLUT_RInf.eq.'lut_RInf_0645_aef_06.dat') then
     pofgLUT = 0.850018000000000000  ! 06µm, 645.5nm
     rgchOutFile(1) = 'solar_lut_0645_06.dat'
    elseif(pochLUT_RInf.eq.'lut_RInf_0645_aef_10.dat') then
     pofgLUT = 0.861757000000000000  ! 10µm, 645.5nm
     rgchOutFile(1) = 'solar_lut_0645_10.dat'
    elseif(pochLUT_RInf.eq.'lut_RInf_0645_aef_16.dat') then
     pofgLUT = 0.869221000000000000  ! 16µm, 645.5nm
     rgchOutFile(1) = 'solar_lut_0645_16.dat'
    elseif(pochLUT_RInf.eq.'lut_RInf_0856_aef_06.dat') then
     pofgLUT = 0.843722000000000000  ! 06µm, 856.5nm
     rgchOutFile(1) = 'solar_lut_0856_06.dat'
    elseif(pochLUT_RInf.eq.'lut_RInf_0856_aef_10.dat') then
     pofgLUT = 0.857995000000000000  ! 10µm, 856.5nm
     rgchOutFile(1) = 'solar_lut_0856_10.dat'
    elseif(pochLUT_RInf.eq.'lut_RInf_0856_aef_16.dat') then
     pofgLUT = 0.867233000000000000  ! 16µm, 856.5nm
     rgchOutFile(1) = 'solar_lut_0856_16.dat'
    elseif(pochLUT_RInf.eq.'lut_RInf_1630_aef_06.dat') then
     pofgLUT = 0.817538000000000000  ! 06µm, 1630nm
     rgchOutFile(1) = 'solar_lut_1630_06.dat'
    elseif(pochLUT_RInf.eq.'lut_RInf_1630_aef_10.dat') then
     pofgLUT = 0.846061000000000000  ! 10µm, 1630nm
     rgchOutFile(1) = 'solar_lut_1630_10.dat'
    elseif(pochLUT_RInf.eq.'lut_RInf_1630_aef_16.dat') then
     pofgLUT = 0.861364000000000000  ! 16µm, 1630nm
     rgchOutFile(1) = 'solar_lut_1630_16.dat'
    elseif(pochLUT_RInf.eq.'lut_RInf_ice.dat') then
     pofgLUT = 0.752400000000000000  ! all aef
     rgchOutFile(1) = 'solar_lut_ice.dat'
    else
     print*, 'No valid LUT for RInf has been selected.'
     print*, 'The program is going to stop...'
     stop
    endif

!   Open output filenames    
    open(550,file=trim(rgchOutFile(1)))
    write(550,*) 'Version: ',trim(chVersion)
    write(550,970)

!*******************************************************************************
!
!   Read LUTs
!
!*******************************************************************************

    call CLOUD_ReadLUT

!*******************************************************************************
!
!   Compute radiative characteristics
!
!*******************************************************************************

!   Program controll settings
    liTauMin = nint(fTauMin*100)
    liTauMax = nint(fTauMax*100)
    liTauStep = nint(fTauStep*100)
    liSSAMin = nint(fSSAMin * 1000000)
    liSSAMax = nint(fSSAMax * 1000000)
    liSSAStep = nint(fSSAStep * 1000000)
    liSZenMin = nint(fSZenMin*100)
    liSZenMax = nint(fSZenMax*100)
    liSZenStep = nint(fSZenStep*100)
    liPZenMin = nint(fPZenMin*100)
    liPZenMax = nint(fPZenMax*100)
    liPZenStep = nint(fPZenStep*100)
    liRelAzmMin = nint(fRelAzmMin*100)
    liRelAzmMax = nint(fRelAzmMax*100)
    liRelAzmStep = nint(fRelAzmStep*100)

!   Convert ssa grid values of the LUTs to similarity parameter (1-fs is used
!   for historical reasons since then, the coloumn with the smalest/largest ssa
!   value is the same).
!   LUTs have been computed for an asymmetry parameter of 0.8500, this value is
!   used for the conversion.

    porgfLUTRInfSSAGrid  = 1.0 - sqrt( (1.0-porgfLUTRInfSSAGrid)/ &
                           (1.0-porgfLUTRInfSSAGrid*pofgLUT) )
    porgfLUTKSSAGrid  = 1.0 - sqrt( (1.0-porgfLUTKSSAGrid)/ &
                           (1.0-porgfLUTKSSAGrid*pofgLUT) )
    porgfLUTRdinfSSAGrid  = 1.0 - sqrt( (1.0-porgfLUTRdinfSSAGrid)/ &
                           (1.0-porgfLUTRdinfSSAGrid*pofgLUT) )

    do 10 liTauCounter = liTauMin, liTauMax, liTauStep
     do 20 liSZenCounter = liSZenMin, liSZenMax, liSZenStep
      do 30 liPZenCounter = liPZenMin, liPZenMax, liPZenStep
       do 40 liRelAzmCounter = liRelAzmMin, liRelAzmMax, liRelAzmStep
        do 50 liSSACounter = liSSAMin, liSSAMax, liSSAStep

        fTau = float(liTauCounter)/100.0
        fSSA = float(liSSACounter)/1000000.0
        fSZen = float(liSZenCounter)/100.0
        fPZen = float(liPZenCounter)/100.0 
        fRelAzm = float(liRelAzmCounter)/100.0

!############# END OF SECTION TO COMMENT IF USED WITHIN SLALOM #################

        fMu0 = cos( fSZen*acos(-1.)/180.)
        fMu  = cos( fPZen*acos(-1.)/180.)

!       Start computation
        fBeta=1.-fSSA
        fs=sqrt( fBeta/(1.-fSSA*fg))
        fastra=sqrt(3.)*fs-(0.985*fs*fs-0.253*fs*fs*fs)/(6.464-5.464*fs)
        fk=(1.-fSSA*fg)*fastra
        fl=(1.-fs)*(1.-0.681*fs)/(1.+0.792*fs)
        fanew=(1.+1.8*fs-7.087*fs*fs+4.74*fs*fs*fs)/(1.-0.819*fs)
        fm=(1.+1.537*fs)*alog(fanew/(1.-fs)/(1.-fs) )
        fn=sqrt(  (1.-fs)*(1.+0.414*fs)/(1.+1.888*fs)  )
        frsinf=(1.-fs)*(1.-0.319*fs)/(1.+1.17*fs)
        fx=fk*fTau
        fy=4*sqrt((1-fSSA)/(3*(1-fg)))

!       Compute actual R_infinite with respect to mu0,mu,fi,ssa
!       Compute R_infinite additionaly for ssa = 1.000 if actual SSA is greater/
!       equal 0.999 (see computation of reflection function later).
!       Use 1 - similarity parameter for the interpolation instead of ssa in
!       order to correct for actual asymmetry parameter.
        if(fSSA.ge.0.999999) then
         fSSALook = 1.000
         fSSALook = 1.0 - sqrt( (1.0-fSSALook)/(1.0-fSSALook*fg) )
         call CLOUD_RInfinite
         dRinf_1000 = dRinf
        else
         fSSALook = fSSA
         fSSALook = 1.0 - sqrt( (1.0-fSSALook)/(1.0-fSSALook*fg) )
         call CLOUD_RInfinite
        endif

        fz=fMu0-0.5
        fa1=-0.9991+3.139*fg-1.874*fg*fg
        fa2=1.435-4.924*fg+2.089*fg*fg
        fa3=0.719-5.801*fg+2.117*fg*fg
        fa4=-0.509+0.418*fg+3.36*fg*fg
        fP=(fa1*fz+fa2*fz*fz)*fs+(fa3*fz+fa4*fz*fz)*fs*fs

!       Compute actual plane albedo for a semi-infinite layer with respect to 
!       mu0 dRDInf=(1.-fs)*exp(fP)/(1.+2.*fs*fMu0)
!       Compute actual plane albedo additionaly for ssa = 1.000 if actual SSA
!       is greater/equal 0.999 (see computation of reflection function later).
!       Use 1 - similarity parameter for the interpolation instead of ssa in
!       order to correct for actual asymmetry parameter.
        if(fSSA.ge.0.999999) then
         fSSALook = 1.000
         fSSALook = 1.0 - sqrt( (1.0-fSSALook)/(1.0-fSSALook*fg) )
         call CLOUD_RdInfinite
         dRDInf_1000 = dRDInf
        else
         fSSALook = fSSA
         fSSALook = 1.0 - sqrt( (1.0-fSSALook)/(1.0-fSSALook*fg) )
         call CLOUD_RdInfinite
        endif

!       Radiative characteristics
!       Spherical
!       This if is necessary since otherwise some environments produce a
!       runtime error
        if(fx.ne.0.00) then
         if(.not.bFinite  ) then
          ft=0.e0
         else
          ft=fm*fn*fn*exp(-fx)/(1.-fl*fl*exp(-2.*fx))
         endif
         frs=frsinf-ft*fl*exp(-fx)
        endif

!       Compute escape function
!       Hemispherical - 2D interpolation of LUT values is required.
!       Compute escape function additionaly for ssa = 1.000 if actual SSA is
!       greater/equal 0.999 (see computation of reflection function later).
!       Use 1 - similarity parameter for the interpolation instead of ssa in
!       order to correct for actual asymmetry parameter.
        if(fSSA.ge.0.999999) then
         fSSALook = 1.000
         fSSALook = 1.0 - sqrt( (1.0-fSSALook)/(1.0-fSSALook*fg) )
         call CLOUD_Escape
         fEscape0_1000 = fEscape0
         fEscape_1000  = fEscape
        else
         fSSALook = fSSA
         fSSALook = 1.0 - sqrt( (1.0-fSSALook)/(1.0-fSSALook*fg) )
         call CLOUD_Escape
        endif

!       Diffuse transmittance
        fTDiff=ft*fEscape0/fn

!       Plane albedo
        fRDiff=dRDInf-ft*fEscape0*fl*exp(-fx)/fn

!       Transmission
        if(fSSA.eq.1.000) then
         fTrans=ft*fEscape*fEscape0/fn/fn
        elseif(fSSA.ge.0.999999) then
         fTrans = sinh(fy) / sinh(1.07*fy+fx)
        else
         fTrans=ft*fEscape*fEscape0/fn/fn
        endif

!       Reflection function and absorptance
        if(fSSA.eq.1.000) then
         fh=1.07079+0.75*(1.-fg)*fTau
         if(.not.bFinite  ) then
          ft=0.0e0
         else
          ft=1./fh
         endif
         frs=1-ft
         fTDiff=ft*fEscape0
         fRDiff=1-fTDiff
         fTrans=ft*fEscape*fEscape0 
         fRefl=dRinf-fTrans
         fADiff=0.0
        elseif(fSSA.ge.0.999999) then
!        fff1= fEscape0 at SSA=1.0
!        fff2= fEscape  at SSA=1.0
!        fff3= dRinf    at SSA=1.0
!        dRinf=fff3-fy*fff1*fff2
!!       dRinf=dRinf_1000-fy*fEscape0_1000*fEscape_1000
       dRinf=dRinf_1000*exp(-fy*fEscape0_1000*fEscape_1000/dRinf_1000) 
       fRefl=dRinf-ft*fl*exp(-fx)*fEscape*fEscape0/fn/fn 
         if(.not.bFinite  ) then
          ft=0.0e0
         else
          ft=sinh(fy) / sinh(1.07079*fy+fx)
         endif
         fRefl=dRInf-ft*exp(-fx-fy)*fEscape0_1000*fEscape_1000 
         fADiff = 1 - fTDiff - fRDiff 
        else
         fRefl=dRinf-ft*fl*exp(-fx)*fEscape*fEscape0/fn/fn
         fADiff = 1 - fTDiff - fRDiff
        endif

!*******************************************************************************
!
!   Print/write output
!
!*******************************************************************************

!############# COMMENT THE FOLLOWING SECTION IF USED WITHIN SLALOM #############

     if(.not.bFinite  ) then
      write(550,971) fSZen, fPZen, fRelAzm, fSSA, fg, &
                     fTau, fRefl, fRDiff, frs, fTrans, fTdiff, ft
     else
      write(550,961) fSZen, fPZen, fRelAzm, fSSA, fg, &
                     fTau, fRefl, fTrans, fRDiff, fTdiff, frs,ft, &
                     dRinf, dRDInf, fEscape0, fEscape
     endif


     if(bPrint) then
         print*,' '
         print*,'Final results from program CLOUD:'
         print*,'tau:           ', fTau
         print*,'rs:            ', frs
         print*,'t:             ', ft
         print*,'RDiff:         ', fRDiff
         print*,'TDiff:         ', fTDiff
         print*,'ADiff:         ', fADiff
         print*,'Trans:         ', fTrans
         print*,'Refl:          ', fRefl
         print*,'Ssa:           ', fSSA
         print*,'fg:            ', fg
         print*,'fSZen:         ', fSZen
         print*,'fPZen:         ', fPZen
         print*,'fRelAzm:       ', fRelAzm
         if(bPause) read(*,*)
        endif

50      continue
40     continue
30    continue
20   continue
10  continue

    close(550)
    close(551)

    stop
    end program CLOUD

!############# END OF SECTION TO COMMENT IF USED WITHIN SLALOM #################
!############# UNCOMMENT THE NEXT TWO LINES IF USED WITHIN SLALOM ##############
!    return
!    end subroutine CLOUD

!###############################################################################
!! 
!!  NAME 
!!  Subroutine "CLOUD_ReadSettings" 
!! 
!! 
!!  VERSION 
!!  $Revision: 1.44 $ 
!! 
!! 
!!  PURPOSE 
!!  Read settings from configuration file.
!! 
!! 
!!  PROCEDURE 
!!  Read settings from configuration file. 
!! 
!!
!!  CONTACT 
!!  Please send any comments, suggestions, criticism, or (for our sake) bug
!!  reports to nauss@lcrs.de
!!
!! 
!!  Development
!!  2003-2008 
!!  Thomas Nauss, Alexander A. Kokhanovsky
!!
!! 
!!  This program is free software without any warranty and without even the
!!  implied warranty of merchantability or fitness for a particular purpose. 
!!  If you use the software you must acknowledge the software and its authors. 
!! 
!###############################################################################
 
    subroutine CLOUD_ReadSettings
 
!*******************************************************************************
! 
!   Declaration of variables 
! 
!******************************************************************************* 

    use CLOUD_Defs
    implicit none

!*******************************************************************************
! 
!   Define namelist
! 
!*******************************************************************************

    namelist /CLOUDControl/ &
    bFinite  ,     &
    chLUTPath,     &
    chLUT_RInf,    &
    fg,            &
    fTauMin,       &
    fTauMax,       &
    fTauStep,      &
    fSSAMin,       &
    fSSAMax,       &
    fSSAStep,      &
    fSZenMin,      &
    fSZenMax,      &
    fSZenStep,     &
    fPZenMin,      &
    fPZenMax,      &
    fPZenStep,     &
    fRelAzmMin,    & 
    fRelAzmMax,    &
    fRelAzmStep,   &
    bPrint,        &
    bPause

!*******************************************************************************
!
!   Read configuration file
!
!*******************************************************************************

    open(501,file=chControlFile)
    rewind(501)
    read(501,nml=CLOUDControl)
    close(501)

    return
    end subroutine CLOUD_ReadSettings 


!###############################################################################
!! 
!!  NAME 
!!  Subroutine "CLOUD_ReadLUT" 
!! 
!! 
!!  VERSION 
!!  $Revision: 1.44 $ 
!! 
!! 
!!  PURPOSE 
!!  Read tabulated values from ASCII file.
!!  1. Read lut_k.dat for 2D interpolation of K(mu[0],ssa)
!!     K/2 equals escape function
!!  2. Read LUT interp_ssa_wwww.dat for 3D interpolation of 
!!     R(szen,pzen,relazm,ssa)
!!     "_ssa"  = single-scattering albedo (ssa=0.95 --> "_0_95")
!!     "_wwww" = wavelenght in microns (wwww=0.865microns --> "_0865")
!!  3. Read LUT lut_RdInf.dat for 1D interpolation of RDinf(szen,ssa)
!! 
!! 
!!  PROCEDURE 
!!  Read tabulated values from ASCII file. 
!! 
!!
!!  CONTACT 
!!  Please send any comments, suggestions, criticism, or (for our sake) bug
!!  reports to nauss@lcrs.de
!!
!! 
!!  Development
!!  2003-2008 
!!  Thomas Nauss, Alexander A. Kokhanovsky
!!
!! 
!!  This program is free software without any warranty and without even the
!!  implied warranty of merchantability or fitness for a particular purpose. 
!!  If you use the software you must acknowledge the software and its authors. 
!! 
!############################################################################### 
 
    subroutine  CLOUD_ReadLUT
 
!*******************************************************************************
 
! 
!   Declaration of variables 
! 
!*******************************************************************************
 
    use CLOUD_Defs
    implicit none
    
    integer(2) iCounter       !! Counter for lut_K.dat
    integer(2) iCounterSZen   !! Counter for sun zenith in lut_RInf.dat
    integer(2) iCounterPZen   !! Counter for pixel zenith in lut_RInf.dat
    integer(2) iCounterRelAzm !! Counter for relative azimuth in lut_RInf.dat

!*******************************************************************************
! 
!   Format 
! 
!*******************************************************************************
 
950 format(f17.8,8f17.8)  !! Read lut_K.dat content
951 format(9f10.6)        !! Print lut_K.dat content
952 format(33i5)          !! Read lut_RInf content
953 format(10f22.14)      !! Read lut_RdInf.dat content

!*******************************************************************************
! 
!   Read lut_K.dat for 2D interpolation of K(mu[0],ssa)
! 
!*******************************************************************************

!   Read LUT ktable 
    iCounter = 0 
    open(501,file=trim(chLUTPath)//'lut_k.dat') 
    read(501,*)
    do iCounter = 1, 30
     read(501,950) rgfLUTKMu0Grid(iCounter), rgfLUTKSSA(iCounter,1:8) 
    enddo 
    close(501)
    rgfLUTKMuGrid = rgfLUTKMu0Grid

!*******************************************************************************
! 
!   Read lut_RInf.dat for 4D interpolation of R(szen,pzen,fi,ssa)
! 
!*******************************************************************************

!   Read LUT values in array
    open(501,file=trim(chLUTPath)//trim(pochLUT_RInf))
    read(501,*)
     do iCounterSZen = 0, 89
      do iCounterRelAzm = 0, 180
       do iCounterPZen = 0, 89
        read(501,952) idummy1, idummy2, idummy3, &
               porgiLUTRInf(iCounterSZen+1,iCounterPZen+1,iCounterRelAzm+1,1:30)
       enddo
      enddo
     enddo
    close(501)

!*******************************************************************************
! 
!   Read lut_RdInf.dat for 1D interpolation of RDinf(szen,ssa)
! 
!*******************************************************************************

!   Read LUT ktable 
    iCounter = 0 
    open(501,file=trim(chLUTPath)//'lut_RdInf.dat') 
    read(501,*)
    do iCounter = 1,100
     read(501,953) rgfLUTRdInfMu0Grid(iCounter), rgfLUTRdInfSSA(iCounter,1:9) 
    enddo 
    close(501)

    return
    end subroutine CLOUD_ReadLUT


!############################################################################### 
!! 
!!  NAME 
!!  Subroutine "CLOUD_RInfinite" 
!! 
!! 
!!  VERSION 
!!  $Revision: 1.44 $ 
!! 
!! 
!!  PURPOSE 
!!  Compute 4D interpolation of tabulated values of 
!!  R_infinite(szen,pzen,relazm,ssa) using lut_RInf.dat.
!! 
!! 
!!  PROCEDURE 
!!  1. Compute 3D interpolation for R_infinite(szen,pzen,relazm) for LUT-grid
!!     ssa value smaler than actual ssa.
!!  2. Compute 3D interpolation for R_infinite(szen,pzen,relazm) for LUT-grid
!!     ssa value larger than actual ssa.
!!  3. Compute 1D interpolation for the two R_infinite values with respect to
!!     acutal ssa.
!! 
!!
!!  CONTACT 
!!  Please send any comments, suggestions, criticism, or (for our sake) bug
!!  reports to nauss@lcrs.de
!!
!! 
!!  Development
!!  2003-2008 
!!  Thomas Nauss, Alexander A. Kokhanovsky
!!
!! 
!!  This program is free software without any warranty and without even the
!!  implied warranty of merchantability or fitness for a particular purpose. 
!!  If you use the software you must acknowledge the software and its authors. 
!!
!############################################################################### 
 
    subroutine  CLOUD_RInfinite
 
!*******************************************************************************
! 
!   Declaration of variables 
! 
!*******************************************************************************
 
    use CLOUD_Defs
    implicit none
 
    integer(2) :: iSmalerPosSZen   !! Position of smaller szen value in array
    integer(2) :: iLargerPosSZen   !! Position of larger szen value in array
    integer(2) :: iSmalerPosPZen   !! Position of smaller pzen value in array
    integer(2) :: iLargerPosPZen   !! Position of larger pzen value in array
    integer(2) :: iSmalerPosRelAzm !! Position of smaller relazm value in array
    integer(2) :: iLargerPosRelAzm !! Position of larger relazm value in array

    real(8) :: dR     !! Parameter for 3D interpolation
    real(8) :: dS     !! Parameter for 3D interpolation
    real(8) :: dT     !! Parameter for 3D interpolation
    real(8) :: dX000  !! Parameter for 3D interpolation
    real(8) :: dX001  !! Parameter for 3D interpolation
    real(8) :: dX010  !! Parameter for 3D interpolation
    real(8) :: dX011  !! Parameter for 3D interpolation
    real(8) :: dX100  !! Parameter for 3D interpolation
    real(8) :: dX101  !! Parameter for 3D interpolation
    real(8) :: dX110  !! Parameter for 3D interpolation
    real(8) :: dX111  !! Parameter for 3D interpolation
    real(8) :: dSmalerRInf !! RInf for smaler LUT-grid ssa value
    real(8) :: dLargerRInf !! RInf for larger LUT-grid ssa value

    real(4),dimension(1) :: rgfMinDiffSSA(30) !! Parameter for 4D interp.
    real(4),dimension(1) :: rgfMinLoc(1)      !! Parameter for 4D interp.

!*******************************************************************************
! 
!   Find position of LUT-grid ssa value closest to actual ssa
! 
!*******************************************************************************
 
    rgfMinDiffSSA = abs(porgfLUTRInfSSAGrid - fSSALook) 
    rgfMinLoc = MinLoc(rgfMinDiffSSA) 
    iMinPosSSA = rgfMinLoc(1)
    if(porgfLUTRInfSSAGrid(iMinPosSSA).eq.fSSALook) then
     iSmalerPosSSA = iMinPosSSA
     iLargerPosSSA = iMinPosSSA
    elseif(porgfLUTRInfSSAGrid(iMinPosSSA).gt.fSSALook) then
     iSmalerPosSSA = iMinPosSSA-1
     iLargerPosSSA = iMinPosSSA
    else
     iSmalerPosSSA = iMinPosSSA
     iLargerPosSSA = iMinPosSSA+1
    endif

    if(iSmalerPosSSA.le.0) iSmalerPosSSA = 1
    if(iLargerPosSSA.gt.30) iLargerPosSSA = 30

!   Set nearest smaler/larger grid value of ssa
    fSmalerSSA = porgfLUTRInfSSAGrid(iSmalerPosSSA) 
    fLargerSSA = porgfLUTRInfSSAGrid(iLargerPosSSA) 

!*******************************************************************************
! 
!   Interpolate value of R_infinite with respect to szen,pzen,relazm for
!   smaler ssa
! 
!*******************************************************************************
 
!    fSZen = acosd(fMu0)
!    fPZen = acosd(fMu)
    dR = fSZen - int(fSZen)
    dS = fPZen - int(fPZen)
    dT = fRelAzm - int(fRelAzm)
    iSmalerPosSZen = int(fSZen)+1
    iLargerPosSZen = int(fSZen)+2
    iSmalerPosPZen = int(fPZen)+1
    iLargerPosPZen = int(fPZen)+2
    iSmalerPosRelAzm = int(fRelAzm)+1
    iLargerPosRelAzm = int(fRelAzm)+2
    if(iLargerPosSZen.ge.91) iLargerPosSZen = 90
    if(iLargerPosPZen.ge.91) iLargerPosSZen = 90
    if(iLargerPosRelAzm.ge.182) iLargerPosRelAzm = 2
    dX000 = &
    float(porgiLUTRInf(iSmalerPosSZen,iSmalerPosPZen, &
          iSmalerPosRelAzm,iSmalerPosSSA))/1000.0
    dX100 = &
    float(porgiLUTRInf(iLargerPosSZen,iSmalerPosPZen, &
          iSmalerPosRelAzm,iSmalerPosSSA))/1000.0
    dX010 = &
    float(porgiLUTRInf(iSmalerPosSZen,iLargerPosPZen, &
          iSmalerPosRelAzm,iSmalerPosSSA))/1000.0
    dX110 = &
    float(porgiLUTRInf(iLargerPosSZen,iLargerPosPZen, &
          iSmalerPosRelAzm,iSmalerPosSSA))/1000.0
    dX001 = &
    float(porgiLUTRInf(iSmalerPosSZen,iSmalerPosPZen, &
          iLargerPosRelAzm,iSmalerPosSSA))/1000.0
    dX101 = &
    float(porgiLUTRInf(iLargerPosSZen,iSmalerPosPZen, &
          iLargerPosRelAzm,iSmalerPosSSA))/1000.0
    dX011 = &
    float(porgiLUTRInf(iSmalerPosSZen,iLargerPosPZen, &
          iLargerPosRelAzm,iSmalerPosSSA))/1000.0
    dX111 = &
    float(porgiLUTRInf(iLargerPosSZen,iLargerPosPZen, &
          iLargerPosRelAzm,iSmalerPosSSA))/1000.0
    call CLOUD_3DInterpolation &
         (dR,dS,dT,dX000,dX001,dX010,dX011,dX100,dX101,dX110,dX111, dSmalerRInf)

!*******************************************************************************
! 
!   Interpolate value of R_infinite with respect to szen,pzen,relazm for
!   larger ssa
! 
!*******************************************************************************
 
    if(iSmalerPosSSA.eq.iLargerPosSSA) then
     dLargerRInf = dSmalerRInf
    else
     dX000 = &
     float(porgiLUTRInf(iSmalerPosSZen,iSmalerPosPZen, &
           iSmalerPosRelAzm,iLargerPosSSA))/1000.0
     dX100 = &
     float(porgiLUTRInf(iLargerPosSZen,iSmalerPosPZen, &
           iSmalerPosRelAzm,iLargerPosSSA))/1000.0
     dX010 = &
     float(porgiLUTRInf(iSmalerPosSZen,iLargerPosPZen, &
           iSmalerPosRelAzm,iLargerPosSSA))/1000.0
     dX110 = &
     float(porgiLUTRInf(iLargerPosSZen,iLargerPosPZen, &
           iSmalerPosRelAzm,iLargerPosSSA))/1000.0
     dX001 = &
     float(porgiLUTRInf(iSmalerPosSZen,iSmalerPosPZen, &
           iLargerPosRelAzm,iLargerPosSSA))/1000.0
     dX101 = &
     float(porgiLUTRInf(iLargerPosSZen,iSmalerPosPZen, &
           iLargerPosRelAzm,iLargerPosSSA))/1000.0
     dX011 = &
     float(porgiLUTRInf(iSmalerPosSZen,iLargerPosPZen, &
           iLargerPosRelAzm,iLargerPosSSA))/1000.0
     dX111 = &
     float(porgiLUTRInf(iLargerPosSZen,iLargerPosPZen, &
           iLargerPosRelAzm,iLargerPosSSA))/1000.0
     call CLOUD_3DInterpolation &
          (dR,dS,dT,dX000,dX001,dX010,dX011,dX100,dX101,dX110,dX111,dLargerRInf)
    endif

!*******************************************************************************
! 
!   Interpolate values of R_infinite for smaler/larger ssa with respect to
!   actual ssa
! 
!*******************************************************************************
 
   if(iSmalerPosSSA.eq.iLargerPosSSA) then
     dRinf = dSmalerRInf
    else
!    Use sqrt(1-ssa) instead of ssa since R=R_no_abs-sqrt(1-ssa)*const
!    Keep in mind: larger SSA and therefore larger Rinf value computed above
!    corresponds to smaler sqrt of beta.
     dRelativeSSA = (fSSALook - fSmalerSSA) / (fLargerSSA - fSmalerSSA) 
     call CLOUD_1DInterpolation(dRelativeSSA,dSmalerRInf,dLargerRInf,dRinf)
    endif

!   Print intermediate results
    if(bPrint) then
     print*,' '
     print*,'Results from subroutine CLOUD_RInfinite'
     print*,'SSA:           ', fSSALook
     print*,'SZen:          ', fSZen
     print*,'PZen:          ', fPZen
     print*,'RelAzm:        ', fRelAzm
     print*,'Smaler SSA:    ', fSmalerSSA
     print*,'Larger SSA:    ', fLargerSSA 
     print*,'Smaler SZen:   ', iSmalerPosSZen-1
     print*,'Larger SZen:   ', iLargerPosSZen-1
     print*,'Smaler PZen:   ', iSmalerPosPZen-1
     print*,'Larger PZen:   ', iLargerPosPZen-1
     print*,'Smaler RelAzm: ', iSmalerPosRelAzm-1
     print*,'Larger RelAzm: ', iLargerPosRelAzm-1
     print*,'dR:            ', dR
     print*,'dS:            ', dS
     print*,'dT:            ', dT
     print*,'Smaler Rinf:   ', dSmalerRInf
     print*,'LargerRinf:    ', dLargerRInf
     print*,'Rinf:          ', dRinf
     if(bPause) read(*,*)
    endif
    
    return
    end subroutine CLOUD_RInfinite


!############################################################################### 
!! 
!!  NAME 
!!  Subroutine "CLOUD_RdInfinite" 
!! 
!! 
!!  VERSION 
!!  $Revision: 1.44 $ 
!! 
!! 
!!  PURPOSE 
!!  Compute 1D interpolation of tabulated values of RDinf(szen,ssa)
!!  using LUT lut_RdInf.dat. 
!! 
!! 
!!  PROCEDURE 
!!  Compute 1D interpolation.
!! 
!!
!!  CONTACT 
!!  Please send any comments, suggestions, criticism, or (for our sake) bug
!!  reports to nauss@lcrs.de
!!
!! 
!!  Development
!!  2003-2008 
!!  Thomas Nauss, Alexander A. Kokhanovsky
!!
!! 
!!  This program is free software without any warranty and without even the
!!  implied warranty of merchantability or fitness for a particular purpose. 
!!  If you use the software you must acknowledge the software and its authors. 
!!
!############################################################################### 
 
    subroutine  CLOUD_RdInfinite
 
!*******************************************************************************
! 
!   Declaration of variables 
! 
!*******************************************************************************
 
    use CLOUD_Defs
    implicit none
    real(4),dimension(1) :: rgfMinDiffMu0(100)  !! Parameter for 2D interp.
    real(4),dimension(1) :: rgfMinDiffSSA(9)    !! Parameter for 2D interp.
    real(4),dimension(1) :: rgfMinLoc(1)        !! Parameter for 2D interp.
 
!*******************************************************************************
! 
!   Interpolate value of RDinf with respect to szen, ssa
! 
!*******************************************************************************
 
!   Cosine sun zenith angle
!   Find LUT-grid mu0 value closest to mu0
    rgfMinDiffMu0 = abs(rgfLUTRdInfMu0Grid - fMu0) 
    rgfMinLoc = MinLoc(rgfMinDiffMu0) 
    iMinPosMu0 = rgfMinLoc(1) 
    if(rgfLUTRdInfMu0Grid(iMinPosMu0).eq.fMu0) then
     iSmalerPosMu0 = iMinPosMu0
     iLargerPosMu0 = iMinPosMu0
    elseif(rgfLUTRdInfMu0Grid(iMinPosMu0).gt.fMu0) then
     iSmalerPosMu0 = iMinPosMu0-1
     iLargerPosMu0 = iMinPosMu0
    else
     iSmalerPosMu0 = iMinPosMu0
     iLargerPosMu0 = iMinPosMu0+1
    endif

    if(iSmalerPosMu0.le.0) iSmalerPosMu0 = 1
    if(iLargerPosMu0.gt.100) iLargerPosMu0 = 100

!   Set nearest smaler/larger grid value of mu0
    fSmalerMu0 = rgfLUTRdInfMu0Grid(iSmalerPosMu0) 
    fLargerMu0 = rgfLUTRdInfMu0Grid(iLargerPosMu0) 

!   Compute relative distance between closest LUT-grid values and
!   the actual mu0 value   
    if(fLargerMu0.eq.fSmalerMu0) then
     dRelativeMu0 = 0.0
    else
     dRelativeMu0 = (fMu0 - fSmalerMu0) / (fLargerMu0 - fSmalerMu0) 
    endif

!   Single-scattering albedo
!   Find LUT-grid ssa value closest to ssa
    rgfMinDiffSSA = abs(porgfLUTRdinfSSAGrid - fSSALook) 
    rgfMinLoc = MinLoc(rgfMinDiffSSA) 
    iMinPosSSA = rgfMinLoc(1) 
    if(porgfLUTRdinfSSAGrid(iMinPosSSA).eq.fSSALook) then
     iSmalerPosSSA = iMinPosSSA
     iLargerPosSSA = iMinPosSSA
    elseif(porgfLUTRdinfSSAGrid(iMinPosSSA).gt.fSSALook) then
     iSmalerPosSSA = iMinPosSSA-1
     iLargerPosSSA = iMinPosSSA
    else
     iSmalerPosSSA = iMinPosSSA
     iLargerPosSSA = iMinPosSSA+1
    endif

    if(iSmalerPosSSA.le.0) iSmalerPosSSA = 1
    if(iLargerPosSSA.gt.9) iLargerPosSSA = 9

!   Set nearest smaler/larger grid value of ssa
    fSmalerSSA = porgfLUTRdinfSSAGrid(iSmalerPosSSA) 
    fLargerSSA = porgfLUTRdinfSSAGrid(iLargerPosSSA) 

!   Compute relative distance between closest LUT-grid values and
!   the actual ssa value   
    if(fLargerSSA.eq.fSmalerSSA) then
     dRelativeSSA = 0.0
    else
     dRelativeSSA = (fSSALook - fSmalerSSA) / (fLargerSSA - fSmalerSSA) 
    endif

!   Read nearest boundary values from LUT-arrays for actual mu0 and ssa
    dx00 = rgfLUTRdInfSSA(iSmalerPosMu0,iSmalerPosSSA) 
    dx01 = rgfLUTRdInfSSA(iSmalerPosMu0,iLargerPosSSA) 
    dx10 = rgfLUTRdInfSSA(iLargerPosMu0,iSmalerPosSSA) 
    dx11 = rgfLUTRdInfSSA(iLargerPosMu0,iLargerPosSSA) 

!   Compute 2D interpolation (column/row-coordinate system)
    call CLOUD_2DInterpolation &
         (dRelativeMu0,dRelativeSSA,dx00,dx01,dx10,dx11,dRDInf)

!   Print intermediate results
    if(bPrint) then
     print*,' '
     print*,'Results from subroutine CLOUD_RdInfinite'
     print*,'SSA:           ', fSSALook
     print*,'Mu0:           ', fMu0
     print*,'Smaler SSA:    ', fSmalerSSA
     print*,'Larger SSA:    ', fLargerSSA 
     print*,'Smaler Mu0:    ', fSmalerMu0
     print*,'Larger Mu0:    ', fLargerMu0
     print*,'Relative SSA:  ', dRelativeSSA
     print*,'Relative Mu0:  ', dRelativeMu0
     print*,'dx00:          ', dx00
     print*,'dx01:          ', dx01
     print*,'dx10:          ', dx10
     print*,'dx11:          ', dx11
     print*,'RDinf:         ', dRDInf
     if(bPause) read(*,*)
    endif

    return
    end subroutine CLOUD_RdInfinite


!############################################################################### 
!! 
!!  NAME 
!!  Subroutine "CLOUD_Escape" 
!! 
!! 
!!  VERSION 
!!  $Revision: 1.44 $ 
!! 
!! 
!!  PURPOSE 
!!  Compute 2D interpolation of tabulated values of K(mu0,mu)
!!  using LUT k_lut.dat. 
!! 
!! 
!!  PROCEDURE 
!!  1. Compute 2D interpolation of tabulated values for fMu0. 
!!  2. Compute 2D interpolation of tabulated values for fMu. 
!!
!!
!!  CONTACT 
!!  Please send any comments, suggestions, criticism, or (for our sake) bug
!!  reports to nauss@lcrs.de
!!
!! 
!!  Development
!!  2003-2008 
!!  Thomas Nauss, Alexander A. Kokhanovsky
!!
!! 
!!  This program is free software without any warranty and without even the
!!  implied warranty of merchantability or fitness for a particular purpose. 
!!  If you use the software you must acknowledge the software and its authors. 
!! 
!############################################################################### 
 
    subroutine  CLOUD_Escape
 
!*******************************************************************************
! 
!   Declaration of variables 
! 
!*******************************************************************************
 
    use CLOUD_Defs
    implicit none
    real(4),dimension(1) :: rgfMinDiffMu0(30)   !! Parameter for 2D interp.
    real(4),dimension(1) :: rgfMinDiffMu(30)    !! Parameter for 2D interp.
    real(4),dimension(1) :: rgfMinDiffSSA(8)    !! Parameter for 2D interp.
    real(4),dimension(1) :: rgfMinLoc(1)        !! Parameter for 2D interp.

!*******************************************************************************
! 
!   Format
! 
!*******************************************************************************

950 format(a16,f15.7)    !! Print settings and results from 2D interpolation
 
!*******************************************************************************
! 
!   Interpolate value for K with respect to mu0, ssa
! 
!*******************************************************************************

!   Cosine sun zenith angle
!   Find LUT-grid mu0 value closest to mu0
    rgfMinDiffMu0 = abs(rgfLUTKMu0Grid - fMu0) 
    rgfMinLoc = MinLoc(rgfMinDiffMu0) 
    iMinPosMu0 = rgfMinLoc(1) 
    if(rgfLUTKMu0Grid(iMinPosMu0).eq.fMu0) then
     iSmalerPosMu0 = iMinPosMu0
     iLargerPosMu0 = iMinPosMu0
    elseif(rgfLUTKMu0Grid(iMinPosMu0).gt.fMu0) then
     iSmalerPosMu0 = iMinPosMu0-1
     iLargerPosMu0 = iMinPosMu0
    else
     iSmalerPosMu0 = iMinPosMu0
     iLargerPosMu0 = iMinPosMu0+1
    endif
     
    if(iSmalerPosMu0.le.0) iSmalerPosMu0 = 1
    if(iLargerPosMu0.gt.30) iLargerPosMu0 = 30

!   Set nearest smaler/larger grid value of mu0
    fSmalerMu0 = rgfLUTKMu0Grid(iSmalerPosMu0) 
    fLargerMu0 = rgfLUTKMu0Grid(iLargerPosMu0) 

!   Compute relative distance between closest LUT-grid values and
!   the actual mu0 value   
    if(fLargerMu0.eq.fSmalerMu0) then
     dRelativeMu0 = 0.0
    else
     dRelativeMu0 = (fMu0 - fSmalerMu0) / (fLargerMu0 - fSmalerMu0) 
    endif


!   Single-scattering albedo
!   Find LUT-grid ssa value closest to ssa with ssa(grid) <= ssa
    rgfMinDiffSSA = abs(porgfLUTKSSAGrid - fSSALook) 
    rgfMinLoc = MinLoc(rgfMinDiffSSA) 
    iMinPosSSA = rgfMinLoc(1) 
    if(porgfLUTKSSAGrid(iMinPosSSA).eq.fSSALook) then
     iSmalerPosSSA = iMinPosSSA
     iLargerPosSSA = iMinPosSSA
    elseif(porgfLUTKSSAGrid(iMinPosSSA).gt.fSSALook) then
     iSmalerPosSSA = iMinPosSSA-1
     iLargerPosSSA = iMinPosSSA
    else
     iSmalerPosSSA = iMinPosSSA
     iLargerPosSSA = iMinPosSSA+1
    endif
     
    if(iSmalerPosSSA.le.0) iSmalerPosSSA = 1
    if(iLargerPosSSA.gt.8) iLargerPosSSA = 8

!   Set nearest smaler/larger grid value of ssa
    fSmalerSSA = porgfLUTKSSAGrid(iSmalerPosSSA) 
    fLargerSSA = porgfLUTKSSAGrid(iLargerPosSSA) 

!   Compute relative distance between closest LUT-grid values and
!   the actual ssa value   
    if(fLargerSSA.eq.fSmalerSSA) then
     dRelativeSSA = 0.0
    else
     dRelativeSSA = (fSSALook - fSmalerSSA) / (fLargerSSA - fSmalerSSA) 
    endif

!   Read nearest boundary values from LUT-arrays for actual mu0 and ssa
    dx00 = rgfLUTKSSA(iSmalerPosMu0,iSmalerPosSSA) 
    dx01 = rgfLUTKSSA(iSmalerPosMu0,iLargerPosSSA) 
    dx10 = rgfLUTKSSA(iLargerPosMu0,iSmalerPosSSA) 
    dx11 = rgfLUTKSSA(iLargerPosMu0,iLargerPosSSA) 
     
!   Compute 2D interpolation (column/row-coordinate system)
    call CLOUD_2DInterpolation(dRelativeMu0,dRelativeSSA,dx00,dx01,dx10,dx11,dK)
    fEscape0 = dK/2.0

!*******************************************************************************
! 
!   Interpolated value for K with respect to mu, ssa
! 
!*******************************************************************************

!   Cosine pixel zenith angle
!   Find LUT-grid mu value closest to mu with mu(grid) <= mu
    rgfMinDiffMu = abs(rgfLUTKMuGrid - fMu) 
    rgfMinLoc = MinLoc(rgfMinDiffMu) 
    iMinPosMu = rgfMinLoc(1) 
    if(rgfLUTKMuGrid(iMinPosMu).eq.fMu) then
     iSmalerPosMu = iMinPosMu
     iLargerPosMu = iMinPosMu
    elseif(rgfLUTKMuGrid(iMinPosMu).gt.fMu) then
     iSmalerPosMu = iMinPosMu-1
     iLargerPosMu = iMinPosMu
    else
     iSmalerPosMu = iMinPosMu
     iLargerPosMu = iMinPosMu+1
    endif

    if(iSmalerPosMu.le.0) iSmalerPosMu = 1
    if(iLargerPosMu.gt.30) iLargerPosMu = 30

!   Set nearest smaler/larger grid value of Mu
    fSmalerMu = rgfLUTKMuGrid(iSmalerPosMu) 
    fLargerMu = rgfLUTKMuGrid(iLargerPosMu) 

!   Set nearest smaler/larger grid value of mu
    fSmalerMu = rgfLUTKMuGrid(iSmalerPosMu) 
    fLargerMu = rgfLUTKMuGrid(iLargerPosMu) 

!   Compute relative distance between closest LUT-grid values and
!   the actual mu value   
    if(fLargerMu.eq.fSmalerMu) then
     dRelativeMu = 0.0
    else
     dRelativeMu = (fMu - fSmalerMu) / (fLargerMu - fSmalerMu) 
    endif

!   Read nearest boundary values from LUT-arrays for actual mu and ssa
    dx00 = rgfLUTKSSA(iSmalerPosMu,iSmalerPosSSA) 
    dx01 = rgfLUTKSSA(iSmalerPosMu,iLargerPosSSA) 
    dx10 = rgfLUTKSSA(iLargerPosMu,iSmalerPosSSA) 
    dx11 = rgfLUTKSSA(iLargerPosMu,iLargerPosSSA) 

!   Compute 2D interpolation (column/row-coordinate system)
    call CLOUD_2DInterpolation(dRelativeMu,dRelativeSSA,dx00,dx01,dx10,dx11,dK)
    fEscape = dK/2.0

    if(bPrint) then
     print*,' '
     print*,'Results from subroutine CLOUD_Escape'
     print*,'SSA:           ', fSSALook
     print*,'Mu0:           ', fMu0
     print*,'Mu:            ', fMu
     print*,'Smaler SSA:    ', fSmalerSSA
     print*,'Larger SSA:    ', fLargerSSA 
     print*,'Smaler Mu0:    ', fSmalerMu0
     print*,'Larger Mu0:    ', fLargerMu0
     print*,'Smaler Mu:     ', fSmalerMu
     print*,'Larger Mu:     ', fLargerMu
     print*,'Relative SSA:  ', dRelativeSSA
     print*,'Relative Mu0:  ', dRelativeMu0
     print*,'Relative Mu:   ', dRelativeMu
     print*,'Escape0:       ', fEscape0
     print*,'Escape:        ', fEscape
     if(bPause) read(*,*)
    endif

    return
    end subroutine CLOUD_Escape


!############################################################################### 
!! 
!!  NAME 
!!  Subroutine "CLOUD_1DInterpolation" 
!! 
!! 
!!  VERSION 
!!  $Revision: 1.44 $ 
!! 
!! 
!!  PURPOSE 
!!  Extends scalar endpoint data to a line. 
!! 
!! 
!!  PROCEDURE 
!!  x0-----r-----x1
!!
!! 
!!  This code is from John Burkardt
!!  http://www.csit.fsu.edu/~burkardt/f_src/f_src.html.
!!
!############################################################################### 

    subroutine CLOUD_1DInterpolation(r,x0,x1,x)
    implicit none

    real ( kind = 8 ) r
    real ( kind = 8 ) x
    real ( kind = 8 ) x0
    real ( kind = 8 ) x1

    x = ( 1.0D+00 - r ) * x0 + r * x1

    return
    end subroutine CLOUD_1DInterpolation


!############################################################################### 
!! 
!!  NAME 
!!  Subroutine "CLOUD_2DInterpolation" 
!! 
!! 
!!  VERSION 
!!  $Revision: 1.44 $ 
!! 
!! 
!!  PURPOSE 
!!  Extends scalar point data into a square. 
!! 
!! 
!!  PROCEDURE 
!!   01------------11
!!    |      .      |
!!    |      .      |
!!    |.....rs......|
!!    |      .      |
!!    |      .      |
!!   00------------10
!!
!! Formula:
!!
!! Written in terms of R and S, the map has the form:
!!
!! X(R,S) =
!!          1     * ( + x00 )
!!        + r     * ( - x00 + x10 )
!!        + s     * ( - x00       + x01 )
!!        + r * s * ( + x00 - x10 - x01 + x11 )
!!
!! Written in terms of the coefficients, the map has the form:
!!
!! X(R,S) =   x00 * ( 1 - r - s + r * s )
!!          + x01 * (         s - r * s )
!!          + x10 * (     r     - r * s )
!!          + x11 * (             r * s )
!!
!!          = x00 * ( 1 - r ) * ( 1 - s )
!!          + x01 * ( 1 - r ) *       s
!!          + x10 *       r   * ( 1 - s )
!!          + x11 *       r           s!! 
!! 
!!  This code is from John Burkardt
!!  http://www.csit.fsu.edu/~burkardt/f_src/f_src.html.
!!
!############################################################################### 

    subroutine CLOUD_2DInterpolation(r,s,x00,x01,x10,x11,x)
    implicit none

    real ( kind = 8 ) r
    real ( kind = 8 ) s
    real ( kind = 8 ) x
    real ( kind = 8 ) x00
    real ( kind = 8 ) x01
    real ( kind = 8 ) x10
    real ( kind = 8 ) x11

    x =             + x00 &
        + r *     ( - x00 + x10 ) &
        + s *     ( - x00       + x01 ) &
        + r * s * ( + x00 - x10 - x01 + x11 )

    return
    end subroutine CLOUD_2DInterpolation


!############################################################################### 
!! 
!!  NAME 
!!  Subroutine "CLOUD_3DInterpolation" 
!! 
!! 
!!  VERSION
!!  $Revision: 1.44 $
!!
!!
!!  PURPOSE
!!  Extends scalar point data into a cube.
!!
!!
!!  PROCEDURE
!!    011--------------111 
!!      |               |
!!      |               | 
!!      |               |
!!      |               |
!!      |               |
!!    001--------------101
!!
!!
!!      *---------------*
!!      |               |
!!      |               |
!!      |      rst      |
!!      |               |
!!      |               |
!!      *---------------*
!!
!!
!!    010--------------110
!!      |               |
!!      |               |
!!      |               |
!!      |               | 
!!      |               |
!!    000--------------100 
!!
!!  Formula:
!!
!!    Written as a polynomial in R, S and T, the interpolation map has the 
!!    form:
!!
!!      X(R,S,T) =
!!        1         * ( + x000 )
!!      + r         * ( - x000 + x100 )
!!      +     s     * ( - x000        + x010 )
!!      +         t * ( - x000               + x001 )
!!      + r * s     * ( + x000 - x100 - x010                       + x110 )
!!      + r     * t * ( + x000 - x100        - x001        + x101 )
!!      +     s * t * ( + x000        - x010 - x001 + x011 )
!!      + r * s * t * ( - x000 + x100 + x010 + x001 - x011 - x101 - x110 + x111 )
!! 
!!  This code is from John Burkardt
!!  http://www.csit.fsu.edu/~burkardt/f_src/f_src.html.
!!
!############################################################################### 

    subroutine CLOUD_3DInterpolation &
              (r, s, t, x000, x001, x010, x011, x100, x101, x110, x111,x)
    implicit none
    
    real(8) :: r
    real(8) :: s
    real(8) :: t
    real(8) :: x
    real(8) :: x000
    real(8) :: x001
    real(8) :: x010
    real(8) :: x011
    real(8) :: x100
    real(8) :: x101
    real(8) :: x110
    real(8) :: x111

    x = &
    1.0E+00     * ( + x000 ) &
    + r         * ( - x000 + x100 ) &
    +     s     * ( - x000        + x010 ) &
    +         t * ( - x000               + x001 ) &
    + r * s     * ( + x000 - x100 - x010                      + x110 ) &
    + r     * t * ( + x000 - x100        - x001        + x101 ) &
    +     s * t * ( + x000        - x010 - x001 + x011 ) &
    + r * s * t * ( - x000 + x100 + x010 + x001 - x011 - x101 - x110 + x111 )
 
    return
    end subroutine CLOUD_3DInterpolation

