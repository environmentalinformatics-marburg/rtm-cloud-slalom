!######################################################################################### 
!! 
!!  NAME 
!!  Program "SLALOM" 
!! 
!!  VERSION 
!!  $Revision: 1.7 $ 
!! 
!!  PURPOSE 
!!  Inverse algorithm for the retrieval of the cloud optical thickness (COT) and the
!!  cloud single scattering albedo (SSA) for optically thick layers. The forward model
!!  used within SLALOM is CLOUD.
!! 
!!  PROCEDURE 
!!  
!!  ADDITIONAL REQUIREMENTS 
!! 
!!  ATTENTION 
!!  This program consists of three basic parts:
!!  1. Read/write and control routines (all in the main sacuraICE program)
!!  2. Main retrieval and necessary functions (starting from sacuraICE_retrieve)
!!  3. Routines from program CLOUD for the forward computation part
!!  For operational or user specified applications, the retrieval algorithm can be
!!  included in other frameworks, as long as all necessary input parameters are given
!!  along with the calling sequence.
!!  For the forward computation (which is necessary for Rmes-Rcomp), all routines of
!!  the CLOUD program are included in SLALOM. The only modifications are:
!!  1. The variables are declared in a special section of the module SLALOM_defs
!!  2. The main routine is included as a subroutine (without control sequence lines) 
!!  Configuration of this version: You don't need a configuration file, even there is a
!!  routine for reading it (but it is not called). This is just for an easier copy/paste
!!  of new versions of CLOUD in this code - so don't worry about that.
!! 
!!  CALLING SEQUENCE 
!!  SLALOM
!! 
!!  COMMENT 
!! 
!!  LANGUAGE 
!!  FORTRAN 90/95 
!! 
!!  COMPILER 
!!  Intel Visual FORTRAN MS Windows
!!  Intel FORTRAN Linux
!! 
!!  CVS INFORMATION 
!!  $Id: slalomICE.f90,v 1.7 2007/07/19 10:11:09 tnauss Exp $ 
!! 
!!  CONTACT 
!!  Please send any comments, suggestions, criticism, or (for our sake) bug reports to 
!!  info@lcrs.de
!! 
!!  Copyright (c) 2003-2006 
!!  Alexander A. Kokhanovsky, Thomas Nauss
!! 
!!  This program is free software without any warranty and without even the implied  
!!  warranty of merchantability or fitness for a particular purpose. 
!!  If you use the software you must acknowledge the software and its author. 
!! 
!!  HISTORY 
!!  $Log: slalomICE.f90,v $
!!  Revision 1.7  2007/07/19 10:11:09  tnauss
!!  Major revision (still test phase)
!!  Include routines for 2nd LUT for abs. channel.
!!
!!  Revision 1.6  2007/03/22 17:08:32  tnauss
!!  Bugfix
!!  1001 things have been changed.
!!
!!  Revision 1.4  2006/10/31 16:52:34  tnauss
!!  Update
!!
!!  Revision 1.2  2006/03/25 11:56:39  tnauss
!!  Minor bug fixes.
!!
!!  Revision 1.1.1.1  2006/03/24 23:19:26  tnauss
!!  Retrieval of ice cloud parameters.
!!
!! 
!######################################################################################### 

!######################################################################################### 
!! 
!!  NAME 
!!  Module "SLALOM_defs" 
!! 
!!  VERSION 
!!  $Revision: 1.7 $ 
!!
!!  PURPOSE 
!!  Define global variables. 
!!
!!  PROCEDURE 
!!  
!!  CONTACT 
!!  Please send any comments, suggestions, criticism, or (for our sake) bug reports to 
!!  info@lcrs.de  
!! 
!!  Copyright (c) 2003-2006 
!!  Alexander A. Kokhanovsky, Thomas Nauss
!! 
!!  This program is free software without any warranty and without even the implied  
!!  warranty of merchantability or fitness for a particular purpose. 
!!  If you use the software you must acknowledge the software and its author. 
!! 
!######################################################################################### 

    module SLALOM_defs
 
!***************************************************************************************** 
! 
!   Declaration of variables for SLALOM
! 
!***************************************************************************************** 

    logical(4) :: bASCII    !! Use ASCII input data file
    logical(4) :: bWaterCloud

    character(100) :: chCommandLineArgument !! Command line argument
    character(300) :: chOutputFolder !! Path to the output folder
    character(300), dimension(5) :: rgchInFile   !! Names of input files

    integer(1),dimension(1) :: prgiCloud(2748620)   !! Input cloud mask
    integer(1),dimension(1) :: prgiMask(2748620)

    integer(2) :: iError    !! Flag for input/output error
    integer(2) :: iCloud    !! Cloud mask flag

    integer(2) :: LUT

    integer(4) :: liRuns     !! Number of input value sets
    integer(4) :: liRecl

    real(8) :: fAefInput
    real(8) :: fDummy
    real(8) :: fNonAbs       !! Actual value of the non-abs. band
    real(8) :: fAbs          !! Actual value of the abs. band
    real(8) :: fNonAbsAlb    !! Actual value of the non-abs. band background albedo
    real(8) :: fAbsAlb       !! Actual value of the abs. band background albedo
    real(8) :: fTau_1000     !! Value of tau for the non-abs. wavelength for ssa = 1.000
    real(8) :: fSSAMinValid  !! Minimum value of SSA for finding root of functionR
    real(8) :: fSSAMaxValid  !! Maximum value of SAA for finding root for functionR
    real(8) :: fAefMinValid  !! Minimum value of aef for finding root of functionR
    real(8) :: fAefMaxValid  !! Maximum value of aef for finding root for functionR
    real(8) :: fdTau         !! Ratio between retrieved and correct tau
    real(8) :: fdBeta        !! Ratio between retrieved and correct beta
    real(8) :: fWavelength1  !! Non-absorbing channel wavelength
    real(8) :: fWavelength2  !! Absorbing channel wavelength
    real(8) :: fMiIce        !! Imaginary part of refractive index for ice clouds
    real(8) :: fPi           !! Pi
    real(8) :: fBetaInf      !! Probability of absorption for an infinite cloud
    real(8) :: fl0           !! Default parameter for particle absorption length    
    real(8) :: fPAL          !! Particle absorption length
    real(8) :: fAef          !! Effective cloud droplet radius
    real(8) :: fIWP          !! Ice water path
    real(8) :: fXef          !! Size parameter
    real(8) :: fLWP          !! Liquid water path
    real(8) :: ft1           !! Global reflectance for wavelength 1
    real(8) :: fRho          !! Density of water
    real(8) :: fSigma_ext    !! Extinction coefficient
    real(8) :: fSigma_abs    !! Absorbtion coefficient
    real(8) :: dR0           !! Reflection of a semi-infinite cloud


    real(4),dimension(1) :: prgfSZen(2748620)   !! Input solar zenith angle
    real(4),dimension(1) :: prgfPZen(2748620)   !! Input pixel zenith angle
    real(4),dimension(1) :: prgfRelAzm(2748620) !! Input relative azimuth angle
    real(4),dimension(1) :: prgfNonAbs(2748620) !! Input non-abs. band reflectance
    real(4),dimension(1) :: prgfAbs(2748620)    !! Input abs. band reflectance
    real(4),dimension(1) :: prgfNonAbsAlb(2748620) !! Input non-abs. background albedo
    real(4),dimension(1) :: prgfAbsAlb(2748620) !! Input abs. background albedo
    real(4),dimension(1) :: prgfSSAIn(2748620) !! Input test ssa value
    real(4),dimension(1) :: prgfTauIn(2748620) !! Input test tau value
    real(4),dimension(1) :: prgfTauRetrieve(2748620)
    real(4),dimension(1) :: prgfSSARetrieve(2748620)
    real(4),dimension(1) :: prgfPAL(2748620)
    real(4),dimension(1) :: prgfAef(2748620)
    real(4),dimension(1) :: prgfLWP(2748620)
    real(4),dimension(1) :: rgfLUTRInfSSAGridIn(30) !! Grid values of lut_RInf.dat for ssa
    real(4),dimension(1) :: rgfLUTKSSAGridIn(8)  !! Grid values of lut_k.dat for ssa
    real(4),dimension(1) :: rgfLUTRdInfSSAGridIn(9) !! Grid values of LUT lut_RdInf.dat for ssa
    real(4),dimension(1) :: rgfAefInput(2748620)

    real(8) :: dSSARetrieve !! Retrieved value of SSA


!***************************************************************************************** 
! 
!   Declaration of variables from program CLOUD
! 
!***************************************************************************************** 

    logical(4) :: bPrint        !! Print intermediate results
    logical(4) :: bPause        !! Hold after printing intermediate results
    logical(4) :: bFinite       !! Compute results for finite (true) or semi-infinite (false) media

    character(20) :: chVersion    !! Version of the program
    character(50) :: chActualFileSelector !! Actual file selector used to extend the wildcard information
    character(300) :: chValue     !! Content of the command line argument at position iArgument
    character(300) :: chControlFile !! Name of the control file
    character(300) :: chLUT_RInf     !! Name of the LUT for R_infinite
    character(300) :: chLUT_RInf1     !! Name of the LUT for R_infinite
    character(300) :: chLUT_RInf2     !! Name of the LUT for R_infinite

    character(50), dimension(10) :: rgchOutFile !! Names of output files

    integer(2) :: iDummy1       !! Dummy
    integer(2) :: iDummy2       !! Dummy
    integer(2) :: iDummy3       !! Dummy
    integer(2) :: iMinPosMu0    !! Temporary variable necessary for 2D interpolation
    integer(2) :: iMinPosMu     !! Temporary variable necessary for 2D interpolation
    integer(2) :: iMinPosSSA    !! Temporary variable necessary for 2D interpolation
    integer(2) :: iSmalerPosMu0 !! Position of next smaler value in mu0 grid
    integer(2) :: iLargerPosMu0 !! Position of next larger value in mu grid
    integer(2) :: iSmalerPosMu  !! Position of next smaler value in mu grid
    integer(2) :: iLargerPosMu  !! Position of next larger value in mu grid
    integer(2) :: iSmalerPosSSA !! Position of next smaler value in ssa grid
    integer(2) :: iLargerPosSSA !! Position of next larger value in ssa grid
    integer(2) :: iStatus       !! Return status of subroutines

    integer(4) :: liArgument     !! Return position of the command line argument
    integer(4) :: liTauMin       !! Same as fTauMin*100 (necessary for do routines)
    integer(4) :: liTauMax       !! Same as fTauMax*100 (necessary for do routines)
    integer(4) :: liTauStep      !! Same as fTauStep*100 (necessary for do routines)
    integer(4) :: liTauCounter   !! Counter for tau
    integer(4) :: liSSAMin       !! Same as fSSAMin*100000 (necessary for do routines)
    integer(4) :: liSSAMax       !! Same as fSSAMax*100000 (necessary for do routines)
    integer(4) :: liSSAStep      !! Same as fSSAStep*100000 (necessary for do routines)
    integer(4) :: liSSACounter   !! Counter for SSA
    integer(4) :: liSZenMin      !! Same as fSZenMin*100 (necessary for do routines)
    integer(4) :: liSZenMax      !! Same as fSZenMax*100 (necessary for do routines)
    integer(4) :: liSZenStep     !! Same as fSZenStep*100 (necessary for do routines)
    integer(4) :: liSZenCounter  !! Counter for sun zenith angles
    integer(4) :: liPZenMin      !! Same as fPZenMin*100 (necessary for do routines)
    integer(4) :: liPZenMax      !! Same as fPZenMax*100 (necessary for do routines)
    integer(4) :: liPZenStep     !! Same as fPZenStep*100 (necessary for do routines)
    integer(4) :: liPZenCounter  !! Counter for pixel zenith angles
    integer(4) :: liRelAzmMin    !! Same as fRelAzmMin*100 (necessary for do routines)
    integer(4) :: liRelAzmMax    !! Same as fRelAzmMax*100 (necessary for do routines)
    integer(4) :: liRelAzmStep   !! Same as fRelAzmStep*100 (necessary for do routines)
    integer(4) :: liRelAzmCounter !! Counter for relative azimuth angles

!    integer(4),dimension(90,90,181,30) :: rgiLUTRInf !! LUT for R_infinite (mu0,mu,fi)
    integer(4),dimension(90,90,181,30) :: rgiLUTRInf1 !! LUT for R_infinite (mu0,mu,fi)
    integer(4),dimension(90,90,181,30) :: rgiLUTRInf2 !! LUT for R_infinite (mu0,mu,fi)

    real(4) :: fMu0          !! Actual value of cos(sun_zenith)=mu0
    real(4) :: fMu           !! Actual value of cos(pixel_zenith)=mu
    real(4) :: fSSA          !! Actual value of single-scattering albedo (ssa)
    real(4) :: fSSALook      !! Actual value of ssa that should be interpolated from LUT
    real(4) :: fSmalerMu0    !! Next smaler mu0 grid value with respect to actual mu
    real(4) :: fLargerMu0    !! Next larger mu0 grid value with respect to actual mu
    real(4) :: fSmalerMu     !! Next smaler mu grid value with respect to actual mu
    real(4) :: fLargerMu     !! Next larger mu grid value with respect to actual mu
    real(4) :: fSmalerSSA    !! Next smaler ssa grid value with respect to actual ssa
    real(4) :: fLargerSSA    !! Next larger ssa grid value with respect to actual ssa
    real(4) :: ftestrun         !! Temporary variable necessary for 2D interpolation
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
    real(4) :: fgLUT         !! Asymmetry parameter used for computation of LUT RInf
    real(4) :: fgLUT1        !! Asymmetry parameter used for computation of LUT RInf
    real(4) :: fgLUT2        !! Asymmetry parameter used for computation of LUT RInf
    real(4) :: fs            !! Similarity parameter
    real(4) :: fk            !! Parameter
    real(4) :: fl            !! Parameter
    real(4) :: fm            !! Parameter
    real(4) :: fn            !! Parameter
    real(4) :: fh            !! Parameter
    real(4) :: fa1           !! Parameter
    real(4) :: fa2           !! Parameter
    real(4) :: fa3           !! Parameter
    real(4) :: fa4           !! Parameter
    real(4) :: fanew         !! Parameter
    real(4) :: fastra        !! Parameter
    real(4) :: fP            !! Parameter
    real(4) :: fx            !! Parameter
    real(4) :: fz            !! Parameter
    real(4) :: fy            !! Paramter

    real(4),dimension(1) :: rgfLUTRInfSSAGrid(30) !! Grid values of lut_RInf.dat for ssa
    real(4),dimension(1) :: rgfLUTRInfSSAGrid1(30) !! Grid values of lut_RInf.dat for ssa
    real(4),dimension(1) :: rgfLUTRInfSSAGrid2(30) !! Grid values of lut_RInf.dat for ssa
    real(4),dimension(1) :: rgfLUTKMu0Grid(30) !! Grid values of lut_k.dat for mu0
    real(4),dimension(1) :: rgfLUTKSSAGrid(8)  !! Grid values of lut_k.dat for ssa
    real(4),dimension(2) :: rgfLUTKSSA(30,8)   !! Values of K(mu0,ssa) from lut_k.dat
    real(4),dimension(1) :: rgfLUTKMuGrid(30)  !! Grid values of lut_k.dat for mu
    real(4),dimension(1) :: rgfLUTRdInfMu0Grid(100) !! Grid values of LUT lut_RdInf.dat for mu0
    real(4),dimension(1) :: rgfLUTRdInfSSAGrid(9) !! Grid values of LUT lut_RdInf.dat for ssa
    real(4),dimension(2) :: rgfLUTRdInfSSA(100,9) !! Values of K(mu0,ssa) from LUT lut_RdInf.dat

    real(8) :: dRelativeMu0 !! Relative distance between actual and gridded mu0 values
    real(8) :: dRelativeMu  !! Relative distance between actual and gridded mu values
    real(8) :: dRelativeSSA !! Relative distance between actual and gridded ssa values
    real(8) :: dx00         !! Temporary variable necessary for 2D interpolation (c/r)
    real(8) :: dx10         !! Temporary variable necessary for 2D interpolation (c/r)
    real(8) :: dx01         !! Temporary variable necessary for 2D interpolation (c/r)
    real(8) :: dx11         !! Temporary variable necessary for 2D interpolation (c/r)
    real(8) :: dK           !! Temporary variable necessary for 2D interpolation (c/r)
    real(8) :: dRinf        !! Actual R_infinite
    real(8) :: dRinf_1000   !! Actual R_infinite for ssa = 1.000
    real(8) :: dRDInf       !! Plane albedo of semi-infinite layer
    real(8) :: dRDInf_1000  !! Plane albedo of semi-infinite layer for ssa = 1.000

!   Grid of single-scattering albedo array in LUT K
    data rgfLUTRInfSSAGrid/0.800,0.810,0.820,0.830,0.840,0.850,0.860,0.870,0.880,0.890, &
                           0.900,0.910,0.920,0.930,0.940,0.950,0.960,0.970,0.980,0.990, &
                           0.991,0.992,0.993,0.994,0.995,0.996,0.997,0.998,0.999,1.000/
    data rgfLUTKSSAGrid/0.800,0.900,0.950,0.980,0.990,0.995,0.999,1.000/ 
    data rgfLUTRdInfSSAGrid/0.500,0.600,0.700,0.800,0.900,0.950,0.990,0.999,0.9999/ 

    end module SLALOM_defs

    !######################################################################################### 
!! 
!!  NAME 
!!  Program "SLALOM" 
!! 
!!  VERSION 
!!  $Revision: 1.7 $ 
!!
!!  PURPOSE 
!!  Main program control. 
!!
!!  PROCEDURE 
!!  
!!  CONTACT 
!!  Please send any comments, suggestions, criticism, or (for our sake) bug reports to 
!!  info@lcrs.de  
!! 
!!  Copyright (c) 2003-2006 
!!  Alexander A. Kokhanovsky, Thomas Nauss
!! 
!!  This program is free software without any warranty and without even the implied  
!!  warranty of merchantability or fitness for a particular purpose. 
!!  If you use the software you must acknowledge the software and its author. 
!! 
!######################################################################################### 


    program SLALOM
 
!***************************************************************************************** 
! 
!   Declaration of variables 
! 
!***************************************************************************************** 
    
    use SLALOM_defs
    implicit none
    
    integer(4) :: liCounter !! Counter
    
!Begin change
    integer(4) :: liErrorCounter
!End change   

!***************************************************************************************** 
! 
!   Format 
! 
!***************************************************************************************** 
  
950 format(3f7.2,f13.9,f8.5,f11.6,2f11.8)  !! Read input data file content
!Begin change
!960 format('SZen   PZen   RelAzm Rna     Ra      Ana      Aa      TauIn   SSAIn   ' &
!           'Tau    SSA     PAL     Aef     LWP     dTau    dBeta')
!970 format(3f7.2,4f8.5,f8.2,f8.5,f8.2,f8.5,3f8.2,f8.4,f8.4) !! Write output data file content
!End change

960 format('SZen   PZen   RelAzm Rna     Ra      Ana      Aa      TauIn   SSAIn   ' &
           'Tau    SSA     PAL     Aef     LWP          dTau         dBeta        error')
970 format(3f7.2,4f8.5,f8.2,f8.5,f8.2,f8.5,3f8.2,f13.4,f13.4,f13.4) !! Write output data file content

!***************************************************************************************** 
! 
!   Say hello...
! 
!***************************************************************************************** 

    chVersion = '$Revision: 1.7 $'
    print*, ' '
    print*, ' '
    print*, 'Program SLALOM'
    print*, chVersion
    print*, ' '
    print*, 'Forward model CLOUD'
    print*, 'Revision: 1.38'
    print*, ' '
    print*, 'Copyright (c) 2006-2007, Thomas Nauss, Alexander A. Kokhanovsky'
    print*, ' '
    print*, ' '
    print*, 'Usage:'
    print*, 'Options:'
    print*, '-1       Wavelength of the non-absorbing channel.'
    print*, '-2       Wavelength of the absorbing channel.'
    print*, '-n       Name of the LUT for Rinf for the non-absorbing channel.'
    print*, '-1       Name of the LUT for Rinf for the absorbing channel.'
    print*, '-o       Name of the output folder.'
    print*, ' '
    print*, 'slalom options'

!***************************************************************************************** 
! 
!   Read input data
! 
!***************************************************************************************** 

!    chControlFile = 'SLALOM.cfg'
!    call CLOUD_readSettings

!    Read command line parameters
!    Non-absorbing wavelength
     call fwargument(1,'-1',chCommandLineArgument)
     read(chCommandLineArgument,'(f5.3)') fWavelength1

!    Absorbing wavelength
     call fwargument(1,'-2',chCommandLineArgument)
     read(chCommandLineArgument,'(f5.3)') fWavelength2

!    LUT Rinf for non-absorbing wavelength
     call fwargument(1,'-n',chLUT_RInf1)

!    LUT Rinf for non-absorbing wavelength
     call fwargument(1,'-a',chLUT_RInf2)

!    Path to the output folder
     call fwargument(1,'-o',chOutputFolder)

     print*, fWavelength1
     print*, fWavelength2
     print*, chLUT_RInf1
     print*, chLUT_RInf2
     print*, chOutputFolder

!   Set parameters
    bWaterCloud    = .TRUE.
!    fWavelength1   = 0.856
!    fWavelength2   = 1.630
    bPause         = .FALSE.
    bPrint         = .FALSE.
    fRho           = 1.0
    bFinite        = .TRUE.
    fSSAMinValid   = 0.8d0
    fSSAMaxValid   = 1.0d0
    fAefMinValid   = 3.0
    fAefMaxValid   = 40.0
    fPi            = acos(-1.)
    bASCII = 0
!    chLUT_RInf1 = 'lut_RInf_0856_aef_06.dat'
!    chLUT_RInf2 = 'lut_RInf_1630_aef_06.dat'

    if(bASCII) then
     liRuns = 0
!!     open(501,file='SLALOM_input.dat')
!!     read(501,*)
     open(501,file='ascii_input/input_0645_1630_aef_16.dat')
!!     read(501,*)
    do
      liRuns = liRuns+1
!!      read(501,950,iostat=iError) prgfSZen(liRuns), prgfPZen(liRuns), prgfRelAzm(liRuns), &
!!                                  prgfSSAIn(liRuns), fDummy, prgfTauIn(liRuns), &
!!                                  prgfAbs(liRuns), prgfNonAbs(liRuns)
      read(501,*,iostat=iError) prgfTauIn(liRuns), prgfNonAbs(liRuns), prgfAbs(liRuns)
!  see below!! fDummy can be deleted later, too.
! prgfNonAbsAlb(liRuns), prgfAbsAlb(liRuns)
      if(iError.ne.0) exit 
     enddo
     close(501)
     liRuns = liRuns-1
     prgfSZen = 60.0d0
     prgfPZen = 0.0d0
     prgfRelAzm = 0.0d0
     prgfSSAIn = 1.0d0
     prgfNonAbsAlb = 0.0d0
     prgfAbsAlb = 0.0d0
     prgiCloud = 2
     rgfAefInput = 6.0d0
     prgiMask = 1

    else
     liRuns = 1354*2030
     liRecl = 1354*2030

     print*, 'read Sat'

     open(501,file='pacific/200107181530_ta01m_ma11danb1_na001_1000_rp0101_001000.rst',access='direct',recl=liRecl)
     read(501,rec=1) prgfSZen
     close(501)
     prgfSZen = abs(prgfSZen)

     open(501,file='pacific/200107181530_ta01m_ma01danb1_na001_1000_rp0101_001000.rst',access='direct',recl=liRecl)
     read(501,rec=1) prgfPZen
     close(501)
     prgfPZen = abs(prgfPZen)

     open(501,file='pacific/200107181530_ta01m_ma22danb1_na001_1000_rp0101_001000.rst',access='direct',recl=liRecl)
     read(501,rec=1) prgfRelAzm
     close(501)

     open(501,file='pacific/200107181530_ta01m_pa20dlnb1_na001_1000_rp0101_001000.rst',access='direct',recl=liRecl)
     read(501,rec=1) prgiCloud
     close(501)

     if(fWavelength1.gt.0.855.and.fWavelength1.lt.0.857) then
      open(501,file='pacific/200107181530_ta01m_ca01p0002_na001_1000_rp0101_001000.rst',access='direct',recl=liRecl)
     elseif(fWavelength1.gt.0.644.and.fWavelength1.lt.0.646) then
      open(501,file='pacific/200107181530_ta01m_ca01p0001_na001_1000_rp0101_001000.rst',access='direct',recl=liRecl)
     else
      print*, 'Wrong channel'
      stop
     endif
     read(501,rec=1) prgfNonAbs
     close(501)
     prgfNonAbs = prgfNonAbs / cosd(prgfSZen)

     open(501,file='pacific/200107181530_ta01m_ca01p0006_na001_1000_rp0101_001000.rst',access='direct',recl=liRecl)
     read(501,rec=1) prgfAbs
     close(501)
     prgfAbs = prgfAbs / cosd(prgfSZen)

     open(501,file='pacific/200107181530_ta01m_py75umnb1_na001_1000_rp0101_001000.rst',access='direct',recl=liRecl)
!     open(501,file='pacific_out/lut_RInf_1630_aef_10_200107181530_ta01m_py81umnb1_na001_1000_rp0101_001000.rst',access='direct',recl=liRecl)
     read(501,rec=1) rgfAefInput
     close(501)


     open(501,file='mask/MASK_TAU_GE_10.rst',access='direct',recl=liRecl)
     read(501,rec=1) prgiMask
     close(501)

          print*, 'read Sat done'

    prgfAef = -99.0
    prgfTauRetrieve = -99.0
    prgfSSARetrieve = -99.0
    prgfPAL = -99.0
    prgfLWP = -99.0
    prgfTauIn = -99.0
    prgfSSAIn = -99.0
    prgfNonAbsAlb = 0.0
    prgfAbsAlb = 0.00
    endif


!***************************************************************************************** 
! 
!   Read LUTs
! 
!***************************************************************************************** 

!   Set imaginary part of refractive index
    call SLALOM_RefractiveIndex

!   Read LUTs necessary for CLOUD code
    call CLOUD_readLUT


    chLUT_RInf = chLUT_RInf1
!   Set asymmetry parameter used for computation of actual LUT for RInf
   if(chLUT_RInf.eq.'lut_RInf_0645_aef_06.dat') then
     fgLUT = 0.850018000000000000  ! 06  microns, 645.5nm
    elseif(chLUT_RInf.eq.'lut_RInf_0645_aef_10.dat') then
     fgLUT = 0.861757000000000000  ! 10  microns, 645.5nm
    elseif(chLUT_RInf.eq.'lut_RInf_0645_aef_16.dat') then
     fgLUT = 0.869221000000000000  ! 16 microns, 645.5nm
    elseif(chLUT_RInf.eq.'lut_RInf_0856_aef_06.dat') then
     fgLUT = 0.843722000000000000  ! 06 microns, 856.5nm
    elseif(chLUT_RInf.eq.'lut_RInf_0856_aef_10.dat') then
     fgLUT = 0.857995000000000000  ! 10 microns, 856.5nm
    elseif(chLUT_RInf.eq.'lut_RInf_0856_aef_16.dat') then
     fgLUT = 0.867233000000000000  ! 16 microns, 856.5nm
    elseif(chLUT_RInf.eq.'lut_RInf_1630_aef_06.dat') then
     fgLUT = 0.817538000000000000  ! 06 microns, 1630nm
    elseif(chLUT_RInf.eq.'lut_RInf_1630_aef_10.dat') then
     fgLUT = 0.846061000000000000  ! 10 microns, 1630nm
    elseif(chLUT_RInf.eq.'lut_RInf_1630_aef_16.dat') then
     fgLUT = 0.861364000000000000  ! 16 microns, 1630nm
    else
     print*, 'No valid LUT for RInf has been selected.'
     print*, 'The program is going to stop...'
     stop
    endif
    fgLUT1= fgLUT

    chLUT_RInf = chLUT_RInf2
    if(chLUT_RInf.eq.'lut_RInf_0645_aef_06.dat') then
     fgLUT = 0.850018000000000000  ! 06 microns, 645.5nm
    elseif(chLUT_RInf.eq.'lut_RInf_0645_aef_10.dat') then
     fgLUT = 0.861757000000000000  ! 10 microns, 645.5nm
    elseif(chLUT_RInf.eq.'lut_RInf_0645_aef_16.dat') then
     fgLUT = 0.869221000000000000  ! 16 microns, 645.5nm
    elseif(chLUT_RInf.eq.'lut_RInf_0856_aef_06.dat') then
     fgLUT = 0.843722000000000000  ! 06 microns, 856.5nm
    elseif(chLUT_RInf.eq.'lut_RInf_0856_aef_10.dat') then
     fgLUT = 0.857995000000000000  ! 10 microns, 856.5nm
    elseif(chLUT_RInf.eq.'lut_RInf_0856_aef_16.dat') then
     fgLUT = 0.867233000000000000  ! 16 microns, 856.5nm
    elseif(chLUT_RInf.eq.'lut_RInf_1630_aef_06.dat') then
     fgLUT = 0.817538000000000000  ! 06 microns, 1630nm
    elseif(chLUT_RInf.eq.'lut_RInf_1630_aef_10.dat') then
     fgLUT = 0.846061000000000000  ! 10 microns, 1630nm
    elseif(chLUT_RInf.eq.'lut_RInf_1630_aef_16.dat') then
     fgLUT = 0.861364000000000000  ! 16 microns, 1630nm
    else
     print*, 'No valid LUT for RInf has been selected.'
     print*, 'The program is going to stop...'
     stop
    endif
    fgLUT2= fgLUT


    rgfLUTRInfSSAGrid1  = 1.0 - sqrt( (1.0-rgfLUTRInfSSAGrid)/(1.0-rgfLUTRInfSSAGrid*fgLUT1) )
    rgfLUTKSSAGrid     = 1.0 - sqrt( (1.0-rgfLUTKSSAGrid)/(1.0-rgfLUTKSSAGrid*fgLUT1) )
    rgfLUTRdInfSSAGrid = 1.0 - sqrt( (1.0-rgfLUTRdInfSSAGrid)/(1.0-rgfLUTRdInfSSAGrid*fgLUT1) )

    rgfLUTRInfSSAGrid2  = 1.0 - sqrt( (1.0-rgfLUTRInfSSAGrid)/(1.0-rgfLUTRInfSSAGrid*fgLUT2) )

!***************************************************************************************** 
! 
!   Compute retrieval
! 
!***************************************************************************************** 

    open(501,file='SLALOM_output.dat')
    write(501,960)
    write(501,*) trim(chVersion), ' (CLOUD: 1.38)'
    
    do liCounter = 1, liRuns
     if(mod(liCounter,25000).eq.0) print*, 'Calculating',liCounter,' Pixel of',liRuns
     fNonAbs = prgfNonAbs(liCounter)
     fAbs    = prgfAbs(liCounter)
     fSZen   = prgfSZen(liCounter)
     fPZen   = prgfPZen(liCounter)
     fRelAzm = prgfRelAzm(liCounter)
     fNonAbsAlb = prgfNonAbsAlb(liCounter)
     fAbsAlb = prgfAbsAlb(liCounter)
     iCloud = int(prgiCloud(liCounter))
     fAefInput = rgfAefInput(liCounter)

 
! szen    pzen  relazm   ssa     g       tau       Refl   Trans  RDiff  TDiff  rs     t       dRinf  dRDInf fEsc0  fEsc
!  45.00  15.00  67.00  0.999998 0.84973  27.00 0.7081680 0.3148 0.6948 0.2528 0.7562 0.2430 1.0122 0.9390 1.0356 1.2397
! szen    pzen  relazm   ssa     g       tau       Refl   Trans  RDiff  TDiff  rs     t       dRinf  dRDInf fEsc0  fEsc
!  45.00  15.00  67.00  0.996506 0.81665  28.58 0.6185527 0.1769 0.6749 0.1431 0.6829 0.1401 0.6545 0.7039 0.8697 1.0521

!     fSZen = 45.00d0
!     fPZen = 15.00d0
!     fRelAzm = 67.00d0
!     fNonAbs = 0.7081680 ! 0.4150 !0.4154  !0.7081680d0
!     fAbs = 0.6185527 !0.4280  !0.6185527


!Begin change
!     do liErrorCounter = 800,1200
!      fNonAbs = prgfNonAbs(liCounter)*float(liErrorCounter)/1000.0
!      fAbs    = prgfAbs(liCounter)*float(liErrorCounter)/1000.0
!      if(fNonAbs.gt.1.000) fNonAbs = 1.000
!      if(fAbs.gt.1.000) fAbs = 1.000
!End change

!    Check if reflection is greater 0.0 (there might be an empty line in the input file)      
     if(fNonAbs.ne.0.000.and.iCloud.gt.1.and.prgiMask(liCounter).eq.1) then
      call SLALOM_Retrieve

      if(bASCII) then
!      Compute some statistics
       fdTau  = ( fTau/prgfTauIn(liCounter)-1.0 )*100.0
       fdBeta = ( (1-dSSARetrieve)/(1-prgfSSAIn(liCounter))-1.0 )*100.0
       write(501,970) fSZen, fPZen, fRelAzm, fNonAbs, fAbs, fNonAbsAlb, fAbsAlb, &
                      prgfTauIn(liCounter), prgfSSAIn(liCounter), &
                      fTau, dSSARetrieve, fPAL, fAef, fLWP, &
                      fdTau, fdBeta, &
!Begin change
                     100.0
!                    (float(liErrorCounter)/1000.0-1.0)*100.0
!      enddo
       if(prgfSSAIn(liCounter).eq.0.80000) then
        write(501,*)
        write(501,*)
       endif
!       write(501,*) ' '
!       write(501,*) ' '
!End change
      else
       prgfTauRetrieve(liCounter) = fTau
!       print*, fTauRetrieve
!       print*, prgfTauRetrieve(liCounter)
       prgfSSARetrieve(liCounter) = dSSARetrieve
       prgfPAL(liCounter) = fPAL
       prgfAef(liCounter) = fAef
       prgfLWP(liCounter) = fLWP
      endif
     endif
    enddo

    close(501)

    if(.not.bASCII) then
     open(501,file=trim(chOutputFolder)//'/200107181530_ta01m_py82dlnb1_na001_1000_rp0101_001000.rst',access='direct',recl=liRecl)
     write(501,rec=1) prgfTauRetrieve
     close(501)

     open(501,file=trim(chOutputFolder)//'/SSARetrieve.rst',access='direct',recl=liRecl)
     write(501,rec=1) prgfSSARetrieve
     close(501)

     open(501,file=trim(chOutputFolder)//'/PAL.rst',access='direct',recl=liRecl)
     write(501,rec=1) prgfPAL
     close(501)

     open(501,file=trim(chOutputFolder)//'/200107181530_ta01m_py81umnb1_na001_1000_rp0101_001000.rst',access='direct',recl=liRecl)
     write(501,rec=1) prgfAef
     close(501)

     open(501,file=trim(chOutputFolder)//'/LWP.rst',access='direct',recl=liRecl)
     write(501,rec=1) prgfLWP
     close(501)
    endif

    print*,'SLALOM finished.'

    end program SLALOM

!######################################################################################### 
!! 
!!  NAME 
!!  Subroutine "SLALOM_RefractiveIndex" 
!! 
!!  VERSION 
!!  $Revision: 1.7 $ 
!! 
!!  PURPOSE 
!!  Read imaginary part of refractive index from LUTs.
!!
!!  CONTACT 
!!  Please send any comments, suggestions, criticism, or (for our sake) bug reports to 
!!  info@lcrs.de  
!! 
!!  Copyright (c) 2003-2006 
!!  Alexander A. Kokhanovsky, Thomas Nauss
!! 
!!  This program is free software without any warranty and without even the implied  
!!  warranty of merchantability or fitness for a particular purpose. 
!!  If you use the software you must acknowledge the software and its author. 
!! 
!######################################################################################### 
 
    subroutine SLALOM_RefractiveIndex
 
!***************************************************************************************** 
! 
!   Declaration of variables 
! 
!***************************************************************************************** 

    use SLALOM_defs
    implicit none

    integer(2) :: iCounter      !! Counter
    integer(2) :: iCheck        !! Check flag for read statement

    real(8) :: dWavelengthCheckSmaler   !! Smaler wavelength in LUT file
    real(8) :: dRealPartSmaler          !! Smaler value of real part of refreactive index
    real(8) :: dImaginaryPartSmaler     !! Larger value of real part of refreactive index
    real(8) :: dWavelengthCheckLarger   !! Larger wavelength in LUT file
    real(8) :: dRealPartLarger          !! Larger value of real part of refreactive index
    real(8) :: dImaginaryPartLarger     !! Larger value of real part of refreactive index
    real(8) :: dRelativeWavelength      !! Relative difference between smaler/larger wavelength
    real(8) :: dTemp                    !! Temporary variable for swithing min/max

!***************************************************************************************** 
! 
!  Format
! 
!*****************************************************************************************     


!***************************************************************************************** 
! 
!  Open LUT file.
! 
!*****************************************************************************************     

    if(bWaterCloud) then
     open(501,file='lut_water_1981.dat')
     do iCounter = 1,4
      read(501,*)
     enddo
    else
     open(501,file='lut_ice_1995.dat')
     do iCounter = 1,5
      read(501,*)
     enddo
    endif

!***************************************************************************************** 
! 
!  Search and read imaginary part of refractive index with respect to wavelength
! 
!*****************************************************************************************     
    
!   Look for the two bounding value of the refractive index
    read(501,*) dWavelengthCheckSmaler, dRealPartSmaler, dImaginaryPartSmaler
    do
     read(501,*,iostat=iCheck) dWavelengthCheckLarger, dRealPartLarger, dImaginaryPartLarger
     if(iCheck.lt.0) then
      print*, 'No appropriate refractive index has been found.'
      print*, 'Program is going to stop...'
      stop
     endif
     if(dWavelengthCheckSmaler.le.fWavelength2.and.dWavelengthCheckLarger.ge.fWavelength2) then
      exit
     else
      dWavelengthCheckSmaler = dWavelengthCheckLarger
      dRealPartSmaler        = dRealPartLarger
      dImaginaryPartSmaler   = dImaginaryPartLarger
     endif
    enddo     

!   Compute refractive index
    dRelativeWavelength = (fWavelength2-dWavelengthCheckSmaler)/(dWavelengthCheckLarger-dWavelengthCheckSmaler)
    call CLOUD_1DInterpolation(dRelativeWavelength,dImaginaryPartSmaler,dImaginaryPartLarger,dTemp)
    fMiIce = dTemp
    return
    end subroutine SLALOM_RefractiveIndex
    
!######################################################################################### 
!! 
!!  NAME 
!!  Subroutine "SLALOM_Retrieve" 
!! 
!!  VERSION 
!!  $Revision: 1.7 $ 
!! 
!!  PURPOSE 
!!  Compute cloud properties.
!!
!!  CONTACT 
!!  Please send any comments, suggestions, criticism, or (for our sake) bug reports to 
!!  info@lcrs.de  
!! 
!!  Copyright (c) 2003-2006 
!!  Alexander A. Kokhanovsky, Thomas Nauss
!! 
!!  This program is free software without any warranty and without even the implied  
!!  warranty of merchantability or fitness for a particular purpose. 
!!  If you use the software you must acknowledge the software and its author. 
!! 
!######################################################################################### 
 
    subroutine SLALOM_Retrieve
     
!***************************************************************************************** 
! 
!   Declaration of variables 
! 
!***************************************************************************************** 

    use SLALOM_defs
    implicit none

!    external   functionR        !! Rabs-Rcomp
!    external   functionTest        !! Rabs-Rcomp
    external functionAef

    integer(2) :: iCounter

    real(8) :: dSumC
    real(8) :: dSumD 
    real(8),dimension(0:4) :: drgC1  !! Parameters for wavelength 0.6457µm
    real(8),dimension(0:4) :: drgD

    data drgC1/0.1121,0.5118,0.8997,0.0,0.0/ ! 0.6457
!    data drgC1/0.1115,0.4513,1.2719,0.0,0.0/  ! 0.8590
!    data drgC2/0.0608,2.465,-32.98,248.94,-636.0/
!    data drgD/1.671,0.0025,-2.365e-4,2.861e-6,-1.05e-8/
    real(8) :: zbrent           !! Find root of functionR using Brent's method

!***************************************************************************************** 
! 
!  Retrieve effective cloud droplet radius.
! 
!*****************************************************************************************     

!   Compute reflection of a semi-infinite cloud
    bFinite = .TRUE.
    fTau    = 5.0        ! Value does not matter.
    fSSA    = 1.0
!!    fg = 0.850018000000000000

!    rgfLUTRInfSSAGrid  = 1.0 - sqrt( (1.0-rgfLUTRInfSSAGridIn)/(1.0-rgfLUTRInfSSAGridIn*fg) )
!    rgfLUTKSSAGrid     = 1.0 - sqrt( (1.0-rgfLUTKSSAGridIn)/(1.0-rgfLUTKSSAGridIn*fg) )
!    rgfLUTRdInfSSAGrid = 1.0 - sqrt( (1.0-rgfLUTRdInfSSAGridIn)/(1.0-rgfLUTRdInfSSAGridIn*fg) )
!    rgfLUTRInfSSAGrid = rgfLUTRInfSSAGrid1
!    rgfLUTKSSAGrid = rgfLUTKSSAGrid1
!    rgfLUTRdInfSSAGrid = rgfLUTRdInfSSAGrid1
!    rgiLUTRInf = rgiLUTRInf1

    LUT = 1
    call cloud

    dR0 = dRinf_1000
    bFinite = .TRUE.
!!    call cloud
!   Compute transmission for non-absorbing measurement
    ft1    = ( (1.0-fNonAbsAlb)*(dRinf_1000-fNonAbs) ) / &
             ( (1.0-fNonAbsAlb)*fEscape0_1000*fEscape_1000 - fNonAbsAlb*(dRInf-fNonAbs) )

!!     fEscape0_1000 = 3.0/7.0*(1.0+2.0*cosd(fSZen))
!!     fEscape_1000  = 3.0/7.0*(1.0+2.0*cosd(fPZen))
!!     ft1    = ( (1.0-fNonAbsAlb)*(dRinf_1000-fNonAbs) ) / &
!!             ( (1.0-fNonAbsAlb)*fEscape0_1000*fEscape_1000 - fNonAbsAlb*(dRInf-fNonAbs) )
!!     print*, fEscape0_1000, fEscape_1000
!!     print*, ft1



!   Retrieve effective radius
    fAef = zbrent(functionAef,fAefMinValid,fAefMaxValid,1d-9)
!!   fAef = fAefInput
!!    fAef = 6.0
!   Compute tau using aef
!!    ft     = ( (1-fNonAbsAlb)*(dRinf_1000-fNonAbs) ) / &
!!             ( (1-fNonAbsAlb)*fEscape0_1000*fEscape_1000 - fNonAbsAlb*(dRInf-fNonAbs) )
!!    fAef = 16.0000

    dSumC = 0.d0
    do iCounter = 0,4
     dSumC = dSumC + drgC1(iCounter) * ( 2.0*fPi/fWavelength1 * fAef )**(-2.0*float(iCounter)/3.0)
    enddo
    fg = 1.d0 - dSumC

    fTau   = (ft1**(-1.0)-1.072) / (0.75*(1.0-fg))
!   Compute LWP from aef for non-absorbing wavelength
     fLWP = fTau * fRho * fAef / &
           1.5 * ( 1.0 + 1.1 / ( 2.0*fPi/fWavelength1 * fAef ) ** (2.0/3.0) ) ** (-1.0)

!    print*, 'R_mes_0650:                   ', fNonAbs
!    print*, 'R_mes_1640:                   ', fAbs
!    print*, 'aef retrieved:                ', fAef
!    print*, 'tau retrieved:                ', fTau
!   print*, 'LWP retrieved:                ', fLWP
!    print*, 't_non_abs used within brent:  ', ft1
!    print*, 'last value of ssa in brent:   ', fSSA
!    print*, 'fg used to compute final tau: ', fg
!    pause

    return
    end subroutine SLALOM_Retrieve

!######################################################################################### 
!! 
!!  NAME 
!!  Function "functionAef" 
!! 
!!  VERSION 
!!  $Revision: 1.7 $ 
!! 
!!  PURPOSE 
!!  Compute Rmes_absorbing - Rcomp_absorbing.
!!
!!  CONTACT 
!!  Please send any comments, suggestions, criticism, or (for our sake) bug reports to 
!!  info@lcrs.de  
!! 
!!  Copyright (c) 2003-2006 
!!  Alexander A. Kokhanovsky, Thomas Nauss
!! 
!!  This program is free software without any warranty and without even the implied  
!!  warranty of merchantability or fitness for a particular purpose. 
!!  If you use the software you must acknowledge the software and its author. 
!! 
!######################################################################################### 
 
    function functionAef(dAef)
 
!***************************************************************************************** 
! 
!   Declaration of variables 
! 
!***************************************************************************************** 

    use SLALOM_defs
    implicit none

    integer(2) :: iCounter

    real(8) :: dAef         !! Actual value of aef within Brent's routine (functionBrent)
    real(8) :: functionAef  !! Rabs-Rcomp
    real(8) :: dSumC
    real(8) :: dSumD 
    real(8),dimension(0:4) :: drgC1,drgC2  !! Parameters for wavelength 0.6457µm
    real(8),dimension(0:4) :: drgD

    data drgC1/0.1121,0.5118,0.8997,0.0,0.0/ !0.6457
!    data drgC1/0.1115,0.4513,1.2719,0.0,0.0/  !0.8590
    data drgC2/0.0608,2.465,-32.98,248.94,-636.0/
    data drgD/1.671,0.0025,-2.365e-4,2.861e-6,-1.05e-8/

!***************************************************************************************** 
! 
!  Compute 
! 
!*****************************************************************************************     

! szen    pzen  relazm   ssa     g       tau       Refl   Trans  RDiff  TDiff  rs     t       dRinf  dRDInf fEsc0  fEsc
!  45.00  15.00  67.00  0.999998 0.84973  27.00 0.7081680 0.3148 0.6948 0.2528 0.7562 0.2430 1.0122 0.9390 1.0356 1.2397
! szen    pzen  relazm   ssa     g       tau       Refl   Trans  RDiff  TDiff  rs     t       dRinf  dRDInf fEsc0  fEsc
!  45.00  15.00  67.00  0.996506 0.81665  28.58 0.6185527 0.1769 0.6749 0.1431 0.6829 0.1401 0.6545 0.7039 0.8697 1.0521



!   Compute asymetry parameter for non-absorbing wavelength
    dSumC = 0.d0
    do iCounter = 0,4
     dSumC = dSumC + drgC1(iCounter) * ( 2.0*fPi/fWavelength1 * dAef )**(-2.0*float(iCounter)/3.0)
    enddo
    fg = 1.d0 - dSumC

!   Compute tau for non-absorbing wavelength
    fTau   = (ft1**(-1.0)-1.072) / (0.75*(1.0-fg))

!   Compute LWP for non-absorbing wavelength
    fLWP = fTau * fRho * dAef / &
           1.5 * ( 1.0 + 1.1 / ( 2.0*fPi/fWavelength1 * dAef ) ** (2.0/3.0) ) ** (-1.0)

!   Compute tau for absorbing wavelength using LWP from non-absorbing wavelength
    fTau = 1.5 * fLWP / ( fRho * dAef) * ( 1.0 + 1.1 / ( 2.0*fPi/fWavelength2 * dAef ) ** (2.0/3.0) )

!   Compute the extinction coefficient for absorbing wavelength
    fSigma_ext = 1.5 / dAef * (1.0 + 1.1 / ( 2.0*fPi/fWavelength2 * dAef ) ** (2.0/3.0) )

!   Compute the absorbtion coefficient for absorbing wavelength
    dSumD = 0.d0
    do iCounter = 0,4
     dSumD = dSumD + drgD(iCounter) * ( 2.0*fPi/fWavelength2 * dAef )**float(iCounter)
    enddo
    fSigma_abs = 4.0 * fPi * fMiIce / fWavelength2 * dSumD

!   Compute ssa for absorbing wavelength based on actual aef value.
    fSSA = 1.d0 - fSigma_abs / fSigma_ext

!   Compute asymetry parameter for absorbing wavelength based on actual aef value.
    dSumC = 0.d0
    do iCounter = 0,4
     dSumC = dSumC + drgC2(iCounter) * ( 2.0*fPi/fWavelength2 * dAef )**(-2.0*float(iCounter)/3.0)
    enddo
    fg = 1.d0 - dSumC

!   Compute reflectance for absorbing channel
!    rgfLUTRInfSSAGrid  = 1.0 - sqrt( (1.0-rgfLUTRInfSSAGridIn)/(1.0-rgfLUTRInfSSAGridIn*fg) )
!    rgfLUTKSSAGrid     = 1.0 - sqrt( (1.0-rgfLUTKSSAGridIn)/(1.0-rgfLUTKSSAGridIn*fg) )
!    rgfLUTRdInfSSAGrid = 1.0 - sqrt( (1.0-rgfLUTRdInfSSAGridIn)/(1.0-rgfLUTRdInfSSAGridIn*fg) )
!    rgfLUTRInfSSAGrid = rgfLUTRInfSSAGrid2
!    rgfLUTKSSAGrid = rgfLUTKSSAGrid2
!    rgfLUTRdInfSSAGrid = rgfLUTRdInfSSAGrid2
!    rgiLUTRInf = rgiLUTRInf2
    LUT = 2
    call cloud

!   Compute difference
    functionAef = fAbs - fRefl


!***************************************************************************************** 
! 
!  Compute R_mes(1640nm)-R_computed
! 
!*****************************************************************************************     

!!    fSSA = dSSATest
!!    fSSALook = fSSA
!!    call cloud
!Begin change?
!!    functionR = fAbs - fRefl
!Only one option - but which?
!    functionR = fAbs - &
!                dRinf + fl * fTrans * fn**(-2) * fEscape0 * fEscape * exp(-fk*fTau_1000)- &
!                ( fAbsAlb * fTrans**2 * fn**(-2) * fEscape0 * fEscape ) / ( 1 - fAbsAlb*frs)
!End change
    
    end function functionAef

!######################################################################################### 
!! 
!!  NAME 
!!  Function "functionTest" 
!! 
!!  VERSION 
!!  $Revision: 1.7 $ 
!! 
!!  PURPOSE 
!!  Compute Rmes-Rcomp.
!!
!!  CONTACT 
!!  Please send any comments, suggestions, criticism, or (for our sake) bug reports to 
!!  info@lcrs.de  
!! 
!!  Copyright (c) 2003-2006 
!!  Alexander A. Kokhanovsky, Thomas Nauss
!! 
!!  This program is free software without any warranty and without even the implied  
!!  warranty of merchantability or fitness for a particular purpose. 
!!  If you use the software you must acknowledge the software and its author. 
!! 
!######################################################################################### 
 
    function functionTest(fAefTest)
 
!***************************************************************************************** 
! 
!   Declaration of variables 
! 
!***************************************************************************************** 

    use SLALOM_defs
    implicit none
    
    integer(2) :: iCounter
    
    real(8) :: dAefTest     !! Actual value of SSA within Brent's routine (functionBrent)
    real(8) :: functionTest !! Rabs-Rcomp
    real(8) :: fAefTest
    real(8),dimension(0:4) :: frgD
    real(8) :: fSumD 

!***************************************************************************************** 
! 
!  Compute R_mes(1640nm)-R_computed
! 
!*****************************************************************************************     

    data frgD/1.671,0.0025,-2.365e-4,2.861e-6,-1.05e-8/
    
    fAef = fAefTest
    fSumD = 0
    do iCounter = 0,4
     fSumD = fSumD + frgD(iCounter)*(2*fPi/fWavelength2 * fAef)**iCounter
    enddo
    
    functionTest = (1-dSSARetrieve) - &
                   ( (4*fPi*fMiIce/fWavelength2)*1.0*fSumD ) / &
                   ( 1.5/fAef * (1+1.1/(2*fPi/fWavelength2*fAef)**(2/3)) )
    
    
    end function functionTest

!######################################################################################### 
!! 
!!  NAME 
!!  Function "zbrent" 
!! 
!!  VERSION 
!!  $Revision: 1.7 $ 
!! 
!!  PURPOSE 
!!  Find root of function using Brent's method (functionR).
!!
!!  CONTACT 
!!  Please send any comments, suggestions, criticism, or (for our sake) bug reports to 
!!  info@lcrs.de  
!! 
!!  Copyright (c) 1986-1992
!!  Numerical Recipes Software ,#23
!! 
!######################################################################################### 
 
    function zbrent(func,x1,x2,tol)
 
!***************************************************************************************** 
! 
!   Declaration of variables 
! 
!***************************************************************************************** 

!    implicit none

!    INTEGER(4) ITMAX
!    real*8 :: x1, x2
!    REAL(8) tol,func,EPS
!    real(8) functionBrent   !! Root of functionR
!    real(8) :: functionR    !! Rabs-Rcomp
!    PARAMETER (ITMAX=100,EPS=3.e-8)
!    INTEGER(4) iter
!    REAL(8) a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
!    a=x1
!    b=x2
!    fa=functionR(a)
!    fb=functionR(b)
      
    INTEGER(4) ::  ITMAX
    REAL(8) ::  zbrent,tol,func,EPS !! Nauss: x1,x2,
    EXTERNAL func
    PARAMETER (ITMAX=100,EPS=3.e-8)
    INTEGER(4) ::  iter
    REAL(8) ::  a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
    real*8 :: x1,x2    !! Nauss
    a=x1
    b=x2
    fa=func(a)
    fb=func(b)
    if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.))then
!     print*, 'root must be bracketed for zbrent'
!     pause
    endif
    c=b
    fc=fb
    do 11 iter=1,ITMAX
     if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
      c=a
      fc=fa
      d=b-a
      e=d
     endif
     if(abs(fc).lt.abs(fb)) then
      a=b
      b=c
      c=a
      fa=fb
      fb=fc
      fc=fa
     endif
     tol1=2.*EPS*abs(b)+0.5*tol
     xm=.5*(c-b)
     if(abs(xm).le.tol1 .or. fb.eq.0.)then
      zbrent=b
      return
     endif
     if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
      s=fb/fa
      if(a.eq.c) then
       p=2.*xm*s
       q=1.-s
      else
       q=fa/fc
       r=fb/fc
       p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
       q=(q-1.)*(r-1.)*(s-1.)
      endif
      if(p.gt.0.) q=-q
      p=abs(p)
      if(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
       e=d
       d=p/q
      else
       d=xm
       e=d
      endif
     else
      d=xm
      e=d
     endif
     a=b
     fa=fb
     if(abs(d) .gt. tol1) then
      b=b+d
     else
      b=b+sign(tol1,xm)
     endif
     fb=func(b)
11  continue
    pause 'zbrent exceeding maximum iterations'
    zbrent=b
    return
!   (C) Copr. 1986-92 Numerical Recipes Software ,#23.
    end function zbrent

!***************************************************************************************** 
!***************************************************************************************** 
! 
!   Include program cloud here (and adjust main routine so it can be used as subroutine)
!
!***************************************************************************************** 
!***************************************************************************************** 

!######################################################################################### 
!! 
!!  NAME 
!!  Program "CLOUD" 
!! 
!! 
!!  VERSION 
!!  $Revision: 1.7 $ 
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
!!  Please send any comments, suggestions, criticism, or (for our sake) bug reports to 
!!  info@lcrs.de  
!! 
!! 
!!  Copyright (c) 2003-2006 
!!  Alexander A. Kokhanovsky, Thomas Nauss
!! 
!!  This program is free software without any warranty and without even the implied  
!!  warranty of merchantability or fitness for a particular purpose. 
!!  If you use the software you must acknowledge the software and its author. 
!! 
!######################################################################################### 

    subroutine CLOUD

!############# UNCOMMENT THE NEXT LINE IF USED WITHIN SLALOM #############################
!    subroutine CLOUD

!***************************************************************************************** 
!
!   Declaration of variables 
!
!***************************************************************************************** 

!    use CLOUD_Defs

!############# UNCOMMENT THE NEXT LINE IF USED WITHIN SLALOM #############################
    use SLALOM_Defs

    implicit none


!***************************************************************************************** 
!
!   Format
!
!***************************************************************************************** 

!#!950 format(a8,f15.7)    !! Print final results
!#!960 format(' szen    pzen  relazm  ssa     g       ' &
!#!           'tau    Refl   Trans  RDiff  TDiff  rs     t       ' &
!#!           'dRinf  dRDInf fEsc0  fEsc')
!#!961 format(3f7.2,f11.7,f8.5,f7.2,6f7.4,4f7.4)
!#!970 format(' szen    pzen  relazm  ssa     g       ' &
!#!           'tau    Refl   RDiff  rs     Trans  TDiff  t')
!#!971 format(3f7.2,f11.7,f8.5,f7.2,6f7.4)

!***************************************************************************************** 
!
!   Say hello...
!
!***************************************************************************************** 

!#!    chVersion = '$Revision: 1.7 $'
!#!    print*, ' '
!#!    print*, ' '
!#!    print*, 'Program CLOUD'
!#!    print*, trim(chVersion)

!***************************************************************************************** 
!
!   Definition of input parameters
!
!***************************************************************************************** 

!############# COMMENT THE FOLLOWING SECTION IF USED WITHIN SLALOM #######################
!#!    chControlFile    = 'cloud.cfg'
!#!    call CLOUD_ReadSettings

!   Set asymmetry parameter used for computation of actual LUT for RInf
!#!    if(chLUT_RInf.eq.'lut_RInf_0645_aef_06.dat') then
!#!     fgLUT = 0.850018000000000000  ! 06�m, 645.5nm
!#!     rgchOutFile(1) = 'solar_lut_0645_06_.dat'
!#!    elseif(chLUT_RInf.eq.'lut_RInf_0645_aef_10.dat') then
!#!     fgLUT = 0.861757000000000000  ! 10�m, 645.5nm
!#!     rgchOutFile(1) = 'solar_lut_0645_10_.dat'
!#!    elseif(chLUT_RInf.eq.'lut_RInf_0645_aef_16.dat') then
!#!     fgLUT = 0.869221000000000000  ! 16�m, 645.5nm
!#!     rgchOutFile(1) = 'solar_lut_0645_16_.dat'
!#!    elseif(chLUT_RInf.eq.'lut_RInf_1630_aef_06.dat') then
!#!     fgLUT = 0.817538000000000000  ! 06�m, 1630nm
!#!     rgchOutFile(1) = 'solar_lut_1630_06_.dat'
!#!    elseif(chLUT_RInf.eq.'lut_RInf_1630_aef_10.dat') then
!#!     fgLUT = 0.846061000000000000  ! 10�m, 1630nm
!#!     rgchOutFile(1) = 'solar_lut_1630_10_.dat'
!#!    elseif(chLUT_RInf.eq.'lut_RInf_1630_aef_16.dat') then
!#!     fgLUT = 0.861364000000000000  ! 16�m, 1630nm
!#!     rgchOutFile(1) = 'solar_lut_1630_16_.dat'
!#!    else
!#!     print*, 'No valid LUT for RInf has been selected.'
!#!     print*, 'The program is going to stop...'
!#!     stop
!#!    endif

!   Open output filenames    
!#!    open(550,file=trim(rgchOutFile(1)))
!#!    write(550,*) 'Version: ',trim(chVersion)
!#!    write(550,970)

!***************************************************************************************** 
!
!   Read LUTs
!
!***************************************************************************************** 

!#!    call CLOUD_ReadLUT

!***************************************************************************************** 
!
!   Compute radiative characteristics
!
!***************************************************************************************** 

!   Program controll settings
!#!    liTauMin = nint(fTauMin*100)
!#!    liTauMax = nint(fTauMax*100)
!#!    liTauStep = nint(fTauStep*100)
!#!    liSSAMin = nint(fSSAMin * 1000000)
!#!    liSSAMax = nint(fSSAMax * 1000000)
!#!    liSSAStep = nint(fSSAStep * 1000000)
!#!    liSZenMin = nint(fSZenMin*100)
!#!    liSZenMax = nint(fSZenMax*100)
!#!    liSZenStep = nint(fSZenStep*100)
!#!    liPZenMin = nint(fPZenMin*100)
!#!    liPZenMax = nint(fPZenMax*100)
!#!    liPZenStep = nint(fPZenStep*100)
!#!    liRelAzmMin = nint(fRelAzmMin*100)
!#!    liRelAzmMax = nint(fRelAzmMax*100)
!#!    liRelAzmStep = nint(fRelAzmStep*100)

!   Convert ssa grid values of the LUTs to similarity parameter (1-fs is used for historical
!   reasons since then, the coloumn with the smalest/largest ssa value is the same).
!   LUTs have been computed for an asymmetry parameter of 0.8500, this value is used for the
!   conversion.

!#!    rgfLUTRInfSSAGrid  = 1.0 - sqrt( (1.0-rgfLUTRInfSSAGrid)/(1.0-rgfLUTRInfSSAGrid*fgLUT) )
!#!    rgfLUTKSSAGrid     = 1.0 - sqrt( (1.0-rgfLUTKSSAGrid)/(1.0-rgfLUTKSSAGrid*fgLUT) )
!#!    rgfLUTRdInfSSAGrid = 1.0 - sqrt( (1.0-rgfLUTRdInfSSAGrid)/(1.0-rgfLUTRdInfSSAGrid*fgLUT) )

!#!    do 10 liTauCounter = liTauMin, liTauMax, liTauStep
!#!     do 20 liSZenCounter = liSZenMin, liSZenMax, liSZenStep
!#!      do 30 liPZenCounter = liPZenMin, liPZenMax, liPZenStep
!#!       do 40 liRelAzmCounter = liRelAzmMin, liRelAzmMax, liRelAzmStep
!#!        do 50 liSSACounter = liSSAMin, liSSAMax, liSSAStep

!#!        fTau = float(liTauCounter)/100.0
!#!        fSSA = float(liSSACounter)/1000000.0
!#!        fSZen = float(liSZenCounter)/100.0
!#!        fPZen = float(liPZenCounter)/100.0 
!#!        fRelAzm = float(liRelAzmCounter)/100.0

!############# END OF SECTION TO COMMENT IF USED WITHIN SLALOM ###########################

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
!       Compute R_infinite additionaly for ssa = 1.000 if actual SSA is greater/equal
!       0.999 (see computation of reflection function later).
!       Use 1 - similarity parameter for the interpolation instead of ssa in order to
!       correct for actual asymmetry parameter.
        if(fSSA.ge.0.999999) then
         fSSALook = 1.000
         fSSALook = 1.0 - sqrt( (1.0-fSSALook)/(1.0-fSSALook*fg) )
         if(LUT.eq.1) then
          call CLOUD_RInfinite1
         elseif(LUT.eq.2) then
          call CLOUD_RInfinite2
         else
          stop
         endif
        dRinf_1000 = dRinf
        else
         fSSALook = fSSA
         fSSALook = 1.0 - sqrt( (1.0-fSSALook)/(1.0-fSSALook*fg) )
         if(LUT.eq.1) then
          call CLOUD_RInfinite1
         elseif(LUT.eq.2) then
          call CLOUD_RInfinite2
         else
          stop
         endif
        endif

        fz=fMu0-0.5
        fa1=-0.9991+3.139*fg-1.874*fg*fg
        fa2=1.435-4.924*fg+2.089*fg*fg
        fa3=0.719-5.801*fg+2.117*fg*fg
        fa4=-0.509+0.418*fg+3.36*fg*fg
        fP=(fa1*fz+fa2*fz*fz)*fs+(fa3*fz+fa4*fz*fz)*fs*fs

!       Compute actual plane albedo for a semi-infinite layer with respect to mu0
!       dRDInf=(1.-fs)*exp(fP)/(1.+2.*fs*fMu0)
!       Compute actual plane albedo additionaly for ssa = 1.000 if actual SSA is
!       greater/equal 0.999 (see computation of reflection function later).
!       Use 1 - similarity parameter for the interpolation instead of ssa in order to
!       correct for actual asymmetry parameter.
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
!       This if is necessary since otherwise some environments produce a runtime error
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
!       Use 1 - similarity parameter for the interpolation instead of ssa in order to
!       correct for actual asymmetry parameter.
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

!*****************************************************************************************
!
!   Print/write output
!
!*****************************************************************************************

!############# COMMENT THE FOLLOWING SECTION IF USED WITHIN SLALOM #######################

!#!     if(.not.bFinite  ) then
!#!      write(550,971) fSZen, fPZen, fRelAzm, fSSA, fg, &
!#!                     fTau, fRefl, fRDiff, frs, fTrans, fTdiff, ft
!#!     else
!#!      write(550,961) fSZen, fPZen, fRelAzm, fSSA, fg, &
!#!                     fTau, fRefl, fTrans, fRDiff, fTdiff, frs,ft, &
!#!                     dRinf, dRDInf, fEscape0, fEscape
!#!     endif


!#!     if(bPrint) then
!#!         print*,' '
!#!         print*,'Final results from program CLOUD:'
!#!         print*,'tau:           ', fTau
!#!         print*,'rs:            ', frs
!#!         print*,'t:             ', ft
!#!         print*,'RDiff:         ', fRDiff
!#!         print*,'TDiff:         ', fTDiff
!#!         print*,'ADiff:         ', fADiff
!#!         print*,'Trans:         ', fTrans
!#!         print*,'Refl:          ', fRefl
!#!         print*,'Ssa:           ', fSSA
!#!         print*,'fg:            ', fg
!#!         print*,'fSZen:         ', fSZen
!#!         print*,'fPZen:         ', fPZen
!#!         print*,'fRelAzm:       ', fRelAzm
!#!        pause
!#!         if(bPause) pause
!#!        endif

!#!50      continue
!#!40     continue
!#!30    continue
!#!20   continue
!#!10  continue

!#!    close(550)
!#!    close(551)

!#!    stop
!#!    end program CLOUD

!############# END OF SECTION TO COMMENT IF USED WITHIN SLALOM ###########################
!############# UNCOMMENT THE NEXT TWO LINES IF USED WITHIN SLALOM ########################
    return
    end subroutine CLOUD

!######################################################################################### 
!! 
!!  NAME 
!!  Subroutine "CLOUD_ReadSettings" 
!! 
!! 
!!  VERSION 
!!  $Revision: 1.7 $ 
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
!!  Please send any comments, suggestions, criticism, or (for our sake) bug reports to 
!!  info@lcrs.de  
!! 
!! 
!!  Copyright (c) 2003-2006 
!!  Alexander A. Kokhanovsky, Thomas Nauss
!! 
!!  This program is free software without any warranty and without even the implied  
!!  warranty of merchantability or fitness for a particular purpose. 
!!  If you use the software you must acknowledge the software and its author. 
!! 
!######################################################################################### 
 
    subroutine CLOUD_ReadSettings
 
!***************************************************************************************** 
! 
!   Declaration of variables 
! 
!***************************************************************************************** 

    use SLALOM_Defs
    implicit none

!***************************************************************************************** 
! 
!   Define namelist
! 
!***************************************************************************************** 

    namelist /CLOUDControl/ &
    bFinite  ,     &
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

!*****************************************************************************************
!
!   Read configuration file
!
!*****************************************************************************************

    open(501,file=chControlFile)
    rewind(501)
    read(501,nml=CLOUDControl)
    close(501)

    return
    end subroutine CLOUD_ReadSettings 


!######################################################################################### 
!! 
!!  NAME 
!!  Subroutine "CLOUD_ReadLUT" 
!! 
!! 
!!  VERSION 
!!  $Revision: 1.7 $ 
!! 
!! 
!!  PURPOSE 
!!  Read tabulated values from ASCII file.
!!  1. Read lut_k.dat for 2D interpolation of K(mu[0],ssa)
!!     K/2 equals escape function
!!  2. Read LUT interp_ssa_wwww.dat for 3D interpolation of R(szen,pzen,relazm,ssa)
!!     "_ssa"  = single-scattering albedo (ssa=0.95 --> "_0_95")
!!     "_wwww" = wavelenght in  microns (wwww=0.865 microns --> "_0865")
!!  3. Read LUT lut_RdInf.dat for 1D interpolation of RDinf(szen,ssa)
!! 
!! 
!!  PROCEDURE 
!!  Read tabulated values from ASCII file. 
!! 
!!
!!  CONTACT 
!!  Please send any comments, suggestions, criticism, or (for our sake) bug reports to 
!!  info@lcrs.de  
!! 
!! 
!!  Copyright (c) 2003-2006 
!!  Alexander A. Kokhanovsky, Thomas Nauss
!! 
!!  This program is free software without any warranty and without even the implied  
!!  warranty of merchantability or fitness for a particular purpose. 
!!  If you use the software you must acknowledge the software and its author. 
!! 
!######################################################################################### 
 
    subroutine  CLOUD_ReadLUT
 
!***************************************************************************************** 
! 
!   Declaration of variables 
! 
!***************************************************************************************** 
 
    use SLALOM_Defs
    implicit none
    
    integer(2) iCounter       !! Counter for lut_K.dat
    integer(2) iCounterSZen   !! Counter for sun zenith in lut_RInf.dat
    integer(2) iCounterPZen   !! Counter for pixel zenith in lut_RInf.dat
    integer(2) iCounterRelAzm !! Counter for relative azimuth in lut_RInf.dat
    integer(2) ios

!***************************************************************************************** 
! 
!   Format 
! 
!***************************************************************************************** 
  
950 format(f17.8,8f17.8)  !! Read lut_K.dat content
951 format(9f10.6)        !! Print lut_K.dat content
952 format(33i5)          !! Read lut_RInf content
953 format(10f22.14)      !! Read lut_RdInf.dat content

!***************************************************************************************** 
! 
!   Read lut_K.dat for 2D interpolation of K(mu[0],ssa)
! 
!***************************************************************************************** 

!   Read LUT ktable 
    iCounter = 0 
    open(501,file='lut_k.dat') 
    read(501,*)
    do iCounter = 1, 30
     read(501,950) rgfLUTKMu0Grid(iCounter), rgfLUTKSSA(iCounter,1:8) 
    enddo 
    close(501)
    rgfLUTKMuGrid = rgfLUTKMu0Grid

!***************************************************************************************** 
! 
!   Read lut_RInf.dat for 4D interpolation of R(szen,pzen,fi,ssa)
! 
!***************************************************************************************** 

!   Read LUT values in array
    open(501,file=trim(chLUT_RInf1))
    read(501,*)
     do iCounterSZen = 0, 89
      do iCounterRelAzm = 0, 180
       do iCounterPZen = 0, 89
        read(501,952) idummy1, idummy2, idummy3, &
                      rgiLUTRInf1(iCounterSZen+1,iCounterPZen+1,iCounterRelAzm+1,1:30)
       enddo
      enddo
     enddo
    close(501)
    print*, 'LUT1 read'
!   Read LUT values in array
    open(501,file=trim(chLUT_RInf2))
    read(501,*)
     do iCounterSZen = 0, 89
      do iCounterRelAzm = 0, 180
       do iCounterPZen = 0, 89
        read(501,952) idummy1, idummy2, idummy3, &
                      rgiLUTRInf2(iCounterSZen+1,iCounterPZen+1,iCounterRelAzm+1,1:30)
       enddo
      enddo
     enddo
    print*, 'LUT2 read'

    close(501)

!***************************************************************************************** 
! 
!   Read lut_RdInf.dat for 1D interpolation of RDinf(szen,ssa)
! 
!***************************************************************************************** 

!   Read LUT ktable 
    iCounter = 0 
    open(501,file='lut_RdInf.dat') 
    read(501,*)
    do iCounter = 1,100
     read(501,953) rgfLUTRdInfMu0Grid(iCounter), rgfLUTRdInfSSA(iCounter,1:9) 
    enddo 
    close(501)
    
    return
    end subroutine CLOUD_ReadLUT


!######################################################################################### 
!! 
!!  NAME 
!!  Subroutine "CLOUD_RInfinite" 
!! 
!! 
!!  VERSION 
!!  $Revision: 1.7 $ 
!! 
!! 
!!  PURPOSE 
!!  Compute 4D interpolation of tabulated values of R_infinite(szen,pzen,relazm,ssa)
!!  using lut_RInf.dat.
!! 
!! 
!!  PROCEDURE 
!!  1. Compute 3D interpolation for R_infinite(szen,pzen,relazm) for LUT-grid ssa value
!!     smaler than actual ssa.
!!  2. Compute 3D interpolation for R_infinite(szen,pzen,relazm) for LUT-grid ssa value
!!     larger than actual ssa.
!!  3. Compute 1D interpolation for the two R_infinite values with respect to acutal ssa
!! 
!!
!!  CONTACT 
!!  Please send any comments, suggestions, criticism, or (for our sake) bug reports to 
!!  info@lcrs.de  
!! 
!! 
!!  Copyright (c) 2003-2006 
!!  Alexander A. Kokhanovsky, Thomas Nauss
!! 
!!  This program is free software without any warranty and without even the implied  
!!  warranty of merchantability or fitness for a particular purpose. 
!!  If you use the software you must acknowledge the software and its author. 
!!
!######################################################################################### 
 
    subroutine  CLOUD_RInfinite1
 
!***************************************************************************************** 
! 
!   Declaration of variables 
! 
!***************************************************************************************** 
 
    use SLALOM_Defs
    implicit none
 
    integer(2) :: iSmalerPosSZen
    integer(2) :: iLargerPosSZen
    integer(2) :: iSmalerPosPZen
    integer(2) :: iLargerPosPZen
    integer(2) :: iSmalerPosRelAzm
    integer(2) :: iLargerPosRelAzm
    
    real(8) :: dR
    real(8) :: dS
    real(8) :: dT
    real(8) :: dX000
    real(8) :: dX001
    real(8) :: dX010
    real(8) :: dX011
    real(8) :: dX100
    real(8) :: dX101
    real(8) :: dX110
    real(8) :: dX111
    real(8) :: dSmalerRInf !! RInf for smaler LUT-grid ssa value
    real(8) :: dLargerRInf  !! RInf for larger LUT-grid ssa value
    real(8) :: dSqrtBeta
    real(8) :: dSmalerSqrtBeta
    real(8) :: dLargerSqrtBeta
    real(8) :: dRelativeSqrtBeta
    
    real(4),dimension(1) :: rgfMinDiffSSA(30)   !! Temp. variable necessary for 4D interp.
    real(4),dimension(1) :: rgfMinLoc(1)        !! Temp. variable necessary for 4D interp.
 
!***************************************************************************************** 
! 
!   Find position of LUT-grid ssa value closest to actual ssa
! 
!***************************************************************************************** 
 
    rgfMinDiffSSA = abs(rgfLUTRInfSSAGrid1 - fSSALook) 
    rgfMinLoc = MinLoc(rgfMinDiffSSA) 
    iMinPosSSA = rgfMinLoc(1)
    if(rgfLUTRInfSSAGrid1(iMinPosSSA).eq.fSSALook) then
     iSmalerPosSSA = iMinPosSSA
     iLargerPosSSA = iMinPosSSA
    elseif(rgfLUTRInfSSAGrid1(iMinPosSSA).gt.fSSALook) then
     iSmalerPosSSA = iMinPosSSA-1
     iLargerPosSSA = iMinPosSSA
    else
     iSmalerPosSSA = iMinPosSSA
     iLargerPosSSA = iMinPosSSA+1
    endif
     
    if(iSmalerPosSSA.le.0) iSmalerPosSSA = 1
    if(iLargerPosSSA.gt.30) iLargerPosSSA = 30
    
!   Set nearest smaler/larger grid value of ssa
    fSmalerSSA = rgfLUTRInfSSAGrid1(iSmalerPosSSA) 
    fLargerSSA = rgfLUTRInfSSAGrid1(iLargerPosSSA) 


!***************************************************************************************** 
! 
!   Interpolate value of R_infinite with respect to szen,pzen,relazm for smaler ssa
! 
!***************************************************************************************** 
   
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
    float(rgiLUTRInf1(iSmalerPosSZen,iSmalerPosPZen,iSmalerPosRelAzm,iSmalerPosSSA))/1000.0
    dX100 = &
    float(rgiLUTRInf1(iLargerPosSZen,iSmalerPosPZen,iSmalerPosRelAzm,iSmalerPosSSA))/1000.0
    dX010 = &
    float(rgiLUTRInf1(iSmalerPosSZen,iLargerPosPZen,iSmalerPosRelAzm,iSmalerPosSSA))/1000.0
    dX110 = &
    float(rgiLUTRInf1(iLargerPosSZen,iLargerPosPZen,iSmalerPosRelAzm,iSmalerPosSSA))/1000.0
    dX001 = &
    float(rgiLUTRInf1(iSmalerPosSZen,iSmalerPosPZen,iLargerPosRelAzm,iSmalerPosSSA))/1000.0
    dX101 = &
    float(rgiLUTRInf1(iLargerPosSZen,iSmalerPosPZen,iLargerPosRelAzm,iSmalerPosSSA))/1000.0
    dX011 = &
    float(rgiLUTRInf1(iSmalerPosSZen,iLargerPosPZen,iLargerPosRelAzm,iSmalerPosSSA))/1000.0
    dX111 = &
    float(rgiLUTRInf1(iLargerPosSZen,iLargerPosPZen,iLargerPosRelAzm,iSmalerPosSSA))/1000.0
    call CLOUD_3DInterpolation &
         (dR,dS,dT,dX000,dX001,dX010,dX011,dX100,dX101,dX110,dX111, dSmalerRInf)


!***************************************************************************************** 
! 
!   Interpolate value of R_infinite with respect to szen,pzen,relazm for larger ssa
! 
!***************************************************************************************** 
   
    if(iSmalerPosSSA.eq.iLargerPosSSA) then
     dLargerRInf = dSmalerRInf
    else
     dX000 = &
     float(rgiLUTRInf1(iSmalerPosSZen,iSmalerPosPZen,iSmalerPosRelAzm,iLargerPosSSA))/1000.0
     dX100 = &
     float(rgiLUTRInf1(iLargerPosSZen,iSmalerPosPZen,iSmalerPosRelAzm,iLargerPosSSA))/1000.0
     dX010 = &
     float(rgiLUTRInf1(iSmalerPosSZen,iLargerPosPZen,iSmalerPosRelAzm,iLargerPosSSA))/1000.0
     dX110 = &
     float(rgiLUTRInf1(iLargerPosSZen,iLargerPosPZen,iSmalerPosRelAzm,iLargerPosSSA))/1000.0
     dX001 = &
     float(rgiLUTRInf1(iSmalerPosSZen,iSmalerPosPZen,iLargerPosRelAzm,iLargerPosSSA))/1000.0
     dX101 = &
     float(rgiLUTRInf1(iLargerPosSZen,iSmalerPosPZen,iLargerPosRelAzm,iLargerPosSSA))/1000.0
     dX011 = &
     float(rgiLUTRInf1(iSmalerPosSZen,iLargerPosPZen,iLargerPosRelAzm,iLargerPosSSA))/1000.0
     dX111 = &
     float(rgiLUTRInf1(iLargerPosSZen,iLargerPosPZen,iLargerPosRelAzm,iLargerPosSSA))/1000.0
     call CLOUD_3DInterpolation &
          (dR,dS,dT,dX000,dX001,dX010,dX011,dX100,dX101,dX110,dX111, dLargerRInf)
    endif

!***************************************************************************************** 
! 
!   Interpolate values of R_infinite for smaler/larger ssa with respect to actual ssa
! 
!***************************************************************************************** 

    if(iSmalerPosSSA.eq.iLargerPosSSA) then
     dRinf = dSmalerRInf
    else
!    Use sqrt(1-ssa) instead of ssa since R=R_no_abs-sqrt(1-ssa)*const
!    Keep in mind: larger SSA and therefore larger Rinf value computed above corresponds
!    to smaler sqrt of beta.
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
     print*,'Sqrt Beta:     ', dSqrtBeta
     print*,'Smaler Beta:   ', dSmalerSqrtBeta
     print*,'Larger Beta:   ', dLargerSqrtBeta
     print*,'Relative Beta: ', dRelativeSqrtBeta
     print*,'Rinf:          ', dRinf
     if(bPause) pause
    endif
    
    return
    end subroutine CLOUD_RInfinite1



    subroutine  CLOUD_RInfinite2
 
!***************************************************************************************** 
! 
!   Declaration of variables 
! 
!***************************************************************************************** 
 
    use SLALOM_Defs
    implicit none
 
    integer(2) :: iSmalerPosSZen
    integer(2) :: iLargerPosSZen
    integer(2) :: iSmalerPosPZen
    integer(2) :: iLargerPosPZen
    integer(2) :: iSmalerPosRelAzm
    integer(2) :: iLargerPosRelAzm
    
    real(8) :: dR
    real(8) :: dS
    real(8) :: dT
    real(8) :: dX000
    real(8) :: dX001
    real(8) :: dX010
    real(8) :: dX011
    real(8) :: dX100
    real(8) :: dX101
    real(8) :: dX110
    real(8) :: dX111
    real(8) :: dSmalerRInf !! RInf for smaler LUT-grid ssa value
    real(8) :: dLargerRInf  !! RInf for larger LUT-grid ssa value
    real(8) :: dSqrtBeta
    real(8) :: dSmalerSqrtBeta
    real(8) :: dLargerSqrtBeta
    real(8) :: dRelativeSqrtBeta
    
    real(4),dimension(1) :: rgfMinDiffSSA(30)   !! Temp. variable necessary for 4D interp.
    real(4),dimension(1) :: rgfMinLoc(1)        !! Temp. variable necessary for 4D interp.
 
!***************************************************************************************** 
! 
!   Find position of LUT-grid ssa value closest to actual ssa
! 
!***************************************************************************************** 
 
    rgfMinDiffSSA = abs(rgfLUTRInfSSAGrid2 - fSSALook) 
    rgfMinLoc = MinLoc(rgfMinDiffSSA) 
    iMinPosSSA = rgfMinLoc(1)
    if(rgfLUTRInfSSAGrid2(iMinPosSSA).eq.fSSALook) then
     iSmalerPosSSA = iMinPosSSA
     iLargerPosSSA = iMinPosSSA
    elseif(rgfLUTRInfSSAGrid2(iMinPosSSA).gt.fSSALook) then
     iSmalerPosSSA = iMinPosSSA-1
     iLargerPosSSA = iMinPosSSA
    else
     iSmalerPosSSA = iMinPosSSA
     iLargerPosSSA = iMinPosSSA+1
    endif
     
    if(iSmalerPosSSA.le.0) iSmalerPosSSA = 1
    if(iLargerPosSSA.gt.30) iLargerPosSSA = 30
    
!   Set nearest smaler/larger grid value of ssa
    fSmalerSSA = rgfLUTRInfSSAGrid2(iSmalerPosSSA) 
    fLargerSSA = rgfLUTRInfSSAGrid2(iLargerPosSSA) 


!***************************************************************************************** 
! 
!   Interpolate value of R_infinite with respect to szen,pzen,relazm for smaler ssa
! 
!***************************************************************************************** 
   
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
    float(rgiLUTRinf2(iSmalerPosSZen,iSmalerPosPZen,iSmalerPosRelAzm,iSmalerPosSSA))/1000.0
    dX100 = &
    float(rgiLUTRinf2(iLargerPosSZen,iSmalerPosPZen,iSmalerPosRelAzm,iSmalerPosSSA))/1000.0
    dX010 = &
    float(rgiLUTRinf2(iSmalerPosSZen,iLargerPosPZen,iSmalerPosRelAzm,iSmalerPosSSA))/1000.0
    dX110 = &
    float(rgiLUTRinf2(iLargerPosSZen,iLargerPosPZen,iSmalerPosRelAzm,iSmalerPosSSA))/1000.0
    dX001 = &
    float(rgiLUTRinf2(iSmalerPosSZen,iSmalerPosPZen,iLargerPosRelAzm,iSmalerPosSSA))/1000.0
    dX101 = &
    float(rgiLUTRinf2(iLargerPosSZen,iSmalerPosPZen,iLargerPosRelAzm,iSmalerPosSSA))/1000.0
    dX011 = &
    float(rgiLUTRinf2(iSmalerPosSZen,iLargerPosPZen,iLargerPosRelAzm,iSmalerPosSSA))/1000.0
    dX111 = &
    float(rgiLUTRinf2(iLargerPosSZen,iLargerPosPZen,iLargerPosRelAzm,iSmalerPosSSA))/1000.0
    call CLOUD_3DInterpolation &
         (dR,dS,dT,dX000,dX001,dX010,dX011,dX100,dX101,dX110,dX111, dSmalerRInf)


!***************************************************************************************** 
! 
!   Interpolate value of R_infinite with respect to szen,pzen,relazm for larger ssa
! 
!***************************************************************************************** 
   
    if(iSmalerPosSSA.eq.iLargerPosSSA) then
     dLargerRInf = dSmalerRInf
    else
     dX000 = &
     float(rgiLUTRinf2(iSmalerPosSZen,iSmalerPosPZen,iSmalerPosRelAzm,iLargerPosSSA))/1000.0
     dX100 = &
     float(rgiLUTRinf2(iLargerPosSZen,iSmalerPosPZen,iSmalerPosRelAzm,iLargerPosSSA))/1000.0
     dX010 = &
     float(rgiLUTRinf2(iSmalerPosSZen,iLargerPosPZen,iSmalerPosRelAzm,iLargerPosSSA))/1000.0
     dX110 = &
     float(rgiLUTRinf2(iLargerPosSZen,iLargerPosPZen,iSmalerPosRelAzm,iLargerPosSSA))/1000.0
     dX001 = &
     float(rgiLUTRinf2(iSmalerPosSZen,iSmalerPosPZen,iLargerPosRelAzm,iLargerPosSSA))/1000.0
     dX101 = &
     float(rgiLUTRinf2(iLargerPosSZen,iSmalerPosPZen,iLargerPosRelAzm,iLargerPosSSA))/1000.0
     dX011 = &
     float(rgiLUTRinf2(iSmalerPosSZen,iLargerPosPZen,iLargerPosRelAzm,iLargerPosSSA))/1000.0
     dX111 = &
     float(rgiLUTRinf2(iLargerPosSZen,iLargerPosPZen,iLargerPosRelAzm,iLargerPosSSA))/1000.0
     call CLOUD_3DInterpolation &
          (dR,dS,dT,dX000,dX001,dX010,dX011,dX100,dX101,dX110,dX111, dLargerRInf)
    endif

!***************************************************************************************** 
! 
!   Interpolate values of R_infinite for smaler/larger ssa with respect to actual ssa
! 
!***************************************************************************************** 

    if(iSmalerPosSSA.eq.iLargerPosSSA) then
     dRinf = dSmalerRInf
    else
!    Use sqrt(1-ssa) instead of ssa since R=R_no_abs-sqrt(1-ssa)*const
!    Keep in mind: larger SSA and therefore larger Rinf value computed above corresponds
!    to smaler sqrt of beta.
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
     print*,'Sqrt Beta:     ', dSqrtBeta
     print*,'Smaler Beta:   ', dSmalerSqrtBeta
     print*,'Larger Beta:   ', dLargerSqrtBeta
     print*,'Relative Beta: ', dRelativeSqrtBeta
     print*,'Rinf:          ', dRinf
     if(bPause) pause
    endif
    
    return
    end subroutine CLOUD_RInfinite2






!######################################################################################### 
!! 
!!  NAME 
!!  Subroutine "CLOUD_RdInfinite" 
!! 
!! 
!!  VERSION 
!!  $Revision: 1.7 $ 
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
!!  Please send any comments, suggestions, criticism, or (for our sake) bug reports to 
!!  info@lcrs.de  
!! 
!! 
!!  Copyright (c) 2003-2006 
!!  Alexander A. Kokhanovsky, Thomas Nauss
!! 
!!  This program is free software without any warranty and without even the implied  
!!  warranty of merchantability or fitness for a particular purpose. 
!!  If you use the software you must acknowledge the software and its author. 
!!
!######################################################################################### 
 
    subroutine  CLOUD_RdInfinite
 
!***************************************************************************************** 
! 
!   Declaration of variables 
! 
!***************************************************************************************** 
 
    use SLALOM_Defs
    implicit none
    real(4),dimension(1) :: rgfMinDiffMu0(100)  !! Temp. variable necessary for 2D interp.
    real(4),dimension(1) :: rgfMinDiffSSA(9)    !! Temp. variable necessary for 2D interp.
    real(4),dimension(1) :: rgfMinLoc(1)        !! Temp. variable necessary for 2D interp.
 
 
!***************************************************************************************** 
! 
!   Interpolate value of RDinf with respect to szen, ssa
! 
!***************************************************************************************** 

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
    rgfMinDiffSSA = abs(rgfLUTRdInfSSAGrid - fSSALook) 
    rgfMinLoc = MinLoc(rgfMinDiffSSA) 
    iMinPosSSA = rgfMinLoc(1) 
    if(rgfLUTRdInfSSAGrid(iMinPosSSA).eq.fSSALook) then
     iSmalerPosSSA = iMinPosSSA
     iLargerPosSSA = iMinPosSSA
    elseif(rgfLUTRdInfSSAGrid(iMinPosSSA).gt.fSSALook) then
     iSmalerPosSSA = iMinPosSSA-1
     iLargerPosSSA = iMinPosSSA
    else
     iSmalerPosSSA = iMinPosSSA
     iLargerPosSSA = iMinPosSSA+1
    endif
     
    if(iSmalerPosSSA.le.0) iSmalerPosSSA = 1
    if(iLargerPosSSA.gt.9) iLargerPosSSA = 9

!   Set nearest smaler/larger grid value of ssa
    fSmalerSSA = rgfLUTRdInfSSAGrid(iSmalerPosSSA) 
    fLargerSSA = rgfLUTRdInfSSAGrid(iLargerPosSSA) 

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
    call CLOUD_2DInterpolation(dRelativeMu0,dRelativeSSA,dx00,dx01,dx10,dx11,dRDInf) 

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
     if(bPause) pause
    endif

    return
    end subroutine CLOUD_RdInfinite


!######################################################################################### 
!! 
!!  NAME 
!!  Subroutine "CLOUD_Escape" 
!! 
!! 
!!  VERSION 
!!  $Revision: 1.7 $ 
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
!!  Please send any comments, suggestions, criticism, or (for our sake) bug reports to 
!!  info@lcrs.de  
!! 
!! 
!!  Copyright (c) 2003-2006 
!!  Alexander A. Kokhanovsky, Thomas Nauss
!! 
!!  This program is free software without any warranty and without even the implied  
!!  warranty of merchantability or fitness for a particular purpose. 
!!  If you use the software you must acknowledge the software and its author. 
!! 
!######################################################################################### 
 
    subroutine  CLOUD_Escape
 
!***************************************************************************************** 
! 
!   Declaration of variables 
! 
!***************************************************************************************** 
 
    use SLALOM_Defs
    implicit none
    real(4),dimension(1) :: rgfMinDiffMu0(30)   !! Temp. variable necessary for 2D interp.
    real(4),dimension(1) :: rgfMinDiffMu(30)    !! Temp. variable necessary for 2D interp.
    real(4),dimension(1) :: rgfMinDiffSSA(8)    !! Temp. variable necessary for 2D interp.
    real(4),dimension(1) :: rgfMinLoc(1)        !! Temp. variable necessary for 2D interp.
 

!***************************************************************************************** 
! 
!   Format
! 
!***************************************************************************************** 

950 format(a16,f15.7)    !! Print settings and results from 2D interpolation
 
!***************************************************************************************** 
! 
!   Interpolate value for K with respect to mu0, ssa
! 
!***************************************************************************************** 

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
    rgfMinDiffSSA = abs(rgfLUTKSSAGrid - fSSALook) 
    rgfMinLoc = MinLoc(rgfMinDiffSSA) 
    iMinPosSSA = rgfMinLoc(1) 
    if(rgfLUTKSSAGrid(iMinPosSSA).eq.fSSALook) then
     iSmalerPosSSA = iMinPosSSA
     iLargerPosSSA = iMinPosSSA
    elseif(rgfLUTKSSAGrid(iMinPosSSA).gt.fSSALook) then
     iSmalerPosSSA = iMinPosSSA-1
     iLargerPosSSA = iMinPosSSA
    else
     iSmalerPosSSA = iMinPosSSA
     iLargerPosSSA = iMinPosSSA+1
    endif
     
    if(iSmalerPosSSA.le.0) iSmalerPosSSA = 1
    if(iLargerPosSSA.gt.8) iLargerPosSSA = 8

!   Set nearest smaler/larger grid value of ssa
    fSmalerSSA = rgfLUTKSSAGrid(iSmalerPosSSA) 
    fLargerSSA = rgfLUTKSSAGrid(iLargerPosSSA) 

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

!***************************************************************************************** 
! 
!   Interpolated value for K with respect to mu, ssa
! 
!***************************************************************************************** 

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
     if(bPause) pause
    endif
    
    return
    end subroutine CLOUD_Escape


!######################################################################################### 
!! 
!!  NAME 
!!  Subroutine "CLOUD_1DInterpolation" 
!! 
!! 
!!  VERSION 
!!  $Revision: 1.7 $ 
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
!!  This code is from John Burkardt, http://www.csit.fsu.edu/~burkardt/f_src/f_src.html.
!!
!######################################################################################### 

    subroutine CLOUD_1DInterpolation(r,x0,x1,x)
    implicit none

    real ( kind = 8 ) r
    real ( kind = 8 ) x
    real ( kind = 8 ) x0
    real ( kind = 8 ) x1

    x = ( 1.0D+00 - r ) * x0 + r * x1

    return
    end subroutine CLOUD_1DInterpolation


!######################################################################################### 
!! 
!!  NAME 
!!  Subroutine "CLOUD_2DInterpolation" 
!! 
!! 
!!  VERSION 
!!  $Revision: 1.7 $ 
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
!!  This code is from John Burkardt, http://www.csit.fsu.edu/~burkardt/f_src/f_src.html.
!!
!######################################################################################### 

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


!######################################################################################### 
!! 
!!  NAME 
!!  Subroutine "CLOUD_3DInterpolation" 
!! 
!! 
!!  VERSION
!!  $Revision: 1.7 $
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
!!  Copyright (c) 2003-2005 
!!  This code is from John Burkardt, http://www.csit.fsu.edu/~burkardt/f_src/f_src.html.
!!
!######################################################################################### 

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


!######################################################################################### 
!! 
!!  NAME 
!!  Subroutine "FWARGUMENT" 
!! 
!! 
!!  VERSION 
!!  $Revision: 1.7 $ 
!! 
!! 
!!  PURPOSE 
!!  Handling of command line arguments.
!! 
!! 
!!  PROCEDURE 
!!  Handle command line arguments.
!! 
!!
!!  CONTACT 
!!  Please send any comments, suggestions, criticism, or (for our sake) bug reports to 
!!  info@lcrs.de  
!! 
!! 
!!  Copyright (c) 2003-2006 
!!  Jan Cermak
!! 
!!  This program is free software without any warranty and without even the implied  
!!  warranty of merchantability or fitness for a particular purpose. 
!!  If you use the software you must acknowledge the software and its author. 
!! 
!######################################################################################### 
    
    SUBROUTINE FWARGUMENT(iMode,chDefString,chArgument)
!    USE DFLIB
    CHARACTER(LEN=*) :: chDefString     !! Definition string for command line argument
    CHARACTER(LEN=*) :: chArgument      !! Command line argument
    INTEGER(2):: iMode                  !! Program mode (see comment)
    INTEGER(KIND=2) :: iStatus          !! Needed for GetArg subroutine
    INTEGER(KIND=4) :: liCArg            !! Counter for arguments
    INTEGER(KIND=2),PARAMETER :: iMaxArg = 50  !! Maximum number of arguments (fixed)
    chArgument = '' ! get rid of any initial values
    IF(iMode.EQ.1) THEN    ! chDefString is the preceeding argument
     find_pre: DO liCArg = 1,iMaxArg
     CALL GETARG(liCArg,chArgument)
    IF (chArgument.EQ.chDefString) EXIT find_pre
    ENDDO find_pre
     liCArg = liCArg + 1
     CALL GETARG(liCArg,chArgument)
     RETURN
    ELSEIF(iMode.EQ.2) THEN  ! chDefString is contained in the argument
     find_arg_substring: DO liCArg = 1,iMaxArg
     CALL GETARG(liCArg,chArgument)
     IF(INDEX(chArgument,chDefString).NE.0) RETURN
     ENDDO find_arg_substring
    ELSEIF(iMode.EQ.3) THEN  ! chDefString is identical with the argument
     find_arg_string: DO liCArg = 1,iMaxArg
     CALL GETARG(liCArg,chArgument)
     IF(chArgument.EQ.chDefString) RETURN
     ENDDO find_arg_string
    ELSEIF(iMode.EQ.4) THEN  ! chDefString contains the argument number
     READ(chDefString,*,ERR=250) liCArg
     CALL GETARG(liCArg,chArgument)
    ENDIF
    RETURN
250 iStatus = -1
    END SUBROUTINE FWARGUMENT
