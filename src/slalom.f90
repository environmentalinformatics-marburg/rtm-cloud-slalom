!###############################################################################
!!
!!  NAME
!!  Program "SLALOM"
!!
!!
!!  VERSION 
!!  $Revision: 1.7 $
!!
!!
!!  PURPOSE 
!!  Inverse algorithm for the retrieval of the cloud optical thickness (COT), 
!!  the effective cloud droplet radius (AEF), the liquid / ice water path (LWP),
!!  the cloud single scattering albedo (SSA) and the particle absorption length
!!  (PAL) for optically thick layers. The forward model used within SLALOM is
!!  CLOUD.
!!
!! 
!!  PROCEDURE 
!!  This program consists of three basic parts:
!!  1. Read/write and control routines (all in the main SLALOM program)
!!  2. Main retrieval and necessary functions (starting from SLALOM_Retrieve)
!!  3. Routines from program CLOUD used for the retrieval (starting from the 
!!     subroutine CLOUD)
!!  For operational or user specified applications, the retrieval algorithm can
!!  be included in other frameworks, as long as all necessary input parameters
!!  are given along with the calling sequence.
!!  For the forward computation (which is necessary for R_meassured-R_computed),
!!  all routines of the CLOUD program are included in SLALOM. The only
!!  modifications are:
!!  1. The variables are declared in a special section of the module SLALOM_defs
!!  2. The main routine is included as a subroutine (without control sequence) 
!!  Configuration of this version: You don't need a configuration file,
!!  even there is a routine for reading it (but it is not called). This is just
!!  for an easier copy/paste of new versions of CLOUD in this code - so don't
!!  worry about that.
!!
!!
!!  ADDITIONAL REQUIREMENTS
!!  none
!!
!!  ATTENTION 
!!  If ASCII input is used, the ASCII file must contain the following parameters
!!  in exactly the same order:
!!  sun zenith, pixel zenith, relative azimuth, 
!!  reflection non-absorbing channel, and reflection absorbing channel.
!!
!!
!!  CALLING SEQUENCE 
!!  Usage if binary input (bi = 1):
!!  slalom <bi> <cp> <w1> <w2> <lut1> <lut2> \
!!  <band1> <band2> <szen> <pzen> <razm> <cmask> <mask> \
!!  <tau> <aef> <lwp> <iwp> <pal> <cols> <rows> \
!!  <alb> <alb1> <alb2>
!!
!!  Usage if ASCII input (bi = 0):
!!  slalom <bi> <cp> <w1> <w2> <lut1> <lut2> <indat> <outdat>
!!
!!  Parameters:
!!  <bi>    Binary input (1 = binary, 0 = ASCII).
!!  <cp>    Cloud phase (1 = water, 0 = ice).
!!  <w1>    Wavelength (microns) of the non-absorbing channel (x.xxx).
!!  <w2>    Wavelength (microns) of the absorbing channel (x.xxx).
!!  <lut1>  Name of the LUT used for Rinf for the non-abs. channel.
!!  <lut2>  Name of the LUT used for Rinf for the absorbing channel.
!!  <band1> Filename for the non-absorbing channel data (real 4)
!!  <band2> Filename for the absorbing channel data (real 4).
!!  <szen>  Filename for the sun zenith angle data (real 4).
!!  <pazm>  Filename for the pixel zenith data (real 4).
!!  <razm>  Filename for the relative azimuth angle data (real 4).
!!  <cmask> Filename for the cloud mask (integer 2).
!!  <mask>  Filename for the computation mask (integer 2).
!!  <tau>   Filename for the optical thickness (output).
!!  <aef>   Filename for the effective droplet radius (output)
!!  <lwp>   Filename for the liquid water path (output)
!!  <iwp>   Filename for the ice water path (output)
!!  <ssa>   Filename for the single scattering albedo (output)
!!  <pal>   Filename for the particle absorption length (output)
!!  <cols>  Number of columns in the input binary datasets.
!!  <rows>  Number of rows in the input binary datasets.
!!  <alb>   Use albedo different from 0.0 (1 = yes, 0 = no).
!!  <alb1>  Filename for the non-absorbing albedo data (real 4)
!!  <alb2>  Filename for the absorbing albedo data (real 4)
!!  <indat> Filename for the input ASCII dataset
!!  <outdat> Filename for the output ASCII dataset
!!
!!
!!  COMMENT 
!! 
!!  LANGUAGE 
!!  FORTRAN 90/95 
!!
!!
!!  COMPILER 
!!  Any intel compiler should do
!!
!!
!!  CVS INFORMATION 
!!  $Id: slalom.f90,v 1.7 2008/08/26 14:21:51 tnauss Exp $ 
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
!!  $Log: slalom.f90,v $
!!  Revision 1.7  2008/08/26 14:21:51  tnauss
!!  Update
!!  Include cloud geometrical thickness and droplet concentration retrieval.
!!
!!  Revision 1.6  2008/08/26 07:06:35  tnauss
!!  Update
!!  Include test function in SLALOM code (i. e. known tau and aef is read from
!!  the ASCII input file in addition to the measurement data and error between
!!  tau and aef input and output is computed in addition).
!!
!!  Revision 1.5  2008/05/20 10:47:06  tnauss
!!  Update
!!  Include ASCII file in- and output.
!!
!!  Revision 1.4  2008/05/19 17:18:01  tnauss
!!  Update
!!  Remove test codelines, include command line calling and include some more
!!  documentation.
!!
!!  Revision 1.3  2008/04/15 18:37:56  tnauss
!!  Major update
!!  Include pointers for a seamless integration of the CLOUD program sources
!!  (starting from CLOUD 1.41).
!!
!!  Revision 1.2  2007/11/29 15:09:31  tnauss
!!  Update
!!
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
!###############################################################################

!###############################################################################
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
!!  Please send any comments, suggestions, criticism, or (for our sake) bug
!!  reports to nauss@lcrs.de
!! 
!!  Development
!!  2003-2008 
!!  Thomas Nauss, Alexander A. Kokhanovsky
!! 
!!  This program is free software without any warranty and without even the
!!  implied warranty of merchantability or fitness for a particular purpose. 
!!  If you use the software you must acknowledge the software and its authors. 
!! 
!############################################################################### 

    module SLALOM_defs
 
!******************************************************************************* 
! 
!   Declaration of variables for SLALOM
! 
!******************************************************************************* 

!   Only necessary for testing - start
    real(4), dimension(1) :: prgfTauIn(30000)
    real(4), dimension(1) :: prgfAefIn(30000)
    real(4), dimension(1) :: prgfTauError(30000)
    real(4), dimension(1) :: prgfAefError(30000)
    real(4), dimension(1) :: prgft1In(30)
    real(4) :: fTauIn
    real(4) :: fAefIn
    real(4) :: ft1In
!   Only necessary for testing - end

    logical(4) :: bError       !! Error indicator
    logical(4) :: bBinary      !! Use binary (1) or ASCII (0) input
    logical(4) :: bWaterCloud  !! Retrieval for water (1) or ice (0) clouds
    logical(4) :: bAlbedo      !! Use background albedo different from 0.0
    logical(4) :: bTest        !! Testrun (known tau and aef are read)

    character(1) :: chBinary       !! Command line argument for bBinary
    character(1) :: chWaterCloud   !! Command line argument for bWaterCloud
    character(1) :: chAlbedo       !! Command line argument for bAlbedo
    character(1) :: chCloudType    !! Command line argument for iCloudType
    character(1) :: chTest         !! Command line argument for bTest

    character(5) :: chWavelength01 !! Command line argument for iWavelength01
    character(5) :: chWavelength02 !! Command line argument for iWavelength02

    character(10) :: chCols    !! Command line argument liCols
    character(10) :: chRows    !! Command line argument for liRows

    character(300) :: chBand01 !! Name of non-absorbing channel dataset
    character(300) :: chBand02 !! Name of absorbing channel dataset
    character(300) :: chSZen   !! Name of sun zenith angle dataset
    character(300) :: chPZen   !! Name of pixel zenith angle dataset
    character(300) :: chRelAzm !! Name of relative azimuth angle dataset
    character(300) :: chCloudmask !! Name of cloud mask dataset
    character(300) :: chComputationMask !! Name of computation mask dataset
    character(300) :: chTau   !! Name of optical thickness dataset
    character(300) :: chAef   !! Name of effective radius dataset
    character(300) :: chLWP   !! Name of liquid water path dataset
    character(300) :: chIWP   !! Name of ice water path dataset
    character(300) :: chCloudAlb !! Name of cloud albedo dataset
    character(300) :: chCGT      !! Name of cloud geometrical thickness dataset
    character(300) :: chCDC      !! Name of columnar droplet concentration dataset
    character(300) :: chDC       !! Name of droplet concentration dataset
    character(300) :: chSSA   !! Name of single scattering albedo dataset
    character(300) :: chPAL   !! Name of particle absorption length dataset
    character(300) :: chAlbedo01 !! Name of non-absorbing albedo dataset
    character(300) :: chAlbedo02 !! Name of absorbing albedo dataset
    character(300) :: chInfile   !! Name of ASCII input dataset
    character(300) :: chOutfile  !! Name of ASCII output dataset

    integer(2) :: iError !! Flag for input/output error
    integer(2) :: iCloud !! Cloud mask flag
    integer(2) :: iCloudType !! Cloud type (Cu/Sc, As/Ac, Ci/Cs/Cc, Ns, Cb,St)

    integer(2),allocatable :: prgiCloud(:) !! Cloud mask (>=1: cloud)
    integer(2),allocatable :: prgiMask(:)  !! Computation mask (1=compute)

    integer(4) :: liCols !! Number of columns in binary input files
    integer(4) :: liRows !! Number of rows in binary input files
    integer(4) :: liRuns !! Number of input value sets
    integer(4) :: liRecl !! Recordlength for satellite data binary input files

    real(8) :: fDummy        !! Dummy
    real(8) :: fNonAbs       !! Actual non-abs. band value
    real(8) :: fAbs          !! Actual abs. band value
    real(8) :: fNonAbsAlb    !! Actual non-abs. band background albedo
    real(8) :: fAbsAlb       !! Actual abs. band background albedo
    real(8) :: fSSAMinValid  !! Minimum of SSA for finding root of functionR
    real(8) :: fSSAMaxValid  !! Maximum of SAA for finding root for functionR
    real(8) :: fAefMinValid  !! Minimum of aef for finding root of functionR
    real(8) :: fAefMaxValid  !! Maximum of aef for finding root for functionR
    real(8) :: fdTau         !! Ratio between retrieved and correct tau
    real(8) :: fdBeta        !! Ratio between retrieved and correct beta
    real(8) :: fWavelength01 !! Non-absorbing channel wavelength
    real(8) :: fWavelength02 !! Absorbing channel wavelength
    real(8) :: fMi           !! Imaginary part of refractive index
    real(8) :: fPi           !! Pi
    real(8) :: fBetaInf      !! Probability of absorption for an infinite cloud
    real(8) :: fPAL          !! Particle absorption length
    real(8) :: fAef          !! Effective cloud droplet radius
    real(8) :: fIWP          !! Ice water path
    real(8) :: fCloudAlb     !! Cloud albedo
    real(8) :: fCGT          !! Cloud geometrical thickness
    real(8) :: fCDC          !! Columnar droplet concentration
    real(8) :: fDC           !! Droplet concentration
    real(8) :: fXef          !! Size parameter
    real(8) :: fLWP          !! Liquid water path
    real(8) :: ft1           !! Global reflectance for wavelength 1
    real(8) :: fRho          !! Density of water
    real(8) :: fSigma_ext    !! Extinction coefficient
    real(8) :: fSigma_abs    !! Absorbtion coefficient
    real(8) :: dR0           !! Reflection of a semi-infinite cloud

    real(4),allocatable :: prgfNonAbs(:) !! Input non-abs. band reflectance
    real(4),allocatable :: prgfAbs(:)    !! Input abs. band reflectance
    real(4),allocatable :: prgfNonAbsAlb(:) !! Input non-abs. background albedo
    real(4),allocatable :: prgfAbsAlb(:) !! Input abs. background albedo
    real(4),allocatable :: prgfSZen(:)   !! Input solar zenith angle
    real(4),allocatable :: prgfPZen(:)   !! Input pixel zenith angle
    real(4),allocatable :: prgfRelAzm(:) !! Input relative azimuth angle
    real(4),allocatable :: prgfTau(:)    !! Output optical thickness
    real(4),allocatable :: prgfAef(:)    !! Output cloud droplet radius
    real(4),allocatable :: prgfLWP(:)    !! Output liquid water path
    real(4),allocatable :: prgfSSA(:)    !! Output single scattering albedo
    real(4),allocatable :: prgfPAL(:)    !! Output particle absorption length
    real(4),allocatable :: prgfCloudAlb(:) !! Output cloud albedo
    real(4),allocatable :: prgfCGT(:)    !! Output cloud geometrical thickness
    real(4),allocatable :: prgfCDC(:)    !! Output columnar droplet concentration
    real(4),allocatable :: prgfDC(:)     !! Output droplet concentration
    real(4),dimension(1) :: rgfLUTRInfSSAGridIn(30) !! Grid of lut_RInf.dat
    real(4),dimension(1) :: rgfLUTKSSAGridIn(8)  !! Grid of lut_k.dat
    real(4),dimension(1) :: rgfLUTRdInfSSAGridIn(9) !! Grid of lut_RdInf.dat
    real(4),dimension(1) :: rgfAefInput(2748620)

    real(8) :: drgC1630(5)  !! Parameters for wavelength 1.6300µm
    real(8) :: drgD(5)      !! Parameter
    real(8) :: drgCWC(6)    !! Cloud water content (Cu/Sc, As/Ac, Ci/Cs/Cc, Ns, Cb,St)
    real(8) :: dSSARetrieve !! Retrieved value of SSA

!   Some pointers and targets
    character(300),target ::  chLUT_RInf02 !! Name of the LUT for Rinf02
    integer(4),target  ::     rgiLUTRInf02(90,90,181,30) !! LUT for Rinf02
    real(4),target ::  fgLUT02 !! Asymmetry parameter used for LUTRInf02
    real(4),target ::  rgfLUTRInfSSAGrid02(30) !! Grid for LUTRInf02
    real(8),target ::  drgC0645(5)  !! Parameters for wavelength 0.6457µm
    real(8),target ::  drgC0859(5)  !! Parameters for wavelength 0.8590µm
    real(8),pointer :: podrgC(:)    !! Pointer to parameters

!   Modified for next version 2008-07-16
!    real(4),dimension(1) :: rgfLUTgNonAbs(28)   !! LUT for g for non-abs. channel
!    real(4),dimension(1) :: rgfLUTgAbs(28)   !! LUT for g for abs. channel
!   Modified for next version 2008-07-16

!***************************************************************************************** 
! 
!   Declaration of variables and data blocks from program CLOUD
! 
!***************************************************************************************** 

    logical(4) :: bPrint        !! Print intermediate results
    logical(4) :: bPause        !! Hold after printing intermediate results
    logical(4) :: bFinite       !! Compute finite (1) or semi-infinite (0) media

    character(20) :: chVersion    !! Version of the program
    character(300) :: chControlFile !! Name of the control file

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
    real(4),dimension(1) :: rgfLUTKSSAGrid(8)  !! Grid (ssa) of lut_k
    real(4),dimension(2) :: rgfLUTKSSA(30,8)   !! Grid (mu0,ssa) from lut_k
    real(4),dimension(1) :: rgfLUTKMuGrid(30)  !! Grid (mu) of lut_k
    real(4),dimension(1) :: rgfLUTRdInfMu0Grid(100) !! Grid (mu0) of lut_RdInf
    real(4),dimension(1) :: rgfLUTRdInfSSAGrid(9) !! Grid (ssa) lut_RdInf
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

!Uncomment the following lines for theoretical error studies
!Start change
    character(8) :: STIME
!!!    real(4) :: fTauTarget, fAefTarget
!End change

!   Grid of single-scattering albedo array in LUT K
    data rgfLUTRInfSSAGrid/0.800,0.810,0.820,0.830,0.840,0.850,0.860,0.870, &
                           0.880,0.890,0.900,0.910,0.920,0.930,0.940,0.950, &
                           0.960,0.970,0.980,0.990,0.991,0.992,0.993,0.994, &
                           0.995,0.996,0.997,0.998,0.999,1.000/
    data rgfLUTKSSAGrid/0.800,0.900,0.950,0.980,0.990,0.995,0.999,1.000/ 
    data rgfLUTRdInfSSAGrid/0.500,0.600,0.700,0.800,0.900,0.950,0.990,0.999, &
                            0.9999/ 

!*******************************************************************************
! 
!   Set data blocks from SLALOM
! 
!*******************************************************************************

    data drgC0645/0.1121,0.5118,0.8997,0.0,0.0/ ! 0.6457 microns
    data drgC0859/0.1115,0.4513,1.2719,0.0,0.0/ ! 0.8590 microns
    data drgC1630/0.0608,2.465,-32.98,248.94,-636.0/ ! 1.630 microns
    data drgD/1.671,0.0025,-2.365e-4,2.861e-6,-1.05e-8/

!   Cloud liquid water content from Kawamoto et al. 2001 after
!   Pruppacher et al. 1978. Ci,Cs,Cc (Pos. 3) and Cb (Pos. 5) are from Nakajima
!   et al. 1995 after Liou 1976.
!   Cu/Sc, As/Ac, Ci/Cs/Cc, Ns, Cb,St
    data drgCWC/0.3,0.25,0.014,0.3,0.3928,0.35/

!   Modified for next version 2008-07-16
!    data rgfLUTgNonAbs/0.8093E+00,0.8285E+00,0.8383E+00,0.8445E+00,0.8491E+00, &
!                       0.8527E+00,0.8556E+00,0.8581E+00,0.8602E+00,0.8620E+00, &
!                       0.8636E+00,0.8649E+00,0.8662E+00,0.8672E+00,0.8682E+00, &
!                       0.8691E+00,0.8699E+00,0.8706E+00,0.8713E+00,0.8719E+00, &
!                       0.8724E+00,0.8730E+00,0.8734E+00,0.8739E+00,0.8743E+00, &
!                       0.8747E+00,0.8751E+00,0.8754E+00/
!    data rgfLUTgAbs/0.7868E+00,0.7861E+00,0.8026E+00,0.8176E+00,0.8285E+00, &
!                       0.8362E+00,0.8419E+00,0.8463E+00,0.8498E+00,0.8528E+00, &
!                       0.8554E+00,0.8576E+00,0.8596E+00,0.8614E+00,0.8631E+00, &
!                       0.8646E+00,0.8659E+00,0.8672E+00,0.8683E+00,0.8694E+00, &
!                       0.8704E+00,0.8713E+00,0.8722E+00,0.8730E+00,0.8737E+00, &
!                       0.8745E+00,0.8752E+00,0.8758E+00/
!   Modified for next version 2008-07-16

    end module SLALOM_defs

!############################################################################### 
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
!!  Please send any comments, suggestions, criticism, or (for our sake) bug
!!  reports to nauss@lcrs.de
!! 
!!  Development
!!  2003-2008 
!!  Thomas Nauss, Alexander A. Kokhanovsky
!! 
!!  This program is free software without any warranty and without even the
!!  implied warranty of merchantability or fitness for a particular purpose. 
!!  If you use the software you must acknowledge the software and its authors. 
!! 
!###############################################################################


    program SLALOM
 
!*******************************************************************************
! 
!   Declaration of variables 
! 
!*******************************************************************************

    use SLALOM_defs
    implicit none

    integer(4) :: liCounter !! Counter
    integer(4) :: liErrorCounter !! Counter used for sensitivity study

!******************************************************************************* 
! 
!   Format 
! 
!******************************************************************************* 

!! Read input data file content
950 format(3f7.2,f13.9,f8.5,f11.6,2f11.8)
960 format('SZen   PZen   RelAzm Rna     Ra      Ana      Aa      TauIn   ' &
           'SSAIn   Tau    SSA     PAL     Aef     LWP          dTau         ' &
           'dBeta        error')

!! Write output data file content
970 format(3f7.2,4f8.5,f8.2,f8.5,f8.2,f8.5,3f8.2,f13.4,f13.4,f13.4)

!   Only necessary for testing - start
980 format('# SZen    PZen      RelAzm    Tau       Aef       TauIn     AefIn     TauError  AefError')
990 format(9f10.3)
!   Only necessary for testing - end

!Uncomment the following lines for theoretical error studies
!Start change
!!!999 format(9f10.3)
!End change

!******************************************************************************* 
! 
!   Say hello...
! 
!*******************************************************************************

    chVersion = '$Revision: 1.7 $'
    print*, ' '
    print*, ' '
    print*, '##################################################################'
    print*, 'Program SLALOM'
    print*, chVersion
    print*, ' '
    print*, 'Forward model CLOUD'
    print*, 'Revision: 1.43'
    print*, ' '
    print*, '2006-2008, Thomas Nauss, Alexander A. Kokhanovsky'
    print*, ' '
    print*, ' '
    print*, 'Usage if binary input (bi = 1):'
    print*, 'slalom <bi> <cp> <w1> <w2> <lut1> <lut2> \'
    print*, '<band1> <band2> <szen> <pzen> <razm> <cmask> <mask> \'
    print*, '<tau> <aef> <lwp> <iwp> <pal> <cols> <rows> \'
    print*, '<ctype> <alb> <alb1> <alb2>'
    print*, ' '
    print*, 'Usage if ASCII input (bi = 0):'
    print*, 'slalom <bi> <cp> <w1> <w2> <lut1> <lut2> <ctype> <indat> <outdat> <test>'
    print*, ' '
    print*, 'Parameters:'
    print*, '<bi>    Binary input (1 = binary, 0 = ASCII).'
    print*, '<cp>    Cloud phase (1 = water, 0 = ice).'
    print*, '<w1>    Wavelength (microns) of the non-absorbing channel (x.xxx).'
    print*, '<w2>    Wavelength (microns) of the absorbing channel (x.xxx).'
    print*, '<lut1>  Name of the LUT used for Rinf for the non-abs. channel.'
    print*, '<lut2>  Name of the LUT used for Rinf for the absorbing channel.'
    print*, '<band1> Filename for the non-absorbing channel data (real 4)'
    print*, '<band2> Filename for the absorbing channel data (real 4).'
    print*, '<szen>  Filename for the sun zenith angle data (real 4).'
    print*, '<pazm>  Filename for the pixel zenith data (real 4).'
    print*, '<razm>  Filename for the relative azimuth angle data (real 4).'
    print*, '<cmask> Filename for the cloud mask (integer 2).'
    print*, '<mask>  Filename for the computation mask (integer 2).'
    print*, '<tau>   Filename for the optical thickness (output).'
    print*, '<aef>   Filename for the effective droplet radius (output)'
    print*, '<lwp>   Filename for the liquid water path (output)'
    print*, '<iwp>   Filename for the ice water path (output)'
    print*, '<ssa>   Filename for the single scattering albedo (output)'
    print*, '<pal>   Filename for the particle absorption length (output)'
    print*, '<calb>  Filename for the cloud albedo (output)'
    print*, '<cgt>   Filename for the cloud geometrical thickness (output)'
    print*, '<cdc>   Filename for the columnar droplet concentration (output)'
    print*, '<dc>    Filename for the droplet concentration (output)'
    print*, '<cols>  Number of columns in the input binary datasets.'
    print*, '<rows>  Number of rows in the input binary datasets.'
    print*, '<ctype> Cloud type'
    print*, '        (1=Cu/Sc, 2=As/Ac, 3=Ci/Cs/Cc, 4=Ns, 5=Cb,6=St).'
    print*, '<alb>   Use albedo different from 0.0 (1 = yes, 0 = no).'
    print*, '<alb1>  Filename for the non-absorbing albedo data (real 4)'
    print*, '<alb2>  Filename for the absorbing albedo data (real 4)'
    print*, '<indat> Filename for the input ASCII dataset'
    print*, '<outdat> Filename for the output ASCII dataset'
    print*, '<test>  Testrun flag (error of retrieved parameters is computed)'
    print*, ' '
    print*, '##################################################################'
    print*, ' '
    print*, ' '

!*******************************************************************************
!
!   Set configuration parameters
!
!*******************************************************************************

!   Set parameters
    bPause         = .FALSE.
    bPrint         = .FALSE.
    bError         = .FALSE.
    bAlbedo        = .FALSE.
    bFinite        = .TRUE.
    fRho           = 1.0
    fSSAMinValid   = 0.8d0
    fSSAMaxValid   = 1.0d0
    fAefMinValid   = 3.0
    fAefMaxValid   = 40.0
    fPi            = acos(-1.)

!******************************************************************************* 
!
!   Read command line parameters
!
!******************************************************************************* 

!   Read command line parameters
!   Binary input (binary = 1, ASCII = 0)
    call GETARG(1, chBinary)
    read(chBinary, '(i1)') bBinary
    if(chBinary.eq.'') then 
     print*, 'No input mode selected <bi>.'
     bError=.TRUE.
    endif

!   Binary input
    if(bBinary) then
!    Cloud phase (water = 1, ice = 0)
     call GETARG(2, chWaterCloud)
     read(chWaterCloud, '(i1)') bWaterCloud
     if(chWaterCloud.eq.'') then
      print*, 'No cloud phase given <cp>.'
      bError=.TRUE.
     endif

!    Non-absorbing wavelength
     call GETARG(3, chWavelength01)
     if(chWavelength01.eq.'') then
      print*, 'No wavelength given <w1>.'
      bError=.TRUE.
     else
      read(chWavelength01, '(f5.3)') fWavelength01
     endif

!    Absorbing wavelength
     call GETARG(4, chWavelength02)
     if(chWavelength02.eq.'') then
      print*, 'No wavelength given <w2>.'
      bError=.TRUE.
     else
      read(chWavelength02, '(f5.3)') fWavelength02
     endif

!    LUT Rinf for non-absorbing wavelength
     call GETARG(5, chLUT_RInf)
     if(chLUT_RInf.eq.'') then
      print*, 'No LUT for Rinf given <lut1>.'
      bError=.TRUE.
     endif

!    LUT Rinf for absorbing wavelength
     call GETARG(6, chLUT_RInf02)
     if(chLUT_RInf02.eq.'') then
      print*, 'No LUT for Rinf given <lut2>.'
      bError=.TRUE.
     endif

!    Filename for the non-absorbing channel data (real 4)
     call GETARG(7, chBand01)
     if(chBand01.eq.'') then
      print*, 'No input channel given <band1>.'
      bError=.TRUE.
     endif

!    Filename for the absorbing channel data (real 4).
     call GETARG(8, chBand02)
     if(chBand02.eq.'') then
      print*, 'No input channel given <band2>.'
      bError=.TRUE.
     endif

!    Filename for the pixel zenith angle data (real 4).
     call GETARG(9, chSZen)
     if(chSZen.eq.'') then
      print*, 'No sun zenith angle given <szen>.'
      bError=.TRUE.
     endif

!    Filename for the pixel zenith data (real 4).
     call GETARG(10, chPZen)
     if(chPZen.eq.'') then
      print*, 'No pixel zenith angle given <pzen>.'
      bError=.TRUE.
     endif

!    Filename for the relative azimuth angle data (real 4).
     call GETARG(11, chRelAzm)
     if(chRelAzm.eq.'') then
      print*, 'No relative azimuth angle given <razm>.'
      bError=.TRUE.
     endif

!    Filename for the cloud mask (integer 2).
     call GETARG(12, chCloudMask)
     if(chCloudMask.eq.'') then
      print*, 'No cloud mask given <cmask>.'
      bError=.TRUE.
     endif

!    Filename for the computation mask (integer 2).band1
     call GETARG(13, chComputationMask)
     if(chComputationMask.eq.'') then
      print*, 'No computaion mask given <mask>.'
      bError=.TRUE.
     endif

!    Filename for the optical thickness (output).
     call GETARG(14, chTau)
     if(chTau.eq.'') then
      print*, 'No tau output file given <tau>.'
      bError=.TRUE.
     endif

!    Filename for the effective droplet radius (output)
     call GETARG(15, chAef)
     if(chAef.eq.'') then
      print*, 'No aef output file given <aef>.'
      bError=.TRUE.
     endif

!    Filename for the liquid water path (output)
     call GETARG(16, chLWP)
     if(chLWP.eq.'') then
      print*, 'No lwp output file given <lwp>.'
      bError=.TRUE.
     endif

!    Filename for the ice water path (output)
     call GETARG(17, chIWP)
     if(chIWP.eq.'') then
      print*, 'No iwp output file given <iwp>.'
      bError=.TRUE.
     endif

!    Filename for the particle absorption length (output)
     call GETARG(18, chSSA)
     if(chSSA.eq.'') then
      print*, 'No ssa output file given <ssa>.'
      bError=.TRUE.
     endif

!    Filename for the particle absorption length (output)
     call GETARG(19, chPAL)
     if(chPAL.eq.'') then
      print*, 'No pal output file given <pal>.'
      bError=.TRUE.
     endif

!    Filename for the particle absorption length (output)
     call GETARG(20, chCloudAlb)
     if(chCloudAlb.eq.'') then
      print*, 'No cloud albedo output file given <calb>.'
      bError=.TRUE.
     endif

!    Filename for the cloud geometrical thickness (output)
     call GETARG(21, chCGT)
     if(chCloudAlb.eq.'') then
      print*, 'No cloud geometrical thickness output file given <calb>.'
      bError=.TRUE.
     endif

!    Filename for the columnar droplet concentration (output)
     call GETARG(22, chCDC)
     if(chCloudAlb.eq.'') then
      print*, 'No columnar droplet concentration output file given <calb>.'
      bError=.TRUE.
     endif

!    Filename for the droplet concentration (output)
     call GETARG(23, chDC)
     if(chCloudAlb.eq.'') then
      print*, 'No droplet concentration output file given <calb>.'
      bError=.TRUE.
     endif

!    Number of columns in the input binary dataset
     call GETARG(24, chCols)
     if(chCols.eq.'') then 
      print*, 'No coloumns given <cols>.'
      bError=.TRUE.
     else
      read(chCols, *) liCols
     endif

!    Number of rows in the input binary dataset
     call GETARG(25, chRows)
     if(chRows.eq.'') then
      print*, 'No rows given <rows>.'
      bError=.TRUE.
     else
     read(chRows, *) liRows
     endif

!    Cloud type
     call GETARG(26, chCloudType)
     if(chCloudType.eq.'') then 
      print*, 'No cloud type given <ctype>.'
      bError=.TRUE.
     else
      read(chCloudType, *) iCloudType
     endif

!    Use background albedo different from 0.0
     call GETARG(27, chAlbedo)
     read(chAlbedo, '(i1)') bAlbedo
     if(chAlbedo.eq.'') then
      print*, 'No albedo mode selected <alb>.'
      bError=.TRUE.
     endif

     if(bAlbedo) then
!     Filename for the non-absorbing albedo data (real 4)
      call GETARG(28, chAlbedo01)
      if(chAlbedo01.eq.'') then
       print*, 'No albedo for band 1 given <alb1>.'
       bError=.TRUE.
      endif

!     Filename for the absorbing channel albedo (real 4).
      call GETARG(29, chAlbedo02)
      if(chAlbedo02.eq.'') then
       print*, 'No albedo for band 2 given <alb2>.'
       bError=.TRUE.
      endif
     endif

!   ASCII input
    elseif(.not.bBinary) then
!    Cloud phase (water = 1, ice = 0)
     call GETARG(2, chWaterCloud)
     read(chWaterCloud, '(i1)') bWaterCloud
     if(chWaterCloud.eq.'') then
      print*, 'No cloud phase given <cp>.'
      bError=.TRUE.
     endif

!    Non-absorbing wavelength
     call GETARG(3, chWavelength01)
     if(chWavelength01.eq.'') then
      print*, 'No wavelength given <w1>.'
      bError=.TRUE.
     else
      read(chWavelength01, '(f5.3)') fWavelength01
     endif

!    Absorbing wavelength
     call GETARG(4, chWavelength02)
     if(chWavelength02.eq.'') then
      print*, 'No wavelength given <w2>.'
      bError=.TRUE.
     else
      read(chWavelength02, '(f5.3)') fWavelength02
     endif

!    LUT Rinf for non-absorbing wavelength
     call GETARG(5, chLUT_RInf)
     if(chLUT_RInf.eq.'') then
      print*, 'No LUT for Rinf given <lut1>.'
      bError=.TRUE.
     endif

!    LUT Rinf for absorbing wavelength
     call GETARG(6, chLUT_RInf02)
     if(chLUT_RInf02.eq.'') then
      print*, 'No LUT for Rinf given <lut2>.'
      bError=.TRUE.
     endif

!    Cloud type
     call GETARG(7, chCloudType)
     if(chCloudType.eq.'') then 
      print*, 'No cloud type given <ctype>.'
      bError=.TRUE.
     else
      read(chCloudType, *) iCloudType
     endif

!    Filename for the ASCII input file
     call GETARG(8, chInfile)
     if(chInfile.eq.'') then
      print*, 'No input file given <indat>.'
      bError=.TRUE.
     endif

!    Filename for the ASCII output file
     call GETARG(9, chOutfile)
     if(chOutfile.eq.'') then
      print*, 'No output file given <outdat>.'
      bError=.TRUE.
     endif

!    Flag for testruns
     call GETARG(10, chTest)
     read(chTest, '(i1)') bTest
     if(chTest.eq.'') then
      print*, 'No test flag setting given <indat>.'
      bError=.TRUE.
     endif
    endif

    if(bError) then
     print*, ' '
     print*, 'The above errors have occured.'
     print*, 'The program will be stoped...'
     stop
    else
     print*, 'The following settings will be used:'
     print*, '<test>  = ', bTest
     print*, '<cp>    = ', bWaterCloud
     print*, '<w1>    = ', fWavelength01
     print*, '<w2>    = ', fWavelength02
     print*, '<lut1>  = ', trim(chLUT_RInf)
     print*, '<lut2>  = ', trim(chLUT_RInf02)
     if(bBinary) then
      print*, '<band1> = ', trim(chBand01)
      print*, '<band2> = ', trim(chBand02)
      print*, '<szen>  = ', trim(chSZen)
      print*, '<pzen>  = ', trim(chPZen)
      print*, '<razm>  = ', trim(chRelAzm)
      print*, '<cmask> = ', trim(chCloudMask)
      print*, '<mask>  = ', trim(chComputationMask)
      print*, '<tau>   = ', trim(chTau)
      print*, '<aef>   = ', trim(chAef)
      print*, '<lwp>   = ', trim(chLWP)
      print*, '<iwp>   = ', trim(chIWP)
      print*, '<ssa>   = ', trim(chSSA)
      print*, '<pal>   = ', trim(chPAL)
      print*, '<calb>   = ', trim(chCloudAlb)
      print*, '<cols>  = ', liCols
      print*, '<rows>  = ', liRows
      print*, '<ctype> = ', iCloudType
      print*, '<alb>   = ', bAlbedo
      if(bAlbedo) then
       print*, '<alb1>  = ', trim(chAlbedo01)
       print*, '<alb2>  = ', trim(chAlbedo02)
      endif
     else
      print*, '<indat> = ', trim(chInfile)
      print*, '<outdat> = ', trim(chOutfile)
     endif
    endif

!*******************************************************************************
!
!   Allocate and initialize arrays
!
!*******************************************************************************

!    For binary runs, set pixel number with respect to dataset dimension
!    For ASCII runs, set liRuns to 10,000 (= 10,000 datasets)
     if(bBinary) then
      liRuns = liCols*liRows
      liRecl = liCols*liRows
     else
      liRuns = 30000
      liRecl = 30000
     endif

!    Allocate arrays
     allocate( &
     prgiCloud(liRecl), &
     prgiMask(liRecl), &
     prgfNonAbs(liRecl), &
     prgfAbs(liRecl), &
     prgfNonAbsAlb(liRecl), &
     prgfAbsAlb(liRecl), &
     prgfSZen(liRecl), &
     prgfPZen(liRecl), &
     prgfRelAzm(liRecl), &
     prgfTau(liRecl), &
     prgfAef(liRecl), &
     prgfLWP(liRecl), &
     prgfSSA(liRecl), &
     prgfPAL(liRecl), &
     prgfCloudAlb(liRecl), &
     prgfCGT(liRecl), &
     prgfCDC(liRecl), &
     prgfDC(liRecl))

!    Initialize arrays
     if(bBinary) then
      prgiCloud     = -99.0
      prgiMask      = -99.0
     else
      prgiCloud     =   2.0
      prgiMask      =   1.0
     endif
     prgfSZen      = -99.0
     prgfPZen      = -99.0
     prgfRelAzm    = -99.0
     prgfNonAbs    = -99.0
     prgfAbs       = -99.0
     prgfTau       = -99.0
     prgfAef       = -99.0
     prgfSSA       = -99.0
     prgfLWP       = -99.0
     prgfPAL       = -99.0
     prgfCloudAlb  = -99.0
     prgfCGT       = -99.0
     prgfCDC       = -99.0
     prgfDC       = -99.0

    if(.not.bAlbedo) then
     prgfNonAbsAlb =  0.00d0
     prgfAbsAlb    =  0.00d0
    endif

!*******************************************************************************
!
!   Read satellite reflection data from idrisi file
!
!*******************************************************************************

!    Binary input
     if(bBinary) then
      print*, ' '
      print*, 'Reading input satellite datasets...'

      open(501,file=chSZen,access='direct',recl=liRecl)
      read(501,rec=1) prgfSZen
      close(501)
      prgfSZen = abs(prgfSZen)

      open(501,file=chPZen,access='direct',recl=liRecl)
      read(501,rec=1) prgfPZen
      close(501)
      prgfPZen = abs(prgfPZen)

      open(501,file=chRelAzm,access='direct',recl=liRecl)
      read(501,rec=1) prgfRelAzm
      close(501)

      open(501,file=chBand01,access='direct',recl=liRecl)
      read(501,rec=1) prgfNonAbs
      close(501)
      prgfNonAbs = prgfNonAbs / cosd(prgfSZen)

      open(501,file=chBand02,access='direct',recl=liRecl)
      read(501,rec=1) prgfAbs
      close(501)
      prgfAbs = prgfAbs / cosd(prgfSZen)

      if(bAlbedo) then
       open(501,file=chAlbedo01,access='direct',recl=liRecl)
       read(501,rec=1) prgfNonAbsAlb
       close(501)

       open(501,file=chAlbedo02,access='direct',recl=liRecl)
       read(501,rec=1) prgfAbsAlb
       close(501)
      endif

      open(501,file=chCloudMask,access='direct',recl=liRecl)
      read(501,rec=1) prgiCloud
      close(501)

      open(501,file=chComputationMask,access='direct',recl=liRecl)
      read(501,rec=1) prgiMask
      close(501)

!    ASCII input
     else
      print*, ' '
      print*, 'Reading input ASCII datasets...'
      open(501,file=chInfile)
      read(501,*)
      liRuns = 0
      do while (.not.EOF(501))
       liRuns = liRuns + 1
!      Only necessary for testing - start
       if(bTest) then
        read(501,*) prgfSZen(liRuns), prgfPZen(liRuns), prgfRelAzm(liRuns), prgfNonAbs(liRuns), prgfAbs(liRuns), prgfTauIn(liRuns), prgfAefIn(liRuns)
!      Only necessary for testing - end
       else
        read(501,*) prgfSZen(liRuns), prgfPZen(liRuns), prgfRelAzm(liRuns), prgfNonAbs(liRuns), prgfAbs(liRuns)
       endif
      enddo
      close(501)
     endif
     print*, 'Ready'
!*******************************************************************************
! 
!   Read LUTs
! 
!*******************************************************************************

!   Set data blocks
    rgfLUTRInfSSAGrid02 = rgfLUTRInfSSAGrid

!   Set pointers
    porgiLUTRInf => rgiLUTRInf
    pochLUT_RInf => chLUT_RInf
    pofgLUT => fgLUT
    porgfLUTRInfSSAGrid => rgfLUTRInfSSAGrid

    if(fWavelength01.gt.0.8) then
     podrgC => drgC0859
    else
     podrgC => drgC0645
    endif

    if(bprint) print*, '001 podrgC: ', podrgC

!   Set imaginary part of refractive index
    call SLALOM_RefractiveIndex

!   Read LUTs necessary for CLOUD code
!   LUT for non-absorbing channel wavelength
    pochLUT_RInf => chLUT_RInf
    porgiLUTRInf => rgiLUTRInf
    if(bprint) print*, '002 pochLUT_RInf (should be band 01): ', pochLUT_RInf
    call CLOUD_readLUT
    if(bprint) print*, '003 rgiLUTRInf  (should be band 01): ', &
                       rgiLUTRInf(10,10,10,10)

!   LUT for absorbing channel wavelength
    pochLUT_RInf => chLUT_RInf02
    porgiLUTRInf => rgiLUTRInf02
    if(bprint) print*, '004 pochLUT_RInf  (should be band 02): ', pochLUT_RInf
    call CLOUD_readLUT
    if(bprint) print*, '005 rgiLUTRInf02  (should be band 02): ', &
                       rgiLUTRInf02(10,10,10,10)

!   Set asymmetry parameter used for computation of actual LUT for RInf
!   Non-absorbing channel settings
    pochLUT_RInf => chLUT_RInf
    pofgLUT=> fgLUT
    if(pochLUT_RInf.eq.'lut_RInf_0645_aef_06.dat') then
     pofgLUT = 0.850018000000000000  ! 06µm, 645.5nm
    elseif(pochLUT_RInf.eq.'lut_RInf_0645_aef_10.dat') then
     pofgLUT = 0.861757000000000000  ! 10µm, 645.5nm
    elseif(pochLUT_RInf.eq.'lut_RInf_0645_aef_16.dat') then
     pofgLUT = 0.869221000000000000  ! 16µm, 645.5nm
    elseif(pochLUT_RInf.eq.'lut_RInf_0856_aef_06.dat') then
     pofgLUT = 0.843722000000000000  ! 06µm, 856.5nm
    elseif(pochLUT_RInf.eq.'lut_RInf_0856_aef_10.dat') then
     pofgLUT = 0.857995000000000000  ! 10µm, 856.5nm
    elseif(pochLUT_RInf.eq.'lut_RInf_0856_aef_16.dat') then
     pofgLUT = 0.867233000000000000  ! 16µm, 856.5nm
    elseif(pochLUT_RInf.eq.'lut_RInf_1630_aef_06.dat') then
     pofgLUT = 0.817538000000000000  ! 06µm, 1630nm
    elseif(pochLUT_RInf.eq.'lut_RInf_1630_aef_10.dat') then
     pofgLUT = 0.846061000000000000  ! 10µm, 1630nm
    elseif(pochLUT_RInf.eq.'lut_RInf_1630_aef_16.dat') then
     pofgLUT = 0.861364000000000000  ! 16µm, 1630nm
    elseif(pochLUT_RInf.eq.'lut_RInf_ice.dat') then
     pofgLUT = 0.752400000000000000  ! all aef
    else
     print*, 'No valid LUT for RInf has been selected.'
     print*, 'The program is going to stop...'
     stop
    endif
    if(bprint) print *, '006 fgLUT  (should be band 01): ', fgLUT

!   Absorbing channel settings
    pochLUT_RInf => chLUT_RInf02
    pofgLUT => fgLUT02
    if(pochLUT_RInf.eq.'lut_RInf_0645_aef_06.dat') then
     pofgLUT = 0.850018000000000000  ! 06µm, 645.5nm
    elseif(pochLUT_RInf.eq.'lut_RInf_0645_aef_10.dat') then
     pofgLUT = 0.861757000000000000  ! 10µm, 645.5nm
    elseif(pochLUT_RInf.eq.'lut_RInf_0645_aef_16.dat') then
     pofgLUT = 0.869221000000000000  ! 16µm, 645.5nm
    elseif(pochLUT_RInf.eq.'lut_RInf_0856_aef_06.dat') then
     pofgLUT = 0.843722000000000000  ! 06µm, 856.5nm
    elseif(pochLUT_RInf.eq.'lut_RInf_0856_aef_10.dat') then
     pofgLUT = 0.857995000000000000  ! 10µm, 856.5nm
    elseif(pochLUT_RInf.eq.'lut_RInf_0856_aef_16.dat') then
     pofgLUT = 0.867233000000000000  ! 16µm, 856.5nm
    elseif(pochLUT_RInf.eq.'lut_RInf_1630_aef_06.dat') then
     pofgLUT = 0.817538000000000000  ! 06µm, 1630nm
    elseif(pochLUT_RInf.eq.'lut_RInf_1630_aef_10.dat') then
     pofgLUT = 0.846061000000000000  ! 10µm, 1630nm
    elseif(pochLUT_RInf.eq.'lut_RInf_1630_aef_16.dat') then
     pofgLUT = 0.861364000000000000  ! 16µm, 1630nm
    elseif(pochLUT_RInf.eq.'lut_RInf_ice.dat') then
     pofgLUT = 0.752400000000000000  ! all aef
    else
     print*, 'No valid LUT for RInf has been selected.'
     print*, 'The program is going to stop...'
     stop
    endif
    if(bprint) print *, '007 fgLUT02  (should be band 02): ', fgLUT02

!   Correct LUT for RInf with respect to asymmetry parameter used for its
!   computation
!   Wavelength of non-absorbing channel
    porgfLUTRInfSSAGrid => rgfLUTRInfSSAGrid
    pofgLUT => fgLUT
    if(bprint) print*, '008 porgfLUTRInfSSAGrid  (should be band 01): ', &
                       porgfLUTRInfSSAGrid(5)
    if(bprint) print*, '009 pofgLUT  (should be band 01): ', pofgLUT
    porgfLUTRInfSSAGrid    = 1.0 - sqrt( (1.0-porgfLUTRInfSSAGrid) / &
                            (1.0-porgfLUTRInfSSAGrid*pofgLUT) )
    rgfLUTRInfSSAGrid = porgfLUTRInfSSAGrid
    if(bprint) print*, '010 rgfLUTRInfSSAGrid  (should be band 01): ', &
                       rgfLUTRInfSSAGrid(5)

!   Wavelength of absorbing channel
    porgfLUTRInfSSAGrid => rgfLUTRInfSSAGrid02
    pofgLUT => fgLUT02
    if(bprint) print*, '011 porgfLUTRInfSSAGrid  (should be band 02): ', &
                       porgfLUTRInfSSAGrid(5)
    if(bprint) print*, '012 pofgLUT  (should be band 02): ', pofgLUT
    porgfLUTRInfSSAGrid = 1.0 - sqrt( (1.0-porgfLUTRInfSSAGrid) / &
                         (1.0-porgfLUTRInfSSAGrid*pofgLUT) )
    rgfLUTRInfSSAGrid02 = porgfLUTRInfSSAGrid
    if(bprint) print*, '013 rgfLUTRInfSSAGrid02  (should be band 02): ', &
                       rgfLUTRInfSSAGrid02(5)

    rgfLUTKSSAGrid     = 1.0 - sqrt( (1.0-rgfLUTKSSAGrid) / &
                        (1.0-rgfLUTKSSAGrid*pofgLUT) )
    rgfLUTRdInfSSAGrid = 1.0 - sqrt( (1.0-rgfLUTRdInfSSAGrid) / &
                        (1.0-rgfLUTRdInfSSAGrid*pofgLUT) )

!*******************************************************************************
! 
!   Compute retrieval pixel by pixel
! 
!*******************************************************************************

!   Open log file
    open(501,file='SLALOM_output.dat')
    write(501,960)
    write(501,*) trim(chVersion), ' (CLOUD: 1.43)'

!   Compute pixel by pixel
    CALL TIME (STIME)
    print*, STIME

    do 10 liCounter = 1, liRuns
     if(mod(liCounter,25000).eq.0) print*, 'Calculating', liCounter, &
                                           ' Pixel of', liRuns
!    Set actual pixel values
     fNonAbs    = prgfNonAbs(liCounter)
     fAbs       = prgfAbs(liCounter)
     fSZen      = prgfSZen(liCounter)
     fPZen      = prgfPZen(liCounter)
     fRelAzm    = prgfRelAzm(liCounter)
     fNonAbsAlb = prgfNonAbsAlb(liCounter)
     fAbsAlb    = prgfAbsAlb(liCounter)
     iCloud     = int(prgiCloud(liCounter))

!    Only necessary for testing - start
     if(bTest) then
      fTauIn = prgfTauIn(liCounter)
      fAefIn = prgfAefIn(liCounter)
     endif
!    Only necessary for testing - end


!Uncomment the following lines for theoretical error studies
!Start change
!!!     do liErrorCounter = 799,1200
!!!      if(liErrorCounter.eq.799) then
!!!       fNonAbs = prgfNonAbs(liCounter)
!!!       fAbs    = prgfAbs(liCounter)
!!!      else
!!!       fNonAbs = prgfNonAbs(liCounter)*float(liErrorCounter)/1000.0
!!!       fAbs    = prgfAbs(liCounter)*float(liErrorCounter)/1000.0
!!!      endif
!!!      if(fNonAbs.gt.1.000) fNonAbs = 1.000
!!!      if(fAbs.gt.1.000) fAbs = 1.000
!End change

!    Check if actual pixel data is valid; if yes: compute
     if(fNonAbs.ne.0.000.and.iCloud.gt.1.and.prgiMask(liCounter).eq.1) then

!     Compute retrieval
      call SLALOM_Retrieve

!     Write results to output arrays
      prgfTau(liCounter) = fTau
      prgfAef(liCounter) = fAef
      prgfLWP(liCounter) = fLWP
      prgfSSA(liCounter) = dSSARetrieve
      prgfPAL(liCounter) = fPAL
      prgfCloudAlb(liCounter) = fCloudAlb
      prgfCGT(liCounter) = fCGT
      prgfCDC(liCounter) = fCDC
      prgfDC(liCounter) = fDC

!     Only necessary for testing - start
      if(bTest) then
       prgfTauError(liCounter) = fTau / prgfTauIn(liCounter) * 100.0 - 100.0
       prgfAefError(liCounter) = fAef / prgfAefIn(liCounter) * 100.0 - 100.0
      endif
!     Only necessary for testing - end

!Uncomment the following lines for theoretical error studies
!Start change
!!!      if(liErrorCounter.lt.800) then
!!!       fTauTarget = fTau
!!!       fAefTarget = fAef
!!!      else
!!!       prgfTauError(liCounter) = fTau / fTauTarget * 100.0 - 100.0
!!!       prgfAefError(liCounter) = fAef / fAefTarget * 100.0 - 100.0
!!!       print 999, prgfTauIn(liCounter), prgfAefIn(liCounter), float(liErrorCounter)/10.0, fNonAbs, fAbs, prgfTau(liCounter), prgfAef(liCounter), prgfTauError(liCounter), prgfAefError(liCounter)
!!!      endif
!End change

     endif

!Uncomment the following lines for theoretical error studies
!Start change
!!!      enddo
!!!      print *, " "
!!!      print *, " "
!End change


10  continue

    close(501)
    CALL TIME (STIME)
    print*, STIME

!   Write retrieval results to output binary files
    if(bBinary) then
     open(501,file=chTau,access='direct',recl=liRecl)
     write(501,rec=1) prgfTau
     close(501)

     open(501,file=chAef,access='direct',recl=liRecl)
     write(501,rec=1) prgfAef
     close(501)

     open(501,file=chLWP,access='direct',recl=liRecl)
     write(501,rec=1) prgfLWP
     close(501)

     open(501,file=chSSA,access='direct',recl=liRecl)
     write(501,rec=1) prgfSSA
     close(501)

     open(501,file=chPAL,access='direct',recl=liRecl)
     write(501,rec=1) prgfPAL
     close(501)

     open(501,file=chCloudAlb,access='direct',recl=liRecl)
     write(501,rec=1) prgfCloudAlb
     close(501)

     open(501,file=chCGT,access='direct',recl=liRecl)
     write(501,rec=1) prgfCGT
     close(501)

     open(501,file=chCDC,access='direct',recl=liRecl)
     write(501,rec=1) prgfCDC
     close(501)

     open(501,file=chDC,access='direct',recl=liRecl)
     write(501,rec=1) prgfDC
     close(501)

!   Write retrieval results to output ASCII files
    else
     open(501,file=chOutfile)
     do liCounter = 1, liRuns

!     Only necessary for testing - start
      if(bTest) then
       write (501,990) prgfSZen(liCounter), prgfPZen(liCounter), &
                       prgfRelAzm(liCounter), &
                       prgfTau(liCounter), prgfAef(liCounter), &
                       prgfTauIn(liCounter), prgfAefIn(liCounter), &
                       prgfTauError(liCounter), prgfAefError(liCounter)
       if(prgfTauIn(liCounter).ge.99.9) then
        write(501,*)
        write(501,*)
       endif
!      Only necessary for testing - end
      else
       write (501,*) prgfSZen(liCounter), prgfPZen(liCounter), &
                     prgfRelAzm(liCounter), &
                     prgfTau(liCounter), prgfAef(liCounter)
      endif

     enddo
     close(501)
    endif

    print*,'SLALOM finished.'

    end program SLALOM

!###############################################################################
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
!!  Please send any comments, suggestions, criticism, or (for our sake) bug
!!  reports to nauss@lcrs.de
!! 
!!  Development
!!  2003-2008 
!!  Thomas Nauss, Alexander A. Kokhanovsky
!! 
!!  This program is free software without any warranty and without even the
!!  implied warranty of merchantability or fitness for a particular purpose. 
!!  If you use the software you must acknowledge the software and its authors. 
!! 
!###############################################################################

    subroutine SLALOM_RefractiveIndex

!*******************************************************************************
!
!   Declaration of variables 
!
!*******************************************************************************

    use SLALOM_defs
    implicit none

    integer(2) :: iCounter      !! Counter
    integer(2) :: iCheck        !! Check flag for read statement

    real(8) :: dWavelengthCheckSmaler   !! Smaler wavelength in LUT file
    real(8) :: dRealPartSmaler          !! Smaler real refreactive index
    real(8) :: dImaginaryPartSmaler     !! Larger imaginary refreactive index
    real(8) :: dWavelengthCheckLarger   !! Larger wavelength in LUT file
    real(8) :: dRealPartLarger          !! Larger real refreactive index
    real(8) :: dImaginaryPartLarger     !! Larger imaginary refreactive index
    real(8) :: dRelativeWavelength      !! Relative wavelength difference
    real(8) :: dTemp                    !! Temporal variable

!*******************************************************************************
!
!  Open LUT file
!
!*******************************************************************************

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

!*******************************************************************************
!
!  Search and read imaginary part of refractive index with respect to wavelength
!
!*******************************************************************************

!   Look for the two bounding values of the refractive index
    read(501,*) dWavelengthCheckSmaler, dRealPartSmaler, dImaginaryPartSmaler
    do
     read(501,*,iostat=iCheck) dWavelengthCheckLarger, dRealPartLarger, &
                               dImaginaryPartLarger
     if(iCheck.lt.0) then
      print*, 'No appropriate refractive index has been found.'
      print*, 'Program is going to stop...'
      stop
     endif
     if(dWavelengthCheckSmaler.le.fWavelength02.and. &
        dWavelengthCheckLarger.ge.fWavelength02) then
      exit
     else
      dWavelengthCheckSmaler = dWavelengthCheckLarger
      dRealPartSmaler        = dRealPartLarger
      dImaginaryPartSmaler   = dImaginaryPartLarger
     endif
    enddo

!   Compute refractive index
    dRelativeWavelength = (fWavelength02-dWavelengthCheckSmaler) / &
                          (dWavelengthCheckLarger-dWavelengthCheckSmaler)
    call CLOUD_1DInterpolation(dRelativeWavelength,dImaginaryPartSmaler, &
                               dImaginaryPartLarger,dTemp)
    fMi = dTemp
    return
    end subroutine SLALOM_RefractiveIndex

!###############################################################################
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
!!  Please send any comments, suggestions, criticism, or (for our sake) bug
!!  reports to nauss@lcrs.de
!! 
!!  Development
!!  2003-2008 
!!  Thomas Nauss, Alexander A. Kokhanovsky
!! 
!!  This program is free software without any warranty and without even the
!!  implied warranty of merchantability or fitness for a particular purpose. 
!!  If you use the software you must acknowledge the software and its authors. 
!! 
!###############################################################################
 
    subroutine SLALOM_Retrieve
     
!*******************************************************************************
! 
!   Declaration of variables 
! 
!******************************************************************************* 

    use SLALOM_defs
    implicit none

    external functionAef  !! Function for the computation of aef

    integer(2) :: iCounter  !! Counter

    real(8) :: dSumC   !! Sum of C parameters (see SLALOM_Defs)
    real(8) :: dSumD   !! Sum of D parameters (see SLALOM_Defs)
    real(8) :: zbrent  !! Find root of functionR using Brent's method

!*******************************************************************************
! 
!  Retrieve effective cloud droplet radius.
! 
!*******************************************************************************

!   Computation for water clouds with iterration using aef
    if(bWaterCloud) then

!    Compute reflection of a semi-infinite cloud
     bFinite = .TRUE.
     fTau    = 5.0   ! Value does not matter.
     fSSA    = 1.0

     porgiLUTRInf => rgiLUTRInf
     pochLUT_RInf => chLUT_RInf
     pofgLUT => fgLUT
     porgfLUTRInfSSAGrid => rgfLUTRInfSSAGrid
     call cloud

     dR0 = dRinf_1000
     if(bprint) print*, '014 dR0: ', dR0
     bFinite = .TRUE.

!    Compute transmission for non-absorbing measurement
     ft1    = ( (1.0-fNonAbsAlb)*(dRinf_1000-fNonAbs) ) / &
              ( (1.0-fNonAbsAlb)*fEscape0_1000*fEscape_1000 - &
              fNonAbsAlb*(dRInf-fNonAbs) )

      if(bprint) print*, '015 fEscape0_1000, fEscape_1000: ', &
                         fEscape0_1000, fEscape_1000
      if(bprint) print*, '016 ft1: ', ft1

!    Retrieve effective radius using Brent's routine
     fAef = zbrent(functionAef,fAefMinValid,fAefMaxValid,1d-9)

!    Compute corresponding optical thickness
     dSumC = 0.d0
     do iCounter = 0,4
      dSumC = dSumC + podrgC(iCounter+1) * &
              ( 2.0*fPi/fWavelength01 * fAef )**(-2.0*float(iCounter)/3.0)
     enddo
     fg = 1.d0 - dSumC

!   Modified for next version 2008-07-16
!    print*, "1a", fg, fAef
!    if(fAef.ge.3.0.and.fAef.le.30.0) then
!     fg = rgfLUTgNonAbs(nint(fAef)-2)
!    elseif(fAef.lt.3.0) then
!     fg = rgfLUTgNonAbs(1)
!    else
!     fg = rgfLUTgNonAbs(28)
!    endif
!    print*, "1b", fg
!   Modified for next version 2008-07-16

     if(bprint) print*, '018 fg: ', fg

     fTau   = (ft1**(-1.0)-1.072) / (0.75*(1.0-fg))

     if(bprint) print*, '019 fTau: ', fTau

!    Compute LWP from aef for non-absorbing wavelength
     fLWP = fTau * fRho * fAef / &
            1.5 * ( 1.0 + 1.1 / &
            ( 2.0*fPi/fWavelength01 * fAef ) ** (2.0/3.0) ) ** (-1.0)

!    Compute cloud albedo
     fCloudAlb = 1 - (1.072 + 0.75 * (1.0 - fg) * fTau)**(-1)

!    Compute cloud geometrical thickness
     fCGT = fLWP / drgCWC(iCloudType)

!    Compute columnar droplet concentration
     fCDC = 81.0 * fTau * 1.0e+8 / (112.0 * fPi * fAef**2)

!    Compute droplet concentration
     fDC = fCDC / fCGT

!!!    print*, 'R_mes_0650:                   ', fNonAbs
!!!    print*, 'R_mes_1640:                   ', fAbs
!!!    print*, 'aef retrieved:                ', fAef
!!!    print*, 'tau retrieved:                ', fTau
!!!    print*, 'LWP retrieved:                ', fLWP
!!!    print*, 't_non_abs used within brent:  ', ft1
!!!    print*, 'last value of ssa in brent:   ', fSSA
!!!    print*, 'fg used to compute final tau: ', fg
!!!    pause
!!!    return

!   Computation of ice clouds with iterration usin SSA
!use the following logic:
!1)g=0.75.the exact number from Table 3 ..as in my e-mail
!2) tau is found immediately from transmittance in the visible for this "g"
!3)single scattering albedo is found using zbrent(func_ice,...)
!FUNCTION funk_ice
!func_ice= funk_1240nm_measured-  (  R_inf(SSA)- fun(tau,g,SSA) 
!SSA is found easily
!-- 
!I advice use 1240nm measurement and not 1640nm measurement for ice
!-- 
!see our JGR paper
!-- 
!after you find SSA - you calculate a_ef from our JGR-paper-formula(1 line) 
    else
!    Compute reflection of a semi-infinite cloud
     bFinite = .TRUE.
     fTau    = 5.0   ! Value does not matter.
     fSSA    = 1.0

     porgiLUTRInf => rgiLUTRInf
     pochLUT_RInf => chLUT_RInf
     pofgLUT => fgLUT
     porgfLUTRInfSSAGrid => rgfLUTRInfSSAGrid
     call cloud

     dR0 = dRinf_1000
     if(bprint) print*, '014 dR0: ', dR0
     bFinite = .TRUE.

!    Compute transmission for non-absorbing measurement
     ft1    = ( (1.0-fNonAbsAlb)*(dRinf_1000-fNonAbs) ) / &
              ( (1.0-fNonAbsAlb)*fEscape0_1000*fEscape_1000 - &
              fNonAbsAlb*(dRInf-fNonAbs) )

      if(bprint) print*, '015 fEscape0_1000, fEscape_1000: ', &
                         fEscape0_1000, fEscape_1000
      if(bprint) print*, '016 ft1: ', ft1

!    Retrieve effective radius using Brent's routine
!     fAef = zbrent(functionAef,fAefMinValid,fAefMaxValid,1d-9)
     fAef = 3.0
     fg = 0.752400000000000000

!   Modified for next version 2008-07-16
!    print*, "1a", fg, fAef
!    if(fAef.ge.3.0.and.fAef.le.30.0) then
!     fg = rgfLUTgNonAbs(nint(fAef)-2)
!    elseif(fAef.lt.3.0) then
!     fg = rgfLUTgNonAbs(1)
!    else
!     fg = rgfLUTgNonAbs(28)
!    endif
!    print*, "1b", fg
!   Modified for next version 2008-07-16

     if(bprint) print*, '018 fg: ', fg

     fTau   = (ft1**(-1.0)-1.072) / (0.75*(1.0-fg))

!    Test Nauss Start
999  format(5f10.5)
     write(*,999) fNonAbs, dRinf_1000, ft1, fg, fTau
!    Test Nauss End

     if(bprint) print*, '019 fTau: ', fTau

!    Compute LWP from aef for non-absorbing wavelength
     fLWP = fTau * fRho * fAef / &
            1.5 * ( 1.0 + 1.1 / &
            ( 2.0*fPi/fWavelength01 * fAef ) ** (2.0/3.0) ) ** (-1.0)

!    Compute cloud albedo
     fCloudAlb = 1 - (1.072 + 0.75 * (1.0 - fg) * fTau)**(-1)

!    Compute cloud geometrical thickness
     fCGT = fLWP / drgCWC(iCloudType)

!    Compute columnar droplet concentration
     fCDC = 81.0 * fTau * 1.0e+8 / (112.0 * fPi * fAef**2)

!    Compute droplet concentration
     fDC = fCDC / fCGT

!!!    print*, 'R_mes_0650:                   ', fNonAbs
!!!    print*, 'R_mes_1640:                   ', fAbs
!!!    print*, 'aef retrieved:                ', fAef
!!!    print*, 'tau retrieved:                ', fTau
!!!    print*, 'LWP retrieved:                ', fLWP
!!!    print*, 't_non_abs used within brent:  ', ft1
!!!    print*, 'last value of ssa in brent:   ', fSSA
!!!    print*, 'fg used to compute final tau: ', fg
!!!    pause
!!!    return
    endif
    end subroutine SLALOM_Retrieve

!###############################################################################
!! 
!!  NAME 
!!  Function "functionAef" 
!! 
!!  VERSION 
!!  $Revision: 1.7 $ 
!! 
!!  PURPOSE 
!!  Compute difference between meassured and computed reflection of the 
!!  absorbing channel in order to retrieve the associated effective droplet
!!  radius (called from Brent's routine).
!!
!!  CONTACT 
!!  Please send any comments, suggestions, criticism, or (for our sake) bug
!!  reports to nauss@lcrs.de
!! 
!!  Development
!!  2003-2008 
!!  Thomas Nauss, Alexander A. Kokhanovsky
!! 
!!  This program is free software without any warranty and without even the
!!  implied warranty of merchantability or fitness for a particular purpose. 
!!  If you use the software you must acknowledge the software and its authors. 
!! 
!###############################################################################
 
    function functionAef(dAef)
 
!*******************************************************************************
! 
!   Declaration of variables 
! 
!*******************************************************************************

    use SLALOM_defs
    implicit none

    integer(2) :: iCounter  !! Counter

    real(8) :: dAef         !! Actual aef within Brent's routine (functionBrent)
    real(8) :: functionAef  !! Rabs-Rcomp
    real(8) :: dSumC        !! Sum of C parameters (see SLALOM_Defs)
    real(8) :: dSumD        !! Sum of D parameters (see SLALOM_Defs)

!******************************************************************************* 
! 
!  Compute 
! 
!*******************************************************************************

!   Compute asymetry parameter for non-absorbing wavelength
    dSumC = 0.d0
    do iCounter = 0,4
     dSumC = dSumC + podrgC(iCounter+1) * &
             ( 2.0*fPi/fWavelength01 * dAef )**(-2.0*float(iCounter)/3.0)
    enddo
    if(bprint) print*, '100 podrgC: ', podrgC
    if(bprint) print*, '100 fWavelength01: ', fWavelength01
    if(bprint) print*, '100 dAef: ', dAef
    if(bprint) print*, '100 dSumC: ', dSumC
 
    fg = 1.d0 - dSumC

    if(bprint) print*, '101 fg: ', fg

!   Modified for next version 2008-07-16
!    print*, "2a", fg, fAef
!    if(dAef.ge.3.0.and.dAef.le.30.0) then
!     fg = rgfLUTgNonAbs(nint(dAef)-2)
!    elseif(dAef.lt.3.0) then
!     fg = rgfLUTgNonAbs(1)
!    else
!     fg = rgfLUTgNonAbs(28)
!    endif
!    print*, "2b", fg
!   Modified for next version 2008-07-16

!   Compute tau for non-absorbing wavelength
    fTau   = (ft1**(-1.0)-1.072) / (0.75*(1.0-fg))

    if(bprint) print*, '102 fTau: ', fTau

!   Compute LWP for non-absorbing wavelength
    fLWP = fTau * fRho * dAef / &
           1.5 * ( 1.0 + 1.1 / &
           ( 2.0*fPi/fWavelength01 * dAef ) ** (2.0/3.0) ) ** (-1.0)

    if(bprint) print*, '103 fLWP: ', fLWP


!   Compute tau for absorbing wavelength using LWP from non-absorbing wavelength
    fTau = 1.5 * fLWP / ( fRho * dAef) * &
           ( 1.0 + 1.1 / ( 2.0*fPi/fWavelength02 * dAef ) ** (2.0/3.0) )

    if(bprint) print*, '104 fTau: ', fTau

!   Compute the extinction coefficient for absorbing wavelength
    fSigma_ext = 1.5 / dAef * &
                 (1.0 + 1.1 / ( 2.0*fPi/fWavelength02 * dAef )**(2.0/3.0) )

    if(bprint) print*, '105 fSigma_ext: ', fSigma_ext

!   Compute the absorbtion coefficient for absorbing wavelength
    dSumD = 0.d0
    do iCounter = 0,4
     dSumD = dSumD + drgD(iCounter+1) * &
             ( 2.0*fPi/fWavelength02 * dAef )**float(iCounter)
    enddo

    if(bprint) print*, '106 dSumD: ', dSumD

    fSigma_abs = 4.0 * fPi * fMi / fWavelength02 * dSumD

    if(bprint) print*, '107 fSigma_abs: ', fSigma_abs

!   Compute ssa for absorbing wavelength based on actual aef value.
    fSSA = 1.d0 - fSigma_abs / fSigma_ext

    if(bprint) print*, '108 fSSA: ', fSSA

!   Compute asymetry parameter for absorbing wavelength based on actual aef.
    dSumC = 0.d0
    do iCounter = 0,4
     dSumC = dSumC + drgC1630(iCounter+1) * &
             ( 2.0*fPi/fWavelength02 * dAef )**(-2.0*float(iCounter)/3.0)
    enddo

    if(bprint) print*, '109 dSumC: ', dSumC

    fg = 1.d0 - dSumC

    if(bprint) print*, '110 fg: ', fg

!   Modified for next version 2008-07-16
!    print*, "2a", fg, fAef
!    if(dAef.ge.3.0.and.dAef.le.30.0) then
!     fg = rgfLUTgAbs(nint(dAef)-2)
!    elseif(dAef.lt.3.0) then
!     fg = rgfLUTgAbs(1)
!    else
!     fg = rgfLUTgAbs(28)
!    endif
!    print*, "2b", fg
!   Modified for next version 2008-07-16

!   Set pointers in order to compute reflection of the absorbing channel with
!   respect to the actual aef using CLOUD
    porgiLUTRInf => rgiLUTRInf02
    pochLUT_RInf => chLUT_RInf02
    pofgLUT => fgLUT02
    porgfLUTRInfSSAGrid => rgfLUTRInfSSAGrid02

    if(bprint) print*, '111 porgiLUTRInf: ', porgiLUTRInf(10,10,10,10)
    if(bprint) print*, '112 pochLUT_RInf: ', pochLUT_RInf
    if(bprint) print*, '113 pofgLUT: ', pofgLUT
    if(bprint) print*, '114 porgfLUTRInfSSAGrid: ', porgfLUTRInfSSAGrid(5)

    call cloud

!   Compute difference
    functionAef = fAbs - fRefl
    if(bprint) print*, '017 fAbs, fRefl: ', fAbs, fRefl

    end function functionAef

!###############################################################################
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
!!  Copyright (c) 1986-1992
!!  Numerical Recipes Software ,#23
!! 
!###############################################################################
 
    function zbrent(func,x1,x2,tol)
 
!*******************************************************************************
! 
!   Declaration of variables 
! 
!*******************************************************************************

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

!*******************************************************************************
!*******************************************************************************
! 
!   Include program cloud here
!   (and adjust main routine so it can be used as subroutine)
!
!******************************************************************************* 
!*******************************************************************************


!###############################################################################
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
!!  ADDITIONAL REQUIREMENTS 
!!  Configuration file "cloud.cfg" for definition of initial parameters.
!!
!!
!!  ATTENTION
!!  If used within SLALOM, several code clocks within the main routine
!!  (program CLOUD) have to be commented out. Look for "SLALOM" to see what to
!!  do. In addition, "use SLALOM_Defs" has to be replaced by "use SLALOM_Defs".
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
!!  $Id: slalom.f90,v 1.7 2008/08/26 14:21:51 tnauss Exp $ 
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
!!  $Log: slalom.f90,v $
!!  Revision 1.7  2008/08/26 14:21:51  tnauss
!!  Update
!!  Include cloud geometrical thickness and droplet concentration retrieval.
!!
!!  Revision 1.6  2008/08/26 07:06:35  tnauss
!!  Update
!!  Include test function in SLALOM code (i. e. known tau and aef is read from
!!  the ASCII input file in addition to the measurement data and error between
!!  tau and aef input and output is computed in addition).
!!
!!  Revision 1.5  2008/05/20 10:47:06  tnauss
!!  Update
!!  Include ASCII file in- and output.
!!
!!  Revision 1.43  2008/05/20 10:35:55  tnauss
!!  Update
!!  Delete some variables no longer necessary.
!!
!!  Revision 1.43  2008/05/19 15:12:18  tnauss
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

!!!    program CLOUD

!############# UNCOMMENT THE NEXT LINE IF USED WITHIN SLALOM ###################
    subroutine CLOUD

!*******************************************************************************
!
!   Declaration of variables 
!
!*******************************************************************************

    use SLALOM_Defs

    implicit none

!############# COMMENT THE FOLLOWING SECTION IF USED WITHIN SLALOM #############

!*******************************************************************************
!
!   Format
!
!*******************************************************************************

!!!950 format(a8,f15.7)    !! Print final results
!!!960 format(' szen    pzen  relazm  ssa     g       ' &
!!!           'tau    Refl   Trans  RDiff  TDiff  rs     t       ' &
!!!           'dRinf  dRDInf fEsc0  fEsc')
!!!961 format(3f7.2,f11.7,f8.5,f7.2,6f7.4,4f7.4)
!!!970 format(' szen    pzen  relazm  ssa     g       ' &
!!!           'tau    Refl   RDiff  rs     Trans  TDiff  t')
!!!971 format(3f7.2,f11.7,f8.5,f7.2,6f7.4)

!*******************************************************************************
!
!   Say hello...
!
!*******************************************************************************

!!!    chVersion = '$Revision: 1.7 $'
!!!    print*, ' '
!!!    print*, ' '
!!!    print*, 'Program CLOUD'
!!!    print*, trim(chVersion)

!*******************************************************************************
!
!   Read settings
!
!*******************************************************************************

!!!    chControlFile    = 'cloud.cfg'
!!!    call CLOUD_ReadSettings

!*******************************************************************************
!
!   Set pointers
!
!*******************************************************************************

!!!    pochLUT_RInf => chLUT_RInf
!!!    porgiLUTRInf => rgiLUTRInf
!!!    pofgLUT => fgLUT
!!!    porgfLUTRInfSSAGrid => rgfLUTRInfSSAGrid

!*******************************************************************************
!
!   Define parameters
!
!*******************************************************************************

!!!    chControlFile    = 'cloud.cfg'
!!!    call CLOUD_ReadSettings

!   Set asymmetry parameter used for computation of actual LUT for RInf
!!!   if(pochLUT_RInf.eq.'lut_RInf_0645_aef_06.dat') then
!!!     pofgLUT = 0.850018000000000000  ! 06µm, 645.5nm
!!!     rgchOutFile(1) = 'solar_lut_0645_06.dat'
!!!    elseif(pochLUT_RInf.eq.'lut_RInf_0645_aef_10.dat') then
!!!     pofgLUT = 0.861757000000000000  ! 10µm, 645.5nm
!!!     rgchOutFile(1) = 'solar_lut_0645_10.dat'
!!!    elseif(pochLUT_RInf.eq.'lut_RInf_0645_aef_16.dat') then
!!!     pofgLUT = 0.869221000000000000  ! 16µm, 645.5nm
!!!     rgchOutFile(1) = 'solar_lut_0645_16.dat'
!!!    elseif(pochLUT_RInf.eq.'lut_RInf_0856_aef_06.dat') then
!!!     pofgLUT = 0.843722000000000000  ! 06µm, 856.5nm
!!!     rgchOutFile(1) = 'solar_lut_0856_06.dat'
!!!    elseif(pochLUT_RInf.eq.'lut_RInf_0856_aef_10.dat') then
!!!     pofgLUT = 0.857995000000000000  ! 10µm, 856.5nm
!!!     rgchOutFile(1) = 'solar_lut_0856_10.dat'
!!!    elseif(pochLUT_RInf.eq.'lut_RInf_0856_aef_16.dat') then
!!!     pofgLUT = 0.867233000000000000  ! 16µm, 856.5nm
!!!     rgchOutFile(1) = 'solar_lut_0856_16.dat'
!!!    elseif(pochLUT_RInf.eq.'lut_RInf_1630_aef_06.dat') then
!!!     pofgLUT = 0.817538000000000000  ! 06µm, 1630nm
!!!     rgchOutFile(1) = 'solar_lut_1630_06.dat'
!!!    elseif(pochLUT_RInf.eq.'lut_RInf_1630_aef_10.dat') then
!!!     pofgLUT = 0.846061000000000000  ! 10µm, 1630nm
!!!     rgchOutFile(1) = 'solar_lut_1630_10.dat'
!!!    elseif(pochLUT_RInf.eq.'lut_RInf_1630_aef_16.dat') then
!!!     pofgLUT = 0.861364000000000000  ! 16µm, 1630nm
!!!     rgchOutFile(1) = 'solar_lut_1630_16.dat'
!!!    else
!!!     print*, 'No valid LUT for RInf has been selected.'
!!!     print*, 'The program is going to stop...'
!!!     stop
!!!    endif


!   Open output filenames    
!!!    open(550,file=trim(rgchOutFile(1)))
!!!    write(550,*) 'Version: ',trim(chVersion)
!!!    write(550,970)

!*******************************************************************************
!
!   Read LUTs
!
!*******************************************************************************

!!!    call CLOUD_ReadLUT

!*******************************************************************************
!
!   Compute radiative characteristics
!
!*******************************************************************************

!   Program controll settings
!!!    liTauMin = nint(fTauMin*100)
!!!    liTauMax = nint(fTauMax*100)
!!!    liTauStep = nint(fTauStep*100)
!!!    liSSAMin = nint(fSSAMin * 1000000)
!!!    liSSAMax = nint(fSSAMax * 1000000)
!!!    liSSAStep = nint(fSSAStep * 1000000)
!!!    liSZenMin = nint(fSZenMin*100)
!!!    liSZenMax = nint(fSZenMax*100)
!!!    liSZenStep = nint(fSZenStep*100)
!!!    liPZenMin = nint(fPZenMin*100)
!!!    liPZenMax = nint(fPZenMax*100)
!!!    liPZenStep = nint(fPZenStep*100)
!!!    liRelAzmMin = nint(fRelAzmMin*100)
!!!    liRelAzmMax = nint(fRelAzmMax*100)
!!!    liRelAzmStep = nint(fRelAzmStep*100)

!   Convert ssa grid values of the LUTs to similarity parameter (1-fs is used
!   for historical reasons since then, the coloumn with the smalest/largest ssa
!   value is the same).
!   LUTs have been computed for an asymmetry parameter of 0.8500, this value is
!   used for the conversion.

!!!    porgfLUTRInfSSAGrid  = 1.0 - sqrt( (1.0-porgfLUTRInfSSAGrid)/ &
!!!                           (1.0-porgfLUTRInfSSAGrid*pofgLUT) )
!!!    rgfLUTKSSAGrid     = 1.0 - sqrt( (1.0-rgfLUTKSSAGrid)/ &
!!!                         (1.0-rgfLUTKSSAGrid*pofgLUT) )
!!!    rgfLUTRdInfSSAGrid = 1.0 - sqrt( (1.0-rgfLUTRdInfSSAGrid)/ &
!!!                         (1.0-rgfLUTRdInfSSAGrid*pofgLUT) )

!!!    do 10 liTauCounter = liTauMin, liTauMax, liTauStep
!!!     do 20 liSZenCounter = liSZenMin, liSZenMax, liSZenStep
!!!      do 30 liPZenCounter = liPZenMin, liPZenMax, liPZenStep
!!!       do 40 liRelAzmCounter = liRelAzmMin, liRelAzmMax, liRelAzmStep
!!!        do 50 liSSACounter = liSSAMin, liSSAMax, liSSAStep

!!!        fTau = float(liTauCounter)/100.0
!!!        fSSA = float(liSSACounter)/1000000.0
!!!        fSZen = float(liSZenCounter)/100.0
!!!        fPZen = float(liPZenCounter)/100.0 
!!!        fRelAzm = float(liRelAzmCounter)/100.0

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

!!!     if(.not.bFinite  ) then
!!!      write(550,971) fSZen, fPZen, fRelAzm, fSSA, fg, &
!!!                     fTau, fRefl, fRDiff, frs, fTrans, fTdiff, ft
!!!     else
!!!      write(550,961) fSZen, fPZen, fRelAzm, fSSA, fg, &
!!!                     fTau, fRefl, fTrans, fRDiff, fTdiff, frs,ft, &
!!!                     dRinf, dRDInf, fEscape0, fEscape
!!!     endif


!!!     if(bPrint) then
!!!         print*,' '
!!!         print*,'Final results from program CLOUD:'
!!!         print*,'tau:           ', fTau
!!!         print*,'rs:            ', frs
!!!         print*,'t:             ', ft
!!!         print*,'RDiff:         ', fRDiff
!!!         print*,'TDiff:         ', fTDiff
!!!         print*,'ADiff:         ', fADiff
!!!         print*,'Trans:         ', fTrans
!!!         print*,'Refl:          ', fRefl
!!!         print*,'Ssa:           ', fSSA
!!!         print*,'fg:            ', fg
!!!         print*,'fSZen:         ', fSZen
!!!         print*,'fPZen:         ', fPZen
!!!         print*,'fRelAzm:       ', fRelAzm
!!!         if(bPause) pause
!!!        endif

!!!50      continue
!!!40     continue
!!!30    continue
!!!20   continue
!!!10  continue

!!!    close(550)
!!!    close(551)

!!!    stop
!!!    end program CLOUD

!############# END OF SECTION TO COMMENT IF USED WITHIN SLALOM #################
!############# UNCOMMENT THE NEXT TWO LINES IF USED WITHIN SLALOM ##############
    return
    end subroutine CLOUD

!###############################################################################
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

    use SLALOM_Defs
    implicit none

!*******************************************************************************
! 
!   Define namelist
! 
!*******************************************************************************

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
!!  $Revision: 1.7 $ 
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
 
    use SLALOM_Defs
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
    open(501,file='lut_k.dat') 
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
    open(501,file=trim(pochLUT_RInf))
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
    open(501,file='lut_RdInf.dat') 
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
!!  $Revision: 1.7 $ 
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
 
    use SLALOM_Defs
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
     if(bPause) pause
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
 
    use SLALOM_Defs
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
     if(bPause) pause
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
 
    use SLALOM_Defs
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
     if(bPause) pause
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


