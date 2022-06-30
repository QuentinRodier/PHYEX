!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
MODULE MODE_TURB_VER_THERMO_FLUX
IMPLICIT NONE
CONTAINS
SUBROUTINE TURB_VER_THERMO_FLUX(D,CST,CSTURB,TURBN,                 &
                      KRR,KRRL,KRRI,KSV,                            &
                      OTURB_FLX,HTURBDIM,HTOM,OOCEAN,ODEEPOC,OHARAT,&
                      OCOUPLES,OLES_CALL, OCOMPUTE_SRC,             &
                      PIMPL,PEXPL,PTSTEP,HPROGRAM,                  &
                      TPFILE,                                       &
                      PDXX,PDYY,PDZZ,PDZX,PDZY,PDIRCOSZW,PZZ,       &
                      PRHODJ,PTHVREF,                               &
                      PSFTHM,PSFRM,PSFTHP,PSFRP,                    &
                      PWM,PTHLM,PRM,PSVM,                           &
                      PTKEM,PLM,PLEPS,                              &
                      PLOCPEXNM,PATHETA,PAMOIST,PSRCM,PFRAC_ICE,    &
                      PBETA, PSQRT_TKE, PDTH_DZ, PDR_DZ, PRED2TH3,  &
                      PRED2R3, PRED2THR3, PBLL_O_E, PETHETA,        &
                      PEMOIST, PREDTH1, PREDR1, PPHI3, PPSI3, PD,   &
                      PFWTH,PFWR,PFTH2,PFR2,PFTHR,MFMOIST,PBL_DEPTH,&
                      PWTHV,PRTHLS,PRRS,PTHLP,PRP,PTP,PWTH,PWRC     )
!     ###############################################################
!
!
!!****  *TURB_VER_THERMO_FLUX* -compute the source terms due to the vertical turbulent
!!       fluxes.
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to compute the vertical turbulent
!     fluxes of the evolutive variables and give back the source 
!     terms to the main program.	In the case of large horizontal meshes,
!     the divergence of these vertical turbulent fluxes represent the whole
!     effect of the turbulence but when the three-dimensionnal version of
!     the turbulence scheme is activated (CTURBDIM="3DIM"), these divergences
!     are completed in the next routine TURB_HOR. 
!		  An arbitrary degree of implicitness has been implemented for the 
!     temporal treatment of these diffusion terms.
!       The vertical boundary conditions are as follows:
!           *  at the bottom, the surface fluxes are prescribed at the same
!              as the other turbulent fluxes
!           *  at the top, the turbulent fluxes are set to 0.
!       It should be noted that the condensation has been implicitely included
!     in this turbulence scheme by using conservative variables and computing
!     the subgrid variance of a statistical variable s indicating the presence 
!     or not of condensation in a given mesh. 
!
!!**  METHOD
!!    ------
!!      1D type calculations are made;
!!      The vertical turbulent fluxes are computed in an off-centered
!!      implicit scheme (a Crank-Nicholson type with coefficients different
!!      than 0.5), which allows to vary the degree of implicitness of the
!!      formulation.
!!      	 The different prognostic variables are treated one by one. 
!!      The contributions of each turbulent fluxes are cumulated into the 
!!      tendency  PRvarS, and into the dynamic and thermal production of 
!!      TKE if necessary.
!!        
!!			 In section 2 and 3, the thermodynamical fields are considered.
!!      Only the turbulent fluxes of the conservative variables
!!      (Thetal and Rnp stored in PRx(:,:,:,1))  are computed. 
!!       Note that the turbulent fluxes at the vertical 
!!      boundaries are given either by the soil scheme for the surface one
!!      ( at the same instant as the others fluxes) and equal to 0 at the 
!!      top of the model. The thermal production is computed by vertically 
!!      averaging the turbulent flux and multiply this flux at the mass point by
!!      a function ETHETA or EMOIST, which preform the transformation from the
!!      conservative variables to the virtual potential temperature. 
!!     
!! 	    In section 4, the variance of the statistical variable
!!      s indicating presence or not of condensation, is determined in function 
!!      of the turbulent moments of the conservative variables and its
!!      squarred root is stored in PSIGS. This information will be completed in 
!!      the horizontal turbulence if the turbulence dimensionality is not 
!!      equal to "1DIM".
!!
!!			 In section 5, the x component of the stress tensor is computed.
!!      The surface flux <u'w'> is computed from the value of the surface
!!      fluxes computed in axes linked to the orography ( i", j" , k"):
!!        i" is parallel to the surface and in the direction of the maximum
!!           slope
!!        j" is also parallel to the surface and in the normal direction of
!!           the maximum slope
!!        k" is the normal to the surface
!!      In order to prevent numerical instability, the implicit scheme has 
!!      been extended to the surface flux regarding to its dependence in 
!!      function of U. The dependence in function of the other components 
!!      introduced by the different rotations is only explicit.
!!      The turbulent fluxes are used to compute the dynamic production of 
!!      TKE. For the last TKE level ( located at PDZZ(:,:,IKB)/2 from the
!!      ground), an harmonic extrapolation from the dynamic production at 
!!      PDZZ(:,:,IKB) is used to avoid an evaluation of the gradient of U
!!      in the surface layer.
!!
!!         In section 6, the same steps are repeated but for the y direction
!!		  and in section 7, a diagnostic computation of the W variance is 
!!      performed.
!!
!!         In section 8, the turbulent fluxes for the scalar variables are 
!!      computed by the same way as the conservative thermodynamical variables
!!
!!            
!!    EXTERNAL
!!    --------
!!      GX_U_M, GY_V_M, GZ_W_M :  cartesian gradient operators 
!!      GX_U_UW,GY_V_VW	         (X,Y,Z) represent the direction of the gradient
!!                               _(M,U,...)_ represent the localization of the 
!!                               field to be derivated
!!                               _(M,UW,...) represent the localization of the 
!!                               field	derivated
!!                               
!!
!!      MXM,MXF,MYM,MYF,MZM,MZF
!!                             :  Shuman functions (mean operators)     
!!      DXF,DYF,DZF,DZM
!!                             :  Shuman functions (difference operators)     
!!                               
!!      SUBROUTINE TRIDIAG     : to compute the split implicit evolution
!!                               of a variable located at a mass point
!!
!!      SUBROUTINE TRIDIAG_WIND: to compute the split implicit evolution
!!                               of a variable located at a wind point
!!
!!      FUNCTIONs ETHETA and EMOIST  :  
!!            allows to compute:
!!            - the coefficients for the turbulent correlation between
!!            any variable and the virtual potential temperature, of its 
!!            correlations with the conservative potential temperature and 
!!            the humidity conservative variable:
!!            -------              -------              -------
!!            A' Thv'  =  ETHETA   A' Thl'  +  EMOIST   A' Rnp'  
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST : contains physical constants
!!
!!           CST%XG         : gravity constant
!!
!!      Module MODD_CTURB: contains the set of constants for
!!                        the turbulence scheme
!!
!!           CSTURB%XCMFS,XCMFB : cts for the momentum flux
!!           CSTURB%XCSHF       : ct for the sensible heat flux
!!           CSTURB%XCHF        : ct for the moisture flux
!!           CSTURB%XCTV,CSTURB%XCHV   : cts for the T and moisture variances
!!
!!      Module MODD_PARAMETERS
!!
!!           JPVEXT_TURB     : number of vertical external points
!!           JPHEXT     : number of horizontal external points
!!
!!
!!    REFERENCE
!!    ---------
!!      Book 1 of documentation (Chapter: Turbulence)
!!
!!    AUTHOR
!!    ------
!!      Joan Cuxart             * INM and Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       August   19, 1994
!!      Modifications: February 14, 1995 (J.Cuxart and J.Stein) 
!!                                  Doctorization and Optimization
!!      Modifications: March 21, 1995 (J.M. Carriere) 
!!                                  Introduction of cloud water
!!      Modifications: June  14, 1995 (J.Cuxart and J. Stein) 
!!                                 Phi3 and Psi3 at w-point + bug in the all
!!                                 or nothing condens. 
!!      Modifications: Sept  15, 1995 (J.Cuxart and J. Stein) 
!!                                 Change the DP computation at the ground
!!      Modifications: October 10, 1995 (J.Cuxart and J. Stein) 
!!                                 Psi for scal var and LES tools
!!      Modifications: November 10, 1995 (J. Stein)
!!                                 change the surface	relations 
!!      Modifications: February 20, 1995 (J. Stein) optimization
!!      Modifications: May 21, 1996 (J. Stein) 
!!                                  bug in the vertical flux of the V wind 
!!                                  component for explicit computation
!!      Modifications: May 21, 1996 (N. wood) 
!!                                  modify the computation of the vertical
!!                                   part or the surface tangential flux
!!      Modifications: May 21, 1996 (P. Jabouille)
!!                                  same modification in the Y direction
!!      
!!      Modifications: Sept 17, 1996 (J. Stein) change the moist case by using
!!                                  Pi instead of Piref + use Atheta and Amoist
!!
!!      Modifications: Nov  24, 1997 (V. Masson) removes the DO loops 
!!      Modifications: Mar  31, 1998 (V. Masson) splits the routine TURB_VER_THERMO_FLUX 
!!      Modifications: Oct  18, 2000 (V. Masson) LES computations
!!      Modifications: Dec  01, 2000 (V. Masson) conservation of energy from
!!                                               surface flux in 1DIM case
!!                                               when slopes are present
!!                     Nov  06, 2002 (V. Masson) LES budgets
!!                     Feb  20, 2003 (JP Pinty)  Add PFRAC_ICE
!!                     May  20, 2003 (JP Pinty)  Correction of ETHETA
!!                                                         and EMOIST calls
!!                     July     2005 (S. Tomas, V. Masson)
!!                                               Add 3rd order moments
!!                                               and implicitation of PHI3 and PSI3
!!                     October 2009 (G. Tanguy) add ILENCH=LEN(YCOMMENT) after
!!                                              change of YCOMMENT
!!                     2012-02 (Y. Seity) add possibility to run with reversed
!!                                             vertical levels
!!      Modifications  July 2015 (Wim de Rooy) OHARAT switch
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!                     2021 (D. Ricard) last version of HGRAD turbulence scheme
!!                                 Leronard terms instead of Reynolds terms
!!                                 applied to vertical fluxes of r_np and Thl
!!                                 for implicit version of turbulence scheme
!!                                 corrections and cleaning
!!                     June 2020 (B. Vie) Patch preventing negative rc and ri in 2.3 and 3.3
!! JL Redelsperger  : 03/2021: Ocean and Autocoupling O-A LES Cases
!!                             Sfc flux shape for LDEEPOC Case
!!--------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!
USE MODD_CST, ONLY: CST_t
USE MODD_CTURB, ONLY: CSTURB_t
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
USE MODD_FIELD,          ONLY: TFIELDDATA, TYPEREAL
USE MODD_GRID_n,         ONLY: XZS, XXHAT, XYHAT
USE MODD_IO,             ONLY: TFILEDATA
USE MODD_METRICS_n,      ONLY: XDXX, XDYY, XDZX, XDZY, XDZZ
USE MODD_PARAMETERS, ONLY: JPVEXT_TURB, JPHEXT
USE MODD_TURB_n,         ONLY: TURB_t
USE MODD_LES
USE MODD_DIM_n, ONLY: NIMAX_ll, NJMAX_ll
USE MODD_OCEANH, ONLY: XSSTFL
USE MODD_TURB_n, ONLY: TURB_t
USE MODD_FRC, ONLY: XCENTX_OC, XCENTY_OC, XRADX_OC,XRADY_OC
!
USE MODI_GRADIENT_U
USE MODI_GRADIENT_V
USE MODI_GRADIENT_W
USE MODI_GRADIENT_M
USE MODI_SHUMAN , ONLY : DZF, DZM, MZF, MZM, MYF, MXF
USE MODI_LES_MEAN_SUBGRID
USE MODE_TRIDIAG_THERMO, ONLY: TRIDIAG_THERMO
USE MODE_TM06_H, ONLY: TM06_H
!
USE MODE_IO_FIELD_WRITE, ONLY: IO_FIELD_WRITE
USE MODE_PRANDTL
USE SHUMAN_PHY, ONLY: MZM_PHY, MZF_PHY, DZM_PHY, DZF_PHY
!
USE MODI_SECOND_MNH
USE MODE_ll
USE MODE_GATHER_ll
!
IMPLICIT NONE
!
!*      0.1  declarations of arguments
!
!
!
TYPE(DIMPHYEX_t),       INTENT(IN)   :: D
TYPE(CST_t),            INTENT(IN)   :: CST
TYPE(CSTURB_t),         INTENT(IN)   :: CSTURB
TYPE(TURB_t),           INTENT(IN)   :: TURBN
INTEGER,                INTENT(IN)   :: KRR           ! number of moist var.
INTEGER,                INTENT(IN)   :: KSV           ! number of scalar var.
INTEGER,                INTENT(IN)   :: KRRL          ! number of liquid water var.
INTEGER,                INTENT(IN)   :: KRRI          ! number of ice water var.
LOGICAL,                INTENT(IN)   ::  OTURB_FLX    ! switch to write the
                                 ! turbulent fluxes in the syncronous FM-file
LOGICAL,                INTENT(IN)   ::  OOCEAN       ! switch for Ocean model version
LOGICAL,                INTENT(IN)   ::  ODEEPOC      ! activates sfc forcing for ideal ocean deep conv
LOGICAL,                INTENT(IN)   ::  OHARAT
LOGICAL,                INTENT(IN)   ::  OCOMPUTE_SRC ! flag to define dimensions of SIGS and
CHARACTER(LEN=4),       INTENT(IN)   ::  HTURBDIM     ! dimensionality of the
                                                      ! turbulence scheme
CHARACTER(LEN=6), INTENT(IN) :: HPROGRAM ! CPROGRAM is the program currently running (modd_conf)
CHARACTER(LEN=4),       INTENT(IN)   ::  HTOM         ! type of Third Order Moment
REAL,                   INTENT(IN)   ::  PIMPL, PEXPL ! Coef. for temporal disc.
REAL,                   INTENT(IN)   ::  PTSTEP       ! Double Time Step
TYPE(TFILEDATA),        INTENT(IN)   ::  TPFILE       ! Output file
LOGICAL,                INTENT(IN)   ::  OLES_CALL    ! compute the LES diagnostics at current time-step
LOGICAL,                INTENT(IN)   ::  OCOUPLES     ! switch to activate atmos-ocean LES version 
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PDZZ, PDXX, PDYY, PDZX, PDZY
                                                      ! Metric coefficients
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)   ::  PDIRCOSZW    ! Director Cosinus of the
                                                      ! normal to the ground surface
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PZZ          ! altitudes
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PRHODJ       ! dry density * grid volum
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  MFMOIST      ! moist mass flux dual scheme
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PTHVREF      ! ref. state Virtual 
                                                      ! Potential Temperature 
!
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)   ::  PSFTHM,PSFRM ! surface fluxes at time
!                                                     ! t - deltat 
!
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)   ::  PSFTHP,PSFRP ! surface fluxes at time
!                                                     ! t + deltat 
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PWM 
! Vertical wind
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PTHLM 
! potential temperature at t-Delta t
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KRR), INTENT(IN) ::  PRM          ! Mixing ratios 
                                                      ! at t-Delta t
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KSV), INTENT(IN) ::  PSVM         ! Mixing ratios 
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PTKEM        ! TKE at time t
!
! In case OHARAT=TRUE, PLM already includes all stability corrections
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PLM          ! Turb. mixing length   
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PLEPS        ! dissipative length   
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PLOCPEXNM    ! Lv(T)/Cp/Exnref at time t-1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PATHETA      ! coefficients between 
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PAMOIST      ! s and Thetal and Rnp
REAL, DIMENSION(MERGE(D%NIT,0,OCOMPUTE_SRC),&
                MERGE(D%NJT,0,OCOMPUTE_SRC),&
                MERGE(D%NKT,0,OCOMPUTE_SRC)), INTENT(IN)   ::  PSRCM        ! normalized 
! 2nd-order flux s'r'c/2Sigma_s2 at t-1 multiplied by Lambda_3
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PFRAC_ICE    ! ri fraction of rc+ri
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PBETA        ! buoyancy coefficient
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PSQRT_TKE    ! sqrt(e)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PDTH_DZ      ! d(th)/dz
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PDR_DZ       ! d(rt)/dz
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PRED2TH3     ! 3D Redeslperger number R*2_th
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PRED2R3      ! 3D Redeslperger number R*2_r
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PRED2THR3    ! 3D Redeslperger number R*2_thr
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PBLL_O_E     ! beta * Lk * Leps / tke
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PETHETA      ! Coefficient for theta in theta_v computation
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PEMOIST      ! Coefficient for r in theta_v computation
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PREDTH1      ! 1D Redelsperger number for Th
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PREDR1       ! 1D Redelsperger number for r
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PPHI3        ! Prandtl number for temperature
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PPSI3        ! Prandtl number for vapor
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PD           ! Denominator in Prandtl numbers
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PFWTH        ! d(w'2th' )/dz (at flux point)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PFWR         ! d(w'2r'  )/dz (at flux point)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PFTH2        ! d(w'th'2 )/dz (at mass point)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PFR2         ! d(w'r'2  )/dz (at mass point)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PFTHR        ! d(w'th'r')/dz (at mass point)
REAL, DIMENSION(MERGE(D%NIT,0,HTOM=='TM06'),&
                MERGE(D%NJT,0,HTOM=='TM06')),   INTENT(INOUT)::  PBL_DEPTH    ! BL depth
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT)  :: PWTHV         ! buoyancy flux
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(INOUT) :: PRTHLS     ! cumulated source for theta
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KRR), INTENT(INOUT) :: PRRS       ! cumulated source for rt
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(OUT)   :: PTHLP      ! guess of thl at t+ deltat
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(OUT)   :: PRP        ! guess of r at t+ deltat
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(OUT)   :: PTP       ! Dynamic and thermal
                                                     ! TKE production terms
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(OUT)   :: PWTH       ! heat flux
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(OUT)   :: PWRC       ! cloud water flux
!
!
!*       0.2  declaration of local variables
!
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT)  ::  &
       ZA,       & ! work variable for wrc or LES computation
       ZFLXZ,    & ! vertical flux of the treated variable
       ZSOURCE,  & ! source of evolution for the treated variable
       ZKEFF,    & ! effectif diffusion coeff = LT * SQRT( TKE )
       ZF,       & ! Flux in dTh/dt =-dF/dz (evaluated at t-1)(or rt instead of Th)
       ZDFDDTDZ, & ! dF/d(dTh/dz)
       ZDFDDRDZ, & ! dF/d(dr/dz)
       Z3RDMOMENT,&  ! 3 order term in flux or variance equation
       ZF_LEONARD,&  ! Leonard terms
       ZRWTHL,    &
       ZRWRNP,    &
       ZCLD_THOLD,&
       ZALT,      &
       ZWORK1,ZWORK2, &
       ZWORK3,ZWORK4 ! working var. for shuman operators (array syntax)
!
INTEGER             :: IKB,IKE      ! I index values for the Beginning and End
                                    ! mass points of the domain in the 3 direct.
INTEGER             :: IKT          ! array size in k direction
INTEGER             :: IKTB,IKTE    ! start, end of k loops in physical domain 
INTEGER             :: JI, JJ, JK ! loop indexes 
!
INTEGER                    :: IIB,IJB       ! Lower bounds of the physical
                                            ! sub-domain in x and y directions
INTEGER                    :: IIE,IJE       ! Upper bounds of the physical
                                            ! sub-domain in x and y directions
!
! NIMPORTE QUOI : TODO TO BE REMOVED OUTSIDE OF TURB ? :
REAL, DIMENSION(1)  :: ZXHAT_ll  !  Position x in the conformal
                                                 ! plane (array on the complete domain)
REAL, DIMENSION(1)  :: ZYHAT_ll  !   Position y in the conformal
                                                 ! plane (array on the complete domain)
!
!
REAL :: ZTIME1, ZTIME2
REAL :: ZDELTAX
REAL :: ZXBEG,ZXEND,ZYBEG,ZYEND ! Forcing size for ocean deep convection
REAL, DIMENSION(D%NIT,D%NJT) :: ZDIST ! distance
                                   ! from the center of the cooling               
REAL :: ZFLPROV
INTEGER           :: JKM          ! vertical index loop
INTEGER           :: JSW
REAL :: ZSWA     ! index for time flux interpolation
!
INTEGER :: IIU, IJU
INTEGER :: IRESP
LOGICAL :: GUSERV   ! flag to use water
LOGICAL :: GFTH2    ! flag to use w'th'2
LOGICAL :: GFWTH    ! flag to use w'2th'
LOGICAL :: GFR2     ! flag to use w'r'2
LOGICAL :: GFWR     ! flag to use w'2r'
LOGICAL :: GFTHR    ! flag to use w'th'r'
TYPE(TFIELDDATA) :: TZFIELD
!----------------------------------------------------------------------------
!
!*       1.   PRELIMINARIES
!             -------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('TURB_VER_THERMO_FLUX',0,ZHOOK_HANDLE)
!
! Size for a given proc & a given model      
IIU=D%NIT
IJU=D%NJT
IIE=D%NIEC
IIB=D%NIBC
IJE=D%NJEC
IJB=D%NJBC
!
!! Compute Shape of sfc flux for Oceanic Deep Conv Case
! 
IF (OOCEAN .AND. ODEEPOC) THEN
  !*       COMPUTES THE PHYSICAL SUBDOMAIN BOUNDS
  !compute ZXHAT_ll = position in the (0:Lx) domain 1 (Lx=Size of domain1 )
  !compute XXHAT_ll = position in the (L0_subproc,Lx_subproc) domain for the current subproc
  !                                     L0_subproc as referenced in the full domain 1
  CALL GATHERALL_FIELD_ll('XX',XXHAT,ZXHAT_ll,IRESP)
  CALL GATHERALL_FIELD_ll('YY',XYHAT,ZYHAT_ll,IRESP)
  CALL GET_DIM_EXT_ll('B',IIU,IJU)
  CALL GET_INDICE_ll(IIB,IJB,IIE,IJE,IIU,IJU)
  DO JJ = IJB,IJE
    DO JI = IIB,IIE
      ZDIST(JI,JJ) = SQRT(                         &
      (( (XXHAT(JI)+XXHAT(JI+1))*0.5 - XCENTX_OC ) / XRADX_OC)**2 + &
      (( (XYHAT(JJ)+XYHAT(JJ+1))*0.5 - XCENTY_OC ) / XRADY_OC)**2   &
                                )
    END DO
  END DO
  DO JJ=IJB,IJE
    DO JI=IIB,IIE
      IF ( ZDIST(JI,JJ) > 1.) XSSTFL(JI,JJ)=0.
    END DO
  END DO
END IF !END DEEP OCEAN CONV CASE
!
IKT=D%NKT  
IKTB=D%NKTB          
IKTE=D%NKTE
IKB=D%NKB
IKE=D%NKE
!
GUSERV = (KRR/=0)
!
!  compute the coefficients for the uncentred gradient computation near the 
!  ground
!
IF (OHARAT) THEN
! OHARAT so TKE and length scales at half levels!
  DO JK=1,D%NKT 
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE 
        ZKEFF(JI,JJ,JK) =  PLM(JI,JJ,JK) * SQRT(PTKEM(JI,JJ,JK)) & 
        +50.*MFMOIST(JI,JJ,JK)
      ENDDO
    ENDDO
  ENDDO
ELSE
  DO JK=1,D%NKT 
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE 
        ZWORK1(JI,JJ,JK) = PLM(JI,JJ,JK) * SQRT(PTKEM(JI,JJ,JK))
      ENDDO
    ENDDO
  ENDDO
  CALL MZM_PHY(D,ZWORK1,ZKEFF)
ENDIF
!
! Define a cloud mask with ri and rc (used after with a threshold) for Leonard terms
!
IF(TURBN%LHGRAD) THEN
  IF ( KRRL >= 1 ) THEN
    IF ( KRRI >= 1 ) THEN
      DO JK=1,D%NKT 
        DO JJ=IJB,IJE 
          DO JI=IIB,IIE 
            ZCLD_THOLD(JI,JJ,JK) = PRM(JI,JJ,JK,2) + PRM(JI,JJ,JK,4)
          ENDDO
        ENDDO
      ENDDO
    ELSE
      DO JK=1,D%NKT 
        DO JJ=IJB,IJE 
          DO JI=IIB,IIE 
            ZCLD_THOLD(JI,JJ,JK) = PRM(JI,JJ,JK,2)
          ENDDO
        ENDDO
      ENDDO
    END IF
  END IF
END IF
!
! Flags for 3rd order quantities
!
GFTH2 = .FALSE.
GFR2  = .FALSE.
GFTHR = .FALSE.
GFWTH = .FALSE.
GFWR  = .FALSE.
!
IF (HTOM/='NONE') THEN
  GFTH2 = ANY(PFTH2/=0.)
  GFR2  = ANY(PFR2 /=0.) .AND. GUSERV
  GFTHR = ANY(PFTHR/=0.) .AND. GUSERV
  GFWTH = ANY(PFWTH/=0.)
  GFWR  = ANY(PFWR /=0.) .AND. GUSERV
END IF
!----------------------------------------------------------------------------
!
!*       2.   SOURCES OF CONSERVATIVE POTENTIAL TEMPERATURE AND 
!                                                  PARTIAL THERMAL PRODUCTION 
!             ---------------------------------------------------------------
!
!*       2.1  Splitted value for cons. potential temperature at t+deltat
!
! Compute the turbulent flux F and F' at time t-dt.
!
CALL DZM_PHY(D,PTHLM,ZWORK1)
ZWORK2 = D_PHI3DTDZ_O_DDTDZ(D,CSTURB,PPHI3,PREDTH1,PREDR1,PRED2TH3,PRED2THR3,HTURBDIM,GUSERV)
IF (OHARAT) THEN
  DO JK=1,D%NKT 
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE 
        ZF(JI,JJ,JK) = -ZKEFF(JI,JJ,JK)*ZWORK1(JI,JJ,JK)/PDZZ(JI,JJ,JK)
        ZDFDDTDZ(JI,JJ,JK) = -ZKEFF(JI,JJ,JK)
      ENDDO
    ENDDO
  ENDDO
ELSE
  DO JK=1,D%NKT 
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE 
        ZF(JI,JJ,JK) = -CSTURB%XCSHF*PPHI3(JI,JJ,JK)*ZKEFF(JI,JJ,JK)& 
        *ZWORK1(JI,JJ,JK)/PDZZ(JI,JJ,JK)
        ZDFDDTDZ(JI,JJ,JK) = -CSTURB%XCSHF*ZKEFF(JI,JJ,JK)*ZWORK2(JI,JJ,JK)
      ENDDO
    ENDDO
  ENDDO
END IF
!
IF (TURBN%LHGRAD) THEN
 ! Compute the Leonard terms for thl
 ZDELTAX= XXHAT(3) - XXHAT(2)
 ZF_LEONARD (:,:,:)= TURBN%XCOEFHGRADTHL*ZDELTAX*ZDELTAX/12.0*(      &
                 MXF(GX_W_UW(PWM(:,:,:), XDXX,XDZZ,XDZX,D%NKA,D%NKU,D%NKL))&
                *MZM(GX_M_M(PTHLM(:,:,:),XDXX,XDZZ,XDZX,D%NKA, D%NKU, D%NKL), D%NKA, D%NKU, D%NKL)  &
              +  MYF(GY_W_VW(PWM(:,:,:), XDYY,XDZZ,XDZY,D%NKA,D%NKU,D%NKL))  &
                *MZM(GY_M_M(PTHLM(:,:,:),XDYY,XDZZ,XDZY,D%NKA, D%NKU, D%NKL), D%NKA, D%NKU, D%NKL) )
END IF
!
! Effect of 3rd order terms in temperature flux (at flux point)
!
! d(w'2th')/dz
IF (GFWTH) THEN
  Z3RDMOMENT= M3_WTH_W2TH(D,CSTURB,PREDTH1,PREDR1,PD,ZKEFF,PTKEM)
  ZWORK1 = D_M3_WTH_W2TH_O_DDTDZ(D,CSTURB,PREDTH1,PREDR1,&
   & PD,PBLL_O_E,PETHETA,ZKEFF,PTKEM)
!
  DO JK=1,D%NKT 
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE 
        ZF(JI,JJ,JK)= ZF(JI,JJ,JK)       + Z3RDMOMENT(JI,JJ,JK) * PFWTH(JI,JJ,JK)
        ZDFDDTDZ(JI,JJ,JK) = ZDFDDTDZ(JI,JJ,JK) + ZWORK1(JI,JJ,JK) * PFWTH(JI,JJ,JK)
      ENDDO
    ENDDO
  ENDDO
END IF
!
! d(w'th'2)/dz
IF (GFTH2) THEN
  Z3RDMOMENT= M3_WTH_WTH2(D,CSTURB,PREDTH1,PREDR1,PD,PBLL_O_E,PETHETA)
  ZWORK1 = D_M3_WTH_WTH2_O_DDTDZ(D,CSTURB,Z3RDMOMENT,PREDTH1,PREDR1,&
    & PD,PBLL_O_E,PETHETA)
  ZWORK2 = MZM(PFTH2, D%NKA, D%NKU, D%NKL)
!
  DO JK=1,D%NKT 
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE 
        ZF(JI,JJ,JK)       = ZF(JI,JJ,JK)       + Z3RDMOMENT(JI,JJ,JK) * ZWORK2(JI,JJ,JK)
        ZDFDDTDZ(JI,JJ,JK) = ZDFDDTDZ(JI,JJ,JK) + ZWORK1(JI,JJ,JK) * ZWORK2(JI,JJ,JK)
      ENDDO
    ENDDO
  ENDDO
END IF
!
! d(w'2r')/dz
IF (GFWR) THEN
  ZWORK1 = M3_WTH_W2R(D,CSTURB,PD,ZKEFF,PTKEM,PBLL_O_E,PEMOIST,PDTH_DZ)
  ZWORK2 = D_M3_WTH_W2R_O_DDTDZ(D,CSTURB,PREDTH1,PREDR1,PD,ZKEFF,PTKEM,PBLL_O_E,PEMOIST)
!
  DO JK=1,D%NKT 
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE 
        ZF(JI,JJ,JK)       = ZF(JI,JJ,JK)       + ZWORK1(JI,JJ,JK) * PFWR(JI,JJ,JK)
        ZDFDDTDZ(JI,JJ,JK) = ZDFDDTDZ(JI,JJ,JK) + ZWORK2(JI,JJ,JK) * PFWR(JI,JJ,JK)
      ENDDO
    ENDDO
  ENDDO
END IF
!
! d(w'r'2)/dz
IF (GFR2) THEN
  ZWORK1 = M3_WTH_WR2(D,CSTURB,PD,ZKEFF,PTKEM,PSQRT_TKE,PBLL_O_E,PBETA,PLEPS,PEMOIST,PDTH_DZ)
  ZWORK2 = MZM(PFR2, D%NKA, D%NKU, D%NKL)
  ZWORK3 = D_M3_WTH_WR2_O_DDTDZ(D,CSTURB,PREDTH1,PREDR1,PD,&
    & ZKEFF,PTKEM,PSQRT_TKE,PBLL_O_E,PBETA,PLEPS,PEMOIST)
!
  DO JK=1,D%NKT 
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE     
        ZF(JI,JJ,JK)       = ZF(JI,JJ,JK)       + ZWORK1(JI,JJ,JK) * ZWORK2(JI,JJ,JK)
        ZDFDDTDZ(JI,JJ,JK) = ZDFDDTDZ(JI,JJ,JK) + ZWORK3(JI,JJ,JK) * ZWORK2(JI,JJ,JK)
      ENDDO
    ENDDO
  ENDDO
END IF
!
! d(w'th'r')/dz
IF (GFTHR) THEN
  Z3RDMOMENT= M3_WTH_WTHR(D,CSTURB,PREDR1,PD,ZKEFF,PTKEM,PSQRT_TKE,PBETA,&
    & PLEPS,PEMOIST)
  ZWORK1 = D_M3_WTH_WTHR_O_DDTDZ(D,CSTURB,Z3RDMOMENT,PREDTH1,PREDR1,PD,PBLL_O_E,PETHETA)
  ZWORK2 = MZM(PFTHR, D%NKA, D%NKU, D%NKL)
!
  DO JK=1,D%NKT 
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE 
        ZF(JI,JJ,JK)       = ZF(JI,JJ,JK)       + Z3RDMOMENT(JI,JJ,JK) * ZWORK2(JI,JJ,JK)
        ZDFDDTDZ(JI,JJ,JK) = ZDFDDTDZ(JI,JJ,JK) + ZWORK1(JI,JJ,JK) * ZWORK2(JI,JJ,JK)
      ENDDO
    ENDDO
  ENDDO
END IF
! compute interface flux
IF (OCOUPLES) THEN   ! Autocoupling O-A LES
  IF (OOCEAN) THEN    ! ocean model in coupled case
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE  
        ZF(JI,JJ,IKE) =  (TURBN%XSSTFL_C(JI,JJ,1)+TURBN%XSSRFL_C(JI,JJ,1)) &
        *0.5* ( 1. + PRHODJ(JI,JJ,D%NKU)/PRHODJ(JI,JJ,IKE) )
      ENDDO
    ENDDO 
  ELSE                ! atmosph model in coupled case
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE  
        ZF(JI,JJ,IKB) =  TURBN%XSSTFL_C(JI,JJ,1) &
        *0.5* ( 1. + PRHODJ(JI,JJ,D%NKA)/PRHODJ(JI,JJ,IKB) )
      ENDDO
    ENDDO 
  ENDIF 
!
ELSE  ! No coupling O and A cases
  ! atmosp bottom
  !*In 3D, a part of the flux goes vertically,
  ! and another goes horizontally (in presence of slopes)
  !*In 1D, part of energy released in horizontal flux is taken into account in the vertical part
  IF (HTURBDIM=='3DIM') THEN
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE  
        ZF(JI,JJ,IKB) = ( PIMPL*PSFTHP(JI,JJ) + PEXPL*PSFTHM(JI,JJ) )   &
        * PDIRCOSZW(JI,JJ)                       &
        * 0.5 * (1. + PRHODJ(JI,JJ,D%NKA) / PRHODJ(JI,JJ,IKB))
      ENDDO
    ENDDO 
  ELSE
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE  
        ZF(JI,JJ,IKB) = ( PIMPL*PSFTHP(JI,JJ) + PEXPL*PSFTHM(JI,JJ) )   &
        / PDIRCOSZW(JI,JJ)                       &
        * 0.5 * (1. + PRHODJ(JI,JJ,D%NKA) / PRHODJ(JI,JJ,IKB))
      ENDDO
    ENDDO 
  END IF
!
  IF (OOCEAN) THEN
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE 
        ZF(JI,JJ,IKE) = XSSTFL(JI,JJ) *0.5*(1. + PRHODJ(JI,JJ,D%NKU) / PRHODJ(JI,JJ,IKE))
      ENDDO
    ENDDO
  ELSE !end ocean case (in nocoupled case)
    ! atmos top
#ifdef REPRO48
#else
      ZF(IIB:IIE,IJB:IJE,IKE)=0.
#endif
  END IF
END IF !end no coupled cases
!
! Compute the split conservative potential temperature at t+deltat
CALL TRIDIAG_THERMO(D,PTHLM,ZF,ZDFDDTDZ,PTSTEP,PIMPL,PDZZ,&
                    PRHODJ,PTHLP)
!
! Compute the equivalent tendency for the conservative potential temperature
!
DO JK=1,D%NKT 
  DO JJ=IJB,IJE 
    DO JI=IIB,IIE     
      ZRWTHL(JI,JJ,JK)= PRHODJ(JI,JJ,JK)*(PTHLP(JI,JJ,JK)-PTHLM(JI,JJ,JK))& 
      /PTSTEP
    ENDDO
  ENDDO
ENDDO    
! replace the flux by the Leonard terms above ZALT and ZCLD_THOLD
IF (TURBN%LHGRAD) THEN
 DO JK=1,D%NKU
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE 
        ZALT(JI,JJ,JK) = PZZ(JI,JJ,JK)-XZS(JI,JJ)
      ENDDO
    ENDDO
 END DO
 ZWORK1 = GZ_W_M(MZM(PRHODJ, D%NKA, D%NKU, D%NKL)*ZF_LEONARD(:,:,:),XDZZ,&
                   D%NKA, D%NKU, D%NKL)
  DO JK=1,D%NKT 
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE 
        IF ( (ZCLD_THOLD(JI,JJ,JK) >= TURBN%XCLDTHOLD) .AND. ( ZALT(JI,JJ,JK) >= TURBN%XALTHGRAD) )THEN
          ZRWTHL(JI,JJ,JK) = -ZWORK1(JI,JJ,JK)
        ENDIF
      ENDDO
    ENDDO
  ENDDO
END IF
!
DO JK=1,D%NKT 
  DO JJ=IJB,IJE 
    DO JI=IIB,IIE 
      ZWORK1(JI,JJ,JK) = PTHLP(JI,JJ,JK) - PTHLM(JI,JJ,JK)
    ENDDO
  ENDDO
ENDDO
CALL DZM_PHY(D,ZWORK1,ZWORK2)
!
DO JK=1,D%NKT 
  DO JJ=IJB,IJE 
    DO JI=IIB,IIE 
      PRTHLS(JI,JJ,JK)= PRTHLS(JI,JJ,JK)  + ZRWTHL(JI,JJ,JK)
!
!*       2.2  Partial Thermal Production
!
!  Conservative potential temperature flux : 
!
!
      ZFLXZ(JI,JJ,JK)   = ZF(JI,JJ,JK) + PIMPL * ZDFDDTDZ(JI,JJ,JK) * & 
      ZWORK2(JI,JJ,JK)/ PDZZ(JI,JJ,JK)
    ENDDO
  ENDDO
ENDDO
!
! replace the flux by the Leonard terms
IF (TURBN%LHGRAD) THEN
  DO JK=1,D%NKT 
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE 
        IF ( (ZCLD_THOLD(JI,JJ,JK) >= TURBN%XCLDTHOLD) .AND. ( ZALT(JI,JJ,JK) >= TURBN%XALTHGRAD) )THEN
          ZFLXZ(JI,JJ,JK) = ZF_LEONARD(JI,JJ,JK)
        ENDIF
      ENDDO
    ENDDO
  ENDDO
END IF
!
DO JJ=IJB,IJE 
  DO JI=IIB,IIE 
    ZFLXZ(JI,JJ,D%NKA) = ZFLXZ(JI,JJ,IKB)
  ENDDO
ENDDO
IF (OOCEAN) THEN
  DO JJ=IJB,IJE 
    DO JI=IIB,IIE 
      ZFLXZ(JI,JJ,D%NKU) = ZFLXZ(JI,JJ,IKE)
    ENDDO
  ENDDO
END IF
!
DO JK=IKTB+1,IKTE-1
  DO JJ=IJB,IJE 
    DO JI=IIB,IIE 
      PWTH(JI,JJ,JK)=0.5*(ZFLXZ(JI,JJ,JK)+ZFLXZ(JI,JJ,JK+D%NKL))
    ENDDO
  ENDDO
END DO
!
DO JJ=IJB,IJE 
  DO JI=IIB,IIE 
    PWTH(JI,JJ,IKB)=0.5*(ZFLXZ(JI,JJ,IKB)+ZFLXZ(JI,JJ,IKB+D%NKL)) 
  ENDDO
ENDDO    
!
IF (OOCEAN) THEN
  DO JJ=IJB,IJE 
    DO JI=IIB,IIE 
      PWTH(JI,JJ,IKE)=0.5*(ZFLXZ(JI,JJ,IKE)+ZFLXZ(JI,JJ,IKE+D%NKL))
      PWTH(JI,JJ,D%NKA)=0. 
      PWTH(JI,JJ,D%NKU)=ZFLXZ(JI,JJ,D%NKU)
    ENDDO
  ENDDO
ELSE
  DO JJ=IJB,IJE 
    DO JI=IIB,IIE 
      PWTH(JI,JJ,D%NKA)=0.5*(ZFLXZ(JI,JJ,D%NKA)+ZFLXZ(JI,JJ,D%NKA+D%NKL))
      PWTH(JI,JJ,IKE)=PWTH(JI,JJ,IKE-D%NKL)
    ENDDO
  ENDDO
END IF
!
IF ( OTURB_FLX .AND. TPFILE%LOPENED ) THEN
  ! stores the conservative potential temperature vertical flux
  TZFIELD%CMNHNAME   = 'THW_FLX'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'THW_FLX'
  TZFIELD%CUNITS     = 'K m s-1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'Conservative potential temperature vertical flux'
  TZFIELD%NGRID      = 4
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZFLXZ)
END IF
!
! Contribution of the conservative temperature flux to the buoyancy flux
IF (OOCEAN) THEN
  CALL MZF_PHY(D,ZFLXZ,ZWORK1)
  DO JK=1,D%NKT 
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE 
        PTP(JI,JJ,JK)= CST%XG*CST%XALPHAOC * ZWORK1(JI,JJ,JK)
      ENDDO
    ENDDO
  ENDDO
ELSE
  IF (KRR /= 0) THEN
    CALL MZM_PHY(D,PETHETA,ZWORK1)
    ZWORK1(IIB:IIE,IJB:IJE,1:D%NKT) = ZWORK1(IIB:IIE,IJB:IJE,1:D%NKT) * ZFLXZ(IIB:IIE,IJB:IJE,1:D%NKT)
    CALL MZF_PHY(D,ZWORK1,ZWORK2)
    !ZWORK1 = MZF( MZM(PETHETA,D%NKA, D%NKU, D%NKL) * ZFLXZ,D%NKA, D%NKU, D%NKL )
    DO JK=1,D%NKT 
      DO JJ=IJB,IJE 
        DO JI=IIB,IIE 
          PTP(JI,JJ,JK)  =  PBETA(JI,JJ,JK) * ZWORK2(JI,JJ,JK)
        ENDDO
      ENDDO
    ENDDO
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE 
        PTP(JI,JJ,IKB)=  PBETA(JI,JJ,IKB) * PETHETA(JI,JJ,IKB) *   &
        0.5 * ( ZFLXZ(JI,JJ,IKB) + ZFLXZ(JI,JJ,IKB+D%NKL) )
      ENDDO
    ENDDO
  ELSE
    CALL MZF_PHY(D,ZFLXZ,ZWORK1)
    DO JK=1,D%NKT 
      DO JJ=IJB,IJE 
        DO JI=IIB,IIE 
          PTP(JI,JJ,JK)=  PBETA(JI,JJ,JK) * ZWORK1(JI,JJ,JK)
        ENDDO
      ENDDO
    ENDDO
  END IF
END IF 
!
! Buoyancy flux at flux points
!
CALL MZM_PHY(D,PETHETA,ZWORK1)
DO JK=1,D%NKT 
  DO JJ=IJB,IJE 
    DO JI=IIB,IIE 
      PWTHV(JI,JJ,JK) = ZWORK1(JI,JJ,JK) * ZFLXZ(JI,JJ,JK)
    ENDDO
  ENDDO
ENDDO
DO JJ=IJB,IJE 
  DO JI=IIB,IIE 
    PWTHV(JI,JJ,IKB) = PETHETA(JI,JJ,IKB) * ZFLXZ(JI,JJ,IKB)
  ENDDO
ENDDO
!
IF (OOCEAN) THEN
  ! temperature contribution to Buy flux
  DO JJ=IJB,IJE 
    DO JI=IIB,IIE  
      PWTHV(JI,JJ,IKE) = PETHETA(JI,JJ,IKE) * ZFLXZ(JI,JJ,IKE)
    ENDDO
  ENDDO
END IF
!*       2.3  Partial vertical divergence of the < Rc w > flux
! Correction for qc and qi negative in AROME 
IF(HPROGRAM/='AROME  ') THEN
 IF ( KRRL >= 1 ) THEN
    DO JK=1,D%NKT 
      DO JJ=IJB,IJE 
        DO JI=IIB,IIE 
          ZWORK1(JI,JJ,JK) = ZFLXZ(JI,JJ,JK)/PDZZ(JI,JJ,JK)
        ENDDO
      ENDDO
    ENDDO
   CALL DZF_PHY(D,ZWORK1,ZWORK2)
   IF ( KRRI >= 1 ) THEN
      DO JK=1,D%NKT 
        DO JJ=IJB,IJE 
          DO JI=IIB,IIE 
            PRRS(JI,JJ,JK,2) = PRRS(JI,JJ,JK,2) -                                        &
            PRHODJ(JI,JJ,JK)*PATHETA(JI,JJ,JK)*2.*PSRCM(JI,JJ,JK)& 
            *ZWORK2(JI,JJ,JK) *(1.0-PFRAC_ICE(JI,JJ,JK))
            PRRS(JI,JJ,JK,4) = PRRS(JI,JJ,JK,4) -                                        &
            PRHODJ(JI,JJ,JK)*PATHETA(JI,JJ,JK)*2.*PSRCM(JI,JJ,JK)&
            * ZWORK2(JI,JJ,JK)*PFRAC_ICE(JI,JJ,JK)
          ENDDO
        ENDDO
      ENDDO
   ELSE
      DO JK=1,D%NKT 
        DO JJ=IJB,IJE 
          DO JI=IIB,IIE 
            PRRS(JI,JJ,JK,2) = PRRS(JI,JJ,JK,2) -                                        &
            PRHODJ(JI,JJ,JK)*PATHETA(JI,JJ,JK)*2.*PSRCM(JI,JJ,JK)&
            *ZWORK2(JI,JJ,JK)
          ENDDO
        ENDDO
      ENDDO
   END IF
 END IF
END IF
!
!*       2.4  Storage in LES configuration
! 
IF (OLES_CALL) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL LES_MEAN_SUBGRID(MZF(ZFLXZ, D%NKA, D%NKU, D%NKL), X_LES_SUBGRID_WThl ) 
  CALL LES_MEAN_SUBGRID(MZF(PWM*ZFLXZ, D%NKA, D%NKU, D%NKL), X_LES_RES_W_SBG_WThl )
  CALL LES_MEAN_SUBGRID(GZ_W_M(PWM,PDZZ, D%NKA, D%NKU, D%NKL)*MZF(ZFLXZ, D%NKA, D%NKU, D%NKL),&
      & X_LES_RES_ddxa_W_SBG_UaThl )
  CALL LES_MEAN_SUBGRID(MZF(PDTH_DZ*ZFLXZ, D%NKA, D%NKU, D%NKL), X_LES_RES_ddxa_Thl_SBG_UaThl )
  CALL LES_MEAN_SUBGRID(-CSTURB%XCTP*PSQRT_TKE/PLM*MZF(ZFLXZ, D%NKA, D%NKU, D%NKL), X_LES_SUBGRID_ThlPz ) 
  CALL LES_MEAN_SUBGRID(MZF(MZM(PETHETA, D%NKA, D%NKU, D%NKL)*ZFLXZ, D%NKA, D%NKU, D%NKL), X_LES_SUBGRID_WThv ) 
  IF (KRR>=1) THEN
    CALL LES_MEAN_SUBGRID(MZF(PDR_DZ*ZFLXZ, D%NKA, D%NKU, D%NKL), X_LES_RES_ddxa_Rt_SBG_UaThl )
  END IF
  !* diagnostic of mixing coefficient for heat
  ZA = DZM(PTHLP, D%NKA, D%NKU, D%NKL)
  DO JK=1,D%NKT 
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE 
        IF (ZA(JI,JJ,JK)==0.) THEN
          ZA(JI,JJ,JK)=1.E-6
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  DO JK=1,D%NKT 
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE 
        ZA(JI,JJ,JK) = - ZFLXZ(JI,JJ,JK) / ZA(JI,JJ,JK) * PDZZ(JI,JJ,JK)
      ENDDO
    ENDDO
  ENDDO
  DO JJ=IJB,IJE 
    DO JI=IIB,IIE 
      ZA(JI,JJ,IKB) = CSTURB%XCSHF*PPHI3(JI,JJ,IKB)*ZKEFF(JI,JJ,IKB)
    ENDDO
  ENDDO
  ZA = MZF(ZA, D%NKA, D%NKU, D%NKL)
  DO JK=1,D%NKT 
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE 
        ZA(JI,JJ,JK) = MIN(MAX(ZA(JI,JJ,JK),-1000.),1000.)
      ENDDO
    ENDDO
  ENDDO
  CALL LES_MEAN_SUBGRID( ZA, X_LES_SUBGRID_Kh   ) 
  !
  CALL SECOND_MNH(ZTIME2)
  XTIME_LES = XTIME_LES + ZTIME2 - ZTIME1
END IF
!
!*       2.5  New boundary layer depth for TOMs
! 
IF (HTOM=='TM06') CALL TM06_H(IKB,IKTB,IKTE,PTSTEP,PZZ,ZFLXZ,PBL_DEPTH)
!
!----------------------------------------------------------------------------
!
!
!*       3.   SOURCES OF CONSERVATIVE AND CLOUD MIXING RATIO AND 
!                                        COMPLETE THERMAL PRODUCTION 
!             ------------------------------------------------------
!
!*       3.1  Splitted value for cons. mixing ratio at t+deltat
!
!
IF (KRR /= 0) THEN
  ! Compute the turbulent flux F and F' at time t-dt.
  !
  CALL DZM_PHY(D,PRM(:,:,:,1),ZWORK1)
 IF (OHARAT) THEN
    DO JK=1,D%NKT 
      DO JJ=IJB,IJE 
        DO JI=IIB,IIE     
          ZF(JI,JJ,JK) = -ZKEFF(JI,JJ,JK)*ZWORK1(JI,JJ,JK)/PDZZ(JI,JJ,JK)
          ZDFDDRDZ(JI,JJ,JK) = -ZKEFF(JI,JJ,JK)
        ENDDO
      ENDDO
    ENDDO    
 ELSE
  ZWORK2 = D_PSI3DRDZ_O_DDRDZ(D,CSTURB,PPSI3,PREDR1,PREDTH1,PRED2R3,PRED2THR3,HTURBDIM,GUSERV)
    DO JK=1,D%NKT 
      DO JJ=IJB,IJE 
        DO JI=IIB,IIE     
          ZF(JI,JJ,JK) = -CSTURB%XCSHF*PPSI3(JI,JJ,JK)*ZKEFF(JI,JJ,JK)& 
          *ZWORK1(JI,JJ,JK)/PDZZ(JI,JJ,JK)
          ZDFDDRDZ(JI,JJ,JK) = -CSTURB%XCSHF*ZKEFF(JI,JJ,JK)*ZWORK2(JI,JJ,JK)
        ENDDO
      ENDDO
    ENDDO    
 ENDIF
  !
  ! Compute Leonard Terms for Cloud mixing ratio
  IF (TURBN%LHGRAD) THEN
    ZDELTAX= XXHAT(3) - XXHAT(2)
    ZF_LEONARD (:,:,:)= TURBN%XCOEFHGRADRM*ZDELTAX*ZDELTAX/12.0*(        &
                MXF(GX_W_UW(PWM(:,:,:),  XDXX,XDZZ,XDZX,D%NKA,D%NKU,D%NKL))       &
                *MZM(GX_M_M(PRM(:,:,:,1),XDXX,XDZZ,XDZX,D%NKA,D%NKU,D%NKL),D%NKA,D%NKU,D%NKL) &
                +MYF(GY_W_VW(PWM(:,:,:), XDYY,XDZZ,XDZY,D%NKA,D%NKU,D%NKL))        &
                *MZM(GY_M_M(PRM(:,:,:,1),XDYY,XDZZ,XDZY,D%NKA,D%NKU,D%NKL),D%NKA,D%NKU,D%NKL) )
   END IF
  !
  ! Effect of 3rd order terms in temperature flux (at flux point)
  !
  ! d(w'2r')/dz
  IF (GFWR) THEN
    Z3RDMOMENT= M3_WR_W2R(D,CSTURB,PREDR1,PREDTH1,PD,ZKEFF,PTKEM)
    ZWORK1 = D_M3_WR_W2R_O_DDRDZ(D,CSTURB,PREDR1,PREDTH1,PD,&
     & PBLL_O_E,PEMOIST,ZKEFF,PTKEM)
  !
    DO JK=1,D%NKT 
      DO JJ=IJB,IJE 
        DO JI=IIB,IIE 
          ZF(JI,JJ,JK)       = ZF(JI,JJ,JK)       + Z3RDMOMENT(JI,JJ,JK) * PFWR(JI,JJ,JK)
          ZDFDDRDZ(JI,JJ,JK) = ZDFDDRDZ(JI,JJ,JK) + ZWORK1(JI,JJ,JK) * PFWR(JI,JJ,JK)
        ENDDO
      ENDDO
    ENDDO
  END IF
  !
  ! d(w'r'2)/dz
  IF (GFR2) THEN
    Z3RDMOMENT= M3_WR_WR2(D,CSTURB,PREDR1,PREDTH1,PD,PBLL_O_E,PEMOIST)
    ZWORK1 = MZM(PFR2, D%NKA, D%NKU, D%NKL)
    ZWORK2 = D_M3_WR_WR2_O_DDRDZ(D,CSTURB,Z3RDMOMENT,PREDR1,&
     & PREDTH1,PD,PBLL_O_E,PEMOIST)
  !
    DO JK=1,D%NKT 
      DO JJ=IJB,IJE 
        DO JI=IIB,IIE 
          ZF(JI,JJ,JK)       = ZF(JI,JJ,JK)       + Z3RDMOMENT(JI,JJ,JK) * ZWORK1(JI,JJ,JK)
          ZDFDDRDZ(JI,JJ,JK) = ZDFDDRDZ(JI,JJ,JK) + ZWORK2(JI,JJ,JK) * ZWORK1(JI,JJ,JK)
        ENDDO
      ENDDO
    ENDDO
  END IF
  !
  ! d(w'2th')/dz
  IF (GFWTH) THEN
    ZWORK1 = M3_WR_W2TH(D,CSTURB,PD,ZKEFF,&
     & PTKEM,PBLL_O_E,PETHETA,PDR_DZ)
    ZWORK2 = D_M3_WR_W2TH_O_DDRDZ(D,CSTURB,PREDR1,PREDTH1,& 
     & PD,ZKEFF,PTKEM,PBLL_O_E,PETHETA)
  !
    DO JK=1,D%NKT 
      DO JJ=IJB,IJE 
        DO JI=IIB,IIE 
          ZF(JI,JJ,JK)       = ZF(JI,JJ,JK)       + ZWORK1(JI,JJ,JK) * PFWTH(JI,JJ,JK)
          ZDFDDRDZ(JI,JJ,JK) = ZDFDDRDZ(JI,JJ,JK) + ZWORK2(JI,JJ,JK) * PFWTH(JI,JJ,JK)
        ENDDO
      ENDDO
    ENDDO
  END IF
  !
  ! d(w'th'2)/dz
  IF (GFTH2) THEN
    ZWORK1 = MZM(PFTH2, D%NKA, D%NKU, D%NKL)
    ZWORK2 = M3_WR_WTH2(D,CSTURB,PD,ZKEFF,PTKEM,&
    & PSQRT_TKE,PBLL_O_E,PBETA,PLEPS,PETHETA,PDR_DZ)
    ZWORK3 = D_M3_WR_WTH2_O_DDRDZ(D,CSTURB,PREDR1,PREDTH1,PD,&
     &ZKEFF,PTKEM,PSQRT_TKE,PBLL_O_E,PBETA,PLEPS,PETHETA)
    !
    DO JK=1,D%NKT 
      DO JJ=IJB,IJE 
        DO JI=IIB,IIE 
          ZF(JI,JJ,JK)       = ZF(JI,JJ,JK)       + ZWORK2(JI,JJ,JK) * ZWORK1(JI,JJ,JK)
          ZDFDDRDZ(JI,JJ,JK) = ZDFDDRDZ(JI,JJ,JK) + ZWORK3(JI,JJ,JK) * ZWORK1(JI,JJ,JK)
        ENDDO
      ENDDO
    ENDDO
  END IF
  !
  ! d(w'th'r')/dz
  IF (GFTHR) THEN
    Z3RDMOMENT= M3_WR_WTHR(D,CSTURB,PREDTH1,PD,ZKEFF,PTKEM,PSQRT_TKE,PBETA,&
     & PLEPS,PETHETA)
    ZWORK1 = MZM(PFTHR, D%NKA, D%NKU, D%NKL)
    ZWORK2 = D_M3_WR_WTHR_O_DDRDZ(D,CSTURB,Z3RDMOMENT,PREDR1, &
     & PREDTH1,PD,PBLL_O_E,PEMOIST)
  !
    DO JK=1,D%NKT 
      DO JJ=IJB,IJE 
        DO JI=IIB,IIE 
          ZF(JI,JJ,JK)       = ZF(JI,JJ,JK)       + Z3RDMOMENT(JI,JJ,JK) * ZWORK1(JI,JJ,JK)
          ZDFDDRDZ(JI,JJ,JK) = ZDFDDRDZ(JI,JJ,JK) + ZWORK2(JI,JJ,JK) * ZWORK1(JI,JJ,JK)
        ENDDO
      ENDDO
    ENDDO
  END IF
  !
  ! compute interface flux
  IF (OCOUPLES) THEN   ! coupling NH O-A
    IF (OOCEAN) THEN    ! ocean model in coupled case
      ! evap effect on salinity to be added later !!!
      ZF(IIB:IIE,IJB:IJE,IKE) =  0.
    ELSE                ! atmosph model in coupled case
      ZF(IIB:IIE,IJB:IJE,IKB) =  0.
      ! AJOUTER FLUX EVAP SUR MODELE ATMOS
    ENDIF
  !
  ELSE  ! No coupling NH OA case
    ! atmosp bottom
    !* in 3DIM case, a part of the flux goes vertically, and another goes horizontally
    ! (in presence of slopes)
    !* in 1DIM case, the part of energy released in horizontal flux
    ! is taken into account in the vertical part
    !
    IF (HTURBDIM=='3DIM') THEN
      DO JJ=IJB,IJE 
        DO JI=IIB,IIE  
          ZF(JI,JJ,IKB) = ( PIMPL*PSFRP(JI,JJ) + PEXPL*PSFRM(JI,JJ) )       &
          * PDIRCOSZW(JI,JJ)                       &
          * 0.5 * (1. + PRHODJ(JI,JJ,D%NKA) / PRHODJ(JI,JJ,IKB))
        ENDDO
      ENDDO 
    ELSE
      DO JJ=IJB,IJE 
        DO JI=IIB,IIE  
          ZF(JI,JJ,IKB) = ( PIMPL*PSFRP(JI,JJ) + PEXPL*PSFRM(JI,JJ) )     &
          / PDIRCOSZW(JI,JJ)                       &
          * 0.5 * (1. + PRHODJ(JI,JJ,D%NKA) / PRHODJ(JI,JJ,IKB))
        ENDDO
      ENDDO 
    END IF
    !
    IF (OOCEAN) THEN
      ! General ocean case
      ! salinity/evap effect to be added later !!!!!
      ZF(IIB:IIE,IJB:IJE,IKE) = 0.
    ELSE !end ocean case (in nocoupled case)
      ! atmos top
#ifdef REPRO48
#else
      ZF(IIB:IIE,IJB:IJE,IKE)=0.
#endif
    END IF
  END IF!end no coupled cases
  ! Compute the split conservative potential temperature at t+deltat
  CALL TRIDIAG_THERMO(D,PRM(:,:,:,1),ZF,ZDFDDRDZ,PTSTEP,PIMPL,&
                      PDZZ,PRHODJ,PRP)
  !
  ! Compute the equivalent tendency for the conservative mixing ratio
  !
  DO JK=1,D%NKT 
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE     
        ZRWRNP(JI,JJ,JK) = PRHODJ(JI,JJ,JK)*(PRP(JI,JJ,JK)-PRM(JI,JJ,JK,1))& 
        /PTSTEP
      ENDDO
    ENDDO
  ENDDO    
  !
  ! replace the flux by the Leonard terms above ZALT and ZCLD_THOLD
  IF (TURBN%LHGRAD) THEN
   DO JK=1,D%NKU
      DO JJ=IJB,IJE 
        DO JI=IIB,IIE 
          ZALT(JI,JJ,JK) = PZZ(JI,JJ,JK)-XZS(JI,JJ)
        ENDDO
      ENDDO
   END DO
   ZWORK1 = GZ_W_M(MZM(PRHODJ(:,:,:),D%NKA,D%NKU,D%NKL)*ZF_LEONARD(:,:,:),XDZZ,D%NKA,D%NKU,D%NKL)
    DO JK=1,D%NKT 
      DO JJ=IJB,IJE 
        DO JI=IIB,IIE 
          IF ( (ZCLD_THOLD(JI,JJ,JK) >= TURBN%XCLDTHOLD ) .AND. ( ZALT(JI,JJ,JK) >= TURBN%XALTHGRAD ) )THEN
            ZRWRNP(JI,JJ,JK) =  -ZWORK1(JI,JJ,JK)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  END IF
  !
  DO JK=1,D%NKT 
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE 
        ZWORK1(JI,JJ,JK) = PRP(JI,JJ,JK) - PRM(JI,JJ,JK,1)
      ENDDO
    ENDDO
  ENDDO
  CALL DZM_PHY(D,ZWORK1,ZWORK2)
  DO JK=1,D%NKT 
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE 
        PRRS(JI,JJ,JK,1) = PRRS(JI,JJ,JK,1) + ZRWRNP(JI,JJ,JK)
  !
  !*       3.2  Complete thermal production
  !
  ! cons. mixing ratio flux :
  !
        ZFLXZ(JI,JJ,JK)   = ZF(JI,JJ,JK)                                                &
        + PIMPL * ZDFDDRDZ(JI,JJ,JK) * ZWORK2(JI,JJ,JK) / PDZZ(JI,JJ,JK) 
      ENDDO
    ENDDO
  ENDDO
  !
  ! replace the flux by the Leonard terms above ZALT and ZCLD_THOLD
  IF (TURBN%LHGRAD) THEN
    DO JK=1,D%NKT 
      DO JJ=IJB,IJE 
        DO JI=IIB,IIE 
          IF ( (ZCLD_THOLD(JI,JJ,JK) >= TURBN%XCLDTHOLD ) .AND. ( ZALT(JI,JJ,JK) >= TURBN%XALTHGRAD ) )THEN
            ZFLXZ(JI,JJ,JK) = ZF_LEONARD(JI,JJ,JK)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  END IF
  !
  DO JJ=IJB,IJE 
    DO JI=IIB,IIE 
      ZFLXZ(JI,JJ,D%NKA) = ZFLXZ(JI,JJ,IKB) 
    ENDDO
  ENDDO
  !
  DO JK=IKTB+1,IKTE-1
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE 
        PWRC(JI,JJ,JK)=0.5*(ZFLXZ(JI,JJ,JK)+ZFLXZ(JI,JJ,JK+D%NKL))
      ENDDO
    ENDDO
  END DO
  DO JJ=IJB,IJE 
    DO JI=IIB,IIE 
      PWRC(JI,JJ,IKB)=0.5*(ZFLXZ(JI,JJ,IKB)+ZFLXZ(JI,JJ,IKB+D%NKL))
      PWRC(JI,JJ,D%NKA)=0.5*(ZFLXZ(JI,JJ,D%NKA)+ZFLXZ(JI,JJ,D%NKA+D%NKL))
      PWRC(JI,JJ,IKE)=PWRC(JI,JJ,IKE-D%NKL)
    ENDDO
  ENDDO
  !
  !
  IF ( OTURB_FLX .AND. TPFILE%LOPENED ) THEN
    ! stores the conservative mixing ratio vertical flux
    TZFIELD%CMNHNAME   = 'RCONSW_FLX'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'RCONSW_FLX'
    TZFIELD%CUNITS     = 'kg m s-1 kg-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'Conservative mixing ratio vertical flux'
    TZFIELD%NGRID      = 4
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZFLXZ)
  END IF
  !
  ! Contribution of the conservative water flux to the Buoyancy flux
  IF (OOCEAN) THEN     
     CALL MZF_PHY(D,ZFLXZ,ZWORK1)
    DO JK=1,D%NKT 
      DO JJ=IJB,IJE 
        DO JI=IIB,IIE 
          ZA(JI,JJ,JK)=  -CST%XG*CST%XBETAOC  * ZWORK1(JI,JJ,JK)
        ENDDO
      ENDDO
    ENDDO
  ELSE
    CALL MZM_PHY(D,PEMOIST,ZWORK1)
    DO JK=1,D%NKT 
      DO JJ=IJB,IJE 
        DO JI=IIB,IIE 
          ZWORK1(JI,JJ,JK) = ZWORK1(JI,JJ,JK) * ZFLXZ(JI,JJ,JK)
        ENDDO
      ENDDO
    ENDDO
    CALL MZF_PHY(D,ZWORK1,ZWORK2)
    !
    DO JK=1,D%NKT 
      DO JJ=IJB,IJE 
        DO JI=IIB,IIE 
          ZA(JI,JJ,JK)   =  PBETA(JI,JJ,JK) * ZWORK2(JI,JJ,JK)
        ENDDO
      ENDDO
    ENDDO
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE 
        ZA(JI,JJ,IKB) =  PBETA(JI,JJ,IKB) * PEMOIST(JI,JJ,IKB) *   &
        0.5 * ( ZFLXZ(JI,JJ,IKB) + ZFLXZ(JI,JJ,IKB+D%NKL) )
      ENDDO
    ENDDO
    DO JK=1,D%NKT 
      DO JJ=IJB,IJE 
        DO JI=IIB,IIE 
          PTP(JI,JJ,JK) = PTP(JI,JJ,JK) + ZA(JI,JJ,JK)
        ENDDO
      ENDDO
    ENDDO
  END IF
  !
  ! Buoyancy flux at flux points
  !
  CALL MZM_PHY(D,PEMOIST,ZWORK1)
  DO JK=1,D%NKT 
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE 
        PWTHV(JI,JJ,JK)=PWTHV(JI,JJ,JK) + ZWORK1(JI,JJ,JK) * ZFLXZ(JI,JJ,JK)
      ENDDO
    ENDDO
  ENDDO
  DO JJ=IJB,IJE 
    DO JI=IIB,IIE 
      PWTHV(JI,JJ,IKB) = PWTHV(JI,JJ,IKB) + PEMOIST(JI,JJ,IKB) * ZFLXZ(JI,JJ,IKB)
    ENDDO
  ENDDO
  IF (OOCEAN) THEN
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE 
        PWTHV(JI,JJ,IKE) = PWTHV(JI,JJ,IKE) + PEMOIST(JI,JJ,IKE)* ZFLXZ(JI,JJ,IKE)
      ENDDO
    ENDDO
  END IF   
!
!*       3.3  Complete vertical divergence of the < Rc w > flux
! Correction of qc and qi negative for AROME
IF(HPROGRAM/='AROME  ') THEN
   IF ( KRRL >= 1 ) THEN
       ZWORK1 = DZF(ZFLXZ/PDZZ,D%NKA,D%NKU,D%NKL )
     IF ( KRRI >= 1 ) THEN
        DO JK=1,D%NKT 
          DO JJ=IJB,IJE 
            DO JI=IIB,IIE 
              PRRS(JI,JJ,JK,2) = PRRS(JI,JJ,JK,2) -                                        &
              PRHODJ(JI,JJ,JK)*PAMOIST(JI,JJ,JK)*2.*PSRCM(JI,JJ,JK)& 
              *ZWORK1(JI,JJ,JK) *(1.0-PFRAC_ICE(JI,JJ,JK))
              PRRS(JI,JJ,JK,4) = PRRS(JI,JJ,JK,4) -                                        &
              PRHODJ(JI,JJ,JK)*PAMOIST(JI,JJ,JK)*2.*PSRCM(JI,JJ,JK)&
              *ZWORK1(JI,JJ,JK) *PFRAC_ICE(JI,JJ,JK)
            ENDDO
          ENDDO
        ENDDO
     ELSE
        DO JK=1,D%NKT 
          DO JJ=IJB,IJE 
            DO JI=IIB,IIE 
              PRRS(JI,JJ,JK,2) = PRRS(JI,JJ,JK,2) -                                        &
              PRHODJ(JI,JJ,JK)*PAMOIST(JI,JJ,JK)*2.*PSRCM(JI,JJ,JK)&
              *ZWORK1(JI,JJ,JK)
            ENDDO
          ENDDO
        ENDDO
     END IF
   END IF
END IF
!
!*       3.4  Storage in LES configuration
! 
  IF (OLES_CALL) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL LES_MEAN_SUBGRID(MZF(ZFLXZ, D%NKA, D%NKU, D%NKL), X_LES_SUBGRID_WRt ) 
    CALL LES_MEAN_SUBGRID(MZF(PWM*ZFLXZ, D%NKA, D%NKU, D%NKL), X_LES_RES_W_SBG_WRt )
    CALL LES_MEAN_SUBGRID(GZ_W_M(PWM,PDZZ, D%NKA, D%NKU, D%NKL)*MZF(ZFLXZ, D%NKA, D%NKU, D%NKL),&
    & X_LES_RES_ddxa_W_SBG_UaRt )
    CALL LES_MEAN_SUBGRID(MZF(PDTH_DZ*ZFLXZ, D%NKA, D%NKU, D%NKL), X_LES_RES_ddxa_Thl_SBG_UaRt )
    CALL LES_MEAN_SUBGRID(MZF(PDR_DZ*ZFLXZ, D%NKA, D%NKU, D%NKL), X_LES_RES_ddxa_Rt_SBG_UaRt )
    CALL LES_MEAN_SUBGRID(MZF(MZM(PEMOIST, D%NKA, D%NKU, D%NKL)*ZFLXZ, D%NKA, D%NKU, D%NKL), X_LES_SUBGRID_WThv , .TRUE. ) 
    CALL LES_MEAN_SUBGRID(-CSTURB%XCTP*PSQRT_TKE/PLM*MZF(ZFLXZ, D%NKA, D%NKU, D%NKL), X_LES_SUBGRID_RtPz ) 
    CALL SECOND_MNH(ZTIME2)
    XTIME_LES = XTIME_LES + ZTIME2 - ZTIME1
  END IF
!
END IF
!
!----------------------------------------------------------------------------
!
!
!*       4.   TURBULENT CORRELATIONS : <w Rc>
!             -------------------------------
!
!
!*       4.1  <w Rc>    
!
IF ( ((OTURB_FLX .AND. TPFILE%LOPENED) .OR. OLES_CALL) .AND. (KRRL > 0) ) THEN
!  
! recover the Conservative potential temperature flux : 
! With OHARAT is true tke and length scales at half levels
! yet modify to use length scale and tke at half levels from vdfexcuhl
  DO JK=1,D%NKT 
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE 
        ZWORK1(JI,JJ,JK) = PIMPL * PTHLP(JI,JJ,JK) + PEXPL * PTHLM(JI,JJ,JK)
      ENDDO
    ENDDO
  ENDDO
 CALL DZM_PHY(D,ZWORK1,ZWORK2)
 IF (OHARAT) THEN
    DO JK=1,D%NKT 
      DO JJ=IJB,IJE 
        DO JI=IIB,IIE 
          ZA(JI,JJ,JK)   = ZWORK2(JI,JJ,JK)/ PDZZ(JI,JJ,JK) * &
          (-PLM(JI,JJ,JK)*PSQRT_TKE(JI,JJ,JK))
        ENDDO
      ENDDO
    ENDDO
 ELSE
    DO JK=1,D%NKT 
      DO JJ=IJB,IJE 
        DO JI=IIB,IIE 
          ZWORK1(JI,JJ,JK) = PLM(JI,JJ,JK)*PSQRT_TKE(JI,JJ,JK)
        ENDDO
      ENDDO
    ENDDO
  CALL MZM_PHY(D,ZWORK1,ZWORK3)
    DO JK=1,D%NKT 
      DO JJ=IJB,IJE 
        DO JI=IIB,IIE 
          ZA(JI,JJ,JK)   = ZWORK2(JI,JJ,JK)/ PDZZ(JI,JJ,JK) * &
          (-PPHI3(JI,JJ,JK)*ZWORK3(JI,JJ,JK)) * CSTURB%XCSHF 
        ENDDO
      ENDDO
    ENDDO
 ENDIF
  DO JJ=IJB,IJE 
    DO JI=IIB,IIE 
      ZA(JI,JJ,IKB) = (PIMPL*PSFTHP(JI,JJ) + PEXPL*PSFTHM(JI,JJ)) * PDIRCOSZW(JI,JJ)
    ENDDO
  ENDDO
  !  
  ! compute <w Rc>
  DO JK=1,D%NKT 
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE 
        ZWORK1(JI,JJ,JK) = PAMOIST(JI,JJ,JK) * 2.* PSRCM(JI,JJ,JK)
        ZWORK2(JI,JJ,JK) = PATHETA(JI,JJ,JK) * 2.* PSRCM(JI,JJ,JK)
      ENDDO
    ENDDO
  ENDDO
  CALL MZM_PHY(D,ZWORK1,ZWORK3)
  CALL MZM_PHY(D,ZWORK2,ZWORK4)
  DO JK=1,D%NKT 
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE 
        ZFLXZ(JI,JJ,JK) = ZWORK3(JI,JJ,JK)* ZFLXZ(JI,JJ,JK) &
        + ZWORK4(JI,JJ,JK)* ZA(JI,JJ,JK)
      ENDDO
    ENDDO
  ENDDO
  DO JJ=IJB,IJE 
    DO JI=IIB,IIE 
      ZFLXZ(JI,JJ,D%NKA) = ZFLXZ(JI,JJ,IKB)
    ENDDO
  ENDDO
  !                 
  ! store the liquid water mixing ratio vertical flux
  IF ( OTURB_FLX .AND. TPFILE%LOPENED ) THEN
    TZFIELD%CMNHNAME   = 'RCW_FLX'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'RCW_FLX'
    TZFIELD%CUNITS     = 'kg m s-1 kg-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'Liquid water mixing ratio vertical flux'
    TZFIELD%NGRID      = 4
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZFLXZ)
  END IF
  !  
! and we store in LES configuration this subgrid flux <w'rc'>
!
  IF (OLES_CALL) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL LES_MEAN_SUBGRID( MZF(ZFLXZ, D%NKA, D%NKU, D%NKL), X_LES_SUBGRID_WRc ) 
    CALL SECOND_MNH(ZTIME2)
    XTIME_LES = XTIME_LES + ZTIME2 - ZTIME1
  END IF
!
END IF !end of <w Rc>
!
!----------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TURB_VER_THERMO_FLUX',1,ZHOOK_HANDLE)
END SUBROUTINE TURB_VER_THERMO_FLUX
END MODULE MODE_TURB_VER_THERMO_FLUX
