!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
MODULE MODE_TURB_VER_DYN_FLUX
IMPLICIT NONE
CONTAINS
SUBROUTINE TURB_VER_DYN_FLUX(D,CST,CSTURB,TURBN,KSV,O2D,OFLAT,&
                      OTURB_FLX,KRR, OOCEAN,OHARAT,OCOUPLES,OLES_CALL,&
                      HTURBDIM,PIMPL,PEXPL,                         &
                      PTSTEP,                                       &
                      TPFILE,                                       &
                      PDXX,PDYY,PDZZ,PDZX,PDZY,PDIRCOSZW,PZZ,       &
                      PCOSSLOPE,PSINSLOPE,                          &
                      PRHODJ,                                       &
                      PCDUEFF,PTAU11M,PTAU12M,PTAU33M,              &
                      PTHLM,PRM,PSVM,PUM,PVM,PWM,PUSLOPEM,PVSLOPEM, &
                      PTKEM,PLM,MFMOIST,PWU,PWV,                    &
                      PRUS,PRVS,PRWS,                               &
                      PDP,PTP                                       )
!     ###############################################################
!
!
!!****  *TURB_VER_DYN_FLUX* -compute the source terms due to the vertical turbulent
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
!!      Modifications: Mar  31, 1998 (V. Masson) splits the routine TURB_VER_DYN_FLUX 
!!      Modifications: Oct  18, 2000 (J. Stein)  Bug in some computations for IKB level
!!      Modifications: Oct  18, 2000 (V. Masson) LES computations + LFLAT switch
!!                     Nov  06, 2002 (V. Masson) LES budgets
!!                     October 2009 (G. Tanguy) add ILENCH=LEN(YCOMMENT) after
!!                                              change of YCOMMENT
!!      2012-02 Y. Seity,  add possibility to run with reversed vertical levels
!!      Modifications  July 2015 (Wim de Rooy) OHARATU switch
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!      Q. Rodier      17/01/2019 : cleaning : remove cyclic conditions on DP and ZA
!! JL Redelsperger 03/2021 : Add Ocean  & O-A Autocoupling LES Cases
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
USE MODD_IO,             ONLY: TFILEDATA
USE MODD_LES
USE MODD_OCEANH, ONLY: XSSUFL, XSSUFL_T,XSSVFL
USE MODD_PARAMETERS, ONLY: JPVEXT_TURB,XUNDEF
USE MODD_TURB_n, ONLY: TURB_t
!
USE SHUMAN_PHY
USE MODE_GRADIENT_U_PHY, ONLY : GZ_U_UW_PHY
USE MODE_GRADIENT_V_PHY, ONLY : GZ_V_VW_PHY
!
USE MODI_GRADIENT_U
USE MODI_GRADIENT_V
USE MODI_GRADIENT_W
USE MODI_GRADIENT_M
USE MODI_SECOND_MNH
USE MODI_SHUMAN , ONLY: MZM, MZF, MXM, MXF, MYM, MYF,&
                      & DZM, DXF, DXM, DYF, DYM
USE MODE_TRIDIAG_WIND, ONLY: TRIDIAG_WIND
USE MODI_LES_MEAN_SUBGRID
!
USE MODE_IO_FIELD_WRITE, only: IO_Field_write
USE MODE_ll
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
INTEGER,                INTENT(IN)   :: KSV           ! number of scalar variables
LOGICAL,                INTENT(IN)   ::  OTURB_FLX    ! switch to write the
                                 ! turbulent fluxes in the syncronous FM-file
LOGICAL,                INTENT(IN)   ::  OOCEAN       ! switch for Ocean model version
LOGICAL,                INTENT(IN)   ::  OHARAT
LOGICAL,                INTENT(IN)   ::  O2D          ! Logical for 2D model version (modd_conf)
LOGICAL,                INTENT(IN)   ::  OFLAT        ! Logical for zero ororography
LOGICAL,                INTENT(IN)   ::  OLES_CALL    ! compute the LES diagnostics at current time-step
LOGICAL,                INTENT(IN)   ::  OCOUPLES     ! switch to activate atmos-ocean LES version 
INTEGER,                INTENT(IN)   ::  KRR          ! number of moist var.
CHARACTER(len=4),       INTENT(IN)   ::  HTURBDIM     ! dimensionality of the
                                                      ! turbulence scheme
REAL,                   INTENT(IN)   ::  PIMPL, PEXPL ! Coef. for temporal disc.
REAL,                   INTENT(IN)   ::  PTSTEP       ! Double Time Step
TYPE(TFILEDATA),        INTENT(IN)   ::  TPFILE       ! Output file
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PDXX, PDYY, PDZZ, PDZX, PDZY 
                                                      ! Metric coefficients
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)   ::  PDIRCOSZW    ! Director Cosinus of the
                                                      ! normal to the ground surface
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PZZ          ! altitude of flux points
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)   ::  PCOSSLOPE    ! cosinus of the angle 
                                      ! between i and the slope vector
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)   ::  PSINSLOPE    ! sinus of the angle 
                                      ! between i and the slope vector
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PRHODJ       ! dry density * grid volum
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  MFMOIST      ! moist mass flux dual scheme

!
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)   ::  PCDUEFF     ! Cd * || u || at time t
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)   ::  PTAU11M      ! <uu> in the axes linked 
       ! to the maximum slope direction and the surface normal and the binormal 
       ! at time t - dt
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)   ::  PTAU12M      ! <uv> in the same axes
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)   ::  PTAU33M      ! <ww> in the same axes
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PUM,PVM,PWM, PTHLM
  ! Wind at t-Delta t
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KRR), INTENT(IN) ::  PRM
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KSV), INTENT(IN) ::  PSVM
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)   ::  PUSLOPEM     ! wind component along the 
                                     ! maximum slope direction
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)   ::  PVSLOPEM     ! wind component along the 
                                     ! direction normal to the maximum slope one
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PTKEM        ! TKE at time t
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PLM          ! Turb. mixing length   
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT)  ::  PWU          ! momentum flux u'w'
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT)  ::  PWV          ! momentum flux v'w'
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(INOUT)   ::  PRUS, PRVS, PRWS
                            ! cumulated sources for the prognostic variables
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT)  ::  PDP          ! Dynamic TKE production term
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PTP          ! Thermal TKE production term
!
!
!
!
!*       0.2  declaration of local variables
!
!
REAL, DIMENSION(D%NIT,D%NJT)  :: ZDIRSINZW ! sinus of the angle
                   ! between the normal and the vertical at the surface
REAL, DIMENSION(D%NIT,D%NJT,1):: ZCOEFS    ! coeff. for the 
                   ! implicit scheme for the wind at the surface
REAL, DIMENSION(D%NIT,D%NJT,D%NKT)  ::  &
       ZA, &       ! under diagonal elements of the tri-diagonal matrix involved
                   ! in the temporal implicit scheme (also used to store coefficient
                   ! J in Section 5)
       ZRES, &     ! guess of the treated variable at t+ deltat when the turbu-
                   ! lence is the only source of evolution added to the ones
                   ! considered in ZSOURCE  
       ZFLXZ,  &   ! vertical flux of the treated variable
       ZSOURCE,  & ! source of evolution for the treated variable
       ZKEFF,    & ! effectif diffusion coeff = LT * SQRT( TKE )
       ZWORK1,ZWORK2,&
       ZWORK3,ZWORK4,&
       ZWORK5,ZWORK6! working var. for shuman operators (array syntax)
!
INTEGER             :: IIE,IIB,IJE,IJB,IKB,IKE      ! index value for the mass points of the domain 
INTEGER             :: IKT          ! array size in k direction
INTEGER             :: IKTB,IKTE    ! start, end of k loops in physical domain
INTEGER             :: JSV,JI,JJ,JK          ! scalar loop counter
REAL, DIMENSION(D%NIT,D%NJT,1)   :: ZCOEFFLXU, &
                                    ZCOEFFLXV, ZUSLOPEM, ZVSLOPEM
                                    ! coefficients for the surface flux
                                    ! evaluation and copy of PUSLOPEM and
                                    ! PVSLOPEM in local 3D arrays
!
REAL :: ZTIME1, ZTIME2, ZCMFS
TYPE(TFIELDDATA) :: TZFIELD
!----------------------------------------------------------------------------
!
!*       1.   PRELIMINARIES
!             -------------
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('TURB_VER_DYN_FLUX',0,ZHOOK_HANDLE)
!
ZA(:,:,:)=XUNDEF
PDP(:,:,:)=XUNDEF
!
IKT=D%NKT  
IKTB=D%NKTB          
IKTE=D%NKTE
IKB=D%NKB
IKE=D%NKE
IIE=D%NIEC
IIB=D%NIBC
IJE=D%NJEC
IJB=D%NJBC
!
ZSOURCE(:,:,:) = 0.
ZFLXZ(:,:,:) = 0.
ZCMFS = CSTURB%XCMFS
IF (OHARAT) ZCMFS=1.
!
DO JJ=IJB,IJE 
  DO JI=IIB,IIE 
    ZDIRSINZW(JI,JJ) = SQRT(1.-PDIRCOSZW(JI,JJ)**2)
  ENDDO
ENDDO
!  compute the coefficients for the uncentred gradient computation near the 
!  ground
!
! With OHARATU length scale and TKE are at half levels so remove MZM
!
IF (OHARAT) THEN
  DO JK=1,D%NKT 
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE 
        ZKEFF(JI,JJ,JK) =  PLM(JI,JJ,JK) * SQRT(PTKEM(JI,JJ,JK)) + & 
        50*MFMOIST(JI,JJ,JK)
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
ZUSLOPEM(IIB:IIE,IJB:IJE,1)=PUSLOPEM(IIB:IIE,IJB:IJE)
ZVSLOPEM(IIB:IIE,IJB:IJE,1)=PVSLOPEM(IIB:IIE,IJB:IJE)
!
!----------------------------------------------------------------------------
!
!
!*       5.   SOURCES OF U,W WIND COMPONENTS AND PARTIAL DYNAMIC PRODUCTION 
!             -------------------------------------------------------------
!
!*       5.1  Source of U wind component
!
! Preparation of the arguments for TRIDIAG_WIND 
!
CALL MXM_PHY(D,ZKEFF,ZWORK1)
CALL MXM_PHY(D,PDZZ,ZWORK2)
CALL MZM_PHY(D,PRHODJ,ZWORK3)
CALL MXM_PHY(D,ZWORK3,ZWORK4)
DO JK=1,D%NKT 
  DO JJ=IJB,IJE 
    DO JI=IIB,IIE 
      ZA(JI,JJ,JK) = -PTSTEP * ZCMFS * ZWORK1(JI,JJ,JK)* ZWORK4(JI,JJ,JK) &
      / ZWORK2(JI,JJ,JK)**2
    ENDDO
  ENDDO
ENDDO
!
!
! Compute the source of U wind component 
!
! compute the coefficient between the vertical flux and the 2 components of the 
! wind following the slope
DO JJ=IJB,IJE 
  DO JI=IIB,IIE 
    ZCOEFFLXU(JI,JJ,1) = PCDUEFF(JI,JJ) * (PDIRCOSZW(JI,JJ)**2 - ZDIRSINZW(JI,JJ)**2) &
    * PCOSSLOPE(JI,JJ)
    ZCOEFFLXV(JI,JJ,1) = PCDUEFF(JI,JJ) * PDIRCOSZW(JI,JJ) * PSINSLOPE(JI,JJ)
    
! prepare the implicit scheme coefficients for the surface flux
    ZCOEFS(JI,JJ,1)=  ZCOEFFLXU(JI,JJ,1) * PCOSSLOPE(JI,JJ) * PDIRCOSZW(JI,JJ)  &
    +ZCOEFFLXV(JI,JJ,1) * PSINSLOPE(JI,JJ)
!
! average this flux to be located at the U,W vorticity point
  ENDDO
ENDDO
ZWORK1(IIB:IIE,IJB:IJE,1:1)=ZCOEFS(IIB:IIE,IJB:IJE,1:1) / PDZZ(IIB:IIE,IJB:IJE,IKB:IKB) 
CALL MXM_PHY(D,ZWORK1(IIB:IIE,IJB:IJE,1:1),ZCOEFS(IIB:IIE,IJB:IJE,1:1))
!
!
! ZSOURCE= FLUX /DZ
CALL MXM_PHY(D,PRHODJ,ZWORK1)
IF (OOCEAN) THEN  ! OCEAN MODEL ONLY
  ! Sfx flux assumed to be in SI & at vorticity point
  DO JJ=IJB,IJE 
    DO JI=IIB,IIE 
      IF (OCOUPLES) THEN
        ZSOURCE(JI,JJ,IKE) = TURBN%XSSUFL_C(JI,JJ,1)/PDZZ(JI,JJ,IKE) &
        *0.5 * ( 1. + ZWORK1(JI,JJ,D%NKU) / ZWORK1(JI,JJ,IKE)) 
      ELSE
        ZSOURCE(JI,JJ,IKE)     = XSSUFL(JI,JJ)
        ZSOURCE(JI,JJ,IKE) = ZSOURCE(JI,JJ,IKE) /PDZZ(JI,JJ,IKE) &
        *0.5 * ( 1. + ZWORK1(JI,JJ,D%NKU) / ZWORK1(JI,JJ,IKE) )
      ENDIF
    ENDDO
  ENDDO
  !No flux at the ocean domain bottom
  ZSOURCE(IIB:IIE,IJB:IJE,IKB)           = 0.
  ZSOURCE(IIB:IIE,IJB:IJE,IKTB+1:IKTE-1) = 0
!
ELSE             !ATMOS MODEL ONLY
  IF (OCOUPLES) THEN
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE 
        ZSOURCE(JI,JJ,IKB) = TURBN%XSSUFL_C(JI,JJ,1)/PDZZ(JI,JJ,IKB) &
        * 0.5 * ( 1. + ZWORK1(JI,JJ,D%NKA) / ZWORK1(JI,JJ,IKB) )
      ENDDO
    ENDDO
  ELSE
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE                
    ! compute the explicit tangential flux at the W point
        ZSOURCE(JI,JJ,IKB)     =                                              &
        PTAU11M(JI,JJ) * PCOSSLOPE(JI,JJ) * PDIRCOSZW(JI,JJ) * ZDIRSINZW(JI,JJ) &
        -PTAU12M(JI,JJ) * PSINSLOPE(JI,JJ) * ZDIRSINZW(JI,JJ)                  &
        -PTAU33M(JI,JJ) * PCOSSLOPE(JI,JJ) * ZDIRSINZW(JI,JJ) * PDIRCOSZW(JI,JJ)  
      ENDDO
    ENDDO
!
    ! add the vertical part or the surface flux at the U,W vorticity point
!
    DO JK=IKB,IKB 
      DO JJ=IJB,IJE 
        DO JI=IIB,IIE 
          ZWORK3(JI,JJ,JK) = ZSOURCE(JI,JJ,JK)/PDZZ(JI,JJ,JK)
        ENDDO
      ENDDO
    ENDDO
    CALL MXM_PHY(D,ZWORK3,ZWORK4)
    ZWORK5(IIB:IIE,IJB:IJE,IKB:IKB)= ZCOEFFLXU(IIB:IIE,IJB:IJE,1:1) / PDZZ(IIB:IIE,IJB:IJE,IKB:IKB)       &
           *ZUSLOPEM(IIB:IIE,IJB:IJE,1:1)                           &
          -ZCOEFFLXV(IIB:IIE,IJB:IJE,1:1) / PDZZ(IIB:IIE,IJB:IJE,IKB:IKB)       &
           *ZVSLOPEM(IIB:IIE,IJB:IJE,1:1)
    CALL MXM_PHY(D,ZWORK5,ZWORK6)
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE 
        ZSOURCE(JI,JJ,IKB) =                                  &
        (   ZWORK4(JI,JJ,IKB) &
        +  ZWORK6(JI,JJ,IKB)   &
        -  ZCOEFS(JI,JJ,1) * PUM(JI,JJ,IKB) * PIMPL        &
        ) * 0.5 * ( 1. + ZWORK1(JI,JJ,D%NKA) / ZWORK1(JI,JJ,IKB) )
      ENDDO
    ENDDO
  ENDIF 
!
  ZSOURCE(IIB:IIE,IJB:IJE,IKTB+1:IKTE-1) = 0.
  ZSOURCE(IIB:IIE,IJB:IJE,IKE) = 0.
ENDIF !end ocean or atmosphere cases
!
! Obtention of the split U at t+ deltat 
!
CALL TRIDIAG_WIND(D,PUM,ZA,ZCOEFS(:,:,1),PTSTEP,PEXPL,PIMPL,   &
                  ZWORK1,ZSOURCE,ZRES)
! 
!  Compute the equivalent tendency for the U wind component
!
CALL MXM_PHY(D,PRHODJ,ZWORK1)
CALL MXM_PHY(D,ZKEFF,ZWORK2)
DO JK=1,D%NKT 
  DO JJ=IJB,IJE 
    DO JI=IIB,IIE 
      ZWORK3(JI,JJ,JK)=PIMPL*ZRES(JI,JJ,JK) + PEXPL*PUM(JI,JJ,JK)
    ENDDO
  ENDDO
ENDDO
CALL DZM_PHY(D,ZWORK3,ZWORK4)
CALL MXM_PHY(D,PDZZ,ZWORK5)
DO JK=1,D%NKT 
  DO JJ=IJB,IJE 
    DO JI=IIB,IIE 
      PRUS(JI,JJ,JK)= PRUS(JI,JJ,JK)+ZWORK1(JI,JJ,JK)*(ZRES(JI,JJ,JK) & 
      - PUM(JI,JJ,JK))/PTSTEP
!
!
!*       5.2  Partial Dynamic Production
!
! vertical flux of the U wind component
!
      ZFLXZ(JI,JJ,JK)     = -ZCMFS * ZWORK2(JI,JJ,JK) * ZWORK4(JI,JJ,JK) &
      / ZWORK5(JI,JJ,JK)
    ENDDO
  ENDDO
ENDDO
!
! surface flux
CALL MXM_PHY(D,PDZZ,ZWORK1)
CALL MXM_PHY(D,PRHODJ,ZWORK2)
DO JJ=IJB,IJE 
  DO JI=IIB,IIE 
    ZFLXZ(JI,JJ,IKB)   =   ZWORK1(JI,JJ,IKB)  *                &
    ( ZSOURCE(JI,JJ,IKB)                                          &
    +ZCOEFS(JI,JJ,1) * ZRES(JI,JJ,IKB) * PIMPL                   &                
    ) / 0.5 / ( 1. + ZWORK2(JI,JJ,D%NKA)/ ZWORK2(JI,JJ,IKB) )
!
    ZFLXZ(JI,JJ,D%NKA) = ZFLXZ(JI,JJ,IKB) 
  ENDDO
ENDDO
!
IF (OOCEAN) THEN !ocean model at phys sfc (ocean domain top)
!
  DO JJ=IJB,IJE 
    DO JI=IIB,IIE 
      ZFLXZ(JI,JJ,IKE)   =  ZWORK1(JI,JJ,IKE) *                &
      ZSOURCE(JI,JJ,IKE)                     &
      / 0.5 / ( 1. + ZWORK2(JI,JJ,D%NKU)/ ZWORK2(JI,JJ,IKE) )
      ZFLXZ(JI,JJ,D%NKU) = ZFLXZ(JI,JJ,IKE)
    ENDDO
  ENDDO
END IF
!
IF ( OTURB_FLX .AND. TPFILE%LOPENED ) THEN
  ! stores the U wind component vertical flux
  TZFIELD%CMNHNAME   = 'UW_VFLX'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'UW_VFLX'
  TZFIELD%CUNITS     = 'm2 s-2'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'U wind component vertical flux'
  TZFIELD%NGRID      = 4
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZFLXZ)
END IF
!
! first part of total momentum flux
!
PWU(IIB:IIE,IJB:IJE,1:D%NKT) = ZFLXZ(IIB:IIE,IJB:IJE,1:D%NKT)
!
! Contribution to the dynamic production of TKE
! compute the dynamic production at the mass point
!
CALL GZ_U_UW_PHY(D,PUM,PDZZ,ZWORK1)
DO JK=1,D%NKT 
  DO JJ=IJB,IJE 
    DO JI=IIB,IIE 
      ZWORK2(JI,JJ,JK) = ZFLXZ(JI,JJ,JK) * ZWORK1(JI,JJ,JK)
    ENDDO
  ENDDO
ENDDO
CALL MXF_PHY(D,ZWORK2,ZWORK3)
CALL MZF_PHY(D,ZWORK3,ZWORK4)
DO JK=1,D%NKT 
  DO JJ=IJB,IJE 
    DO JI=IIB,IIE 
      PDP(JI,JJ,JK) = -ZWORK4(JI,JJ,JK)
    ENDDO
  ENDDO
ENDDO
!
! evaluate the dynamic production at w(IKB+D%NKL) in PDP(IKB)
CALL MXM_PHY(D,PDZZ,ZWORK1)
ZWORK2(IIB:IIE,IJB:IJE,IKB) = ZFLXZ(IIB:IIE,IJB:IJE,IKB+D%NKL) * (PUM(IIB:IIE,IJB:IJE,IKB+D%NKL)-PUM(IIB:IIE,IJB:IJE,IKB))  &
                         / ZWORK1(IIB:IIE,IJB:IJE,IKB+D%NKL)
CALL MXF_PHY(D,ZWORK2,ZWORK3)
DO JJ=IJB,IJE 
  DO JI=IIB,IIE 
    PDP(JI,JJ,IKB) = -ZWORK3(JI,JJ,IKB)
  ENDDO
ENDDO
!
IF (OOCEAN) THEN
  ZWORK2(IIB:IIE,IJB:IJE,IKE) = ZFLXZ(IIB:IIE,IJB:IJE,IKE-D%NKL) * (PUM(IIB:IIE,IJB:IJE,IKE)-PUM(IIB:IIE,IJB:IJE,IKE-D%NKL))  &
                         / ZWORK1(IIB:IIE,IJB:IJE,IKE-D%NKL)
  CALL MXF_PHY(D,ZWORK2,ZWORK3)
  ! evaluate the dynamic production at w(IKE-D%NKL) in PDP(IKE)
  DO JJ=IJB,IJE 
    DO JI=IIB,IIE 
      PDP(JI,JJ,IKE) = -ZWORK3(JI,JJ,IKE)
    ENDDO
  ENDDO
END IF
!
! Storage in the LES configuration
! 
IF (OLES_CALL) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL LES_MEAN_SUBGRID(MZF(MXF(ZFLXZ), D%NKA, D%NKU, D%NKL), X_LES_SUBGRID_WU ) 
  CALL LES_MEAN_SUBGRID(MZF(MXF(GZ_U_UW(PUM,PDZZ, D%NKA, D%NKU, D%NKL) &
                       & *ZFLXZ), D%NKA, D%NKU, D%NKL), X_LES_RES_ddxa_U_SBG_UaU )
  CALL LES_MEAN_SUBGRID( ZCMFS * ZKEFF, X_LES_SUBGRID_Km )
  CALL SECOND_MNH(ZTIME2)
  XTIME_LES = XTIME_LES + ZTIME2 - ZTIME1
END IF
!
!*       5.3  Source of W wind component
!
!
IF(HTURBDIM=='3DIM') THEN
  ! Compute the source for the W wind component
                ! used to compute the W source at the ground
  DO JJ=1,D%NJT 
    DO JI=1,D%NIT 
      ZFLXZ(JI,JJ,D%NKA) = 2 * ZFLXZ(JI,JJ,IKB) - ZFLXZ(JI,JJ,IKB+D%NKL) ! extrapolation 
    ENDDO
  ENDDO
 IF (OOCEAN) THEN
    DO JJ=1,D%NJT 
      DO JI=1,D%NIT 
        ZFLXZ(JI,JJ,D%NKU) = 2 * ZFLXZ(JI,JJ,IKE) - ZFLXZ(JI,JJ,IKE-D%NKL) ! extrapolation
      ENDDO
    ENDDO
 END IF     
     
  !
  IF (.NOT. OFLAT) THEN
    ZWORK1 = DXF( MZM(MXM(PRHODJ) /PDXX, D%NKA, D%NKU, D%NKL)  * ZFLXZ )
    ZWORK2 = DZM(PRHODJ / MZF(PDZZ, D%NKA, D%NKU, D%NKL) *                &
                      MXF(MZF(ZFLXZ*PDZX, D%NKA, D%NKU, D%NKL) / PDXX ),      &
                     D%NKA, D%NKU, D%NKL)
    DO JK=1,D%NKT 
      DO JJ=1,D%NJT 
        DO JI=1,D%NIT 
          PRWS(JI,JJ,JK)= PRWS(JI,JJ,JK) - ZWORK1(JI,JJ,JK) + ZWORK2(JI,JJ,JK)
        ENDDO
      ENDDO
    ENDDO
  ELSE
    ZWORK1 = DXF(MZM(MXM(PRHODJ) /PDXX, D%NKA, D%NKU, D%NKL)  * ZFLXZ )
    DO JK=1,D%NKT 
      DO JJ=1,D%NJT 
        DO JI=1,D%NIT 
          PRWS(JI,JJ,JK)= PRWS(JI,JJ,JK) - ZWORK1(JI,JJ,JK)
        ENDDO
      ENDDO
    ENDDO
  END IF
  !
  ! Complete the Dynamical production with the W wind component 
  !
  ZA(:,:,:)=-MZF(MXF(ZFLXZ * GX_W_UW(PWM,PDXX,PDZZ,PDZX, D%NKA, D%NKU, D%NKL)), D%NKA, D%NKU, D%NKL)
  !
  !
  ! evaluate the dynamic production at w(IKB+D%NKL) in PDP(IKB)
  ZA(:,:,IKB:IKB) = - MXF (                                                  &
   ZFLXZ(:,:,IKB+D%NKL:IKB+D%NKL) *                                              &
     ( DXM( PWM(:,:,IKB+D%NKL:IKB+D%NKL) )                                       &
      -MXM(  (PWM(:,:,IKB+2*D%NKL:IKB+2*D%NKL   )-PWM(:,:,IKB+D%NKL:IKB+D%NKL))      &
              /(PDZZ(:,:,IKB+2*D%NKL:IKB+2*D%NKL)+PDZZ(:,:,IKB+D%NKL:IKB+D%NKL))     &
            +(PWM(:,:,IKB+D%NKL:IKB+D%NKL)-PWM(:,:,IKB:IKB  ))                   &
              /(PDZZ(:,:,IKB+D%NKL:IKB+D%NKL)+PDZZ(:,:,IKB:IKB  ))               &
          )                                                                  &
        * PDZX(:,:,IKB+D%NKL:IKB+D%NKL)                                          &
     ) / (0.5*(PDXX(:,:,IKB+D%NKL:IKB+D%NKL)+PDXX(:,:,IKB:IKB)))                 &
                          )
  !
IF (OOCEAN) THEN
  ! evaluate the dynamic production at w(IKE-D%NKL) in PDP(IKE)
  ZA(:,:,IKE:IKE) = - MXF (                                                  &
   ZFLXZ(:,:,IKE-D%NKL:IKE-D%NKL) *                                              &
     ( DXM( PWM(:,:,IKE-D%NKL:IKE-D%NKL) )                                       &
      -MXM(  (PWM(:,:,IKE-2*D%NKL:IKE-2*D%NKL   )-PWM(:,:,IKE-D%NKL:IKE-D%NKL))      &
              /(PDZZ(:,:,IKE-2*D%NKL:IKE-2*D%NKL)+PDZZ(:,:,IKE-D%NKL:IKE-D%NKL))     &
            +(PWM(:,:,IKE-D%NKL:IKE-D%NKL)-PWM(:,:,IKE:IKE  ))                   &
              /(PDZZ(:,:,IKE-D%NKL:IKE-D%NKL)+PDZZ(:,:,IKE:IKE  ))               &
          )                                                                  &
         * PDZX(:,:,IKE-D%NKL:IKE-D%NKL)                                         &
     ) / (0.5*(PDXX(:,:,IKE-D%NKL:IKE-D%NKL)+PDXX(:,:,IKE:IKE)))                 &
                          )
END IF
  !
  DO JK=1,D%NKT 
    DO JJ=1,D%NJT 
      DO JI=1,D%NIT 
        PDP(JI,JJ,JK)=PDP(JI,JJ,JK)+ZA(JI,JJ,JK)
      ENDDO
    ENDDO
  ENDDO
  !
  ! Storage in the LES configuration
  ! 
  IF (OLES_CALL) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL LES_MEAN_SUBGRID(MZF(MXF(GX_W_UW(PWM,PDXX,&
      PDZZ,PDZX, D%NKA, D%NKU, D%NKL)*ZFLXZ), D%NKA, D%NKU, D%NKL), X_LES_RES_ddxa_W_SBG_UaW )
    CALL LES_MEAN_SUBGRID(MXF(GX_M_U(D%NKA, D%NKU, D%NKL,PTHLM,PDXX,PDZZ,PDZX)&
      * MZF(ZFLXZ, D%NKA, D%NKU, D%NKL)), X_LES_RES_ddxa_Thl_SBG_UaW )
    IF (KRR>=1) THEN
      CALL LES_MEAN_SUBGRID(MXF(GX_U_M(PRM(:,:,:,1),PDXX,PDZZ,PDZX, D%NKA, D%NKU, D%NKL)&
      *MZF(ZFLXZ, D%NKA, D%NKU, D%NKL)),X_LES_RES_ddxa_Rt_SBG_UaW )
    END IF
    DO JSV=1,KSV
      CALL LES_MEAN_SUBGRID( MXF(GX_U_M(PSVM(:,:,:,JSV),PDXX,PDZZ,&
      PDZX, D%NKA, D%NKU, D%NKL)*MZF(ZFLXZ, D%NKA, D%NKU, D%NKL)),X_LES_RES_ddxa_Sv_SBG_UaW(:,:,:,JSV) )
    END DO
    CALL SECOND_MNH(ZTIME2)
    XTIME_LES = XTIME_LES + ZTIME2 - ZTIME1
  END IF
END IF
!
!----------------------------------------------------------------------------
!
!
!*       6.   SOURCES OF V,W WIND COMPONENTS AND COMPLETE 1D DYNAMIC PRODUCTION 
!             -----------------------------------------------------------------
!
!*       6.1  Source of V wind component
!
! Preparation of the arguments for TRIDIAG_WIND
!!
CALL MYM_PHY(D,ZKEFF,ZWORK1)
CALL MYM_PHY(D,PDZZ,ZWORK2)
CALL MZM_PHY(D,PRHODJ,ZWORK3)
CALL MYM_PHY(D,ZWORK3,ZWORK4)
DO JK=1,D%NKT 
  DO JJ=IJB,IJE 
    DO JI=IIB,IIE 
      ZA(JI,JJ,JK) = -PTSTEP * ZCMFS * ZWORK1(JI,JJ,JK)* ZWORK4(JI,JJ,JK) & 
      / ZWORK2(JI,JJ,JK)**2
    ENDDO
  ENDDO
ENDDO
!
!
!
! Compute the source of V wind component
! compute the coefficient between the vertical flux and the 2 components of the 
! wind following the slope
DO JJ=IJB,IJE 
  DO JI=IIB,IIE 
    ZCOEFFLXU(JI,JJ,1) = PCDUEFF(JI,JJ) * (PDIRCOSZW(JI,JJ)**2 - ZDIRSINZW(JI,JJ)**2) &
    * PSINSLOPE(JI,JJ)
    ZCOEFFLXV(JI,JJ,1) = PCDUEFF(JI,JJ) * PDIRCOSZW(JI,JJ) * PCOSSLOPE(JI,JJ)
    
! prepare the implicit scheme coefficients for the surface flux
    ZCOEFS(JI,JJ,1)=  ZCOEFFLXU(JI,JJ,1) * PSINSLOPE(JI,JJ) * PDIRCOSZW(JI,JJ)  &
    +ZCOEFFLXV(JI,JJ,1) * PCOSSLOPE(JI,JJ)
  ENDDO
ENDDO
!
! average this flux to be located at the V,W vorticity point
ZWORK1(IIB:IIE,IJB:IJE,1:1)=ZCOEFS(IIB:IIE,IJB:IJE,1:1) / PDZZ(IIB:IIE,IJB:IJE,IKB:IKB) 
CALL MYM_PHY(D,ZWORK1(IIB:IIE,IJB:IJE,1:1),ZCOEFS(IIB:IIE,IJB:IJE,1:1))
!
CALL MYM_PHY(D,PRHODJ,ZWORK1) ! it was MXM(PRHODJ) : bug corrected now ? REPRO55 OCEAN
IF (OOCEAN) THEN ! Ocean case
  IF (OCOUPLES) THEN
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE 
        ZSOURCE(JI,JJ,IKE) =  TURBN%XSSVFL_C(JI,JJ,1)/PDZZ(JI,JJ,IKE) &
        *0.5 * ( 1. + ZWORK1(JI,JJ,D%NKU) / ZWORK1(JI,JJ,IKE))
      ENDDO
    ENDDO
  ELSE 
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE 
        ZSOURCE(JI,JJ,IKE) = XSSVFL(JI,JJ)
        ZSOURCE(JI,JJ,IKE) = ZSOURCE(JI,JJ,IKE)/PDZZ(JI,JJ,IKE) &
        *0.5 * ( 1. + ZWORK1(JI,JJ,D%NKU) / ZWORK1(JI,JJ,IKE))
      ENDDO
    ENDDO
  END IF
  !No flux at the ocean domain bottom
  ZSOURCE(IIB:IIE,IJB:IJE,IKB) = 0.
ELSE ! Atmos case
  ZWORK3(IIB:IIE,IJB:IJE,IKB:IKB) = ZCOEFFLXU(IIB:IIE,IJB:IJE,1:1) / PDZZ(IIB:IIE,IJB:IJE,IKB:IKB)           &
            *ZUSLOPEM(IIB:IIE,IJB:IJE,1:1)                                &
            +ZCOEFFLXV(IIB:IIE,IJB:IJE,1:1) / PDZZ(IIB:IIE,IJB:IJE,IKB:IKB)           &
            *ZVSLOPEM(IIB:IIE,IJB:IJE,1:1) 
    CALL MYM_PHY(D,ZWORK3,ZWORK6)
!
  IF (.NOT.OCOUPLES) THEN !  only atmosp without coupling
  ! compute the explicit tangential flux at the W point
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE 
        ZSOURCE(JI,JJ,IKB)       =                                                  &
        PTAU11M(JI,JJ) * PSINSLOPE(JI,JJ) * PDIRCOSZW(JI,JJ) * ZDIRSINZW(JI,JJ)         &
        +PTAU12M(JI,JJ) * PCOSSLOPE(JI,JJ) * ZDIRSINZW(JI,JJ)                          &
        -PTAU33M(JI,JJ) * PSINSLOPE(JI,JJ) * ZDIRSINZW(JI,JJ) * PDIRCOSZW(JI,JJ) 
      ENDDO
    ENDDO
    DO JK=IKB,IKB 
      DO JJ=IJB,IJE 
        DO JI=IIB,IIE 
          ZWORK3(JI,JJ,JK) = ZSOURCE(JI,JJ,JK)/PDZZ(JI,JJ,JK)
        ENDDO
      ENDDO
    ENDDO
    CALL MYM_PHY(D,ZWORK3,ZWORK5)
!
  ! add the vertical part or the surface flux at the V,W vorticity point
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE 
        ZSOURCE(JI,JJ,IKB) =                                      &
        (   ZWORK5(JI,JJ,IKB)     &
        +  ZWORK6(JI,JJ,IKB)        &
        - ZCOEFS(JI,JJ,1) * PVM(JI,JJ,IKB) * PIMPL             &
        ) * 0.5 * ( 1. + ZWORK1(JI,JJ,D%NKA) / ZWORK1(JI,JJ,IKB) )
      ENDDO
    ENDDO
!
  ELSE   !atmosphere when coupling
    ! input flux assumed to be in SI and at vorticity point
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE 
        ZSOURCE(JI,JJ,IKB) =     -TURBN%XSSVFL_C(JI,JJ,1)/(1.*PDZZ(JI,JJ,IKB)) &
        * 0.5 * ( 1. + ZWORK1(JI,JJ,D%NKA) / ZWORK1(JI,JJ,IKB)  )
      ENDDO
    ENDDO
  ENDIF
  !No flux at the atmosphere top
  ZSOURCE(IIB:IIE,IJB:IJE,IKE) = 0.
ENDIF ! End of Ocean or Atmospher Cases
ZSOURCE(IIB:IIE,IJB:IJE,IKTB+1:IKTE-1) = 0.
! 
!  Obtention of the split V at t+ deltat 
CALL TRIDIAG_WIND(D,PVM,ZA,ZCOEFS(:,:,1),PTSTEP,PEXPL,PIMPL,  &
                  ZWORK1,ZSOURCE,ZRES)
!
! Compute the equivalent tendency for the V wind component
!
CALL MYM_PHY(D,PRHODJ,ZWORK1)
CALL MYM_PHY(D,ZKEFF,ZWORK2)
DO JK=1,D%NKT 
  DO JJ=IJB,IJE 
    DO JI=IIB,IIE 
      ZWORK3(JI,JJ,JK)=PIMPL*ZRES(JI,JJ,JK) + PEXPL*PVM(JI,JJ,JK)
    ENDDO
  ENDDO
ENDDO
CALL DZM_PHY(D,ZWORK3,ZWORK4)
CALL MYM_PHY(D,PDZZ,ZWORK5)
DO JK=1,D%NKT 
  DO JJ=IJB,IJE 
    DO JI=IIB,IIE 
      PRVS(JI,JJ,JK) = PRVS(JI,JJ,JK)+ZWORK1(JI,JJ,JK)*(ZRES(JI,JJ,JK)& 
      - PVM(JI,JJ,JK))/PTSTEP
!
!
!*       6.2  Complete 1D dynamic Production
!
!  vertical flux of the V wind component
!
      ZFLXZ(JI,JJ,JK)   = -ZCMFS * ZWORK2(JI,JJ,JK) * ZWORK4(JI,JJ,JK) &
      / ZWORK5(JI,JJ,JK)
    ENDDO
  ENDDO
ENDDO
!
DO JJ=IJB,IJE 
  DO JI=IIB,IIE 
    ZFLXZ(JI,JJ,IKB)   =   ZWORK5(JI,JJ,IKB)  *                       &
    ( ZSOURCE(JI,JJ,IKB)                                                 &
    +ZCOEFS(JI,JJ,1) * ZRES(JI,JJ,IKB) * PIMPL                      &      
    ) / 0.5 / ( 1. + ZWORK1(JI,JJ,D%NKA) / ZWORK1(JI,JJ,IKB) )
!  
!
    ZFLXZ(JI,JJ,D%NKA) = ZFLXZ(JI,JJ,IKB)
  ENDDO
ENDDO
!
IF (OOCEAN) THEN
  DO JJ=IJB,IJE 
    DO JI=IIB,IIE 
      ZFLXZ(JI,JJ,IKE)   =   ZWORK5(JI,JJ,IKE)  *                &
      ZSOURCE(JI,JJ,IKE)                                          &
      / 0.5 / ( 1. + ZWORK1(JI,JJ,D%NKU) / ZWORK1(JI,JJ,IKE) )
      ZFLXZ(JI,JJ,D%NKU) = ZFLXZ(JI,JJ,IKE)
    ENDDO
  ENDDO
END IF
!
IF ( OTURB_FLX .AND. TPFILE%LOPENED ) THEN
  ! stores the V wind component vertical flux
  TZFIELD%CMNHNAME   = 'VW_VFLX'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'VW_VFLX'
  TZFIELD%CUNITS     = 'm2 s-2'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'V wind component vertical flux'
  TZFIELD%NGRID      = 4
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZFLXZ)
END IF
!
! second part of total momentum flux
!
PWV(IIB:IIE,IJB:IJE,1:D%NKT) = ZFLXZ(IIB:IIE,IJB:IJE,1:D%NKT)
!
!  Contribution to the dynamic production of TKE
! compute the dynamic production contribution at the mass point
!
CALL GZ_V_VW_PHY(D,PVM,PDZZ,ZWORK1)
DO JK=1,D%NKT 
  DO JJ=IJB,IJE 
    DO JI=IIB,IIE 
      ZWORK2(JI,JJ,JK) = ZFLXZ(JI,JJ,JK) * ZWORK1(JI,JJ,JK)
    ENDDO
  ENDDO
ENDDO
CALL MYF_PHY(D,ZWORK2,ZWORK3)
CALL MZF_PHY(D,ZWORK3,ZWORK4)
DO JK=1,D%NKT 
  DO JJ=IJB,IJE 
    DO JI=IIB,IIE 
      ZA(JI,JJ,JK) = -ZWORK4(JI,JJ,JK)
    ENDDO
  ENDDO
ENDDO
!
! evaluate the dynamic production at w(IKB+D%NKL) in PDP(IKB)
CALL MYM_PHY(D,PDZZ,ZWORK1)
ZWORK2(IIB:IIE,IJB:IJE,IKB) = ZFLXZ(IIB:IIE,IJB:IJE,IKB+D%NKL) * (PVM(IIB:IIE,IJB:IJE,IKB+D%NKL)-PVM(IIB:IIE,IJB:IJE,IKB))  &
                         / ZWORK1(IIB:IIE,IJB:IJE,IKB+D%NKL)
CALL MYF_PHY(D,ZWORK2,ZWORK3)
DO JJ=IJB,IJE 
  DO JI=IIB,IIE 
    ZA(JI,JJ,IKB) = -ZWORK3(JI,JJ,IKB)
  ENDDO
ENDDO
!
IF (OOCEAN) THEN
  ! evaluate the dynamic production at w(IKE-D%NKL) in PDP(IKE)
  ZWORK2(IIB:IIE,IJB:IJE,IKE) = ZFLXZ(IIB:IIE,IJB:IJE,IKE-D%NKL) * (PVM(IIB:IIE,IJB:IJE,IKE)-PVM(IIB:IIE,IJB:IJE,IKE-D%NKL))  &
                         / ZWORK1(IIB:IIE,IJB:IJE,IKE-D%NKL)
  CALL MYF_PHY(D,ZWORK2,ZWORK3)
  DO JJ=IJB,IJE 
    DO JI=IIB,IIE 
      ZA(JI,JJ,IKE) = -ZWORK3(JI,JJ,IKE)
    ENDDO
  ENDDO
END IF
!
DO JK=1,D%NKT 
  DO JJ=IJB,IJE 
    DO JI=IIB,IIE 
      PDP(JI,JJ,JK)=PDP(JI,JJ,JK)+ZA(JI,JJ,JK)
    ENDDO
  ENDDO
ENDDO
!
! Storage in the LES configuration
!
IF (OLES_CALL) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL LES_MEAN_SUBGRID(MZF(MYF(ZFLXZ), D%NKA, D%NKU, D%NKL), X_LES_SUBGRID_WV ) 
  CALL LES_MEAN_SUBGRID(MZF(MYF(GZ_V_VW(PVM,PDZZ, D%NKA, D%NKU, D%NKL)*&
                    & ZFLXZ), D%NKA, D%NKU, D%NKL), X_LES_RES_ddxa_V_SBG_UaV )
  CALL SECOND_MNH(ZTIME2)
  XTIME_LES = XTIME_LES + ZTIME2 - ZTIME1
END IF
!
!
!*       6.3  Source of W wind component 
!
IF(HTURBDIM=='3DIM') THEN
  ! Compute the source for the W wind component
  DO JJ=1,D%NJT 
    DO JI=1,D%NIT 
      ZFLXZ(JI,JJ,D%NKA) = 2 * ZFLXZ(JI,JJ,IKB) - ZFLXZ(JI,JJ,IKB+D%NKL) ! extrapolation
    ENDDO
  ENDDO
  IF (OOCEAN) THEN
    DO JJ=1,D%NJT 
      DO JI=1,D%NIT 
        ZFLXZ(JI,JJ,D%NKU) = 2 * ZFLXZ(JI,JJ,IKE) - ZFLXZ(JI,JJ,IKE-D%NKL) ! extrapolation 
      ENDDO
    ENDDO
  END IF
  !
  IF (.NOT. O2D) THEN
    ZWORK1 = DYF( MZM(MYM(PRHODJ) /PDYY, D%NKA, D%NKU, D%NKL) * ZFLXZ ) 
    IF (.NOT. OFLAT) THEN
      ZWORK2 = DZM(PRHODJ / MZF(PDZZ, D%NKA, D%NKU, D%NKL) *                &
                        MYF(MZF(ZFLXZ*PDZY, D%NKA, D%NKU, D%NKL) / PDYY ),      &
                       D%NKA, D%NKU, D%NKL)
      DO JK=1,D%NKT 
        DO JJ=1,D%NJT 
          DO JI=1,D%NIT 
            PRWS(JI,JJ,JK)= PRWS(JI,JJ,JK) - ZWORK1(JI,JJ,JK) + ZWORK2(JI,JJ,JK)
          ENDDO
        ENDDO
      ENDDO
    ELSE
      DO JK=1,D%NKT 
        DO JJ=1,D%NJT 
          DO JI=1,D%NIT 
            PRWS(JI,JJ,JK)= PRWS(JI,JJ,JK) - ZWORK1(JI,JJ,JK)
          ENDDO
        ENDDO
      ENDDO
    END IF
  END IF
  ! 
  ! Complete the Dynamical production with the W wind component 
  IF (.NOT. O2D) THEN
    ZA(:,:,:) = - MZF(MYF(ZFLXZ * GY_W_VW(PWM,PDYY,PDZZ,PDZY, D%NKA, D%NKU, D%NKL)), D%NKA, D%NKU, D%NKL)
  !
  ! evaluate the dynamic production at w(IKB+D%NKL) in PDP(IKB)
    ZA(:,:,IKB:IKB) = - MYF (                                              &
     ZFLXZ(:,:,IKB+D%NKL:IKB+D%NKL) *                                          &
       ( DYM( PWM(:,:,IKB+D%NKL:IKB+D%NKL) )                                   &
        -MYM(  (PWM(:,:,IKB+2*D%NKL:IKB+2*D%NKL)-PWM(:,:,IKB+D%NKL:IKB+D%NKL))     &
                /(PDZZ(:,:,IKB+2*D%NKL:IKB+2*D%NKL)+PDZZ(:,:,IKB+D%NKL:IKB+D%NKL)) &
              +(PWM(:,:,IKB+D%NKL:IKB+D%NKL)-PWM(:,:,IKB:IKB  ))               &
                /(PDZZ(:,:,IKB+D%NKL:IKB+D%NKL)+PDZZ(:,:,IKB:IKB  ))           &
            )                                                              &
          * PDZY(:,:,IKB+D%NKL:IKB+D%NKL)                                      &
       ) / (0.5*(PDYY(:,:,IKB+D%NKL:IKB+D%NKL)+PDYY(:,:,IKB:IKB)))             &
                            )
  !
    IF (OOCEAN) THEN
     ZA(:,:,IKE:IKE) = - MYF (                                              &
      ZFLXZ(:,:,IKE-D%NKL:IKE-D%NKL) *                                          &
        ( DYM( PWM(:,:,IKE-D%NKL:IKE-D%NKL) )                                   &
         -MYM(  (PWM(:,:,IKE-2*D%NKL:IKE-2*D%NKL)-PWM(:,:,IKE-D%NKL:IKE-D%NKL))     &
                 /(PDZZ(:,:,IKE-2*D%NKL:IKE-2*D%NKL)+PDZZ(:,:,IKE-D%NKL:IKE-D%NKL)) &
               +(PWM(:,:,IKE-D%NKL:IKE-D%NKL)-PWM(:,:,IKE:IKE  ))               &
                 /(PDZZ(:,:,IKE-D%NKL:IKE-D%NKL)+PDZZ(:,:,IKE:IKE  ))           &
             )                                                              &
           * PDZY(:,:,IKE-D%NKL:IKE-D%NKL)                                      &
        ) / (0.5*(PDYY(:,:,IKE-D%NKL:IKE-D%NKL)+PDYY(:,:,IKE:IKE)))             &
                            )
    END IF
!    
    DO JK=1,D%NKT 
      DO JJ=1,D%NJT 
        DO JI=1,D%NIT 
          PDP(JI,JJ,JK)=PDP(JI,JJ,JK)+ZA(JI,JJ,JK)
        ENDDO
      ENDDO
    ENDDO
  !
  END IF
  !
  ! Storage in the LES configuration
  !
  IF (OLES_CALL) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL LES_MEAN_SUBGRID(MZF(MYF(GY_W_VW(PWM,PDYY,&
                         &PDZZ,PDZY, D%NKA, D%NKU, D%NKL)*ZFLXZ), D%NKA, D%NKU, D%NKL), &
                         &X_LES_RES_ddxa_W_SBG_UaW , .TRUE. )
    CALL LES_MEAN_SUBGRID(MYF(GY_M_V(D%NKA, D%NKU, D%NKL,PTHLM,PDYY,PDZZ,PDZY)*&
                         &MZF(ZFLXZ, D%NKA, D%NKU, D%NKL)), &
                         &X_LES_RES_ddxa_Thl_SBG_UaW , .TRUE. )
    IF (KRR>=1) THEN
      CALL LES_MEAN_SUBGRID(MYF(GY_V_M(PRM(:,:,:,1),PDYY,PDZZ,&
                           &PDZY, D%NKA, D%NKU, D%NKL)*MZF(ZFLXZ, D%NKA, D%NKU, D%NKL)),&
                           &X_LES_RES_ddxa_Rt_SBG_UaW , .TRUE. )
    END IF
    CALL SECOND_MNH(ZTIME2)
    XTIME_LES = XTIME_LES + ZTIME2 - ZTIME1
  END IF
  !
END IF
!
!
!----------------------------------------------------------------------------
!
!*       7.   DIAGNOSTIC COMPUTATION OF THE 1D <W W> VARIANCE
!             -----------------------------------------------
!
IF ( OTURB_FLX .AND. TPFILE%LOPENED .AND. HTURBDIM == '1DIM') THEN
  ZWORK1 = GZ_W_M(PWM,PDZZ, D%NKA, D%NKU, D%NKL)
  DO JK=1,D%NKT 
    DO JJ=1,D%NJT 
      DO JI=1,D%NIT 
        ZFLXZ(JI,JJ,JK)= (2./3.) * PTKEM(JI,JJ,JK)                     &
        -ZCMFS*PLM(JI,JJ,JK)*SQRT(PTKEM(JI,JJ,JK))*ZWORK1(JI,JJ,JK)
      ENDDO
    ENDDO
  ENDDO
  ! to be tested &
  !   +XCMFB*(4./3.)*PLM(:,:,:)/SQRT(PTKEM(:,:,:))*PTP(:,:,:) 
  ! stores the W variance
  TZFIELD%CMNHNAME   = 'W_VVAR'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'W_VVAR'
  TZFIELD%CUNITS     = 'm2 s-2'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_W_VVAR'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZFLXZ)
END IF
!
!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('TURB_VER_DYN_FLUX',1,ZHOOK_HANDLE)
END SUBROUTINE TURB_VER_DYN_FLUX
END MODULE MODE_TURB_VER_DYN_FLUX
