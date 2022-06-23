!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
MODULE MODE_TRIDIAG
IMPLICIT NONE
CONTAINS
SUBROUTINE TRIDIAG(D,PVARM,PA,PTSTEP,PEXPL,PIMPL, &
                                  PRHODJ,PSOURCE,PVARP )
       USE PARKIND1, ONLY : JPRB
       USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!      #################################################
!
!
!!****   *TRIDIAG* - routine to solve a time implicit scheme
!!
!!
!!     PURPOSE
!!     -------
!        The purpose of this routine is to give a field PVARP at t+1, by 
!      solving an implicit tridiagonal system obtained by the 
!      discretization of the vertical turbulent diffusion. It should be noted 
!      that the degree of implicitness can be varied (PIMPL parameter) and the
!      sources of evolution other than the turbulent diffusion can be taken
!      into account through the PSOURCE field. PVARP is localized at a mass 
!      point.
!
!!**   METHOD
!!     ------
!!        First, the Right Hand Side of the implicit equation is computed.
!!     It is build as follows:
!!        ZY = PVARM + PTSTEP*PSOURCE + DIFF_EXPLI
!!     where PVARM is the variable at t-dt, PSOURCE the supplementary sources of
!!     PVAR ( and not PVAR * PRHODJ !!) and  DIFF_EXPLI is the explicit part
!!     of the vertical turbulent diffusion. This operator is spatially 
!!     discretized as the implicit one, thus:
!!        DIFF_EXPLI(k) = - PEXPL / PRHODJ(k) * 
!!                       ( PA(k+1) * (PVARM(k+1) - PVARM(k)  )
!!                        -PA(k)   * (PVARM(k)   - PVARM(k-1)) )
!!     For the first level, only the upper part is considered, the lower one 
!!     is replaced by the turbulent surface flux (taken into account in the 
!!     PSOURCE(ikb) term).
!!        DIFF_EXPLI(ikb) = - PEXPL / PRHODJ(ikb) * 
!!                       ( PA(ikb+1) * (PVARM(ikb+1) - PVARM(ikb))  )
!!     For the last level, only the lower part is considered, the upper one 
!!     is replaced by the turbulent flux which is taken equal to 0 
!!     (taken into account in the PSOURCE(ike) term).
!!
!!        DIFF_EXPLI(ike) = + PEXPL / PRHODJ(ike) * 
!!                       ( PA(ike) * (PVARM(ike) - PVARM(ike-1))  )
!!                      
!!        Then, the classical tridiagonal algorithm is used to invert the 
!!     implicit operator. Its matrix is given by:
!!
!!     ( b(ikb)   c(ikb)      0        0        0         0        0        0  )
!!     (   0      a(ikb+1) b(ikb+1) c(ikb+1)    0  ...    0        0        0  ) 
!!     (   0         0     a(ikb+2) b(ikb+2) c(ikb+2).    0        0        0  ) 
!!      .......................................................................
!!     (   0   ...   0     a(k)     b(k)     c(k)         0   ...  0        0  ) 
!!      .......................................................................
!!     (   0         0        0        0        0 ...a(ike-1) b(ike-1) c(ike-1))
!!     (   0         0        0        0        0 ...     0   a(ike)   b(ike)  )
!!
!!     ikb and ike represent the first and the last inner mass levels of the
!!     model. The coefficients are:
!!         
!!          a(k) = PIMPL * PA(k)/PRHODJ(k)   
!!          b(k) = 1 - PIMPL * PA(k)/PRHODJ(k) - PIMPL * PA(k+1)/PRHODJ(k)
!!          c(k) = PIMPL * PA(k+1)/PRHODJ(k)
!!
!!          for all k /= ikb or ike
!!
!!          b(ikb) = 1 - PIMPL * PA(ikb+1)/PRHODJ(ikb)
!!          c(ikb) = PIMPL * PA(ikb+1)/PRHODJ(ikb)
!!             (discretization of the upper part of the implicit operator)
!!          b(ike) = 1 - PIMPL * PA(ike)/PRHODJ(ike)
!!          a(ike) = PIMPL * PA(ike)/PRHODJ(ike)
!!             (discretization of the lower part of the implicit operator)
!!       Finally, the marginal points are prescribed.
!!
!!       All these computations are purely vertical and vectorizations are 
!!     easely achieved by processing all the verticals in parallel.
!!
!!     EXTERNAL
!!     --------
!!
!!       NONE
!!
!!     IMPLICIT ARGUMENTS
!!     ------------------
!!       MODD_PARAMETERS
!!            JPVEXT_TURB: number of vertical external points
!!
!!     REFERENCE
!!     ---------
!!       Book 1 of Meso-NH documentation (chapter Turbulence)
!!       Press et al: Numerical recipes (1986) Cambridge Univ. Press
!!
!!     AUTHOR
!!     ------
!!       Joan Cuxart       * INM and Meteo-France *   
!! 
!!     MODIFICATIONS
!!     -------------
!!       Original        August 29, 1994
!!       Modification :  January 29, 1995 Algorithm written with two 
!!                                        local variables less
!!      (Cuxart, Stein)  August 21, 1995  Bug correction for PRHODJ
!!      (Stein)         November 16, 1995 new version
!!      (Stein)         February 28, 1995 no inversion in the explicit case
!!      (Seity)         February 2012 add possibility to run with reversed 
!!                            vertical levels
!! ---------------------------------------------------------------------
!
!*       0. DECLARATIONS
!
USE MODD_PARAMETERS, ONLY: JPVEXT_TURB
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
!
IMPLICIT NONE
!
!
!*       0.1 declarations of arguments
!
TYPE(DIMPHYEX_t),          INTENT(IN)   :: D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),    INTENT(IN)  :: PVARM       ! variable at t-1  
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),    INTENT(IN)  :: PA          ! upper diag. elements
REAL,                      INTENT(IN)  :: PTSTEP      ! Double time step
REAL,                      INTENT(IN)  :: PEXPL,PIMPL ! weights of the temporal scheme
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),    INTENT(IN)  :: PRHODJ      ! (dry rho)*J
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),    INTENT(IN)  :: PSOURCE     ! source term of PVAR    
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),    INTENT(OUT) :: PVARP       ! variable at t+1        
!
!*       0.2 declarations of local variables
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT)  :: ZY ,ZGAM 
                                         ! RHS of the equation, 3D work array
REAL, DIMENSION(D%NIT,D%NJT)         :: ZBET
                                         ! 2D work array
INTEGER                              :: JI,JJ,JK            ! loop counter
INTEGER                              :: IKB,IKE,IIB,IIE,IJB,IJE       ! inner vertical limits
INTEGER                              :: IKT           ! array size in k direction
INTEGER                              :: IKTB,IKTE     ! start, end of k loops in physical domain 

!
! ---------------------------------------------------------------------------
!                                              
!*      1.  COMPUTE THE RIGHT HAND SIDE
!           ---------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('TRIDIAG',0,ZHOOK_HANDLE)
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
DO JJ=IJB,IJE 
  DO JI=IIB,IIE 
    ZY(JI,JJ,IKB) = PVARM(JI,JJ,IKB)  + PTSTEP*PSOURCE(JI,JJ,IKB) -   &
    PEXPL / PRHODJ(JI,JJ,IKB) * PA(JI,JJ,IKB+D%NKL) * &
    (PVARM(JI,JJ,IKB+D%NKL) - PVARM(JI,JJ,IKB))
  ENDDO
ENDDO
!
DO JK=IKTB+1,IKTE-1
  DO JJ=IJB,IJE 
    DO JI=IIB,IIE 
      ZY(JI,JJ,JK)= PVARM(JI,JJ,JK)  + PTSTEP*PSOURCE(JI,JJ,JK) -               &
      PEXPL / PRHODJ(JI,JJ,JK) *                                           &
      ( PVARM(JI,JJ,JK-D%NKL)*PA(JI,JJ,JK)                &
      -PVARM(JI,JJ,JK)*(PA(JI,JJ,JK)+PA(JI,JJ,JK+D%NKL))   &
      +PVARM(JI,JJ,JK+D%NKL)*PA(JI,JJ,JK+D%NKL)              &
      ) 
    ENDDO
  ENDDO
END DO
! 
DO JJ=IJB,IJE 
  DO JI=IIB,IIE 
    ZY(JI,JJ,IKE)= PVARM(JI,JJ,IKE) + PTSTEP*PSOURCE(JI,JJ,IKE) +               &
    PEXPL / PRHODJ(JI,JJ,IKE) * PA(JI,JJ,IKE) * (PVARM(JI,JJ,IKE)-PVARM(JI,JJ,IKE-D%NKL))
  ENDDO
ENDDO
!
!
!*       2.  INVERSION OF THE TRIDIAGONAL SYSTEM
!            -----------------------------------
!
IF ( PIMPL > 1.E-10 ) THEN
!
  !
  !  going up
  !
  DO JJ=IJB,IJE 
    DO JI=IIB,IIE 
      ZBET(JI,JJ) = 1. - PIMPL * PA(JI,JJ,IKB+D%NKL) / PRHODJ(JI,JJ,IKB)  ! bet = b(ikb)
      PVARP(JI,JJ,IKB) = ZY(JI,JJ,IKB) / ZBET(JI,JJ)
    ENDDO
  ENDDO               
  !
  DO JK = IKB+D%NKL,IKE-D%NKL,D%NKL
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE 
        ZGAM(JI,JJ,JK) = PIMPL * PA(JI,JJ,JK) / PRHODJ(JI,JJ,JK-D%NKL) / ZBET(JI,JJ)  
                                                    ! gam(k) = c(k-1) / bet
        ZBET(JI,JJ)    = 1. - PIMPL * (  PA(JI,JJ,JK) * (1. + ZGAM(JI,JJ,JK))  &
        + PA(JI,JJ,JK+D%NKL)                      &
        ) / PRHODJ(JI,JJ,JK)  
                                                    ! bet = b(k) - a(k)* gam(k)  
        PVARP(JI,JJ,JK)= ( ZY(JI,JJ,JK) - PIMPL * PA(JI,JJ,JK) / PRHODJ(JI,JJ,JK) &
        * PVARP(JI,JJ,JK-D%NKL)                                 &
        ) / ZBET(JI,JJ)
                                        ! res(k) = (y(k) -a(k)*res(k-1))/ bet 
      ENDDO
    ENDDO
  END DO
  DO JJ=IJB,IJE 
    DO JI=IIB,IIE 
  ! special treatment for the last level
      ZGAM(JI,JJ,IKE) = PIMPL * PA(JI,JJ,IKE) / PRHODJ(JI,JJ,IKE-D%NKL) / ZBET(JI,JJ) 
                                                    ! gam(k) = c(k-1) / bet
      ZBET(JI,JJ)    = 1. - PIMPL * (  PA(JI,JJ,IKE) * (1. + ZGAM(JI,JJ,IKE))  &
      ) / PRHODJ(JI,JJ,IKE)  
                                                    ! bet = b(k) - a(k)* gam(k)  
      PVARP(JI,JJ,IKE)= ( ZY(JI,JJ,IKE) - PIMPL * PA(JI,JJ,IKE) / PRHODJ(JI,JJ,IKE) &
      * PVARP(JI,JJ,IKE-D%NKL)                      &
      ) / ZBET(JI,JJ)
                                       ! res(k) = (y(k) -a(k)*res(k-1))/ bet 
  !
  !  going down
  !
    ENDDO
  ENDDO
  DO JK = IKE-D%NKL,IKB,-1*D%NKL
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE 
        PVARP(JI,JJ,JK) = PVARP(JI,JJ,JK) - ZGAM(JI,JJ,JK+D%NKL) * PVARP(JI,JJ,JK+D%NKL) 
      ENDDO
    ENDDO
  END DO
!
ELSE
! 
  DO JK=IKTB,IKTE
    DO JJ=IJB,IJE 
      DO JI=IIB,IIE 
        PVARP(JI,JJ,JK) = ZY(JI,JJ,JK)
      ENDDO
    ENDDO
  END DO
!
END IF 
!
!
!*       3.  FILL THE UPPER AND LOWER EXTERNAL VALUES
!            ----------------------------------------
!
DO JJ=IJB,IJE 
  DO JI=IIB,IIE 
    PVARP(JI,JJ,D%NKA)=PVARP(JI,JJ,IKB)
    PVARP(JI,JJ,D%NKU)=PVARP(JI,JJ,IKE)
  ENDDO
ENDDO
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('TRIDIAG',1,ZHOOK_HANDLE)
END SUBROUTINE TRIDIAG
END MODULE MODE_TRIDIAG
