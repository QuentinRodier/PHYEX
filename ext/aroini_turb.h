INTERFACE
SUBROUTINE AROINI_TURB(PLINI,OHARATU,OSTATNW,OSUBG_COND)
USE PARKIND1  ,ONLY : JPRB
REAL(KIND=JPRB),INTENT(IN)::PLINI
LOGICAL, INTENT(IN) ::OHARATU
LOGICAL, INTENT(IN) ::OSTATNW
LOGICAL, INTENT(IN) :: OSUBG_COND
END SUBROUTINE AROINI_TURB
END INTERFACE