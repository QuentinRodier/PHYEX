INTERFACE
SUBROUTINE AROINI_NSV(KSV,KSV_CHEMBEG, KSV_CHEMEND, KSV_AERBEG, KSV_AEREND, &
                      KSV_DSTBEG, KSV_DSTEND, KSV_DSTDEPBEG, KSV_DSTDEPEND,&
                      KSV_CO2)
USE PARKIND1  ,ONLY : JPIM

INTEGER(KIND=JPIM), INTENT(IN)::KSV,KSV_CHEMBEG, KSV_CHEMEND, KSV_AERBEG,&
                              & KSV_AEREND , KSV_DSTBEG , KSV_DSTEND,&
                              & KSV_DSTDEPBEG, KSV_DSTDEPEND, KSV_CO2

END SUBROUTINE AROINI_NSV
END INTERFACE