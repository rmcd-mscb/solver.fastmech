    module RivStagrMod_bmi
    USE fastmech


    IMPLICIT none
    type (fastmech_model) :: model

    CONTAINS

    SUBROUTINE STAGRBMI(STR_IN )
    IMPLICIT NONE
    INCLUDE "iriclib_f.h"
    CHARACTER(LEN=*), INTENT(IN) :: STR_IN
    
    call initialize_from_file(model, str_in)
    
    
    END SUBROUTINE STAGRBMI

    End module RivStagrMod_bmi