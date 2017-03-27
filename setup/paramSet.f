MODULE paramSet
    integer, save :: &
    TIPCALC = 0, &
    SETFROZEN = 1, &
    SETCORECA = 1, &
    IRIG = 0
    real, save :: &
    DINT0 = 8., &
    DINT1 = 12., &
    BONDCUTOFF = 8., &
    OFFSETX = 0., &
    OFFSETY = 0., &
    OFFSETZ = 0., &
    FSOLV = 0.5, &
    FVDW = 0.6, &
    EPS = 1., &
    DPS = 7.0, &
    DHF = 5.0

CONTAINS
!===============================================================================
    subroutine readInputParamSet(unit_i)
        integer, intent(IN) :: unit_i

        namelist /setup/tipcalc, dint0, dint1, setfrozen, setcoreca, &
        offsetx, offsety, offsetz, irig, bondcutoff, fsolv, fvdw, eps, dps, dhf
        read(unit_i, SETUP)

    end subroutine readInputParamSet
!===============================================================================
    subroutine printInputParamSet(unit_o)
        integer, intent(IN) :: unit_o
        !TODO pendent format, provisional namelist
        namelist /setup/tipcalc, dint0, dint1, setfrozen, setcoreca, &
        offsetx, offsety, offsetz, irig, bondcutoff, fsolv, fvdw, eps, dps, dhf
        write (unit_o,*) "INPUT PARAMETERS"
        write (unit_o, SETUP)
    end subroutine printInputParamSet
        
END MODULE paramSet

