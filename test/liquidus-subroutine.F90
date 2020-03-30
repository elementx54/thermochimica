program liquidusSubroutine

    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModuleParseCS

    implicit none

    character(25) :: cLiqName
    real(8)       :: dLiqTol, dDriveScale, dLiquidus

    ! Specify units:
    cInputUnitTemperature = 'K'
    cInputUnitPressure    = 'atm'
    cInputUnitMass        = 'moles'
    cThermoFileName       = DATA_DIRECTORY // 'Kaye_NobleMetals.dat'

    ! Specify values:
    dTemperature          = 2000D0
    dPressure             = 1.0D0
    dElementMass          = 0D0
    dElementMass(46)      = 0.9D0                              ! Pd
    dElementMass(44)      = 0.1D0                              ! Ru
    dElementMass(43)      = 0.0D0                              ! Tc
    dElementMass(42)      = 1.9D0                              ! Mo

    ! Specify output and debug modes:
    iPrintResultsMode     = 2
    lDebugMode            = .FALSE.
    cLiqName              = 'LiqN'
    dLiqTol               = 1D-3
    dDriveScale           = 1D2

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    if (INFOThermo == 0)        call Thermochimica

    if (iPrintResultsMode > 0)  call PrintResults

    if (INFOThermo == 0)        call FindLiquidus(cLiqName, dLiquidus, dLiqTol, dDriveScale)

    if (iPrintResultsMode > 0)  call PrintResults

    print *, 'Liquidus Temperature is: ', dLiquidus

end program liquidusSubroutine
