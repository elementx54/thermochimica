program dmu

    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModuleParseCS

    implicit none

    integer :: iElem, iDerMode
    real(8) :: dDelta, dDChemPot


    ! Specify units:
    cInputUnitTemperature = 'K'
    cInputUnitPressure    = 'atm'
    cInputUnitMass        = 'moles'
    cThermoFileName       = DATA_DIRECTORY // 'Ted_Data_files/UO2-FP_no_Gd_hcp_A3.dat'

    ! Specify values:
    dTemperature          = 1500D0
    dPressure             = 1.0D0
    dElementMass          = 0D0

    dElementMass(8)        =      8.4030E+03
    dElementMass(92)       =      3.6820E+03
    dElementMass(54)       =      1.3150E+02
    dElementMass(42)       =      1.0080E+02
    dElementMass(40)       =      9.8930E+01
    dElementMass(44)       =      8.9420E+01
    dElementMass(60)       =      8.0156E+01
    dElementMass(46)       =      6.7630E+01
    dElementMass(94)       =      6.5847E+01
    dElementMass(56)       =      5.7420E+01
    dElementMass(55)       =      5.7380E+01
    dElementMass(58)       =      5.4720E+01
    dElementMass(64)       =      3.0632E+01
    dElementMass(57)       =      2.5200E+01
    dElementMass(59)       =      2.1986E+01
    dElementMass(43)       =      1.8530E+01
    dElementMass(52)       =      1.2290E+01
    dElementMass(39)       =      1.1550E+01
    dElementMass(36)       =      1.0110E+01
    dElementMass(37)       =      9.3000E+00
    dElementMass(45)       =      7.5310E+00
    dElementMass(48)       =      7.3460E+00
    dElementMass(53)       =      5.3120E+00
    dElementMass(2)        =      3.4670E+00

    ! Specify output and debug modes:
    iPrintResultsMode     = 2
    lDebugMode            = .FALSE.

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    if (INFOThermo == 0)        call Thermochimica

    iElem = 8
    dDelta = 1D-5
    iDerMode = 4
    INFO = 0

    call CompDChemicalPotentialDC(iElem, dDelta, iDerMode, dDChemPot, INFO)
    print *, dDChemPot

    ! Destruct everything:
    call ResetThermoAll

    ! Call the debugger:
    call ThermoDebug

end program dmu
