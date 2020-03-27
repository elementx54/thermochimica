program liquidus

    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModuleParseCS

    implicit none

    character(25) :: cLiqName
    logical       :: liqFound
    integer       :: i, j, iLiq
    real(8)       :: dLiqGibbs, dEquilGibbs, dLiqDrive, dLiqTol

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
    cLiqName = 'LiqN'
    liqFound = .FALSE.
    dLiqTol  = 1D-8
    lReinitRequested = .TRUE.

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    if (INFOThermo == 0)        call Thermochimica

    do i = 1, nSolnPhasesSys
      if (cSolnPhaseName(i)== cLiqName) iLiq = i
    end do

    dLiqGibbs = 0
    dEquilGibbs = 0
    do i = nSpeciesPhase(iLiq - 1) + 1, nSpeciesPhase(iLiq)
      dLiqGibbs = dLiqGibbs + (dMolFraction(i) * dChemicalPotential(i))

      do j = 1, nElements
        dEquilGibbs = dEquilGibbs + (dMolFraction(i) * dStoichSpecies(i,j) * dElementPotential(j))
      end do

    end do

    dLiqDrive = dLiqGibbs - dEquilGibbs

    print *, dTemperature, dLiqGibbs, dEquilGibbs, dLiqDrive

    if (DABS(dLiqDrive) < dLiqTol) liqFound = .TRUE.

    do while ((.NOT. liqFound) .AND. (INFOThermo == 0))
      call SaveReinitData
      call ResetThermo
      dTemperature = dTemperature + 1D2 * dLiqDrive
      call Thermochimica

      dLiqGibbs = 0
      dEquilGibbs = 0
      do i = nSpeciesPhase(iLiq - 1) + 1, nSpeciesPhase(iLiq)
        dLiqGibbs = dLiqGibbs + (dMolFraction(i) * dChemicalPotential(i))

        do j = 1, nElements
          dEquilGibbs = dEquilGibbs + (dMolFraction(i) * dStoichSpecies(i,j) * dElementPotential(j))
        end do

      end do

      dLiqDrive = dLiqGibbs - dEquilGibbs

      print *, dTemperature, dLiqGibbs, dEquilGibbs, dLiqDrive

      if (DABS(dLiqDrive) < dLiqTol) liqFound = .TRUE.

    end do

    if (iPrintResultsMode > 0)  call PrintResults


end program liquidus
