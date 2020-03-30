subroutine FindLiquidus(cLiqName, dLiquidus, dLiqTol, dDriveScale)

    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModuleParseCS

    implicit none

    character(25), intent(in)  :: cLiqName
    logical       :: liqFound, tolReached
    integer       :: i, j, iLiq
    real(8)       :: dLiqGibbs, dEquilGibbs, dLiqDrive, dStep
    real(8),       intent(in)  :: dLiqTol, dDriveScale
    real(8),       intent(out) :: dLiquidus

    liqFound = .FALSE.
    lReinitRequested = .TRUE.

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

    if (DABS(dLiqDrive) < dLiqTol) tolReached = .TRUE.

    do while ((.NOT. tolReached) .AND. (INFOThermo == 0))
      call SaveReinitData
      call ResetThermo
      dStep = dDriveScale * dLiqDrive
      dTemperature = dTemperature + dStep
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

      if (DABS(dLiqDrive * dDriveScale) < dLiqTol) tolReached = .TRUE.

    end do

    print *, '-------------------------------------------------------------'

    do while ((.NOT. liqFound) .AND. (INFOThermo == 0))
      call SaveReinitData
      call ResetThermo
      dStep = dLiqTol
      dTemperature = dTemperature + dStep
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

      do i = 1, nSolnPhases
        if (cSolnPhaseName(-iAssemblage(nElements - i + 1)) == cLiqName) liqFound = .TRUE.
      end do

    end do

    dLiquidus = dTemperature

end subroutine FindLiquidus
