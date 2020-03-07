subroutine CompDChemicalPotentialDC(iElem, dDelta, iDerMode, dDChemPot, INFO)

    USE ModuleThermo
    USE ModuleThermoIO

    implicit none

    integer,       intent(in)    :: iElem, iDerMode
    integer,       intent(out)   :: INFO
    integer                      :: iElemLoc
    real(8),       intent(in)    :: dDelta
    real(8),       intent(out)   :: dDChemPot
    real(8)                      :: dOriginalElemMass, dElemChemPotPP, dElemChemPotP, dElemChemPotM, dElemChemPotMM

    ! Save reinit data
    call SaveReinitData

    LOOP_findElem: do iElemLoc = 1, nElements
      if (iElementSystem(iElemLoc) == iElem) EXIT LOOP_findElem
    end do LOOP_findElem
    dOriginalElemMass = dElementMass(iElem)

    call ResetThermo

    ! Load reinit data
    lReinitRequested = .TRUE.

    ! Call Thermochimica again with new concentration:
    dElementMass = dElementMass * dNormalizeInput
    dElementMass(iElem) = dOriginalElemMass * (1D0 + dDelta)
    if (INFOThermo == 0)        call Thermochimica
    ! if (iPrintResultsMode > 0)  call PrintResults
    dElemChemPotP = dElementPotential(iElemLoc)

    call ResetThermo

    ! Call Thermochimica again with new concentration:
    dElementMass = dElementMass * dNormalizeInput
    dElementMass(iElem) = dOriginalElemMass * (1D0 - dDelta)
    if (INFOThermo == 0)        call Thermochimica
    ! if (iPrintResultsMode > 0)  call PrintResults
    dElemChemPotM = dElementPotential(iElemLoc)

    if (iDerMode == 4) then
        ! Call Thermochimica again with new concentration:
        dElementMass = dElementMass * dNormalizeInput
        dElementMass(iElem) = dOriginalElemMass * (1D0 + 2D0*dDelta)
        if (INFOThermo == 0)        call Thermochimica
        ! if (iPrintResultsMode > 0)  call PrintResults
        dElemChemPotPP = dElementPotential(iElemLoc)

        ! Call Thermochimica again with new concentration:
        dElementMass = dElementMass * dNormalizeInput
        dElementMass(iElem) = dOriginalElemMass * (1D0 - 2D0*dDelta)
        if (INFOThermo == 0)        call Thermochimica
        ! if (iPrintResultsMode > 0)  call PrintResults
        dElemChemPotMM = dElementPotential(iElemLoc)

        ! 4-point first derivative
        dDChemPot = (-dElemChemPotPP + 8D0*dElemChemPotP - 8D0*dElemChemPotM + dElemChemPotMM) / (12D0 * dDelta * dOriginalElemMass)
    else
        ! 2-point first derivative
        dDChemPot = (dElemChemPotP - dElemChemPotM) / (2D0 * dDelta * dOriginalElemMass)
    end if
    dElementMass(iElem) = dOriginalElemMass

    dDChemPot = dDChemPot * dIdealConstant * dTemperature

    INFO = 0

end subroutine CompDChemicalPotentialDC
