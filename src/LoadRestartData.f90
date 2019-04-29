
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    LoadRestartData.f90
    !> \brief   Load data for restarting calculation from previous results.
    !> \author  M. Poschmann
    !> \sa      SaveRestartData.f90
    !
    !
    ! References:
    ! ===========
    !
    ! Revisions:
    ! ==========
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   30/11/2018      M. Poschmann         Create file.
    !
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this subroutine is to save all pertinent data such that a new call to Thermochimica
    !! may be restarted from that data.
    !
    ! Pertinent variables:
    ! ====================
    !> \param   lRestart        A logical indicating whether restart data is available.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine LoadRestartData

  USE ModuleThermo
  USE ModuleRestart
  USE ModuleThermoIO
  USE ModuleGEMSolver

  implicit none

  integer                              :: i, j, k
  integer,     dimension(0:118)        :: iElementsUsed
  real(8)                              :: dChemPotSpecies
  logical      :: lAddPhase

  ! Initialize storage variables if not allocated already
  if (.NOT. lRestartAvailable) then
    print *, 'Restart requested but data not available'
    return
  endif

  ! Check that the number of elements hasn't changed. If it has, don't load data.
  i = SIZE(iAssemblage_Old)
  if (i /= nElements) return
  ! Also check that the elements involved themselves are the same.
  iElementsUsed = min(ceiling(dElementMass),1)
  do j = 0, (nElementsPT)
    if (iElementsUsed(j) /= iElementsUsed_Old(j)) return
  enddo

  ! Save old chemical potential data
  if (.NOT. allocated(dChemicalPotential)) allocate(dChemicalPotential(nSpecies))
  if (.NOT. allocated(dElementPotential)) allocate(dElementPotential(nElements))
  ! dChemicalPotential  = dChemicalPotential_Old
  dElementPotential   = dElementPotential_Old

  dChemPotSpecies = 0D0
  do k = 1, nSpecies
      dChemPotSpecies = 0D0
      do j = 1, nElements
          dChemPotSpecies = dChemPotSpecies + dElementPotential(j) * dStoichSpecies(k,j)
      end do
      dChemPotSpecies = dChemPotSpecies / iParticlesPerMole(k)
      print *, dChemPotSpecies, dChemicalPotential_Old(k), iParticlesPerMole(k)
  end do

  do j = 1, nSolnPhasesSys
      lAddPhase         = .FALSE.
      ! Only compute the mole fractions of constituents belonging to unstable phases:
      if (lSolnPhases(j) .EQV. .FALSE.) then
          if (lMiscibility(j) .EQV. .TRUE.) then
              ! Check for a miscibility gap:
              call CheckMiscibilityGap(j,lAddPhase)
          else
              ! This is phase does not contain a miscibility gap.  Compute the mole fractions:
              call CompMolFraction(j)
          end if
      end if
  end do

  ! Save old phase data
  if (.NOT. allocated(iAssemblage)) allocate(iAssemblage(nElements))
  if (.NOT. allocated(dMolesPhase)) allocate(dMolesPhase(nElements))
  if (.NOT. allocated(dMolFraction)) allocate(dMolFraction(nSpecies))
  iAssemblage         = iAssemblage_Old
  dMolesPhase         = dMolesPhase_Old
  dMolFraction        = dMolFraction_Old

  lRestartLoaded = .TRUE.

end subroutine LoadRestartData
