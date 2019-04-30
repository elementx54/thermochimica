
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

  integer                              :: i, j, k, l, iFirst, iLast, nVar, nTotalPhases, nMaxSolutionPhases
  integer,     dimension(0:118)        :: iElementsUsed
  real(8)                              :: dChemPotSpecies, dMaxDrivingForce, dStoichSum, dDrivingForceCutoff
  real(8), dimension(:),   allocatable :: dChemicalPotentialStar, dDrivingForce, dDrivingForceSum
  real(8), dimension(:),   allocatable :: dMinDrivingForces
  logical      :: lAddPhase

  ! Initialize storage variables if not allocated already
  if (.NOT. lRestartAvailable) then
    print *, 'Restart requested but data not available'
    return
  endif

  ! Save old chemical potential data
  if (.NOT. allocated(dChemicalPotential)) allocate(dChemicalPotential(nSpecies))
  if (.NOT. allocated(dElementPotential)) allocate(dElementPotential(nElements))
  dElementPotential   = dElementPotential_Old

  allocate(dChemicalPotentialStar(nSpecies),dDrivingForce(nSpecies))
  dChemPotSpecies = 0D0
  do k = 1, nSpecies
      dChemPotSpecies = 0D0
      dStoichSum = 0
      do j = 1, nElements
          dChemPotSpecies = dChemPotSpecies + dElementPotential(j) * dStoichSpecies(k,j)
          dStoichSum = dStoichSum + dStoichSpecies(k,j)
      end do
      if (k < MAXVAL(nSpeciesPhase)) then
          dChemPotSpecies = dChemPotSpecies / iParticlesPerMole(j)
      else
          dChemPotSpecies = dChemPotSpecies / dStoichSum
      end if
      ! print *, dChemPotSpecies, cSpeciesName(k), iParticlesPerMole(j)
      dChemicalPotentialStar(k) = dChemPotSpecies
  end do

  l = MAX(1,nSolnPhasesSys)
  allocate(dMolesSpecies(nSpecies))
  allocate(dPartialExcessGibbs(nSpecies),dPartialExcessGibbsLast(nSpecies))
  allocate(dUpdateVar(nElements*2))
  allocate(iterHistory(nElements,iterGlobalMax))
  allocate(dSumMolFractionSoln(l))
  allocate(dDrivingForceSoln(l))
  allocate(dEffStoichSolnPhase(l,nElements))
  allocate(iAssemblage(nElements))
  iAssemblage = 0
  allocate(dMolesPhase(nElements))
  dMolesPhase = 0D0
  allocate(dMolFraction(nSpecies))

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
      call CompExcessGibbsEnergy(j)
  end do


  do j = 1, nSpecies
      ! Update the driving force:
      if (j < MAXVAL(nSpeciesPhase)) then
          dDrivingForce(j) = dMolFraction(j) * (dChemicalPotential(j) - dChemicalPotentialStar(j))
      else
          dDrivingForce(j) = (dChemicalPotential(j) - dChemicalPotentialStar(j))
      end if
      ! print *, dDrivingForce(j), dChemicalPotential(j), dChemicalPotentialStar(j), iParticlesPerMole(j), &
      !          cSpeciesName(j), dMolFraction(j)
  end do

  nTotalPhases = nSolnPhasesSys + (nSpecies - MAXVAL(nSpeciesPhase))
  allocate(dDrivingForceSum(nTotalPhases))

  nConPhases = 0
  do j = MAXVAL(nSpeciesPhase) + 1, nSpecies
      i = j + nSolnPhasesSys - MAXVAL(nSpeciesPhase)
      dDrivingForceSum(i) = dDrivingForce(j)
      ! print *, dDrivingForceSum(i), cSpeciesName(j)
      if (dDrivingForce(j) < 1D-12) then
          nConPhases = nConPhases + 1
          iAssemblage(nConPhases) = j
      end if
  end do

  nMaxSolutionPhases = nElements - nConPhases
  allocate(dMinDrivingForces(nMaxSolutionPhases))
  dDrivingForceCutoff = 1D0
  dMinDrivingForces = dDrivingForceCutoff
  nSolnPhases = 0
  do j = 1, nSolnPhasesSys
      iFirst = nSpeciesPhase(j-1) + 1
      iLast  = nSpeciesPhase(j)
      nVar   = iLast - iFirst + 1
      dDrivingForceSum(j) = 0D0
      do i = iFirst, iLast
          dDrivingForceSum(j) = dDrivingForceSum(j) + dDrivingForce(i)
      end do
      ! print *, dDrivingForceSum(j), cSolnPhaseName(j)
      if (dDrivingForceSum(j) < MAXVAL(dMinDrivingForces)) then
          iAssemblage(nElements - nSolnPhases) = -j
          nSolnPhases = nSolnPhases + 1
          dMinDrivingForces(MAXLOC(dMinDrivingForces)) = dDrivingForceSum(j)
      end if
  end do

  call CompMolAllSolnPhases

  lRestartLoaded = .TRUE.

  deallocate(dDrivingForceSum,dDrivingForce,dChemicalPotentialStar,dMinDrivingForces)

end subroutine LoadRestartData
