

    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    ParseCSDataBlockSUBG.f90
    !> \brief   Parse the data block section corresponding to a SUBG phase of a ChemSage data-file.
    !> \author  M.H.A. Piro
    !> \date    Mar. 4, 2018
    !> \sa      ParseCSDataFile.f90
    !> \sa      ParseCSDataBlock.f90
    !> \sa      ParseCSDataBlockGibbs.f90
    !> \todo    There are a number of lines in SUBG phases that I do not yet understand.
    !!           I've asked some experts and they don't know either, which tells me that
    !!           they're not important. Once I
    !!           gain more experience with these models, this will likely become more clear.
    !
    !
    ! DISCLAIMER
    ! ==========
    !
    ! All of the programming herein is original unless otherwise specified and is completely
    ! independent of ChemApp and related products, including Solgas, Solgasmix, Fact, FactSage
    ! and ChemSage.
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer      Description of change
    !   ----            ----------      ---------------------
    !   03/04/2018      M.H.A. Piro     Original code
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to parse the "data block" section of a ChemSage data-file
    !! containing a "SUBG" phase, which represents the modified quasichemical model. This phase differs
    !! from many other types of thermodynamic models in that it attempts to capture Short Range Order (SRO)
    !! in liquid or solid solutions. This is achieved by focusing on pairs of species, rather than the species
    !! themselves. For more information, see the following paper:
    !!
    !! A.D. Pelton, S.A. Degterov, G. Eriksson, C. Roberlin, Y. Dessureault, "The Modified Quasichemical
    !! Model I -- Binary Solutions", Metallurgical and Materials Transactions B, 31B (2000) 651-659.
    !!
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! INFO                      A scalar integer that indicates a successful exit or identifies an error.
    ! nSpeciesCS                Number of species in the system (combined solution species and pure
    !                            separate phases).
    ! nGibbsEqSpecies           Number of Gibbs energy equations for a particular species.
    ! iSpeciesAtomsCS           Integer matrix representing the number of atoms of a particular
    !                            elements in a species.
    ! iParticlesPerMoleCS       An integer vector containing the number of particles per mole of the
    !                            constituent species formula mass.  The default value is 1.
    ! cSolnPhaseNameCS          The name of a solution phase.
    ! cSolnPhaseTypeCS          The type of a solution phase.
    ! cSolnPhaseTypeSupport     A character array representing solution phase types that are supported.
    ! iRegularParamCS           An integer matrix representing the parameter index for the first dimension
    !                            and the mixing terms on the second dimension.  For the second dimension, the
    !                            first coefficient indicates whether the parameter is a binary or ternary term (n),
    !                            the next n coefficients correspond to the constituent indices, and the last
    !                            coefficient corresponds to the exponent.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine ParseCSDataBlockSUBM( i )

    USE ModuleParseCS

    implicit none

    integer                     :: i, j, k, l, n, x, y, p, a, b, nPairs, nCSCS, nTotalConst
    integer                     :: xa, ya, nA2X2, iax, iay, ibx, iby, ia2x2, ia2y2, ib2x2, ib2y2
    integer                     :: iaaxy, ibbxy, iabxx, iabyy
    integer,     dimension(10)  :: iTempVec
    real(8)                     :: qa, qb, qx, qy, za, zb, zx, zy, dF, dCoax, dCobx, dCoay, dCoby
    real(8),     dimension(20)  :: dTempVec
    character(8),dimension(20)  :: cTempVec
    logical, dimension(:), allocatable :: lPairSet

    real(8), dimension(nSpeciesCS,nElementsCS) :: dStoichSpeciesOld

    nCSCS = nCountSublatticeCS

    ! Initialize variables:
    dTempVec = 0D0
    iTempVec = 0

    ! This line contains N integers (where N is the number of sublattices)
    ! where each integer represents the number of constituents on the respective
    ! sublattice. There are always two sublattices for SUBG phases.
    read (1,*,IOSTAT = INFO) nSublatticeElementsCS(nCSCS,1:2)
    nConstituentSublatticeCS(nCSCS,1:2) = nSublatticeElementsCS(nCSCS,1:2)
    nSublatticePhaseCS(nCSCS) = 2
    nTotalConst = nConstituentSublatticeCS(nCSCS,1)+nConstituentSublatticeCS(nCSCS,2)
    allocate(dStoichConstituentCS(nTotalConst,nElementsCS))
    dStoichConstituentCS = 0D0

    nPairs = nSublatticeElementsCS(nCSCS,1) * nSublatticeElementsCS(nCSCS,2)

    ! Read in names of constituents on first sublattice:
    read (1,*,IOSTAT = INFO) cConstituentNameSUBCS(nCSCS,1,1:nSublatticeElementsCS(nCSCS,1))

    ! Read in names of constituents on second sublattice: (ignore for now): 143 1
    read (1,*,IOSTAT = INFO) cConstituentNameSUBCS(nCSCS,2,1:nSublatticeElementsCS(nCSCS,2))

    ! Read in the charge of each constituent on the first sublattice.
    read (1,*,IOSTAT = INFO) dSublatticeChargeCS(nCSCS,1,1:nSublatticeElementsCS(nCSCS,1))

    ! Chemical groups on sublattice 1:
    read (1,*,IOSTAT = INFO) iChemicalGroupCS(nCSCS,1,1:nSublatticeElementsCS(nCSCS,1))

    ! Read in the charge of each constituent on the second sublattice.
    read (1,*,IOSTAT = INFO) dSublatticeChargeCS(nCSCS,2,1:nSublatticeElementsCS(nCSCS,2))

    ! Chemical groups on sublattice 2:
    read (1,*,IOSTAT = INFO) iChemicalGroupCS(nCSCS,2,1:nSublatticeElementsCS(nCSCS,2))

    ! This entry appears to represent the IDs matching constituents on the first sublattice to species:
    read (1,*,IOSTAT = INFO) iConstituentSublatticeCS(nCSCS, 1, 1:nPairs)

    ! This entry appears to represent the IDs matching constituents on the second sublattice to species:
    read (1,*,IOSTAT = INFO) iConstituentSublatticeCS(nCSCS, 2, 1:nPairs)

    ! Set up default pair IDs and coordination numbers


    ! Loop through excess mixing parameters:
    j = 0
    LOOP_ExcessMixingSUBG: do
        j = j + 1
        ! Read in number of constituents involved in parameter:
        read (1,*,IOSTAT = INFO) iRegularParamCS(nParamCS+1,1)

        ! The end of the parameter listing is marked by "0"
        ! or a negative number indicating the number of extra parameter lines.
        ! These lines indicate interpolation schemes, but I don't understand
        ! what these add, given that we can already generate interpolation
        ! schemes based on the chemical groups.
        if (iRegularParamCS(nParamCS+1,1) <= 0) then
            do k = 1, -iRegularParamCS(nParamCS+1,1)
                read (1,*,IOSTAT = INFO) cTempVec(1:10)
            end do
            exit LOOP_ExcessMixingSUBG
        end if

        ! Check if the parameter is binary or ternary:
        if ((iRegularParamCS(nParamCS+1,1) == 3) .OR. (iRegularParamCS(nParamCS+1,1) == 4)) then

            ! Count the number of parameters:
            nParamCS = nParamCS + 1

            ! Mixing terms:
            read (1,*,IOSTAT = INFO) iRegularParamCS(nParamCS,1:6), dRegularParamCS(nParamCS,1:6)

        else
            !! This parameter is not recognized; record an error.
            INFO = 10000 + 1000*j + i
            return
        end if

    end do LOOP_ExcessMixingSUBG

    ! Report an error if necessary:
    if (INFO /= 0) INFO = 1600 + i

    !deallocate(lPairSet,dStoichConstituentCS)

    return

end subroutine ParseCSDataBlockSUBM
