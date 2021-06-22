

    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    ParseCSDataBlockSUBM.f90
    !> \brief   Parse the data block section corresponding to a SUBM phase of a ChemSage data-file.
    !> \author  M.H.A. Piro
    !> \date    Mar. 4, 2018
    !> \sa      ParseCSDataFile.f90
    !> \sa      ParseCSDataBlock.f90
    !> \sa      ParseCSDataBlockGibbs.f90
    !> \todo    There are a number of lines in SUBM phases that I do not yet understand.
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
    !! containing a "SUBM" phase, which represents the modified quasichemical model. This phase differs
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
    ! sublattice. There are always two sublattices for SUBM phases.
    read (1,*,IOSTAT = INFO) nSublatticeElementsCS(nCSCS,1:2)
    nConstituentSublatticeCS(nCSCS,1:2) = nSublatticeElementsCS(nCSCS,1:2)
    nSublatticePhaseCS(nCSCS) = 2
    nTotalConst = nConstituentSublatticeCS(nCSCS,1)+nConstituentSublatticeCS(nCSCS,2)
    allocate(dStoichConstituentCS(nTotalConst,nElementsCS))
    dStoichConstituentCS = 0D0

    nPairs = nSublatticeElementsCS(nCSCS,1) * nSublatticeElementsCS(nCSCS,2)

    ! Read in names of constituents on first sublattice:
    read (1,*,IOSTAT = INFO) cConstituentNameSUBCS(nCSCS,1,1:nSublatticeElementsCS(nCSCS,1))

    ! Read in names of constituents on second sublattice: (ignore for now):
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
    dCoordinationNumberCS(nCSCS,1:nMaxSpeciesPhaseCS,1:4) = 0D0
    do y = 1, nSublatticeElementsCS(nCSCS,2)
        LOOP_sroPairsOuter: do x = 1, nSublatticeElementsCS(nCSCS,2)
            if (x == y) then
                p = (x - 1) * (nSublatticeElementsCS(nCSCS,1) * (nSublatticeElementsCS(nCSCS,1) + 1) / 2)
            else if (x > y) then
                cycle LOOP_sroPairsOuter
            else
                p = (nSublatticeElementsCS(nCSCS,2) + (x - 1) + ((y-2)*(y-1)/2)) &
                  * (nSublatticeElementsCS(nCSCS,1) * (nSublatticeElementsCS(nCSCS,1) + 1) / 2)
            end if
            do k = 1, nSublatticeElementsCS(nCSCS,1)
                LOOP_sroPairsInner: do j = 1, nSublatticeElementsCS(nCSCS,1)
                    if (j == k) then
                        l = j
                    else if (j > k) then
                        cycle LOOP_sroPairsInner
                    else
                        l = nSublatticeElementsCS(nCSCS,1) + j + ((k-2)*(k-1)/2)
                    end if
                    iPairIDCS(nCSCS, l + p, 1) = j
                    iPairIDCS(nCSCS, l + p, 2) = k
                    iPairIDCS(nCSCS, l + p, 3) = x + nSublatticeElementsCS(nCSCS,1)
                    iPairIDCS(nCSCS, l + p, 4) = y + nSublatticeElementsCS(nCSCS,1)
                    end do LOOP_sroPairsInner
            end do
        end do LOOP_sroPairsOuter
    end do



    ! Copy previously-read end member info into appropriate variables before it gets overwritten by
    ! quadruplet data calculated below.
    cPairNameCS(nCSCS,1:nPairsSROCS(nCSCS,1)) = &
                cSpeciesNameCS((nSpeciesPhaseCS(i-1)+1):(nSpeciesPhaseCS(i-1)+nPairsSROCS(nCSCS,1)))
    dStoichSpeciesOld = dStoichSpeciesCS(1:nSpeciesCS,1:nElementsCS)
    dStoichPairsCS(nCSCS,1:nPairsSROCS(nCSCS,2),1:nElementsCS) &
                  = dStoichSpeciesCS((nSpeciesPhaseCS(i-1) + 1):nSpeciesPhaseCS(i),1:nElementsCS)
    dStoichSpeciesCS((nSpeciesPhaseCS(i-1) + 1):nSpeciesPhaseCS(i),1:nElementsCS) = 0D0

    ! Loop through all pairs to calculate stoichiometry entries for quadruplets:
    do j = 1, nPairsSROCS(nCSCS,2)
        a = iPairIDCS(nCSCS, j, 1)
        b = iPairIDCS(nCSCS, j, 2)
        x = iPairIDCS(nCSCS, j, 3)
        y = iPairIDCS(nCSCS, j, 4)

        xa = x - nSublatticeElementsCS(nCSCS,1)
        ya = y - nSublatticeElementsCS(nCSCS,1)

        nA2X2 = nSublatticeElementsCS(nCSCS,1) * nSublatticeElementsCS(nCSCS,2)
        do k = 1, nA2X2
            if   ((iConstituentSublatticeCS(nCSCS,1,k) == a) &
            .AND. (iConstituentSublatticeCS(nCSCS,2,k) == xa)) then
                iax = k
            end if
            if   ((iConstituentSublatticeCS(nCSCS,1,k) == b) &
            .AND. (iConstituentSublatticeCS(nCSCS,2,k) == xa)) then
                ibx = k
            end if
            if   ((iConstituentSublatticeCS(nCSCS,1,k) == a) &
            .AND. (iConstituentSublatticeCS(nCSCS,2,k) == ya)) then
                iay = k
            end if
            if   ((iConstituentSublatticeCS(nCSCS,1,k) == b) &
            .AND. (iConstituentSublatticeCS(nCSCS,2,k) == ya)) then
                iby = k
            end if
        end do

        ia2x2 = a + ((xa - 1) * (nSublatticeElementsCS(nCSCS,1) &
                                * (nSublatticeElementsCS(nCSCS,1) + 1) / 2))
        ib2x2 = b + ((xa - 1) * (nSublatticeElementsCS(nCSCS,1) &
                                * (nSublatticeElementsCS(nCSCS,1) + 1) / 2))
        ia2y2 = a + ((ya - 1) * (nSublatticeElementsCS(nCSCS,1) &
                                * (nSublatticeElementsCS(nCSCS,1) + 1) / 2))
        ib2y2 = b + ((ya - 1) * (nSublatticeElementsCS(nCSCS,1) &
                                * (nSublatticeElementsCS(nCSCS,1) + 1) / 2))

        l = j + nSpeciesPhaseCS(i-1)

        ! Create quadruplet names
        cSpeciesNameCS(l) = TRIM(cConstituentNameSUBCS(nCSCS,1,a)) // '-' &
                         // TRIM(cConstituentNameSUBCS(nCSCS,1,b)) // '-' &
                         // TRIM(cConstituentNameSUBCS(nCSCS,2,x - nSublatticeElementsCS(nCSCS,1))) // '-' &
                         // TRIM(cConstituentNameSUBCS(nCSCS,2,y - nSublatticeElementsCS(nCSCS,1)))

        dCoax = dConstituentCoefficientsCS(nCSCS,iax,1)
        dCobx = dConstituentCoefficientsCS(nCSCS,ibx,1)
        dCoay = dConstituentCoefficientsCS(nCSCS,iay,1)
        dCoby = dConstituentCoefficientsCS(nCSCS,iby,1)
        do k = 1, nElementsCS
            dStoichSpeciesCS(l,k) = dStoichSpeciesCS(l,k) + &
                                   (dStoichPairsCS(nCSCS,iax,k) / dCoordinationNumberCS(nCSCS, j, 1)) / (2D0 * dCoax)
            dStoichSpeciesCS(l,k) = dStoichSpeciesCS(l,k) + &
                                   (dStoichPairsCS(nCSCS,ibx,k) / dCoordinationNumberCS(nCSCS, j, 2)) / (2D0 * dCobx)
            dStoichSpeciesCS(l,k) = dStoichSpeciesCS(l,k) + &
                                   (dStoichPairsCS(nCSCS,iay,k) / dCoordinationNumberCS(nCSCS, j, 1)) / (2D0 * dCoay)
            dStoichSpeciesCS(l,k) = dStoichSpeciesCS(l,k) + &
                                   (dStoichPairsCS(nCSCS,iby,k) / dCoordinationNumberCS(nCSCS, j, 2)) / (2D0 * dCoby)
        end do
    end do

    ! Loop through excess mixing parameters:
    j = 0
    LOOP_ExcessMixingSUBM: do
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
            exit LOOP_ExcessMixingSUBM
        end if

        ! Check if the parameter is binary or ternary:
        if ((iRegularParamCS(nParamCS+1,1) == 3) .OR. (iRegularParamCS(nParamCS+1,1) == 4)) then

            ! Count the number of parameters:
            nParamCS = nParamCS + 1

            ! Mixing terms:
            read (1,*,IOSTAT = INFO) cRegularParamCS(nParamCS), iRegularParamCS(nParamCS,2:9)
            if (.NOT.((cRegularParamCS(nParamCS) == 'G') &
                 .OR. (cRegularParamCS(nParamCS) == 'Q') .OR. (cRegularParamCS(nParamCS) == 'R') &
                 .OR. (cRegularParamCS(nParamCS) == 'B'))) then
                INFO = 10000 + 1000*j + i
                return
            end if

            ! According to Patrice Chartrand, he has no idea what these two lines mean. Ignore.
            read (1,*,IOSTAT = INFO) dTempVec(1:6)
            read (1,*,IOSTAT = INFO) dTempVec(1:6)

            ! Read in the excess gibbs energy of mixing terms.
            read (1,*,IOSTAT = INFO) iRegularParamCS(nParamCS,10:11), dRegularParamCS(nParamCS,1:6)

        else
            !! This parameter is not recognized; record an error.
            INFO = 10000 + 1000*j + i
            return
        end if

    end do LOOP_ExcessMixingSUBM

    ! Report an error if necessary:
    if (INFO /= 0) INFO = 1600 + i

    deallocate(lPairSet,dStoichConstituentCS)

    return

end subroutine ParseCSDataBlockSUBM
