subroutine get_ipair(ipair_ret, iatm, jatm, bonded_pairs, npair)
    implicit none

    integer, intent(in) :: iatm, jatm, npair
    integer, intent(in) :: bonded_pairs(npair, 2)
    integer, intent(out):: ipair_ret
    integer :: ipair, iatm_, jatm_

    ipair_ret = 0
    do ipair=1, npair
        iatm_ = bonded_pairs(ipair, 1)
        jatm_ = bonded_pairs(ipair, 2)

        if (iatm_ == iatm .and. jatm_ == jatm) then
            ipair_ret = ipair
            exit
        else if (iatm_ == jatm .and. jatm_ == iatm) then
            ipair_ret = - ipair
            exit
        end if

    end do

end subroutine

subroutine get_ibnd_to_ipair(ibnd_to_ipair, two_atoms, bonded_pairs, &
                            & nbond, npair)
    implicit none

    integer, intent(in) :: nbond, npair
    integer, intent(in) :: two_atoms(nbond, 2), bonded_pairs(npair, 2)
    integer, intent(out):: ibnd_to_ipair(nbond)

    integer :: ibnd, ipair, iatm, jatm

    do ibnd=1, nbond
        iatm = two_atoms(ibnd, 1)
        jatm = two_atoms(ibnd, 2)

        ! iatm, jatm
        call get_ipair(ipair, iatm, jatm, bonded_pairs, npair)
        ibnd_to_ipair(ibnd) = ipair

    end do

end subroutine

subroutine get_iang_to_ipair(iang_to_ipair, three_atoms, bonded_pairs, &
                            & nangle, npair)
    implicit none

    integer, intent(in) :: nangle, npair
    integer, intent(in) :: three_atoms(nangle, 3), bonded_pairs(npair, 2)
    integer, intent(out):: iang_to_ipair(nangle, 3)

    integer :: iang, ipair, iatm, jatm, katm

    do iang=1, nangle
        iatm = three_atoms(iang, 1)
        jatm = three_atoms(iang, 2)
        katm = three_atoms(iang, 3)

        ! iatm, jatm
        call get_ipair(ipair, iatm, jatm, bonded_pairs, npair)
        iang_to_ipair(iang, 1) = ipair
        ! iatm, katm
        call get_ipair(ipair, iatm, katm, bonded_pairs, npair)
        iang_to_ipair(iang, 2) = ipair
        ! jatm, katm
        call get_ipair(ipair, jatm, katm, bonded_pairs, npair)
        iang_to_ipair(iang, 3) = ipair

    end do

end subroutine

subroutine get_itor_to_ipair(itor_to_ipair, four_atoms, bonded_pairs, &
                            & ntorsion, npair)
    implicit none

    integer, intent(in) :: ntorsion, npair
    integer, intent(in) :: four_atoms(ntorsion, 4), bonded_pairs(npair, 2)
    integer, intent(out):: itor_to_ipair(ntorsion, 6)

    integer :: itor, ipair, iatm, jatm, katm, latm

    do itor=1, ntorsion
        iatm = four_atoms(itor, 1)
        jatm = four_atoms(itor, 2)
        katm = four_atoms(itor, 3)
        latm = four_atoms(itor, 4)

        ! iatm, jatm
        call get_ipair(ipair, iatm, jatm, bonded_pairs, npair)
        itor_to_ipair(itor, 1) = ipair
        ! iatm, katm
        call get_ipair(ipair, iatm, katm, bonded_pairs, npair)
        itor_to_ipair(itor, 2) = ipair
        ! iatm, latm
        call get_ipair(ipair, iatm, latm, bonded_pairs, npair)
        itor_to_ipair(itor, 3) = ipair

        ! jatm, katm
        call get_ipair(ipair, jatm, katm, bonded_pairs, npair)
        itor_to_ipair(itor, 4) = ipair
        ! jatm, latm
        call get_ipair(ipair, jatm, latm, bonded_pairs, npair)
        itor_to_ipair(itor, 5) = ipair
        ! katm, latm
        call get_ipair(ipair, katm, latm, bonded_pairs, npair)
        itor_to_ipair(itor, 6) = ipair

    end do

end subroutine

