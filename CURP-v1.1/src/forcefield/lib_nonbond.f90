
module without_bonded

    implicit none
    integer :: npair
    integer, allocatable :: bonded_pairs(:, :) ! (npair, 2)
    integer, allocatable :: begins(:) ! (natom)
    integer, allocatable :: ends(:)   ! (natom)

    integer :: iatm, jatm
    integer :: i, ipair

contains

    subroutine setup(bonded_pairs_in, natom, npair_in)
        implicit none

        integer, intent(in) :: npair_in, natom
        integer, intent(in) :: bonded_pairs_in(npair_in, 2)

        npair = npair_in

        allocate(bonded_pairs(npair_in, 2))
        bonded_pairs = bonded_pairs_in
        allocate(begins(natom))
        allocate(ends(natom))

        call get_begins_ends(begins, ends, bonded_pairs, npair, natom)

    end subroutine


    logical function is_bonded_pair(iatm_in, jatm_in)
        implicit none
        integer, intent(in) :: iatm_in, jatm_in
        integer :: iatm_, jatm_

        is_bonded_pair = .false.
        do ipair=begins(iatm_in), ends(iatm_in)
            iatm_ = bonded_pairs(ipair, 1)
            jatm_ = bonded_pairs(ipair, 2)

            if ((iatm_in==iatm_) .and. (jatm_in==jatm_)) then
                is_bonded_pair = .true.
                exit
            end if
        end do

    end function


    subroutine get_nonbonded_table(nonbonded_table, ntable_new, &
                            & base_table, ntable)

        implicit none

        integer, intent(in) :: ntable
        integer, intent(in) :: base_table(ntable, 3)

        integer, intent(out) :: ntable_new
        integer, intent(out) :: nonbonded_table(ntable*10, 3)

        integer :: itab, jatm_beg, jatm_end, jatm_beg_new, jatm_end_new
        integer :: itab_new

        itab_new = 0
        do itab=1, ntable
            iatm     = base_table(itab, 1)
            jatm_beg = base_table(itab, 2)
            jatm_end = base_table(itab, 3)

            jatm_beg_new = 0
            jatm_end_new = jatm_beg

            do jatm=jatm_beg, jatm_end

                if (is_bonded_pair(iatm, jatm)) then

                    if (jatm_beg_new > 0) then
                        itab_new = itab_new + 1
                        nonbonded_table(itab_new, 1) = iatm
                        nonbonded_table(itab_new, 2) = jatm_beg_new
                        nonbonded_table(itab_new, 3) = jatm_end_new

                        jatm_beg_new = 0

                    end if

                else
                    if (jatm_beg_new == 0) jatm_beg_new = jatm
                    jatm_end_new = jatm

                end if

            end do

            if (jatm_beg_new > 0) then
                itab_new = itab_new + 1
                nonbonded_table(itab_new, 1) = iatm
                nonbonded_table(itab_new, 2) = jatm_beg_new
                nonbonded_table(itab_new, 3) = jatm_end_new
            
            end if

        end do

        ntable_new = itab_new

    end subroutine

    subroutine get_begins_ends(begins, ends, bonded_pairs, npair, natom)

        implicit none
        integer, intent(in) :: npair, natom
        integer, intent(in) :: bonded_pairs(npair, 2)

        integer, intent(out):: begins(natom)
        integer, intent(out):: ends(natom)

        integer :: nums(natom)

        integer :: ipair, iatm, jatm, iatm_prev

        begins(:) = 0
        nums(:)   = 0

        iatm_prev = 0
        do ipair=1, npair
            iatm = bonded_pairs(ipair, 1)
            jatm = bonded_pairs(ipair, 2)
            nums(iatm) = nums(iatm) + 1

            if (iatm /= iatm_prev) begins(iatm) = ipair
            iatm_prev = iatm
        end do

        do iatm=1, natom
            if (begins(iatm) == 0) then
                ends(iatm) = 0
            else
                ends(iatm) = begins(iatm) + nums(iatm) - 1
            end if
        end do

    end subroutine 

end module

