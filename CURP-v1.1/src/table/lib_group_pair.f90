module within_gpair

    implicit none
    integer :: npair
    integer, allocatable :: gpair_table(:, :) ! (npair, 3)
    integer, allocatable :: begins(:) ! (natom)
    integer, allocatable :: ends(:)   ! (natom)

    integer :: iatm, jatm
    integer :: i, ipair

contains

    subroutine setup(gpair_table_in, natom, npair_in)
        implicit none

        integer, intent(in) :: npair_in, natom
        integer, intent(in) :: gpair_table_in(npair_in, 3)

        npair = npair_in

        allocate(gpair_table(npair_in, 3))
        gpair_table = gpair_table_in
        allocate(begins(natom))
        allocate(ends(natom))

        call get_begins_ends(begins, ends, gpair_table, npair, natom)

    end subroutine

    logical function is_within_gpair(iatm, jatm)
        implicit none
        integer, intent(in) :: iatm, jatm
        integer :: jatm_beg, jatm_end

        is_within_gpair = .false.
        do ipair=begins(iatm), ends(iatm)
            jatm_beg = gpair_table(ipair, 2)
            jatm_end = gpair_table(ipair, 3)

            if ((jatm_beg<=jatm) .and. (jatm<=jatm_end)) then
                is_within_gpair = .true.
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

                if (.not.is_within_gpair(iatm, jatm)) then

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

    subroutine get_begins_ends(begins, ends, gpair_table, npair, natom)

        implicit none
        integer, intent(in) :: npair, natom
        integer, intent(in) :: gpair_table(npair, 3)

        integer, intent(out):: begins(natom)
        integer, intent(out):: ends(natom)

        integer :: nums(natom)

        integer :: ipair, iatm, jatm, iatm_prev

        begins(:) = 0
        nums(:)   = 0

        iatm_prev = 0
        do ipair=1, npair
            iatm     = gpair_table(ipair, 1)
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


