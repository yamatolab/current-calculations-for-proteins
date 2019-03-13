module utils
    implicit none

    integer :: iatm, jatm, igrp, jgrp, itar, jtar
    real(8) :: f_ij(3), r_ij(3), f_ij_r_ij(3, 3)

contains

    function tensor(r1, r2)
        implicit none
        real(8), intent(in) :: r1(3)
        real(8), intent(in) :: r2(3)
        real(8) :: tensor(3, 3)
        integer :: i, j

        do i=1, 3
        do j=1, 3
            tensor(i, j) = r1(i) * r2(j)
        end do
        end do
    end function

end module

module bonded

    implicit none
    ! input
    integer :: natom, ntarget, ngroup, ntbf
    integer, allocatable :: target_atoms(:)  ! (ntarget)
    integer, allocatable :: iatm_to_itars(:) ! (natom)
    integer, allocatable :: iatm_to_igrps(:) ! (natom)
    integer, allocatable :: bonded_pairs(:, :)  ! (ntbf, 2)

    real(8), allocatable :: current(:, :, :)     ! (ntarget, 3, 3)
    real(8), allocatable :: current_inn(:, :, :) ! (ngroup, 3, 3)
    real(8), allocatable :: current_out(:, :, :) ! (ngroup, 3, 3)

contains

    subroutine initialize(target_atoms_in, iatm_to_igrps_in, &
                        & bonded_pairs_in, &
                        & ntarget_in, natom_in, ntbf_in)

        use utils

        implicit none

        integer, intent(in) :: ntarget_in, natom_in, ntbf_in
        integer, intent(in) :: target_atoms_in(ntarget_in)
        integer, intent(in) :: bonded_pairs_in(ntbf_in, 2)
        integer, intent(in) :: iatm_to_igrps_in(natom_in)

        natom   = natom_in
        ntarget = ntarget_in
        ntbf    = ntbf_in

        if (allocated(target_atoms)) deallocate(target_atoms)
        allocate(target_atoms(ntarget))
        target_atoms = target_atoms_in

        if (allocated(bonded_pairs)) deallocate(bonded_pairs)
        allocate(bonded_pairs(ntbf,2))
        bonded_pairs = bonded_pairs_in

        ! iatm_to_itars(iatm) => itar
        if (allocated(iatm_to_itars)) deallocate(iatm_to_itars)
        allocate(iatm_to_itars(natom))

        iatm_to_itars = 0
        do itar=1, ntarget
            iatm = target_atoms(itar)
            iatm_to_itars(iatm) = itar
        end do

        ! for group
        ngroup = maxval(iatm_to_igrps_in)
        if (allocated(iatm_to_igrps)) deallocate(iatm_to_igrps)
        allocate(iatm_to_igrps(natom))
        iatm_to_igrps = iatm_to_igrps_in

        ! for current
        ! atomic
        if (allocated(current)) deallocate(current)
        allocate(current(ntarget, 3, 3))
        current = 0.0d0

        ! group inner
        if (allocated(current_inn)) deallocate(current_inn)
        allocate(current_inn(ngroup, 3, 3))
        current_inn = 0.0d0

        ! group outer
        if (allocated(current_out)) deallocate(current_out)
        allocate(current_out(ngroup, 3, 3))
        current_out = 0.0d0

    end subroutine

    subroutine cal_bonded(crd, tbforces, volumes, gvolumes, &
                & natom_in, ntbf, ntarget_in, ngroup_in)

        use utils
        implicit none
        integer, intent(in) :: natom_in, ntarget_in, ngroup_in
        integer, intent(in) :: ntbf
        real(8), intent(in) :: crd(natom_in, 3)
        real(8), intent(in) :: tbforces(ntbf, 3)
        real(8), intent(in) :: volumes(ntarget_in), gvolumes(ngroup_in)

        integer :: itbf

        current = 0.0d0
        current_inn = 0.0d0
        current_out = 0.0d0

        do itbf=1, ntbf
            iatm = bonded_pairs(itbf, 1)
            jatm = bonded_pairs(itbf, 2)

            itar = iatm_to_itars(iatm)
            jtar = iatm_to_itars(jatm)

            if ( (itar == 0) .and. (jtar == 0) ) cycle

            f_ij(:) = tbforces(itbf,:)
            r_ij(:) = crd(iatm,:) - crd(jatm,:)
            f_ij_r_ij = tensor(f_ij, r_ij)

            ! sum up for iatm
            if ( itar /= 0 ) then
                current(itar, :, :) = current(itar, :, :) & 
                    & + 0.5 * f_ij_r_ij(:, :) / volumes(itar)
            end if

            ! sum up for jatm
            if ( jtar /= 0 ) then
                current(jtar, :, :) = current(jtar, :, :) & 
                    & + 0.5 * f_ij_r_ij(:, :) / volumes(jtar)
            end if

            ! if group calculation is not applied
            if (ngroup == 0) cycle

            ! sum up atomic currents in region of the groups
            igrp = iatm_to_igrps(iatm)
            jgrp = iatm_to_igrps(jatm)

            if ( (igrp == 0) .and. (jgrp == 0) ) cycle

            ! sum up for inside region of the group
            if ( igrp == jgrp ) then
                ! plus the atomic current two times for one group
                current_inn(igrp, :, :) = current_inn(igrp, :, :) & 
                    & + f_ij_r_ij(:,:) / gvolumes(igrp)
                cycle
            end if

            ! contributions from outside of region for group
            ! sum up for ith group
            if ( igrp /= 0 ) then
                current_out(igrp, :, :) = current_out(igrp, :, :)  & 
                    & + 0.5 * f_ij_r_ij(:, :) / gvolumes(igrp)
            end if

            ! sum up for jth group
            if ( jgrp /= 0 ) then
                current_out(jgrp, :, :) = current_out(jgrp, :, :)  & 
                    & + 0.5 * f_ij_r_ij(:, :) / gvolumes(jgrp)
            end if

        end do

    end subroutine

end module 


module nonbonded

    implicit none
    ! input
    integer :: natom, ntarget, ngroup
    integer, allocatable :: target_atoms(:)  ! (ntarget)
    integer, allocatable :: iatm_to_itars(:) ! (natom)
    integer, allocatable :: iatm_to_igrps(:) ! (natom)

    real(8), allocatable :: current(:, :, :)     ! (ntarget, 3, 3)
    real(8), allocatable :: current_inn(:, :, :) ! (ngroup, 3, 3)
    real(8), allocatable :: current_out(:, :, :) ! (ngroup, 3, 3)

    real(8), allocatable :: crd(:, :) ! (natom, 3)
    real(8), allocatable :: volumes(:), gvolumes(:) ! (ntarget), (ngroup)

contains

    subroutine initialize(target_atoms_in, iatm_to_igrps_in, &
                        & ntarget_in, natom_in)
        use utils

        implicit none

        integer, intent(in) :: ntarget_in, natom_in
        integer, intent(in) :: target_atoms_in(ntarget_in)
        integer, intent(in) :: iatm_to_igrps_in(natom_in)

        natom   = natom_in
        ntarget = ntarget_in

        if (allocated(target_atoms)) deallocate(target_atoms)
        allocate(target_atoms(ntarget))
        target_atoms = target_atoms_in

        ! iatm_to_itars(iatm) => itar
        if (allocated(iatm_to_itars)) deallocate(iatm_to_itars)
        allocate(iatm_to_itars(natom))

        iatm_to_itars = 0
        do itar=1, ntarget
            iatm = target_atoms(itar)
            iatm_to_itars(iatm) = itar
        end do

        ! for group
        ngroup = maxval(iatm_to_igrps_in)
        if (allocated(iatm_to_igrps)) deallocate(iatm_to_igrps)
        allocate(iatm_to_igrps(natom))
        iatm_to_igrps = iatm_to_igrps_in

        ! for current
        ! atomic
        if (allocated(current)) deallocate(current)
        allocate(current(ntarget, 3, 3))
        current = 0.0d0

        ! group inner
        if (allocated(current_inn)) deallocate(current_inn)
        allocate(current_inn(ngroup, 3, 3))
        current_inn = 0.0d0

        ! group outer
        if (allocated(current_out)) deallocate(current_out)
        allocate(current_out(ngroup, 3, 3))
        current_out = 0.0d0

    end subroutine

    subroutine init_cal(crd_in, volumes_in, gvolumes_in, &
            & ntarget_in, natom_in, ngroup_in)

        implicit none

        integer, intent(in) :: ntarget_in, natom_in, ngroup_in
        real(8), intent(in) :: crd_in(natom_in, 3)
        real(8), intent(in) :: volumes_in(ntarget_in), gvolumes_in(ngroup_in)

        ! coordinate
        if (allocated(crd)) deallocate(crd)
        allocate(crd(natom, 3))
        crd = crd_in

        ! volumes
        if (allocated(volumes)) deallocate(volumes)
        allocate(volumes(ntarget))
        volumes = volumes_in

        ! group volumes
        if (allocated(gvolumes)) deallocate(gvolumes)
        allocate(gvolumes(ngroup))
        gvolumes = gvolumes_in

        if ( natom_in /= natom ) then
            print*, " The number of atoms is invalid: ", natom_in, natom
            stop
        end if

        if ( ntarget_in /= ntarget ) then
            print*, " The number of targets is invalid: ", ntarget_in, ntarget
            stop
        end if

        if ( ngroup_in /= ngroup ) then
            print*, " The number of groups is invalid: ", ngroup_in, ngroup
            stop
        end if

        ! current
        current = 0.0d0
        current_inn = 0.0d0
        current_out = 0.0d0

    end subroutine

    subroutine cal_nonbonded(tbforces, interact_table, ninteract, ntbf)

        use utils
        implicit none
        integer, intent(in) :: ninteract, ntbf
        integer, intent(in) :: interact_table(ninteract, 3)
        ! (iatm_beg:iatm_end, iatm_beg:natom, :)
        real(8), intent(in) :: tbforces(ntbf, 3)
        integer :: iatm_beg, iatm_end, jatm_beg, jatm_end, jint, itbf

        itbf = 0
        do jint=1, ninteract
            iatm = interact_table(jint, 1)
            itar = iatm_to_itars(iatm)

            jatm_beg = interact_table(jint, 2)
            jatm_end = interact_table(jint, 3)

            do jatm=jatm_beg, jatm_end
                itbf = itbf + 1
                jtar = iatm_to_itars(jatm)
                if ( (itar == 0) .and. (jtar == 0) ) cycle

                f_ij(:) = tbforces(itbf, :)
                r_ij(:) = crd(iatm, :) - crd(jatm, :)
                f_ij_r_ij = tensor(f_ij, r_ij)

                ! sum up for iatm
                if ( itar /= 0 ) then
                    current(itar, :, :) = current(itar, :, :) &
                        & + 0.5 * f_ij_r_ij(:, :) / volumes(itar)
                end if

                ! sum up for jatm
                if ( jtar /= 0 ) then
                    current(jtar, :, :) = current(jtar, :, :) & 
                        & + 0.5 * f_ij_r_ij(:, :) / volumes(jtar)
                end if

                ! if group calculation is not applied
                if (ngroup == 0) cycle

                ! sum up atomic currents in region of the groups
                igrp = iatm_to_igrps(iatm)
                jgrp = iatm_to_igrps(jatm)

                if ( (igrp == 0) .and. (jgrp == 0) ) cycle

                ! sum up for inner region of the group
                if ( igrp == jgrp ) then
                    ! plus the atomic current two times for one group
                    current_inn(igrp, :, :) = current_inn(igrp, :, :) & 
                        & + f_ij_r_ij(:,:) / gvolumes(igrp)
                    cycle
                end if

                ! contributions from outside of region for group
                ! sum up for ith group
                if ( igrp /= 0 ) then
                    current_out(igrp, :, :) = current_out(igrp, :, :)  & 
                        & + 0.5 * f_ij_r_ij(:, :) / gvolumes(igrp)
                end if

                ! sum up for jth group
                if ( jgrp /= 0 ) then
                    current_out(jgrp, :, :) = current_out(jgrp, :, :)  & 
                        & + 0.5 * f_ij_r_ij(:, :) / gvolumes(jgrp)
                end if

            end do

        end do

    end subroutine

end module


module kinetic

    implicit none
    ! input
    integer :: natom, ntarget, ngroup
    integer, allocatable :: target_atoms(:)  ! (ntarget)
    integer, allocatable :: iatm_to_itars(:) ! (natom)
    integer, allocatable :: iatm_to_igrps(:) ! (natom)

    real(8), allocatable :: current(:, :, :)     ! (ntarget, 3, 3)
    real(8), allocatable :: current_inn(:, :, :) ! (ngroup, 3, 3)
    real(8), allocatable :: current_out(:, :, :) ! (ngroup, 3, 3)

contains

    subroutine initialize(target_atoms_in, iatm_to_igrps_in, &
                        & ntarget_in, natom_in)

        use utils

        implicit none

        integer, intent(in) :: ntarget_in, natom_in
        integer, intent(in) :: target_atoms_in(ntarget_in)
        integer, intent(in) :: iatm_to_igrps_in(natom_in)

        natom   = natom_in
        ntarget = ntarget_in

        if (allocated(target_atoms)) deallocate(target_atoms)
        allocate(target_atoms(ntarget))
        target_atoms = target_atoms_in

        ! iatm_to_itars(iatm) => itar
        if (allocated(iatm_to_itars)) deallocate(iatm_to_itars)
        allocate(iatm_to_itars(natom))

        iatm_to_itars = 0
        do itar=1, ntarget
            iatm = target_atoms(itar)
            iatm_to_itars(iatm) = itar
        end do

        ! for group
        ngroup = maxval(iatm_to_igrps_in)
        if (allocated(iatm_to_igrps)) deallocate(iatm_to_igrps)
        allocate(iatm_to_igrps(natom))
        iatm_to_igrps = iatm_to_igrps_in

        ! for current
        ! atomic
        if (allocated(current)) deallocate(current)
        allocate(current(ntarget, 3, 3))
        current = 0.0d0

        ! group inner
        if (allocated(current_inn)) deallocate(current_inn)
        allocate(current_inn(ngroup, 3, 3))
        current_inn = 0.0d0

        ! group outer
        if (allocated(current_out)) deallocate(current_out)
        allocate(current_out(ngroup, 3, 3))
        current_out = 0.0d0

    end subroutine

    subroutine cal_kinetic(vel, masses, volumes, gvolumes, &
                & natom_in, ntarget_in, ngroup_in)

        use utils

        implicit none
        integer, intent(in) :: natom_in, ntarget_in, ngroup_in
        real(8), intent(in) :: vel(natom_in, 3), masses(natom_in)
        real(8), intent(in) :: volumes(ntarget_in), gvolumes(ngroup_in)

        real(8) :: vel2_mass(3, 3)

        ! initialization
        current = 0.0d0
        current_inn = 0.0d0
        current_out = 0.0d0

        do itar=1, ntarget
            iatm = target_atoms(itar)

            vel2_mass(:, :) = tensor( vel(iatm,:), vel(iatm,:)*masses(iatm) )

            current(itar, :, :) = vel2_mass(:, :) / volumes(itar)

            if (ngroup == 0) cycle

            ! sum up atomic currents in region of the groups
            igrp = iatm_to_igrps(iatm)

            ! sum up for ith group
            if ( igrp /= 0 ) then
                current_inn(igrp, :, :) = current_inn(igrp, :, :)  & 
                    & + vel2_mass(:, :) / gvolumes(igrp)
            end if

        end do

    end subroutine

end module 

