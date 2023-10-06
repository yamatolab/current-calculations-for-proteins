!
! CURP 1.2: Yamato, 2021. Minor modification for intra-residue heat flux.
!

module utils
    implicit none

    integer :: iatm, jatm, igrp, jgrp, itar, jtar
    real(8) :: f_ij(3), v_ij(3), d_ij(3), hflux_ij(3)

end module

module bonded

    implicit none
    ! input
    integer :: natom, ntarget, ngroup, ntbf
    integer, allocatable :: target_atoms(:)  ! (ntarget)
    integer, allocatable :: iatm_to_itar(:) ! (natom)
    integer, allocatable :: iatm_to_igrp(:) ! (natom)
    integer, allocatable :: bonded_pairs(:, :)  ! (ntbf, 2)
    logical :: flag_atom, flag_group

    real(8), allocatable :: hflux_atm(:, :, :) ! (ntarget, ntarget, 3)
    real(8), allocatable :: hflux_grp(:, :, :) ! (ngroup, ngroup, 3)

contains

    subroutine initialize(target_atoms_in, iatm_to_igrp_in, & 
            & bonded_pairs_in, &
            & flag_atom_in, flag_group_in, ntarget_in, natom_in, ntbf_in)

        use utils

        implicit none

        integer, intent(in) :: ntarget_in, natom_in, ntbf_in
        integer, intent(in) :: target_atoms_in(ntarget_in)
        integer, intent(in) :: bonded_pairs_in(ntbf_in, 2)
        integer, intent(in) :: iatm_to_igrp_in(natom_in)
        logical, intent(in) :: flag_atom_in, flag_group_in

        natom   = natom_in
        ntarget = ntarget_in
        ntbf    = ntbf_in

        flag_atom  = flag_atom_in
        flag_group = flag_group_in

        if (allocated(target_atoms)) deallocate(target_atoms)
        allocate(target_atoms(ntarget))
        target_atoms = target_atoms_in

        if (allocated(bonded_pairs)) deallocate(bonded_pairs)
        allocate(bonded_pairs(ntbf,2))
        bonded_pairs = bonded_pairs_in
        ! iatm_to_itar(iatm) => itar
        if (allocated(iatm_to_itar)) deallocate(iatm_to_itar)
        allocate(iatm_to_itar(natom))

        iatm_to_itar = 0
        do itar=1, ntarget
            iatm = target_atoms(itar)
            iatm_to_itar(iatm) = itar
        end do

        ! for group
        ngroup = maxval(iatm_to_igrp_in)
        if (allocated(iatm_to_igrp)) deallocate(iatm_to_igrp)
        allocate(iatm_to_igrp(natom))
        iatm_to_igrp = iatm_to_igrp_in

        ! for hflux
        ! atomic
        if ( flag_atom ) then
            if (allocated(hflux_atm)) deallocate(hflux_atm)
            allocate(hflux_atm(ntarget, ntarget, 3))
            hflux_atm = 0.0d0
        end if 

        ! group
        if ( flag_group ) then
            if (allocated(hflux_grp)) deallocate(hflux_grp)
            allocate(hflux_grp(ngroup, ngroup, 3))
            hflux_grp = 0.0d0
        end if
    end subroutine

    subroutine cal_bonded(vel, tbforces, displacement, natom_in, ntbf)

        use utils
        implicit none
        integer, intent(in) :: natom_in
        integer, intent(in) :: ntbf
        real(8), intent(in) :: vel(natom_in, 3)
        real(8), intent(in) :: tbforces(ntbf, 3)
        real(8), intent(in) :: displacement(ntbf, 3)

        integer :: itbf
        if ( flag_atom  ) hflux_atm = 0.0d0
        if ( flag_group ) hflux_grp = 0.0d0
        do itbf=1, ntbf
            iatm = bonded_pairs(itbf, 1)
            jatm = bonded_pairs(itbf, 2)
            itar = iatm_to_itar(iatm)
            jtar = iatm_to_itar(jatm)
            if ( (itar==0) .or. (jtar==0) ) cycle

            f_ij(:) = tbforces(itbf, :)

            v_ij(:) = 0.5 * ( vel(iatm,:) + vel(jatm,:) )
            d_ij(:) = displacement(itbf, :)

            hflux_ij(:) = d_ij(:) * (f_ij(1)*v_ij(1) + f_ij(2)*v_ij(2) + f_ij(3)*v_ij(3))

            ! sum up for ith target and jth target
            if ( flag_atom ) then
                hflux_atm(itar, jtar, :) = hflux_atm(itar, jtar, :) + hflux_ij(:)
                hflux_atm(jtar, itar, :) = hflux_atm(jtar, itar, :) + hflux_ij(:)
            end if

            ! if group calculation is not applied
            if (ngroup == 0) cycle

            ! sum up atomic currents in region of the groups
            if ( flag_group ) then
                igrp = iatm_to_igrp(iatm)
                jgrp = iatm_to_igrp(jatm)

                if ( (igrp == 0) .or. (jgrp==0) ) cycle

                if ( igrp /= jgrp ) then
                  hflux_grp(igrp, jgrp, :) = hflux_grp(igrp, jgrp, :) + hflux_ij(:)
                  hflux_grp(jgrp, igrp, :) = hflux_grp(jgrp, igrp, :) + hflux_ij(:)
                else 
                  hflux_grp(igrp, jgrp, :) = hflux_grp(igrp, jgrp, :) + hflux_ij(:)
                end if

            end if

        end do

    end subroutine

end module


module nonbonded

    implicit none
    ! input
    integer :: natom, ntarget, ngroup
    integer, allocatable :: target_atoms(:)  ! (ntarget)
    integer, allocatable :: iatm_to_itar(:) ! (natom)
    integer, allocatable :: iatm_to_igrp(:) ! (natom)
    logical :: flag_atom, flag_group

    real(8), allocatable :: hflux_atm(:, :, :) ! (ntarget, ntarget)
    real(8), allocatable :: hflux_grp(:, :, :) ! (ngroup, ngroup)

    real(8), allocatable :: vel(:, :) ! (natom, 3)

contains

    subroutine initialize(target_atoms_in, iatm_to_igrp_in, &
                        & flag_atom_in, flag_group_in, &
                        & ntarget_in, natom_in)
        use utils

        implicit none

        integer, intent(in) :: ntarget_in, natom_in
        integer, intent(in) :: target_atoms_in(ntarget_in)
        integer, intent(in) :: iatm_to_igrp_in(natom_in)
        logical, intent(in) :: flag_atom_in, flag_group_in

        natom   = natom_in
        ntarget = ntarget_in

        flag_atom  = flag_atom_in
        flag_group = flag_group_in

        if (allocated(target_atoms)) deallocate(target_atoms)
        allocate(target_atoms(ntarget))
        target_atoms = target_atoms_in

        ! iatm_to_itar(iatm) => itar
        if (allocated(iatm_to_itar)) deallocate(iatm_to_itar)
        allocate(iatm_to_itar(natom))

        iatm_to_itar = 0
        do itar=1, ntarget
            iatm = target_atoms(itar)
            iatm_to_itar(iatm) = itar
        end do

        ! velocity, allcation only
        if (allocated(vel)) deallocate(vel)
        allocate(vel(natom, 3))

        ! for group
        ngroup = maxval(iatm_to_igrp_in)
        if (allocated(iatm_to_igrp)) deallocate(iatm_to_igrp)
        allocate(iatm_to_igrp(natom))
        iatm_to_igrp = iatm_to_igrp_in

        ! for current
        ! atomic
        if ( flag_atom ) then
            if (allocated(hflux_atm)) deallocate(hflux_atm)
            allocate(hflux_atm(ntarget, ntarget, 3))
            hflux_atm = 0.0d0
        end if 

        ! group inner
        if ( flag_group ) then
            if (allocated(hflux_grp)) deallocate(hflux_grp)
            allocate(hflux_grp(ngroup, ngroup, 3))
            hflux_grp = 0.0d0
        end if

    end subroutine

    subroutine init_cal(vel_in, natom_in)

        implicit none

        integer, intent(in) :: natom_in
        real(8), intent(in) :: vel_in(natom_in, 3)

        ! velocity
        vel = vel_in

        ! initialize hflux values
        if ( flag_atom )  hflux_atm = 0.0d0
        if ( flag_group ) hflux_grp = 0.0d0

    end subroutine

    subroutine cal_nonbonded(tbforces, interact_table, displacement, ninteract, ntbf)

        use utils
        implicit none
        integer, intent(in) :: ninteract, ntbf
        integer, intent(in) :: interact_table(ninteract, 3)
        ! (iatm_beg:iatm_end, iatm_beg:natom, :)
        real(8), intent(in) :: tbforces(ntbf, 3)
        real(8), intent(in) :: displacement(ntbf, 3)
        integer :: jatm_beg, jatm_end, jint, itbf

        itbf = 0
        do jint=1, ninteract
            iatm = interact_table(jint, 1)
            itar = iatm_to_itar(iatm)
            if ( itar == 0 ) cycle

            jatm_beg = interact_table(jint, 2)
            jatm_end = interact_table(jint, 3)

            do jatm=jatm_beg, jatm_end
                itbf = itbf + 1 ! increment tbf pointer
                jtar = iatm_to_itar(jatm)
                if ( jtar == 0 ) cycle

                f_ij(:) = tbforces(itbf, :)
                v_ij(:) = 0.5 * ( vel(iatm,:) + vel(jatm,:) )
                d_ij(:) = displacement(itbf, :)

                hflux_ij = d_ij * (f_ij(1)*v_ij(1) + f_ij(2)*v_ij(2) + f_ij(3)*v_ij(3))

                ! sum up for ith target and jth target
                if ( flag_atom ) then
                    hflux_atm(itar, jtar, :) = hflux_atm(itar, jtar, :) + hflux_ij(:)
                    hflux_atm(jtar, itar, :) = hflux_atm(jtar, itar, :) + hflux_ij(:)
                end if

                ! if group calculation is not applied
                if (ngroup == 0) cycle

                ! sum up atomic currents in region of the groups
                if ( flag_group ) then
                    igrp = iatm_to_igrp(iatm)
                    jgrp = iatm_to_igrp(jatm)

                    if ( (igrp == 0) .or. (jgrp==0) ) cycle

                    if ( igrp /= jgrp ) then
                      hflux_grp(igrp, jgrp, :) = hflux_grp(igrp, jgrp, :) + hflux_ij(:)
                      hflux_grp(jgrp, igrp, :) = hflux_grp(jgrp, igrp, :) + hflux_ij(:)
                    else
                      hflux_grp(igrp, jgrp, :) = hflux_grp(igrp, jgrp, :) + hflux_ij(:)
                    end if

                end if

            end do

        end do

    end subroutine

end module

