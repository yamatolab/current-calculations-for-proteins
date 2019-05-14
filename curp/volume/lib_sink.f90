subroutine sink_nohyd(solv_crd, natom_solv, base_crd, shell_crd, &
    & probe_length, natom_base, natom_shell)

    implicit none

    ! arguments variables
    integer, intent(in) :: natom_base, natom_shell
    real(8), intent(in) :: probe_length
    real(8), intent(in) :: base_crd(natom_base, 3)
    real(8), intent(in) :: shell_crd(natom_shell, 3)

    ! return variables
    integer, intent(out):: natom_solv
    real(8), intent(out):: solv_crd(natom_shell+100, 3)

    ! internal variables
    logical :: remove_list(natom_shell)
    real(8) :: r_ij(3), l_ij2, lcut2
    integer :: iatm, jatm


    ! initialization
    solv_crd = 0.0
    remove_list = .false. 
    lcut2 = probe_length*probe_length

    ! pick up the water atoms to be removed by the probe probe length
    do iatm=1, natom_shell

        do jatm=1, natom_base

            r_ij(:) = shell_crd(iatm,:) - base_crd(jatm,:)
            l_ij2 = r_ij(1)**2 + r_ij(2)**2 + r_ij(3)**2

            if ( l_ij2 < lcut2 ) then
                remove_list(iatm) = .true.
                exit
            end if

        end do

    enddo

    ! make new coordinate for solvation
    do iatm=1, natom_base
        solv_crd(iatm, :) = base_crd(iatm, :)
    end do

    iatm = 0
    do jatm=1, natom_shell
        if ( remove_list(jatm) ) cycle

        iatm = iatm + 1
        solv_crd(natom_base+iatm, :) = shell_crd(jatm, :)

    end do

    natom_solv = natom_base + iatm

end subroutine

subroutine sink(solv_crd, natom_solv, base_crd, shell_crd, &
    & probe_length, natom_base, natom_shell)

    implicit none

    ! arguments variables
    integer, intent(in) :: natom_base, natom_shell
    real(8), intent(in) :: probe_length
    real(8), intent(in) :: base_crd(natom_base, 3)
    real(8), intent(in) :: shell_crd(natom_shell, 3)

    ! return variables
    integer, intent(out):: natom_solv
    real(8), intent(out):: solv_crd(natom_shell+100, 3)

    ! internal variables
    logical :: remove_list(natom_shell)
    real(8) :: r_ij(3), l_ij2, lcut2
    integer :: iatm, jatm, ires_1


    ! initialization
    solv_crd = 0.0
    remove_list = .false. 
    lcut2 = probe_length*probe_length

    ! pick up the water atoms to be removed by the probe probe length

    do iatm=1, natom_shell

        do jatm=1, natom_base

            r_ij(:) = shell_crd(iatm,:) - base_crd(jatm,:)
            l_ij2 = r_ij(1)**2 + r_ij(2)**2 + r_ij(3)**2

            if ( l_ij2 < lcut2 ) then
                ires_1 = (iatm-1) / 3

                remove_list(3*(ires_1)+1) = .true.
                remove_list(3*(ires_1)+2) = .true.
                remove_list(3*(ires_1)+3) = .true.
                exit
            end if

        end do

    enddo

    ! make new coordinate for solvation
    do iatm=1, natom_base
        solv_crd(iatm, :) = base_crd(iatm, :)
    end do

    iatm = 0
    do jatm=1, natom_shell
        if ( remove_list(jatm) ) cycle

        iatm = iatm + 1
        solv_crd(natom_base+iatm, :) = shell_crd(jatm, :)

    end do

    natom_solv = natom_base + iatm

end subroutine
