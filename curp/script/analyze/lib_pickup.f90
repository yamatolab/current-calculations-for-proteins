module pickup

    implicit none
    integer, allocatable :: res_beg_end_pairs(:, :) ! (nres, 2)
    real(8) :: cutoff2
    integer :: cutoff_method ! 1: com, 2: nearest, 3: farthest

contains

    subroutine get_candidate_pairs(candidate, crd, nres, natom)
        
        implicit none

        integer, intent(in) :: natom, nres
        real(8), intent(in) :: crd(natom, 3)
        logical, intent(out):: candidate(nres, nres)

        integer :: ires, jres, ibeg, iend, jbeg, jend


        candidate = .false.

        do ires=1, nres-1
            ibeg = res_beg_end_pairs(ires, 1)
            iend = res_beg_end_pairs(ires, 2)

            do jres=ires+1, nres
                jbeg = res_beg_end_pairs(jres, 1)
                jend = res_beg_end_pairs(jres, 2)

                if (cutoff_method == 1) then
                    if ( is_com(crd(ibeg:iend,:), crd(jbeg:jend,:)) ) then
                        candidate(ires, jres) = .true.
                    end if
                else if (cutoff_method == 2) then
                    if ( is_nearest(crd(ibeg:iend,:), crd(jbeg:jend,:)) ) then
                        candidate(ires, jres) = .true.
                    end if
                else if (cutoff_method == 3) then
                    if ( is_farthest(crd(ibeg:iend,:), crd(jbeg:jend,:)) ) then
                        candidate(ires, jres) = .true.
                    end if
                end if

            end do

        end do

    end subroutine

    logical function is_com(crd_i, crd_j)

        implicit none
        real(8), intent(in) :: crd_i(:,:), crd_j(:,:)

        real(8) :: rcom_ij(3)

        rcom_ij(:) = sum(crd_i,1)/size(crd_i,1) - sum(crd_j,1)/size(crd_j,1)
        is_com =  dot_product(rcom_ij, rcom_ij) <= cutoff2

    end function

    logical function is_nearest(crd_i, crd_j)

        implicit none
        real(8), intent(in) :: crd_i(:,:), crd_j(:,:)

        integer :: iatm, jatm
        real(8) :: r_ij(3)

        is_nearest = .false.
        do iatm=1, size(crd_i, 1)
            do jatm=1, size(crd_j, 1)

                r_ij(:) = crd_i(iatm,:) - crd_j(jatm,:)
                if (dot_product(r_ij, r_ij) <= cutoff2) then
                    is_nearest = .true.
                    exit
                end if

            end do
            if (is_nearest) exit

        end do

    end function

    logical function is_farthest(crd_i, crd_j)

        implicit none
        real(8), intent(in) :: crd_i(:,:), crd_j(:,:)

        integer :: iatm, jatm
        real(8) :: r_ij(3)

        is_farthest = .true.
        do iatm=1, ubound(crd_i, 1)
            do jatm=1, ubound(crd_j, 1)

                r_ij(:) = crd_i(iatm,:) - crd_j(jatm,:)
                if (dot_product(r_ij, r_ij) >= cutoff2) then
                    is_farthest = .false.
                    exit
                end if

            end do
            if (.not.is_farthest) exit

        end do

    end function 

end module


subroutine is_com(ok, crd_i, crd_j, cutoff2, natom_i, natom_j)

    implicit none
    integer, intent(in) :: natom_i, natom_j
    real(8), intent(in) :: crd_i(natom_i,3), crd_j(natom_j,3), cutoff2
    logical, intent(out):: ok

    real(8) :: rcom_ij(3)

    rcom_ij(:) = sum(crd_i,1)/size(crd_i,1) - sum(crd_j,1)/size(crd_j,1)
    ok =  dot_product(rcom_ij, rcom_ij) <= cutoff2

end subroutine 

subroutine is_nearest(ok, crd_i, crd_j, cutoff2, natom_i, natom_j)

    implicit none
    integer, intent(in) :: natom_i, natom_j
    real(8), intent(in) :: crd_i(natom_i,3), crd_j(natom_j,3), cutoff2
    logical, intent(out):: ok

    integer :: iatm, jatm
    real(8) :: r_ij(3)

    ok = .false.
    do iatm=1, natom_i
        do jatm=1, natom_j

            r_ij(:) = crd_i(iatm,:) - crd_j(jatm,:)
            if (dot_product(r_ij, r_ij) <= cutoff2) then
                ok = .true.
                exit
            end if

        end do
        if (ok) exit

    end do

end subroutine 

subroutine is_farthest(ok, crd_i, crd_j, cutoff2, natom_i, natom_j)

    implicit none
    integer, intent(in) :: natom_i, natom_j
    real(8), intent(in) :: crd_i(natom_i,3), crd_j(natom_j,3), cutoff2
    logical, intent(out):: ok

    integer :: iatm, jatm
    real(8) :: r_ij(3)

    ok = .true.
    do iatm=1, natom_i
        do jatm=1, natom_j

            r_ij(:) = crd_i(iatm,:) - crd_j(jatm,:)
            if (dot_product(r_ij, r_ij) >= cutoff2) then
                ok = .false.
                exit
            end if

        end do
        if (.not.ok) exit

    end do

end subroutine 

