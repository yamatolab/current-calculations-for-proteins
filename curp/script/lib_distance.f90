module distance

    implicit none
    integer, allocatable :: res_beg_end_pairs(:, :) ! (nres, 2)
    real(8) :: cutoff2
    integer :: method ! 1: cog, 2: nearest, 3: farthest

    ! common variables
    integer :: iatm, jatm
    real(8) :: r_ij(3), l2

contains

    subroutine cal_dist2s(dist2s, crd, nres, natom)
        
        implicit none

        integer, intent(in) :: natom, nres
        real(8), intent(in) :: crd(natom, 3)
        real(8), intent(out):: dist2s(nres, nres)

        integer :: ires, jres, ibeg, iend, jbeg, jend
        
        dist2s = 0.0d0

        do ires=1, nres-1
            ibeg = res_beg_end_pairs(ires, 1)
            iend = res_beg_end_pairs(ires, 2)

            do jres=ires+1, nres
                jbeg = res_beg_end_pairs(jres, 1)
                jend = res_beg_end_pairs(jres, 2)

                if (method == 1) then
                    dist2s(ires,jres) = & 
                        & cog_dist2(crd(ibeg:iend,:), crd(jbeg:jend,:))

                else if (method == 2) then
                    dist2s(ires,jres) = & 
                        & nearest_dist2(crd(ibeg:iend,:), crd(jbeg:jend,:))

                else if (method == 3) then
                    dist2s(ires,jres) = & 
                        & farthest_dist2(crd(ibeg:iend,:), crd(jbeg:jend,:))
                end if

            end do

        end do

    end subroutine

    real(8) function cog_dist2(crd_i, crd_j)

        implicit none
        real(8), intent(in) :: crd_i(:,:), crd_j(:,:)

        r_ij(:) = sum(crd_i,1)/size(crd_i,1) - sum(crd_j,1)/size(crd_j,1)
        l2 = dot_product(r_ij, r_ij)
        cog_dist2 = l2

    end function

    real(8) function nearest_dist2(crd_i, crd_j)

        implicit none
        real(8), intent(in) :: crd_i(:,:), crd_j(:,:)

        nearest_dist2 = cutoff2

        do iatm=1, size(crd_i, 1)
            do jatm=1, size(crd_j, 1)

                r_ij(:) = crd_i(iatm,:) - crd_j(jatm,:)
                l2 = dot_product(r_ij, r_ij)
                if (l2 <= nearest_dist2) then
                    nearest_dist2 = l2
                end if

            end do
        end do

    end function

    real(8) function farthest_dist2(crd_i, crd_j)

        implicit none
        real(8), intent(in) :: crd_i(:,:), crd_j(:,:)

        farthest_dist2 = cutoff2

        do iatm=1, size(crd_i, 1)
            do jatm=1, size(crd_j, 1)

                r_ij(:) = crd_i(iatm,:) - crd_j(jatm,:)
                l2 = dot_product(r_ij, r_ij)
                if (l2 >= farthest_dist2) then
                    farthest_dist2 = l2
                end if

            end do
        end do

    end function 

end module
