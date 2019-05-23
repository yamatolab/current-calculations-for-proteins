
subroutine calvolume_org(radii, volumes, crd, well2s, natom)
    implicit none
    integer, intent(in) :: natom
    real(8), intent(in) :: crd(natom, 3), well2s(natom)
    real(8), intent(out):: radii(natom), volumes(natom)

    real(8), parameter :: PI = 3.141592653589793
    integer :: iatm, jatm
    real(8) :: r_ij(3), l_ij, l_ij2, well2, num, den, nearest_l_ij2

    ! initialization
    radii = 0.0d0
    volumes = 0.0d0

    ! calculate radius and volume
    do iatm=1, natom
        well2 = well2s(iatm)
        num = 0.0d0 ! numerator
        den = 0.0d0 ! denominator

        ! judge
        nearest_l_ij2 = 100.0d0

        do jatm=1, natom
            if (iatm == jatm) cycle
            r_ij  = crd(iatm,:) - crd(jatm,:)
            l_ij2 = dot_product(r_ij, r_ij)
            if (l_ij2 < nearest_l_ij2) nearest_l_ij2 = l_ij2

            if (l_ij2 > well2) cycle

            l_ij = sqrt(l_ij2)

            num = num + 1.0d0/l_ij
            den = den + 1.0d0/l_ij2
        enddo

        if (nearest_l_ij2 <= well2) then
            radii(iatm) = 0.5d0 * num / den
        else
            radii(iatm) = 0.5d0 * sqrt(nearest_l_ij2)
        end if

    enddo

    volumes = 4.0d0 * PI * radii**3 / 3.0d0

end subroutine

! the volume of atom that include all atoms
subroutine calvolume(radii, volumes, crd, well2s, natom)
    implicit none
    integer, intent(in) :: natom
    real(8), intent(in) :: crd(natom, 3), well2s(natom)
    real(8), intent(out):: radii(natom), volumes(natom)

    real(8), parameter :: PI = 3.141592653589793
    integer :: iatm, jatm
    real(8) :: r_ij(3), l_ij, l_ij2, well2, num, den, nearest_l_ij2

    ! initialization
    radii = 0.0d0
    volumes = 0.0d0

    ! calculate radius and volume
    do iatm=1, natom
        well2 = well2s(iatm)
        num = 0.0d0 ! numerator
        den = 0.0d0 ! denominator

        ! judge
        nearest_l_ij2 = 100.0d0

        do jatm=1, natom
            if (iatm == jatm) cycle
            r_ij  = crd(iatm,:) - crd(jatm,:)
            l_ij2 = dot_product(r_ij, r_ij)
            if (l_ij2 < nearest_l_ij2) nearest_l_ij2 = l_ij2

            ! if (l_ij2 > well2) cycle

            l_ij = sqrt(l_ij2)

            num = num + 1.0d0/l_ij
            den = den + 1.0d0/l_ij2
        enddo

        radii(iatm) = 0.5d0 * num / den

    enddo

    volumes = 4.0d0 * PI * radii**3 / 3.0d0

end subroutine
