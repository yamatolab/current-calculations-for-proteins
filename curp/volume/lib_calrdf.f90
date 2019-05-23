
subroutine calrdf(rdfs, crd, num_rs, rmax, dr, per_area, natom)
    implicit none
    integer, intent(in) :: num_rs, natom
    real(8), intent(in) :: crd(natom, 3), rmax, dr
    logical, intent(in) :: per_area
    real(8), intent(out):: rdfs(natom, num_rs)

    integer :: iatm, jatm, ir
    real(8) :: r_ij(3), l_ij, r

    ! initialization
    rdfs = 0.0d0

    ! from here !
    do iatm=1, natom
        do jatm=1, natom
            if (iatm == jatm) cycle
            r_ij = crd(iatm,:) - crd(jatm,:)
            l_ij = sqrt(dot_product(r_ij, r_ij))
            if (l_ij >= rmax) cycle

            ! get r index
            ir = nint(l_ij/dr)

            ! increment for rdf on iatm
            rdfs(iatm, ir) = rdfs(iatm, ir) + 1.0

        enddo
    enddo

    if (per_area) then
        do ir=1, num_rs
            r = dr * ir
            rdfs(:, ir) = rdfs(:, ir) / (r*r)
        enddo
    endif

end subroutine
