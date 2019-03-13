module limit
    implicit none

    integer, parameter :: maxcan=1000, maxver=500

end module

module voronoi_one

    use limit

    implicit none

    ! input variables
    logical, allocatable :: is_enable_atoms(:) ! (natom)
    real(8), allocatable :: crd(:, :)
    integer :: natom
    real(8) :: cutoff, box

    ! voronoi analysis informations obtained
    integer :: nface, ncan, nedge, nver
    integer :: can_to_iatm(maxcan), edges(maxcan)
    real(8) :: can_crd(maxcan, 3), can_len2(maxcan)
    real(8) :: ver_crd(maxver, 3)
    integer :: ver_to_3cans(maxver,3)

    ! neighbour list
    integer :: nnab, nab_list(maxver)

contains

    subroutine initialize(is_enable_atoms_in, cutoff_in, box_in, natom_in)
        implicit none
        integer, intent(in) :: natom_in   ! the number of atoms
        logical, intent(in) :: is_enable_atoms_in(natom_in)
        real(8), intent(in) :: cutoff_in, box_in

        ! initialize
        natom = natom_in
        cutoff = cutoff_in
        box = box_in

        ! allocate and initialize
        if ( allocated(is_enable_atoms) ) deallocate(is_enable_atoms)
        allocate(is_enable_atoms(natom))
        is_enable_atoms = is_enable_atoms_in

        ! allocate only
        if ( allocated(crd) ) deallocate(crd)
        allocate(crd(natom, 3))

    end subroutine

end module


subroutine set_crd(crd_in, natom_in)

    use voronoi_one

    integer, intent(in) :: natom_in
    real(8), intent(in) :: crd_in(natom_in, 3)

    crd = crd_in

end subroutine


subroutine cal_voronoi_one(volume, iatm)

! *******************************************************************
! ** one of all main loop
! *******************************************************************

    use limit
    use voronoi_one

    implicit none

    ! return variables
    real(8), intent(out):: volume

    ! arguments variables
    integer :: iatm

    integer :: itar, jatm
    real(8) :: r_ij(3), l_ij2, lcut2
    integer :: ican, jcan, kcan, iver

    ! **************************************************

    nnab = 0
    nab_list = 0

    ! select candidates
    call get_candidates( iatm )

    ! sort into ascending order of distance 
    ! TODO: should improve efficiency
    call sort_dist(can_crd, can_len2, can_to_iatm, ncan)

    ! perform voronoi analysis
    call ana_voronoi( volume )

    ! summary
    nnab = 0
    nab_list = 0
    do ican=1, ncan
        if ( edges(ican) == 0 ) cycle

        nnab = nnab + 1
        nab_list(nnab) = can_to_iatm(ican)

    end do

end subroutine

subroutine get_candidates(iatm)

    use limit
    use voronoi_one

    implicit none

    ! argument variables
    integer :: iatm

    ! other variables
    integer :: jatm, ican, jcan
    real(8) :: pos_i(3), r_ij(3), l_ij2, lcut2

    lcut2 = cutoff*cutoff

    pos_i(:) = crd(iatm, :)
    ican = 0
    do jatm=1, natom
        if ( iatm == jatm ) cycle
        if ( .not. is_enable_atoms(jatm) ) cycle

        r_ij(:) = crd(iatm, :) - crd(jatm, :)
        r_ij(:) = r_ij(:) - anint(r_ij(:) / box) * box
        l_ij2 = r_ij(1)**2 + r_ij(2)**2 + r_ij(3)**2

        if (l_ij2 < lcut2) then

            ican = ican + 1

            if (ican > maxcan) then
                print*, ' too many candidates, ', ican
                stop
            endif

            can_crd(ican,:)   = r_ij(:)
            can_len2(ican)    = l_ij2
            can_to_iatm(ican) = jatm

        endif

    enddo

    ! candidates have been selected
    ncan = ican

end subroutine

subroutine ana_voronoi(volume)

! *******************************************************************
! ** routine to perform voronoi analysis                           **
! ** and calculate volume per polyhedron                           **
! **                                                               **
! ** We work initially on double the correct scale,                **
! ** i.e. the faces of the polyhedron go through the points.       **
! *******************************************************************

    use limit
    use voronoi_one

    implicit none

    ! output arguments
    real(8) :: volume

    logical :: ok
    integer :: iver, ican, jcan, kcan, lcan
    real(8) :: ai, bi, ci, di, aj, bj, cj, dj, ak, bk, ck, dk
    real(8) :: ab, bc, ca, da, db, dc, det
    real(8) :: vxijk, vyijk, vzijk
    real(8), parameter :: tol = 1.e-6

    real(8) :: face_volume
    integer :: faces_to_vertices(maxcan, 15)

    ! *******************************************************************

    ! ** If there are less than 4 points given, 
    ! ** We cannot construct a polyhedron

    if ( ncan < 4 ) then
        print*, ' less than 4 points given to voronoi analysis ', ncan
        stop
    endif

    iver = 0

    ! ** We aim to examine each possible vertex defined by the intersection
    ! ** of 3 planes each plane is specified by can_crd, can_len2

    do ican=1, ncan-2 ! 400

        ai = can_crd(ican, 1)
        bi = can_crd(ican, 2)
        ci = can_crd(ican, 3)
        di = -can_len2(ican)

        do jcan=ican+1, ncan-1 ! 300

            aj = can_crd(jcan, 1)
            bj = can_crd(jcan, 2)
            cj = can_crd(jcan, 3)
            dj = -can_len2(jcan)

            ab = ai * bj - aj * bi
            bc = bi * cj - bj * ci
            ca = ci * aj - cj * ai
            da = di * aj - dj * ai
            db = di * bj - dj * bi
            dc = di * cj - dj * ci

            do kcan=jcan+1, ncan ! 200

                ak = can_crd(kcan, 1)
                bk = can_crd(kcan, 2)
                ck = can_crd(kcan, 3)
                dk = -can_len2(kcan)

                det = ak * bc + bk * ca + ck * ab

                if ( abs(det) > tol ) then

                    ! the planes intersect

                    vxijk = ( - dk * bc + bk * dc - ck * db ) / det
                    vyijk = ( - ak * dc - dk * ca + ck * da ) / det
                    vzijk = (   ak * db - bk * da - dk * ab ) / det

                    ! now we take shots at the vertex
                    ! using the remaining planes ...

                    ok = .true.
                    do lcan=1, ncan
                        if ( .not. ok) exit

                        if ( ( lcan /= ican ) .and. &
                           & ( lcan /= jcan ) .and. &
                           & ( lcan /= kcan )       ) then

                             oK = ( ( can_crd(lcan,1) * vxijk &
                                  & + can_crd(lcan,2) * vyijk &
                                  & + can_crd(lcan,3) * vzijk ) & 
                                  & <= can_len2(lcan) )

                        endif

                    enddo

                    ! if the vertex made it add it to the hall of fame
                    ! convert to correct scale

                    if ( ok ) then
                        iver = iver + 1
                        if ( iver > maxver ) stop 'too many vertices'

                        ver_to_3cans(iver, 1) = ican
                        ver_to_3cans(iver, 2) = jcan
                        ver_to_3cans(iver, 3) = kcan
                        ver_crd(iver, 1) = 0.5 * vxijk
                        ver_crd(iver, 2) = 0.5 * vyijk
                        ver_crd(iver, 3) = 0.5 * vzijk

                    endif

                endif

            enddo

        enddo

    enddo

    nver = iver

    if ( nver  < 4 ) then
        print*, ' less than 4 vertices fount in work ', nver
        stop
    endif

    ! ** Identify neighbouring points and
    ! ** gather edges(vertices) index for each faces(polygon or plane). 

    faces_to_vertices = 0
    edges = 0

    do iver = 1, nver

        ican = ver_to_3cans(iver, 1)
        jcan = ver_to_3cans(iver, 2)
        kcan = ver_to_3cans(iver, 3)

        edges(ican) = edges(ican) + 1
        edges(jcan) = edges(jcan) + 1
        edges(kcan) = edges(kcan) + 1

        faces_to_vertices( ican, edges(ican) ) = iver
        faces_to_vertices( jcan, edges(jcan) ) = iver
        faces_to_vertices( kcan, edges(kcan) ) = iver

    enddo

    ! sort vertices in order.
    do ican=1, ncan
        call sort_vertices(faces_to_vertices(ican,:), & 
                & edges(ican), ver_crd, nver)
    enddo

    ! accumulate the volumes for all polygons
    volume = 0.0
    do ican=1, ncan
        call cal_volume(face_volume, & 
             & faces_to_vertices(ican,:), edges(ican), ver_crd, nver)
        volume = volume + face_volume
    enddo

    ! points with nonzero edges are neighbours

    ! check euler relation

    nface = 0
    nedge = 0

    do ican=1, ncan
       if ( edges(ican) > 0 ) nface = nface + 1
       nedge = nedge + edges(ican)
    enddo

    if ( mod(nedge, 2) /= 0 ) then
        print*, 'The number of edges in non-integer ', nedge
        stop
    endif

    if ( ( nver - nedge/2 + nface ) /= 2 ) then
        print*,' **** EULER ERROR: DEGENERACY ? **** '
    endif

end subroutine



subroutine sort_dist(pcrd, plen2, tag, ncan)

!    *******************************************************************
!    ** ROUTINE TO SORT NEIGHBOURS INTO INCREASING ORDER OF DISTANCE  **
!    **                                                               **
!    ** FOR SIMPLICITY WE USE A BUBBLE SORT - OK FOR MAXCAN SMALL.    **
!    *******************************************************************

    use limit

    ! input and output
    integer :: ncan
    real(8), intent(inout) :: pcrd(maxcan, 3)
    real(8), intent(inout) :: plen2(maxcan)
    integer, intent(inout) :: tag(maxcan)

    logical :: change
    integer :: itop, ican, tag_i
    real(8) :: r_i(3), l_i

!*******************************************************************

    change = .true.

    do itop=ncan-1, 1, -1
        if (.not. change) exit

        change = .false.

        do ican=1, itop
            if ( plen2(ican) <= plen2(ican+1) ) cycle

            r_i(:) = pcrd(ican,:)
            l_i = plen2(ican)
            tag_i = tag(ican)

            pcrd(ican,:) = pcrd(ican+1,:)
            plen2(ican) = plen2(ican+1)
            tag(ican) = tag(ican+1)

            pcrd(ican+1,:) = r_i(:)
            plen2(ican+1) = l_i
            tag(ican+1) = tag_i

            change = .true.

        enddo

    enddo
    ! return
end subroutine


subroutine sort_vertices(face_to_vertices, nedge, ver_crd, nver)

! *******************************************************************
! ** Routine to sort vertice in one face into increasing order     **
! ** of distance                                                   **
! *******************************************************************

    use limit

    implicit none

    ! input arguments
    integer :: nedge, nver
    real(8) :: ver_crd(maxver, 3)

    ! input output arguments
    integer :: face_to_vertices(15)

    ! other variables
    integer :: ordered_table(15)
    logical :: defined_edges(15)
    integer :: next_edge, iver, jver, iedge, jedge, can_ver, can_edge
    real(8) :: r_ij(3), dist2, dist2_min

! *******************************************************************

    ! prepare
    ordered_table = 0
    defined_edges = .false.

    ! define first edge
    next_edge = 1
    defined_edges(next_edge) = .true.
    ordered_table(1) = face_to_vertices(next_edge)

    do iedge=2, nedge

        ! define the base edge
        iver = face_to_vertices(next_edge)

        dist2_min = 10000.0d0
        do jedge=2, nedge

            ! eliminate the connection defined
            if (defined_edges(jedge)) cycle

            jver = face_to_vertices(jedge)
            r_ij(:) = ver_crd(iver,:) - ver_crd(jver,:)

            dist2 = r_ij(1)**2 + r_ij(2)**2 + r_ij(3)**2

            if (dist2 < dist2_min) then
                dist2_min = dist2
                can_ver = jver
                can_edge = jedge
            endif
        enddo

        ordered_table(iedge) = can_ver
        next_edge = can_edge
        defined_edges(next_edge) = .true.

    enddo

    face_to_vertices = ordered_table

end subroutine


subroutine cal_volume(face_volume, face_to_vertices, nedge, ver_crd, nver)

! *******************************************************************
! ** Routine to calculate the volume of one face                   **
! *******************************************************************
    use limit

    implicit none

    ! output arguments
    real(8) :: face_volume

    ! input arguments
    integer :: nedge, nver
    integer :: face_to_vertices(15)
    real(8) :: ver_crd(maxver,3)

    ! other variables
    integer :: iver, jver, kver, jedge, kedge
    real(8) :: r_i(3), r_j(3), r_k(3), pyra_volume

    ! *******************************************************************

    face_volume = 0.0d0
    pyra_volume = 0.0d0

    ! define edge 1 as base edge 
    iver = face_to_vertices(1)       
    r_i(:) = ver_crd(iver, :)

    do jedge=2, nedge-1

        jver = face_to_vertices(jedge)
        r_j(:) = ver_crd(jver, :)

        kedge = jedge + 1
        kver = face_to_vertices(kedge)
        r_k(:) = ver_crd(kver, :)

        ! calculate the volume of pyramid
        pyra_volume = r_i(1)*r_j(2)*r_k(3) + r_i(2)*r_j(3)*r_k(1) &
                  & + r_i(3)*r_j(1)*r_k(2) - r_i(3)*r_j(2)*r_k(1) &
                  & - r_i(2)*r_j(1)*r_k(3) - r_i(1)*r_j(3)*r_k(2)

        if (pyra_volume >= 0.0d0) then
            face_volume = face_volume + pyra_volume/6.0d0
        else if (pyra_volume < 0.0d0) then
            face_volume = face_volume - pyra_volume/6.0d0
        endif

    enddo

end subroutine



subroutine cal_voronoi_all(radii, volumes, crd, target_atoms, & 
                         & is_enable_atoms, natom, ntar)

    use limit

    implicit none

    ! arguments
    integer, intent(in) :: natom, ntar
    integer, intent(in) :: target_atoms(ntar)
    logical, intent(in) :: is_enable_atoms(natom)
    real(8), intent(in) :: crd(natom, 3)

    ! return values
    real(8), intent(out):: radii(ntar), volumes(ntar)

    integer :: itar, iatm

    real(8) :: volume

    real(8), parameter :: hyd_volume = 8.0

    ! *******************************************************************
    ! ** MAIN LOOP STARTS                                              **
    ! *******************************************************************

    call set_crd(crd, natom)

    do itar=1, ntar
        iatm = target_atoms(itar)

        if ( is_enable_atoms(iatm) ) then

            call cal_voronoi_one(volume, iatm)
            volumes(itar) = volume

        else
            volumes(itar) = hyd_volume

        end if

    end do

    call cal_radii_sphere(radii, volumes, ntar)

end subroutine

subroutine cal_radii_gauss(radii, volumes, ntar_in)
    implicit none
    integer, intent(in) :: ntar_in
    real(8), intent(in) :: volumes(ntar_in)

    real(8), intent(out) :: radii(ntar_in)

    real(8) :: coef
    real(8), parameter :: pi = 3.141592653589793

    coef = 1.0/pi**(1.0/6.0)
    radii = coef * volumes**(1.0/3.0)

end subroutine

subroutine cal_radii_sphere(radii, volumes, ntar_in)
    implicit none
    integer, intent(in) :: ntar_in
    real(8), intent(in) :: volumes(ntar_in)

    real(8), intent(out) :: radii(ntar_in)

    real(8) :: coef
    real(8), parameter :: pi = 3.141592653589793

    coef = (3.0/(4.0*pi))**(1.0/6.0)
    radii = coef * volumes**(1.0/3.0)

end subroutine
