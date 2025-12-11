
!###############################################################################
module common_vars
    implicit none

    real(8), parameter :: PI = 3.141592653589793
    integer :: jint, itbf
    integer :: iatm, jatm, katm, latm
    integer :: iatm_beg, iatm_end, jatm_beg, jatm_end
    real(8), dimension(3) :: r_ij, r_ji, r_ik, r_ki, r_il, r_li
    real(8), dimension(3) :: r_jk, r_kj, r_jl, r_lj, r_kl, r_lk
    real(8), dimension(3) :: r_im, r_mi, r_jm, r_mj, r_km, r_mk, r_lm, r_ml
    real(8) :: l_ij, l_ji, l_ik, l_ki, l_il, l_li
    real(8) :: l_jk, l_kj, l_jl, l_lj, l_kl, l_lk
    real(8) :: l_im, l_mi, l_jm, l_mj, l_km, l_mk, l_lm, l_ml
    real(8) :: l_ij_inv
    real(8), dimension(3) :: f_ij, f_ik, f_il, f_im, f_jk, f_jl, f_jm
    real(8), dimension(3) :: f_kl, f_km, f_lm
    real(8), dimension(3) :: f_i, f_j, f_k, f_l, f_m
    real(8) :: ene

contains

    function outer_prod(r1, r2)
        implicit none
        real(8), intent(in) :: r1(3)
        real(8), intent(in) :: r2(3)
        real(8) :: outer_prod(3)

        outer_prod(1) = r1(2)*r2(3) - r1(3)*r2(2)
        outer_prod(2) = r1(3)*r2(1) - r1(1)*r2(3)
        outer_prod(3) = r1(1)*r2(2) - r1(2)*r2(1)
        
    end function

    function check_tbforce(f1, f2)
        implicit none
        real(8), intent(in) :: f1(3)
        real(8), intent(in) :: f2(3)
        real(8) :: eps = 0.001d0, diff
        logical :: check_tbforce
        logical :: bool(3)
        integer :: d

        diff = sqrt(dot_product(f1, f1)) * eps

        do d=1, 3
            if (-eps < f1(d) .and. f1(d) < eps) then
                bool(d) = .true.
            else if ((f1(d)-diff < f2(d)) .and. (f2(d) < f1(d)+diff)) then
                bool(d) = .true.
            else
                bool(d) = .false.
            end if
        enddo

        if (bool(1) .and. bool(2) .and. bool(3)) then
            check_tbforce = .true.
        else
            check_tbforce = .false.
        endif

    end function

end module

!###############################################################################
module total
    implicit none
    logical :: check = .false. ! parameter
    integer :: natom   ! the number of atoms
    integer :: nbonded ! the number of bonded pairs
    real(8), allocatable :: crd(:, :) ! (natom, 3)
    ! itbf => (iatm, jatm)
    integer, allocatable :: bonded_pairs(:, :)    !(nbonded, 2)
    real(8), allocatable :: bonded_tbforces(:, :) !(nbonded, 3)
end module

!###############################################################################
module bond
    implicit none
    ! input
    integer, allocatable :: two_atoms(:,:)  ! (nbond, 2)
    real(8), allocatable :: force_consts(:) ! (nbond)
    real(8), allocatable :: length_eqs(:)   ! (nbond)
    integer, allocatable :: ibnd_to_itbf(:) ! (nbond)
    ! output
    real(8) :: energy
    real(8), allocatable :: forces(:, :)   ! (natom, 3)
    real(8), allocatable :: tbforces(:, :) ! (ntbf,  3)
    real(8), allocatable :: displacement(:, :)! (ntbf, 3)

    ! main variales
    integer :: ibnd, nbond, itbf_ij
    integer :: iatm, jatm

contains
    subroutine calculate()
        use total
        use common_vars
        implicit none
        real(8) :: kb, l_eq

        ! for each bond,
        ! E_bond = K_b (l_ij - l_eq)^2
        
        ! F_i = - 2 K_b (l_ij - l_eq) r_ij/l_ij
        ! F_j = - F_i

        ! F_ij = - 2 K_b (l_ij - l_eq) r_ij/l_ij
        ! F_ji = - F_ij

        ! initialization
        nbond    = size(force_consts)
        energy   = 0.0d0
        forces   = 0.0d0
        tbforces = 0.0d0
        displacement = 0.0d0

        do ibnd=1, nbond
            itbf_ij = ibnd_to_itbf(ibnd) ! => ij
            iatm    = two_atoms(ibnd, 1)
            jatm    = two_atoms(ibnd, 2)

            kb   = force_consts(ibnd)
            l_eq = length_eqs(ibnd)

            r_ij = crd(iatm, :) - crd(jatm, :)
            l_ij = sqrt( dot_product(r_ij, r_ij) )

            ! calculate energy
            ene = kb * (l_ij - l_eq) ** 2
            energy = energy + ene

            ! calculate force
            f_i = - 2.0 * kb * (l_ij - l_eq)/l_ij * r_ij

            ! store force
            forces(iatm, :) = forces(iatm, :) + f_i(:)
            forces(jatm, :) = forces(jatm, :) - f_i(:)

            ! calculate two-body force
            f_ij = f_i

            ! store two-body force and two-body distance vector
            if (itbf_ij > 0) then
                tbforces(itbf_ij,:) = tbforces(itbf_ij,:) + f_ij(:)
                displacement(itbf_ij,:) = r_ij(:)
            else
                tbforces(-itbf_ij,:) = tbforces(-itbf_ij,:) - f_ij(:)
                displacement(-itbf_ij,:) = -r_ij(:)
            end if

            ! check
            if (check) then
                print*, 'TB_CHECK: i, j =', iatm, jatm
                print*, 'TB_CHECK: f_i vs. f_ij =>', check_tbforce(f_i, f_ij)
                print*, 'TB_CHECK:',f_i
                print*, 'TB_CHECK:',f_ij
                print*, 'TB_CHECK: f_j vs. f_ji =>', check_tbforce(-f_i, -f_ij)
                print*, 'TB_CHECK:',-f_i
                print*, 'TB_CHECK:',-f_ij
            end if

        enddo
        
    end subroutine

    subroutine print_tbforce()
        use total, only: bonded_pairs, nbonded
        use common_vars
        implicit none
        real(8) :: tbf_abs

        do itbf=1, nbonded
            iatm = bonded_pairs(itbf, 1)
            jatm = bonded_pairs(itbf, 2)

            tbf_abs = dot_product(tbforces(itbf,:), tbforces(itbf,:))
            if (-0.0001 < tbf_abs .and. tbf_abs < 0.0001) cycle
            write(*, '(2I5,3F15.7)'), iatm, jatm, tbforces(itbf_ij, :)
         end do
            
    end subroutine

    subroutine print_atom_order()
        use total, only: bonded_pairs
        use common_vars
        implicit none
        integer :: iatm_, jatm_

        print*, '** print atom order for bond **'
        nbond = size(force_consts)
        do ibnd=1, nbond
            iatm = two_atoms(ibnd,1)
            jatm = two_atoms(ibnd,2)

            itbf_ij = ibnd_to_itbf(ibnd)
            if (itbf_ij > 0) then
                iatm_ = bonded_pairs(itbf_ij, 1)
                jatm_ = bonded_pairs(itbf_ij, 2)
            else
                iatm_ = bonded_pairs(-itbf_ij, 2)
                jatm_ = bonded_pairs(-itbf_ij, 1)
            end if

            print*, ibnd
            print*, iatm, iatm_
            print*, jatm, jatm_

         end do
            
    end subroutine

end module

!###############################################################################
module angle
    implicit none
    ! input
    integer, allocatable :: three_atoms(:,:)  ! (nangle, 3)
    real(8), allocatable :: force_consts(:)   ! (nangle)
    real(8), allocatable :: theta_eqs(:)      ! (nangle)
    integer, allocatable :: iang_to_itbf(:,:) ! (nangle, 3)
    ! output
    real(8) :: energy
    real(8), allocatable :: forces(:, :)   ! (natom, 3)
    real(8), allocatable :: tbforces(:, :) ! (ntbf, 3)
    real(8), allocatable :: displacement(:, :) ! (ntbf, 3)

    ! main variales
    integer :: iang, nangle, itbf_ij, itbf_ik, itbf_jk
    integer :: iatm, jatm, katm

contains
    subroutine calculate()
        use total
        use common_vars
        implicit none
        real(8) :: ka, theta_eq, cos_theta, theta, coeff, sin_theta
        real(8), parameter :: eps = 1.0d-10

        ! Energy:
        ! E_angle = K_a (theta_ijk - theta_eq)^2

        ! Force:
        ! A = 2 K_a (theta - theta_eq) / sin(theta)
        ! F_i = A/ (r_kj/l_kj )

        ! two-body force:
        ! A = 2 K_a (theta - theta_eq) / sin(theta)
        ! F_ij = A/l_ij (1/l_kj - cos(theta)/l_ij) r_ij
        ! F_ik = -A/(l_ij l_kj) r_ik
        ! F_jk = A/l_kj (1/l_ij - cos(theta)/l_kj) (-r_kj)

        ! initialization
        nangle   = size(force_consts)
        energy   = 0.0d0
        forces   = 0.0d0
        tbforces = 0.0d0
        displacement = 0.0d0

        do iang=1, nangle
            itbf_ij = iang_to_itbf(iang, 1) ! => ij
            itbf_ik = iang_to_itbf(iang, 2) ! => ik
            itbf_jk = iang_to_itbf(iang, 3) ! => jk

            iatm     = three_atoms(iang, 1)
            jatm     = three_atoms(iang, 2)
            katm     = three_atoms(iang, 3)

            ka       = force_consts(iang)
            theta_eq = theta_eqs(iang)

            r_ij = crd(iatm, :) - crd(jatm, :)
            r_kj = crd(katm, :) - crd(jatm, :)
            r_ik = crd(iatm, :) - crd(katm, :)
            l_ij = sqrt( dot_product(r_ij, r_ij) )
            l_kj = sqrt( dot_product(r_kj, r_kj) )

            ! calculate theta
            cos_theta = dot_product(r_ij, r_kj) / (l_ij*l_kj)
            if (cos_theta >  1.0d0) cos_theta =  1.0d0
            if (cos_theta < -1.0d0) cos_theta = -1.0d0
            theta = acos(cos_theta)

            ! calculate energy
            ene = ka * (theta - theta_eq) ** 2
            energy = energy + ene

            ! calculate force
            sin_theta = sin(theta)
            if (sin_theta < eps) sin_theta = eps ! for sin(theta) == 0
            coeff = 2 * ka * (theta-theta_eq)/sin_theta

            f_i = coeff/l_ij * (r_kj/l_kj - cos_theta*r_ij/l_ij)
            f_k = coeff/l_kj * (r_ij/l_ij - cos_theta*r_kj/l_kj)

            ! store force
            forces(iatm, :) = forces(iatm, :) + f_i(:)
            forces(jatm, :) = forces(jatm, :) - f_i(:) - f_k(:)
            forces(katm, :) = forces(katm, :) + f_k(:)

            ! calculate two-body-force
            f_ij =   coeff/l_ij * (1.0d0/l_kj - cos_theta/l_ij) * r_ij
            f_ik = - coeff/(l_ij*l_kj) * r_ik
            f_jk =   coeff/l_kj * (1.0d0/l_ij - cos_theta/l_kj) * (-r_kj)

            ! store two-body force and two-body distance vector
            if (itbf_ij > 0) then
                tbforces(itbf_ij,:) = tbforces(itbf_ij,:) + f_ij(:)
                displacement(itbf_ij,:) = r_ij(:)
            else
                tbforces(-itbf_ij,:) = tbforces(-itbf_ij,:) - f_ij(:)
                displacement(-itbf_ij,:) = -r_ij(:)
            end if

            if (itbf_ik > 0) then
                tbforces(itbf_ik,:) = tbforces(itbf_ik,:) + f_ik(:)
                displacement(itbf_ik,:) = r_ik(:)
            else
                tbforces(-itbf_ik,:) = tbforces(-itbf_ik,:) - f_ik(:)
                displacement(-itbf_ik,:) = -r_ik(:)
            end if

            if (itbf_jk > 0) then
                tbforces(itbf_jk,:) = tbforces(itbf_jk,:) + f_jk(:)
                displacement(itbf_jk,:) = -r_kj(:)
            else
                tbforces(-itbf_jk,:) = tbforces(-itbf_jk,:) - f_jk(:)
                displacement(-itbf_jk,:) = r_kj(:)
            end if

            if (check) then
                print*, 'TB_CHECK:','i, j, k =', iatm, jatm, katm
                print*, 'TB_CHECK: f_i vs. f_ij + f_ik =>', &
                    & check_tbforce(f_i, f_ij+f_ik)
                print*, 'TB_CHECK:',f_i
                print*, 'TB_CHECK:',f_ij + f_ik
                print*, 'TB_CHECK: f_j vs. f_ji + f_jk =>', & 
                    & check_tbforce(-f_i-f_k, -f_ij+f_jk)
                print*, 'TB_CHECK:',-f_i - f_k
                print*, 'TB_CHECK:',-f_ij + f_jk
                print*, 'TB_CHECK: f_k vs. f_ki + f_kj =>', &
                    & check_tbforce(f_k, -f_ik-f_jk)
                print*, 'TB_CHECK:',f_k
                print*, 'TB_CHECK:',-f_ik - f_jk
            end if

        enddo

    end subroutine

    subroutine print_tbforce()
        use total, only: bonded_pairs, nbonded
        use common_vars
        implicit none
        real(8) :: tbf_abs

        do itbf=1, nbonded
            iatm = bonded_pairs(itbf, 1)
            jatm = bonded_pairs(itbf, 2)

            tbf_abs = dot_product(tbforces(itbf,:), tbforces(itbf,:))
            if (-0.0001 < tbf_abs .and. tbf_abs < 0.0001) cycle
            write(*, '(2I5,3F15.7)'), iatm, jatm, tbforces(itbf_ij, :)
         end do
            
    end subroutine

    subroutine print_atom_order()
        use total, only: bonded_pairs
        use common_vars
        implicit none
        integer :: iatoms(2), jatoms(2), katoms(2)

        print*, '** print atom order for angle **'
        nangle   = size(force_consts)
        do iang=1, nangle
            iatm = three_atoms(iang,1)
            jatm = three_atoms(iang,2)
            katm = three_atoms(iang,3)

            itbf_ij = iang_to_itbf(iang,1)
            if (itbf_ij > 0) then
                iatoms(1) = bonded_pairs(itbf_ij, 1)
                jatoms(1) = bonded_pairs(itbf_ij, 2)
            else
                iatoms(1) = bonded_pairs(-itbf_ij, 2)
                jatoms(1) = bonded_pairs(-itbf_ij, 1)
            end if

            itbf_ik = iang_to_itbf(iang, 2)
            if (itbf_ik > 0) then
                iatoms(2) = bonded_pairs(itbf_ik, 1)
                katoms(1) = bonded_pairs(itbf_ik, 2)
            else
                iatoms(2) = bonded_pairs(-itbf_ik, 2)
                katoms(1) = bonded_pairs(-itbf_ik, 1)
            end if

            itbf_jk = iang_to_itbf(iang, 3)
            if (itbf_jk > 0) then
                jatoms(2) = bonded_pairs(itbf_jk, 1)
                katoms(2) = bonded_pairs(itbf_jk, 2)
            else
                jatoms(2) = bonded_pairs(-itbf_jk, 2)
                katoms(2) = bonded_pairs(-itbf_jk, 1)
            end if

            print*, iang
            print*, iatm, iatoms
            print*, jatm, jatoms
            print*, katm, katoms

         end do
            
    end subroutine

end module

module torsion
    implicit none
    ! input
    integer, allocatable :: four_atoms(:,:)   ! (ntorsion, 4)
    integer, allocatable :: num_torsions(:)   ! (ntorsion)
    integer, allocatable :: num_freqs(:)      ! (ntorsion)
    real(8), allocatable :: force_consts(:)   ! (ntorsion)
    real(8), allocatable :: initial_phases(:) ! (ntorsion)
    integer, allocatable :: itor_to_itbf(:,:) ! (ntorsion, 6)
    ! output
    real(8) :: energy
    real(8), allocatable :: forces(:, :)   ! (natom, 3)
    real(8), allocatable :: tbforces(:, :) ! (ntfb, 3)
    real(8), allocatable :: displacement(:, :) ! (ntbf, 3)

    ! main variales
    integer :: itor, ntorsion
    integer :: iatm, jatm, katm, latm
    integer :: itbf_ij, itbf_ik, itbf_il, itbf_jk, itbf_jl, itbf_kl

contains
    subroutine calculate()
        use total
        use common_vars
        implicit none
        real(8) :: kt, coeff, cos_phi, sin_phi, phi, new_phi, judge
        real(8) :: f_a, f_b(3), f_c(3), n_1(3), n_2(3), l_1, l_2
        integer :: nfreq
        real(8), parameter :: eps = 1.0d-10

        ! Energy:
        ! E_torsion = K_t [1 + cos(eta phi - gamma)]

        ! Force:
        ! f_a = - K_t*ita/2 sin(ita*phi - gamma) / sin(phi)
        ! f_b = (n_2/l_2 - cos(phi)*n_1/l_1) / l_1
        ! f_c = (n_1/l_1 - cos(phi)*n_2/l_2) / l_2
        !
        ! when sin(phi) is nearly equal to 0,
        ! f_a = 
        ! if eta == 1 then  cos(gamma)
        ! if eta == 2 then  2cos(phi)*cos(gamma)
        ! if eta == 3 then  (-4sin(phi)^2 + 3) cos(gamma)
        ! if eta == 4 then  (4(cos(phi)(2cos(phi)^2 -1))) cos(gamma)
        ! if eta == 5 then  (1 + 4*cos(phi)^2 (1-4*sin(phi)^2))) cos(gamma)
        ! if eta == 6 then 

        ! F_i = f_a (f_b x r_jk)
        ! F_j = f_a (f_c x r_kl - f_b x r_ik)
        ! F_k = f_a (f_b x r_ij - f_c x r_jl)
        ! F_l = f_a (f_c x r_jk)

        ! initialization
        ntorsion = size(force_consts)
        energy   = 0.0d0
        forces   = 0.0d0
        tbforces = 0.0d0
        displacement = 0.0d0

        ! print*, num_torsions

        do itor=1, ntorsion

            ! If force constant of torsional energy is small value that 
            ! less than eps, then do not calculate it
            if( dabs( force_consts(itor) ) <= eps ) cycle

            itbf_ij = itor_to_itbf(itor, 1) ! => ij
            itbf_ik = itor_to_itbf(itor, 2) ! => ik
            itbf_il = itor_to_itbf(itor, 3) ! => il
            itbf_jk = itor_to_itbf(itor, 4) ! => jk
            itbf_jl = itor_to_itbf(itor, 5) ! => jl
            itbf_kl = itor_to_itbf(itor, 6) ! => kl

            iatm  = four_atoms(itor, 1)
            jatm  = four_atoms(itor, 2)
            katm  = four_atoms(itor, 3)
            latm  = four_atoms(itor, 4)
            kt    = force_consts(itor)
            nfreq = num_freqs(itor)

            r_ij = crd(iatm,:) - crd(jatm,:)
            r_ik = crd(iatm,:) - crd(katm,:)
            r_il = crd(iatm,:) - crd(latm,:)
            r_jk = crd(jatm,:) - crd(katm,:)
            r_jl = crd(jatm,:) - crd(latm,:)
            r_kl = crd(katm,:) - crd(latm,:)
            l_ij = sqrt( dot_product(r_ij, r_ij) )
            l_ik = sqrt( dot_product(r_ik, r_ik) )
            l_il = sqrt( dot_product(r_il, r_il) )
            l_jk = sqrt( dot_product(r_jk, r_jk) )
            l_jl = sqrt( dot_product(r_jl, r_jl) )
            l_kl = sqrt( dot_product(r_kl, r_kl) )

            ! make normal vector
            n_1 = outer_prod(r_ij, -r_jk)
            n_2 = outer_prod(-r_jk, r_kl)
            l_1 = sqrt( dot_product(n_1, n_1) )
            l_2 = sqrt( dot_product(n_2, n_2) )
            ! l_1 = sqrt( dot_product(r_kl, r_kl) )
            ! l_2 = sqrt( dot_product(r_kl, r_kl) )

            cos_phi = dot_product(n_1, n_2)/(l_1*l_2)
            if (cos_phi >  1.0d0) cos_phi =  1.0d0
            if (cos_phi < -1.0d0) cos_phi = -1.0d0

            ! judge sign and calculate phi
            phi = acos(cos_phi)
            judge = -dot_product(-r_jk, outer_prod(n_1, n_2))
            if (judge < 0.0d0) phi = 2.0d0*PI - phi

            ! calculate energy
            new_phi = nfreq*phi + initial_phases(itor)
            ene = kt * (1.0d0 + cos(new_phi)) / num_torsions(itor)
            energy = energy + ene

            ! calculate force
            ! when phi = 0.0d0
            sin_phi = sin(phi)
            if (nfreq == 1) then
                coeff = 1.0d0
            else if (nfreq == 2) then
                coeff = 2.0d0 * cos_phi
            else if (nfreq == 3) then
                coeff = 3.0d0 - 4.0d0 * sin_phi**2
            else if (nfreq == 4) then
                coeff = (4.0d0 - 8.0d0*sin_phi**2) * cos_phi
            else if (nfreq == 5) then
                coeff = (16.0d0*sin_phi**2 - 20.0d0) * (sin_phi**2) + 5.0d0
            else if (nfreq == 6) then
                coeff = 32.0d0*(cos_phi**2)*(cos_phi**2 - 1.0d0)
                coeff = cos_phi * (6.0d0 + coeff)
            end if
            coeff = coeff * cos(initial_phases(itor))

            f_a = - kt*nfreq*coeff/num_torsions(itor)

            f_b = (n_2/l_2 - cos(phi)*n_1/l_1) / l_1
            f_c = (n_1/l_1 - cos(phi)*n_2/l_2) / l_2

            f_i = f_a * outer_prod(f_b, r_jk)
            f_j = f_a * ( outer_prod(f_c, r_kl) - outer_prod(f_b, r_ik) )
            f_k = f_a * ( outer_prod(f_b, r_ij) - outer_prod(f_c, r_jl) )
            f_l = f_a * outer_prod(f_c, r_jk)

            ! store force
            forces(iatm, :) = forces(iatm, :) + f_i(:)
            forces(jatm, :) = forces(jatm, :) + f_j(:)
            forces(katm, :) = forces(katm, :) + f_k(:)
            forces(latm, :) = forces(latm, :) + f_l(:)

            ! calculate two-body-force
            f_ij = f_a/l_1 * ( dot_product(-r_jk, r_kl)/l_2 &
            &                + dot_product(-r_jk, r_ik)*cos_phi/l_1 ) * r_ij

            f_ik = f_a/l_1 * ( dot_product(-r_jk, -r_jl)/l_2 &
            &                - dot_product(r_ij, -r_jk)*cos_phi/l_1 ) * r_ik

            f_il = - f_a * dot_product(r_jk, r_jk) / (l_1*l_2) * r_il

            f_jk = - f_a * ( dot_product(r_kl, r_ij)/(l_1*l_2) &
            &              + dot_product(r_ik, r_jl)/(l_1*l_2) &
            &              + dot_product(r_ij, r_ik)*cos_phi/(l_1**2) &
            &              + dot_product(r_kl, r_jl)*cos_phi/(l_2**2) ) * r_jk

            f_jl = f_a/l_2 * ( dot_product(r_ik, r_jk)/l_1 &
            &                - dot_product(-r_jk, r_kl)*cos_phi/l_2 ) * r_jl

            f_kl = f_a/l_2 * ( dot_product(r_ij, -r_jk)/l_1 &
            &                - dot_product(-r_jk, -r_jl)*cos_phi/l_2) * r_kl

            ! store two-body force and two-body distance vector
            if (itbf_ij > 0) then
                tbforces(itbf_ij,:) = tbforces(itbf_ij,:) + f_ij(:)
                displacement(itbf_ij,:) = r_ij(:)
            else
                tbforces(-itbf_ij,:) = tbforces(-itbf_ij,:) - f_ij(:)
                displacement(-itbf_ij,:) = -r_ij(:)
            end if

            if (itbf_ik > 0) then
                tbforces(itbf_ik,:) = tbforces(itbf_ik,:) + f_ik(:)
                displacement(itbf_ik,:) = r_ik(:)
            else
                tbforces(-itbf_ik,:) = tbforces(-itbf_ik,:) - f_ik(:)
                displacement(-itbf_ik,:) = -r_ik(:)
            end if

            if (itbf_il > 0) then
                tbforces(itbf_il,:) = tbforces(itbf_il,:) + f_il(:)
                displacement(itbf_il,:) = r_il(:)
            else
                tbforces(-itbf_il,:) = tbforces(-itbf_il,:) - f_il(:)
                displacement(-itbf_il,:) = -r_il(:)
            end if

            if (itbf_jk > 0) then
                tbforces(itbf_jk,:) = tbforces(itbf_jk,:) + f_jk(:)
                displacement(itbf_jk,:) = r_jk(:)
            else
                tbforces(-itbf_jk,:) = tbforces(-itbf_jk,:) - f_jk(:)
                displacement(-itbf_jk,:) = -r_jk(:)
            end if

            if (itbf_jl > 0) then
                tbforces(itbf_jl,:) = tbforces(itbf_jl,:) + f_jl(:)
                displacement(itbf_jl,:) = r_jl(:)
            else
                tbforces(-itbf_jl,:) = tbforces(-itbf_jl,:) - f_jl(:)
                displacement(-itbf_jl,:) = -r_jl(:)
            end if

            if (itbf_kl > 0) then
                tbforces(itbf_kl,:) = tbforces(itbf_kl,:) + f_kl(:)
                displacement(itbf_kl,:) = r_kl(:)
            else
                tbforces(-itbf_kl,:) = tbforces(-itbf_kl,:) - f_kl(:)
                displacement(-itbf_kl,:) = -r_kl(:)
            end if

            if (check) then
                print*, 'TB_CHECK: i, j, k, l=', iatm, jatm, katm, latm
                print*, 'TB_CHECK: f_i vs. f_ij + f_ik + f_il =>', &
                    & check_tbforce(f_i, f_ij+f_ik+f_il)
                print*, 'TB_CHECK:',f_i
                print*, 'TB_CHECK:',f_ij + f_ik + f_il
                print*, 'TB_CHECK: f_j vs. f_ji + f_jk + f_jl =>', & 
                    & check_tbforce(f_j, -f_ij+f_jk+f_jl)
                print*, 'TB_CHECK:',f_j
                print*, 'TB_CHECK:',- f_ij + f_jk + f_jl
                print*, 'TB_CHECK: f_k vs. f_ki + f_kj + f_kl =>', &
                    & check_tbforce(f_k, -f_ik-f_jk+f_kl)
                print*, 'TB_CHECK:',f_k
                print*, 'TB_CHECK:',- f_ik - f_jk + f_kl
                print*, 'TB_CHECK: f_l vs. f_li + f_lj + f_lk =>', &
                    & check_tbforce(f_l, -f_il-f_jl-f_kl)
                print*, 'TB_CHECK:',f_l
                print*, 'TB_CHECK:',- f_il - f_jl - f_kl
            end if

        enddo

    end subroutine

    subroutine print_tbforce()
        use total, only: bonded_pairs, nbonded
        use common_vars
        implicit none
        real(8) :: tbf_abs

        do itbf=1, nbonded
            iatm = bonded_pairs(itbf, 1)
            jatm = bonded_pairs(itbf, 2)

            tbf_abs = dot_product(tbforces(itbf,:), tbforces(itbf,:))
            if (-0.0001 < tbf_abs .and. tbf_abs < 0.0001) cycle
            write(*, '(2I5,3F15.7)'), iatm, jatm, tbforces(itbf_ij, :)
         end do
            
    end subroutine

    subroutine print_atom_order()
        use total, only: bonded_pairs
        use common_vars
        implicit none
        integer :: iatoms(3), jatoms(3), katoms(3), latoms(3)

        print*, '** print atom order for torsion **'
        ntorsion = size(force_consts)
        do itor=1, ntorsion
            iatm = four_atoms(itor,1)
            jatm = four_atoms(itor,2)
            katm = four_atoms(itor,3)
            latm = four_atoms(itor,4)

            itbf_ij = itor_to_itbf(itor, 1)
            if (itbf_ij > 0) then
                iatoms(1) = bonded_pairs(itbf_ij, 1)
                jatoms(1) = bonded_pairs(itbf_ij, 2)
            else
                iatoms(1) = bonded_pairs(-itbf_ij, 2)
                jatoms(1) = bonded_pairs(-itbf_ij, 1)
            end if

            itbf_ik = itor_to_itbf(itor, 2)
            if (itbf_ik > 0) then
                iatoms(2) = bonded_pairs(itbf_ik, 1)
                katoms(1) = bonded_pairs(itbf_ik, 2)
            else
                iatoms(2) = bonded_pairs(-itbf_ik, 2)
                katoms(1) = bonded_pairs(-itbf_ik, 1)
            end if

            itbf_il = itor_to_itbf(itor, 3)
            if (itbf_il > 0) then
                iatoms(3) = bonded_pairs(itbf_il, 1)
                latoms(1) = bonded_pairs(itbf_il, 2)
            else
                iatoms(3) = bonded_pairs(-itbf_il, 2)
                latoms(1) = bonded_pairs(-itbf_il, 1)
            end if

            itbf_jk = itor_to_itbf(itor, 4)
            if (itbf_jk > 0) then
                jatoms(2) = bonded_pairs(itbf_jk, 1)
                katoms(2) = bonded_pairs(itbf_jk, 2)
            else
                jatoms(2) = bonded_pairs(-itbf_jk, 2)
                katoms(2) = bonded_pairs(-itbf_jk, 1)
            end if

            itbf_jl = itor_to_itbf(itor, 5)
            if (itbf_jl > 0) then
                jatoms(3) = bonded_pairs(itbf_jl, 1)
                latoms(2) = bonded_pairs(itbf_jl, 2)
            else
                jatoms(3) = bonded_pairs(-itbf_jl, 2)
                latoms(2) = bonded_pairs(-itbf_jl, 1)
            end if

            itbf_kl = itor_to_itbf(itor, 6)
            if (itbf_kl > 0) then
                katoms(3) = bonded_pairs(itbf_kl, 1)
                latoms(3) = bonded_pairs(itbf_kl, 2)
            else
                katoms(3) = bonded_pairs(-itbf_kl, 2)
                latoms(3) = bonded_pairs(-itbf_kl, 1)
            end if

            print*, itor
            print*, iatm, iatoms
            print*, jatm, jatoms
            print*, katm, katoms
            print*, latm, latoms

         end do
            
    end subroutine

end module

module improper
    implicit none
    ! input
    integer, allocatable :: four_atoms(:,:)   ! (ntorsion, 4)
    integer, allocatable :: num_torsions(:)   ! (ntorsion)
    integer, allocatable :: num_freqs(:)      ! (ntorsion)
    real(8), allocatable :: force_consts(:)   ! (ntorsion)
    real(8), allocatable :: initial_phases(:) ! (ntorsion)
    integer, allocatable :: itor_to_itbf(:,:) ! (ntorsion, 6)
    ! output
    real(8) :: energy
    real(8), allocatable :: forces(:, :)   ! (natom, 3)
    real(8), allocatable :: tbforces(:, :) ! (ntfb, 3)
    real(8), allocatable :: displacement(:, :)! (ntbf, 3)

    ! main variales
    integer :: itor, ntorsion
    integer :: iatm, jatm, katm, latm
    integer :: itbf_ij, itbf_ik, itbf_il, itbf_jk, itbf_jl, itbf_kl

contains
    subroutine calculate()
        use total
        use common_vars
        implicit none
        real(8) :: kt, coeff, cos_phi, sin_phi, phi, new_phi, judge
        real(8) :: f_a, f_b(3), f_c(3), n_1(3), n_2(3), l_1, l_2
        integer :: nfreq
        real(8), parameter :: eps = 1.0d-10

        ! Energy:
        ! E_torsion = K_t [1 + cos(eta phi - gamma)]

        ! Force:
        ! f_a = - K_t*ita/2 sin(ita*phi - gamma) / sin(phi)
        ! f_b = (n_2/l_2 - cos(phi)*n_1/l_1) / l_1
        ! f_c = (n_1/l_1 - cos(phi)*n_2/l_2) / l_2
        !
        ! when sin(phi) is nearly equal to 0,
        ! f_a = 
        ! if eta == 1 then  cos(gamma)
        ! if eta == 2 then  2cos(phi)*cos(gamma)
        ! if eta == 3 then  (-4sin(phi)^2 + 3) cos(gamma)
        ! if eta == 4 then  (4(cos(phi)(2cos(phi)^2 -1))) cos(gamma)
        ! if eta == 5 then  (1 + 4*cos(phi)^2 (1-4*sin(phi)^2))) cos(gamma)
        ! if eta == 6 then 

        ! F_i = f_a (f_b x r_jk)
        ! F_j = f_a (f_c x r_kl - f_b x r_ik)
        ! F_k = f_a (f_b x r_ij - f_c x r_jl)
        ! F_l = f_a (f_c x r_jk)

        ! initialization
        ntorsion = size(force_consts)
        energy   = 0.0d0
        forces   = 0.0d0
        tbforces = 0.0d0
        displacement = 0.0d0

        ! print*, num_torsions

        do itor=1, ntorsion

            ! If force constant of torsional energy is small value that 
            ! less than eps, then do not calculate it
            if( dabs( force_consts(itor) ) <= eps ) cycle

            itbf_ij = itor_to_itbf(itor, 1) ! => ij
            itbf_ik = itor_to_itbf(itor, 2) ! => ik
            itbf_il = itor_to_itbf(itor, 3) ! => il
            itbf_jk = itor_to_itbf(itor, 4) ! => jk
            itbf_jl = itor_to_itbf(itor, 5) ! => jl
            itbf_kl = itor_to_itbf(itor, 6) ! => kl

            iatm  = four_atoms(itor, 1)
            jatm  = four_atoms(itor, 2)
            katm  = four_atoms(itor, 3)
            latm  = four_atoms(itor, 4)
            kt    = force_consts(itor)
            nfreq = num_freqs(itor)

            r_ij = crd(iatm,:) - crd(jatm,:)
            r_ik = crd(iatm,:) - crd(katm,:)
            r_il = crd(iatm,:) - crd(latm,:)
            r_jk = crd(jatm,:) - crd(katm,:)
            r_jl = crd(jatm,:) - crd(latm,:)
            r_kl = crd(katm,:) - crd(latm,:)
            l_ij = sqrt( dot_product(r_ij, r_ij) )
            l_ik = sqrt( dot_product(r_ik, r_ik) )
            l_il = sqrt( dot_product(r_il, r_il) )
            l_jk = sqrt( dot_product(r_jk, r_jk) )
            l_jl = sqrt( dot_product(r_jl, r_jl) )
            l_kl = sqrt( dot_product(r_kl, r_kl) )

            ! make normal vector
            n_1 = outer_prod(r_ij, -r_jk)
            n_2 = outer_prod(-r_jk, r_kl)
            l_1 = sqrt( dot_product(n_1, n_1) )
            l_2 = sqrt( dot_product(n_2, n_2) )
            ! l_1 = sqrt( dot_product(r_kl, r_kl) )
            ! l_2 = sqrt( dot_product(r_kl, r_kl) )

            cos_phi = dot_product(n_1, n_2)/(l_1*l_2)
            if (cos_phi >  1.0d0) cos_phi =  1.0d0
            if (cos_phi < -1.0d0) cos_phi = -1.0d0

            ! judge sign and calculate phi
            phi = acos(cos_phi)
            judge = -dot_product(-r_jk, outer_prod(n_1, n_2))
            if (judge < 0.0d0) phi = 2.0d0*PI - phi

            ! calculate energy
            new_phi = nfreq*phi + initial_phases(itor)
            ene = kt * (1.0d0 + cos(new_phi)) / num_torsions(itor)
            energy = energy + ene

            ! calculate force
            ! when phi = 0.0d0
            sin_phi = sin(phi)
            if (nfreq == 1) then
                coeff = 1.0d0
            else if (nfreq == 2) then
                coeff = 2.0d0 * cos_phi
            else if (nfreq == 3) then
                coeff = 3.0d0 - 4.0d0 * sin_phi**2
            else if (nfreq == 4) then
                coeff = (4.0d0 - 8.0d0*sin_phi**2) * cos_phi
            else if (nfreq == 5) then
                coeff = (16.0d0*sin_phi**2 - 20.0d0) * (sin_phi**2) + 5.0d0
            else if (nfreq == 6) then
                coeff = 32.0d0*(cos_phi**2)*(cos_phi**2 - 1.0d0)
                coeff = cos_phi * (6.0d0 + coeff)
            end if
            coeff = coeff * cos(initial_phases(itor))

            f_a = - kt*nfreq*coeff/num_torsions(itor)

            f_b = (n_2/l_2 - cos(phi)*n_1/l_1) / l_1
            f_c = (n_1/l_1 - cos(phi)*n_2/l_2) / l_2

            f_i = f_a * outer_prod(f_b, r_jk)
            f_j = f_a * ( outer_prod(f_c, r_kl) - outer_prod(f_b, r_ik) )
            f_k = f_a * ( outer_prod(f_b, r_ij) - outer_prod(f_c, r_jl) )
            f_l = f_a * outer_prod(f_c, r_jk)

            ! store force
            forces(iatm, :) = forces(iatm, :) + f_i(:)
            forces(jatm, :) = forces(jatm, :) + f_j(:)
            forces(katm, :) = forces(katm, :) + f_k(:)
            forces(latm, :) = forces(latm, :) + f_l(:)

            ! calculate two-body-force
            f_ij = f_a/l_1 * ( dot_product(-r_jk, r_kl)/l_2 &
            &                + dot_product(-r_jk, r_ik)*cos_phi/l_1 ) * r_ij

            f_ik = f_a/l_1 * ( dot_product(-r_jk, -r_jl)/l_2 &
            &                - dot_product(r_ij, -r_jk)*cos_phi/l_1 ) * r_ik

            f_il = - f_a * dot_product(r_jk, r_jk) / (l_1*l_2) * r_il

            f_jk = - f_a * ( dot_product(r_kl, r_ij)/(l_1*l_2) &
            &              + dot_product(r_ik, r_jl)/(l_1*l_2) &
            &              + dot_product(r_ij, r_ik)*cos_phi/(l_1**2) &
            &              + dot_product(r_kl, r_jl)*cos_phi/(l_2**2) ) * r_jk

            f_jl = f_a/l_2 * ( dot_product(r_ik, r_jk)/l_1 &
            &                - dot_product(-r_jk, r_kl)*cos_phi/l_2 ) * r_jl

            f_kl = f_a/l_2 * ( dot_product(r_ij, -r_jk)/l_1 &
            &                - dot_product(-r_jk, -r_jl)*cos_phi/l_2) * r_kl

            ! store two-body force and two-body distance vector
            if (itbf_ij > 0) then
                tbforces(itbf_ij,:) = tbforces(itbf_ij,:) + f_ij(:)
                displacement(itbf_ij,:) = r_ij(:)
            else
                tbforces(-itbf_ij,:) = tbforces(-itbf_ij,:) - f_ij(:)
                displacement(-itbf_ij,:) = -r_ij(:)
            end if

            if (itbf_ik > 0) then
                tbforces(itbf_ik,:) = tbforces(itbf_ik,:) + f_ik(:)
                displacement(itbf_ik,:) = r_ik(:)
            else
                tbforces(-itbf_ik,:) = tbforces(-itbf_ik,:) - f_ik(:)
                displacement(-itbf_ik,:) = -r_ik(:)
            end if

            if (itbf_il > 0) then
                tbforces(itbf_il,:) = tbforces(itbf_il,:) + f_il(:)
                displacement(itbf_il,:) = r_il(:)
            else
                tbforces(-itbf_il,:) = tbforces(-itbf_il,:) - f_il(:)
                displacement(-itbf_il,:) = -r_il(:)
            end if

            if (itbf_jk > 0) then
                tbforces(itbf_jk,:) = tbforces(itbf_jk,:) + f_jk(:)
                displacement(itbf_jk,:) = r_jk(:)
            else
                tbforces(-itbf_jk,:) = tbforces(-itbf_jk,:) - f_jk(:)
                displacement(-itbf_jk,:) = -r_jk(:)
            end if

            if (itbf_jl > 0) then
                tbforces(itbf_jl,:) = tbforces(itbf_jl,:) + f_jl(:)
                displacement(itbf_jl,:) = r_jl(:)
            else
                tbforces(-itbf_jl,:) = tbforces(-itbf_jl,:) - f_jl(:)
                displacement(-itbf_jl,:) = -r_jl(:)
            end if

            if (itbf_kl > 0) then
                tbforces(itbf_kl,:) = tbforces(itbf_kl,:) + f_kl(:)
                displacement(itbf_kl,:) = r_kl(:)
            else
                tbforces(-itbf_kl,:) = tbforces(-itbf_kl,:) - f_kl(:)
                displacement(-itbf_kl,:) = -r_kl(:)
            end if

            if (check) then
                print*, 'TB_CHECK: i, j, k, l=', iatm, jatm, katm, latm
                print*, 'TB_CHECK: f_i vs. f_ij + f_ik + f_il =>', &
                    & check_tbforce(f_i, f_ij+f_ik+f_il)
                print*, 'TB_CHECK:',f_i
                print*, 'TB_CHECK:',f_ij + f_ik + f_il
                print*, 'TB_CHECK: f_j vs. f_ji + f_jk + f_jl =>', & 
                    & check_tbforce(f_j, -f_ij+f_jk+f_jl)
                print*, 'TB_CHECK:',f_j
                print*, 'TB_CHECK:',- f_ij + f_jk + f_jl
                print*, 'TB_CHECK: f_k vs. f_ki + f_kj + f_kl =>', &
                    & check_tbforce(f_k, -f_ik-f_jk+f_kl)
                print*, 'TB_CHECK:',f_k
                print*, 'TB_CHECK:',- f_ik - f_jk + f_kl
                print*, 'TB_CHECK: f_l vs. f_li + f_lj + f_lk =>', &
                    & check_tbforce(f_l, -f_il-f_jl-f_kl)
                print*, 'TB_CHECK:',f_l
                print*, 'TB_CHECK:',- f_il - f_jl - f_kl
            end if

        enddo

    end subroutine

    subroutine print_tbforce()
        use total, only: bonded_pairs, nbonded
        use common_vars
        implicit none
        real(8) :: tbf_abs

        do itbf=1, nbonded
            iatm = bonded_pairs(itbf, 1)
            jatm = bonded_pairs(itbf, 2)

            tbf_abs = dot_product(tbforces(itbf,:), tbforces(itbf,:))
            if (-0.0001 < tbf_abs .and. tbf_abs < 0.0001) cycle
            write(*, '(2I5,3F15.7)'), iatm, jatm, tbforces(itbf_ij, :)
         end do
            
    end subroutine

    subroutine print_atom_order()
        use total, only: bonded_pairs
        use common_vars
        implicit none
        integer :: iatoms(3), jatoms(3), katoms(3), latoms(3)

        print*, '** print atom order for torsion **'
        ntorsion = size(force_consts)
        do itor=1, ntorsion
            iatm = four_atoms(itor,1)
            jatm = four_atoms(itor,2)
            katm = four_atoms(itor,3)
            latm = four_atoms(itor,4)

            itbf_ij = itor_to_itbf(itor, 1)
            if (itbf_ij > 0) then
                iatoms(1) = bonded_pairs(itbf_ij, 1)
                jatoms(1) = bonded_pairs(itbf_ij, 2)
            else
                iatoms(1) = bonded_pairs(-itbf_ij, 2)
                jatoms(1) = bonded_pairs(-itbf_ij, 1)
            end if

            itbf_ik = itor_to_itbf(itor, 2)
            if (itbf_ik > 0) then
                iatoms(2) = bonded_pairs(itbf_ik, 1)
                katoms(1) = bonded_pairs(itbf_ik, 2)
            else
                iatoms(2) = bonded_pairs(-itbf_ik, 2)
                katoms(1) = bonded_pairs(-itbf_ik, 1)
            end if

            itbf_il = itor_to_itbf(itor, 3)
            if (itbf_il > 0) then
                iatoms(3) = bonded_pairs(itbf_il, 1)
                latoms(1) = bonded_pairs(itbf_il, 2)
            else
                iatoms(3) = bonded_pairs(-itbf_il, 2)
                latoms(1) = bonded_pairs(-itbf_il, 1)
            end if

            itbf_jk = itor_to_itbf(itor, 4)
            if (itbf_jk > 0) then
                jatoms(2) = bonded_pairs(itbf_jk, 1)
                katoms(2) = bonded_pairs(itbf_jk, 2)
            else
                jatoms(2) = bonded_pairs(-itbf_jk, 2)
                katoms(2) = bonded_pairs(-itbf_jk, 1)
            end if

            itbf_jl = itor_to_itbf(itor, 5)
            if (itbf_jl > 0) then
                jatoms(3) = bonded_pairs(itbf_jl, 1)
                latoms(2) = bonded_pairs(itbf_jl, 2)
            else
                jatoms(3) = bonded_pairs(-itbf_jl, 2)
                latoms(2) = bonded_pairs(-itbf_jl, 1)
            end if

            itbf_kl = itor_to_itbf(itor, 6)
            if (itbf_kl > 0) then
                katoms(3) = bonded_pairs(itbf_kl, 1)
                latoms(3) = bonded_pairs(itbf_kl, 2)
            else
                katoms(3) = bonded_pairs(-itbf_kl, 2)
                latoms(3) = bonded_pairs(-itbf_kl, 1)
            end if

            print*, itor
            print*, iatm, iatoms
            print*, jatm, jatoms
            print*, katm, katoms
            print*, latm, latoms

         end do
            
    end subroutine

end module

!###############################################################################
module cmap
    implicit none
    ! input
    integer, allocatable :: five_atoms(:, :)            ! (ncmap, 5)
    integer, allocatable :: cmap_types(:)               ! (ncmap)
    integer, allocatable :: cmap_grid_step_size(:)      ! (ntypes)
    integer, allocatable :: cmap_resolutions(:)          ! (ntypes)
    real(8), allocatable :: cmap_grid_energy(:, :, :)   ! (ntypes, ngrid_x(24), ngrid_y(24))
    real(8), allocatable :: cmap_grid_dphi(:, :, :)     ! (ntypes, ngrid_x(24), ngrid_y(24))
    real(8), allocatable :: cmap_grid_dpsi(:, :, :)     ! (ntypes, ngrid_x(24), ngrid_y(24))
    real(8), allocatable :: cmap_grid_dphi_dpsi(:, :, :)! (ntypes, ngrid_x(24), ngrid_y(24))
    integer, allocatable :: icmp_to_itbf(:, :)         ! (ncmap, 10)
    ! output
    real(8) :: energy
    real(8), allocatable :: forces(:, :)            ! (natom, 3)
    real(8), allocatable :: tbforces(:, :)          ! (ntfb, 3)
    real(8), allocatable :: displacement(:, :)      ! (ntbf, 3)

    ! main values
    integer :: icmap, ncmap
    integer :: iatm, jatm, katm, latm, matm, cmap_type
    integer :: itbf_ij, itbf_ik, itbf_il, itbf_im
    integer :: itbf_jk, itbf_jl, itbf_jm, itbf_kl, itbf_km, itbf_lm
    real(8), dimension(2,2) :: E_grid, dphi_grid, dpsi_grid, dphidpsi_grid
    real(8), dimension(4)   :: E_grid_1D, dphi_grid_1D, dpsi_grid_1D, dphidpsi_grid_1D

contains
    subroutine calculate()

        use total
        use common_vars
        implicit none
        real(8) :: phi, psi, cos_phi, cos_psi, deg_phi, deg_psi
        real(8) :: sin_phi, sin_psi
        real(8) :: cosphi, sinphi, cospsi, sinpsi, phii, psii
        real(8) :: dphi_dra(4,3), dpsi_dra(4,3)
        real(8) :: n_1(3), n_2(3), n_3(3), l_1, l_2, l_3
        real(8) :: judge_phi, judge_psi
        real(8) :: dx, dy
        real(8) :: E, dE_dphi, dE_dpsi
        real(8) :: f_phi, f_psi, f_a(3), f_b(3), f_c(3), f_d(3)
        real(8) :: d_1, d_2, d_3, d_4, d_phi, d_psi
        real(8), dimension(4,4) :: aij
        real(8), parameter :: RAD_TO_DEG = 180.d0 / PI
        integer :: gridstart = -180
        integer :: x, y, grid_step_size
        integer :: i, j

        ! Energy:
        ! E_cmap = 

        ! initialization
        ncmap = size(cmap_types)
        energy = 0.0d0
        forces = 0.0d0
        tbforces = 0.0d0
        displacement = 0.0d0


        do icmap=1, ncmap

            ! get the index of bonded_pairs of each two atoms
            itbf_ij = icmp_to_itbf(icmap, 1)
            itbf_ik = icmp_to_itbf(icmap, 2)
            itbf_il = icmp_to_itbf(icmap, 3)
            itbf_im = icmp_to_itbf(icmap, 4)
            itbf_jk = icmp_to_itbf(icmap, 5)
            itbf_jl = icmp_to_itbf(icmap, 6)
            itbf_jm = icmp_to_itbf(icmap, 7)
            itbf_kl = icmp_to_itbf(icmap, 8)
            itbf_km = icmp_to_itbf(icmap, 9)
            itbf_lm = icmp_to_itbf(icmap, 10)

            ! get the five atoms for this cmap
            iatm = five_atoms(icmap, 1)
            jatm = five_atoms(icmap, 2)
            katm = five_atoms(icmap, 3)
            latm = five_atoms(icmap, 4)
            matm = five_atoms(icmap, 5)

            ! get cmap type
            cmap_type = cmap_types(icmap)

            ! get the vectors and lengths between atoms
            r_ij = crd(iatm, :) - crd(jatm, :)
            r_ik = crd(iatm, :) - crd(katm, :)
            r_il = crd(iatm, :) - crd(latm, :)
            r_im = crd(iatm, :) - crd(matm, :)
            r_jk = crd(jatm, :) - crd(katm, :)
            r_jl = crd(jatm, :) - crd(latm, :)
            r_jm = crd(jatm, :) - crd(matm, :)
            r_kl = crd(katm, :) - crd(latm, :)
            r_km = crd(katm, :) - crd(matm, :)
            r_lm = crd(latm, :) - crd(matm, :)
            l_ij = sqrt( dot_product(r_ij, r_ij) )
            l_ik = sqrt( dot_product(r_ik, r_ik) )
            l_il = sqrt( dot_product(r_il, r_il) )
            l_im = sqrt( dot_product(r_im, r_im) )
            l_jk = sqrt( dot_product(r_jk, r_jk) )
            l_jl = sqrt( dot_product(r_jl, r_jl) )
            l_jm = sqrt( dot_product(r_jm, r_jm) )
            l_kl = sqrt( dot_product(r_kl, r_kl) )
            l_km = sqrt( dot_product(r_km, r_km) )
            l_lm = sqrt( dot_product(r_lm, r_lm) )

            ! make normal vectors
            n_1 = outer_prod(r_ij, -r_jk)
            n_2 = outer_prod(-r_jk, r_kl)
            n_3 = outer_prod(-r_kl, r_lm)
            l_1 = sqrt( dot_product(n_1, n_1) )
            l_2 = sqrt( dot_product(n_2, n_2) )
            l_3 = sqrt( dot_product(n_3, n_3) )

            call cal_torsion(r_ij, -r_jk, -r_kl, dphi_dra, cosphi, sinphi)
            call cal_torsion(r_jk, -r_kl, -r_lm, dpsi_dra, cospsi, sinpsi)

            phii = sign(acos(cosphi), sinphi) * RAD_TO_DEG
            psii = sign(acos(cospsi), sinpsi) * RAD_TO_DEG

            ! calculate dihedral angle phi
            cos_phi = dot_product(n_1, n_2) / (l_1 * l_2)
            if (cos_phi >  1.0d0) cos_phi =  1.0d0
            if (cos_phi < -1.0d0) cos_phi = -1.0d0

            phi = acos(cos_phi)
            judge_phi = dot_product(-r_jk, outer_prod(n_1, n_2))
            if (judge_phi < 0.0d0) phi = 2.0d0 * PI - phi
            deg_phi = phi * RAD_TO_DEG

            ! calculate dihedral angle psi
            cos_psi = dot_product(n_2, n_3) / (l_2 * l_3)
            if (cos_psi >  1.0d0) cos_psi =  1.0d0
            if (cos_psi < -1.0d0) cos_psi = -1.0d0

            psi = acos(cos_psi)
            judge_psi = dot_product(-r_kl, outer_prod(n_2, n_3))
            if (judge_psi < 0.0d0) psi = 2.0d0 * PI - psi
            deg_psi = psi * RAD_TO_DEG

            ! determine grid indices for phi and psi
            grid_step_size = cmap_grid_step_size(cmap_type)
            x = int( (phii - gridstart) / grid_step_size ) + 1
            y = int( (psii - gridstart) / grid_step_size ) + 1

            ! compute fractional positions inside grid cell
            dx = modulo(phii - dble(gridstart), dble(grid_step_size)) / dble(grid_step_size)

            dy = modulo(psii - dble(gridstart), dble(grid_step_size)) / dble(grid_step_size)

            ! get the 2x2 grid energies and derivatives
            do i=1,2
                do j=1,2
                    E_grid(i,j)        = cmap_grid_energy(cmap_type, calculate_cmap_grid(cmap_type, y+i-1), &
                                                            calculate_cmap_grid (cmap_type, x+j-1) )
                    dphi_grid(i,j)     = cmap_grid_dphi(cmap_type, calculate_cmap_grid(cmap_type, y+i-1), &
                                                            calculate_cmap_grid (cmap_type, x+j-1) )
                    dpsi_grid(i,j)     = cmap_grid_dpsi(cmap_type, calculate_cmap_grid(cmap_type, y+i-1), &
                                                            calculate_cmap_grid (cmap_type, x+j-1) )
                    dphidpsi_grid(i,j) = cmap_grid_dphi_dpsi(cmap_type, calculate_cmap_grid(cmap_type, y+i-1), &
                                                            calculate_cmap_grid (cmap_type, x+j-1) )
                end do
            end do

            ! get 1d energy and derivatives
            E_grid_1D(1) = E_grid(1,1)
            E_grid_1D(2) = E_grid(1,2)
            E_grid_1D(3) = E_grid(2,2)
            E_grid_1D(4) = E_grid(2,1)

            dphi_grid_1D(1) = dphi_grid(1,1)
            dphi_grid_1D(2) = dphi_grid(1,2)
            dphi_grid_1D(3) = dphi_grid(2,2)
            dphi_grid_1D(4) = dphi_grid(2,1)

            dpsi_grid_1D(1) = dpsi_grid(1,1)
            dpsi_grid_1D(2) = dpsi_grid(1,2)
            dpsi_grid_1D(3) = dpsi_grid(2,2)
            dpsi_grid_1D(4) = dpsi_grid(2,1)

            dphidpsi_grid_1D(1) = dphidpsi_grid(1,1)
            dphidpsi_grid_1D(2) = dphidpsi_grid(1,2)
            dphidpsi_grid_1D(3) = dphidpsi_grid(2,2)
            dphidpsi_grid_1D(4) = dphidpsi_grid(2,1)

            ! calculate cmap coefficients
            call cmap_coeff(grid_step_size, &
                            E_grid_1D, &
                            dphi_grid_1D, &
                            dpsi_grid_1D, &
                            dphidpsi_grid_1D, &
                            aij)

            ! calculate energy and derivatives
            call calculate_energy(dx, dy, aij, E, dE_dphi, dE_dpsi)
            energy = energy + E


            ! calculate forces
            sin_phi = sin(phi)
            sin_psi = sin(psi)
            f_phi = dE_dphi/ sinphi
            f_psi = dE_dpsi/ sinpsi

            d_1   = sqrt( 4.0d0 * l_jk**2 * l_ij**2 - ( l_ij**2 + l_jk**2 - l_ik**2 )**2 )
            d_2   = sqrt( 4.0d0 * l_jk**2 * l_kl**2 - ( l_kl**2 + l_jk**2 - l_jl**2 )**2 )
            d_3   = sqrt( 4.0d0 * l_kl**2 * l_lm**2 - ( l_lm**2 + l_kl**2 - l_km**2 )**2 )

            d_phi = 2.0d0 * l_jk**2 * ( l_jl**2 + l_ik**2 - l_il**2 - l_jk**2 ) &
                     + ( l_ij**2 + l_jk**2 - l_ik**2 ) * ( l_jk**2 + l_kl**2 - l_jl**2 )
            d_psi = 2.0d0 * l_kl**2 * ( l_km**2 + l_jl**2 - l_jm**2 - l_kl**2 ) &
                     + ( l_jk**2 + l_kl**2 - l_jl**2 ) * ( l_kl**2 + l_lm**2 - l_km**2 )
            f_phi = f_phi * 4.0d0 / (d_1 * d_2) * RAD_TO_DEG / grid_step_size
            f_psi = f_psi * 4.0d0 / (d_2 * d_3) * RAD_TO_DEG / grid_step_size
            
            dphi_dra(:,:) = dphi_dra(:, :) * RAD_TO_DEG / grid_step_size
            dpsi_dra(:,:) = dpsi_dra(:, :) * RAD_TO_DEG / grid_step_size


            f_i =   -dE_dphi * dphi_dra(1, :)
            f_j =   -dE_dphi * dphi_dra(2, :) &
                  -  dE_dpsi * dpsi_dra(1, :)
            f_k =   -dE_dphi * dphi_dra(3, :) &
                  -  dE_dpsi * dpsi_dra(2, :)
            f_l =   -dE_dphi * dphi_dra(4, :) &
                  -  dE_dpsi * dpsi_dra(3, :)
            f_m =   -dE_dpsi * dpsi_dra(4, :)


            forces(iatm, :) = forces(iatm, :) + f_i(:)
            forces(jatm, :) = forces(jatm, :) + f_j(:)
            forces(katm, :) = forces(katm, :) + f_k(:)
            forces(latm, :) = forces(latm, :) + f_l(:)
            forces(matm, :) = forces(matm, :) + f_m(:)

            ! calculate two-body forces
            f_ij = f_phi/l_1 * ( dot_product(-r_jk, r_kl)/ l_2 &
                               - dot_product(r_ik, r_jk)*cos_phi/l_1 ) * r_ij
            
            f_ik = f_phi/l_1 * ( dot_product(r_jl, r_jk)/ l_2 &
                               - dot_product(r_ij, -r_jk)/l_2 ) * r_ik

            f_il = f_phi * -dot_product(r_jk, r_jk) / (l_1 * l_2) * r_il

            f_im = 0.0d0

            f_jk = ( f_phi * ( ( -dot_product(r_ik, r_jl) &
                                 -dot_product(r_ij, r_kl) &
                               )/(l_1*l_2) &
                             - dot_product(r_ij, r_ik)*cos_phi/(l_1*l_1) &
                             - dot_product(r_kl, r_jl)*cos_phi/(l_2*l_2) ) &

                   + f_psi/l_2 * ( dot_product(-r_kl, r_lm)/ l_3 &
                                 - dot_product(r_jk, -r_kl)*cos_psi/l_2 ) ) * r_jk

            f_jl = ( f_phi/l_2 * ( dot_product(r_jk, r_ik)/l_1 &
                                 - dot_product(r_jk, -r_kl)*cos_phi/l_2 ) &

                   + f_psi/l_2 * ( dot_product(r_km, r_kl)/l_3 &
                                 - dot_product(r_jk, -r_kl)*cos_psi/l_2 ) ) * r_jl

            f_jm = f_psi * -dot_product(r_kl, r_kl) / (l_2 * l_3) * r_jm

            f_kl = ( f_phi/l_2 * ( dot_product(r_ij, -r_jk)/l_1 &
                                 - dot_product(r_jk, r_jl)*cos_phi/l_2 ) &
                   
                   + f_psi * ( -dot_product(r_lm, r_jk) &
                               -dot_product(r_jl, r_km) &
                             )/(l_2*l_3) &
                             - dot_product(r_jk, r_jl)*cos_psi/(l_2*l_2) &
                             - dot_product(-r_kl, r_lm)*cos_psi/(l_3*l_3) ) * r_kl

            f_km = f_psi/l_3 * ( dot_product(r_kl, r_jl)/ l_2 &
                               - dot_product(-r_kl, r_lm)*cos_psi/l_3 ) * r_km

            f_lm = f_psi/l_3 * ( dot_product(r_jk, -r_kl)/ l_2 &
                               - dot_product(r_kl, r_km)*cos_psi/l_3 ) * r_lm

            f_ij = f_phi * ( dot_product(-r_jk, r_kl)  &
                            - (d_phi / d_1**2) * dot_product(-r_jk, -r_ik) ) * r_ij

            f_ik = f_phi * ( dot_product(-r_jk, -r_jl)  &
                            - (d_phi / d_1**2) * dot_product(r_ij, -r_jk) ) * r_ik

            f_il = -f_phi * dot_product(-r_jk, -r_jk) * r_il

            f_im = 0.0d0

            f_jk = ( f_phi * ( ( dot_product(-r_jl, r_ij) &
                                + dot_product(r_ik, -r_kl) - dot_product(-r_jk, -r_jk) ) &
                              - (d_phi / d_1**2) * dot_product(r_ij, r_ik) &
                              - (d_phi / d_2**2) * dot_product(-r_kl, -r_jl) ) &
                    + f_psi * ( dot_product(-r_kl, r_lm)  &
                               - ( d_psi / d_2**2 ) * dot_product(-r_kl, -r_jl) ) ) * r_jk
            
            f_jl = ( f_phi * ( dot_product(r_jk, r_ik)  &
                              - ( d_phi / d_2**2 ) * dot_product(-r_kl, r_jk) ) &
                    + f_psi * ( dot_product(-r_kl, -r_km)  &
                            - ( d_psi / d_2**2 ) * dot_product(r_jk, -r_kl) ) ) * r_jl

            f_jm = -f_psi * dot_product(-r_kl, -r_kl) * r_jm

            f_kl = ( f_phi * ( dot_product(r_ij, -r_jk) &
                              - ( d_phi / d_2**2 ) * dot_product(-r_jk, -r_jl) ) &
                    + f_psi * ( ( dot_product(-r_km, r_jk) &
                                + dot_product(r_jl, -r_lm) - dot_product(-r_kl, -r_kl) ) &
                               - (d_psi / d_2**2) * dot_product(r_jk, r_jl) &
                               - (d_psi / d_3**2) * dot_product(-r_lm, -r_km) ) ) * r_kl

            f_km = f_psi * ( dot_product(r_kl, r_jl)  &
                             - ( d_psi / d_3**2 ) * dot_product(-r_lm, r_kl) ) * r_km

            f_lm = f_psi * ( dot_product(r_jk, -r_kl) &
                             - ( d_psi / d_3**2 ) * dot_product(-r_kl, -r_km) ) * r_lm

            if (itbf_ij > 0) then
                tbforces(itbf_ij, :) = tbforces(itbf_ij, :) + f_ij(:)
                displacement(itbf_ij, :) = r_ij(:)
            else
                tbforces(-itbf_ij, :) = tbforces(-itbf_ij, :) - f_ij(:)
                displacement(-itbf_ij, :) = -r_ij(:)
            end if

            if (itbf_ik > 0) then
                tbforces(itbf_ik, :) = tbforces(itbf_ik, :) + f_ik(:)
                displacement(itbf_ik, :) = r_ik(:)
            else
                tbforces(-itbf_ik, :) = tbforces(-itbf_ik, :) - f_ik(:)
                displacement(-itbf_ik, :) = -r_ik(:)
            end if

            if (itbf_il > 0) then
                tbforces(itbf_il, :) = tbforces(itbf_il, :) + f_il(:)
                displacement(itbf_il, :) = r_il(:)
            else
                tbforces(-itbf_il, :) = tbforces(-itbf_il, :) - f_il(:)
                displacement(-itbf_il, :) = -r_il(:)
            end if

            if (itbf_im > 0) then
                tbforces(itbf_im,:) = tbforces(itbf_im,:) + f_im(:)
                displacement(itbf_im,:) = r_im(:)
            else
                tbforces(-itbf_im,:) = tbforces(-itbf_im,:) - f_im(:)
                displacement(-itbf_im,:) = -r_im(:)
            end if
            
            if (itbf_jk > 0) then
                tbforces(itbf_jk,:) = tbforces(itbf_jk,:) + f_jk(:)
                displacement(itbf_jk,:) = r_jk(:)
            else
                tbforces(-itbf_jk,:) = tbforces(-itbf_jk,:) - f_jk(:)
                displacement(-itbf_jk,:) = -r_jk(:)
            end if

            if (itbf_jl > 0) then
                tbforces(itbf_jl,:) = tbforces(itbf_jl,:) + f_jl(:)
                displacement(itbf_jl,:) = r_jl(:)
            else
                tbforces(-itbf_jl,:) = tbforces(-itbf_jl,:) - f_jl(:)
                displacement(-itbf_jl,:) = -r_jl(:)
            end if

            if (itbf_jm > 0) then
                tbforces(itbf_jm,:) = tbforces(itbf_jm,:) + f_jm(:)
                displacement(itbf_jm,:) = r_jm(:)
            else
                tbforces(-itbf_jm,:) = tbforces(-itbf_jm,:) - f_jm(:)
                displacement(-itbf_jm,:) = -r_jm(:)
            end if

            if (itbf_kl > 0) then
                tbforces(itbf_kl,:) = tbforces(itbf_kl,:) + f_kl(:)
                displacement(itbf_kl,:) = r_kl(:)
            else
                tbforces(-itbf_kl,:) = tbforces(-itbf_kl,:) - f_kl(:)
                displacement(-itbf_kl,:) = -r_kl(:)
            end if

            if (itbf_km > 0) then
                tbforces(itbf_km,:) = tbforces(itbf_km,:) + f_km(:)
                displacement(itbf_km,:) = r_km(:)
            else
                tbforces(-itbf_km,:) = tbforces(-itbf_km,:) - f_km(:)
                displacement(-itbf_km,:) = -r_km(:)
            end if

            if (itbf_lm > 0) then
                tbforces(itbf_lm,:) = tbforces(itbf_lm,:) + f_lm(:)
                displacement(itbf_lm,:) = r_lm(:)
            else
                tbforces(-itbf_lm,:) = tbforces(-itbf_lm,:) - f_lm(:)
                displacement(-itbf_lm,:) = -r_lm(:)
            end if

            if (check) then
                print*, 'TB_CHECK: i, j, k, l, m=', iatm, jatm, katm, latm, matm
                print*, 'TB_CHECK: f_i vs. f_ij + f_ik + f_il + f_im =>', &
                    & check_tbforce(f_i, f_ij+f_ik+f_il+f_im)
                print*, 'TB_CHECK:',f_i
                print*, 'TB_CHECK:',f_ij + f_ik + f_il + f_im
                print*, 'TB_CHECK: f_j vs. f_ji + f_jk + f_jl + f_jm =>', & 
                    & check_tbforce(f_j, -f_ij+f_jk+f_jl+f_jm)
                print*, 'TB_CHECK:',f_j
                print*, 'TB_CHECK:',- f_ij + f_jk + f_jl + f_jm
                print*, 'TB_CHECK: f_k vs. f_ki + f_kj + f_kl + f_km =>', &
                    & check_tbforce(f_k, -f_ik-f_jk+f_kl+f_km)
                print*, 'TB_CHECK:',f_k
                print*, 'TB_CHECK:',- f_ik - f_jk + f_kl + f_km
                print*, 'TB_CHECK: f_l vs. f_li + f_lj + f_lk + f_lm =>', &
                    & check_tbforce(f_l, -f_il-f_jl-f_kl+f_lm)
                print*, 'TB_CHECK:',f_l
                print*, 'TB_CHECK:',- f_il - f_jl - f_kl + f_lm
                print*, 'TB_CHECK: f_m vs. f_mi + f_mj + f_mk + f_ml =>', &
                    & check_tbforce(f_m, -f_im-f_jm-f_km-f_lm)
                print*, 'TB_CHECK:',f_m
                print*, 'TB_CHECK:',- f_im - f_jm - f_km - f_lm

            end if


        end do

    end subroutine

    subroutine cal_torsion(rab, rcb, rdc, dphi_dra, cosphi, sinphi)

        use common_vars

        implicit none
        real(8), intent(in)  :: rab(3), rcb(3), rdc(3)
        real(8), intent(out) :: dphi_dra(4,3), cosphi, sinphi

        real(8) :: e_cb(3), up_ab(3), up_dc(3), up_abc(3), up_bcd(3)
        real(8) :: tmp(3)
        real(8) :: d_cb, d_pab, d_pdc, dot_ab_cb, dot_dc_cb
        real(8) :: cosphi_pre, sinphi_pre

        d_cb = 1.0d0 / sqrt(dot_product(rcb, rcb))
        e_cb(1:3) = rcb(1:3) * d_cb
        
        dot_ab_cb = dot_product(rab, e_cb)
        up_ab(1:3) = rab(1:3) - dot_ab_cb * e_cb(1:3)
        d_pab = 1.0d0 / sqrt(dot_product(up_ab, up_ab))
        up_ab(1:3) = up_ab(1:3) * d_pab

        dot_dc_cb = dot_product(rdc, e_cb)
        up_dc(1:3) = rdc(1:3) - dot_dc_cb * e_cb(1:3)
        d_pdc = 1.0d0 / sqrt(dot_product(up_dc, up_dc))
        up_dc(1:3) = up_dc(1:3) * d_pdc

        cosphi_pre = dot_product(up_ab, up_dc)
        cosphi = min( max(cosphi_pre, -1.d0), 1.d0 )

        sinphi_pre = dot_product( outer_prod(up_ab, up_dc), e_cb)
        sinphi = min( max(sinphi_pre, -1.d0), 1.d0 )

        up_abc = outer_prod(up_ab, e_cb)
        up_bcd = outer_prod(e_cb, up_dc)

        dphi_dra(1,1:3) = up_abc(1:3) * d_pab
        dphi_dra(4,1:3) = up_bcd(1:3) * d_pdc

        tmp(1:3) = (dot_ab_cb * d_cb) * dphi_dra(1,1:3) + &
                   (dot_dc_cb * d_cb) * dphi_dra(4,1:3)
        dphi_dra(2,1:3) = tmp(1:3) - dphi_dra(1,1:3)
        dphi_dra(3,1:3) = -tmp(1:3) - dphi_dra(4,1:3)

    end subroutine cal_torsion

    integer function calculate_cmap_grid(cmap_type, value)
        implicit none
        integer, intent(in) :: cmap_type
        integer, intent(in) :: value
        integer :: resolution

        resolution = cmap_resolutions(cmap_type)
        calculate_cmap_grid = modulo(value-1, resolution) + 1

    end function calculate_cmap_grid

    subroutine cmap_coeff(grid_step_size, &
                          E_grid_1D, &
                          dphi_grid_1D, &
                          dpsi_grid_1D, &
                          dphidpsi_grid_1D, &
                          aij)
        implicit none
        integer, intent(in) :: grid_step_size
        !inputs
        real(8), dimension(4), intent(in)    :: E_grid_1D
        real(8), dimension(4), intent(in)    :: dphi_grid_1D
        real(8), dimension(4), intent(in)    :: dpsi_grid_1D
        real(8), dimension(4), intent(in)    :: dphidpsi_grid_1D
        !outputs
        real(8), dimension(4,4), intent(out) :: aij

        real(8), dimension(16)    :: a_tmp
        real(8), dimension(16)    :: E_grids
        real(8), dimension(16,16) :: weight
        integer, dimension(1:2)   :: shapes = (/4, 4/)
        integer, dimension(1:2)   :: orders = (/2, 1/)

        weight(1 ,1:16) =   (/ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        weight(2 ,1:16) =   (/ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0/)
        weight(3 ,1:16) =   (/-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1, 0, 0, 0, 0/)
        weight(4 ,1:16) =   (/ 2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0/)
        weight(5 ,1:16) =   (/ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        weight(6 ,1:16) =   (/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0/)
        weight(7 ,1:16) =   (/ 0, 0, 0, 0,-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1/)
        weight(8 ,1:16) =   (/ 0, 0, 0, 0, 2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1/)
        weight(9 ,1:16) =   (/-3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        weight(10,1:16) =   (/ 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0/)
        weight(11,1:16) =   (/ 9,-9, 9,-9, 6, 3,-3,-6, 6,-6,-3, 3, 4, 2, 1, 2/)
        weight(12,1:16) =   (/-6, 6,-6, 6,-4,-2, 2, 4,-3, 3, 3,-3,-2,-1,-1,-2/)
        weight(13,1:16) =   (/ 2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        weight(14,1:16) =   (/ 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0/)
        weight(15,1:16) =   (/-6, 6,-6, 6,-3,-3, 3, 3,-4, 4, 2,-2,-2,-2,-1,-1/)
        weight(16,1:16) =   (/ 4,-4, 4,-4, 2, 2,-2,-2, 2,-2,-2, 2, 1, 1, 1, 1/)

        E_grids(1:4) = E_grid_1D(1:4)
        E_grids(5:8) = dphi_grid_1D(1:4) * grid_step_size
        E_grids(9:12) = dpsi_grid_1D(1:4) * grid_step_size
        E_grids(13:16) = dphidpsi_grid_1D(1:4) * grid_step_size * grid_step_size

        a_tmp = matmul(weight, E_grids)
        aij = reshape(a_tmp, shapes, order=orders)

    end subroutine cmap_coeff

    subroutine calculate_energy(dx, dy, aij, E, dE_dphi, dE_dpsi)
        implicit none
        real(8), intent(in)                 :: dx, dy
        real(8), dimension(4,4), intent(in) :: aij
        real(8), intent(out) :: E
        real(8), intent(out) :: dE_dphi
        real(8), intent(out) :: dE_dpsi

        integer :: i
        E       = 0.0d0
        dE_dphi  = 0.0d0
        dE_dpsi  = 0.0d0

        do i=4, 1, -1

            E       = E*dx + ( ( aij(i,4)*dy + aij(i,3) )*dy + aij(i,2) )* dy + aij(i,1)
            dE_dphi = dE_dphi*dy + ( 3.0d0*aij(4,i)*dx + 2.0d0*aij(3,i) )*dx + aij(2,i)
            dE_dpsi = dE_dpsi*dx + ( 3.0d0*aij(i,4)*dy + 2.0d0*aij(i,3) )*dy + aij(i,2)

        end do

    end subroutine calculate_energy


end module cmap


!###############################################################################
module coulomb
    implicit none
    ! input
    real(8), allocatable :: charges(:) ! (natom)
    real(8) :: cutoff_length           ! cutoff length
    ! output
    real(8) :: energy                      ! enegry
    real(8), allocatable :: forces(:, :)   ! (natom, 3)
    real(8), allocatable :: tbforces(:, :) ! (max_tbf, 3)
    !(iatm_beg:iatm_end, iatm_beg:natom, 3)
    real(8), allocatable :: displacement(:, :) ! (ntbf, 3)

contains
    subroutine calculate(interact_table, ninteract)
        use total
        use common_vars
        implicit none
        integer, intent(in) :: ninteract
        integer, intent(in) :: interact_table(ninteract, 3)
                ! 1:iatm, 2:begin of jatm, 3:end of jatm
        real(8) :: coeff, cutoff_inv
        ! gained from PRESTO .  e^2 / (4 pi eps_0) [kcal/mol * A/eV^2]
        ! real(8), parameter :: coeff15 = 332.06378d0
        ! gained from Amebr difinition. ! e^2 / (4 pi eps_0) [kcal/mol * A/eV^2]
        real(8), parameter :: coeff15 = 332.05221729d0

        ! dielectric constant = 8.85418782... x 10^-12 [C^2/J M]

        ! for each nonbonded atom pairs,
        ! E_coulomb = 1/(4 pi eps_0)  q_i q_j / l_ij
        
        ! F_i = 1/(4 pi eps_0) q_i q_j / l_ij^3   r_ij
        ! F_j = - F_i

        ! F_ij = F_i
        ! F_ji = F_j

        ! initialization
        energy = 0.0d0
        forces = 0.0d0
        tbforces = 0.0d0
        displacement = 0.0d0

        coeff = coeff15
        cutoff_inv = 1.0d0/cutoff_length

        itbf = 0

        do jint=1, ninteract
            iatm     = interact_table(jint, 1)
            jatm_beg = interact_table(jint, 2)
            jatm_end = interact_table(jint, 3)

            do jatm=jatm_beg, jatm_end
                itbf = itbf + 1

                r_ij = crd(iatm, :) - crd(jatm, :)
                l_ij_inv = 1.0d0/sqrt( dot_product(r_ij, r_ij) )

                ! cutoff
                if (l_ij_inv < cutoff_inv) cycle

                ! calculate energy
                ! ene = coeff*charges(iatm)*charges(jatm) / l_ij
                ene = coeff*charges(iatm)*charges(jatm) * l_ij_inv
                energy = energy + ene

                ! calculate force
                ! f_i = ene * r_ij / (l_ij**2)
                f_i = ene * r_ij *l_ij_inv*l_ij_inv

                ! store force
                forces(iatm, :) = forces(iatm, :) + f_i(:)
                forces(jatm, :) = forces(jatm, :) - f_i(:)

                ! calculate and store two-body force and two-body distance vector
                f_ij = f_i
                tbforces(itbf, :) = f_ij(:)
                displacement(itbf, :) = r_ij(:)

                ! check
                if (check) then
                    print*, 'TB_CHECK: i, j =', iatm, jatm
                    print*, 'TB_CHECK: f_i vs. f_ij =>', check_tbforce(f_i, f_ij)
                    print*, 'TB_CHECK',f_i
                    print*, 'TB_CHECK',f_ij
                    print*, 'TB_CHECK: f_j vs. f_ji =>', check_tbforce(-f_i, -f_ij)
                    print*, 'TB_CHECK',-f_i
                    print*, 'TB_CHECK',-f_ij
                end if

            enddo

        enddo

    end subroutine

end module

module coulomb14
    implicit none
    ! input
    real(8), allocatable :: charges(:)     ! (natom)
    integer, allocatable :: i14_to_itbf(:) ! (n14)
    ! output
    real(8) :: energy                      ! enegry
    real(8), allocatable :: forces(:, :)   ! (natom, 3)
    real(8), allocatable :: tbforces(:, :) ! (ntbf, 3)
    real(8), allocatable :: displacement(:, :)! (ntbf, 3)

    ! parameter
    real(8) :: scale_factor = 1.2d0    ! scaling factor for 1-4 interaction
    
    ! main variales
    integer :: i14, n14, itbf_ij
    integer :: iatm, jatm

contains
    subroutine calculate()
        use total
        use common_vars
        implicit none
        real(8) :: coeff
        ! gained from PRESTO .  e^2 / (4 pi eps_0) [kcal/mol * A/eV^2]
        ! real(8), parameter :: coeff15 = 332.06378d0
        ! gained from Amebr difinition. ! e^2 / (4 pi eps_0) [kcal/mol * A/eV^2]
        real(8), parameter :: coeff15 = 332.05221729d0
        ! dielectric constant = 8.85418782... x 10^-12 [C^2/J M]

        ! for each nonbonded atom pairs,
        ! E_coulomb = 1/(4 pi eps_0)  q_i q_j / l_ij
        
        ! F_i = 1/(4 pi eps_0) q_i q_j / l_ij^3   r_ij
        ! F_j = - F_i

        ! F_ij = F_i
        ! F_ji = F_j

        ! initialization
        energy = 0.0d0
        forces = 0.0d0
        tbforces = 0.0d0
        displacement = 0.0d0

        coeff = coeff15/scale_factor
        n14   = size(i14_to_itbf)

        do i14=1, n14
            itbf_ij = i14_to_itbf(i14)
            iatm = bonded_pairs(abs(itbf_ij), 1)
            jatm = bonded_pairs(abs(itbf_ij), 2)

            r_ij = crd(iatm, :) - crd(jatm, :)
            l_ij = sqrt( dot_product(r_ij, r_ij) )

            ! calculate energy
            ene = coeff*charges(iatm)*charges(jatm) / l_ij
            energy = energy + ene

            ! calculate force
            f_i = ene * r_ij / (l_ij**2)

            ! store force
            forces(iatm, :) = forces(iatm, :) + f_i(:)
            forces(jatm, :) = forces(jatm, :) - f_i(:)

            ! calculate two-body force
            f_ij = f_i

            ! store two-body force and two-body distance vector
            if (itbf_ij > 0) then
                tbforces(itbf_ij,:) = tbforces(itbf_ij,:) + f_ij(:)
                displacement(itbf_ij,:) = r_ij(:)
            else
                tbforces(-itbf_ij,:) = tbforces(-itbf_ij,:) - f_ij(:)
                displacement(-itbf_ij,:) = -r_ij(:)
            end if

            ! check
            if (check) then
                print*, 'TB_CHECK: i, j =', iatm, jatm
                print*, 'TB_CHECK: f_i vs. f_ij =>', check_tbforce(f_i, f_ij)
                print*, 'TB_CHECK',f_i
                print*, 'TB_CHECK',f_ij
                print*, 'TB_CHECK: f_j vs. f_ji =>', check_tbforce(-f_i, -f_ij)
                print*, 'TB_CHECK',-f_i
                print*, 'TB_CHECK',-f_ij
            end if

        enddo

        ! print*, 'c14'
        ! do iatm=1, natom
            ! print*, iatm, forces(iatm,:)
        ! end do
        
    end subroutine

    subroutine print_tbforce()
        use total, only: bonded_pairs, nbonded
        use common_vars
        implicit none
        real(8) :: tbf_abs

        do itbf=1, nbonded
            iatm = bonded_pairs(itbf, 1)
            jatm = bonded_pairs(itbf, 2)

            tbf_abs = dot_product(tbforces(itbf,:), tbforces(itbf,:))
            if (-0.0001 < tbf_abs .and. tbf_abs < 0.0001) cycle
            write(*, '(2I5,3F15.7)'), iatm, jatm, tbforces(itbf_ij, :)
         end do
            
    end subroutine

end module

!###############################################################################
module vdw
    implicit none
    ! input
    ! integer, allocatable :: interact_table(:,:) ! (nint, 3)
    !         ! 1:iatm, 2:begin of jatm, 3:end of jatm
    integer, allocatable :: atom_types(:)  ! (natom)
    real(8), allocatable :: c6s(:,:)       ! (num_atomtypes, num_atomtypes)
    real(8), allocatable :: c12s(:,:)      ! (num_atomtypes, num_atomtypes)
    real(8) :: cutoff_length               ! cutoff length
    ! output
    real(8) :: energy
    real(8), allocatable :: forces(:, :)   ! (natom, 3)
    real(8), allocatable :: tbforces(:, :) ! (max_tbf, 3)
    !(iatm_beg:iatm_end, iatm_beg:natom, 3)
    real(8), allocatable :: displacement(:, :)! (ntbf, 3)

contains
    subroutine calculate(interact_table, ninteract)
        use total
        use common_vars
        implicit none
        integer, intent(in) :: ninteract
        integer, intent(in) :: interact_table(ninteract, 3)
                ! 1:iatm, 2:begin of jatm, 3:end of jatm
        real(8) :: c6, c12, cutoff_len2, coeff
        real(8) :: l_ij2
        real(8), parameter :: coeff15 = 1.0

        ! for each nonbonded atom pairs,
        ! E_vdw = eps_ij ( A_ij/l_ij^12 - B_ij/l_ij^6)
        !       = C12/l_ij^12 - C6/l_ij^6
        
        ! F_i = (-12*C12/l_ij^13 + 6*C6/l_ij^7 ) r_ij/l_ij
        ! F_j = - F_i

        ! F_ij = F_i
        ! F_ji = F_j

        ! initialization
        energy = 0.0d0
        forces = 0.0d0
        tbforces = 0.0d0
        displacement = 0.0d0
        cutoff_len2 = cutoff_length * cutoff_length
        coeff = coeff15

        itbf = 0
        do jint=1, ninteract
            iatm     = interact_table(jint, 1)
            jatm_beg = interact_table(jint, 2)
            jatm_end = interact_table(jint, 3)

            do jatm=jatm_beg, jatm_end
                itbf = itbf + 1

                r_ij = crd(iatm, :) - crd(jatm, :)
                l_ij2 = dot_product(r_ij, r_ij)

                ! cutoff
                if (l_ij2 > cutoff_len2) cycle

                c6   = c6s( atom_types(iatm), atom_types(jatm) )
                c12  = c12s( atom_types(iatm), atom_types(jatm) )

                ! calculate energy
                ene = coeff*(c12/l_ij2**6 - c6/l_ij2**3)
                energy = energy + ene

                ! calculate force
                f_i = - coeff*(-12.0d0*c12/l_ij2**7 + 6.0d0*c6/l_ij2**4 ) * r_ij

                ! store force
                forces(iatm, :) = forces(iatm, :) + f_i(:)
                forces(jatm, :) = forces(jatm, :) - f_i(:)

                ! calculate and store two-body force and two-body distance vector
                f_ij = f_i
                tbforces(itbf, :) = f_ij(:)
                displacement(itbf, :) = r_ij(:)

                ! check
                if (check) then
                    print*, 'TB_CHECK: i, j =', iatm, jatm
                    print*, 'TB_CHECK: f_i vs. f_ij =>', check_tbforce(f_i, f_ij)
                    print*, 'TB_CHECK',f_i
                    print*, 'TB_CHECK',f_ij
                    print*, 'TB_CHECK: f_j vs. f_ji =>', check_tbforce(-f_i, -f_ij)
                    print*, 'TB_CHECK',-f_i
                    print*, 'TB_CHECK',-f_ij
                end if

            enddo

        enddo

    end subroutine

end module


!###############################################################################
module vdw14
    implicit none
    ! input
    integer, allocatable :: atom_types(:)  ! (natom)
    real(8), allocatable :: c6s(:,:)       ! (num_atomtypes, num_atomtypes)
    real(8), allocatable :: c12s(:,:)      ! (num_atomtypes, num_atomtypes)
    integer, allocatable :: i14_to_itbf(:) ! (n14)
    ! output
    real(8) :: energy
    real(8), allocatable :: forces(:, :)   ! (natom, 3)
    real(8), allocatable :: tbforces(:, :) ! (ntbf, 3)
    real(8), allocatable :: displacement(:, :)! (ntbf, 3)

    ! parameter
    real(8) :: scale_factor = 2.0d0 ! scaling factor for 1-4 interaction

    ! main variales
    integer :: i14, n14, itbf_ij
    integer :: iatm, jatm

contains
    subroutine calculate()
        use total
        use common_vars
        implicit none
        real(8) :: c6, c12, l_ij2, coeff
        real(8), parameter :: coeff15 = 1.0

        ! for each nonbonded atom pairs,
        ! E_vdw = eps_ij ( A_ij/l_ij^12 - B_ij/l_ij^6)
        !       = C12/l_ij^12 - C6/l_ij^6
        
        ! F_i = (-12*C12/l_ij^13 + 6*C6/l_ij^7 ) r_ij/l_ij
        ! F_j = - F_i

        ! F_ij = F_i
        ! F_ji = F_j

        ! initialization
        energy = 0.0d0
        forces = 0.0d0
        tbforces = 0.0d0
        displacement = 0.0d0

        coeff = coeff15/scale_factor
        n14 = size(i14_to_itbf)

        do i14=1, n14
            itbf_ij = i14_to_itbf(i14)
            iatm = bonded_pairs(abs(itbf_ij), 1)
            jatm = bonded_pairs(abs(itbf_ij), 2)

            r_ij = crd(iatm, :) - crd(jatm, :)
            l_ij2 = dot_product(r_ij, r_ij)

            c6   = c6s( atom_types(iatm), atom_types(jatm) )
            c12  = c12s( atom_types(iatm), atom_types(jatm) )

            ! calculate energy
            ene = coeff*(c12/l_ij2**6 - c6/l_ij2**3)
            energy = energy + ene

            ! calculate force
            f_i =  - coeff*(- 12.0d0*c12/l_ij2**7 + 6.0d0*c6/l_ij2**4 ) * r_ij

            ! store force
            forces(iatm, :) = forces(iatm, :) + f_i(:)
            forces(jatm, :) = forces(jatm, :) - f_i(:)

            ! calculate two-body force
            f_ij = f_i

            ! store two-body force and two-body distance vector
            if (itbf_ij > 0) then
                tbforces(itbf_ij,:) = tbforces(itbf_ij,:) + f_ij(:)
                displacement(itbf_ij,:) = r_ij(:)
            else
                tbforces(-itbf_ij,:) = tbforces(-itbf_ij,:) - f_ij(:)
                displacement(-itbf_ij,:) = -r_ij(:)
            end if

            ! check
            if (check) then
                print*, 'TB_CHECK: i, j =', iatm, jatm
                print*, 'TB_CHECK: f_i vs. f_ij =>', check_tbforce(f_i, f_ij)
                print*, 'TB_CHECK',f_i
                print*, 'TB_CHECK',f_ij
                print*, 'TB_CHECK: f_j vs. f_ji =>', check_tbforce(-f_i, -f_ij)
                print*, 'TB_CHECK',-f_i
                print*, 'TB_CHECK',-f_ij
            end if

        enddo

    end subroutine

    subroutine print_tbforce()
        use total, only: bonded_pairs, nbonded
        use common_vars
        implicit none
        real(8) :: tbf_abs

        do itbf=1, nbonded
            iatm = bonded_pairs(itbf, 1)
            jatm = bonded_pairs(itbf, 2)

            tbf_abs = dot_product(tbforces(itbf,:), tbforces(itbf,:))
            if (-0.0001 < tbf_abs .and. tbf_abs < 0.0001) cycle
            write(*, '(2I5,3F15.7)'), iatm, jatm, tbforces(itbf_ij, :)
         end do
            
    end subroutine

end module

! subroutine cal_distance_restraint()
! end

!###############################################################################
subroutine setup(natom, check, bonded_pairs, nbonded, max_tbf)
    ! call at the beginning of the all snap shots.
    use total, &
            & t_natom             => natom           ,  & 
            & t_crd               => crd             ,  & 
            & t_nbonded           => nbonded         ,  & 
            & t_bonded_pairs      => bonded_pairs    ,  & 
            & t_bonded_tbforces   => bonded_tbforces ,  & 
            & t_check             => check


    use common_vars
    use bond      , bnd_forces => forces, bnd_tbforces => tbforces, bnd_disp => displacement
    use angle     , ang_forces => forces, ang_tbforces => tbforces, ang_disp => displacement
    use torsion   , tor_forces => forces, tor_tbforces => tbforces, tor_disp => displacement
    use improper  , imp_forces => forces, imp_tbforces => tbforces, imp_disp => displacement
    use cmap      , cmp_forces => forces, cmp_tbforces => tbforces, cmp_disp => displacement
    use coulomb   , cou_forces => forces, cou_tbforces => tbforces, cou_disp => displacement
    use vdw       , vdw_forces => forces, vdw_tbforces => tbforces, vdw_disp => displacement
    use coulomb14 , cou14_forces => forces, cou14_tbforces => tbforces, cou14_disp => displacement
    use vdw14     , vdw14_forces => forces, vdw14_tbforces => tbforces, vdw14_disp => displacement

    implicit none
    integer, intent(in) :: natom
    integer, intent(in) :: nbonded
    integer, intent(in) :: max_tbf
    logical, intent(in) :: check
    integer, intent(in) :: bonded_pairs(nbonded, 2)

    ! allocate
    t_natom = natom
    t_check = check
    t_nbonded = size(bonded_pairs, 1)
    allocate(t_crd(t_natom, 3))
    allocate(t_bonded_pairs(t_nbonded, 2))
    allocate(t_bonded_tbforces(t_nbonded, 3))
    t_bonded_pairs   = bonded_pairs

    ! allocate each force
    allocate(bnd_forces(t_natom, 3))
    allocate(ang_forces(t_natom, 3))
    allocate(tor_forces(t_natom, 3))
    allocate(imp_forces(t_natom, 3))
    allocate(cmp_forces(t_natom, 3))
    allocate(cou_forces(t_natom, 3))
    allocate(vdw_forces(t_natom, 3))
    allocate(cou14_forces(t_natom, 3))
    allocate(vdw14_forces(t_natom, 3))

    ! allocate each two-body force for bonded
    allocate(bnd_tbforces(t_nbonded, 3))
    allocate(ang_tbforces(t_nbonded, 3))
    allocate(tor_tbforces(t_nbonded, 3))
    allocate(imp_tbforces(t_nbonded, 3))
    allocate(cmp_tbforces(t_nbonded, 3))
    allocate(cou14_tbforces(t_nbonded, 3))
    allocate(vdw14_tbforces(t_nbonded, 3))

    ! allocate each two-body force for nonbonded
    allocate(cou_tbforces(max_tbf, 3))
    allocate(vdw_tbforces(max_tbf, 3))

    ! allocate each two-body distance (crd_i - crd_j) matrix
    allocate(bnd_disp(t_nbonded, 3))
    allocate(ang_disp(t_nbonded, 3))
    allocate(tor_disp(t_nbonded, 3))
    allocate(imp_disp(t_nbonded, 3))
    allocate(cmp_disp(t_nbonded, 3))
    allocate(cou_disp(max_tbf, 3))
    allocate(vdw_disp(max_tbf, 3))
    allocate(cou14_disp(t_nbonded, 3))
    allocate(vdw14_disp(t_nbonded, 3))

end subroutine

subroutine initialize(crd, natom)
    ! call each the snapshot
    use total, &
         & t_natom             => natom             , &
         & t_crd               => crd               , & 
         & t_bonded_tbforces   => bonded_tbforces
    implicit none
    integer, intent(in) :: natom
    real(8), intent(in) :: crd(natom, 3)

    t_crd = crd
    t_bonded_tbforces = 0.0d0

end subroutine

! subroutine get_bonded_tbforce()
    ! use total, &
          ! & t_natom             => natom             ,  & 
          ! & t_nbonded      => nbonded      ,  & 
          ! & t_bonded_pairs => bonded_pairs ,  & 
          ! & t_bonded_tbforces   => bonded_tbforces
    ! use bond      , bnd_tbforces => tbforces
    ! use angle     , ang_tbforces => tbforces
    ! use torsion   , tor_tbforces => tbforces
    ! use improper  , imp_tbforces => tbforces
    ! use common_vars
    ! implicit none
    ! real(8) :: total_tbforces(nbonded, 3)

    ! ! allocate(total_tbforces(t_natom, t_nextent, 3))
    ! total_tbforces = bnd_tbforces + ang_tbforces + tor_tbforces + imp_tbforces

! end subroutine

