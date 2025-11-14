module cmap
    implicit none
    ! input
    integer, allocatable :: five_atoms(:, :)            ! (ncmap, 5)
    integer, allocatable :: cmap_types(:)               ! (ncmap)
    integer, allocatable :: cmap_grid_step_size(:)      ! (ncmap)
    integer, allocatable :: cmap_resolution(:)          ! (ntypes)
    real(8), allocatable :: cmap_grid_energy(:, :, :)   ! (ncmap, ngrid_x(24), ngrid_y(24))
    real(8), allocatable :: cmap_grid_dphi(:, :, :)     ! (ncmap, ngrid_x(24), ngrid_y(24))
    real(8), allocatable :: cmap_grid_dpsi(:, :, :)     ! (ncmap, ngrid_x(24), ngrid_y(24))
    real(8), allocatable :: cmap_grid_dphi_dpsi(:, :, :)! (natom, ngrid_x(24), ngrid_y(24))
    integer, allocatable :: icmap_to_itbf(:, :)         ! (ncmap, 10)
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
        real(8) :: n_1(3), n_2(3), n_3(3), l_1, l_2, l_3
        real(8) :: judge_phi, judge_psi
        real(8) :: dx, dy
        real(8) :: E, dE_dphi, dE_dpsi
        real(8) :: f_phi, f_psi
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
            itbf_ij = icmap_to_itbf(icmap, 1)
            itbf_ik = icmap_to_itbf(icmap, 2)
            itbf_il = icmap_to_itbf(icmap, 3)
            itbf_im = icmap_to_itbf(icmap, 4)
            itbf_jk = icmap_to_itbf(icmap, 5)
            itbf_jl = icmap_to_itbf(icmap, 6)
            itbf_jm = icmap_to_itbf(icmap, 7)
            itbf_kl = icmap_to_itbf(icmap, 8)
            itbf_km = icmap_to_itbf(icmap, 9)
            itbf_lm = icmap_to_itbf(icmap, 10)

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

            ! calculate dihedral angle phi
            cos_phi = dot_product(n_1, n_2) / (l_1 * l_2)
            if (cos_phi >  1.0d0) cos_phi =  1.0d0
            if (cos_phi < -1.0d0) cos_phi = -1.0d0

            phi = acos(cos_phi)
            judge_phi = -dot_product(-r_jk, outer_prod(n_1, n_2))
            if (judge_phi < 0.0d0) phi = 2.0d0 * PI - phi
            deg_phi = phi * RAD_TO_DEG

            ! calculate dihedral angle psi
            cos_psi = dot_product(n_2, n_3) / (l_2 * l_3)
            if (cos_psi >  1.0d0) cos_psi =  1.0d0
            if (cos_psi < -1.0d0) cos_psi = -1.0d0

            psi = acos(cos_psi)
            judge_psi = -dot_product(-r_kl, outer_prod(n_2, n_3))
            if (judge_psi < 0.0d0) psi = 2.0d0 * PI - psi
            deg_psi = psi * RAD_TO_DEG

            ! determine grid indices for phi and psi
            grid_step_size = cmap_grid_step_size(cmap_type)
            x = int( (deg_phi - gridstart) / grid_step_size ) + 1
            y = int( (deg_psi - gridstart) / grid_step_size ) + 1

            ! compute fractional positions inside grid cell
            dx = modulo(deg_phi - dble(gridstart), dble(grid_step_size)) / dble(grid_step_size)

            dy = modulo(deg_psi - dble(gridstart), dble(grid_step_size)) / dble(grid_step_size)

            ! get the 2x2 grid energies and derivatives
            do i=1,2
                do j=1,2
                    E_grid(i,j)        = cmap_grid_energy(cmap_type, calculate_cmap_grid(cmap_type, x+i-1), &
                                                            calculate_cmap_grid (cmap_type, y+j-1) )
                    dphi_grid(i,j)     = cmap_grid_dphi(cmap_type, calculate_cmap_grid(cmap_type, x+i-1), &
                                                            calculate_cmap_grid (cmap_type, y+j-1) )
                    dpsi_grid(i,j)     = cmap_grid_dpsi(cmap_type, calculate_cmap_grid(cmap_type, x+i-1), &
                                                            calculate_cmap_grid (cmap_type, y+j-1) )
                    dphidpsi_grid(i,j) = cmap_grid_dphi_dpsi(cmap_type, calculate_cmap_grid(cmap_type, x+i-1), &
                                                            calculate_cmap_grid (cmap_type, y+j-1) )
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
            f_phi = -dE_dphi / sin_phi
            f_psi = -dE_dpsi / sin_psi

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

        end do

    end subroutine

    function calculate_cmap_grid(cmap_type, value)
        implicit none
        integer, intent(in) :: cmap_type
        integer, intent(in) :: value
        integer :: resolution
        integer :: cmap_grid

        resolution = cmap_resolution(cmap_type)
        cmap_grid = modulo(value-1, resolution) + 1
        calculate_cmap_grid = cmap_grid

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
