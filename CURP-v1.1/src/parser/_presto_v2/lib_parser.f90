module io_manager
    implicit none
    integer :: iounit_int = 10

contains

     integer function io_unit()
        implicit none

        iounit_int = iounit_int + 1
        io_unit = iounit_int
        
    end function
    
end module

!###############################################################################
module unformatted_velocity
    implicit none
    integer :: fileno

contains
    subroutine initialize(filename)
        use io_manager
        implicit none
        character(*), intent(in) :: filename
        integer :: ierr

        fileno = io_unit()

        open(fileno, file=filename, status='old', form='unformatted', & 
            & iostat=ierr)

        ! check
        if (ierr/= 0) then
            write(*,*) 'Could not create a file: ', filename
            stop
        end if

    end subroutine

    subroutine finalize()
        implicit none
        close(fileno)
    end subroutine

    subroutine parse_next(natom, vel)
        implicit none
        integer, intent(in) :: natom
        real(8), intent(out):: vel(natom, 3)
        integer :: iloop
        real(8) :: simtime, cputime, energy, kinetic, temperature, potential, rmsf
        integer :: num_vdw_pair, num_hydrogen_pair
        real(8) :: rmsd

        read(fileno) iloop, simtime, cputime, &
            & energy, kinetic, temperature, potential, rmsf, &
            & num_vdw_pair, num_hydrogen_pair, rmsd

        read(fileno) vel
        
    end subroutine

end module

!###############################################################################
module unformatted_coordinate
    implicit none
    integer :: fileno

contains
    subroutine initialize(filename)
        use io_manager
        implicit none
        character(*), intent(in) :: filename
        integer :: ierr

        fileno = io_unit()

        open(fileno, file=filename, status='old', form='unformatted', & 
            & iostat=ierr)

        ! check
        if (ierr/= 0) then
            write(*,*) 'Could not create a file: ', filename
            stop
        end if

    end subroutine

    subroutine finalize()
        implicit none
        close(fileno)
    end subroutine

    subroutine parse_next(natom, crd)
        implicit none
        integer, intent(in) :: natom
        real(8), intent(out):: crd(natom, 3)
        integer :: iloop
        real(8) :: simtime, cputime, energy, kinetic, temperature, potential, rmsf
        integer :: num_vdw_pair, num_hydrogen_pair
        real(8) :: rmsd

        read(fileno) iloop, simtime, cputime, &
            & energy, kinetic, temperature, potential, rmsf, &
            & num_vdw_pair, num_hydrogen_pair, rmsd

        read(fileno) crd
        
    end subroutine

end module

!###############################################################################
module unformatted_restart
    implicit none
    integer :: fileno

contains
    subroutine initialize(filename)
        use io_manager
        implicit none
        character(*), intent(in) :: filename
        integer :: ierr

        fileno = io_unit()

        open(fileno, file=filename, status='old', form='unformatted', & 
            & iostat=ierr)

        ! check
        if (ierr/= 0) then
            write(*,*) 'Could not create a file: ', filename
            stop
        end if

    end subroutine

    subroutine finalize()
        implicit none
        close(fileno)
    end subroutine

    subroutine parse1(title, iloop, simtime, energy, kinetic, potential, ierr)
        implicit none
        character(*), intent(out) :: title
        integer, intent(out) :: iloop
        real(8), intent(out) :: simtime, energy, kinetic, potential
        integer, intent(out) :: ierr

        integer :: num_atom, num_var

        read(fileno, iostat=ierr) title
        read(fileno, iostat=ierr) num_atom, num_var
        read(fileno, iostat=ierr) iloop, simtime, energy, kinetic, potential

        ! if (ierr/= 0) then
        !     write(*,*) 'Could not create a file: ', filename
        !     stop
        ! end if

    end subroutine

    subroutine parse_crd(natom, crd, ierr)
        implicit none
        ! real(8), allocatable, intent(out):: crd(:,:)
        integer, intent(in) :: natom
        real(8), intent(out):: crd(natom,3)
        integer, intent(out):: ierr

        read(fileno, iostat=ierr) crd

    end subroutine

    subroutine parse_vel(natom, vel, ierr)
        implicit none
        integer, intent(in) :: natom
        real(8), intent(out):: vel(natom,3)
        integer, intent(out):: ierr

        read(fileno, iostat=ierr) vel
    end subroutine

end module


! subroutine parse_restart()
! end subroutine
