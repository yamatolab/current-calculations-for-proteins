
subroutine initialize__(filename, formlabel, fileno)
    implicit none
    character(*), intent(in) :: filename
    character(*), intent(in) :: formlabel
    integer, intent(out) :: fileno
    integer :: ierr, ifile=10
    logical :: op

    ! get unused file number
    do
        inquire(ifile, opened=op)
        if (.not. op) exit
        ifile = ifile + 1
    end do
    fileno = ifile

    open(fileno, file=filename, status='old', form=formlabel, iostat=ierr)

    ! check
    if (ierr/= 0) then
        write(*,*) 'Could not create a file: ', filename
        stop
    end if
end subroutine

!###############################################################################
module formatted_coordinate
    implicit none
    integer :: fileno

contains
    subroutine initialize(filename)
        implicit none
        character(*), intent(in) :: filename
        call initialize__(filename, 'formatted', fileno)
    end subroutine

    subroutine finalize()
        implicit none
        close(fileno)
    end subroutine

    subroutine parse_next(natom, crd, ierr)
        implicit none
        integer, intent(in) :: natom
        real(8), intent(out):: crd(natom, 3)
        integer, intent(out):: ierr 
        integer :: iatm
        read(fileno, '(10F8.3)', iostat=ierr) &
            & (crd(iatm,1), crd(iatm,2), crd(iatm,3) ,iatm=1, natom)
    end subroutine

    subroutine parse_pbc(box, ierr)
        implicit none
        real(8), intent(out) :: box(3)
        integer, intent(out) :: ierr
        read(fileno, '(10F8.3)', iostat=ierr) box(1), box(2), box(3)
    end subroutine
end module

!###############################################################################
module unformatted_coordinate
    implicit none
    integer :: fileno

contains
    subroutine initialize(filename)
        implicit none
        character(*), intent(in) :: filename
        call initialize__(filename, 'unformatted', fileno)
    end subroutine

    subroutine finalize()
        implicit none
        close(fileno)
    end subroutine

    subroutine parse_next(natom, crd, ierr)
        implicit none
        integer, intent(in) :: natom
        real(8), intent(out):: crd(natom, 3)
        integer, intent(out):: ierr 
        integer :: iatm
        read(fileno, iostat=ierr) &
            & (crd(iatm,1), crd(iatm,2), crd(iatm,3) ,iatm=1, natom)
    end subroutine

    subroutine parse_pbc(box, ierr)
        implicit none
        real(8), intent(out) :: box(3)
        integer, intent(out) :: ierr
        read(fileno, iostat=ierr) box(1), box(2), box(3)
    end subroutine
end module

!###############################################################################
module formatted_velocity
    implicit none
    integer :: fileno

contains
    subroutine initialize(filename)
        implicit none
        character(*), intent(in) :: filename
        call initialize__(filename, 'formatted', fileno)
    end subroutine

    subroutine finalize()
        implicit none
        close(fileno)
    end subroutine

    subroutine parse_next(natom, vel, ierr)
        implicit none
        integer, intent(in) :: natom
        real(8), intent(out):: vel(natom, 3)
        integer, intent(out):: ierr 
        integer :: iatm
        read(fileno, '(10F8.3)', iostat=ierr) &
            & (vel(iatm,1), vel(iatm,2), vel(iatm,3) ,iatm=1, natom)
    end subroutine

    subroutine parse_pbc(box, ierr)
        implicit none
        real(8), intent(out) :: box(3)
        integer, intent(out) :: ierr
        read(fileno, '(10F8.3)', iostat=ierr) box(1), box(2), box(3)
    end subroutine
end module

!###############################################################################
module unformatted_velocity
    implicit none
    integer :: fileno

contains
    subroutine initialize(filename)
        implicit none
        character(*), intent(in) :: filename
        call initialize__(filename, 'unformatted', fileno)
    end subroutine

    subroutine finalize()
        implicit none
        close(fileno)
    end subroutine

    subroutine parse_next(natom, vel, ierr)
        implicit none
        integer, intent(in) :: natom
        real(8), intent(out):: vel(natom, 3)
        integer, intent(out):: ierr 
        integer :: iatm
        read(fileno, iostat=ierr) &
            & (vel(iatm,1), vel(iatm,2), vel(iatm,3) ,iatm=1, natom)
    end subroutine

    subroutine parse_pbc(box, ierr)
        implicit none
        real(8), intent(out) :: box(3)
        integer, intent(out) :: ierr
        read(fileno, iostat=ierr) box(1), box(2), box(3)
    end subroutine
end module

