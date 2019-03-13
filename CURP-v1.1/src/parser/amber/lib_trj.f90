subroutine write_trajectory(snap, trj_fn, natom)
    implicit none

    integer, intent(in) :: natom
    real(8), intent(in) :: snap(natom, 3)
    character(*), intent(in) :: trj_fn

    integer :: iatm, ixyz
    integer, parameter :: iunit = 11

    ! coorditate trajectory
    open(unit=iunit, file=trj_fn, position='append', form='formatted')
    write(iunit,'(10F8.3)') ((snap(iatm,ixyz),ixyz=1,3), iatm=1,natom)
    close(iunit)

end subroutine
