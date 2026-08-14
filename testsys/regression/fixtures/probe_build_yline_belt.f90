! Standalone driver for src/func_lib.f90's build_yline_belt, compiled and
! linked against src/globalvar.f90 + src/func_lib.f90 only -- no MPI, no
! MUMPS, no case setup. Reads a tiny geometry from stdin so one compiled
! binary can drive both the refuse case and the commensurate (working) case.
! See testsys/regression/test_fault_normal_offset_meshing.py.
!
! stdin format:
!   ntotft dx
!   y(1)
!   y(2)
!   ...
!   y(ntotft)
!
! On success, prints "NYT <n>" then one "Y <i> <value>" line per node line.
! On refusal, build_yline_belt itself prints the boxed MESH ERROR and calls
! STOP with a non-zero code, which is the thing under test.
program probe_build_yline_belt

    use globalvar
    implicit none

    real (kind = dp), allocatable :: ylinet(:)
    integer (kind = 4) :: nyt_local, i

    interface
        subroutine build_yline_belt(ylinet, nyt)
            use globalvar, only: dp
            real (kind = dp), allocatable, intent(out) :: ylinet(:)
            integer (kind = 4), intent(out) :: nyt
        end subroutine build_yline_belt
    end interface

    read(*,*) ntotft, dx
    allocate(fltxyz(2,4,ntotft))
    do i = 1, ntotft
        read(*,*) fltxyz(1,2,i)
        fltxyz(2,2,i) = fltxyz(1,2,i)   ! degenerate box: one y-plane per fault
    enddo

    dis4uniF = 5
    dis4uniB = 5
    rat = 1.3d0
    dymax = min(12.0d0*dx, 3.0d3)
    ymin = -35000.0d0
    ymax = 30000.0d0

    call build_yline_belt(ylinet, nyt_local)

    write(*,'(A,I0)') 'NYT ', nyt_local
    do i = 1, nyt_local
        write(*,'(A,I0,A,F14.4)') 'Y ', i, ' ', ylinet(i)
    enddo

end program probe_build_yline_belt
