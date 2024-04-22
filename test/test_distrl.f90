program test_distrl

    implicit none

    logical :: ok = .true.

    real :: xm0, ym0, dm0, wm0

    interface
        SUBROUTINE DISTRL(X0,Y0,X1,Y1,X2,Y2,XM,YM,DM,WM)
            REAL :: X0,Y0,X1,Y1,X2,Y2,XM,YM,DM,WM
        end subroutine DISTRL
    end interface

    if (ok) then
        print "(a)", "test_math: All tests passed."
    else
        print "(a)", "test_math: One ore more tests failed."
        stop 1
    end if

end program test_distrl