module time_functions
    !! Various procedures for manipulating time

    use uemep_constants, only: dp, secphour, secpday

    implicit none
    private

    public :: get_sun_angles, number_to_date, date_to_number, &
        date_to_datestr_bracket, date_to_datestr_squarebracket, date_to_datestr, &
        datestr_to_date, day_of_week, summer_time_europe, date_to_julian

contains

    subroutine number_to_date(date_num, date_array, ref_year)
        !! Returns array with date and time from a date number
        !!
        !! 'date-num' can be any non-negative number.
        real(dp), intent(in) :: date_num !! Date number in seconds since ref_year
        integer, intent(out) :: date_array(6) !! Datetime [y,m,d,h,m,s]
        integer, intent(in) :: ref_year !! Reference year

        ! Local variables
        integer :: y, m, d
        real(dp) :: day_fraction
        integer :: day_int
        integer :: day_count, days_in_year
        integer :: rest_seconds
        integer :: daysinmonth(12) = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

        ! Check that date number if positive
        if (date_num < 0) then
            print "(a)", "ERROR! In number_to_date, 'date_num' must be >= 0.0"
            stop 1
        end if

        ! Set day fraction to the nearest second. Avoiding round off errors
        day_int = idint(date_num)
        day_fraction = date_num - day_int

        ! Determine hours, minutes and seconds
        date_array = 0
        rest_seconds = int(day_fraction*24.0*3600.0 + 0.5) ! Rounded off
        date_array(4) = int(rest_seconds/3600.0)
        date_array(5) = int((rest_seconds/60.0 - date_array(4)*60.0))
        date_array(6) = int((rest_seconds - date_array(4)*3600.0 - date_array(5)*60.0))

        ! Count up days keeping track of the year, month and day of month
        
        ! Determine year
        y = ref_year - 1
        day_count = 0
        do while (day_count .le. day_int)
            y = y + 1
            days_in_year = 365
            if (((mod(y,4) .eq. 0) .and. (mod(y,100) .ne. 0)) .or. (mod(y,400) .eq. 0)) days_in_year = 366 ! Leap year
            day_count = day_count + days_in_year
        end do
        date_array(1) = y
        day_count = day_count - days_in_year

        ! Determine month given the year
        daysinmonth(2) = 28
        if (((mod(date_array(1),4) .eq. 0) .and. (mod(date_array(1),100) .ne. 0)) .or. (mod(date_array(1),400) .eq. 0)) daysinmonth(2) = 29 ! Leap month
        m = 0

        do while (day_count .le. day_int)
            m = m + 1
            day_count = day_count + daysinmonth(m)
        end do
        date_array(2) = m
        day_count = day_count - daysinmonth(m)

        ! Determine day
        d = 0
        do while (day_count .le. day_int)
            d = d + 1
            day_count = day_count + 1
        end do
        date_array(3) = d
    end subroutine number_to_date

    function date_to_number(date_array, ref_year) result(res)
        !! Returns a date number from an array with date and time
        !!
        !! Note: Some wrong dates (e.g., 2023-02-29) will still return a number which will be wrong!
        integer, intent(in) :: date_array(6) !! Datetime [y,m,d,h,m,s]
        integer, intent(in) :: ref_year !! Reference year
        real(dp) :: res !! Date number in seconds since ref_year

        ! Local variables
        integer :: y, m
        integer :: daysinmonth(12) = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

        res = 0.0
        if (date_array(1) .gt. ref_year) then
            ! Add up days in the year
            do y = ref_year, date_array(1) - 1
                if (((mod(y,4) .eq. 0) .and. (mod(y,100) .ne. 0)) .or. (mod(y,400) .eq. 0)) then
                    daysinmonth(2) = 29
                else
                    daysinmonth(2) = 28
                end if
                do m = 1, 12
                    res = res + sngl(daysinmonth(m))
                end do
            end do
        end if

        ! Add up days in the remaining months
        if (((mod(date_array(1),4) .eq. 0) .and. (mod(date_array(1),100) .ne. 0)) .or. (mod(date_array(1),400) .eq. 0)) then
            daysinmonth(2) = 29
        else
            daysinmonth(2) = 28
        end if
        if (date_array(2) .gt. 1) then
            do m = 1, date_array(2) - 1
                res = res + sngl(daysinmonth(m))
            end do
        end if

        res = res + sngl(date_array(3)) - 1.0
        res = res + sngl(date_array(4))/24.0 ! starts at 0
        res = res + sngl(date_array(5))/24.0/60.0 ! starts at 0
        res = res + sngl(date_array(6))/24.0/60.0/60.0 ! starts at 0
    end function date_to_number

    function date_to_julian(date_array, ref_year) result(res)
        !! Returns Julian day from an array with date and time
        integer, intent(in) :: date_array(6) !! Datetime [y,m,d,h,m,s]
        integer, intent(in) :: ref_year !! Reference year
        integer :: res !! Julian day

        ! Local variables
        integer :: b(6)

        b(1) = date_array(1)
        b(2) = 1
        b(3) = 1
        b(4) = 0
        b(5) = 0
        b(6) = 0

        res = int(date_to_number(date_array, ref_year) - date_to_number(b, ref_year) + 1)
    end function date_to_julian

    subroutine datestr_to_date(date_str, format_str, date_array)
        !! Converts date string to a date array based on the format string
        !! 
        !! Example:
        !!     call datestr_to_date("20190429123005", "yyyymmdd", date_array_out)
        !! Returns:
        !!     date_arrat_out = [2019, 4, 29, 0, 0, 0]
        character(*), intent(in) :: date_str !! Date string [yyyymmddHHMMSS]
        character(*), intent(in) :: format_str !! Format string
        integer, intent(out) :: date_array(6) !! Datetime [y,m,d,h,m,s]
        
        ! Local variables
        integer :: pos

        ! Extract elements of the date string
        pos = index(format_str, 'yyyy')
        if (pos .gt. 0) then
            read(date_str(pos:pos+3),*) date_array(1)
        else
            date_array(1) = 0
        end if
        pos = index(format_str, 'mm')
        if (pos .gt. 0) then
            read(date_str(pos:pos+1),*) date_array(2)
        else
            date_array(2) = 0
        end if
        pos = index(format_str, 'dd')
        if (pos .gt. 0) then
            read(date_str(pos:pos+1),*) date_array(3)
        else
            date_array(3) = 0
        end if
        pos = index(format_str, 'HH')
        if (pos .gt. 0) then
            read(date_str(pos:pos+1),*) date_array(4)
        else
            date_array(4) = 0
        end if
        pos = index(format_str, 'MM')
        if (pos .gt. 0) then
            read(date_str(pos:pos+1),*) date_array(5)
        else
            date_array(5) = 0
        end if
        pos = index(format_str, 'SS')
        if (pos .gt. 0) then
            read(date_str(pos:pos+1),*) date_array(6)
        else
            date_array(6) = 0
        end if
    end subroutine datestr_to_date

    subroutine date_to_datestr(date_array, format_str, date_str)
        !! Returns a date string based on (yyyy.mm.dd HH:MM:SS) from a date array
        !!
        !! Currently only accepts 'yyyymmddHH', 'yyymmdd' and 'yyyymm' as the format string.
        !!
        !! Example:
        !!     call date_to_datestr([2019, 4, 29, 30, 0, 0], "yyyymmdd", date_str)
        !! Returns:
        !!     date_str = "20190429"
        integer, intent(in) :: date_array(6) !! Datetime [y,m,d,h,m,s]
        character(len=*), intent(in) :: format_str !! Format string
        character(len=*), intent(out) :: date_str !! Date string

        ! Local variables
        integer :: pos

        date_str = format_str

        ! To avoid just putting in date parts e.g. mm or dd that might occurr in a string then it is required that at least two of the date
        ! strings are present, i.e. yyyy, mm and dd or HH, MM and SS
        if ((index(format_str, 'yyyy') .gt. 0 .and. index(format_str, 'mm') .gt. 0) .or. (index(format_str, 'yyyy') .gt. 0 &
            .and. index(format_str, 'dd') .gt. 0) .or. (index(format_str, 'dd') .gt. 0 .and.index(format_str, 'mm').gt.0) &
            .or. (index(format_str, 'HH') .gt. 0 .and. index(format_str,'MM') .gt. 0) .or. (index(format_str,'HH') .gt. 0 &
            .and. index(format_str,'SS') .gt. 0) .or. (index(format_str,'MM') .gt. 0 .and. index(format_str,'SS') .gt. 0)) then
            
            ! Do nothing but continue with routine as this is a valid format for date string substitution
        else
            ! Leave the routine
            return
        end if

        ! Now it only accepts the strings 'yyyymmddHH', 'yyyymmdd' and 'yyyymm' for replacement
        pos = index(format_str, 'yyyymmddHH')
        if (pos .gt. 0) then
            write(date_str(pos:pos+3),'(i4)') date_array(1)
            if (date_array(2) .gt. 9) then
                write(date_str(pos+4:pos+5),'(i2)') date_array(2)
            else
                write(date_str(pos+4:pos+5),'(a1,i1)') '0', date_array(2)
            end if
            if (date_array(3) .gt. 9) then
                write(date_str(pos+6:pos+7),'(i2)') date_array(3)
            else
                write(date_str(pos+6:pos+7),'(a1,i1)') '0', date_array(3)
            end if
            if (date_array(4) .gt. 9) then
                write(date_str(pos+8:pos+9),'(i2)') date_array(4)
            else
                write(date_str(pos+8:pos+9),'(a1,i1)') '0', date_array(4)
            end if
            return
        end if

        pos = index(format_str, 'yyyymmdd')
        if (pos .gt. 0) then
            write(date_str(pos:pos+3),'(i4)') date_array(1)
            if (date_array(2) .gt. 9) then
                write(date_str(pos+4:pos+5),'(i2)') date_array(2)
            else
                write(date_str(pos+4:pos+5),'(a1,i1)') '0', date_array(2)
            end if
            if (date_array(3) .gt. 9) then
                write(date_str(pos+6:pos+7),'(i2)') date_array(3)
            else
                write(date_str(pos+6:pos+7),'(a1,i1)') '0', date_array(3)
            end if
            return
        end if

        pos = index(format_str, 'yyyymm')
        if (pos .gt. 0) then
            write(date_str(pos:pos+3),'(i4)') date_array(1)
            if (date_array(2) .gt. 9) then
                write(date_str(pos+4:pos+5),'(i2)') date_array(2)
            else
                write(date_str(pos+4:pos+5),'(a1,i1)') '0', date_array(2)
            end if
            return
        end if
    end subroutine date_to_datestr

    subroutine date_to_datestr_new(date_array, in_format_str, out_a_str)
        integer, intent(in) :: date_array(6) !! Datetime [y,m,d,h,m,s]
        character(*), intent(in) :: in_format_str !! Format string
        character(*), intent(out) ::  out_a_str !! Date string

        ! Local variables
        integer :: pos
        
        ! Extract elements of the date string
        pos = index(in_format_str, 'yyyy')
        if (pos .gt. 0) then
            write(out_a_str(pos:pos+3),'(i4)') date_array(1)
        end if

        pos = index(in_format_str, 'mm')
        if (pos .gt. 0) then
            if (date_array(2) .gt. 9) then
                write(out_a_str(pos:pos+1),'(i2)') date_array(2)
            else
                write(out_a_str(pos:pos+1),'(a1,i1)') '0', date_array(2)
            end if
        end if

        pos = index(in_format_str, 'dd')
        if (pos.gt.0) then
            if (date_array(3) .gt. 9) then
                write(out_a_str(pos:pos+1),'(i2)') date_array(3)
            else
                write(out_a_str(pos:pos+1),'(a1,i1)') '0', date_array(3)
            end if
        end if

        pos = index(in_format_str, 'HH')
        if (pos .gt. 0) then
            if (date_array(4) .gt. 9) then
                write(out_a_str(pos:pos+1),'(i2)') date_array(4)
            else
                write(out_a_str(pos:pos+1),'(a1,i1)') '0', date_array(4)
            end if
        end if

        pos = index(in_format_str,'MM')
        if (pos .gt. 0) then
            if (date_array(5) .gt. 9) then
                write(out_a_str(pos:pos+1),'(i2)') date_array(5)
            else
                write(out_a_str(pos:pos+1),'(a1,i1)') '0', date_array(5)
            end if
        end if

        pos = index(in_format_str,'SS')
        if (pos .gt. 0) then
            if (date_array(6) .gt. 9) then
                write(out_a_str(pos:pos+1),'(i2)') date_array(6)
            else
                write(out_a_str(pos:pos+1),'(a1,i1)') '0', date_array(6)
            end if
        end if
    end subroutine date_to_datestr_new

    subroutine date_to_datestr_bracket(date_array, in_format_str, out_a_str)
        !! Converts a date array to a date string with brackets in format string
        !!
        !! format string brackets are '<' and '>'.
        integer, intent(in) :: date_array(6) !! Datetime [y,m,d,h,m,s]
        character(*), intent(in) :: in_format_str !! Format string
        character(*), intent(out) ::  out_a_str !! Date string

        ! Local variables
        character(256) :: format_str,a_str
        integer :: pos1, pos2

        ! Only changes dates when they are inside '<.....>'
        ! Removes these once changed
        pos1 = index(in_format_str, '<')
        pos2 = index(in_format_str, '>')

        if (pos1 .le. 0 .or. pos2 .le. 0 .or. pos1+1 .gt. pos2-1) then
            out_a_str = in_format_str
            return
        end if

        ! Reassign format_str to be just the text between <..>
        format_str = in_format_str(pos1+1:pos2-1)
        a_str = format_str

        ! Get the date string
        call date_to_datestr_new(date_array, format_str, a_str)

        ! Insert the a_str into out_a_str, removing the '<>' text
        if (len_trim(in_format_str) .gt. pos2) then
            out_a_str = in_format_str(1:pos1-1)//trim(a_str)//in_format_str(pos2+1:)
        else
            out_a_str = in_format_str(1:pos1-1)//trim(a_str)
        end if
    end subroutine date_to_datestr_bracket

    subroutine date_to_datestr_squarebracket(date_array, in_format_str, out_a_str)
        !! Converts a date array to a date string with brackets in format string
        !!
        !! format string brackets are '<' and '>'.
        integer, intent(in) :: date_array(6)
        character(*), intent(in) :: in_format_str
        character(*), intent(out) :: out_a_str

        ! Local variables
        character(256) :: format_str, a_str
        integer :: pos1, pos2

        ! Only changes dates when they are inside '[.....]'
        ! Removes these once changed
        pos1 = index(in_format_str, '[')
        pos2 = index(in_format_str, ']')

        if (pos1 .le. 0 .or. pos2 .le. 0 .or. pos1+1 .gt. pos2-1) then
            out_a_str = in_format_str
            return
        end if

        ! Reassign format_str to be just the text between <..>
        format_str = in_format_str(pos1+1:pos2-1)
        a_str = format_str

        ! Get date string
        call date_to_datestr_new(date_array, format_str, a_str)

        ! Insert the a_str into out_a_str, removing the '[]' text
        if (len_trim(in_format_str) .gt. pos2) then
            out_a_str = in_format_str(1:pos1-1)//trim(a_str)//in_format_str(pos2+1:)
        else
            out_a_str = in_format_str(1:pos1-1)//trim(a_str)
        end if
    end subroutine date_to_datestr_squarebracket

    function day_of_week (date_array) result(res)
        !! The subroutine calculates the day of week given current datetime,
        !! where DAYW = 1 corresponds to Monday and DAYW = 7 to Sunday. The
        !! algorithm is based on the tables in "Hvem Hva Hvor 1971" (p. 121)
        !! and is valid for all years from 1800 to infinity
        !!
        !! Adapted from EPISODE code
        integer, intent(in) :: date_array(6) !! Datetime [y,m,d,h,m,s]
        integer :: res !! Day of week [1-7]

        ! Local variables
        integer :: jm(12) = [1, 5, 5, 2, 7, 4, 2, 6, 3, 1, 5, 3] ! Column number for each month
        integer :: ir ! Row in HHH table for day of month
        integer :: jc ! Column in HHH table for month
        integer :: nt ! Number in HHH table for row ir and column jc
        integer :: jk ! Column in HHH table for year
        integer :: j4, j100, j400 ! Adjustment values for leap year
        logical :: leap ! If leap year then true else false
        integer :: daym, mnth, year

        ! Extract values from array
        daym = date_array(3)
        mnth = date_array(2)
        year = date_array(1)

        ! Calculate leap year or not
        leap = .false.
        if (mod(year,  4) .eq. 0 .and. .not. (mod(year, 100) .eq. 0 .and. mod(year, 400) .ne. 0)) leap = .true.

        ! Calculate row number for day of month
        ir = mod(daym - 1, 7) + 1

        ! Calculate column number for month
        jc = jm(mnth)
        if (leap .and. (mnth .eq. 1 .or. mnth .eq. 2)) jc = jc + 1

        ! Calculate "number" in HHH table with row IR and column JC
        nt = mod(ir + 7 - jc, 7) + 1

        ! Calculate column number for year (adjusting for leap years)
        j4 = (year - 1800)/4
        j100 = (year - 1800)/100
        j400 = (year - 1600)/400
        jk = mod(year - 1800 + j4 - j100 + j400 + 3 - 1, 7) + 1

        ! Calculate day of week
        res = mod(jk - 1 + nt - 1, 7) + 1
    end function day_of_week

    function summer_time_europe(date_array) result(summer_time)
        !! Returns true if supplied date is during summer time in Europe, else false
        integer, intent(in) :: date_array(6) !! Datetime [y,m,d,h,m,s]
        logical :: summer_time !! True if summer time, else false

        ! Local variables
        integer :: a(6), b_start(6), b_end(6)
        integer :: ref_year, year
        real(dp) :: datenum_start, datenum_end, datenum

        a = date_array
        ref_year = 2000
        b_start = 0
        b_end = 0
        year = a(1)
        b_start(1) = a(1)
        b_start(2) = 3
        b_start(3) = (31 - mod((((5 * year)/4) + 4), 7))
        b_start(4) = 1
        b_end(1) = a(1)
        b_end(2) = 10
        b_end(3) = (31 - mod((((5 * year)/4) + 1), 7))
        b_end(4) = 1

        datenum_start = date_to_number(b_start, ref_year)
        datenum_end = date_to_number(b_end, ref_year)
        datenum = date_to_number(a, ref_year)

        summer_time = .false.
        if (datenum .ge. datenum_start .and. datenum .lt. datenum_end) summer_time = .true.
    end function summer_time_europe
    
    subroutine get_sun_angles(lat, lon, date_a, date_num, difutc_h, azimuth_ang, zenith_ang)
        !! Returns the azimuth and zenith angles
        !!
        !! Routine modified from NORTRIP
        real, intent(in) :: lat
        real, intent(in) :: lon
        integer, intent(in) :: date_a(6)
        real(dp), intent(in) :: date_num
        real, intent(in) :: difutc_h
        real, intent(out) :: azimuth_ang
        real, intent(out) :: zenith_ang

        ! Local variables
        real :: julian_day, time_s, dayang, dec, eqtime, solartime, hourang, azt, az
        real, parameter :: s0 = 1367.0
        real, parameter :: pi = 3.14159/180.0
        integer :: ref_year = 2000

        if (date_a(1) .eq. 0) then
            julian_day = real(date_num)
        else
            julian_day = date_to_julian(date_a, ref_year)
        end if

        time_s = (julian_day - 1)*24.0*3600.0
        dayang = 360.0/365.0*(julian_day - 1.0)
        dec = 0.396 - 22.91*cos(pi*dayang) + 4.025*sin(pi*dayang)
        eqtime = (1.03 + 25.7*cos(pi*dayang) - 440.0*sin(pi*dayang) - 201.0*cos(2.0*pi*dayang) - 562.0*sin(2.0*pi*dayang))/secphour
        solartime = mod(time_s + secpday + secphour*(lon/15.0 + difutc_h + eqtime), secpday)
        hourang = 15.0*(12.0 - solartime/secphour)

        ! Set zenith angle for atmospheric corrections
        azt = sin(pi*dec)*sin(pi*lat) + cos(pi*dec)*cos(pi*lat)*cos(pi*hourang)
        if (abs(azt) .lt. 1.0) then
            az = acos(azt)/pi
        else
            az = 0.0
        end if

        azimuth_ang = 180.0 - hourang
        zenith_ang = az
    end subroutine get_sun_angles

end module time_functions

