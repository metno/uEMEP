module time_functions

    implicit none
    private

    public :: global_radiation_sub, number_to_date, date_to_number, &
        date_to_datestr_bracket, date_to_datestr_squarebracket, date_to_datestr, &
        datestr_to_date, day_of_week, summer_time_europe

contains
!----------------------------------------------------------------------
! Various functions for manipulating time
!----------------------------------------------------------------------

!----------------------------------------------------------------------
    subroutine number_to_date(date_num,date_array,ref_year)
    
    implicit none
    
    double precision date_num
    integer ref_year
    integer y,m,d,i
    integer date_array(6)
    double precision day_fraction
    integer day_int
    integer day_count,days_in_year
    integer rest_seconds
    integer daysinmonth(12)
    data (daysinmonth(i),i=1,12) /31,28,31,30,31,30,31,31,30,31,30,31/ 
    
    !ref_year=1900
    !Set day fraction to the nearest second. Avoiding round off errors
    day_int=idint(date_num)
    day_fraction=(date_num-day_int)
      
    !Determine hours, minutes and seconds
    date_array=0
    rest_seconds=int(day_fraction*24.*3600.+.5) !Rounded off
    date_array(4)=int(rest_seconds/3600.)
    date_array(5)=int((rest_seconds/60.-date_array(4)*60.))
    date_array(6)=int((rest_seconds-date_array(4)*3600.-date_array(5)*60.))
    
    !Count up days keeping track of the year, month and day of month
    
    !Determine year
    y=ref_year-1
    day_count=0
    do while (day_count.le.day_int)
        y=y+1
        days_in_year=365
        if (((mod(y,4).eq.0).and.(mod(y,100).ne.0)).or.(mod(y,400).eq.0)) days_in_year=366
        day_count=day_count+days_in_year     
    enddo
    date_array(1)=y
    day_count=day_count-days_in_year

    !Determine month given the year
    daysinmonth(2)=28
    if (((mod(date_array(1),4).eq.0).and.(mod(date_array(1),100).ne.0)).or.(mod(date_array(1),400).eq.0)) daysinmonth(2)=29
    m=0
    !day_count=0

    do while (day_count.le.day_int)
        m=m+1
        day_count=day_count+daysinmonth(m)
    enddo
    date_array(2)=m
    day_count=day_count-daysinmonth(m)
    
    !Determine day
    d=0
    do while (day_count.le.day_int)
        d=d+1
        day_count=day_count+1   
    enddo
    date_array(3)=d

    end subroutine number_to_date
!----------------------------------------------------------------------

!----------------------------------------------------------------------
    function date_to_number(a,ref_year)
    
    implicit none
    
    double precision date_to_number
    integer ref_year
    integer y,m,i
    integer a(6)
    
    integer daysinmonth(12)
    data (daysinmonth(i),i=1,12) /31,28,31,30,31,30,31,31,30,31,30,31/ 
    
    !ref_year=1900
    date_to_number=0.
    daysinmonth(2)=28
    if (a(1).gt.ref_year) then
    !Add up days in the year
        do y=ref_year,a(1)-1
            if (((mod(y,4).eq.0).and.(mod(y,100).ne.0)).or.(mod(y,400).eq.0)) then
                daysinmonth(2)=29
            else
                daysinmonth(2)=28
            endif
            do m=1,12            
                date_to_number=date_to_number+sngl(daysinmonth(m))
            end do     
        end do
    endif
    !Add up days in the remaining months
    if (((mod(a(1),4).eq.0).and.(mod(a(1),100).ne.0)).or.(mod(a(1),400).eq.0)) then
        daysinmonth(2)=29
    else
        daysinmonth(2)=28
    endif
    if (a(2).gt.1) then
        do m=1,a(2)-1
            date_to_number=date_to_number+sngl(daysinmonth(m))
        enddo
    endif
    
    date_to_number=date_to_number+sngl(a(3))-1.
    date_to_number=date_to_number+sngl(a(4))/24. !starts at 0
    date_to_number=date_to_number+sngl(a(5))/24./60. !starts at 0  
    date_to_number=date_to_number+sngl(a(6))/24./60./60. !starts at 0  
    !write(*,*) date_to_number

    end function date_to_number
!----------------------------------------------------------------------

!----------------------------------------------------------------------
    function date_to_julian(a,ref_year)
    
    implicit none
    
    real date_to_julian
    integer a(6),b(6)
    integer ref_year
    
    b(1)=a(1)
    b(2)=1
    b(3)=1
    b(4)=0
    b(5)=0
    b(6)=0
    
    date_to_julian=date_to_number(a,ref_year)-date_to_number(b,ref_year)+1
    
    end function date_to_julian
!----------------------------------------------------------------------

!----------------------------------------------------------------------
    subroutine datestr_to_date(a_str,format_str,a)
    
    implicit none
    
    character(256) a_str,format_str
    integer a(6)
    integer pos
    
    !based on (yyyy.mm.dd HH:MM:SS)
    
    !extract year
    pos=index(format_str,'yyyy')
    if (pos.gt.0) then
        read(a_str(pos:pos+3),*) a(1)
    else
        a(1)=0
    endif
    pos=index(format_str,'mm')
    if (pos.gt.0) then
        read(a_str(pos:pos+1),*) a(2)
    else
        a(2)=0
    endif
    pos=index(format_str,'dd')
    if (pos.gt.0) then
        read(a_str(pos:pos+1),*) a(3)
    else
        a(3)=0
    endif
    pos=index(format_str,'HH')
    if (pos.gt.0) then
        read(a_str(pos:pos+1),*) a(4)
    else
        a(4)=0
    endif
    pos=index(format_str,'MM')
    if (pos.gt.0) then
        read(a_str(pos:pos+1),*) a(5)
    else
        a(5)=0
    endif
    pos=index(format_str,'SS')
    if (pos.gt.0) then
        read(a_str(pos:pos+1),*) a(6)
    else
        a(6)=0
    endif
    
    end subroutine datestr_to_date
!----------------------------------------------------------------------

!----------------------------------------------------------------------
    subroutine date_to_datestr(a,format_str,a_str)
    
    implicit none
    
    character(*) a_str,format_str
    integer a(6)
    integer pos
    
    !based on (yyyy.mm.dd HH:MM:SS)
    
    a_str=format_str
    
    !To avoid just putting in date parts e.g. mm or dd that might occurr in a string then it is required that at least two of the date
    !strings are present, i.e. yyyy, mm and dd or HH, MM and SS
    
    if ((index(format_str,'yyyy').gt.0.and.index(format_str,'mm').gt.0).or.(index(format_str,'yyyy').gt.0.and.index(format_str,'dd').gt.0).or.(index(format_str,'dd').gt.0.and.index(format_str,'mm').gt.0).or. &
        (index(format_str,'HH').gt.0.and.index(format_str,'MM').gt.0).or.(index(format_str,'HH').gt.0.and.index(format_str,'SS').gt.0).or.(index(format_str,'MM').gt.0.and.index(format_str,'SS').gt.0)) then
        !Do nothing but continue with routine as this is a valid format for date string substitution
    else
        !Leave the routine
        return
    endif
    

    !Now it only accepts the two strings 'yyyymmdd' and 'yyyymmddHH' for replacement
    
    pos=index(format_str,'yyyymmddHH')
    if (pos.gt.0) then
        write(a_str(pos:pos+3),'(i4)') a(1)
        if (a(2).gt.9) then
            write(a_str(pos+4:pos+5),'(i2)') a(2)
        else
            write(a_str(pos+4:pos+5),'(a1,i1)') '0',a(2)
        endif     
        if (a(3).gt.9) then
            write(a_str(pos+6:pos+7),'(i2)') a(3)
        else
            write(a_str(pos+6:pos+7),'(a1,i1)') '0',a(3)
        endif     
        if (a(4).gt.9) then
            write(a_str(pos+8:pos+9),'(i2)') a(4)
        else
            write(a_str(pos+8:pos+9),'(a1,i1)') '0',a(4)
        endif     
        return
    else
        !a_str(pos:pos+3)='0000'
    endif

    pos=index(format_str,'yyyymmdd')
    if (pos.gt.0) then
        write(a_str(pos:pos+3),'(i4)') a(1)
        if (a(2).gt.9) then
            write(a_str(pos+4:pos+5),'(i2)') a(2)
        else
            write(a_str(pos+4:pos+5),'(a1,i1)') '0',a(2)
        endif     
        if (a(3).gt.9) then
            write(a_str(pos+6:pos+7),'(i2)') a(3)
        else
            write(a_str(pos+6:pos+7),'(a1,i1)') '0',a(3)
        endif     
        return
    else
        !a_str(pos:pos+3)='0000'
    endif
    
    pos=index(format_str,'yyyymm')
    if (pos.gt.0) then
        write(a_str(pos:pos+3),'(i4)') a(1)
        if (a(2).gt.9) then
            write(a_str(pos+4:pos+5),'(i2)') a(2)
        else
            write(a_str(pos+4:pos+5),'(a1,i1)') '0',a(2)
        endif         
        return
    else
        !a_str(pos:pos+3)='0000'
    endif
    
    !return
    !Do not do the rest
    
    !extract year
    pos=index(format_str,'yyyy')
    if (pos.gt.0) then
        write(a_str(pos:pos+3),'(i4)') a(1)
    else
        !a_str(pos:pos+3)='0000'
    endif
    
    pos=index(format_str,'mm')
    if (pos.gt.0) then
        if (a(2).gt.9) then
            write(a_str(pos:pos+1),'(i2)') a(2)
        else
            write(a_str(pos:pos+1),'(a1,i1)') '0',a(2)
        endif     
    else
        !a_str(pos:pos+1)='00'
    endif

    pos=index(format_str,'dd')
    if (pos.gt.0) then
        if (a(3).gt.9) then
            write(a_str(pos:pos+1),'(i2)') a(3)
        else
            write(a_str(pos:pos+1),'(a1,i1)') '0',a(3)
        endif     
    else
        !a_str(pos:pos+1)='00'
    endif
    
    pos=index(format_str,'HH')
    if (pos.gt.0) then
        if (a(4).gt.9) then
            write(a_str(pos:pos+1),'(i2)') a(4)
        else
            write(a_str(pos:pos+1),'(a1,i1)') '0',a(4)
        endif     
    else
        !a_str(pos:pos+1)='00'
    endif

    pos=index(format_str,'MM')
    if (pos.gt.0) then
        if (a(5).gt.9) then
            write(a_str(pos:pos+1),'(i2)') a(5)
        else
            write(a_str(pos:pos+1),'(a1,i1)') '0',a(5)
        endif     
    else
    !    !a_str(pos:pos+1)='00'
    endif
    
    pos=index(format_str,'SS')
    if (pos.gt.0) then
        if (a(6).gt.9) then
            write(a_str(pos:pos+1),'(i2)') a(6)
        else
            write(a_str(pos:pos+1),'(a1,i1)') '0',a(6)
        endif     
    else
        !a_str(pos:pos+1)='00'
    endif
    
    end subroutine date_to_datestr
!----------------------------------------------------------------------

!----------------------------------------------------------------------
    subroutine date_to_datestr_bracket(a,in_format_str,out_a_str)
    
    implicit none
    
    character(*), intent(out) ::  out_a_str
    character(*), intent(in) :: in_format_str
    integer, intent(in) :: a(6)
    character(256) format_str,a_str
    integer pos
    integer pos1,pos2
    
    !based on (yyyy.mm.dd HH:MM:SS)
    
    !a_str=format_str
    
    !To avoid just putting in date parts e.g. mm or dd that might occurr in a string then it is required that at least two of the date
    !strings are present, i.e. yyyy, mm and dd or HH, MM and SS
    
    !Only changes dates when they are inside '<.....>'
    !Removes these once changed
    pos1=index(in_format_str,'<')
    pos2=index(in_format_str,'>')
    
    if (pos1.le.0.or.pos2.le.0.or.pos1+1.gt.pos2-1) then
        out_a_str=in_format_str
        return
    endif
    
    !Reassign format_str to be just the text between <..>
    format_str=in_format_str(pos1+1:pos2-1)
    a_str=format_str
    
    !extract year
    pos=index(format_str,'yyyy')
    if (pos.gt.0) then
        write(a_str(pos:pos+3),'(i4)') a(1)
    endif
    
    pos=index(format_str,'mm')
    if (pos.gt.0) then
        if (a(2).gt.9) then
            write(a_str(pos:pos+1),'(i2)') a(2)
        else
            write(a_str(pos:pos+1),'(a1,i1)') '0',a(2)
        endif     
    endif

    pos=index(format_str,'dd')
    if (pos.gt.0) then
        if (a(3).gt.9) then
            write(a_str(pos:pos+1),'(i2)') a(3)
        else
            write(a_str(pos:pos+1),'(a1,i1)') '0',a(3)
        endif     
    endif
    
    pos=index(format_str,'HH')
    if (pos.gt.0) then
        if (a(4).gt.9) then
            write(a_str(pos:pos+1),'(i2)') a(4)
        else
            write(a_str(pos:pos+1),'(a1,i1)') '0',a(4)
        endif     
    endif

    pos=index(format_str,'MM')
    if (pos.gt.0) then
        if (a(5).gt.9) then
            write(a_str(pos:pos+1),'(i2)') a(5)
        else
            write(a_str(pos:pos+1),'(a1,i1)') '0',a(5)
        endif     
    endif
    
    pos=index(format_str,'SS')
    if (pos.gt.0) then
        if (a(6).gt.9) then
            write(a_str(pos:pos+1),'(i2)') a(6)
        else
            write(a_str(pos:pos+1),'(a1,i1)') '0',a(6)
        endif     
    endif
    
    !insert the a_str into out_a_str, removing the '<>' text
    if (len_trim(in_format_str).gt.pos2) then
        out_a_str=in_format_str(1:pos1-1)//trim(a_str)//in_format_str(pos2+1:)
    else
        out_a_str=in_format_str(1:pos1-1)//trim(a_str)
    endif
    
    !write(*,*) trim(in_format_str),trim(out_a_str)
    !stop
    
    end subroutine date_to_datestr_bracket
!----------------------------------------------------------------------

!----------------------------------------------------------------------
    subroutine date_to_datestr_squarebracket(a,in_format_str,out_a_str)
    
    implicit none
    
    character(*), intent(out) ::  out_a_str
    character(*), intent(in) :: in_format_str
    integer, intent(in) :: a(6)
    character(256) format_str,a_str
    integer pos
    integer pos1,pos2
    
    !based on (yyyy.mm.dd HH:MM:SS)
    
    !a_str=format_str
    
    !To avoid just putting in date parts e.g. mm or dd that might occurr in a string then it is required that at least two of the date
    !strings are present, i.e. yyyy, mm and dd or HH, MM and SS
    
    !Only changes dates when they are inside '[.....]'
    !Removes these once changed
    pos1=index(in_format_str,'[')
    pos2=index(in_format_str,']')
    
    if (pos1.le.0.or.pos2.le.0.or.pos1+1.gt.pos2-1) then
        out_a_str=in_format_str
        return
    endif
    
    !Reassign format_str to be just the text between <..>
    format_str=in_format_str(pos1+1:pos2-1)
    a_str=format_str
    
    !extract year
    pos=index(format_str,'yyyy')
    if (pos.gt.0) then
        write(a_str(pos:pos+3),'(i4)') a(1)
    endif
    
    pos=index(format_str,'mm')
    if (pos.gt.0) then
        if (a(2).gt.9) then
            write(a_str(pos:pos+1),'(i2)') a(2)
        else
            write(a_str(pos:pos+1),'(a1,i1)') '0',a(2)
        endif     
    endif

    pos=index(format_str,'dd')
    if (pos.gt.0) then
        if (a(3).gt.9) then
            write(a_str(pos:pos+1),'(i2)') a(3)
        else
            write(a_str(pos:pos+1),'(a1,i1)') '0',a(3)
        endif     
    endif
    
    pos=index(format_str,'HH')
    if (pos.gt.0) then
        if (a(4).gt.9) then
            write(a_str(pos:pos+1),'(i2)') a(4)
        else
            write(a_str(pos:pos+1),'(a1,i1)') '0',a(4)
        endif     
    endif

    pos=index(format_str,'MM')
    if (pos.gt.0) then
        if (a(5).gt.9) then
            write(a_str(pos:pos+1),'(i2)') a(5)
        else
            write(a_str(pos:pos+1),'(a1,i1)') '0',a(5)
        endif     
    endif
    
    pos=index(format_str,'SS')
    if (pos.gt.0) then
        if (a(6).gt.9) then
            write(a_str(pos:pos+1),'(i2)') a(6)
        else
            write(a_str(pos:pos+1),'(a1,i1)') '0',a(6)
        endif     
    endif
    
    !insert the a_str into out_a_str, removing the '[]' text
    if (len_trim(in_format_str).gt.pos2) then
        out_a_str=in_format_str(1:pos1-1)//trim(a_str)//in_format_str(pos2+1:)
    else
        out_a_str=in_format_str(1:pos1-1)//trim(a_str)
    endif
    
    !write(*,*) trim(in_format_str),trim(out_a_str)
    !stop
    
    end subroutine date_to_datestr_squarebracket
!----------------------------------------------------------------------

!----------------------------------------------------------------------
    function day_of_week (a)
    !Adapted from EPISODE code
    
    implicit none
    
!The subroutine calculates the day of week given current datetime,
!where DAYW = 1 corresponds to Monday and DAYW = 7 to Sunday. The
!algorithm is based on the tables in "Hvem Hva Hvor 1971" (p. 121)
!and is valid for all years from 1800 to infinity!

      !USE mod_time

! Local variables

      INTEGER JM(12)
      INTEGER IR
      INTEGER JC
      INTEGER NT
      INTEGER JK
      INTEGER J4
      INTEGER J100
      INTEGER J400
      LOGICAL LEAP

! JM   - Column number for each month
! IR   - Row    in HHH table for day of month
! JC   - Column in HHH table for month
! NT   - Number in HHH table for row IR and column JC
! JK   - Column in HHH table for year
! J4   - Adjustment value for leap year
! J100 - Adjustment value for leap year
! J400 - Adjustment value for leap year
! LEAP - If leap year then true else false
      
      integer DAYM,MNTH,YEAR
      integer day_of_week
      integer a(6)
      
      DAYM=a(3)
      MNTH=a(2)
      YEAR=a(1)

!Calculate leap year or not
      LEAP = .FALSE.
      IF (MOD(YEAR,  4) .EQ. 0 .AND. .NOT. (MOD(YEAR,100) .EQ. 0 .AND. MOD(YEAR,400) .NE. 0)) LEAP = .TRUE. 

      ! Set data in table JM

      DATA JM/1,5,5,2,7,4,2,6,3,1,5,3/

! Calculate row    number for day of month

      IR = MOD(DAYM - 1,7) + 1

! Calculate column number for month

      JC = JM(MNTH)
      IF (LEAP .AND. (MNTH .EQ. 1 .OR. MNTH .EQ. 2)) JC = JC + 1

! Calculate "number" in HHH table with row IR and column JC

      NT = MOD(IR + 7 - JC,7) + 1

! Calculate column number for year (adjusting for leap years)

      J4   = (YEAR - 1800)/4
      J100 = (YEAR - 1800)/100
      J400 = (YEAR - 1600)/400

      JK = MOD(YEAR - 1800 + J4 - J100 + J400 + 3 - 1,7) + 1

! Calculate day of week

      day_of_week = MOD(JK - 1 + NT - 1,7) + 1

      RETURN

! End of subroutine CDAYW

    end function day_of_week
!----------------------------------------------------------------------

!----------------------------------------------------------------------
    function summer_time_europe(a_in)
    
    implicit none
    
    logical summer_time_europe
    integer a(6),a_in(6)
    integer b_start(6),b_end(6)
    integer ref_year
    integer year
    double precision datenum_start,datenum_end,datenum

    a=a_in
    ref_year=2000
    b_start=0
    b_end=0
    year=a(1)
    b_start(1)=a(1)
    b_start(2)=3
    b_start(3)=(31 - mod((((5 * year)/4) + 4),7))
    b_start(4)=1
    b_end(1)=a(1)
    b_end(2)=10
    b_end(3)=(31 - mod((((5 * year)/4) + 1),7))
    b_end(4)=1
    
    datenum_start=date_to_number(b_start,ref_year)
    datenum_end=date_to_number(b_end,ref_year)
    datenum=date_to_number(a,ref_year)
    
    summer_time_europe=.false.
    if (datenum.ge.datenum_start.and.datenum.lt.datenum_end) summer_time_europe=.true.
    
    !write(*,*) b_start(3),b_end(3),summer_time_europe
   
    end function summer_time_europe
!----------------------------------------------------------------------
    
!==========================================================================
    !Routine taken from NORTRIP
    subroutine global_radiation_sub(LAT,LON,date_a,date_num,DIFUTC_H,Z_SURF,N_CLOUD,ALBEDO,SOLAR_NET,azimuth_ang,zenith_ang)
    !RETURNS THE NET SHORT WAVE RADIATION
    !DETERMINES SHORT WAVE FLUXES ON A HORIZONTAL SURFACE
    
    implicit none
    !INPUT
    real, intent(in) :: LAT,LON,DIFUTC_H,Z_SURF,N_CLOUD,ALBEDO
    integer, intent(in) ::  date_a(6)
    double precision, intent(in) ::  date_num
    !OUTPUT
    real, intent(out) ::  SOLAR_NET,azimuth_ang,zenith_ang
    !INTERNAL
    real JULIAN_DAY,TIME_S,DAYANG,DEC,EQTIME,SOLARTIME,HOURANG,AZT,AZ
    real TAU_A,TAU_C,DAY_BIG,DAY_END,SOLAR_IN
    real SECPHOUR,SECPDAY,PI,S0
    parameter (SECPHOUR=3600.,SECPDAY=86400.,PI=3.14159/180.,S0=1367.)
    integer :: ref_year=2000
    logical :: calculate_solar=.false.
    
    !FUNCTIONS
    !double precision date_to_number
    
    if (date_a(1).eq.0) then
        JULIAN_DAY=real(date_num)
    else
        JULIAN_DAY=date_to_julian(date_a,ref_year)
    endif
    
    TIME_S=(JULIAN_DAY-1)*24.*3600.
    ![Y, M, D, H, MN, S] = datevec(date_num)
    !JULIAN_DAY=floor(date_num(i)-datenum(Y, 0, 0, 0, 0, 0)+1)
    !TIME_S=(date_num(i)-datenum(Y, M, D, 0, 0, 0))*24*3600

	!DAYANG=1
	!DEC=1
	!EQTIME=1
	!SOLARTIME=1
	!HOURANG=1

    DAYANG=360./365*(JULIAN_DAY-1.)
	DEC=0.396-22.91*cos(PI*DAYANG)+4.025*sin(PI*DAYANG)
	EQTIME=(1.03+25.7*cos(PI*DAYANG)-440.*sin(PI*DAYANG)-201.*cos(2.*PI*DAYANG)-562.*sin(2.*PI*DAYANG))/SECPHOUR
	SOLARTIME=mod(TIME_S+SECPDAY+SECPHOUR*(LON/15.+DIFUTC_H+EQTIME),SECPDAY)
	HOURANG=15.*(12.-SOLARTIME/SECPHOUR)
    
!	SET ZENITH ANGLE FOR ATMOSPHERIC CORRECTIONS
    !AZT=0.5
	AZT=sin(PI*DEC)*sin(PI*LAT)+cos(PI*DEC)*cos(PI*LAT)*cos(PI*HOURANG)
	if (abs(AZT).lt.1.) then
	  AZ=acos(AZT)/PI
    else
	  AZ=0.
    endif
    
    !write(*,*) AZT,AZ

    if (calculate_solar) then
        
!	CORRECTIONS FOR ATMOSPHERE AND CLOUD FROM OERLEMANS (GREENLAND)
    !These need to be updated
    !Have included a correction of 1.1 to match the Stockholm data
    !THe cloud cover transmission is still not assessed
	TAU_A=1.1*(0.75+6.8E-5*Z_SURF-7.1E-9*Z_SURF**2)*(1-.001*AZ)
    TAU_C=1-0.78*N_CLOUD**2*exp(-8.5E-4*Z_SURF)
    

!	SET DAY BEGINNING AND END
	if (abs(tan(PI*DEC)*tan(PI*LAT)).lt.1.) then
	  DAY_BIG=(12.-acos(-tan(PI*DEC)*tan(PI*LAT))/PI/15.)*SECPHOUR
	  DAY_END=(12.+acos(-tan(PI*DEC)*tan(PI*LAT))/PI/15.)*SECPHOUR
    else
	  DAY_BIG=0.
	  DAY_END=24.*SECPHOUR
    endif
    
!	DETERMINE SOLAR RADIATION AT SURFACE DURING DAY
	if ((SOLARTIME.gt.DAY_BIG).and.(SOLARTIME.lt.DAY_END)) then
	  SOLAR_IN=S0*TAU_A*TAU_C*cos(AZ*PI)
    else
	  SOLAR_IN=0.
    endif
    
    SOLAR_NET=SOLAR_IN*(1-ALBEDO)
	!if (SOLARNEW.lt.0.) then
    !    SOLARNEW=0.
    !endif
    
    else
        SOLAR_NET=0
    endif
    
    azimuth_ang=180-HOURANG
    zenith_ang=AZ
    
    end subroutine global_radiation_sub
    
 !==========================================================================    

end module time_functions

