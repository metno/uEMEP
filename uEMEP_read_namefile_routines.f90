!Subroutines and functions for reading in name files
    
    function read_name_real(name_str,default_val,unit_in,unit_out)
    
    use uEMEP_definitions

    implicit none   
    
    integer i,j,k
    real read_name_real
    real default_val
    character(*) name_str
    integer index_val
    integer unit_in,unit_out
    character(256) temp_str,temp_str1,temp_str2
    
    
    read_name_real=default_val
    
    rewind(unit_in)
    
    do while (.not.eof(unit_in))
    
        read(unit_in,'(A)',ERR=30) temp_str
        
        !Remove tabs
        index_val=0
        do i=1,len(temp_str)
            if (ichar(temp_str(i:i)).ne.9) then
                index_val=index_val+1
                temp_str1(index_val:index_val)=temp_str(i:i)
            endif
        enddo
        
        temp_str=ADJUSTL(temp_str1)
        temp_str1=''
        
        !If not a comment
        if (trim(temp_str(1:1)).ne.'!') then
            
            !write(*,*) trim(temp_str)
            
            !Find the position of the equals sign if there is one
            index_val=index(temp_str,'=',back=.false.)
    
            if (index_val.gt.1) then
            
                !Create the string before the equals sign
                temp_str1=trim(temp_str(1:index_val-1))
                
                !Check to see if it is a matching string
                !write(*,*) trim(temp_str1),' ',trim(name_str)
                if (trim(temp_str1).eq.trim(name_str)) then

                    !Create the string after the equals sign
                    temp_str2=temp_str(index_val+1:)
                    temp_str2=adjustl(temp_str2)
                    
                    if (len(trim(temp_str2)).ge.1) then
                        read(temp_str2,*,ERR=20) read_name_real
                        write(unit_out,'(A,es12.4)') 'Setting: '//trim(name_str)//' = ',read_name_real
                    endif
                
                endif
    
            endif
        endif
    
    
20  enddo
    
30  end function read_name_real

    
    
    function read_name_integer(name_str,default_val,unit_in,unit_out)
    
    use uEMEP_definitions

    implicit none   
    
    integer i,j,k
    integer read_name_integer
    integer default_val
    character(*) name_str
    integer index_val
    integer unit_in,unit_out
    character(256) temp_str,temp_str1,temp_str2
    
    
    read_name_integer=default_val
    
    rewind(unit_in)
    
    do while (.not.eof(unit_in))
    
        read(unit_in,'(A)',ERR=30) temp_str
        
        !Remove tabs
        index_val=0
        do i=1,len(temp_str)
            if (ichar(temp_str(i:i)).ne.9) then
                index_val=index_val+1
                temp_str1(index_val:index_val)=temp_str(i:i)
            endif
        enddo
        
        temp_str=ADJUSTL(temp_str1)
        temp_str1=''
        
        !If not a comment
        if (trim(temp_str(1:1)).ne.'!') then
            
            !write(*,*) trim(temp_str)
            
            !Find the position of the equals sign if there is one
            index_val=index(temp_str,'=',back=.false.)
    
            if (index_val.gt.1) then
            
                !Create the string before the equals sign
                temp_str1=trim(temp_str(1:index_val-1))
                
                !Check to see if it is a matching string
                !write(*,*) trim(temp_str1),' ',trim(name_str)
                if (trim(temp_str1).eq.trim(name_str)) then

                    !Create the string after the equals sign
                    temp_str2=temp_str(index_val+1:)
                    temp_str2=adjustl(temp_str2)
                    
                    if (len(trim(temp_str2)).ge.1) then
                        read(temp_str2,*,ERR=20) read_name_integer
                        write(unit_out,'(A,i12)') 'Setting: '//trim(name_str)//' = ',read_name_integer
                    endif
                
                endif
    
            endif
        endif
    
    
20  enddo
    
    
30  end function read_name_integer
    
    
    function read_name_char(name_str,default_val,unit_in,unit_out)
    
    use uEMEP_definitions

    implicit none   
    
    integer i,j,k
    character(256) read_name_char
    character(*) default_val
    character(*) name_str
    integer index_val
    integer unit_in,unit_out
    character(256) temp_str,temp_str1,temp_str2
    character(256) :: call_str='read_name_char'
    
    
    read_name_char=default_val
    
    rewind(unit_in)
    
    do while (.not.eof(unit_in))
    
        read(unit_in,'(A)',ERR=30) temp_str

        !Remove tabs
        index_val=0
        do i=1,len(temp_str)
            if (ichar(temp_str(i:i)).ne.9) then
                index_val=index_val+1
                temp_str1(index_val:index_val)=temp_str(i:i)
            endif
        enddo
        
        temp_str=ADJUSTL(temp_str1)
        temp_str1=''
        
        !If not a comment
        if (trim(temp_str(1:1)).ne.'!') then
            
            !write(*,'(2a)') 'CHECK: ',trim(temp_str)
            
            !Find the position of the equals sign if there is one
            index_val=index(temp_str,'=',back=.false.)
    
            if (index_val.gt.1) then
            
                !Create the string before the equals sign
                temp_str1=trim(temp_str(1:index_val-1))
                
                !Check to see if it is a matching string
                !write(*,*) trim(temp_str1),' ',trim(name_str)
                if (trim(temp_str1).eq.trim(name_str)) then

                    !Create the string after the equals sign
                    temp_str2=temp_str(index_val+1:)
                    temp_str2=adjustl(temp_str2)
    
                    if (len(trim(temp_str2)).ge.1) then
                        !Special for characters so it doesn't read it's own call
                        if (trim(temp_str2(1:min(len(trim(temp_str2)),len(trim(call_str))))).ne.trim(call_str)) then
                            read(temp_str2,*,ERR=20) read_name_char
                            read_name_char=adjustl(read_name_char)
                            write(unit_out,'(A,A)') 'Setting: '//trim(name_str)//' = ',trim(read_name_char)
                        endif
                    endif
                
                endif
    
            endif
        endif
    
    
20  enddo
    
    
30  end function read_name_char
 
    function read_name_logical(name_str,default_val,unit_in,unit_out)
    
    use uEMEP_definitions

    implicit none   
    
    integer i,j,k
    logical read_name_logical
    logical default_val
    character(*) name_str
    integer index_val
    integer unit_in,unit_out
    character(256) temp_str,temp_str1,temp_str2
    
    
    read_name_logical=default_val
    
    rewind(unit_in)
    
    do while (.not.eof(unit_in))
    
        read(unit_in,'(A)',ERR=30) temp_str
        !write(*,*) trim(temp_str)
        
        !Remove tabs
        index_val=0
        do i=1,len(temp_str)
            if (ichar(temp_str(i:i)).ne.9) then
                index_val=index_val+1
                temp_str1(index_val:index_val)=temp_str(i:i)
            endif
        enddo
        
        temp_str=ADJUSTL(temp_str1)
        temp_str1=''

        !If not a comment
        if (trim(temp_str(1:1)).ne.'!') then
            
            !write(*,*) trim(temp_str)
            
            !Find the position of the equals sign if there is one
            index_val=index(temp_str,'=',back=.false.)
    
            if (index_val.gt.1) then
        
                !Create the string before the equals sign
                temp_str1=trim(temp_str(1:index_val-1))
                
                !Check to see if it is a matching string
                !write(*,*) trim(temp_str1),' ',trim(name_str)
                if (trim(temp_str1).eq.trim(name_str)) then

                    !Create the string after the equals sign
                    temp_str2=temp_str(index_val+1:)
                    temp_str2=adjustl(temp_str2)

                    if (len(trim(temp_str2)).ge.1) then
                        read(temp_str2,*,ERR=20) read_name_logical
                        write(unit_out,'(A,L)') 'Setting: '//trim(name_str)//' = ',read_name_logical
                    endif
                
                endif
    
            endif
        endif
    
    
20  enddo    
    
30    end function read_name_logical
