module mod_read_esri_ascii_file

    implicit none
    private

    public :: write_esri_ascii_file, read_esri_ascii_file, read_esri_ascii_header

contains

!read_esri_ascii_file.f90
    
    subroutine read_esri_ascii_header(unit_logfile,filename_ascii_sub,ncols_sub,nrows_sub,cellsize_sub,xllcorner,yllcorner,read_nodata_flag)
    
    implicit none
    character(*) filename_ascii_sub
    character(256) temp_str
    integer i,j,ii,jj
    integer ncols_sub,nrows_sub
    integer ncols_sub_temp,nrows_sub_temp
    real cellsize_sub
    real xllcorner
    real yllcorner
    real :: NODATA_value=-999.
    integer :: unit_in=20
    integer unit_logfile
    logical read_nodata_flag

    !write(unit_logfile,'(a)') ' Opening ascii file: '//trim(filename_ascii_sub)
    open(unit_in,file=filename_ascii_sub,access='sequential',form='formatted',status='old',readonly)
    
        rewind(unit_in)

        read(unit_in,*)temp_str,ncols_sub
        !write(*,*)trim(temp_str),ncols
        read(unit_in,*)temp_str,nrows_sub
        !write(*,*)trim(temp_str),nrows
        read(unit_in,*)temp_str,xllcorner
        !write(*,*)trim(temp_str),xllcorner
        read(unit_in,*)temp_str,yllcorner
        !write(*,*)trim(temp_str),yllcorner
        read(unit_in,*)temp_str,cellsize_sub
        !write(*,*)trim(temp_str),cellsize
        if (read_nodata_flag) read(unit_in,*)temp_str,NODATA_value
        !write(*,*)trim(temp_str),NODATA_value
        write(unit_logfile,'(2a10,4a12)')'ncols','nrows','xllcorner','yllcorner','cellsize','NODATA_val'
        write(unit_logfile,'(2i10,4f12.1)')ncols_sub,nrows_sub,xllcorner,yllcorner,cellsize_sub,NODATA_value
    
        close(unit_in)

    end subroutine read_esri_ascii_header

    
    subroutine read_esri_ascii_file(unit_logfile,filename_ascii_sub,ncols_sub,nrows_sub,cellsize_sub,val_array,x_array,y_array,read_nodata_flag)
    
    implicit none
    character(*) filename_ascii_sub
    character(256) temp_str
    integer i,j,ii,jj
    integer ncols_sub,nrows_sub
    integer ncols_sub_temp,nrows_sub_temp
    real cellsize_sub
    real xllcorner
    real yllcorner
    real NODATA_value
    integer :: unit_in=20
    real val_array(ncols_sub,nrows_sub)
    real x_array(ncols_sub,nrows_sub)
    real y_array(ncols_sub,nrows_sub)
    integer unit_logfile
    logical read_nodata_flag

    !write(unit_logfile,'(a)') ' Opening ascii file: '//trim(filename_ascii_sub)
    open(unit_in,file=filename_ascii_sub,access='sequential',form='formatted',status='old',readonly)
    
        rewind(unit_in)

        read(unit_in,*)temp_str,ncols_sub_temp
        !write(*,*)trim(temp_str),ncols
        read(unit_in,*)temp_str,nrows_sub_temp
        !write(*,*)trim(temp_str),nrows
        read(unit_in,*)temp_str,xllcorner
        !write(*,*)trim(temp_str),xllcorner
        read(unit_in,*)temp_str,yllcorner
        !write(*,*)trim(temp_str),yllcorner
        read(unit_in,*)temp_str,cellsize_sub
        !write(*,*)trim(temp_str),cellsize
        if (read_nodata_flag) read(unit_in,*)temp_str,NODATA_value
        !write(*,*)trim(temp_str),NODATA_value
        !write(unit_logfile,'(2a10,4a12)')'ncols','nrows','xllcorner','yllcorner','cellsize','NODATA_val'
        !write(unit_logfile,'(2i10,4f12.1)')ncols_sub,nrows_sub,xllcorner,yllcorner,cellsize_sub,NODATA_value
    
        !Reduce for testing
        !nrows=1000
    
        !if (.not.allocated(val_array)) allocate(val_array(ncols_sub,nrows_sub))
        !if (.not.allocated(x_array)) allocate(x_array(ncols_sub,nrows_sub))
        !if (.not.allocated(y_array)) allocate(y_array(ncols_sub,nrows_sub))
        
        if (ncols_sub.ne.ncols_sub_temp.or.nrows_sub.ne.nrows_sub_temp) then
            write(unit_logfile,*)'ERROR: Subgrid columns or rows do not match. Columns: ',ncols_sub,ncols_sub_temp,'Rows: ',nrows_sub,nrows_sub_temp
            stop
        endif
    
        read(unit_in,*) ((val_array(ii,jj),ii=1,ncols_sub),jj=nrows_sub,1,-1)
    
        close(unit_in)

        j=1
        !write(*,'(<ncols>f10.1)') (array(ii,j),ii=1,ncols)
        !write(*,*)
        !j=nrows
        !write(*,'(<ncols>f10.1)') (array(ii,j),ii=1,ncols)
        !write(unit_logfile,'(a,f6.1,a,f6.1)') ' Min array val: ',minval(array),' Max array val: ',maxval(array)

        !Set position arrays
        do i=1,ncols_sub
         do j=1,nrows_sub
            x_array(i,j)=xllcorner+cellsize_sub/2.+(i-1)*cellsize_sub
            y_array(i,j)=yllcorner+cellsize_sub/2.+(j-1)*cellsize_sub
         enddo
       enddo

    end subroutine read_esri_ascii_file

    subroutine write_esri_ascii_file(unit_logfile,filename_ascii_sub,ncols_sub,nrows_sub,cellsize_sub,val_array,x_array,y_array)
    
    implicit none
    character(*) filename_ascii_sub
    character(256) temp_str
    integer i,j,ii,jj
    integer ncols_sub,nrows_sub
    real cellsize_sub
    real xllcorner
    real yllcorner
    real :: NODATA_value=-999.
    integer :: unit_in=20
    real val_array(ncols_sub,nrows_sub)
    real x_array(ncols_sub,nrows_sub)
    real y_array(ncols_sub,nrows_sub)
    integer unit_logfile

    xllcorner=x_array(1,1)-cellsize_sub/2.
    yllcorner=y_array(1,1)-cellsize_sub/2.
    !write(unit_logfile,'(a)') ' Writing to ascii file: '//trim(filename_ascii_sub)
    open(unit_in,file=filename_ascii_sub,access='sequential',form='formatted',status='unknown')
    
        write(unit_in,*)'ncols',ncols_sub
        write(unit_in,*)'nrows',nrows_sub
        write(unit_in,*)'xllcorner',xllcorner
        write(unit_in,*)'yllcorner',yllcorner
        write(unit_in,*)'cellsize',cellsize_sub
        write(unit_in,*)'NODATA_value',NODATA_value
        !write(unit_logfile,'(2a10,4a12)')'ncols','nrows','xllcorner','yllcorner','cellsize','NODATA_val'
        !write(unit_logfile,'(2i10,4f12.1)')ncols_sub,nrows_sub,xllcorner,yllcorner,cellsize_sub,NODATA_value  
    
        do jj=nrows_sub,1,-1
            write(unit_in,'(<ncols_sub>es12.3)') (val_array(ii,jj),ii=1,ncols_sub)
        enddo
          
    close(unit_in)

    end subroutine write_esri_ascii_file

    
    subroutine read_esri_ascii_3d_file(unit_logfile,filename_ascii_sub,ncols_sub,nrows_sub,nblocks_sub,cellsize_sub,val_array,x_array,y_array)
    
    implicit none
    character(*) filename_ascii_sub
    character(256) temp_str
    integer i,j,ii,jj,tt
    integer ncols_sub,nrows_sub,nblocks_sub
    integer ncols_sub_temp,nrows_sub_temp,nblocks_sub_temp
    real nrows_sub_temp_real
    real cellsize_sub
    real xllcorner
    real yllcorner
    real NODATA_value
    integer :: unit_in=20
    real val_array(ncols_sub,nrows_sub,nblocks_sub)
    real x_array(ncols_sub,nrows_sub)
    real y_array(ncols_sub,nrows_sub)
    integer unit_logfile

    !write(unit_logfile,'(a)') ' Opening ascii file: '//trim(filename_ascii_sub)
    open(unit_in,file=filename_ascii_sub,access='sequential',form='formatted',status='old',readonly)
    
        rewind(unit_in)

        read(unit_in,*)temp_str,ncols_sub_temp
        !write(*,*)trim(temp_str),ncols
        read(unit_in,*)temp_str,nrows_sub_temp
        !write(*,*)trim(temp_str),nrows
        read(unit_in,*)temp_str,nrows_sub_temp_real
        !write(*,*)trim(temp_str),ncols
        if (temp_str.eq.'xllcorner') then
            nblocks_sub_temp=1
            xllcorner=nrows_sub_temp_real
        else
            nrows_sub_temp=int(nrows_sub_temp_real)
            read(unit_in,*)temp_str,xllcorner
            !write(*,*)trim(temp_str),xllcorner
       endif
        read(unit_in,*)temp_str,yllcorner
        !write(*,*)trim(temp_str),yllcorner
        read(unit_in,*)temp_str,cellsize_sub
        !write(*,*)trim(temp_str),cellsize
        read(unit_in,*)temp_str,NODATA_value
        !write(*,*)trim(temp_str),NODATA_value
        !write(unit_logfile,'(2a10,4a12)')'ncols','nrows','xllcorner','yllcorner','cellsize','NODATA_val'
        !write(unit_logfile,'(2i10,4f12.1)')ncols_sub,nrows_sub,xllcorner,yllcorner,cellsize_sub,NODATA_value
    
        
        if (ncols_sub.ne.ncols_sub_temp.or.nrows_sub.ne.nrows_sub_temp.or.nblocks_sub.ne.nblocks_sub_temp) then
            write(unit_logfile,'(A,2I,A,2I,A,2I)')'ERROR: Subgrid columns or rows do not match. Columns: ',ncols_sub,ncols_sub_temp,'  Rows: ',nrows_sub,nrows_sub_temp,'  Blocks: ',nblocks_sub,nblocks_sub_temp
            stop
        endif
    
        read(unit_in,*) (((val_array(ii,jj,tt),ii=1,ncols_sub),jj=nrows_sub,1,-1),tt=1,nblocks_sub)
    
        close(unit_in)

        j=1
        !write(*,'(<ncols>f10.1)') (array(ii,j),ii=1,ncols)
        !write(*,*)
        !j=nrows
        !write(*,'(<ncols>f10.1)') (array(ii,j),ii=1,ncols)
        !write(unit_logfile,'(a,f6.1,a,f6.1)') ' Min array val: ',minval(array),' Max array val: ',maxval(array)

        !Set position arrays
        do i=1,ncols_sub
         do j=1,nrows_sub
            x_array(i,j)=xllcorner+cellsize_sub/2.+(i-1)*cellsize_sub
            y_array(i,j)=yllcorner+cellsize_sub/2.+(j-1)*cellsize_sub
         enddo
       enddo

    end subroutine read_esri_ascii_3d_file

    subroutine write_esri_ascii_3d_file(unit_logfile,filename_ascii_sub,ncols_sub,nrows_sub,nblocks_sub,cellsize_sub,val_array,x_array,y_array)
    
    implicit none
    character(*) filename_ascii_sub
    character(256) temp_str
    integer i,j,ii,jj,tt
    integer ncols_sub,nrows_sub,nblocks_sub
    real cellsize_sub
    real xllcorner
    real yllcorner
    real :: NODATA_value=-999.
    integer :: unit_in=20
    real val_array(ncols_sub,nrows_sub,nblocks_sub)
    real x_array(ncols_sub,nrows_sub)
    real y_array(ncols_sub,nrows_sub)
    integer unit_logfile

    xllcorner=x_array(1,1)-cellsize_sub/2.
    yllcorner=y_array(1,1)-cellsize_sub/2.
    !write(unit_logfile,'(a)') ' Writing to ascii file: '//trim(filename_ascii_sub)
    open(unit_in,file=filename_ascii_sub,access='sequential',form='formatted',status='unknown')
    
        write(unit_in,*)'ncols',ncols_sub
        write(unit_in,*)'nrows',nrows_sub
        if (nblocks_sub.gt.1) write(unit_in,*)'nblocks',nblocks_sub
        write(unit_in,*)'xllcorner',xllcorner
        write(unit_in,*)'yllcorner',yllcorner
        write(unit_in,*)'cellsize',cellsize_sub
        write(unit_in,*)'NODATA_val',NODATA_value
        !write(unit_logfile,'(2a10,4a12)')'ncols','nrows','xllcorner','yllcorner','cellsize','NODATA_val'
        !write(unit_logfile,'(3i10,4f12.1)')ncols_sub,nrows_sub,nblocks_sub,xllcorner,yllcorner,cellsize_sub,NODATA_value  
    
        do tt=1,nblocks_sub
        do jj=nrows_sub,1,-1
            write(unit_in,'(<ncols_sub>es12.3)') (val_array(ii,jj,tt),ii=1,ncols_sub)
        enddo
        enddo
          
    close(unit_in)

    end subroutine write_esri_ascii_3d_file

end module mod_read_esri_ascii_file

