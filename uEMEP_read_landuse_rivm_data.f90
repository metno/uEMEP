!uEMEP_read_landuse_rivm_data.f90
    
    subroutine uEMEP_read_landuse_rivm_data
    
    use uEMEP_definitions

    implicit none

    integer i,j
    integer ncols_sub,nrows_sub
    real cellsize_sub,xll_corner_sub,yll_corner_sub
    real, allocatable :: landuse_array(:,:)
    integer depac_index(9)
    integer emep_landuse_index
    integer exists
    
    !landuse_subgrid=0
    !landuse_subgrid(:,:,temp_decid_index)=1.
    
    pathfilename_landuse=trim(pathname_landuse)//trim(filename_landuse)
    inquire(file=trim(pathfilename_landuse),exist=exists)
    if (.not.exists) then
        write(unit_logfile,'(A,A)') ' ERROR: Landuse file does not exist: ', trim(pathfilename_landuse)
        stop
    endif
    
    call read_esri_ascii_header(unit_logfile,pathfilename_landuse,ncols_sub,nrows_sub,cellsize_sub,xll_corner_sub,yll_corner_sub,.false.)

    landuse_subgrid_dim(x_dim_index)=ncols_sub
    landuse_subgrid_dim(y_dim_index)=nrows_sub
    landuse_subgrid_min(x_dim_index)=xll_corner_sub
    landuse_subgrid_min(y_dim_index)=yll_corner_sub
    landuse_subgrid_delta(x_dim_index)=cellsize_sub
    landuse_subgrid_delta(y_dim_index)=cellsize_sub
    
    !Deallocate grids if they are already allocated.
    if (allocated(landuse_subgrid)) deallocate (landuse_subgrid)
    if (allocated(x_landuse_subgrid)) deallocate (x_landuse_subgrid)
    if (allocated(y_landuse_subgrid)) deallocate (y_landuse_subgrid)
    if (allocated(lon_landuse_subgrid)) deallocate (lon_landuse_subgrid)
    if (allocated(lat_landuse_subgrid)) deallocate (lat_landuse_subgrid)
    if (allocated(xproj_landuse_subgrid)) deallocate (xproj_landuse_subgrid)
    if (allocated(yproj_landuse_subgrid)) deallocate (yproj_landuse_subgrid)
    
    !Reefine landuse grid
    if (.not.allocated(landuse_subgrid)) allocate (landuse_subgrid(landuse_subgrid_dim(x_dim_index),landuse_subgrid_dim(y_dim_index),n_landuse_index))
    if (.not.allocated(x_landuse_subgrid)) allocate (x_landuse_subgrid(landuse_subgrid_dim(x_dim_index),landuse_subgrid_dim(y_dim_index)))
    if (.not.allocated(y_landuse_subgrid)) allocate (y_landuse_subgrid(landuse_subgrid_dim(x_dim_index),landuse_subgrid_dim(y_dim_index)))
    if (.not.allocated(lon_landuse_subgrid)) allocate (lon_landuse_subgrid(landuse_subgrid_dim(x_dim_index),landuse_subgrid_dim(y_dim_index)))
    if (.not.allocated(lat_landuse_subgrid)) allocate (lat_landuse_subgrid(landuse_subgrid_dim(x_dim_index),landuse_subgrid_dim(y_dim_index)))
    if (.not.allocated(xproj_landuse_subgrid)) allocate (xproj_landuse_subgrid(landuse_subgrid_dim(x_dim_index),landuse_subgrid_dim(y_dim_index)))
    if (.not.allocated(yproj_landuse_subgrid)) allocate (yproj_landuse_subgrid(landuse_subgrid_dim(x_dim_index),landuse_subgrid_dim(y_dim_index)))

    !Set the depsoition grid to be the same dimensions as the landuse grid
    deposition_subgrid_dim(x_dim_index)=ncols_sub
    deposition_subgrid_dim(y_dim_index)=nrows_sub
    deposition_subgrid_min(x_dim_index)=xll_corner_sub
    deposition_subgrid_min(y_dim_index)=yll_corner_sub
    deposition_subgrid_delta(x_dim_index)=cellsize_sub
    deposition_subgrid_delta(y_dim_index)=cellsize_sub
    
    
    !Deallocate grids if they are already allocated.
    if (allocated(deposition_subgrid)) deallocate (deposition_subgrid)
    if (allocated(x_deposition_subgrid)) deallocate (x_deposition_subgrid)
    if (allocated(y_deposition_subgrid)) deallocate (y_deposition_subgrid)
    if (allocated(lon_deposition_subgrid)) deallocate (lon_deposition_subgrid)
    if (allocated(lat_deposition_subgrid)) deallocate (lat_deposition_subgrid)
    if (allocated(xproj_deposition_subgrid)) deallocate (xproj_deposition_subgrid)
    if (allocated(yproj_deposition_subgrid)) deallocate (yproj_deposition_subgrid)
    
    !Reefine deposition grid
    if (.not.allocated(deposition_subgrid)) allocate (deposition_subgrid(deposition_subgrid_dim(x_dim_index),deposition_subgrid_dim(y_dim_index),deposition_subgrid_dim(t_dim_index),n_deposition_index,n_pollutant_loop))
    if (.not.allocated(x_deposition_subgrid)) allocate (x_deposition_subgrid(deposition_subgrid_dim(x_dim_index),deposition_subgrid_dim(y_dim_index)))
    if (.not.allocated(y_deposition_subgrid)) allocate (y_deposition_subgrid(deposition_subgrid_dim(x_dim_index),deposition_subgrid_dim(y_dim_index)))
    if (.not.allocated(lon_deposition_subgrid)) allocate (lon_deposition_subgrid(deposition_subgrid_dim(x_dim_index),deposition_subgrid_dim(y_dim_index)))
    if (.not.allocated(lat_deposition_subgrid)) allocate (lat_deposition_subgrid(deposition_subgrid_dim(x_dim_index),deposition_subgrid_dim(y_dim_index)))
    if (.not.allocated(xproj_deposition_subgrid)) allocate (xproj_deposition_subgrid(deposition_subgrid_dim(x_dim_index),deposition_subgrid_dim(y_dim_index)))
    if (.not.allocated(yproj_deposition_subgrid)) allocate (yproj_deposition_subgrid(deposition_subgrid_dim(x_dim_index),deposition_subgrid_dim(y_dim_index)))
    
    
    do j=1,landuse_subgrid_dim(y_dim_index)
    do i=1,landuse_subgrid_dim(x_dim_index)                 

        x_landuse_subgrid(i,j)=landuse_subgrid_min(x_dim_index)+landuse_subgrid_delta(x_dim_index)*(i-0.5)
        y_landuse_subgrid(i,j)=landuse_subgrid_min(y_dim_index)+landuse_subgrid_delta(y_dim_index)*(j-0.5)
        
        !Set the lat-lon coordinates of the landuse
        call PROJ2LL(x_landuse_subgrid(i,j),y_landuse_subgrid(i,j),lon_landuse_subgrid(i,j),lat_landuse_subgrid(i,j),projection_attributes,projection_type)
        !if (projection_type.eq.RDM_projection_index) then
        !    call RDM2LL(y_landuse_subgrid(i,j),x_landuse_subgrid(i,j),lat_landuse_subgrid(i,j),lon_landuse_subgrid(i,j))
        !elseif (projection_type.eq.UTM_projection_index) then
        !    call UTM2LL(utm_zone,y_landuse_subgrid(i,j),x_landuse_subgrid(i,j),lat_landuse_subgrid(i,j),lon_landuse_subgrid(i,j))
        !endif
        
        !If the EMEP projection is lambert then set the proj coordinates to lambert, otherwise to lat-lon
        if (EMEP_projection_type.eq.LCC_projection_index) then
            call lb2lambert2_uEMEP(xproj_landuse_subgrid(i,j),yproj_landuse_subgrid(i,j),lon_landuse_subgrid(i,j),lat_landuse_subgrid(i,j),EMEP_projection_attributes)
        else
            xproj_landuse_subgrid(i,j)=lon_landuse_subgrid(i,j)
            yproj_landuse_subgrid(i,j)=lat_landuse_subgrid(i,j)            
        endif

    enddo
    enddo
    
    !Set the deposition x and y to be the same
    x_deposition_subgrid=x_landuse_subgrid
    y_deposition_subgrid=y_landuse_subgrid
    xproj_deposition_subgrid=xproj_landuse_subgrid
    yproj_deposition_subgrid=yproj_landuse_subgrid
    lon_deposition_subgrid=lon_landuse_subgrid
    lat_deposition_subgrid=lat_landuse_subgrid
   
   
    !Recalculate the cross references again since these could have changed
    call uEMEP_crossreference_grids
    
    !Read the landuse index into the temporary landuse array
    if (.not.allocated(landuse_array)) allocate (landuse_array(landuse_subgrid_dim(x_dim_index),landuse_subgrid_dim(y_dim_index)))

    call read_esri_ascii_file(unit_logfile,pathfilename_landuse,ncols_sub,nrows_sub,cellsize_sub,landuse_array,x_landuse_subgrid,y_landuse_subgrid,.false.)

    !set the depac indicies to the matching EMEP ones
    depac_index=0
    depac_index(4)=temp_conif_index
    depac_index(5)=temp_decid_index
    depac_index(2)=temp_crop_index
    depac_index(3)=temp_crop_index
    depac_index(8)=moorland_index
    depac_index(1)=grass_index
    depac_index(9)=desert_index
    depac_index(6)=water_index
    depac_index(7)=urban_index
    
    !Distribute to the depac indexes to the EMEP ones
    landuse_subgrid=0
    do j=1,landuse_subgrid_dim(y_dim_index)
    do i=1,landuse_subgrid_dim(x_dim_index) 
        emep_landuse_index=depac_index(int(landuse_array(i,j)))
        landuse_subgrid(i,j,emep_landuse_index)=1.
        !Put the landuse index in the last array
        landuse_subgrid(i,j,grid_index)=emep_landuse_index
        !write(*,*) i,j,int(landuse_array(i,j)), emep_landuse_index
    enddo
    enddo
    
    if (allocated(landuse_array)) deallocate (landuse_array)
    
    end subroutine uEMEP_read_landuse_rivm_data
