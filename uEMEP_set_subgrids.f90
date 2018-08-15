!uEMEP_set_subgrids.f90
    
    subroutine uEMEP_set_subgrids

    use uEMEP_definitions

    implicit none

    integer i,j,k    
    
    !Reset min and max with the buffer and calculate dimensions
    !subgrid_min(x_dim_index)=subgrid_min(x_dim_index)-buffer(x_dim_index);subgrid_min(y_dim_index)=subgrid_min(y_dim_index)-buffer(y_dim_index)
    !subgrid_max(x_dim_index)=subgrid_max(x_dim_index)+buffer(x_dim_index);subgrid_max(y_dim_index)=subgrid_max(y_dim_index)+buffer(y_dim_index)
    subgrid_dim(x_dim_index)=floor((subgrid_max(x_dim_index)-subgrid_min(x_dim_index))/subgrid_delta(x_dim_index))
    subgrid_dim(y_dim_index)=floor((subgrid_max(y_dim_index)-subgrid_min(y_dim_index))/subgrid_delta(y_dim_index))
 
    !write(*,*) subgrid_dim(x_dim_index),subgrid_max(x_dim_index),subgrid_min(x_dim_index),subgrid_delta(x_dim_index)
    !write(*,*) subgrid_dim(y_dim_index),subgrid_max(y_dim_index),subgrid_min(y_dim_index),subgrid_delta(y_dim_index)
    
    !Set all integral subgrids relative to the target subgrid
    if (integral_subgrid_delta_ref.eq.0.) then
        integral_subgrid_delta=subgrid_delta*integral_subgrid_step
    else
        integral_subgrid_delta(x_dim_index)=max(integral_subgrid_delta_ref,subgrid_delta(x_dim_index))
        integral_subgrid_delta(y_dim_index)=max(integral_subgrid_delta_ref,subgrid_delta(y_dim_index))
        integral_subgrid_step=floor(integral_subgrid_delta_ref/subgrid_delta(x_dim_index)+.5)
    endif
    
    integral_subgrid_min=subgrid_min
    integral_subgrid_max=subgrid_max
    integral_subgrid_dim(x_dim_index)=floor((subgrid_max(x_dim_index)-subgrid_min(x_dim_index))/integral_subgrid_delta(x_dim_index))
    integral_subgrid_dim(y_dim_index)=floor((subgrid_max(y_dim_index)-subgrid_min(y_dim_index))/integral_subgrid_delta(y_dim_index))
    integral_subgrid_dim(t_dim_index)=subgrid_dim(t_dim_index)
    !Set the integral subgrid dimensions so they cannot be larger than the target subgrid
    integral_subgrid_dim(x_dim_index)=min(integral_subgrid_dim(x_dim_index),subgrid_dim(x_dim_index))
    integral_subgrid_dim(y_dim_index)=min(integral_subgrid_dim(y_dim_index),subgrid_dim(y_dim_index))

    !Set all population subgrids relative to the target subgrid
    if (population_data_type.eq.population_index) then
        !When 250 m population data is used then set this as a limit
        population_subgrid_delta(x_dim_index)=max(subgrid_delta(x_dim_index),limit_population_delta)
        population_subgrid_delta(y_dim_index)=max(subgrid_delta(y_dim_index),limit_population_delta)
    else
        !Allow the point population data to have the same grid as the target grid
        population_subgrid_delta(x_dim_index)=subgrid_delta(x_dim_index)
        population_subgrid_delta(y_dim_index)=subgrid_delta(y_dim_index)
    endif        
    
    population_subgrid_min=subgrid_min
    population_subgrid_max=subgrid_max
    population_subgrid_dim(x_dim_index)=floor((population_subgrid_max(x_dim_index)-population_subgrid_min(x_dim_index))/population_subgrid_delta(x_dim_index))
    population_subgrid_dim(y_dim_index)=floor((population_subgrid_max(y_dim_index)-population_subgrid_min(y_dim_index))/population_subgrid_delta(y_dim_index))
    !Set the population subgrid dimensions so they cannot be larger than the target subgrid. Not certain why I do this.
    population_subgrid_dim(x_dim_index)=min(population_subgrid_dim(x_dim_index),subgrid_dim(x_dim_index))
    population_subgrid_dim(y_dim_index)=min(population_subgrid_dim(y_dim_index),subgrid_dim(y_dim_index))
    !Set population subgrid so it has a minimum of 1 dimensions, to avoid problems when running receptor calculations
    population_subgrid_dim(x_dim_index)=max(population_subgrid_dim(x_dim_index),1)
    population_subgrid_dim(y_dim_index)=max(population_subgrid_dim(y_dim_index),1)
    

    !Set all emission subgrids to be the same as the target subgrid
    emission_max_subgrid_dim=subgrid_dim
    do i=1,n_source_index
        emission_subgrid_delta(:,i)=subgrid_delta
        emission_subgrid_min(:,i)=subgrid_min
        emission_subgrid_max(:,i)=subgrid_max
        emission_subgrid_dim(:,i)=subgrid_dim
    enddo

    !Set shipping data to a minimum vale for all sources (Cannot be smaller than the target subgrid)
    emission_subgrid_delta(x_dim_index,shipping_index)=max(subgrid_delta(x_dim_index),limit_shipping_delta)
    emission_subgrid_delta(y_dim_index,shipping_index)=max(subgrid_delta(y_dim_index),limit_shipping_delta)
    emission_subgrid_delta(x_dim_index,heating_index)=max(subgrid_delta(x_dim_index),limit_heating_delta)
    emission_subgrid_delta(y_dim_index,heating_index)=max(subgrid_delta(y_dim_index),limit_heating_delta)
    
    !Set all the emission subgrid dimmensions after changes
    do i=1,n_source_index
        emission_subgrid_dim(x_dim_index,i)=floor((emission_subgrid_max(x_dim_index,i)-emission_subgrid_min(x_dim_index,i))/emission_subgrid_delta(x_dim_index,i))
        emission_subgrid_dim(y_dim_index,i)=floor((emission_subgrid_max(y_dim_index,i)-emission_subgrid_min(y_dim_index,i))/emission_subgrid_delta(y_dim_index,i))
        write(unit_logfile,'(A,I6,A5,2I6)') 'Emission grid dimensions for source ',i,': ',emission_subgrid_dim(1:2,i)
    enddo
    
    write(unit_logfile,'(A,2I6)') 'Target grid dimensions: ',subgrid_dim(1:2)
    write(unit_logfile,'(A,2I6)') 'Integral grid dimensions: ',integral_subgrid_dim(1:2) 
    write(unit_logfile,'(A,2f10.1)') 'Target subgrid grid sizes: ',subgrid_delta
    write(unit_logfile,'(A,I6,2f10.1)') 'Integral subgrid step and grid sizes: ',integral_subgrid_step,integral_subgrid_delta
    
    end subroutine uEMEP_set_subgrids