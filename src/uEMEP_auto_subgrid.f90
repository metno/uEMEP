module auto_subgrid

    use mod_area_interpolation, only: area_weighted_interpolation_function
    use uEMEP_definitions, only: allsource_index, use_receptor_positions_for_auto_subgrid_flag, &
        use_emission_positions_for_auto_subgrid_flag, use_population_positions_for_auto_subgrid_flag, &
        use_subgrid, unit_logfile, use_subgrid_val, use_subgrid_interpolation_index, &
        max_interpolation_subgrid_size, traffic_index, shipping_index, n_use_subgrid_levels, &
        x_dim_index, subgrid_delta, use_subgrid_step_delta, calculate_source, source_file_str, &
        emission_subgrid_delta, y_dim_index, subgrid_dim, crossreference_target_to_emission_subgrid, &
        emission_subgrid_dim, proxy_emission_subgrid, interpolate_subgrids_flag, n_source_index, &
        x_subgrid, y_subgrid, outside_interpolation_region_index, t_dim_index, n_pollutant_loop, &
        proxy_subgrid_index, subgrid, traveltime_subgrid, outside_region_index, municipality_index, &
        filename_population, region_index, pathfilename_population, pathname_population, region_id, &
        subgrid_min, trace_emissions_from_in_region, emission_subgrid_min, x_emission_subgrid, &
        y_emission_subgrid, use_subgrid_region

    implicit none
    private

    public :: uEMEP_auto_subgrid, uEMEP_region_mask, uEMEP_interpolate_auto_subgrid

contains

!uEMEP_auto_subgrid.f90

!==========================================================================
!   uEMEP model uEMEP_auto_subgrid
!   Automatically creates a grid dependent on the distance to source
!
!==========================================================================
    subroutine uEMEP_auto_subgrid


        implicit none

        integer i,j,k
        integer         i_source
        integer         t

        real :: max_use_subgrid_size(n_source_index)
        real :: use_subgrid_range(n_source_index)
        integer use_subgrid_step
        integer use_emission_subgrid_space
        integer i_cross,j_cross
        integer i_start,i_end,j_start,j_end
        integer i_start_k,i_end_k,j_start_k,j_end_k
        real sum_emission

        !Exit if emission positions are not to be used and if one of the other auto positions is to be used
        if (.not.use_emission_positions_for_auto_subgrid_flag(allsource_index).and.(use_population_positions_for_auto_subgrid_flag.or.use_receptor_positions_for_auto_subgrid_flag)) then
            return
        endif
        !Exit and fill subgrid if none of the auto subgrids are to be used
        if (.not.use_emission_positions_for_auto_subgrid_flag(allsource_index).and..not.use_population_positions_for_auto_subgrid_flag.and..not.use_receptor_positions_for_auto_subgrid_flag) then
            use_subgrid=.true.
            return
        endif


        write(unit_logfile,'(A)') ''
        write(unit_logfile,'(A)') '================================================================'
        write(unit_logfile,'(A)') 'Setting subgrids automatically (uEMEP_auto_subgrid)'
        write(unit_logfile,'(A)') '================================================================'

        !Set all subgrids to do not use
        use_subgrid_val=0
        use_subgrid_interpolation_index=-1

        !Set time index used for emissions to 1, so only tests for the first hour if there are emissions
        t=1
        !Set the maximum grid size
        max_use_subgrid_size=max_interpolation_subgrid_size
        !max_use_subgrid_size(traffic_index)=1000.
        !max_use_subgrid_size(shipping_index)=1000.
        use_subgrid_range=8.
        use_subgrid_range(traffic_index)=4.
        use_subgrid_range(shipping_index)=12.

        n_use_subgrid_levels=0
        !Sets the auto gridding
        if (subgrid_delta(x_dim_index).eq.100.) then
            use_subgrid_step_delta(0)=1
            use_subgrid_step_delta(1)=2
            use_subgrid_step_delta(2)=5
            use_subgrid_step_delta(3)=10
            use_subgrid_step_delta(4)=20
            if (max_interpolation_subgrid_size.eq.250.) n_use_subgrid_levels=1
            if (max_interpolation_subgrid_size.eq.500.) n_use_subgrid_levels=2
            if (max_interpolation_subgrid_size.eq.1000.) n_use_subgrid_levels=3
            if (max_interpolation_subgrid_size.eq.2000.) n_use_subgrid_levels=4
            use_subgrid_range(traffic_index)=6.
        elseif (subgrid_delta(x_dim_index).eq.25.) then
            use_subgrid_step_delta(0)=1
            use_subgrid_step_delta(1)=2
            use_subgrid_step_delta(2)=4
            use_subgrid_step_delta(3)=10
            use_subgrid_step_delta(4)=20
            use_subgrid_step_delta(5)=40
            use_subgrid_step_delta(5)=80
            if (max_interpolation_subgrid_size.eq.250.) n_use_subgrid_levels=2
            if (max_interpolation_subgrid_size.eq.500.) n_use_subgrid_levels=3
            if (max_interpolation_subgrid_size.eq.1000.) n_use_subgrid_levels=4
            if (max_interpolation_subgrid_size.eq.2000.) n_use_subgrid_levels=5
            use_subgrid_range(traffic_index)=4.
        elseif (subgrid_delta(x_dim_index).eq.50.) then
            use_subgrid_step_delta(0)=1
            use_subgrid_step_delta(1)=2
            use_subgrid_step_delta(2)=4
            use_subgrid_step_delta(3)=10
            use_subgrid_step_delta(4)=20
            use_subgrid_step_delta(5)=40
            if (max_interpolation_subgrid_size.eq.250.) n_use_subgrid_levels=2
            if (max_interpolation_subgrid_size.eq.500.) n_use_subgrid_levels=3
            if (max_interpolation_subgrid_size.eq.1000.) n_use_subgrid_levels=4
            if (max_interpolation_subgrid_size.eq.2000.) n_use_subgrid_levels=5
            use_subgrid_range(traffic_index)=6.
        elseif (subgrid_delta(x_dim_index).eq.250.) then
            use_subgrid_step_delta(0)=1
            use_subgrid_step_delta(1)=2
            use_subgrid_step_delta(2)=4
            use_subgrid_step_delta(3)=8
            if (max_interpolation_subgrid_size.eq.250.) n_use_subgrid_levels=0
            if (max_interpolation_subgrid_size.eq.500.) n_use_subgrid_levels=1
            if (max_interpolation_subgrid_size.eq.1000.) n_use_subgrid_levels=2
            if (max_interpolation_subgrid_size.eq.2000.) n_use_subgrid_levels=3
            use_subgrid_range(traffic_index)=8.
        else
            write(unit_logfile,'(a,f12.2)') 'When auto gridding target grid sizes must be either 25, 50, 100 or 250 metres. Stopping because target grid size is:',subgrid_delta(x_dim_index)
            stop
        endif

        if (n_use_subgrid_levels(allsource_index).eq.0) then
            write(unit_logfile,'(a,f12.2)') 'When auto gridding maximum grid sizes must be either 250, 500, 1000 or 2000 metres. Stopping because max grid size is:',max_interpolation_subgrid_size
            stop
        endif


        !Set the number of levels to match this
        !n_use_subgrid_levels=floor(log(max_use_subgrid_size/sqrt(subgrid_delta(x_dim_index)*subgrid_delta(y_dim_index)))/log(2.)+.5)
        do i_source=1,n_source_index
            if (calculate_source(i_source)) then
                write(*,*) 'Using auto subgrid for source ',trim(source_file_str(i_source)),use_emission_positions_for_auto_subgrid_flag(i_source)
            endif
        enddo

        do i_source=1,n_source_index
            if (calculate_source(i_source).and.use_emission_positions_for_auto_subgrid_flag(i_source)) then
                write(unit_logfile,'(a,2f10.1,i10)') trim(source_file_str(i_source))//': maximum use subgrid size (m), grid range and number of levels: ',max_use_subgrid_size(i_source),use_subgrid_range(i_source),n_use_subgrid_levels(i_source)
                !Fill the interpolation index with the highest level
                !use_subgrid_interpolation_index(:,:,i_source)=n_use_subgrid_levels(i_source)

                do k=n_use_subgrid_levels(i_source),0,-1
                    use_subgrid_step=2**k
                    use_subgrid_step=use_subgrid_step_delta(k)
                    !write(*,*) k,use_subgrid_step
                    use_emission_subgrid_space=floor(use_subgrid_step/sqrt(emission_subgrid_delta(x_dim_index,i_source)*emission_subgrid_delta(y_dim_index,i_source)/subgrid_delta(x_dim_index)/subgrid_delta(y_dim_index))/2.*use_subgrid_range(i_source))
                    !use_emission_subgrid_space=floor(use_subgrid_step/2.*use_subgrid_range(i_source))
                    use_emission_subgrid_space=max(1,use_emission_subgrid_space)
                    !write(*,*) k,use_subgrid_step,use_emission_subgrid_space,use_subgrid_step*subgrid_delta(y_dim_index),use_emission_subgrid_space*emission_subgrid_delta(y_dim_index,i_source)

                    !do j=use_subgrid_step,subgrid_dim(y_dim_index),use_subgrid_step
                    !do i=use_subgrid_step,subgrid_dim(x_dim_index),use_subgrid_step
                    do j=1,subgrid_dim(y_dim_index),use_subgrid_step
                        do i=1,subgrid_dim(x_dim_index),use_subgrid_step

                            i_cross=crossreference_target_to_emission_subgrid(i,j,x_dim_index,i_source)
                            j_cross=crossreference_target_to_emission_subgrid(i,j,y_dim_index,i_source)

                            !If the target grid is inside an emission grid with non zero emissions then include it
                            !if (sum(emission_subgrid(i_cross,j_cross,t,i_source,:)).gt.0) use_subgrid_val(i,j,i_source)=1

                            !Search in the neighbourhood and find the sum of the emissions
                            i_start=max(1,i_cross-use_emission_subgrid_space)
                            i_end=min(emission_subgrid_dim(x_dim_index,i_source),i_cross+use_emission_subgrid_space)
                            j_start=max(1,j_cross-use_emission_subgrid_space)
                            j_end=min(emission_subgrid_dim(y_dim_index,i_source),j_cross+use_emission_subgrid_space)
                            sum_emission=sum(proxy_emission_subgrid(i_start:i_end,j_start:j_end,i_source,:))

                            !Select those with emission sums > 0 but always include the coarsest level everywhere
                            if (sum_emission.gt.0..or.k.eq.n_use_subgrid_levels(i_source)) then
                                !if (sum_emission.gt.0..and.k.ne.n_use_subgrid_levels(i_source)) then
                                use_subgrid_val(i,j,i_source)=1
                                !Label the grids for interpolation later
                                i_start_k=max(1,i-use_subgrid_step/2)
                                i_end_k=min(subgrid_dim(x_dim_index),i+use_subgrid_step/2)
                                j_start_k=max(1,j-use_subgrid_step/2)
                                j_end_k=min(subgrid_dim(y_dim_index),j+use_subgrid_step/2)

                                use_subgrid_interpolation_index(i_start_k:i_end_k,j_start_k:j_end_k,i_source)=k
                                !write(*,'(6i)') k,i_start_k,i_end_k,j_start_k,j_end_k,sum(use_subgrid_interpolation_index(i_start_k:i_end_k,j_start_k:j_end_k,i_source))/(i_end_k-i_start_k+1)/(j_end_k-j_start_k+1)
                                !write(*,*) k,use_subgrid_step/2,use_subgrid_interpolation_index(i_start_k:i_end_k,j_start_k:j_end_k,i_source)
                            endif

                            if (i_source.eq.traffic_index) then
                                !write(*,'(7i,f)') k,i,j,i_start,i_end,j_start,j_end,sum_emission
                            endif


                            !Always include the coarsest level everywhere
                            !if (k.eq.n_use_subgrid_levels(i_source)) use_subgrid_val(i,j,i_source)=1

                        enddo
                    enddo
                enddo

                !Put points around the domain edge at lowest resolution for proper interpolation
                !k=n_use_subgrid_levels(i_source)
                !use_subgrid_step=2**k
                !do j=1,subgrid_dim(y_dim_index),use_subgrid_step
                !use_subgrid_val(1,j,i_source)=1
                !use_subgrid_val(subgrid_dim(x_dim_index),j,i_source)=1
                !enddo
                !do i=1,subgrid_dim(x_dim_index),use_subgrid_step
                !use_subgrid_val(i,1,i_source)=1
                !use_subgrid_val(i,subgrid_dim(y_dim_index),i_source)=1
                !enddo
                !use_subgrid_val(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),i_source)=1


                use_subgrid_val(:,:,allsource_index)=use_subgrid_val(:,:,allsource_index)+use_subgrid_val(:,:,i_source)


                !Check
                do j=1,subgrid_dim(y_dim_index)
                    do i=1,subgrid_dim(x_dim_index)
                        if (use_subgrid_interpolation_index(i,j,i_source).lt.0) write(*,*) i,j,use_subgrid_interpolation_index(i,j,i_source)
                    enddo
                enddo

            endif
        enddo

        use_subgrid_val=min(1,use_subgrid_val)
        !stop

        !Convert values to logical
        do i_source=1,n_source_index
            if (use_emission_positions_for_auto_subgrid_flag(i_source)) then
                if (calculate_source(i_source).or.i_source.eq.allsource_index) then
                    do j=1,subgrid_dim(y_dim_index)
                        do i=1,subgrid_dim(x_dim_index)
                            if (use_subgrid_val(i,j,i_source).gt.0) then
                                use_subgrid(i,j,i_source)=.true.
                            else
                                use_subgrid(i,j,i_source)=.false.
                            endif
                            !write(*,'(5i)') n_use_subgrid_levels(i_source),i_source,i,j,use_subgrid_interpolation_index(i,j,i_source)
                        enddo
                    enddo

                endif
            else
                !If not to be auto gridded then set use to true
                use_subgrid(:,:,i_source)=.true.
            endif

        enddo

        !If the grids are not to be interpolated then they must all be calculated at the same place for the different sources for interpolation later
        if (.not.interpolate_subgrids_flag) then
            do i_source=1,n_source_index
                !use_subgrid(:,:,i_source)=use_subgrid(:,:,allsource_index)
                !use_subgrid_val(:,:,i_source)=use_subgrid_val(:,:,allsource_index)
            enddo
        endif

        do i_source=1,n_source_index
            if (calculate_source(i_source).and.use_emission_positions_for_auto_subgrid_flag(i_source)) then
                write(unit_logfile,'(a,2i10,f6.1)') 'Number of calculation subgrids for '//trim(source_file_str(i_source))//' (number, total, percent):',sum(use_subgrid_val(:,:,i_source)),subgrid_dim(1)*subgrid_dim(2),sum(use_subgrid_val(:,:,i_source))*100./(subgrid_dim(1)*subgrid_dim(2))
            endif
        enddo

    end subroutine uEMEP_auto_subgrid

    subroutine uEMEP_interpolate_auto_subgrid
        !This is the corresponding routine for interpolating the auto selected data

        implicit none

        integer xdim,ydim
        real delta(2)
        real xval,yval
        real xgrid(3,3),ygrid(3,3),zgrid(3,3)

        integer i,j,k,ii,jj
        integer i_source,i_pollutant,t
        integer i_in(3),j_in(3)
        integer use_subgrid_step

        write(unit_logfile,'(A)') ''
        write(unit_logfile,'(A)') '================================================================'
        write(unit_logfile,'(A)') 'Interpolating auto grid (uEMEP_interpolate_auto_subgrid)'
        write(unit_logfile,'(A)') '================================================================'


        xdim=3
        ydim=3

        do i_source=1,n_source_index
            if (calculate_source(i_source).and.i_source.ne.allsource_index.and.use_emission_positions_for_auto_subgrid_flag(i_source)) then
                !Only interpolate for chosen sources that have been auto gridded based on emissions

                do k=n_use_subgrid_levels(i_source),1,-1
                    !use_subgrid_step=2**k
                    use_subgrid_step=use_subgrid_step_delta(k)
                    delta=use_subgrid_step*subgrid_delta


                    do j=1,subgrid_dim(y_dim_index)
                        do i=1,subgrid_dim(x_dim_index)
                            xval=x_subgrid(i,j)
                            yval=y_subgrid(i,j)
                            !Only do the interpolation if it is the right interpolation_index, it is inside the region and it is not a valid subgrid
                            !if (use_subgrid_interpolation_index(i,j,i_source).eq.k.and.use_subgrid_val(i,j,i_source).ne.outside_interpolation_region_index.and..not.use_subgrid_val(i,j,i_source)) then
                            !Do it everywhere in the grid
                            if (use_subgrid_interpolation_index(i,j,i_source).eq.k.and. &
                            !(use_subgrid_val(i,j,i_source).ne.outside_interpolation_region_index.or.use_subgrid_val(i,j,i_source).ne.outside_region_index).and. &
                                use_subgrid_val(i,j,i_source).ne.outside_interpolation_region_index.and. &
                                .not.use_subgrid(i,j,i_source)) then
                                !i_in(2)=((i-1)/use_subgrid_step)*use_subgrid_step+1
                                i_in(2)=floor(real(i-1)/use_subgrid_step+0.5)*use_subgrid_step+1
                                i_in(2)=min(subgrid_dim(x_dim_index),max(1,i_in(2)))
                                i_in(1)=i_in(2)-use_subgrid_step
                                i_in(3)=i_in(2)+use_subgrid_step
                                !j_in(2)=((j-1)/use_subgrid_step)*use_subgrid_step+1
                                j_in(2)=floor(real(j-1)/use_subgrid_step+0.5)*use_subgrid_step+1
                                j_in(2)=min(subgrid_dim(y_dim_index),max(1,j_in(2)))
                                j_in(1)=j_in(2)-use_subgrid_step
                                j_in(3)=j_in(2)+use_subgrid_step
                                if (i_in(1).lt.1) i_in(1)=i_in(2)
                                if (j_in(1).lt.1) j_in(1)=j_in(2)
                                if (i_in(3).gt.subgrid_dim(x_dim_index)) i_in(3)=i_in(2)
                                if (j_in(3).gt.subgrid_dim(y_dim_index)) j_in(3)=j_in(2)


                                do t=1,subgrid_dim(t_dim_index)
                                    do i_pollutant=1,n_pollutant_loop
                                        do jj=1,3
                                            do ii=1,3
                                                !xgrid(ii,jj)=x_subgrid(i_in(ii),j_in(jj))
                                                !ygrid(ii,jj)=y_subgrid(i_in(ii),j_in(jj))
                                                xgrid(ii,jj)=x_subgrid(i_in(2),j_in(2))+(ii-2)*delta(1)
                                                ygrid(ii,jj)=y_subgrid(i_in(2),j_in(2))+(jj-2)*delta(2)
                                                zgrid(ii,jj)=subgrid(i_in(ii),j_in(jj),t,proxy_subgrid_index,i_source,i_pollutant)
                                                if (i_in(ii).lt.1.or.i_in(ii).gt.subgrid_dim(x_dim_index).or.j_in(jj).lt.1.or.j_in(jj).gt.subgrid_dim(y_dim_index)) then
                                                    zgrid(ii,jj)=subgrid(i_in(2),j_in(2),t,proxy_subgrid_index,i_source,i_pollutant)
                                                elseif (use_subgrid_val(i_in(ii),j_in(jj),i_source).eq.outside_interpolation_region_index) then
                                                    zgrid(ii,jj)=subgrid(i_in(2),j_in(2),t,proxy_subgrid_index,i_source,i_pollutant)
                                                endif

                                                !If the middle grid is not valid (on the edge) then give it an average of the valid values surrounding it
                                                !count=0;val_temp=0.
                                                !if (use_subgrid_val(i_in(2),j_in(2),i_source).eq.outside_interpolation_region_index.and.ii.ne.2.and.jj.ne.2.and.use_subgrid_val(i_in(ii),j_in(jj),i_source).ne.outside_interpolation_region_index) then
                                                !    count=count+1
                                                !    val_temp=val_temp+zgrid(ii,jj)
                                                !    write(*,*) 'Middle interpolation grid outside'
                                                !endif
                                                !if (count.gt.0) then
                                                !    zgrid(ii,jj)=val_temp/count
                                                !endif

                                                !write(*,'(4i,4f12.2)') k,use_subgrid_step,ii,jj,xgrid(ii,jj)-xgrid(1,1),ygrid(ii,jj)-ygrid(1,1),xval-xgrid(1,1),yval-ygrid(1,1)
                                                !write(*,'(4i,4f12.2)') k,use_subgrid_step,ii,jj,xgrid(ii,jj),ygrid(ii,jj),xval,yval
                                            enddo
                                        enddo

                                        subgrid(i,j,t,proxy_subgrid_index,i_source,i_pollutant)=area_weighted_interpolation_function(xgrid,ygrid,zgrid,xdim,ydim,delta,xval,yval)
                                        !write(*,'(7i,f12.2)') k,t,i_pollutant,i,j,i_in(2),j_in(2),subgrid(i,j,t,proxy_subgrid_index,i_source,i_pollutant)

                                        !Travel time interpolation as well
                                        do jj=1,3
                                            do ii=1,3
                                                zgrid(ii,jj)=traveltime_subgrid(i_in(ii),j_in(jj),t,1,i_pollutant)
                                                if (i_in(ii).lt.1.or.i_in(ii).gt.subgrid_dim(x_dim_index).or.j_in(jj).lt.1.or.j_in(jj).gt.subgrid_dim(y_dim_index)) then
                                                    zgrid(ii,jj)=traveltime_subgrid(i_in(2),j_in(2),t,1,i_pollutant)
                                                elseif (use_subgrid_val(i_in(ii),j_in(jj),i_source).eq.outside_interpolation_region_index) then
                                                    zgrid(ii,jj)=traveltime_subgrid(i_in(2),j_in(2),t,1,i_pollutant)
                                                endif
                                            enddo
                                        enddo
                                        traveltime_subgrid(i,j,t,1,i_pollutant)=area_weighted_interpolation_function(xgrid,ygrid,zgrid,xdim,ydim,delta,xval,yval)

                                        do jj=1,3
                                            do ii=1,3
                                                zgrid(ii,jj)=traveltime_subgrid(i_in(ii),j_in(jj),t,2,i_pollutant)
                                                if (i_in(ii).lt.1.or.i_in(ii).gt.subgrid_dim(x_dim_index).or.j_in(jj).lt.1.or.j_in(jj).gt.subgrid_dim(y_dim_index)) then
                                                    zgrid(ii,jj)=traveltime_subgrid(i_in(2),j_in(2),t,2,i_pollutant)
                                                elseif (use_subgrid_val(i_in(ii),j_in(jj),i_source).eq.outside_interpolation_region_index) then
                                                    zgrid(ii,jj)=traveltime_subgrid(i_in(2),j_in(2),t,2,i_pollutant)
                                                endif
                                            enddo
                                        enddo
                                        traveltime_subgrid(i,j,t,2,i_pollutant)=area_weighted_interpolation_function(xgrid,ygrid,zgrid,xdim,ydim,delta,xval,yval)

                                    enddo
                                enddo

                            endif
                        enddo
                    enddo

                enddo
            endif
        enddo


        !Reset the use_subgrid values so chemistry and exposure happens everywhere but not outside the region
        use_subgrid(:,:,allsource_index)=.true.
        where (use_subgrid_val(:,:,allsource_index).eq.outside_region_index) use_subgrid(:,:,allsource_index)=.false.

    end subroutine uEMEP_interpolate_auto_subgrid

    subroutine uEMEP_region_mask
        !This routine defines use_subgrid_val=2 for regions outside the selected region. This can be used to control the interpolation routine
        !It also sets use_subgrid=.false. outside the region so that no calculations are made there either

        implicit none

        character(256) temp_name,temp_str,temp_str1
        integer, allocatable :: tile_municipality_subgrid(:,:,:)

        integer municipality_id
        real x_ssb,f_easting,ssb_dx,y_ssb,ssb_dy
        integer count
        integer*8 ssb_id
        integer i,j,k,i_range,j_range,i_range_interp,j_range_interp,i_tile,j_tile
        logical         exists
        integer unit_in
        integer index_val
        integer i_source
        character(256) region_number_str
        integer n_search
        parameter (n_search=5)
        character(16) search_str(n_search)
        real search_delta(n_search)
        integer temp_search
        integer i_tile_region,j_tile_region
        integer i_range_region,j_range_region
        integer :: total_grids=0
        integer :: source_count=0

        data search_str /'_1000m','_500m','_250m','_100m','_50m'/
        data search_delta /1000.,500.,250.,100.,50./

        write(unit_logfile,'(A)') ''
        write(unit_logfile,'(A)') '================================================================'
        write(unit_logfile,'(A)') 'Masking region (uEMEP_region_mask)'
        write(unit_logfile,'(A)') '================================================================'

        if (.not.allocated(tile_municipality_subgrid)) allocate (tile_municipality_subgrid(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),2))
        tile_municipality_subgrid=0

        i_source=allsource_index

        !Search file name to define the grid size
        ssb_dx=0.;ssb_dy=0.
        do k=1,n_search
            temp_search=index(filename_population(municipality_index),trim(adjustl(search_str(k))))
            if (temp_search.ne.0) then
                ssb_dx=search_delta(k)
                ssb_dy=search_delta(k)
                write(unit_logfile,'(i,A)') temp_search,' Reading municipality data with resolution '//trim(adjustl(search_str(k)))
            endif
        enddo

        if (ssb_dx.eq.0) then
            write(unit_logfile,'(A)') 'Cannot find a valid SSB grid size. Stopping. '//trim(filename_population(municipality_index))
            stop
        endif

        region_number_str=''
        write(region_number_str,*) region_index
        region_number_str=trim(region_number_str)//'_'

        !Read in SSB file containing gridded municipality ids
        f_easting=2.e6
        pathfilename_population(municipality_index)=trim(pathname_population(municipality_index))//trim(adjustl(region_number_str))//trim(filename_population(municipality_index))
        !Test existence of the heating filename. If does not exist then use default
        inquire(file=trim(pathfilename_population(municipality_index)),exist=exists)
        if (.not.exists) then
            write(unit_logfile,'(A,A)') ' ERROR: Masking region SSB file with municipality IDs does not exist: ', trim(pathfilename_population(municipality_index))
            stop
        endif
        temp_name=pathfilename_population(municipality_index)
        !Open the file for reading
        unit_in=20
        open(unit_in,file=temp_name,access='sequential',status='old',readonly)
        write(unit_logfile,'(a)') ' Opening SSB municipality file '//trim(temp_name)
        rewind(unit_in)
        !Read header SSBID0250M;kommunenum
        read(unit_in,'(A)') temp_str
        write(unit_logfile,'(A)') 'Header: '//trim(temp_str)
        count=0
        i_range=ceiling(ssb_dx/subgrid_delta(x_dim_index)/2.)+1
        j_range=ceiling(ssb_dy/subgrid_delta(y_dim_index)/2.)+1
        i_range_interp=ceiling(ssb_dx/subgrid_delta(x_dim_index)/2.+max_interpolation_subgrid_size/subgrid_delta(x_dim_index))+1
        j_range_interp=ceiling(ssb_dy/subgrid_delta(y_dim_index)/2.+max_interpolation_subgrid_size/subgrid_delta(y_dim_index))+1
        do while(.not.eof(unit_in))
            ssb_id=0;municipality_id=0
            !Read in file string
            read(unit_in,'(A)') temp_str
            index_val=index(temp_str,';',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:);if (index_val.gt.1) read(temp_str1,*) ssb_id
            read(temp_str,*) municipality_id
            !If this ssb municipality grid has the correct ID then find the target grid that matches it
            if (municipality_id.eq.region_id) then
                count=count+1
                !if (mod(count,100000).eq.0) write(*,*) count,ssb_id,municipality_id
                !Convert id to grid centre coordinates that are already in UTM33 for SSB data
                x_ssb=floor(ssb_id/10000000.)-f_easting+ssb_dx/2.
                y_ssb=mod(ssb_id,10000000)+ssb_dy/2.

                !Find the tile this ssb grid is in
                i_tile=1+floor((x_ssb-subgrid_min(x_dim_index))/subgrid_delta(x_dim_index)) !New definition
                j_tile=1+floor((y_ssb-subgrid_min(y_dim_index))/subgrid_delta(y_dim_index)) !New definition

                !write(*,*) i_tile,j_tile,i_range,j_range,subgrid_dim(x_dim_index),subgrid_dim(y_dim_index)

                do j=j_tile-j_range,j_tile+j_range
                    do i=i_tile-i_range,i_tile+i_range
                        !Make sure i and j are inside a valid range
                        if (i.ge.1.and.i.le.subgrid_dim(x_dim_index).and.j.ge.1.and.j.le.subgrid_dim(y_dim_index)) then
                            !Find the target subgrid within the ssb grid and the region
                            !write(*,*) ssb_dx,ssb_dy
                            !write(*,*) i,j,x_subgrid(i,j)-(x_ssb-ssb_dx/2.0),x_subgrid(i,j)-(x_ssb+ssb_dx/2.0),y_subgrid(i,j)-(y_ssb-ssb_dy/2.0),y_subgrid(i,j)-(y_ssb+ssb_dy/2.0)
                            if (x_subgrid(i,j).ge.x_ssb-ssb_dx/2.0.and.x_subgrid(i,j).lt.x_ssb+ssb_dx/2.0.and. &
                                y_subgrid(i,j).ge.y_ssb-ssb_dy/2.0.and.y_subgrid(i,j).lt.y_ssb+ssb_dy/2.0) then
                                !If there is a target subgrid within this range then allocate the mask value if it is not the correct region_id
                                tile_municipality_subgrid(i,j,1)=1
                            endif
                            !write(*,*) tile_municipality_subgrid(i,j,1)
                        endif
                    enddo
                enddo

                do j=j_tile-j_range_interp,j_tile+j_range_interp
                    do i=i_tile-i_range_interp,i_tile+i_range_interp
                        !Make sure i and j are inside a valid range
                        if (i.ge.1.and.i.le.subgrid_dim(x_dim_index).and.j.ge.1.and.j.le.subgrid_dim(y_dim_index)) then
                            !Find the target subgrid within the ssb grid and the region
                            if (x_subgrid(i,j).ge.x_ssb-ssb_dx/2.0-max_interpolation_subgrid_size.and.x_subgrid(i,j).le.x_ssb+ssb_dx/2.0+max_interpolation_subgrid_size.and. &
                                y_subgrid(i,j).ge.y_ssb-ssb_dy/2.0-max_interpolation_subgrid_size.and.y_subgrid(i,j).le.y_ssb+ssb_dy/2.0+max_interpolation_subgrid_size) then
                                !If there is a target subgrid within this range then allocate the mask value if it is not the correct region_id
                                tile_municipality_subgrid(i,j,2)=1
                            endif
                        endif
                    enddo
                enddo

                !Specify the emissions grids that will be used for regional calculations
                if (trace_emissions_from_in_region) then
                    source_count=0
                    do i_source=1,n_source_index
                        if (calculate_source(i_source)) then
                            i_tile_region=1+floor((x_ssb-emission_subgrid_min(x_dim_index,i_source))/emission_subgrid_delta(x_dim_index,i_source)) !New definition
                            j_tile_region=1+floor((y_ssb-emission_subgrid_min(y_dim_index,i_source))/emission_subgrid_delta(y_dim_index,i_source)) !New definition
                            i_range_region=ceiling(ssb_dx/emission_subgrid_delta(x_dim_index,i_source)/2.)+1
                            j_range_region=ceiling(ssb_dy/emission_subgrid_delta(y_dim_index,i_source)/2.)+1

                            !write(*,*) i_tile_region,j_tile_region,i_source,emission_subgrid_dim(x_dim_index,i_source),emission_subgrid_dim(y_dim_index,i_source)
                            do j=j_tile_region-j_range_region,j_tile_region+j_range_region
                                do i=i_tile_region-i_range_region,i_tile_region+i_range_region
                                    !Make sure i and j are inside a valid range
                                    if (i.ge.1.and.i.le.emission_subgrid_dim(x_dim_index,i_source).and.j.ge.1.and.j.le.emission_subgrid_dim(y_dim_index,i_source)) then
                                        !Find the target subgrid within the ssb grid and the region
                                        if (x_emission_subgrid(i,j,i_source).ge.x_ssb-ssb_dx/2.0.and.x_emission_subgrid(i,j,i_source).lt.x_ssb+ssb_dx/2.0.and. &
                                            y_emission_subgrid(i,j,i_source).ge.y_ssb-ssb_dy/2.0.and.y_emission_subgrid(i,j,i_source).lt.y_ssb+ssb_dy/2.0) then
                                            !Default is false
                                            use_subgrid_region(i,j,i_source)=.true.
                                            total_grids=total_grids+1
                                        endif
                                        !write(*,*) use_subgrid_region(i,j,i_source)
                                    endif
                                enddo
                            enddo
                            source_count=source_count+1
                        endif
                    enddo
                    !write(unit_logfile,'(A,i,i)')'Number of grids found (municipality, region) = ',sum(tile_municipality_subgrid(:,:,1)),total_grids/source_count

                endif

            endif
        enddo
        close(unit_in)


        do j=1,subgrid_dim(y_dim_index)
            do i=1,subgrid_dim(x_dim_index)
                !write(*,*) i,j,tile_municipality_subgrid(i,j,1),tile_municipality_subgrid(i,j,2)
                if (tile_municipality_subgrid(i,j,2).eq.0) then
                    use_subgrid_val(i,j,:)=outside_interpolation_region_index
                    use_subgrid(i,j,:)=.false.
                endif
                if (tile_municipality_subgrid(i,j,1).eq.0) then
                    use_subgrid_val(i,j,:)=outside_region_index
                    !use_subgrid(i,j,:)=.false.
                endif
            enddo
        enddo


        if (allocated(tile_municipality_subgrid)) deallocate (tile_municipality_subgrid)

    end subroutine uEMEP_region_mask

end module auto_subgrid

