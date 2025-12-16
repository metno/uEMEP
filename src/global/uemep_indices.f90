module uemep_indices

    use index_tables

    implicit none

    type(index_table_t) :: source_idx
    type(index_table_t) :: nc_idx

    ! Legacy indices
    integer :: allsource_index, traffic_index, shipping_index, heating_index, agriculture_index, &
        industry_index, publicpower_index, fugitive_index, solvents_index, aviation_index, &
        offroad_index, waste_index, livestock_index, other_index, traffic_exhaust_index, traffic_nonexhaust_index, &
        n_source_index

    integer :: allsource_nc_index, traffic_nc_index, shipping_nc_index, heating_nc_index, agriculture_nc_index, &
        industry_nc_index, publicpower_nc_index, fugitive_nc_index, solvents_nc_index, aviation_nc_index, &
        offroad_nc_index, waste_nc_index, livestock_nc_index, other_nc_index, traffic_exhaust_nc_index, &
        traffic_nonexhaust_nc_index, traffic_gasoline_nc_index, traffic_diesel_nc_index, traffic_gas_nc_index, &
        publicpower_point_nc_index, publicpower_area_nc_index, extrasource_nc_index, n_source_nc_index
    

contains

    subroutine set_indices()

        call source_idx%init()
        call source_idx%add("allsource")
        call source_idx%add("traffic")
        call source_idx%add("shipping")
        call source_idx%add("heating")
        call source_idx%add("agriculture")
        call source_idx%add("industry")
        call source_idx%add("publicpower")
        call source_idx%add("fugitive")
        call source_idx%add("solvents")
        call source_idx%add("aviation")
        call source_idx%add("offroad")
        call source_idx%add("waste")
        call source_idx%add("livestock")
        call source_idx%add("other")
        call source_idx%add("traffic_exhaust")
        call source_idx%add("traffic_nonexhaust")

        call nc_idx%init()
        call nc_idx%add("allsource")
        call nc_idx%add("traffic")
        call nc_idx%add("shipping")
        call nc_idx%add("heating")
        call nc_idx%add("agriculture")
        call nc_idx%add("industry")
        call nc_idx%add("publicpower")
        call nc_idx%add("fugitive")
        call nc_idx%add("solvents")
        call nc_idx%add("aviation")
        call nc_idx%add("offroad")
        call nc_idx%add("waste")
        call nc_idx%add("livestock")
        call nc_idx%add("other")
        call nc_idx%add("traffic_exhaust")
        call nc_idx%add("traffic_nonexhaust")
        call nc_idx%add("traffic_gasoline")
        call nc_idx%add("traffic_diesel")
        call nc_idx%add("traffic_gas")
        call nc_idx%add("publicpower_point")
        call nc_idx%add("publicpower_area")
        call nc_idx%add("extrasource")

        call set_legacy_indices()

    end subroutine set_indices

    subroutine set_legacy_indices()
        
        allsource_index = source_idx%get("allsource")
        traffic_index = source_idx%get("traffic")
        shipping_index = source_idx%get("shipping")
        heating_index = source_idx%get("heating")
        agriculture_index = source_idx%get("agriculture")
        industry_index = source_idx%get("industry")
        publicpower_index = source_idx%get("publicpower")
        fugitive_index = source_idx%get("fugitive")
        solvents_index = source_idx%get("solvents")
        aviation_index = source_idx%get("aviation")
        offroad_index = source_idx%get("offroad")
        waste_index = source_idx%get("waste_index")
        livestock_index = source_idx%get("livestock")
        other_index = source_idx%get("other")
        traffic_exhaust_index = source_idx%get("traffic_exhaust")
        traffic_nonexhaust_index = source_idx%get("traffic_nonexhaust")

        n_source_index = source_idx%n

        allsource_nc_index = nc_idx%get("allsource")
        traffic_nc_index = nc_idx%get("traffic")
        shipping_nc_index = nc_idx%get("shipping")
        heating_nc_index = nc_idx%get("heating")
        agriculture_nc_index = nc_idx%get("agriculture")
        industry_nc_index = nc_idx%get("industry")
        publicpower_nc_index = nc_idx%get("publicpower")
        fugitive_nc_index = nc_idx%get("fugitive")
        solvents_nc_index = nc_idx%get("solvents")
        aviation_nc_index = nc_idx%get("aviation")
        offroad_nc_index = nc_idx%get("offroad")
        waste_nc_index = nc_idx%get("waste")
        livestock_nc_index = nc_idx%get("livestock")
        other_nc_index = nc_idx%get("other")
        traffic_exhaust_nc_index = nc_idx%get("traffic_exhaust")
        traffic_nonexhaust_nc_index = nc_idx%get("traffic_nonexhaust")
        traffic_gasoline_nc_index = nc_idx%get("traffic_gasoline")
        traffic_diesel_nc_index = nc_idx%get("traffic_diesel")
        traffic_gas_nc_index = nc_idx%get("traffic_gas")
        publicpower_point_nc_index = nc_idx%get("publicpower_point")
        publicpower_area_nc_index = nc_idx%get("publicpower_area")
        extrasource_nc_index = nc_idx%get("extrasource")

        n_source_nc_index = nc_idx%n

    end subroutine set_legacy_indices

end module uemep_indices