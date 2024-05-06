module set_filenames

    use uEMEP_definitions !, only: pathname_grid, filename_grid, source_file_str
    use uemep_configuration !, only: pathname_output_grid, save_netcdf_fraction_as_contribution_flag

    implicit none
    private

    public :: uEMEP_set_filenames

contains

    subroutine uEMEP_set_filenames()
        !! Set filenames for all gridded data to be saved
        !! These are the names now given in the netcdf files
        
        ! Local variables
        integer :: i

        ! Set pathname for all gridded data to be saved
        pathname_grid(:) = pathname_output_grid
        
        do i = 1, n_source_index
            filename_grid(proxy_emission_file_index(i)) = trim('proxy_emission_subgrid')//'_'//trim(source_file_str(i))
            filename_grid(emission_file_index(i)) = trim('emission_subgrid')//'_'//trim(source_file_str(i))
            filename_grid(proxy_file_index(i)) = trim('proxy_subgrid')//'_'//trim(source_file_str(i))
            filename_grid(proxy_integral_file_index(i)) = trim('proxy_integral_subgrid')//'_'//trim(source_file_str(i))
            filename_grid(emep_subgrid_file_index(i)) = trim('EMEP_subgrid')//'_'//trim(source_file_str(i))
            filename_grid(emep_subgrid_nonlocal_file_index(i)) = trim('EMEP_nonlocal_subgrid')//'_'//trim(source_file_str(i))
            filename_grid(emep_subgrid_local_file_index(i)) = trim('EMEP_local_subgrid')//'_'//trim(source_file_str(i))
            filename_grid(emep_additional_subgrid_nonlocal_file_index(i)) = trim('EMEP_additional_nonlocal_subgrid')//'_'//trim(source_file_str(i))
            filename_grid(emep_additional_subgrid_local_file_index(i)) = trim('EMEP_additional_local_subgrid')//'_'//trim(source_file_str(i))
            filename_grid(emep_subgrid_frac_file_index(i)) = trim('EMEP_frac_subgrid')//'_'//trim(source_file_str(i))
            filename_grid(subgrid_local_file_index(i)) = trim('local_subgrid')//'_'//trim(source_file_str(i))
            filename_grid(subgrid_total_file_index(i)) = trim('total_subgrid')//'_'//trim(source_file_str(i))
            filename_grid(use_subgrid_file_index(i)) = trim('use_subgrid')//'_'//trim(source_file_str(i))
            filename_grid(emep_emission_subgrid_file_index(i)) = trim('EMEP_emission_subgrid')//'_'//trim(source_file_str(i))
        end do
        
        ! Alternative set of names for outputs to netcdf
        do i=1,n_source_index
            filename_grid(proxy_emission_file_index(i)) = trim('proxy_emission')//'_'//trim(source_file_str(i))
            filename_grid(emission_file_index(i)) = trim('emission')//'_'//trim(source_file_str(i))
            filename_grid(proxy_file_index(i)) = trim('proxy')//'_'//trim(source_file_str(i))
            filename_grid(proxy_integral_file_index(i)) = trim('proxy_integral')//'_'//trim(source_file_str(i))
            filename_grid(emep_subgrid_file_index(i)) = trim('EMEP')//'_'//trim(source_file_str(i))
            if (save_netcdf_fraction_as_contribution_flag) then
                filename_grid(emep_subgrid_nonlocal_file_index(i)) = trim('EMEP_nonlocal_contribution')
                filename_grid(emep_subgrid_local_file_index(i)) = trim('EMEP_local_contribution')//'_'//trim(source_file_str(i))
                filename_grid(emep_additional_subgrid_nonlocal_file_index(i)) = trim('EMEP_additional_nonlocal_contribution')
                filename_grid(emep_additional_subgrid_local_file_index(i)) = trim('EMEP_additional_local_contribution')//'_'//trim(source_file_str(i))
                filename_grid(emep_subgrid_frac_file_index(i)) = trim('EMEP_contribution')//'_'//trim(source_file_str(i))
                filename_grid(subgrid_local_file_index(i)) = trim('local_contribution')//'_'//trim(source_file_str(i))
            else
                filename_grid(emep_subgrid_nonlocal_file_index(i)) = trim('EMEP_nonlocal_fraction')
                filename_grid(emep_subgrid_local_file_index(i)) = trim('EMEP_local_fraction')//'_'//trim(source_file_str(i))
                filename_grid(emep_additional_subgrid_nonlocal_file_index(i)) = trim('EMEP_additional_nonlocal_fraction')
                filename_grid(emep_additional_subgrid_local_file_index(i)) = trim('EMEP_additional_local_fraction')//'_'//trim(source_file_str(i))
                filename_grid(emep_subgrid_frac_file_index(i)) = trim('EMEP_fraction')//'_'//trim(source_file_str(i))
                filename_grid(subgrid_local_file_index(i)) = trim('local_fraction')//'_'//trim(source_file_str(i))
            end if
            filename_grid(subgrid_total_file_index(i)) = trim('total')//'_'//trim(source_file_str(i))
            filename_grid(use_subgrid_file_index(i)) = trim('use_subgrid')//'_'//trim(source_file_str(i))
            filename_grid(emep_emission_subgrid_file_index(i)) = trim('EMEP_emission')//'_'//trim(source_file_str(i))
        end do

        filename_grid(population_file_index(dwelling_index)) = trim('dwelling')
        filename_grid(population_file_index(population_index)) = trim('population')
        filename_grid(population_file_index(school_index)) = trim('school')
        filename_grid(population_file_index(establishment_index)) = trim('establishment')
        filename_grid(population_file_index(kindergaten_index)) = trim('kindergaten')
        filename_grid(population_file_index(home_index)) = trim('home')

        ! Meteo file names
        filename_grid(subgrid_ugrid_file_index) = 'xgrid_wind'
        filename_grid(subgrid_vgrid_file_index) = 'ygrid_wind'
        filename_grid(subgrid_u10_file_index) = 'x10_wind'
        filename_grid(subgrid_v10_file_index) = 'y10_wind'
        filename_grid(subgrid_hmix_file_index) = 'hmix'
        filename_grid(subgrid_kz_file_index) = 'kz'
        filename_grid(subgrid_logz0_file_index) = 'logz0'
        filename_grid(subgrid_invL_file_index) = 'inv_L'
        filename_grid(subgrid_FFgrid_file_index) = 'wind_speed_grid'
        filename_grid(subgrid_DDgrid_file_index) = 'wind_direction_grid'
        filename_grid(subgrid_FF10_file_index) = 'wind_speed_10m'
        filename_grid(subgrid_DD10_file_index) = 'wind_direction_10m'
        filename_grid(subgrid_invFFgrid_file_index) = 'inv_FFgrid'
        filename_grid(subgrid_invFF10_file_index) = 'inv_FF10'
        filename_grid(subgrid_ustar_file_index) = 'ustar'
        filename_grid(subgrid_J_file_index) = 'J_photo'
        filename_grid(subgrid_t2m_file_index) = 'air_temperature_2m'

    end subroutine uEMEP_set_filenames

end module set_filenames
