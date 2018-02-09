!uEMEP_set_filenames.f90
    
    subroutine uEMEP_set_filenames

    use uEMEP_definitions

    implicit none
    
    !Set all the input data here that is to be read
    !pathname_rl(1)='C:\uEMEP\Data\Road_data\'
    !filename_rl(1)='Road_data_ascii_ERFK_large_Oslo.txt'
    !filename_rl(1)='Road_data_ascii_ERFK_Bergen.txt'
    !filename_rl(1)='Road_data_ascii_ERFK_Stavanger.txt'
    !pathname_ship(1)='C:\uEMEP\Data\Shipping_data\'
    !filename_ship(1)='Bruce_201501.csv'
    !pathname_agriculture(1)='C:\uEMEP\Data\Nederland\'
    !filename_agriculture(1)='NH3_agriculture_animal_housings_2014.csv'
    !pathname_agriculture(2)='C:\uEMEP\Data\Nederland\'
    !filename_agriculture(2)='NH3_agriculture_other_2014_1000x1000.csv'
    !pathname_heating(dwelling_index)='C:\uEMEP\Data\SSB\NOR0250M_DWE_2017\'
    !filename_heating(dwelling_index)='SSB0250M_DWE_2017.csv'
    !pathname_heating(population_index)='C:\uEMEP\Data\SSB\r250m_web_2017\'
    !filename_heating(population_index)='r250m_web_2017.csv'
    
    !filename_receptor='Station_positions_Norway.txt'
    !pathname_receptor='C:\uEMEP\Data\Airbase\'

    
    !EMEP input files
    !if (hourly_calculations) then
    !    pathname_EMEP='C:\uEMEP\EMEP_data\Hourly\'
    !    filename_EMEP='Save_Base_hour.nc'
    ! endif
    !if (annual_calculations) then
     !   pathname_EMEP='C:\uEMEP\EMEP_data\nh3\'
     !   filename_EMEP='uEMEP_nh3_RIVM_EMIS_fullrun.nc'
    !endif
    
    !Set pathname for all gridded data to be saved
    pathname_grid(:)=pathname_output_grid

    !Set filenames for all gridded data to be saved
    do i=1,n_source_index
        filename_grid(proxy_emission_file_index(i))=trim('proxy_emission_subgrid')//'_'//trim(source_file_str(i))
        filename_grid(emission_file_index(i))=trim('emission_subgrid')//'_'//trim(source_file_str(i))
        filename_grid(proxy_file_index(i))=trim('proxy_subgrid')//'_'//trim(source_file_str(i))
        filename_grid(proxy_integral_file_index(i))=trim('proxy_integral_subgrid')//'_'//trim(source_file_str(i))
        filename_grid(emep_subgrid_file_index(i))=trim('EMEP_subgrid')//'_'//trim(source_file_str(i))
        filename_grid(emep_subgrid_nonlocal_file_index(i))=trim('EMEP_nonlocal_subgrid')//'_'//trim(source_file_str(i))
        filename_grid(emep_subgrid_local_file_index(i))=trim('EMEP_local_subgrid')//'_'//trim(source_file_str(i))
        filename_grid(emep_subgrid_frac_file_index(i))=trim('EMEP_frac_subgrid')//'_'//trim(source_file_str(i))
        filename_grid(subgrid_local_file_index(i))=trim('local_subgrid')//'_'//trim(source_file_str(i))
        filename_grid(subgrid_total_file_index(i))=trim('total_subgrid')//'_'//trim(source_file_str(i))
        filename_grid(use_subgrid_file_index(i))=trim('use_subgrid')//'_'//trim(source_file_str(i))
        filename_grid(emep_emission_subgrid_file_index(i))=trim('EMEP_emission_subgrid')//'_'//trim(source_file_str(i))       
    enddo
    
    filename_grid(population_file_index(dwelling_index))=trim('dwelling_subgrid')
    filename_grid(population_file_index(population_index))=trim('population_subgrid')
    filename_grid(population_file_index(school_index))=trim('school_subgrid')
    filename_grid(population_file_index(establishment_index))=trim('establishment_subgrid')
    filename_grid(population_file_index(kindergaten_index))=trim('kindergaten_subgrid')
    filename_grid(population_file_index(home_index))=trim('home_subgrid')

    !Meteo files
    filename_grid(subgrid_ugrid_file_index)='ugrid_subgrid'
    filename_grid(subgrid_vgrid_file_index)='vgrid_subgrid'
    filename_grid(subgrid_hmix_file_index)='hmix_subgrid'
    filename_grid(subgrid_kz_file_index)='kz_subgrid'
    filename_grid(subgrid_logz0_file_index)='logz0_subgrid'
    filename_grid(subgrid_invL_file_index)='invL_subgrid'
    filename_grid(subgrid_FFgrid_file_index)='FFgrid_subgrid'
    filename_grid(subgrid_FF10_file_index)='FF10_subgrid'
    filename_grid(subgrid_invFFgrid_file_index)='invFFgrid_subgrid'
    filename_grid(subgrid_invFF10_file_index)='invFF10_subgrid'
    filename_grid(subgrid_ustar_file_index)='ustar_subgrid'
    filename_grid(subgrid_J_file_index)='J_subgrid'
    
    end subroutine uEMEP_set_filenames