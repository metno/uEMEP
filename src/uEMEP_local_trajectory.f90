module local_trajectory

    implicit none
    private

    public :: uEMEP_calculate_all_trajectory, uEMEP_minimum_distance_trajectory_fast

contains

!uEMEP_local_trajectory

    !subroutine uEMEP_local_trajectory(i_rec,j_rec,i_emis,j_emis,t,traj_max_index,dr_traj,i_source,x_loc,y_loc,valid_traj)
    subroutine uEMEP_local_trajectory(x_r,y_r,x_emis,y_emis,t,traj_max_index,dr_traj,x_loc,y_loc,valid_traj)

        use uEMEP_definitions

        implicit none

        integer k
        real x_loc,y_loc

        !integer i_rec,j_rec,i_emis,j_emis
        integer t,traj_max_index
        !integer i_source
        real dr_traj
        real x_r,y_r,x_emis,y_emis
        real x_traj(traj_max_index),y_traj(traj_max_index)
        real distance_traj(traj_max_index),distance_intercept_traj(traj_max_index)
        real x_intercept_traj(traj_max_index),y_intercept_traj(traj_max_index),frac_length_traj(traj_max_index)
        logical exit_traj,valid_traj
        integer i_integral,j_integral
        !real x_r,y_r

        k=1
        !Set the initial trajectory position to the emission source
        !x_traj(k)=x_emission_subgrid(i_emis,j_emis,i_source)
        !y_traj(k)=y_emission_subgrid(i_emis,j_emis,i_source)
        x_traj(k)=x_emis
        y_traj(k)=y_emis

        !Set the position of the receptor
        !x_r=x_subgrid(i_rec,j_rec)
        !y_r=y_subgrid(i_rec,j_rec)

        !Set the distances for the initial emission grid
        !distance_traj(k)=sqrt((x_traj(k)-x_r)*(x_traj(k)-x_r)+(y_traj(k)-y_r)*(y_traj(k)-y_r))
        !distance_intercept_traj(k)=distance_traj(k)
        distance_traj(k)=(x_traj(k)-x_r)*(x_traj(k)-x_r)+(y_traj(k)-y_r)*(y_traj(k)-y_r)
        distance_intercept_traj(k)=distance_traj(k)*distance_traj(k)
        y_loc=0.
        x_loc=0.
        exit_traj=.false.
        valid_traj=.false.

        do while (.not.exit_traj.and.k.lt.traj_max_index)

            !Find what meteorological grid the trajectory point is in
            !i_integral=1+floor((x_traj(k)-integral_subgrid_min(x_dim_index))/(integral_subgrid_max(x_dim_index)-integral_subgrid_min(x_dim_index))*integral_subgrid_dim(x_dim_index))
            !j_integral=1+floor((y_traj(k)-integral_subgrid_min(y_dim_index))/(integral_subgrid_max(y_dim_index)-integral_subgrid_min(y_dim_index))*integral_subgrid_dim(y_dim_index))
            i_integral=1+floor((x_traj(k)-integral_subgrid_min(x_dim_index))/integral_subgrid_delta(x_dim_index))
            j_integral=1+floor((y_traj(k)-integral_subgrid_min(y_dim_index))/integral_subgrid_delta(y_dim_index))


            !Test to see if it is still in the grid
            if (i_integral.ge.1.and.i_integral.le.integral_subgrid_dim(x_dim_index).and.j_integral.ge.1.and.j_integral.le.integral_subgrid_dim(y_dim_index)) then

                !Account for at grid value
                if (distance_traj(k).eq.0) then
                    exit_traj=.true.
                    y_loc=0.
                    x_loc=0.
                    valid_traj=.true.
                endif

                k=k+1

                x_traj(k)=x_traj(k-1)+dr_traj*meteo_subgrid(i_integral,j_integral,t,cos_subgrid_index)
                y_traj(k)=y_traj(k-1)+dr_traj*meteo_subgrid(i_integral,j_integral,t,sin_subgrid_index)

                !distance_traj(k)=sqrt((x_traj(k)-x_r)*(x_traj(k)-x_r)+(y_traj(k)-y_r)*(y_traj(k)-y_r))
                distance_traj(k)=(x_traj(k)-x_r)*(x_traj(k)-x_r)+(y_traj(k)-y_r)*(y_traj(k)-y_r)

                call DISTRL_SQR(x_r,y_r,x_traj(k-1),y_traj(k-1),x_traj(k),y_traj(k),x_intercept_traj(k),y_intercept_traj(k),distance_intercept_traj(k),frac_length_traj(k))

                if (distance_intercept_traj(k).lt.distance_traj(k).and.distance_intercept_traj(k).le.distance_traj(k-1)) then
                    exit_traj=.true.
                    y_loc=distance_intercept_traj(k)
                    x_loc=dr_traj*(k-2)+frac_length_traj(k)*dr_traj
                    valid_traj=.true.
                    !elseif (distance_intercept_traj(k).gt.distance_intercept_traj(k-1)) then
                    !    valid_traj=.false.
                    !    exit_traj=.true.
                    !    k=k-1
                endif


            else
                exit_traj=.true.
                valid_traj=.false.
            endif

        enddo

        y_loc=sqrt(y_loc)

    end subroutine uEMEP_local_trajectory


    subroutine uEMEP_calculate_all_trajectory(x_emis,y_emis,t,traj_max_index,dr_traj,x_traj,y_traj)

        use uEMEP_definitions

        implicit none

        integer k
        !integer i_rec,j_rec,i_emis,j_emis
        integer t,traj_max_index
        !integer i_source
        real dr_traj
        real x_emis,y_emis
        real x_traj(traj_max_index),y_traj(traj_max_index)
        integer i_integral,j_integral
        logical exit_traj

        k=1
        !Set the initial trajectory position to the emission source
        x_traj(k)=x_emis
        y_traj(k)=y_emis

        exit_traj=.false.

        do while (.not.exit_traj.and.k.lt.traj_max_index)

            !Find what meteorological grid the trajectory point is in
            i_integral=1+floor((x_traj(k)-integral_subgrid_min(x_dim_index))/integral_subgrid_delta(x_dim_index))
            j_integral=1+floor((y_traj(k)-integral_subgrid_min(y_dim_index))/integral_subgrid_delta(y_dim_index))

            !Test to see if it is still in the grid
            if (i_integral.ge.1.and.i_integral.le.integral_subgrid_dim(x_dim_index).and.j_integral.ge.1.and.j_integral.le.integral_subgrid_dim(y_dim_index)) then

                k=k+1
                x_traj(k)=x_traj(k-1)+dr_traj*meteo_subgrid(i_integral,j_integral,t,cos_subgrid_index)
                y_traj(k)=y_traj(k-1)+dr_traj*meteo_subgrid(i_integral,j_integral,t,sin_subgrid_index)

            else
                exit_traj=.true.
            endif

        enddo

    end subroutine uEMEP_calculate_all_trajectory


    subroutine uEMEP_minimum_distance_trajectory(x_r,y_r,x_emis,y_emis,traj_max_index_in,dr_traj,x_traj,y_traj,x_loc,y_loc,valid_traj)

        use uEMEP_definitions

        implicit none

        real, intent(out) ::  x_loc,y_loc
        logical, intent(out) ::  valid_traj
        real, intent(in) ::  dr_traj
        real, intent(in) ::  x_r,y_r,x_emis,y_emis
        integer, intent(in) :: traj_max_index_in
        !real, intent(in) ::  x_traj(traj_max_index_in),y_traj(traj_max_index_in)
        real, intent(in) ::  x_traj(*),y_traj(*)

        integer k
        real distance_traj(traj_max_index_in),distance_intercept_traj(traj_max_index_in)
        real x_intercept_traj(traj_max_index_in),y_intercept_traj(traj_max_index_in),frac_length_traj(traj_max_index_in)

        real distance_intercept_min

        !write(*,*) traj_max_index_in

        k=1
        !Set the initial trajectory position to the emission source
        !x_traj(k)=x_emission_subgrid(i_emis,j_emis,i_source)
        !y_traj(k)=y_emission_subgrid(i_emis,j_emis,i_source)
        !x_traj(k)=x_emis
        !y_traj(k)=y_emis

        !Set the position of the receptor
        !x_r=x_subgrid(i_rec,j_rec)
        !y_r=y_subgrid(i_rec,j_rec)
        return

        !Set the distances for the initial emission grid
        !distance_traj(k)=sqrt((x_traj(k)-x_r)*(x_traj(k)-x_r)+(y_traj(k)-y_r)*(y_traj(k)-y_r))
        distance_traj(k)=(x_traj(k)-x_r)*(x_traj(k)-x_r)+(y_traj(k)-y_r)*(y_traj(k)-y_r)

        !Leave the routine because the receptor is the same as the emission grid
        if (distance_traj(k).eq.0) then
            y_loc=0.
            x_loc=0.
            valid_traj=.true.
            return
        endif

        !distance_intercept_traj(k)=distance_traj(k)
        y_loc=0.
        x_loc=0.
        valid_traj=.false.

        distance_intercept_min=distance_traj(k)

        do k=2,traj_max_index_in

            if (x_traj(k).ne.NODATA_value) then
                call DISTRL_SQR(x_r,y_r,x_traj(k-1),y_traj(k-1),x_traj(k),y_traj(k),x_intercept_traj(k),y_intercept_traj(k),distance_intercept_traj(k),frac_length_traj(k))

                if (distance_intercept_traj(k).lt.distance_intercept_min) then
                    distance_intercept_min=distance_intercept_traj(k)
                    y_loc=distance_intercept_traj(k)
                    x_loc=dr_traj*(k-2)+frac_length_traj(k)*dr_traj
                    valid_traj=.true.
                endif
            endif
        enddo

        !Remove most of the results because they are upwind
        if (x_loc.eq.0.and.y_loc.gt.dr_traj) then
            valid_traj=.false.
        endif

        y_loc=sqrt(y_loc)

    end subroutine uEMEP_minimum_distance_trajectory

    subroutine uEMEP_minimum_distance_trajectory_fast(x_r,y_r,x_emis,y_emis,traj_max_index_in,dr_traj,x_traj,y_traj,x_loc,y_loc,valid_traj)

        use uEMEP_definitions

        implicit none

        real, intent(out) ::  x_loc,y_loc
        logical, intent(out) ::  valid_traj
        real, intent(in) ::  dr_traj
        real, intent(in) ::  x_r,y_r,x_emis,y_emis
        integer, intent(in) :: traj_max_index_in
        !real, intent(in) ::  x_traj(traj_max_index_in),y_traj(traj_max_index_in)
        real, intent(in) ::  x_traj(*),y_traj(*)

        integer k
        real distance_traj,distance_intercept_traj
        real x_intercept_traj,y_intercept_traj,frac_length_traj

        real distance_intercept_min

        k=1
        !return

        !Set the distances for the initial emission grid
        !distance_traj(k)=sqrt((x_traj(k)-x_r)*(x_traj(k)-x_r)+(y_traj(k)-y_r)*(y_traj(k)-y_r))
        distance_traj=(x_traj(k)-x_r)*(x_traj(k)-x_r)+(y_traj(k)-y_r)*(y_traj(k)-y_r)

        !Leave the routine because the receptor is the same as the emission grid
        if (distance_traj.eq.0) then
            y_loc=0.
            x_loc=0.
            valid_traj=.true.
            return
        endif

        !distance_intercept_traj(k)=distance_traj(k)
        y_loc=0.
        x_loc=0.
        valid_traj=.false.

        distance_intercept_min=distance_traj

        do k=2,traj_max_index_in
            if (x_traj(k).ne.NODATA_value) then
                call DISTRL_SQR(x_r,y_r,x_traj(k-1),y_traj(k-1),x_traj(k),y_traj(k),x_intercept_traj,y_intercept_traj,distance_intercept_traj,frac_length_traj)

                if (distance_intercept_traj.lt.distance_intercept_min) then
                    distance_intercept_min=distance_intercept_traj
                    y_loc=distance_intercept_traj
                    x_loc=dr_traj*(k-2)+frac_length_traj*dr_traj
                    valid_traj=.true.
                endif
            endif
        enddo

        !Remove most of the results because they are upwind
        if (x_loc.eq.0.and.y_loc.gt.dr_traj) then
            valid_traj=.false.
        endif

        y_loc=sqrt(y_loc)

    end subroutine uEMEP_minimum_distance_trajectory_fast

end module local_trajectory

