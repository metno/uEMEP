!Area_weighted_interpolation_function
      
function area_weighted_interpolation_function(xgrid,ygrid,zgrid,xdim,ydim,delta,xval,yval)  
    
    implicit none
    
    integer, intent(in) :: xdim,ydim
    real, intent(in) :: delta(2)
    real, intent(in) :: xgrid(xdim,ydim),ygrid(xdim,ydim),zgrid(xdim,ydim)
    real, intent(in) :: xval,yval
    real zval
    real area_weighted_interpolation_function
    real sum_weight
    
    real weighting
    integer i,j,ii,jj
    real  xpos_area_max,xpos_area_min,ypos_area_max,ypos_area_min
    real  xpos_max,xpos_min,ypos_max,ypos_min
        
    !Find grid index for position val
    i=1+floor((xval-xgrid(1,1))/delta(1)+0.5)
    j=1+floor((yval-ygrid(1,1))/delta(2)+0.5)     

    if (i.lt.1.or.j.lt.1.or.i.gt.xdim.or.j.gt.ydim) then
        write(*,'(A,4i6)') 'Interpolation out of range. Stopping. (i,j,xdim,ydim)',i,j,xdim,ydim
        write(*,'(4f12.2)') xval,yval,xgrid(1,1),ygrid(1,1)
        stop
    else
        
        xpos_area_max=xval+delta(1)/2.
        xpos_area_min=xval-delta(1)/2.
        ypos_area_max=yval+delta(2)/2.
        ypos_area_min=yval-delta(2)/2.
                  
        zval=0.
        sum_weight=0.
        
        do jj=j-1,j+1
        do ii=i-1,i+1         
            
            xpos_min=max(xpos_area_min,xgrid(ii,jj)-delta(1)/2.)
            xpos_max=min(xpos_area_max,xgrid(ii,jj)+delta(1)/2.)
            ypos_min=max(ypos_area_min,ygrid(ii,jj)-delta(2)/2.)
            ypos_max=min(ypos_area_max,ygrid(ii,jj)+delta(2)/2.)
               
            if (xpos_max.gt.xpos_min.and.ypos_max.gt.ypos_min) then
                weighting=(ypos_max-ypos_min)*(xpos_max-xpos_min)/delta(1)/delta(2)
            else
                weighting=0.
            endif                

            zval=zval+zgrid(ii,jj)*weighting

            !write(*,'(4i6,f12.4,f12.4)') i,j,ii,jj,weighting,zgrid(ii,jj)
            sum_weight=sum_weight+weighting
        enddo
        enddo
            !write(*,'(f12.4,f12.4)') zval,sum_weight
          
    endif
    
    area_weighted_interpolation_function=zval
        
    
end function area_weighted_interpolation_function
    
