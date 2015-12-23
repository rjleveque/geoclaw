
    real(kind=8) function sea_level_fcn(x,y,t) result(vsea_level)

        ! Sample sea_level_fcn that defines a varying sea level.
        ! This will be called only if variabel_sea_level == .true. in
        ! geoclaw_module.f90 (should be a setrun parameter eventually).

        use geoclaw_module, only: sea_level

        implicit none

        ! Arguments
        real(kind=8), intent(in) :: x,y,t

        real(kind=8) :: x1,x2,y1,y2,dz

        x1 = -124.4d0
        x2 = -123.6d0
        y1 = 46.8d0
        y2 = 47.1d0

        if ((t < 1.d0) .or. (x <= x1) .or. (x >= x2) &
                       .or. (y <= y1) .or. (y >= y2)) &
          then
            ! use default sea_level parameter
            vsea_level = sea_level
          else
            ! approximate subsidence of L1 near Gray's Harbor:
            dz = -3.d0 + 1.2d0*((x+123.9d0)/0.3d0)**2
            vsea_level = sea_level + dz
          endif

        end function sea_level_fcn
            
        

