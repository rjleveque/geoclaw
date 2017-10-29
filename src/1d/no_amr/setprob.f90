subroutine setprob

   !set gauges and other parameters for 1D geoclaw
   ! for other problem set-up features, copy to your directory and modify

   use gauges_module, only: set_gauges
   use geoclaw_module, only: set_geo, coordinate_system
   use grid_module, only: read_grid

   implicit none

   call set_gauges()
   call set_geo()

   if (coordinate_system == 2) then
        ! mapped grid, assume file grid.data exists...
        read_grid('grid.data')
      endif

end subroutine setprob
