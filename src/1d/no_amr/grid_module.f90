module grid_module

    use geoclaw_module, only: bouss
    
    implicit none
    save

    integer, parameter :: mx_grid_max=100000
    integer :: mx_grid
    logical :: radial
    real(kind=8), dimension(0:mx_grid_max+1) :: dxm, dxc, cm, cp, c0

    real(kind=8), dimension(0:mx_grid_max+1) ::  xgrid, zgrid, xcell

contains

subroutine read_grid(fname)

   implicit none
   character(len=*), intent(in) :: fname
   integer :: i,j
   real(kind=8) :: rim,rip,ric,c0i,cmi,cpi,r

   open(unit=58, file=fname, status='old',form='formatted')
   read(58,*) mx_grid
   if (mx_grid+1 > mx_grid_max) then
      write(6,*) '*** too many topo values'
      write(6,*) '*** increase mx_grid_max in grid_module.f90'
      stop
      endif

   ! read in grid cell edges and topo value at each edge:
   do i=1,mx_grid+1
      read(58,*) xgrid(i),zgrid(i)
      enddo

      
   ! extend to ghost cells, to have same width as first interior cell:
   xgrid(0) = 2.d0*xgrid(1) - xgrid(2)
   xgrid(mx_grid+2) = 2.d0*xgrid(mx_grid+1) - xgrid(mx_grid)
   zgrid(0) = zgrid(1)
   zgrid(mx_grid+2) = zgrid(mx_grid+1)

   do i=0,mx_grid+1
      ! cell centers based on grid edges:
      xcell(i) = 0.5d0*(xgrid(i) + xgrid(i+1))
      enddo

   do i=1,mx_grid+1
      dxm(i) = xcell(i) - xcell(i-1)
      enddo


   do i=1,mx_grid
      ! define coefficients cm(i), cp(i), c0(i) so that 
      !   q_{xx} \approx cm(i)*Q(i-1) - c0(i)*Q(i) + cp(i)*Q(i+1)
      ! Note that approximation is not centered at xcell(i), so only 1st order
      ! in general, but should be 2nd order if grid is smoothly varying (?).
      
      dxc(i) = (xcell(i+1) - xcell(i-1))  ! = 2*dx for uniform grid!
      cm(i) = 2.d0 / (dxm(i)*dxc(i))
      cp(i) = 2.d0 / (dxm(i+1)*dxc(i))
      c0(i) = cm(i) + cp(i)

      if (radial .and. bouss) then
         ! include factors for ((1/r)*(r*q)_r)_r 
         ! using form q_{rr} + (1/r)q_r - (1/r**2)q
         r = xcell(i) - xgrid(1)  ! assuming radial about left edge!
         cm(i) = cm(i) - 1.d0/(r*dxc(i))
         cp(i) = cp(i) + 1.d0/(r*dxc(i))
         c0(i) = c0(i) + 1.d0/(r**2)
         endif
      enddo
    
end subroutine read_grid

end module grid_module
