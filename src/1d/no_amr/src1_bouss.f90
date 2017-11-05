
  subroutine src1(meqn,mbc,mx,xlower,dx,q,maux,aux,t,dt)
  
  ! This version includes:
  !   friction
  !   radial source terms (if radial==.true.)
  !   Boussinesq source terms (if bouss==.true.)
  
  ! Must be compiled with LAPACK for linear solvers.

      use grid_module, only: mx_grid, xgrid, radial
      use geoclaw_module, only: dry_tolerance, grav, DEG2RAD
      use geoclaw_module, only: ifrictiontype => friction_forcing
      use geoclaw_module, only: frictioncoeff => friction_coefficient
      use geoclaw_module, only: radial, bouss


      implicit none
            
      integer, intent(in) :: meqn,mbc,mx,maux
      real(kind=8), intent(inout) ::   q(meqn,1-mbc:mx+mbc)
      real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc)
      real(kind=8), intent(in) :: xlower, dx, t, dt

      real(kind=8) gamma, rcell, u

      real(kind=8)  q0(meqn,1-mbc:mx+mbc)
      real(kind=8) psi(mx+2)
     
      real(kind=8)   rk_stage(1:mx,4)
      real(kind=8) :: delt
      integer ::  i,k,ii,nstep,rk_order
      real(kind=8) tol,x
      integer(kind=4) INFO
     
      INTEGER            LDB
      INTEGER            IPIV( mx+2 )
      DOUBLE PRECISION   D(mx+2), DL(mx+1), DU(mx+1), DU2(1:mx)             
      
      
          ! Friction:

          if (frictioncoeff.eq.0.d0 .or. ifrictiontype.eq.0) return
            ! integrate source term based on Manning formula
          if (ifrictiontype.eq.1) then
            do i=1,mx
               if (q(1,i)<=dry_tolerance) then
                  q(2,i) = 0.0
               else
                  gamma= dsqrt(q(2,i)**2)*(grav*frictioncoeff**2)/(q(1,i)**(7.0/3.0))
                  q(2,i)= q(2,i)/(1.d0 + dt*gamma)
              endif
            enddo
          elseif (ifrictiontype.eq.2) then
            do i=1,mx
               if (q(1,i)<=dry_tolerance) then
                  q(2,i) = 0.0
               else
                  gamma= q(1,i)*grav*dtan(frictioncoeff*DEG2RAD)
                  gamma = max(0.d0, abs(q(2,i)) - dt*abs(gamma))
                  q(2,i) = sign(gamma, q(2,i))
              endif
            enddo
          endif

      !    ----------------------------------------------------------------
          ! Radial source terms:

          if (radial) then
                ! radial source term:
                do i=1,mx
                   ! assume radial about left edge!
                   rcell = 0.5*(xgrid(i) + xgrid(i+1)) - xgrid(1) 
                   q(1,i) = q(1,i) - dt/rcell * q(2,i)
                   if (q(1,i) .gt. dry_tolerance) then
                       u = q(2,i)/q(1,i)
                     else
                       u = 0.d0
                     endif
                   q(2,i) = q(2,i) - dt/rcell * q(1,i)*u**2
                   enddo
               endif

      !    ----------------------------------------------------------------

      ! Boussinesq terms:

       if (bouss) then
  !
         nstep = 1  ! > 1 to sub-step the source terms

         k = 1

         delt=dt/nstep

         LDB = mx+2
         D  =0.d0
         DU =0.d0
         DL =0.d0
         DU2=0.d0

         rk_order = 2  ! 2 or 4
        
         do ii=1,nstep

              call set_diag(mx,meqn,mbc,dx,q,maux,aux, DL, D, DU)

              call DGTTRF( mx+2, DL, D, DU, DU2, IPIV, INFO )
        
              psi = 0.d0
              rk_stage = 0.d0
              q0  = q

              !-----------------------
              if (rk_order == 2) then
              ! RK2   

              ! First Stage
          
              call set_psi(mx,meqn,mbc,dx,q0,maux,aux,psi)


              call DGTTRS( 'N', mx+2, 1, DL, D, DU, DU2 &
                          , IPIV, psi, LDB,INFO )


              rk_stage(1:mx,1) = psi(2:mx+1)
            
              q0(2,1:mx)=q(2,1:mx)-delt/2.d0*rk_stage(1:mx,1)
            
              ! Second Stage

              call set_psi(mx,meqn,mbc,dx,q0,maux,aux,psi)

              call DGTTRS( 'N', mx+2, 1, DL, D, DU, DU2 &
                          , IPIV, psi, LDB,INFO )

              rk_stage(1:mx,2)=psi(2:mx+1)
            
              q(2,1:mx) = q(2,1:mx)- delt*rk_stage(1:mx,2)

              endif ! rk_order==2

              !-----------------------
              if (rk_order == 4) then
              ! RK4   

              ! First Stage
          
              call set_psi(mx,meqn,mbc,dx,q0,maux,aux,psi)


              call DGTTRS( 'N', mx+2, 1, DL, D, DU, DU2 &
                          , IPIV, psi, LDB,INFO )


              rk_stage(1:mx,1) = psi(2:mx+1)
            
              q0(2,1:mx)=q(2,1:mx)-delt/2.d0*rk_stage(1:mx,1)
            
              ! Second Stage

              call set_psi(mx,meqn,mbc,dx,q0,maux,aux,psi)

              call DGTTRS( 'N', mx+2, 1, DL, D, DU, DU2 &
                          , IPIV, psi, LDB,INFO )

              rk_stage(1:mx,2)=psi(2:mx+1)
            
              q0(2,1:mx)=q(2,1:mx)-delt/2.d0*rk_stage(1:mx,2)
            
              ! Third Stage
            
              call set_psi(mx,meqn,mbc,dx,q0,maux,aux,psi)

              call DGTTRS( 'N', mx+2, 1, DL, D, DU, DU2 &
                          , IPIV, psi, LDB,INFO )
                    
              rk_stage(1:mx,3)=psi(2:mx+1)
            
              q0(2,1:mx)=q(2,1:mx)-delt*rk_stage(1:mx,3)
            
              ! Fourth Stage
            
              call set_psi(mx,meqn,mbc,dx,q0,maux,aux,psi)

              call DGTTRS( 'N', mx+2, 1, DL, D, DU, DU2 &
                          , IPIV, psi, LDB,INFO )

              rk_stage(1:mx,4)=psi(2:mx+1)

                                  
              q(2,1:mx) = q(2,1:mx)- delt/6.d0*(rk_stage(1:mx,1) &
                   + 2.d0*rk_stage(1:mx,2) &
                   + 2.d0*rk_stage(1:mx,3) + rk_stage(1:mx,4))

              endif ! rk_order==4

          enddo
          
        endif ! end of Bouss terms
      
  end subroutine src1

  
! =========================================================
      subroutine set_diag(mx,meqn,mbc,dx,q,maux,aux, DL, D, DU)
! =========================================================
      use bous_module
      use geoclaw_module, only: sea_level
      use grid_module, only: xcell, dxm, dxc, cm, cp, c0

      implicit none
     
      integer(kind=4) mx,meqn,mbc,maux,mxy
      integer(kind=4) i,k
      real(kind=8) dx

      real(kind=8)   q(meqn,1-mbc:mx+mbc)
      real(kind=8) aux(maux,1-mbc:mx+mbc)

      DOUBLE PRECISION   D(mx+2), DL(mx+1), DU(mx+1)                
      
      real(kind=8), dimension(1-mbc:mx+mbc) :: h0, hh, hhh
      
      h0 = sea_level - aux(1,:)
     
        do i=1-mbc,mx+mbc
           hh(i)=(max(0.,h0(i)))**2
           hhh(i)=(max(0.,h0(i)))**3
        enddo
      
        D = 1.d0
        DU= 0.d0
        DL= 0.d0

   
        do i=1,mx
        
            if ((minval(q(1,i-2:i+2)) < 0.1d0) &
                .or. (minval(h0(i-2:i+2)) < 0.1d0)) then

            !if ((q(1,i)<0.1d0) .or. (h0(i)<0.1d0) .or. &
            !    (minval(h0(i-2:i+2))<0.1d0) .or. &
            !    (maxval(q(1,i-2:i+2))<0.1d0) .or. &
            !    (minval(h0(i-1:i+1))<=0.1d0)) then

              ! do nothing
        
            else


              D(i+1) = 1.d0 + c0(i)*((B_param+.5d0)*hh(i) &
                    - 1.d0/6.d0*hhh(i)/h0(i))
     
              DU(i+1)= -cp(i)*((B_param+.5d0)*hh(i) &
                    - 1.d0/6.d0*hhh(i)/h0(i+1))
     
              DL(i)= -cm(i)*((B_param+.5d0)*hh(i) &
                    - 1.d0/6.d0*hhh(i)/h0(i-1))
     
            endif
            
        enddo

      if (.true.) then
         ! for wall-reflecting BC at left: impose Q_0 = -Q_1
         D(2) = D(2) - DL(1)   
         DL(1) = 0.d0     
         endif
             
      return
      end subroutine set_diag
   
!======================================================================
      subroutine set_psi(mx,meqn,mbc,dx,q,maux,aux,psi)

      use geoclaw_module, only: g => grav, sea_level
      use bous_module
      use grid_module, only: xcell, dxm, dxc, cm, cp, c0
     
      implicit none

      integer, intent(in) :: mx,meqn,mbc,maux
      real(kind=8), intent(in) :: dx
     
      real(kind=8), intent(in) ::  q(meqn,1-mbc:mx+mbc)
      real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc)
      real(kind=8), intent(out) :: psi(mx+2)

      real(kind=8)  tol
      real(kind=8)  depth
      real(kind=8), dimension(1-mbc:mx+mbc) :: h0, hh, hhh, eta, hu2
      real(kind=8), dimension(1-mbc:mx+mbc) :: hetax,detax,s1,s1_h
      real(kind=8)  detaxxx,s1xx,s1_hxx
      integer :: i,j,k,kk,iL,iR

      ! phase out dispersive term between sw_depth0 and sw_depth1:
      ! parameters now in bous_module
      !sw_depth0 = 180.d0  ! pure SWE in shallower water than this
      !sw_depth1 = 190.d0  ! full Bouss in deeper water than this
    
      h0 = sea_level - aux(1,:)
     
         do i=1-mbc,mx+mbc
           if (q(1,i).gt.1d-4) then
              hu2(i)= q(2,i)**2/q(1,i)
             else
              hu2(i)= 0.d0
             endif
           eta(i)= q(1,i) - h0(i)
           hh(i)=(max(0.d0,h0(i)))**2
           hhh(i)=(max(0.d0,h0(i)))**3
        enddo
     
     hetax = 0.d0
     detax = 0.d0
     !write(6,*) '+++ eta(0:1): ',eta(0),eta(1)

         do i= 1,mx
         if (minval(h0(i-1:i+1)) > 0.d0) then
            if (i<1) then
                !write(6,*) '+++ dxm = ', (dxm(j), j=1,4)
                hetax(i)= q(1,i)*(eta(i+1)-eta(i))/dxm(i+1)
                detax(i)= h0(i)*(eta(i+1)-eta(i))/dxm(i+1)
            elseif (i>mx) then
                hetax(i)= q(1,i)*(eta(i)-eta(i-1))/dxm(i)
                detax(i)= h0(i)*(eta(i)-eta(i-1))/dxm(i)
            else
                hetax(i)= q(1,i)*(eta(i+1)-eta(i-1))/dxc(i)
                detax(i)= h0(i)*(eta(i+1)-eta(i-1))/dxc(i)
            endif
         endif
        enddo
     
     s1=0.d0
     s1_h=0.d0

     !write(6,*) '+++ hu2(0:1): ',hu2(0),hu2(1)
     do i=1,mx
          if (i<1) then
             s1(i)= (hu2(i+1)-hu2(i))/dxm(i+1)
          elseif (i>mx) then
             s1(i)= (hu2(i)-hu2(i-1))/dxm(i)
          else
             s1(i)= (hu2(i+1)-hu2(i-1))/dxc(i)
          endif
         
          s1(i)= s1(i)+g*hetax(i)
      
           if (h0(i) > 0.1d0) then 
              s1_h(i)=s1(i)/h0(i)
           else
              s1_h(i)=0.d0
           endif
     enddo
     
      tol = 1d-8

     !write(6,*) '+++ s1(0:1): ',s1(0),s1(1)
     !write(6,*) '+++ s1_h(0:1): ',s1_h(0),s1_h(1)
     !write(6,*) '+++ detax(0:1): ',detax(0),detax(1)
     s1(0) = s1(1)
     s1_h(0) = s1_h(1)
     detax(0) = detax(1)

      psi = 0.d0

        do i=1,mx

           k = i+1

           if ((h0(i) < sw_depth0) .or. (q(1,i) < sw_depth0)) then
              ! no dispersive term:
              psi(k) = 0.d0

           else
            
              ! check if we need to use one-sided differences:
              if (i<1) then
                 j = i+1
              elseif (i>mx) then
                 j = i-1
              else
                 j = i
              endif

              s1xx = cp(j)*s1(j+1) - c0(j)*s1(j) + cm(j)*s1(j-1)
              detaxxx = cp(j)*detax(j+1) - c0(j)*detax(j) + &
                        cm(j)*detax(j-1)
              s1_hxx = cp(j)*s1_h(j+1) - c0(j)*s1_h(j) + cm(j)*s1_h(j-1)

              psi(k) = (B_param+.5d0)*hh(i)*s1xx - B_param*g*hh(i)*detaxxx &
                  -hhh(i)/6.d0 * s1_hxx
                 
           endif


        if (h0(i) .lt. sw_depth1) then
            ! reduce psi linearly between sw_depth0 and sw_depth1:
            psi(k) = psi(k) * (h0(i)-sw_depth0)/(sw_depth1-sw_depth0)
            endif
                
        enddo

      return
      end subroutine set_psi
