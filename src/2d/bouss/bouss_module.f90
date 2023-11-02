
module bouss_module

    use amr_module, only: t0, maxlv

#ifdef HAVE_PETSC
#include <petsc/finclude/petscksp.h>
    use petscksp
#endif

    implicit none
    save


    integer :: bc_xlo, bc_xhi, bc_ylo, bc_yhi

    real(kind=8) :: Bparam
    real(kind=8) :: alpha

    logical :: startWithBouss
    real(kind=8) :: startBoussTime

    !
!   Added for Boussinesq hybrid solver
! ================================
!
    type matrix_patchIndex
       integer, allocatable, dimension(:,:) :: mindex
       logical*1, allocatable, dimension(:,:) :: isBouss
       integer, allocatable, dimension(:,:) :: mindexSolo
       !!! TEST FOR 2 LEVEL CODE, REALLY WANT LISTS
       integer, allocatable, dimension(:,:) :: mfine
       integer, allocatable, dimension(:,:) :: ifine
       integer, allocatable, dimension(:,:) :: jfine
       integer, allocatable, dimension(:,:) :: mcoarse
       integer, allocatable, dimension(:,:) :: icoarse
       integer, allocatable, dimension(:,:) :: jcoarse
    end type

    type matrix_levInfo
    !! matrix_indices is indexed into by patch number, only Bouss Grids at that level
       type (matrix_patchIndex), allocatable, dimension(:) :: matrix_indices
       integer, allocatable, dimension(:) :: matrix_ia, matrix_ja
       real(kind=8), allocatable, dimension(:) :: matrix_sa
       integer :: numBoussGrids, numBoussCells, numBoussCellsSolo, numGhostCount, numUnset
       ! intfCountc are the # equations (cells)  added for you as a coarse grid
       ! these are the equations under a fine grid on the first interior cells
       ! that say the coarse cell is the conservative average of the refined cells

       ! intfCountf are the # equations (cells)  if you are the fine grid
       ! these are the equations to interpolate a ghost cell from the
       ! coarser patch
       integer :: intfCount
       integer :: matrix_nelt
    end type

    double precision :: boussMinDepth
    integer :: minLevelBouss, maxLevelBouss
    integer :: isolver
    integer :: ibouss

    type (matrix_levInfo), target :: matrix_info_allLevs(maxlv)

#ifdef HAVE_PARDISO
    ! vars for pardiso solver
    integer :: mtype,solver,maxfct,mnum,phase,nn,nrhs,error,msglvl
    integer iparm(64), idum
    real*8 dparm(64),ddum
    integer*8 pt(64,maxlv)
#endif

    ! if new grids have to redo factorization
    ! either first time in setgrd of in regrid
    logical newGrids(maxlv), newGridsTower(maxlv)
    integer*8  itcount(maxlv), numTimes(maxlv)

#ifdef HAVE_PETSC
    Mat Jr(maxlv)
    KSP ksp(maxlv)
    Mat JrTower(maxlv)
    KSP kspTower(maxlv)
#endif




contains

    subroutine set_bouss(rest)

    ! Set Bparam and bc choices for Boussinesq implicit solver

    ! Eventually allow more things to be specified in setrun.py
    ! and read them in here.

    use amr_module, only: mthbc, outunit
    implicit none
    logical rest
    integer iunit,i
    character(len=25) fname

#ifdef WHERE_AM_I
    write(*,*) 'starting set_bouss'
#endif

    iunit = 7
    fname = 'bouss.data'
!   # open the unit with new routine from Clawpack 4.4 to skip over
!   # comment lines starting with #:
    call opendatafile(iunit, fname)

    ! set constants for whichever bouss solver being used
    ! MS equation parameter B:
    Bparam = 1.d0/15.d0  

    ! SGN equation parameter alpha
    !alpha = 1.d0
    alpha = 1.153d0

    read(7,*) boussEquations

    ! modify write statements to say what value of Bparam or alpha is used:

    read(7,*) ibouss
    if (ibouss .eq. 0) then
       write(*,*)" Using Shallow Water equations"
       write(outunit,*)" Using Shallow Water equations"
    else if (ibouss .eq. 1) then
       write(*,*)" Using Madsen Sorensen equations"
       write(outunit,*)" Using Madsen Sorensen equations"
    else if (ibouss .eq. 2) then
       write(*,*)" Using SGN equations"
       write(outunit,*)" Using SGN equations"
    else
       write(*,*)" Unrecognized option for equation set"
       write(outunit,*)" Unrecognized option for equation set"
       stop
    endif

    ! CHECK ORDER!
    
    read(7,*) boussMinDepth
    read(7,*) minLevelBouss
    read(7,*) maxLevelBouss
    read(7,*) isolver
    read(7,*) startBoussTime


 99   write(*,900) minLevelBouss, maxLevelBouss
      write(outunit,900) minLevelBouss, maxLevelBouss
 900  format("==> Applying Bouss equations to selected grids between levels ",i3," and ",i3)

      write(*,*)"==> Use Bouss. in water deeper than ",boussMinDepth
      write(outunit,*)"==> Use Bouss. in water deeper than ",boussMinDepth

      if (isolver .eq. 1) then
         write(*,*)" No longer supporting GMRES solver"
         stop
      else if (isolver .eq. 2) then
#ifdef HAVE_PARDISO
         !write(*,*)" Using Pardiso solver"
         !write(outunit,*)" Using Pardiso solver"
         write(*,*)"Cannot use expired Pardiso solver"
         write(outunit,*)"Cannot use expired Pardiso solver"
         !stop
#else
         write(*,*)"need to install Pardiso for this option"
         stop
#endif
        else if (isolver .eq. 3) then
#ifdef HAVE_PETSC
         write(*,*)"Using a PETSc solver"
         write(outunit,*)"Using PETSc solver"
#else
         write(*,*)"need to install PETSc for this option"
         stop
#endif
      else
         write(*,*)"Unknown solver",isolver," choose 1,2 or 3"
         write(outunit,*)"Unknown solver",isolver," choose 1,2 or 3"
         stop
      endif

      if (startBoussTime .le. t0) then
         write(*,*)"Using Bouss equations from the start"
         write(outunit,*)"Using Bouss equations from the start"
         startWithBouss = .true.
      else
         write(*,*)"==> Wait until time ",startBoussTime," for before starting Bouss"
         write(*,*)"==> Using SWE until then."
         startWithBouss = .false.
      endif
     endif

    close(unit=iunit)



    ! Boundary conditions to impose in computing Boussinesq update:
    if ((mthbc(1)==2) .or. (mthbc(2)==2) &
        .or. (mthbc(3)==2) .or. (mthbc(4)==2)) then
        write(6,*) '*** Periodic BCs not supported in bouss_module'
        stop
    endif

    if ((mthbc(1)==4) .or. (mthbc(2)==4) &
        .or. (mthbc(3)==4) .or. (mthbc(4)==4)) then
        write(6,*) '*** Sphere BCs not supported in bouss_module'
        stop
    endif
    
    ! Dirichlet BCs using ghost cell values saved from coarser level:
    ! Requires num_eqn == 5.
    ! This is used in buildSparseMatrix at patch boundaries that are not
    ! domain boundaries, so the values of bc_xlo etc. only affect domain bdries
    ! Default (only if mthbc=0?) unless set otherwise by mthbc tests below

    bc_xlo = 2
    bc_xhi = 2
    bc_ylo = 2
    bc_yhi = 2
    

    ! For extrapolation in x:
    if (mthbc(1)==1) bc_xlo = 1
    if (mthbc(2)==1) bc_xhi = 1

    ! For extrapolation in y:
    if (mthbc(3)==1) bc_ylo = 1
    if (mthbc(4)==1) bc_yhi = 1


    ! For solid wall in x:
    if (mthbc(1)==3) bc_xlo = 3
    if (mthbc(2)==3) bc_xhi = 3

    ! For solid wall in y:
    if (mthbc(3)==3) bc_ylo = 3
    if (mthbc(4)==3) bc_yhi = 3
    
    
    ! FOR A PLANE WAVE CASE, SOLID WALL SHOULD WORK - Don't need what's below
    
    ! For plane wave in x-direction, should instead use Neumann conditions in y:
    ! Can use wall BC in SWE step so assume 
    !  clawdata.bc_lower[1] = clawdata.bc_upper[1] = 3
    !if (mthbc(3)==3) bc_ylo = 1
    !if (mthbc(4)==3) bc_yhi = 1
    !write(6,*) 'For plane wave in x-direction, using Neumann at top,bottom'


    ! For plane wave in y-direction, should instead use Neumann conditions in x:
    !if (mthbc(1)==3) bc_xlo = 1
    !if (mthbc(2)==3) bc_xhi = 1
    !write(6,*) 'For plane wave in y-direction, using Neumann at left,right'

#ifdef WHERE_AM_I
    write(*,*) 'ending   set_bouss'
#endif

    end subroutine set_bouss

end module bouss_module