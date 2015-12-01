!=============================================================================
!   PDE - ODE Coupled Solver
!-----------------------------------------------------------------------------
!     This fortran code solves the PDE - ODE coupled system that describes 
!   the growth of Clostridium Thermocellum and its consumption of the carbon
!   substrait.
!   
!   The system is based off the following system:
!     M_t = nabla (D(M) nabla M ) + f(C,M) 
!     C_t = - g(C,M)
!   where M is the biomass of C. Thermocellum and C is the concentration of
!   Carbon. D(M) is the diffusion coeffient for the biomass movement. f(C,M) 
!   is the growth and death term for the biomass. g(C,M) is the consumption
!   of carbon substrait.
!-----------------------------------------------------------------------------
!     The method of solving is by trapzidral rule for C, and by finite
!   difference for M 
!==============================================================================

program cThermoPDEODE
!    use omp_lib
    implicit none
    INTERFACE
        subroutine solveLSDIAG(n,ndiag,ioff,a,sol,rhs,nit,err1,err2,stopcrit)
        !---------------------------------------------------------------------
        ! input:  n       problem size
        !         ndiag:  number of diagonals
        !         ioff:   offsets (distance of sub diagonals to main diagonal)
        !         Mnew:      matrix values
        !         sol:    initial guess for iteration ( overwritten by result)
        !         rhs:    the righ hand side of the linear system
        !         nit:    max num of iter to be carried out (overwritten)
        !         err1:   tol for 1st stopping criterion (will be overwritten)
        !         err2:   tol for 2nd stopping criterion (will be overwritten)
        ! output: sol:    solution of Ax=b
        !         nit:    number of iterations taken
        !         err1:   computed value for 1st stopping criterion
        !         err2:   computed value for 2nd stopping criterion
        !         stopcrit: value that tells which stopping criti activated
        !---------------------------------------------------------------------
        implicit none
        integer, intent(in)               :: n,ndiag
        real,dimension(n,ndiag),intent(in):: a
        integer, dimension(ndiag),intent(in)::ioff
        real,dimension(n),intent(in)      :: rhs
        real,dimension(n),intent(inout)   :: sol
        real,intent(inout)                :: err1,err2
        integer,intent(inout)             :: nit
        integer,intent(out)               :: stopcrit
        end subroutine
    END INTERFACE
  
    !=======================
    ! Variable Descriptions
    !=======================
    
    ! filename Variables
    character(100) :: filename

    ! Problem Parameters
    real :: kappa,delta,nu,gama
    integer :: alpha,beta
    real :: depth,height
    integer :: fSelect, dSelect, gSelect, MinitialCond
    integer :: num_innocu_points

    ! Numerical Method Parameters
    integer :: pSize,row,col,n,ndiag,nit,nOuts
    real :: tEnd,tDel,xDel,e1,e2,eSoln,eTrav
    integer :: true2D, checkTrav

    ! Solution variables
    real,dimension(:),allocatable  :: C, M
    
    ! Reporting variables
    real :: avgIters, maxIters
    real :: avgNit, maxNit
    integer :: startTime, endTime, clock_rate, clock_max
    real :: elapsedTime
    
    write(*,*) "Enter parameter file name: ex. 'parameter.txt' "
    write(*,'(A)', advance="no") "    "
    read(*,*) filename

    write(*,*) "Setting Parameters that shouldn't change"
    ndiag = 5
    write(*,'(A, I5)') "     ndiag = ", ndiag

    write(*,*) "Getting problem size from file"
    call getProbSize(pSize, true2D, checkTrav, filename, len(filename)) 
    if (true2D == 1) then 
      write(*,*) "    Problem is 2D"
      col = 4
    else if (true2D == 0) then 
      write(*,*) "    Problem is 3D"
      col = pSize
    end if
    write(*,*) "    pSize = ", pSize
    row = pSize
    n = row * col
    write(*,*) "    row   = ", row
    write(*,*) "    col   = ", col        
    write(*,*) "    n     = ", n

    write(*,*) "Opening Parameter file"
    call paramSet(depth, height, num_innocu_points, alpha, beta, kappa, &
                  gama, nu, delta, nOuts, tEnd, tDel, e1, e2, eTrav, &
                  eSoln, fSelect, dSelect, gSelect, MinitialCond, filename, &
                  len(filename))
    xDel = 1/real(row)
    nit = n * 100
    write(*,*) "Parameters set:"
    write(*,*) "    depth       = ", depth
    write(*,*) "    height      = ", height
    write(*,*) "    alpha       = ", alpha
    write(*,*) "    beta        = ", beta
    write(*,*) "    gama        = ", gama
    write(*,*) "    nu          = ", nu
    write(*,*) "    delta       = ", delta
    write(*,*) "    nOuts       = ", nOuts
    write(*,*) "    tEnd        = ", tEnd
    write(*,*) "    tDel        = ", tDel
    write(*,*) "    e1          = ", e1
    write(*,*) "    e2          = ", e2
    write(*,*) "    eSoln       = ", eSoln
    write(*,*) "    nit         = ", nit
    write(*,*) "    xDel        = ", xDel
    write(*,*) "    fSelect     = ", fSelect
    write(*,*) "    dSelect     = ", dSelect
    write(*,*) "    gSelect     = ", gSelect
    write(*,*) "    MinitialCond= ", MinitialCond
   write(*,*) "Allocating the size of arrays"
    allocate(C(n),M(n))
    write(*,'(A,I12,A)') "    C and M are now dimension(", n, ") arrays"

    write(*,*) "Setting Initial Conditions"
    call setInitialConditions(C,M,row,col,n,depth,height,xDel,MinitialCond, &
                              num_innocu_points)

    write(*,*) "Starting Solver"
    call system_clock(COUNT_RATE=clock_rate, COUNT_MAX=clock_max)
    call system_clock(startTime)
    call solveOrder2(tEnd,nOuts,tDel,n,row,col,M,C,xDel,&
       gama,nu,alpha,beta,kappa,delta,ndiag,e1,e2,nit,eSoln,dSelect,&
       fSelect,gSelect,avgIters,maxIters,avgNit,maxNit,true2D,checkTrav,eTrav)
    call system_clock(endTime)
   
    ! Convert times into seconds
    elapsedTime = real(endTime - startTime)/real(clock_rate)


    ! Report Statistics
    write(*,*) "-------------------------------------------------------------"
    write(*,*) "Statsitcs:"
    write(*,*) "-------------------------------------------------------------"
    write(*,*) "Time to compute = ", elapsedTime
    write(*,*) ""
    write(*,*) "Avg Iters for iterating betn. soln. =", avgIters
    write(*,*) "Max Iters for iterating betn. soln. =", maxIters
    write(*,*) "Avg Iters for linear solver =", avgNit
    write(*,*) "Max Iters for linear solver =", maxNit  
    write(*,*) "Writing statistics to file"
    call reportStats(avgIters,maxIters,avgNit,maxNit,elapsedTime)
    
    write(*,*) "Execution completed"
  
end program


!==============================================================================
!   Sets the problem size and checks for auxillary settings
!-----------------------------------------------------------------------------
!     This must be done first since the problem size is used in the variable 
!   declaration of the arrays.
!     The auxillary settings are if the problem is 2D in nature; this means 
!   that the computation time can be drastically reduced by letting col = 4
!   instead of col = pSize.
!     Also, the option for checking for the existance of travelling wave
!   solutions is done here.
!==============================================================================
subroutine getProbSize(pSize, true2D, checkTrav, filename, nameLen) 
    implicit none
    integer,intent(out) :: pSize, true2D, checkTrav
    integer,intent(in)  :: nameLen
    character(nameLen),intent(in) :: filename
    character :: dum   ! Dummy variable
    write(*,*) filename
    open(UNIT = 19, FILE = filename, STATUS = "old", ACTION = "read")
    read(19,*); read(19,*); read(19,*)  ! Skip first 3 lines
    read(19,*) dum, dum, pSize
    read(19,*) dum, dum, true2D
    read(19,*) dum, dum, checkTrav
    close(19)
end subroutine getprobSize


!==============================================================================
!   Sets the all the parameter values
!-----------------------------------------------------------------------------
!     Opens the parameter.txt file, which holds all the parameter values and 
!   function uses, and sets the variables accordingly.
!     Takes in all the parameters on entry and outputs them with the appropriate
!   value.
!==============================================================================
subroutine paramSet(depth, height, numInnocu, alpha, beta, kappa, &
                  gama, nu, delta, nOuts, tEnd, tDel, e1, e2, eTrav, &
                  eSoln, fSelect, dSelect, gSelect, MinitialCond, filename, &
                  nameLen)
    implicit none
    real, intent(out) :: depth, height, kappa, delta, nu 
    real, intent(out) :: tEnd, tDel, e1, e2, eSoln, eTrav, gama
    integer, intent(out) :: alpha, beta,  numInnocu, nOuts
    integer, intent(out) :: fSelect, dSelect, gSelect, MinitialCond
    integer, intent(in) :: nameLen
    character(nameLen), intent(in) :: filename

    character :: dum, dum2      ! Dummy variable
    
    open(UNIT = 19, FILE = filename, STATUS = "old", ACTION = "read")
    ! Skip first 6 lines
    read(19,*); read(19,*); read(19,*); read(19,*); read(19,*); read(19,*)

    read(19,*) dum, dum2, depth
    read(19,*) dum, dum2, height
    read(19,*) dum, dum2, numInnocu
    read(19,*) dum, dum2, alpha
    read(19,*) dum, dum2, beta
    read(19,*) dum, dum2, nu
    read(19,*) dum, dum2, kappa 
    read(19,*) dum, dum2, gama 
    read(19,*) dum, dum2, delta
    read(19,*) dum, dum2, nOuts
    read(19,*) dum, dum2, tEnd
    read(19,*) dum, dum2, tDel
    read(19,*)                          ! Skip xDel 
    read(19,*) dum, dum2, e1
    read(19,*) dum, dum2, e2
    read(19,*) dum, dum2, eSoln
    read(19,*) dum, dum2, eTrav
    read(19,*)                          ! Skip nit
    read(19,*); read(19,*); read(19,*)  ! Skip 3 lines
    read(19,*) dum, dum2, fSelect
    read(19,*) dum, dum2, dSelect
    read(19,*) dum, dum2, gSelect
    read(19,*) dum, dum2, MinitialCond

    close(19)
    
end subroutine paramSet


!==============================================================================
!   Sets the initial conditions for the system
!-----------------------------------------------------------------------------
!   For C, we have homogenous initial conditions, trivial to do
!   For M, the innoculation point has a smooth curve to avoid sharp numerical 
!     artifacts in the system. This is done with a polynomial f(x) = a*x^8+b, 
!     this function is computed based on the depth and height parameters, 
!     calculating b = height and a = b/(depth)^8
!==============================================================================
subroutine setInitialConditions(C,M,row,col,n,depth,height,xDel,MinitialCond, &
                                num_innocu_points)
    implicit none
    integer,intent(in) :: row,col,n,MinitialCond, num_innocu_points
    real,intent(in) :: depth, height, xDel
    real,dimension(n),intent(out) :: C,M

    integer :: i,j,x,y, xi, yi, p
    real :: f,a       ! function for IC curve
    real :: innoc     ! Placeholder for new IC value at point
    real :: total     ! Counts total innoculation height
    real :: mIC, cIC      ! Used to store previous simulation values while reading
    real :: xr, yr
    
    C = 1.; j = 0; M = 0

    ! 2D trav wave smooth curve
    if (MinitialCond == 1) then
        a = -height/(depth)**4
        !$omp parallel do private(f) shared(height,a) 
        do i = 1, n
          x = (i-1)/col
          f = a*(x*xDel)**4 + height
          M(i) = f
          if (f .LE. 0) M(i) = 0
        enddo 
        !$omp end parallel do

    ! homogenous everywhere
    else if (MinitialCond == 2) then
        M = height

    ! 3D initial conditions
    else if (MinitialCond == 3) then
        a = -height/(depth)**2
        !$omp parallel do private(f,x,y) shared(height,a) 
        do i = 1, n
          x = MOD(i-1, col)
          y = i / row
          f = a*((x*xDel-0.5)*(x*xDel-0.5) + (y*xDel-0.5)*(y*xDel-0.5)) + height
          M(i) = f
          if (M(i) .LE. 0) M(i) = 0
        enddo
        !$omp end parallel do

    !(Random innoclation points over whole domain [3D])
    else if (MinitialCond == 4) then
      a = -height/(depth*depth)
      call init_random_seed()
      do i = 1, num_innocu_points
        call random_number(xr)
        call random_number(yr)
        do j = 1, n
          if (M(j) .LE. 0) then
            x = MOD(j-1, col)
            y = j / row
            M(j) = M(j) + a*((x*xDel-xr)*(x*xDel-xr) + (y*xDel-yr)*(y*xDel-yr)) + height
            if (M(j) .LE. 0) M(j) = 0
          endif
        enddo
      enddo 

    !(Random innoculation points on y=0 side [3D])
    else if (MinitialCond == 5) then
      a = -height/(depth*depth)
      call init_random_seed()
      do i = 1, num_innocu_points
        call random_number(xr)
        call random_number(yr)
        yr = yr*0.1
        do j = 1, n
            x = MOD(j-1, col)
            y = j / row
            innoc = a*((x*xDel-xr)*(x*xDel-xr) + (y*xDel-yr)*(y*xDel-yr)) + height
            if (innoc .GE. 0) M(j) = M(j) + innoc 
        enddo
      enddo
        ! Divide innoculation height by average non-zero value
      i = 0
      do j = 1, n
        if (M(j) .GT. 0) then
          i = i + 1
          total = total + M(j)
        endif
      enddo
      do j = 1, n
        M(j) = M(j) * height/(total/real(i))
      enddo


    ! sharp IC. homogenous in y-dir
    else if (MinitialCond == 6) then
      do i = 1, n
        x = (i-1) / col
        if ( x*xDel .LE. depth) then
          M(i) = height
        endif
      enddo 

    ! pertubations in y-dir; check trav.wave. stability
    else if (MinitialCond == 7) then
      call init_random_seed()
      do i = 1, n
        if ( i*xDel/col .LE. depth) then
          call random_number(xr)
          M(i) = xr*0.2
        endif
      enddo

    ! Test the grid ordering/ printing
    else if (MinitialCond == 8) then
     ! do i = 1, row
     !   do j = 1, col
     !     if (i*xDel .LE. depth) then
     !       p = j + (i-1)*col
     !       M(i) = real(p)/real(n)/4.0
     !     endif
     !   enddo
     ! enddo
     do i = 1, n
       M(i) = real(i)/real(n)/10
     enddo
  
    !(Evenly spaced innoculation points on y=0 side [3D])
    else if (MinitialCond == 9) then
      a = -height/(depth*depth)
      do i = 1, num_innocu_points
        xr = depth + (i-1)*2*depth + real((i-1)*(1-num_innocu_points*depth*2))/real(num_innocu_points-1) 
        do j = 1, n
          if (M(j) .LE. 0) then
            x = MOD(j-1, col)
            y = j / row
            M(j) = M(j) + a*((x*xDel-xr)*(x*xDel-xr) + (y*xDel)*(y*xDel)) + height
            if (M(j) .LE. 0) M(j) = 0
          endif
        enddo
      enddo

    ! Random innoculation points on y=0 and y = 1
    else if (MinitialCond == 10) then
      a = -height/(depth*depth)
      call init_random_seed()
      do i = 1, 2*num_innocu_points
        call random_number(xr)
        call random_number(yr)
        yr = yr*0.1
        if (i .GT. num_innocu_points) then 
          yr = yr + 0.9 ! For y = 1 side
        endif
        do j = 1, n
            x = MOD(j-1, col)
            y = j / row
            innoc = a*((x*xDel-xr)*(x*xDel-xr) + (y*xDel-yr)*(y*xDel-yr)) + height
            if (innoc .GE. 0) M(j) = M(j) + innoc 
        enddo
      enddo

    ! Read IC from output file... (Continues a previous sim) CANNOT GET WORKING>>>>
!    else if (MinitialCond == 11) then
!      open(UNIT = 119, FILE = "initialCond.dat", STATUS = "old", ACTION = "read")
!      read(119,*) ! Skip first line
!      do i = 1, n/2+1
!        read(119,*) innoc, innoc, innoc, innoc
!      enddo
!      do i = n/2+1, n      <<<<< THIS IS THE PROBLEM AREA, Can't divide the region in two
!        read(119,*) innoc, innoc, mIC, cIC ! innoc is used as a dummy variable here
!        M(i) = mIC 
!        C(i) = cIC 
!      enddo
!      close(119)
!

    ! Clumped up innoculation at origin. Same total as MinitCond == 9
    else if (MinitialCond == 12) then
      innoc = (2*height*depth*depth*num_innocu_points) ** (1.0/3.0)
      a = -1.0/(innoc)
      do j = 1, n
        x = MOD(j-1, col)
        y = j / row
        M(j) = a*((x*xDel)*(x*xDel) + (y*xDel)*(y*xDel)) + innoc 
        if (M(j) .LE. 0) M(j) = 0
      enddo

    ! Evenly Spaced innocu point in cubes (exact total)
    else if (MinitialCond == 13) then
      xr = real(col)/real(num_innocu_points)
      do j = 1, n
        x = MOD(j-1, col)
        y = j / row
        if (y*xDel .LE. depth .AND. MOD(x+xr/4.0, xr) >= xr/2.0) M(j) = height
      enddo

    ! Clumped up like ==12, but exact cube
    else if (MinitialCond == 14) then
      xr = (num_innocu_points * depth * height * real(col)/(2*real(num_innocu_points))) ** (1.0/3.0)
      do j = 1, n
        x = MOD(j-1, col)
        y = j / row
        if (y*xDel .LE. xr .AND. x*xDel .LE. xr) M(j) = xr
      enddo
      
    endif

end subroutine setInitialConditions


!==============================================================================
!   Function f(C,M)
!==============================================================================
subroutine fFunc(M,C,f,kappa,nu,fSelect)
    implicit none    
    integer, intent(in) :: fSelect
    real, intent(in) :: M,C,kappa,nu
    real, intent(out) :: f
    real :: eps
    eps = C*0.00000001
    if (fSelect == 1) f = C/(kappa+C)-nu 
    if (fSelect == 2) f = C/(kappa+C)
    if (fSelect == 3) f = C/(kappa*M+C+eps)-nu
    if (fSelect == 4) f = C/(kappa*M+C+eps)
    if (fSelect == 5) f = 1
    if (fSelect == 6) f = C/(kappa+C)*(1-M/(C+eps))-nu
end subroutine fFunc


!==============================================================================
!   Function d(M)
!==============================================================================
subroutine dFunc(M,d,delta,alpha,beta,dSelect)
    implicit none    
    integer,intent(in) :: alpha, beta, dSelect
    real, intent(in) :: M, delta
    real, intent(out) :: d
    if (dSelect == 1) d = delta
    if (dSelect == 2) d = delta * M** alpha
    if (dSelect == 3) d = delta * M**alpha / ((1 - M)**beta)
end subroutine dFunc


!==============================================================================
!   Solves the solution, C
!-----------------------------------------------------------------------------
!   Uses the trapizoidal rule for an ODE and solves for C.
!   Different gSelects give different values for b and c. 
!   gSelect = 1 --> g = gama*C/(kappa+C)
!==============================================================================
subroutine solveC(M,Mnew,Csol,Cnew,n,kappa,gama,tDel,gSelect)
    implicit none 
    integer,intent(in) :: n,gSelect
    real,intent(in) :: kappa,tDel,gama
    real,dimension(n),intent(in) :: M,Mnew,Csol
    real,dimension(n),intent(out) :: Cnew

    integer :: i,j    ! grid index
    integer :: g      ! current grid point
    real :: b,c       ! quadratic equation terms
    real :: f         ! f(C,M) value at gridpoint
    real :: aux 
    
    real :: r,s
    r = 1
    s = 0.5

    aux = tDel*0.5*gama
    
    if (gSelect == 1) then
     !$omp parallel do private(b,c) shared(kappa,Csol,aux,Mnew,M)
      do i = 1,n 
        b = kappa - Csol(i) + aux*(Mnew(i) + Csol(i)*M(i)/(kappa+Csol(i)) )
        c = -kappa*Csol(i) + aux*kappa*Csol(i)*M(i)/(kappa + Csol(i))
        Cnew(i) = 0.5 * (-b + SQRT(b*b - 4 * c))
      enddo
     !$omp end parallel do
    end if    
      
    end subroutine solveC


!====================================================================
!   Generates matrix M in diagonal format for the current timestep
!====================================================================
subroutine GenMatrixM(M,C,MatrixM,Mioff,Mrhs,row,col,n,ndiag,delta,nu,&
                      alpha,beta,kappa,tDel,xDel,dSelect,fSelect)
    implicit none
    integer,intent(in) :: row,col,n,ndiag,alpha,beta,dSelect,fSelect
    real,intent(in) :: delta,kappa,xDel,tDel,nu
    real,dimension(n),intent(in) :: M,C

    real,dimension(n,ndiag),intent(out) :: MatrixM
    real,dimension(n),intent(out) :: Mrhs
    integer,dimension(ndiag),intent(out) :: Mioff

    real :: xCof
    real :: f
    real,dimension(n) :: diff
    !integer :: x,y,g
    integer :: i

    real :: tDela
    tDela = tDel      ! Used for testing purposes

    xCof = 1/(xDel*xDel)

    Mioff = (/ -col, -1, 0, 1, col/)
    MatrixM(:,:) = 0

    ! Compute all the diffusion coefficients
    !$omp parallel do shared(M,diff,delta,alpha,beta,dSelect)
    do i = 1,n
        call dFunc(M(i), diff(i), delta, alpha, beta, dSelect)
    enddo
    !$omp end parallel do

    !!$omp parallel do shared(MatrixM,xCof,diff,Mrhs,tDel,M,C,kappa,nu,fSelect) private(i, f)
    !$omp parallel do private(f)
    do i = 1,n
        if (i .LE. col) then
            MatrixM(i,5) = MatrixM(i,5) - xCof*0.5*(diff(i+col)+diff(i))
            MatrixM(i,3) = MatrixM(i,3) + xCof*0.5*(diff(i+col)+diff(i))
        else
            MatrixM(i,1) = MatrixM(i,1) - xCof*0.5*(diff(i-col)+diff(i))
            MatrixM(i,3) = MatrixM(i,3) + xCof*0.5*(diff(i-col)+diff(i))         
        endif
          
        if (MOD(i,col) == 1) then
            MatrixM(i,4) = MatrixM(i,4) - xCof*0.5*(diff(i+1)+diff(i))
            MatrixM(i,3) = MatrixM(i,3) + xCof*0.5*(diff(i+1)+diff(i))
        else
            MatrixM(i,2) = MatrixM(i,2) - xCof*0.5*(diff(i-1)+diff(i))
            MatrixM(i,3) = MatrixM(i,3) + xCof*0.5*(diff(i-1)+diff(i))         
        endif

        if (MOD(i,col) == 0) then
            MatrixM(i,2) = MatrixM(i,2) - xCof*0.5*(diff(i-1)+diff(i))
            MatrixM(i,3) = MatrixM(i,3) + xCof*0.5*(diff(i-1)+diff(i))
        else
            MatrixM(i,4) = MatrixM(i,4) - xCof*0.5*(diff(i+1)+diff(i))
            MatrixM(i,3) = MatrixM(i,3) + xCof*0.5*(diff(i+1)+diff(i))         
        endif
          
        if  (i .GE. n-col) then
            MatrixM(i,1) = MatrixM(i,1) - xCof*0.5*(diff(i-col)+diff(i))
            MatrixM(i,3) = MatrixM(i,3) + xCof*0.5*(diff(i-col)+diff(i))
        else
            MatrixM(i,5) = MatrixM(i,5) - xCof*0.5*(diff(i+col)+diff(i))
            MatrixM(i,3) = MatrixM(i,3) + xCof*0.5*(diff(i+col)+diff(i))         
        endif
        
        call fFunc(M(i),C(i),f,kappa,nu,fSelect)
        MatrixM(i,3) = MatrixM(i,3) - f + (1/tDel)
        Mrhs(i) = M(i)/tDel
    enddo
    !$omp end parallel do

end subroutine GenMatrixM


!==============================================================================
!   Solve Order 2
!-----------------------------------------------------------------------------
!   1. Solves for M_{i+1} using C_i and M_i
!   2. Solves for C_{i+1} using C_i, M_i, and M_{i+1}
!   ... repeat until convergence
!==============================================================================
subroutine solveOrder2(tEnd,nOuts,tDel,n,row,col,M,C,xDel,&
    gama,nu,alpha,beta,kappa,delta,ndiag,e1,e2,nit,eSoln,dSelect,fSelect,&
    gSelect,avgIters,maxIters,avgNit,maxNit,true2D,checkTrav,eTrav)
  implicit none
  integer,intent(in) :: nOuts,n,row,col,alpha,beta,ndiag
  integer,intent(in) :: dSelect,fSelect,gSelect,true2D,checkTrav
  real,intent(in) :: tEnd,tDel,xDel,kappa,delta,nu,gama
  real,intent(in) :: eSoln,eTrav
  real,intent(inout) :: e1,e2
  integer,intent(inout) :: nit
  real,dimension(n),intent(out) :: M,C
  real,intent(out) :: avgIters,maxIters,avgNit,maxNit
  
  integer :: counter        ! Counts the num. of iter. in system solver
  integer :: endLoop        ! endLoop = 1 -> solving loop can stop
  integer :: filter         ! Controls frequency of outputs written
  
  real :: stored_e1, stored_e2
 
  real :: totalMassM        ! Total Biomass over region
  real :: totalMassC        ! Total Substrait over region
  real :: prevMassC         ! Total Substrait from X timesteps ago
 
  real,dimension(n) :: Mnew
  real,dimension(n) :: Cnew
  real,dimension(n,ndiag) :: MatrixM
  real,dimension(n) :: Mrhs
  real,dimension(n) :: Cprev, Mprev
  integer,dimension(ndiag) :: Mioff 
  integer :: stopcritria 
  integer :: stat 
  real :: diffC, diffM 
  integer :: countIters 
  real :: peak, height, intfac  ! Track Interface and wave peak
  integer :: travExist          ! 1 if trav wave exist at this timestep
  real :: waveSpeed             ! calculated wavespeed at current timestpe
  real,dimension(n) :: Mprev_10 ! M from 10-ish timesteps ago, for travCheck

  stored_e1 = e1
  stored_e2 = e2 
 
  filter = int(tEnd/(nOuts*tDel))
  endLoop = 0
  counter = 0
  countIters = 0
  avgIters = 0
  avgNit = 0
  Mprev_10 = M  ! To initiatize for travCheck
  travExist = 0
  
  prevMassC = 1

  open(UNIT = 124, IOSTAT = stat, FILE = "total.dat", STATUS = "old")
  if (stat .EQ. 0) close(124, STATUS = "delete")
  open(UNIT = 120, FILE = "total.dat", POSITION = "append", ACTION = "write")

  open(UNIT = 124, IOSTAT = stat, FILE = "COprod.dat", STATUS = "old")
  if (stat .EQ. 0) close(124, STATUS = "delete")
  open(UNIT = 126, FILE = "COprod.dat", POSITION = "append", ACTION = "write")
  
  if (true2D == 1) then
    open(UNIT = 124, IOSTAT = stat, FILE = "peakInfo.dat", STATUS = "old")
    if (stat .EQ. 0) close(124, STATUS = "delete")
    open(UNIT = 121, FILE = "peakInfo.dat", POSITION = "append", ACTION = "write")
  endif    

  if (checkTrav == 1) then
    open(UNIT = 124, IOSTAT = stat, FILE = "travCheck.dat", STATUS = "old")
    if (stat .EQ. 0) close(124, STATUS = "delete")
    open(UNIT = 138, FILE = "travCheck.dat", POSITION = "append", ACTION = "write")
  endif
  
  write(*,*) "   time    avgIter  maxIter      avgNit  maxNit        avgM        avgC"

  do while((counter) * tDel <= tEnd)
    ! Output every 100 times more then nOuts for the peak info
    if (MOD(counter, int(filter/100)+1) == 0) then
      ! Get total M and C
      call calcMass(M,totalMassM,n,row,col)
      call calcMass(C,totalMassC,n,row,col)
      write(120,*) counter*tDel, totalMassM, totalMassC

      ! CO_2 production: t, current produced CO2, total produced CO2
      write(126,*) counter*tDel, prevMassC - totalMassC, 1 - totalMassC
      prevMassC = totalMassC

      ! peak-interface, only availbale for 2D graphs
      if (true2D == 1) then
        call calcPeakInterface(M, row, col, peak, height, intfac)
        write (121,*) tDel*counter, peak, height, intfac
      end if
    endif 
    if (true2D == 1 .AND. checkTrav == 1 .AND. MOD(counter, int(filter))==0) then
      call checkTravWave(M, Mprev_10, row, col, travExist, wavespeed, height, eTrav)
      wavespeed = wavespeed/filter ! Makes wavespeed independent of number of outputs
      write (138,*) tDel*counter, travExist, wavespeed
      Mprev_10 = M
    endif
  
    ! Write to file / report Total Mass
    if (MOD(counter, filter) == 0) then
      if (true2D == 1) then
        call printToFile2D(n,row,col,M,C)
      else
        call printToFile(n,row,col,M,C)
      end if

      if (counter == 0) then
        write(*,'(F8.2,F12.2,I8,F12.2,I8,F12.6,F12.6)') 0.0, 0.0, &
          int(maxIters), 0.0, int(maxNit), totalMassM, totalMassC
      else
        write(*,'(F8.2,F12.2,I8,F12.2,I8,F12.6,F12.6)') tDel*counter, &
          real(avgIters/counter), int(maxIters), real(avgNit/avgIters), &
          int(maxNit),totalMassM, totalMassC
      endif
    endif

    diffC = 1
    diffM = 1
    countIters = 0
    Cprev = C
    Mprev = M

    do while(diffC + diffM > eSoln)
        nit = 100*n
        e1 = stored_e1 
        e2 = stored_e2
        
        ! Solve M
        call GenMatrixM(Mprev,C,MatrixM,Mioff,Mrhs,row,col,n,ndiag,delta,nu,&
                    alpha,beta,kappa,tDel,xDel,dSelect,fSelect)
        call solveLSDIAG(n,ndiag,Mioff,MatrixM,Mnew,Mrhs,nit,e1,e2,stopcritria)
  
        ! Solve C
        call solveC(Mprev,Mnew,Cprev,Cnew,n,kappa,gama,tDel,gSelect)

        ! Solve CO_2 production
        

        if(nit > maxNit) maxNit = nit
        avgNit = avgNit + nit

        call calcDiff(diffC, C, Cnew, row, col)
        call calcDiff(diffM, M, Mnew, row, col)

        C = Cnew
        M = Mnew
        countIters = countIters+1

        if (countIters > 10000) then
          write(*,*) "[!] Over 10000 iterations in one timestep"
          write(*,*) "[!] Solutions not converging. Exit!"
          stop
        end if
!        write(*,*) countIters, diffC, diffM, C(12),M(12)
    enddo
    if(countIters > maxIters) maxIters = countIters
    avgIters = avgIters + countIters
    
    counter = counter + 1
  enddo
  avgNit = avgNit/(avgIters) ! avgIters right now is the total
  avgIters = avgIters/counter

!  ! Print final solution
!  if (true2D == 1) then
!    call printToFile2D(n,row,col,M,C)
!  else
!    call printToFile(n,row,col,M,C)
!  end if
!  write(*,'(F8.2,F12.2,I8,F12.2,I8,F12.6,F12.6)') tDel*counter, &
!          real(avgIters), int(maxIters), real(avgNit), &
!          int(maxNit),totalMassM, totalMassC

  close(120)
  close(138)
  close(126)
  if (true2D == 1) then
    close(121)
  endif
  
end subroutine solveOrder2


!==============================================================================
!   Calculate the Total Mass
!-----------------------------------------------------------------------------
!   Uses Riemann Sums kind of concept to check if the program runs correctly
!   since the total mass of the region can be computed and checked.
!==============================================================================
subroutine calcMass(X,totalMass,n,row,col)
    implicit none
    integer,intent(in) :: n,row,col
    real,dimension(n),intent(in) :: X
  
    real,intent(out) :: totalMass
  
    integer :: i,j,g
  
    totalMass = 0

    !$omp parallel do reduction(+:totalMass) 
    do i=1,n
      totalMass = totalMass + X(i)
    enddo
    !$omp end parallel do 
    
    totalMass = totalMass / n
  
end subroutine calcMass


!==============================================================================
!   Prints out the solution in a format accepted by gnuplot, used for graphing
!-----------------------------------------------------------------------------
!   Runs through the grid, row-by-row. The MOD and filter act to reduce the 
!     number of grid points written, useful when comparing different grid 
!     sizes.
!-----------------------------------------------------------------------------
!   Input:
!     n      =  The problem size
!     row    =  The number of rows
!     col    =  The number of columns
!     M      =  The solution vector for biomass
!     C      =  The solution for substrait
!   Output:
!     none
!==============================================================================
subroutine printToFile(n,row,col,M,C)
    implicit none

    integer,intent(in) :: n, row, col
    real,dimension(n),intent(in) :: M,C

    integer :: p
    integer :: i,j
    integer :: stat
    integer :: filter
  
    !-------------------------------------------
    ! Deletes the old output file if it exist
    !-------------------------------------------
    open(UNIT = 123, IOSTAT = stat, FILE = "output.dat", STATUS = "old")
    if (stat .EQ. 0) close(123, STATUS = "delete")

    open(UNIT = 11, FILE = "output.dat", POSITION = "append", ACTION = "write")

    filter = 1
    if (row .gt. 257) then
      filter = (row-1)/256
    endif

    do i=1,row
      if (MOD(i-1, filter) == 0) then
        do j=1,col
          if (MOD(j-1, filter) == 0) then
            p = (j + (i-1)*col)
            write(11,*) real(j-1)/real(col-1), &
                        real(i-1)/real(row-1), M(p), C(p)
          endif
        enddo
        write(11,*) ' '
      endif
    enddo
    write(11,*) 
    
      
end subroutine printToFile


!==============================================================================
!   Prints out the solution in a 2D format, used for graphing the Travelling
!     Wave Example.
!-----------------------------------------------------------------------------
!   Runs through the grid, row-by-row. The MOD and filter act to reduce the 
!     number of grid points written
!   Unique here is that the average of the x-axis is taken so that the system 
!     can be reduced to just y. Also written are the max and min for each y 
!     value; this is used for showing that the system can be reduced.
!-----------------------------------------------------------------------------
!   Input:
!     n      =  The problem size
!     row    =  The number of rows
!     col    =  The number of columns
!     M      =  The solution vector for biomass
!     C      =  The solution for substrait
!   Output:
!     none
!==============================================================================
subroutine printToFile2D(n,row,col,M,C)
   implicit none

   integer,intent(in) :: n, row, col
   real,dimension(n),intent(in) :: M,C

   integer :: p
   integer :: i,j
   integer :: stat
   integer :: filter
   real :: averageM, averageC
   real :: maxM, minM, maxC, minC
   real :: y
  
   !-------------------------------------------
   ! Deletes the old output file if it exist
   !-------------------------------------------
   open(UNIT = 124, IOSTAT = stat, FILE = "2D_output.dat", STATUS = "old")
   if (stat .EQ. 0) close(124, STATUS = "delete")

   open(UNIT = 12, FILE = "2D_output.dat", POSITION = "append", ACTION = "write")

   filter = 1
   if (row .gt. 257) then
     filter = (row-1)/256
   endif
 
   do, i=1,row
     if (MOD(i-1, filter) == 0) then
       averageM = 0
       averageC = 0
       maxM = 0
       maxC = 0
       minM = 1
       minC = 1
       do, j=1,col
           p = j + (i-1)*col
           averageM = averageM + M(p)
           averageC = averageC + C(p)
           if(M(p) .ge. maxM) then 
               maxM = M(p)
           endif
           if(M(p) .le. minM) then 
               minM = M(p)
           endif
           if(C(p) .ge. maxC) then
               maxC = C(p)
           endif
           if(C(p) .le. minC) then
               minC = C(p)
           endif
       enddo
       averageM = averageM/(col)
       averageC = averageC/(col)
       y = real(i-1)/real(row-1)
       write(12,'(f20.12,f20.12,f20.12)') y,averageM,averageC
!wrte(12,'(f14.10,f14.10,f14.10,f14.10,f14.10,f14.10,f14.10)') y, averageM, averageC, minM, maxM, minC, maxC
     endif
   enddo
   write(12,*) 
      
      
end subroutine printToFile2D


!==============================================================================
!   Calculates the difference between each grid point for the solutions at
!       different iterations
!-----------------------------------------------------------------------------
!   Input:
!     row    =  The number of rows
!     col    =  The number of columns
!     C      =  The previous solution for substrait
!     Cnew   =  The current solution for substrait
!   Output:
!     diff   =  The average difference between the two C's
!==============================================================================
subroutine calcDiff(diff, C, Cnew, row, col)
    integer, intent(in) :: row, col
    real, dimension(row*col), intent(in) :: C, Cnew
    real, intent(out) :: diff

    integer :: i

    diff = 0
    !$omp parallel do reduction(+:diff) 
    do i=1,row*col
        diff = diff + abs(C(i) - Cnew(i))
    enddo
    !$omp end parallel do
    diff = diff/real(row*col)
    
end subroutine calcDiff


!==============================================================================
!   Calculates the peak and interface info at a single timestep
!-----------------------------------------------------------------------------
!   Input:
!     row    =  The number of rows
!     col    =  The number of columns
!     M      =  The solution for biomass, array size (n)
!   Output:
!     peak   =  Peak location
!     height =  Peak height
!     intfac =  Interface location
!==============================================================================
subroutine calcPeakInterface(M, row, col, peak, height, intfac)
    implicit none
    integer, intent(in):: row,col
    real, dimension(row*col), intent(in) :: M
    real, intent(out) :: peak, height, intfac

    real :: hei
    real :: y
    integer :: i,j,p

    hei = 0
    do, i=1,row
      y = real(i-1)/real(row-1)
      do, j=1,col
          p = (j + (i-1)*col)
          if (M(p) >= hei) then 
              peak = y
              hei = M(p)
          endif
          if (M(p) > 0.1) then
              intfac = y
          endif
      enddo
    enddo
    height = hei

end subroutine 



!=======================================================================
!   Check if there is evidence of a travelling wave solution
!----------------------------------------------------------------------
!   Makes an educated guess for the wavespeed (which would be incorrect
!   by, at worst, 2*eTrav) and then uses this wavespeed to check if the
!   difference between M and Mprev is consistently wavespeed +- eTrav
!     Reports 1 or 0 for travExists based on if traveilling exists (1)
!   or not (0). Also reports the approximated wavespeed.
!=======================================================================
subroutine checkTravWave(M, Mprev, row, col, travExist, wavespeed, height, eTrav)
  implicit none
  integer, intent(in) :: row, col
  real, intent(in) :: height, eTrav
  real, dimension(row * col), intent(in) :: M, Mprev
  
  integer, intent(out) :: travExist
  real, intent(out) :: wavespeed

  integer :: i, wavePoint1, wavePoint2
  real :: diff  ! placeholder for difference between M and Mprev
  wavePoint1 = -1 
  wavePoint2 = -1 

  do, i=row,1,-1
    ! -3 because each column is 4 and I want to start at 1
    if (M(i*col-3) > 0.09 - eTrav .AND. M(i*col-3) < 0.09 + eTrav) then
      wavePoint1 = i 
      exit
    endif
  enddo

  do, i=row,1,-1
    if (Mprev(i*col-3) > 0.09 - eTrav .AND. Mprev(i*col-3) < 0.09 + eTrav) then
      wavePoint2 = i
      exit
    endif
  enddo

  ! Here the wavespeed is in 'i' units
  wavespeed = abs(wavePoint1 - wavePoint2)    

  
  travExist = 1
  do i = 1, row-int(wavespeed)
   !write(*,*) wavespeed, i,  M((i - int(wavespeed))* col - 3 ), Mprev(i*col - 3)
    diff = abs(M((i+int(wavespeed)) * col - 3) - Mprev(i * col - 3))
    if ( diff > eTrav*10 ) travExist = 0
  enddo 
  if (wavespeed == 0) travExist = 0 ! Can't have trav wave with 0 speed

  ! Here wavespeed is converted to dimensionless units over X time
  wavespeed = wavespeed / float(row)


end subroutine checkTravWave



!==============================================================================
!   Writes a bunch of statistics to file
!-----------------------------------------------------------------------------
!   Input:
!       avgIters = average number of iterations from between solutions
!       maxIters = maximum number of iterations from between solutions
!       avgNit   = average number of iterations from linear solver
!       maxNit   = maximum number of iterations from linear solver
!       time     = time to complete solveOrder
!   Output: 
!       write everything to the file.
!==============================================================================
subroutine reportStats(avgIters,maxIters,avgNit,maxNit,time)
  implicit none
  real,intent(in)::avgIters,maxIters,avgNit,maxNit
  real,intent(in)::time
  integer :: stat

  open(UNIT = 125, IOSTAT = stat, FILE = "statReport.dat", STATUS = "old")
  if (stat .EQ. 0) close(125, STATUS = "delete")
  open(UNIT = 128, FILE = "statReport.dat", POSITION = "append", ACTION = "write")
  
  write(128,*) "Statsitcs:"
  write(128,*) "-------------------------------------------------------------"
  write(128,*) "Time to compute = ", time
  write(128,*) ""
  write(128,*) "Avg Iters for iterating betn. soln. =", avgIters
  write(128,*) "Max Iters for iterating betn. soln. =", maxIters
  write(128,*) "Avg Iters for linear solver =", avgNit
  write(128,*) "Max Iters for linear solver =", maxNit  
  
  close(128)
  
end subroutine



subroutine amuxd (n,x,y,diag,idiag,ioff) 
!-----------------------------------------------------------------------
!        Mnew times a vector in Diagonal storage format (DIA) 
!        f90/f95 version of the sparskit f77 subroutine
!----------------------------------------------------------------------- 
! multiplies a matrix by a vector when the original matrix is stored 
! in the diagonal storage format.
!-----------------------------------------------------------------------
!
! on entry:
!----------
! n     = row dimension of Mnew
! x     = real array of length equal to the column dimension of
!         the Mnew matrix.
! ndiag  = integer. The first dimension of array adiag as declared in
!         the calling program.
!         (obscolete (=n always)
! idiag  = integer. The number of diagonals in the matrix.
! diag   = real array containing the diagonals stored of Mnew.
! idiag  = number of diagonals in matrix.
! diag   = real array of size (ndiag x idiag) containing the diagonals
!          
! ioff   = integer array of length idiag, containing the offsets of the
!   	   diagonals of the matrix:
!          diag(i,j) contains the element a(i,i+ioff(j)) of the matrix.
!
! on return:
!-----------
! y     = real array of length n, containing the product y=Mnew*x
!
!-----------------------------------------------------------------------
implicit none
    integer, intent(in)::  n, idiag
    integer, intent(in),dimension(idiag) :: ioff
    real, dimension(n), intent(in) :: x    
    real, dimension(n,idiag), intent(in) :: diag
    real, dimension(n), intent(out) :: y
    integer :: j, io, i1, i2, i       

    !$omp parallel shared(y,diag,x,n) private(j,io,i1,i2)

    !!$omp workshare
    !$omp do
    do i=1,n
        y(i)=0.  
    enddo
    !$omp enddo
    !!$omp end workshare    

    do j=1, idiag
        io = ioff(j)
        i1 = max0(1,1-io)
        i2 = min0(n,n-io)
        !$omp do
        do i=i1,i2
            y(i) = y(i)+diag(i,j)*x(i+io)
        enddo
        !$omp end do
    enddo  
    !$omp end parallel
    
end subroutine amuxd





