subroutine solveLSDIAG(n,ndiag,ioff,A,sol,rhs,nit,err1,err2,stopcrit)
!----------------------------------------------------------------------------
! HJE FOR M6020
!----------------------------------------------------------------------------
! uses the conjugate gradient method to solve Ax=b in diagonal format
!----------------------------------------------------------------------------
! input:  n       problem size
!         ndiag:  number of diagonals
!         ioff:   offsets (distance of sub diagonals to main diagonal)
!         A:      matrix values
!         sol:    initial guess for iteration (will be overwritten with result)
!         rhs:    the righ hand side of the linear system
!         nit:    maximum number of iterations to be carried out (overwritten)
!         err1:   tolerance for 1st stopping criterion (will be overwritten)
!         err2:   tolerance for 2nd stopping criterion (will be overwritten)
! output: sol:    solution of Ax=b
!         nit:    number of iterations taken
!         err1:   computed value for 1st stopping criterion
!         err2:   computed value for 2nd stopping criterion
!         stopcrit:  int value that depict which stop criterions became active
!----------------------------------------------------------------------------
implicit none
 integer, intent(in)               :: n,ndiag
 real,dimension(n,ndiag),intent(in):: a
 integer, dimension(ndiag),intent(in)::ioff
 real,dimension(n),intent(in)      :: rhs
 real,dimension(n),intent(inout)   :: sol
 real,intent(inout)                :: err1,err2
 integer,intent(inout)             :: nit
 integer,intent(out)               :: stopcrit

 real,dimension(n) :: r, z, p, q
 real :: rho, beta, rhoold, alfa, dot

 integer           :: maxit,i
 real              :: tol1,tol2,normx,normA,normb

 ! initialization: set tolerances, max number of iterations
 ! and rough estimates of matrix and rhs norms for stopping criteria
 maxit=nit; tol1=err1; tol2=err2
 normA=maxval(abs(a)); normb=maxval(abs(rhs))


 rho=0.

 call amuxd(n,sol,r,a,ndiag,ioff)
 !$omp parallel do
 do i=1,n
    r(i)=rhs(i)-r(i)
 enddo
 !$omp end parallel do

 nit=0; stopcrit=0
 do while(stopcrit==0)
    nit=nit+1
    rhoold=rho
    normx=sqrt(sum(sol*sol))/n

    call dotProd(n,r,r,rho)
    
    if (nit==1) then
        !$omp parallel do shared (p,r)
        do i=1,n
            p(i)=r(i)
        enddo
        !$omp end parallel do
    else
        beta=rho/rhoold
        !$omp parallel do shared (p,r,beta)
        do i=1,n
            p(i)=r(i)+beta*p(i)        
        enddo
        !$omp end parallel do
    endif
    call amuxd(n,p,q,a,ndiag,ioff)
    call dotProd(n,p,q,dot)
    alfa=rho/dot
    !$omp parallel do shared (sol,r,alfa,p) 
    do i = 1, n
        sol(i)=sol(i)+alfa*p(i)
        r(i)=r(i)-alfa*q(i)
        ! test for convergence
    enddo
    !$omp end parallel do

    err1=sqrt(sum(r*r))/n
    err2=abs(alfa)*sqrt(sum(p*p))/n
    if (nit>maxit) stopcrit=-1
    !if (err1<tol1*(norma*normx+normb)) stopcrit=stopcrit-10
    if (err2<tol2*normx) stopcrit=stopcrit-100

    ! uncomment the next line to monitor convergence progress
    !write(*,'(I6,4(E14.7,X))') nit,err1,err2,maxval(sol),minval(sol)
 enddo
end subroutine



subroutine dotProd(n,u,v,sol)
implicit none
    integer, intent(in) :: n
    real, dimension(n), intent(in) :: u,v
    integer :: i

    real, intent(out) :: sol

    sol = 0.
    !$omp parallel do reduction(+:sol)
    do i=1,n
        sol = sol + u(i) * v(i)
    enddo
    !$omp end parallel do

end subroutine







