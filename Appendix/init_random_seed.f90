! initialize a random seed from the system clock at every run (fortran
! 95 code

subroutine init_random_seed
  integer :: i, n, clock
  integer, dimension(:), allocatable :: seed

  call random_seed(size = n)
  allocate(seed(n))

  call system_clock(count = clock)

  seed = clock + 37 * (/ (i -1, i = 1, n) /)
  call random_seed(put = seed)

  deallocate(seed)
end subroutine 

