program total_prob_and_gf
  implicit none
  double precision :: j_lo, f_lo, sum, total_gf, total_prob, sum2
  integer :: i, n, j, m
  double precision, dimension(:), allocatable :: gf_hfs, a
  
  print*, ' Give the angular momentum J of the lower level'
  read(*,'(f3.1)') j_lo
  print*, ' Give the F value of the lower level'
  read(*,'(f3.1)') f_lo

  print*, ' Number of gf values'
  read(*,'(i2)') n
  print*, ' Number of A values'
  read(*,'(i2)') m
  allocate(gf_hfs(n),a(m))

  sum = 0
  sum2 = 0
  
  do i=1,n
     print*, ' Give gf value ', i
     read(*,'(f8.5)') gf_hfs(i)
     sum = sum + gf_hfs(i)
      
  enddo
  
  do j = 1, m
     print*, ' Give A value ', j
     read(*,'(d11.5)') a(j)
     sum2 = sum2 + a(j)
  enddo

  
  total_gf = sum*(2*j_lo + 1)/(2*f_lo + 1)
  
  print*, total_gf, sum2
  
end program total_prob_and_gf
