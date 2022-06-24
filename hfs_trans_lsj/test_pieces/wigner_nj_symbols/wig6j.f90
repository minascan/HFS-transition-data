program sixj_symbol
  implicit none

  integer :: step
  double precision :: j1, j2, j3, l1, l2, l3, a, b, c, d1, d2, d3, d4, w6j
  double precision :: k, kmax, kmin, numerator, denominator
  double precision, parameter :: modtwo = 2

  print *, 'Give the values of the six symbols '
  read(*,'(f3.1)') j1
  read(*,'(f3.1)') j2 
  read(*,'(f3.1)') j3
  read(*,'(f3.1)') l1
  read(*,'(f3.1)') l2
  read(*,'(f3.1)') l3

  if (j1 .ge. 0 .and. j2 .ge. 0 .and. j3 .ge. 0 .and. l1 .ge. 0 .and.  &
       l2 .ge. 0 .and. l3 .ge. 0 .and.  &                                                              ! quantities .ge. 0
       mod(4*j1,modtwo) .eq. 0 .and. mod(4*j2,modtwo) .eq. 0 .and. mod(4*j3,modtwo) .eq. 0 .and.  &    ! ji,li integer 
       mod(4*l1,modtwo) .eq. 0 .and. mod(4*l2,modtwo) .eq. 0 .and. mod(4*l3,modtwo) .eq. 0 .and.  &    ! or half integer
       mod(2*(j1+j2+j3),modtwo) .eq. 0 .and.                                   &                       ! j1+j2+j3 integer 
       j1 + j2 .ge. j3 .and. j2 + j3 .ge. j1 .and. j3 + j1 .ge. j2 .and.  &                            ! triangle relations
       mod(2*(j1+l2+l3),modtwo) .eq. 0 .and.                                   &                       ! j1+l2+l3 integer 
       j1 + l2 .ge. l3 .and. l2 + l3 .ge. j1 .and. l3 + j1 .ge. l2 .and.  &                            ! triangle relations
       mod(2*(l1+j2+l3),modtwo) .eq. 0 .and.                                   &                       ! l1+j2+l3 integer 
       l1 + j2 .ge. l3 .and. j2 + l3 .ge. l1 .and. l3 + l1 .ge. j2 .and.  &                            ! triangle relations
       mod(2*(l1+l2+j3),modtwo) .eq. 0 .and.                                   &                       ! l1+l2+j3 integer 
       l1 + l2 .ge. j3 .and. l2 + j3 .ge. l1 .and. j3 + l1 .ge. l2 )     then                          ! triangle relations
     
     a = j1
     b = j2
     c = j3
     d1 = sqrt(factorial(a+b-c)*factorial(a-b+c)*factorial(-a+b+c)/factorial(a+b+c+1))
     a = j1
     b = l2
     c = l3
     d2 = sqrt(factorial(a+b-c)*factorial(a-b+c)*factorial(-a+b+c)/factorial(a+b+c+1))
     a = l1
     b = j2
     c = l3
     d3 = sqrt(factorial(a+b-c)*factorial(a-b+c)*factorial(-a+b+c)/factorial(a+b+c+1))
     a = l1
     b = l2
     c = j3
     d4 = sqrt(factorial(a+b-c)*factorial(a-b+c)*factorial(-a+b+c)/factorial(a+b+c+1))
     
     kmin = max(j1+j2+j3, j1+l2+l3, l1+j2+l3, l1+l2+j3)
     kmax = min(j1+j2+l1+l2, j2+j3+l2+l3, j3+j1+l3+l1)
     
     w6j = 0
     step = 0
     k = kmin
     do while (k.le.kmax)
        numerator = factorial(k+1)
        denominator = factorial(k-j1-j2-j3)*factorial(k-j1-l2-l3)*factorial(k-l1-j2-l3)* &
             factorial(k-l1-l2-j3)*factorial(j1+j2+l1+l2-k)*factorial(j2+j3+l2+l3-k)*factorial(j3+j1+l3+l1-k)
        print*, numerator, denominator
        w6j = w6j + (-1)**k*numerator/denominator
        step = step +1
        k = kmin + step
     end do
     w6j = w6j*d1*d2*d3*d4
  else
     w6j = 0
  endif

  print *, w6j

contains

  double precision function factorial(n)
    implicit none
    integer :: i
    double precision, intent(in) :: n
    double precision :: fac
    double precision, parameter :: modtwo = 2
    
    if (n.eq.0 .or. n.eq.1) then
       fac = 1
    elseif (mod(n,modtwo).eq.0 .or. mod(n,modtwo).eq.1) then
       fac = 1
       
       do i = 1, int(n)
          fac = fac*i
       enddo      
    else
       fac = gamma(n)
    endif

    factorial = fac
    
  end function factorial

  
end program sixj_symbol
