program trans_prop

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
   !!!             input file: name1.name2.ct.lsj                           !!!
  
   !!!             output file: hfs.ct.lsj                                  !!!
   !!!                                                                      !!!
   !!!     Written by Asimina Papoulia ,   September 2018                   !!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  implicit none

  double precision, parameter :: nuc_spin = 2.5, l2 = 1.0
  character(100) :: filename1, dummy, conf_lo, conf_up, E_lo, E_up
  character(180) :: string
  integer :: count1, j, i, l_up, l_lo, test, nom1, den1, nom2, den2, k, l, blank
  integer :: twoj_up, twoj_lo
  double precision :: j_up, j_lo, A_c, Ahfs_c, A_b, Ahfs_b, f_up, f_up_max, f_up_min 
  double precision :: f_lo, f_lo_max, f_lo_min, gf_c, gf_b, w6j, gf_hfs_c, gf_hfs_b
  double precision :: uncertainty
  
  count1 = 0
  
  print *, 'Full name of the transition data file'
  read (*, '(a)') filename1
  !print*, 'Give the nuclear spin I'
  !read (*,'(f3.1)') nuc_spin
  
  open(10, file=filename1, status='old', form='formatted', &  !Open the input transition data file
       action='read')
  open(40, file='hfs.ct.lsj', status='old', form='formatted', &   !Open the output hfs transition data file
       action='write', position='append')
!Writing the headers in the output file
  write(40,'(a51)') '         UPPER                     LOWER                     '
  write(40,'(a98)') '    Conf        J    F       Conf        J    F       A (s-1)         gf         log(gf)    +/-    '
!First reading of the input file to count the number of lines
  do
     read(10, '(a)', end=19) string
     count1 = count1 + 1
  end do
19 continue
  rewind(10)
  
  write(*,'(a42,i3)') ' The number of lines in the file is ', count1

  do j = 1, 3
     read(10, '(a)') dummy
  enddo
!------------------------------------------------------------------------------------------------------
! main loop starts here
  do i = 1, (count1-3)/7        !header lines are 3 in the .ct.lsj input file
     ! every loop reads the two blank lines first + 5 more lines that contain all the info
     do j = 1, 2
        read(10, '(a)') dummy
     enddo
     ! and then reads the other five lines of the input file and assign the useful values to variables
     read(10, '(i4, a14, a20, a11)') twoj_lo, E_lo, dummy, conf_lo
     read(10, '(i4, a14, a20, a11)') twoj_up, E_up, dummy, conf_up
     read(10, '(a)') dummy
     read(10, '(a30, d12.7, a9, d12.1, a17)') dummy, gf_b, dummy, A_b, dummy
     read(10, '(a30, d12.7, a8, d13.3)') dummy, gf_c, dummy, A_c
     !print*, E_lo, twoj_lo, conf_lo
     !print*, E_up, twoj_up, conf_up
     !print*, gf_b, A_b
     !print*, gf_c, A_c
     
     j_up = real(twoj_up)/2
     j_lo = real(twoj_lo)/2
     !print*, j_up, j_lo 
     ! estimating the possible F quantum numbers for upper and lower levels
     f_up_max = j_up + nuc_spin
     f_up_min = abs(j_up - nuc_spin)    
     f_lo_max = j_lo + nuc_spin
     f_lo_min = abs(j_lo - nuc_spin)
     !print*, f_up_max, f_up_min
     k = 0  ! step counter of possible F values for the UPPER level !***
     l = 0  ! step counter of possible F values for the LOWER level !***
     f_up = f_up_min
     f_lo = f_lo_min
     
     do while (f_up.le.f_up_max)  ! loop for all possible F values of the UPPER level
        blank = 0  ! when zero there is a new combination of quantum nubers and therefore a
                   ! complete line of variables in the output file  !%%%
        l = 0      ! step counter of possible F values for the LOWER level !***
        do while (f_lo.le.f_lo_max)  ! loop for all possible F values of the LOWER level !%%%
           ! checking whether the transition is allowed i.e. DF = 0 or +/-1
           ! and whether it is a transition for a new F value of the upper level i.e. blank = 0 
           if ((0.le.abs(f_up-f_lo) .and. abs(f_up-f_lo).le.1 .and. blank.eq.0) .and. &
                (f_up.ne.0 .or. f_lo.ne.0)) then
                  
              ! calculate the transition data in the COULOBM gauge for this hfs transition
              w6j = wig6j(j_lo, nuc_spin, f_lo, f_up, l2, j_up)
              call probabilities(f_lo, j_up, w6j, A_c, Ahfs_c)
              call weightedf(f_up, f_lo, j_lo, w6j, gf_c, gf_hfs_c)           
              ! write data in the output file
              if (conf_up.eq.'3s(2).4d_2D' .and. (E_up.eq.' -242.14686040' .or. E_up.eq.' -242.14683776')) then
                 write(40,'(1x,a11,a1,2x,f3.1,2x,f3.1,3x,a11,3x,f3.1,2x,f3.1,2x,a1,1x,d11.5,2x,d12.6)') &
                      conf_up,'b', j_up, f_up, conf_lo, j_lo, f_lo, 'C', Ahfs_c, gf_hfs_c
              elseif (conf_lo.eq.'3s(2).4d_2D' .and. (E_lo.eq.' -242.14686040' .or. E_lo .eq.' -242.14683776')) then
                 write(40,'(1x,a11,3x,f3.1,2x,f3.1,3x,a11,a1,2x,f3.1,2x,f3.1,2x,a1,1x,d11.5,2x,d12.6)') &
                      conf_up, j_up, f_up, conf_lo,'b', j_lo, f_lo, 'C', Ahfs_c, gf_hfs_c
              else
                 write(40,'(1x,a11,3x,f3.1,2x,f3.1,3x,a11,3x,f3.1,2x,f3.1,2x,a1,1x,d11.5,2x,d12.6)') &
                      conf_up, j_up, f_up, conf_lo, j_lo, f_lo, 'C', Ahfs_c, gf_hfs_c
              end if
              print*, Ahfs_c, gf_hfs_c
              ! calculate the transition data in the BABUSHKIN gauge for this hfs transition
              call probabilities(f_lo, j_up, w6j, A_b, Ahfs_b)
              call weightedf(f_up, f_lo, j_lo, w6j, gf_b, gf_hfs_b)
              uncertainty = uncert(real(gf_hfs_c), real(gf_hfs_b))          
              write(40,'(a50,a1,1x,d11.5,2x,d12.6,3x,f10.7,2x,f9.7)') '                                                  ', &
                   'B', Ahfs_b, gf_hfs_b, log10(gf_hfs_b), uncertainty
              print*, Ahfs_b, gf_hfs_b, log10(gf_hfs_b), uncertainty
              blank = blank + 1  !%%%
              ! checking whether the transition is allowed i.e. DF = 0 or +/-1
              ! and whether it is a transition for the same F value of the upper level i.e. blank >=1 
           elseif (0.le.abs(f_up-f_lo) .and. abs(f_up-f_lo).le.1 .and. blank.ne.0) then
              
              ! calculate the transition data in the COULOBM gauge for this hfs transition         
              w6j = wig6j(j_lo, nuc_spin, f_lo, f_up, l2, j_up)
              call probabilities(f_lo, j_up, w6j, A_c, Ahfs_c)
              call weightedf(f_up, f_lo, j_lo, w6j, gf_c, gf_hfs_c)                    
              ! write data in the output file
              write(40, '(a45,f3.1,2x,a1,1x,d11.5,2x,d12.6)') '                                            ', &
                   f_lo, 'C', Ahfs_c, gf_hfs_c
              ! calculate the transition data in the BABUSHKIN gauge for this hfs transition
              call probabilities(f_lo, j_up, w6j, A_b, Ahfs_b)
              call weightedf(f_up, f_lo, j_lo, w6j, gf_b, gf_hfs_b)
              uncertainty = uncert(real(gf_hfs_c), real(gf_hfs_b))
              write(40,'(a50,a1,1x,d11.5,2x,d12.6,3x,f10.7,2x,f9.7)') '                                                  ', &
                    'B', Ahfs_b, gf_hfs_b, log10(gf_hfs_b), uncertainty
              
           endif
           
           l = l + 1           !***
           f_lo = f_lo_min + l !***       
        enddo
        !l = 0
        !blank = 0  !%%%
        f_lo = f_lo_min
        k = k + 1             !***
        f_up = f_up_min + k   !***
     enddo
       
     !print*, j_up, j_lo, A_c, A_b, gf_c, gf_b       !test print
     
  enddo
!------------------------------------------------------------------------------------------------------  
write(40,'(a50)') '--------------------------------------------------'


contains
!-----------------------------------------------------------------------------------------------------
  subroutine probabilities(f_lower, j_upper, w6j, prop, outp)

    implicit none
    double precision, intent(in) :: f_lower, j_upper, w6j, prop
    double precision :: brack
    double precision, intent(out) :: outp

    brack = (2*f_lower + 1)*(2*j_upper + 1) 
    outp = brack*(w6j**2)*prop
    
  end subroutine probabilities
!-----------------------------------------------------------------------------------------------------
  subroutine weightedf(f_upper, f_lower, j_lower, w6j, val, outc)

    implicit none
    double precision, intent(in) :: f_upper, f_lower, j_lower, w6j, val
    double precision :: brack, weight
    double precision, intent(out) :: outc

    brack = (2*f_upper + 1)*(2*j_lower + 1)
    weight = (2*f_lower + 1)/(2*j_lower + 1)
    outc = weight*brack*(w6j**2)*val
    
  end subroutine weightedf
!-----------------------------------------------------------------------------------------------------
  double precision function uncert(gf_hfs_coul, gf_hfs_babus)
    implicit none
    real, intent(in) :: gf_hfs_coul, gf_hfs_babus
    double precision :: delta_gf, delta_log_gf

    delta_gf = gf_hfs_babus*abs(gf_hfs_coul - gf_hfs_babus)/max(gf_hfs_coul,gf_hfs_babus)
    delta_log_gf = ( abs( log10(gf_hfs_babus + delta_gf) - &
         log10(gf_hfs_babus - delta_gf) ) )/ 2d0

    uncert = delta_log_gf
    
  end function uncert
!-----------------------------------------------------------------------------------------------------  
  double precision function wig6j(j1,j2,j3,l1,l2,l3)

    !Computes the Wigner 6j-symbol {j1 j2 j3}
    !                              {l1 l2 l3}   

    implicit none
    double precision, intent(in) :: j1,j2,j3,l1,l2,l3
    double precision :: a, b, c, d1, d2, d3, d4, w6j
    double precision :: k, kmax, kmin, numerator, denominator
    double precision, parameter :: modtwo = 2 !l2=1
    integer :: step
    
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
             factorial(k-l1-l2-j3)*factorial(j1+j2+l1+l2-k)*factorial(j2+j3+l2+l3-k)* &
             factorial(j3+j1+l3+l1-k)
        w6j = w6j + (-1)**k*numerator/denominator
        step = step +1
        k = kmin + step
     end do
     w6j = w6j*d1*d2*d3*d4
  else
     w6j = 0
  endif

  wig6j = w6j
  
  end function wig6j
!-----------------------------------------------------------------------------------------------------  
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
!-----------------------------------------------------------------------------------------------------
end program trans_prop

