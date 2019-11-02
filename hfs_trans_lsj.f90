program hyperfine_transition_data

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
   !!!             input files: name1.name2.ct.lsj                          !!!
   !!!                                 name1.chlsj                          !!!
   !!!                                 name2.chlsj                          !!!
  
   !!!             output file: hfs.ct.lsj                                  !!!
   !!!                                                                      !!!
   !!!     Written by Asimina Papoulia ,   September 2018                   !!!
   !!!                     Last Update ,    November 2019                   !!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  implicit none

  integer, parameter :: dp = selected_real_kind(15, 307)
  real(kind=dp),parameter :: MHz_to_invcm = 3.3356410D-05, au_to_invcm = 21947463068D-06
  real(kind=dp),parameter :: l2 = 1.0  ! Wigner 6j symbol parameter
  
  character(100) :: filename1, filename2, filename3, dummy, conf_lo, conf_up, E_lo, E_up
  character(100) :: Bab_data
  character(20), dimension(:), allocatable :: E_levels

  integer :: i, j, levels, count1, count2, count3, l_up, l_lo, test, nom1, den1, nom2, den2 
  integer :: m, ll, uu, k, l, blank, twoj_up, twoj_lo
  
  double precision, dimension(:,:), allocatable :: hfs_results, ab_constants 
  double precision :: nuc_spin, j_up, j_lo, A_b, Ahfs_b, A_c, Ahfs_c, f_up, gf_b, gf_c, gf_hfs_c   
  double precision :: gf_hfs_b, f_up_max, f_up_min, f_lo, f_lo_max, f_lo_min, w6j, uncertainty
  double precision :: E_l, E_u, A_lo, A_up, B_lo, B_up

  logical :: low_level, upper_level
  
  print *, 'Full name of the transition data file'
  read (*, '(a)') filename1
  print *, 'Full name of the first hfs data file'
  read (*, '(a)') filename2
  print *, 'Full name of the second hfs data file'
  read (*, '(a)') filename3
  print *, 'Total number for the levels of the computed A and B hfs constants'
  read (*, '(i3)') levels      ! needed for ALLOCATION of E_levels, hfs_results, ab_constants
  print*, 'Give the nuclear spin I'
  read (*,'(f3.1)') nuc_spin   ! For 27Al is 2.5

  !print*, 'NEW !!!'
  allocate(E_levels(levels), hfs_results(levels,3), ab_constants(levels,3))
  
  !Open and First Reading of the input files to count the number of lines
  call line_counting(10,filename1,count1)
  call line_counting(20,filename2,count2)
  call line_counting(30,filename3,count3)
  !================================================================================================================
  !================================================================================================================
  !Actual Reading of the files containing the hfs constants. Assign them to variables, convert them to cm-1
  !and save them into the 2-D matrix called 'ab_constants'
  !--------------- FILE 1 - EVEN states -----------------
  do j = 1, 6
     read(20, '(a)') dummy ! headers and blank lines
  end do
  do i = 1, count2-6
     read(20, '(d14.6, a52, d12.7, d16.7, a16)') hfs_results(i,1), dummy, hfs_results(i,2), hfs_results(i,3), dummy
     ab_constants(i,1) = hfs_results(i,1)*au_to_invcm
     ab_constants(i,2) = hfs_results(i,2)*MHz_to_invcm
     ab_constants(i,3) = hfs_results(i,3)*MHz_to_invcm
     !print*, hfs_results(i,1), hfs_results(i,2), hfs_results(i,3)
     print*, ab_constants(i,1), ab_constants(i,2), ab_constants(i,3)
  enddo
  rewind(20)
  do j = 1, 6
     read(20, '(a)') dummy
  end do
  do i = 1, count2-6
     read(20, '(a2, a12, a96)') dummy, E_levels(i), dummy
     print*, E_levels(i)
  end do
  
  !--------------- FILE 2 - ODD states -----------------
  do j = 1, 6
     read(30, '(a)') dummy
  end do
  do i = count2-6+1, count2+count3-12
     read(30, '(d14.6, a46, d12.7, d16.7, a16)') hfs_results(i,1), dummy, hfs_results(i,2), hfs_results(i,3), dummy
     ab_constants(i,1) = hfs_results(i,1)*au_to_invcm
     ab_constants(i,2) = hfs_results(i,2)*MHz_to_invcm
     ab_constants(i,3) = hfs_results(i,3)*MHz_to_invcm
     !print *, hfs_results(i,1), hfs_results(i,2), hfs_results(i,3)
     print*, ab_constants(i,1), ab_constants(i,2), ab_constants(i,3)
  enddo
  rewind(30)
  do j = 1, 6
     read(30, '(a)') dummy
  end do
  do i = count2-6+1, count2+count3-12
     read(30, '(a2, a12, a96)') dummy, E_levels(i), dummy
     print*, E_levels(i)
  end do
  !================================================================================================================
  !================================================================================================================
  !writing the headers of the output file
  open(40, file='hfs.ct.lsj', status='old', form='formatted', &   !Open the output hfs transition data file
       action='write', position='append')
  !Writing the headers in the output file
  write(40,'(a51)') '         UPPER                     LOWER                     '
  write(40,'(a113)') '    Conf        J    F       Conf        J    F     E (cm-1)     A (s-1)         gf          log(gf) &
       &   +/-      '

  !writing the rest of the data is done simultaneously with the reading of the transition ct.lsj file
  do j = 1, 3
     read(10, '(a)') dummy
  enddo
!------------------------------------------------------------------------------------------------------
! main loop starts here
  do i = 1, (count1-3)/7        !header lines are 3 in the .ct.lsj input file
     low_level = .false.
     upper_level = .false.
     ! every loop reads the two blank lines and the next 5 contain the data for each transition 
     do j = 1, 2
        read(10, '(a)') dummy
     enddo
     ! and then reads the other five lines of the input file and assign the useful values to variables
     read(10, '(i4, a12, a22, a11)') twoj_lo, E_lo, dummy, conf_lo
     read(10, '(i4, a12, a22, a11)') twoj_up, E_up, dummy, conf_up
     read(10, '(a)') dummy
     read(10, '(a30, d12.7, a9, d12.1, a17)') dummy, gf_b, dummy, A_b, dummy
     read(10, '(a30, d12.7, a8, d13.3)') dummy, gf_c, dummy, A_c
     print*, E_lo, twoj_lo, conf_lo
     print*, E_up, twoj_up, conf_up
     !print*, gf_b, A_b
     !print*, gf_c, A_c
     !-----------------------------------------------------------------------------------------
     ! based on the E_lo & E_up levels we get the corresponding values of A and B hfs constants
     do m = 1, levels
        if (E_lo.eq.E_levels(m)) then
           ll = m
           low_level = .true.
           print*, 'The lower level has been identified', ll
        elseif (E_up.eq.E_levels(m)) then
           uu = m
           upper_level = .true.
           print*, 'The upper level has been identified', uu
        endif
        if (low_level .and. upper_level) then
           print*, 'Both levels have been identified', m
           exit
        end if     
     end do
     E_l  = ab_constants(ll,1)
     A_lo = ab_constants(ll,2)
     B_lo = ab_constants(ll,3)
     E_u  = ab_constants(uu,1)
     A_up = ab_constants(uu,2)
     B_up = ab_constants(uu,3)
     !-----------------------------------------------------------------------------------------
     
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

     !The strings of the format for writing data into the output file
     Bab_data = '(a61,a1,1x,d11.5,2x,d12.6,3x,f8.5,2x,f7.5)'
     
     do while (f_up.le.f_up_max)  ! loop for all possible F values of the UPPER level
        blank = 0  ! when zero there is a new combination of quantum nubers and therefore a
                   ! complete line of variables in the output file  !%%%
        l = 0      ! step counter of possible F values for the LOWER level !***
        do while (f_lo.le.f_lo_max)  ! loop for all possible F values of the LOWER level !%%%
           ! checking whether the transition is allowed i.e. DF = 0 or +/-1
           ! and whether it is a transition for a new F value of the upper level i.e. blank = 0 
           if ((0.le.abs(f_up-f_lo) .and. abs(f_up-f_lo).le.1 .and. blank.eq.0) .and. &
                (f_up.ne.0 .or. f_lo.ne.0)) then
              !===============NEW======================
              ! calculate the hfs transition energies
              !call hfs_trans_energies()
              ! calculate the transition data in the COULOBM gauge for this hfs transition
              w6j = wig6j(j_lo, nuc_spin, f_lo, f_up, l2, j_up)
              call probabilities(f_lo, j_up, w6j, A_c, Ahfs_c)
              call weightedf(f_up, f_lo, j_lo, w6j, gf_c, gf_hfs_c)           
              ! write data in the output file
              !---------------------------------------------------------------------------------------------------------------------------------
              ! special case for the Al_I computations since two states have the same label and we need to distinguish the second one with a "b"
              ! this state can be either the lower or the upper and thus we need to check that for both cases
              if (conf_up.eq.'3s(2).4d_2D' .and. (E_up.eq.' -242.146860' .or. E_up.eq.' -242.146837')) then ! we need to check for both j=3/2,5/2 
                 write(40,'(1x,a11,a1,2x,f3.1,2x,f3.1,3x,a11,3x,f3.1,2x,f3.1,2x,f9.2,2x,a1,1x,d11.5,2x,d12.6)') &
                      conf_up,'b', j_up, f_up, conf_lo, j_lo, f_lo, hfs_results(1,1), 'C', Ahfs_c, gf_hfs_c
              elseif (conf_lo.eq.'3s(2).4d_2D' .and. (E_lo.eq.' -242.146860' .or. E_lo .eq.' -242.146837')) then
                 write(40,'(1x,a11,3x,f3.1,2x,f3.1,3x,a11,a1,2x,f3.1,2x,f3.1,2x,f9.2,2x,a1,1x,d11.5,2x,d12.6)') &
                      conf_up, j_up, f_up, conf_lo,'b', j_lo, f_lo, hfs_results(1,1), 'C', Ahfs_c, gf_hfs_c
              !---------------------------------------------------------------------------------------------------------------------------------   
              else
                 write(40,'(1x,a11,3x,f3.1,2x,f3.1,3x,a11,3x,f3.1,2x,f3.1,2x,f9.2,2x,a1,1x,d11.5,2x,d12.6)') &
                      conf_up, j_up, f_up, conf_lo, j_lo, f_lo, hfs_results(1,1), 'C', Ahfs_c, gf_hfs_c
              end if
              !print*, Ahfs_c, gf_hfs_c
              ! calculate the transition data in the BABUSHKIN gauge for this hfs transition
              call probabilities(f_lo, j_up, w6j, A_b, Ahfs_b)
              call weightedf(f_up, f_lo, j_lo, w6j, gf_b, gf_hfs_b)
              uncertainty = uncert(gf_hfs_c, gf_hfs_b)          
              write(40, FMT=Bab_data) '                                                  ', &
                   'B', Ahfs_b, gf_hfs_b, log10(gf_hfs_b), uncertainty
              !print*, Ahfs_b, gf_hfs_b, log10(gf_hfs_b), uncertainty
              blank = blank + 1  !%%%
              ! checking whether the transition is allowed i.e. DF = 0 or +/-1
              ! and whether it is a transition for the same F value of the upper level i.e. blank >=1 
           elseif (0.le.abs(f_up-f_lo) .and. abs(f_up-f_lo).le.1 .and. blank.ne.0) then
              !===============NEW======================
              ! calculate the hfs transition energies
              !call hfs_trans_energies()
              ! calculate the transition data in the COULOBM gauge for this hfs transition         
              w6j = wig6j(j_lo, nuc_spin, f_lo, f_up, l2, j_up)
              call probabilities(f_lo, j_up, w6j, A_c, Ahfs_c)
              call weightedf(f_up, f_lo, j_lo, w6j, gf_c, gf_hfs_c)                    
              ! write data in the output file
              write(40, '(a45,f3.1,2x,f9.2,2x,a1,1x,d11.5,2x,d12.6)') '                                            ', &
                   f_lo, hfs_results(1,1), 'C', Ahfs_c, gf_hfs_c
              ! calculate the transition data in the BABUSHKIN gauge for this hfs transition
              call probabilities(f_lo, j_up, w6j, A_b, Ahfs_b)
              call weightedf(f_up, f_lo, j_lo, w6j, gf_b, gf_hfs_b)
              uncertainty = uncert(gf_hfs_c, gf_hfs_b)
              write(40, FMT=Bab_data) '                                                  ', &
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
     write(40,'(a100)') '----------------------------------------------------------------------------------------------------'
  enddo 

contains
  !-----------------------------------------------------------------------------------------------------
  subroutine line_counting(fnum,filename, counter)
    implicit none
    
    integer, intent(in) :: fnum
    character(100), intent(in) :: filename
    integer, intent(out) :: counter
    character(180) :: string
    
    open(unit=fnum, file=filename, status='old', form='formatted', action='read')

    counter = 0    
    do
       read(fnum, '(a)', end=19) string
       counter = counter + 1
    end do
19  continue
    rewind(fnum)

    write(*,'(a42,i3)') ' The number of lines in the file is ', counter
    
  end subroutine line_counting
  !-----------------------------------------------------------------------------------------------------
  subroutine hfs_transitionE(DE_hfs)
    implicit none

    double precision, intent(in) ::
    double precision, intent(out) :: DE_hfs
    double precision :: C_lower, E_lower, Ehfs_lower, C_upper, E_upper, Ehfs_upper
    
    C_lower =
    Ehfs_lower = E_lower +
    C_upper =
    Ehfs_upper = E_upper
    DE_hfs = Ehfs_upper - Ehfs_lower
    
  end subroutine hfs_transitionE
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
    double precision, intent(in) :: gf_hfs_coul, gf_hfs_babus
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
end program hyperfine_transition_data
