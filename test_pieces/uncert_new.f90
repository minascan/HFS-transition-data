program uncertainty
    implicit none
    double precision :: gf_hfs_coul, gf_hfs_babus, delta_gf, delta_log_gf, uncert

    print*, 'Give the gf value in the Coulomb gauge'
    read(*,'(f9.7)') gf_hfs_coul
    print*, 'Give the gf value in the Babushkin gauge'
    read(*,'(f9.7)') gf_hfs_babus

    delta_gf = gf_hfs_babus*abs(gf_hfs_coul - gf_hfs_babus)/max(gf_hfs_coul,gf_hfs_babus)

    delta_log_gf = ( abs( log10(gf_hfs_babus + delta_gf) - &
         log10(gf_hfs_babus - delta_gf) ) )/ 2

    uncert = delta_log_gf

    write(*,*) gf_hfs_coul, gf_hfs_babus
    write(*,*) max(gf_hfs_coul,gf_hfs_babus)
    write(*,*) gf_hfs_babus/max(gf_hfs_coul,gf_hfs_babus)
    write(*,*) delta_gf
    write(*,'(f9.7)') uncert
    
  end program uncertainty
