# hfs_trans
Program that computes parameters for transitions between hyperfine structure levels 

In this repository we monitor the progress of the source code that is written to
provide additional output, such as the transition energies and wavelengths that were
missing from the initial version.

The hfs transition energies will be given by:
DE = E_upper + corr - E_lower + corr = E_up + (E_M1 + E_E2) - E_lo + (E_M1 + E_E2),
where
E_M1 = 0.5*A*C and E_E2 = B*...
where
C = F(F+1) - J(J+1) - I(I+1)

Steps:
1. Two extra input files "even....chlsj" and "odd....chlsj" need to be read
   so that we get the values of the A and B constants.
2. These values, together with the E_upper and E_lower, can be used as input 
   to a subroutine that will compute the hfs transition energies. 
   
