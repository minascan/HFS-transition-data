# hfs_trans
Program that computes parameters for transitions between hyperfine structure (HFS) levels 

The program assumes that a .ct.lsj file containing transition data for fine structure levels and .chlsj
files containing the HFS constants for each of these levels exist. The number of the targeted levels and
the nuclear spin of the system of interest are also taken as an input.

For each of the two levels involved in a fine structure transition, the possible F quantum numbers are
evaluated and the HFS transitions are defined. For each HFS transition, the energy, (vacuum) wavelength,
transition rate and oscillator strength gf (in both Babhuskin and Coulomb gauges) are computed.
The transition rates and oscillator strengths are evaluated based on 6j-symbols. The expressions can be
found in Sec. 17-9 of Cowan's book. The log(gf) values together with their uncertainties are also given.


--------------------------------------------------------------------------------------------------------
In the latest version of the code, contained in this repository, the output of the HFS transition
energies and wavelengths was added.

In short, the HFS transition energies are given by:
DE = E_upper + corr - E_lower + corr = E_up + (E_M1 + E_E2) - E_lo + (E_M1 + E_E2),
where
E_M1 = 0.5*A*C and E_E2 = B*...
where
C = F(F+1) - J(J+1) - I(I+1)

Steps:
1. Two extra input files "even....chlsj" and "odd....chlsj" need to be read
   so that we get the values of the A and B constants.
2. These values, together with the E_upper and E_lower, are used as input 
   to a subroutine that will compute the hfs transition energies. 
--------------------------------------------------------------------------------------------------------


