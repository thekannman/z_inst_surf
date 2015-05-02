# z_inst_surf
This program allows calculation of instantaneous liquid interfaces based on the work of Willard and Chandler (J. Phys. Chem. B, 2010, 114 (5), pp 1954–1958)

THIS PROGRAM HAS NOT BEEN THOROUGHLY TESTED!!!

Required libraries:
* GROMACS XTC library for reading position/velocity files. 
    http://www.gromacs.org/Developer_Zone/Programming_Guide/XTC_Library
* The Armadillo library for matrix calculations.
    http://arma.sourceforge.net/
* The Boost libraries for input option parsing, 4D tensor calculations, lexical casting, and output formatting.
    http://www.boost.org/


This program is still in the testing phase, and updates will be forthcoming. The current status is:
* The surface calculations seem to work.
* The average, deviation correlation, and density/temperature profiles are still being tested.
* Only GROMACS input is accepted, with positions in a .xtc file and velocities in a .trr file.
* Convolution parameters are hard-coded to match those in the Willard/Chandler paper.
