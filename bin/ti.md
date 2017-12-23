######## Ti bulk ########
NAME	Ti bulk (HCP)
ZNUM	22
# atomic mass
MNUM	22.000
# lattice constant 2.951
LATT	2.951

# lattice vectors
LVEC
 1.00000000 0.00000000 0.00000000
 0.00000000 1.73205015 0.00000000
 0.00000000 0.00000000 1.58787529
 
# number of atoms in the supercell
NATM	4

# atoms position in the supercell
APOS
 0.00000000 0.00000000 0.00000000
 0.50000000 0.50000000 0.00000000
 0.50000000 0.16666667 0.50000000
 0.00000000 0.66666667 0.50000000
# shift for the entire supercell
SHFT
0.0 0.0 0.0
# offset for the entire supercell
OSET
0.0 0.0 0.0
########  MD parameters ##########
TIME	0.00001
# number of steps
NSTP	10000
# cut radius 6.72
CUTR	5.1
# shell thickness
SHLL	0.3
# temperature
TEMP	1.0
# unit cells in each direction
CLLS	4 4 4