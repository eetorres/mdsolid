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
 0.00000000 1.73205143 0.00000000
 0.00000000 0.00000000 1.58500000
 
# number of atoms in the supercell
NATM	4

# atoms position in the supercell
APOS
 0.00000000 0.00000000 0.00000000
 0.50000000 0.50000000 0.00000000
 0.00000000 0.33333333 0.50000000
 0.50000000 0.83333333 0.50000000
  
# offset for the entire supercell
OSET
0.25 0.25 0.25
########  MD parameters ##########
TIME	0.0005
# number of steps
NSTP	10000
# cut radius 6.72
CUTR    5.1
# shell thickness
SHLL	0.0
# temperature
TEMP	1.0
# unit cells in each direction
CLLS	6 10 6