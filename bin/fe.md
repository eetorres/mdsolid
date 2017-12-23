#################################
#       BCC Fe supercell        #
#################################
NAME	BCC Fe bulk
ZNUM	26
# atomic mass
MNUM	55.845
# lattice constant
LATT	2.8664
# lattice vectors
LVEC
1.000000 0.000000 0.000000
0.000000 1.000000 0.000000
0.000000 0.000000 1.000000
# number of basis atoms
NATM	2
# basis atoms positon
APOS
0.000000 0.000000 0.000000
0.500000 0.500000 0.500000
# shift for the entire supercell
SHFT
0.25 0.25 0.25
#################################
#        MD parameters          #
#################################
# time step
TIME	0.0001
# number of steps
NSTP	10000
# cut radius
CUTR	4.095
# shell thickness
SHLL	0.3
# temperature
TEMP	1.0
# unit cells in each direction
CLLS	10 10 10