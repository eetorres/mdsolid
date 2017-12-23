#################################
#       FCC Al supercell        #
#################################
NAME	Al Fe bulk
ZNUM	13
# atomic mass
MNUM	26.982
# lattice constant
#LATT	2.698
LATT	4.05
# lattice vectors
LVEC
1.000000 0.000000 0.000000
0.000000 1.000000 0.000000
0.000000 0.000000 1.000000
# number of basis atoms
NATM	4
# basis atoms positon
APOS
0.000000 0.000000 0.000000
0.500000 0.500000 0.000000
0.000000 0.500000 0.500000
0.500000 0.000000 0.500000
# shift for the entire supercell
SHFT
0.25 0.25 0.25

#################################
#        MD parameters          #
#################################
# timestep
TIME	0.0001
# number of steps
NSTP	10000
# cut radius
CUTR	6.724
# shell thickness
SHLL	0.3
# temperature
TEMP	0.001
# unit cells in each direction
CLLS	6 6 6