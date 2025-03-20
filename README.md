# Occupancy.py
This code will examine the OUTCAR, POSCAR, and POTCAR files to determine the occupancy of a system optimized in VASP with DFT.
It is designed for use with VASP.6.3.1, and uses numpy, re, and PyMatGen
Optionally, the CHGCAR file can be used to determine the shape of the total spin and spin difference arrays, as well as integrate for total magnetism
Outside of this, magnetism reported in the OUTCAR file when ISPIN = 2 is output by the code
