# Count the occurrences of distances less than 8 Angstroms

# Import the libraries
import glob
import pandas as pd
import numpy as np
# Read dat files from directory
filenames = glob.glob('./*.dat') 



num_charges = 0
# Determine how many charged residues there are
interactions = len(filenames)
coeff = [1, -1,-2*interactions ]
num_charges = np.roots(coeff)
num_charges = num_charges[num_charges > 0]


# Count opposite charge interactions w/ distance less than 8A
edge = 0
for iter in range(0, interactions):
    file = filenames[iter].split('_')
    holder = file[1]
    primary_resi = holder[0:3]
    holder = file[2]
    secondary_resi = holder[0:3]
    if (primary_resi == 'arg') and (secondary_resi == 'glu' or secondary_resi == 'asp'):
        f = pd.read_table(filenames[iter], delim_whitespace=True)
        distance = f.iloc[0,1]
        if (float(distance) <= 8.0):
            edge+=1
    if (primary_resi == 'lys') and (secondary_resi == 'glu' or secondary_resi == 'asp'):
        f = pd.read_table(filenames[iter], delim_whitespace=True)
        f = pd.read_table(filenames[iter], delim_whitespace=True)
        distance = f.iloc[0,1]
        if (float(distance) <= 8.0):
            edge+=1
    if (primary_resi == 'hip') and (secondary_resi == 'glu' or secondary_resi == 'asp'):
        f = pd.read_table(filenames[iter], delim_whitespace=True)
        distance = f.iloc[0,1]
        if (float(distance) <= 8.0):
            edge+=1
    if (primary_resi == 'asp') and (secondary_resi == 'arg' or secondary_resi == 'lys'
                                                                or secondary_resi == 'hip'):
        f = pd.read_table(filenames[iter], delim_whitespace=True)
        distance = f.iloc[0,1]
        if (float(distance) <= 8.0):
            edge+=1
    if (primary_resi == 'glu') and (secondary_resi == 'arg' or secondary_resi == 'lys'
                                                                or secondary_resi == 'hip'):
        f = pd.read_table(filenames[iter], delim_whitespace=True)
        distance = f.iloc[0,1]
        if (float(distance) <= 8.0):
            edge+=1
print num_charges
print edge
print float(edge/num_charges)





