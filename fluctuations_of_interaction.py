# Determine the fluctuations of energy 
# Trying to make this generalizable from the input of a "protein_info.txt" file
# that should contain the names of all the proteins. 
# Importing the libraries
import numpy as np
import pandas as pd
import sys 
import math
from scipy import stats
import matplotlib.pyplot as plt

# Import the "protein_info.txt" file for name information 
# protInfo_filename = sys.argv[1]
protInfo_filename = 'protein_info.txt'
protInfo_file = pd.read_table(protInfo_filename, header = None , delim_whitespace=True)
protnames = protInfo_file.iloc[:,0]

# Define the subset for output file
subset = 'folded_state_fifty_clusters'

# Initialize number of clusters
numclusters_all = 50

# Could also determine number of clusters from a column in "protein_info.txt"
# numclusters_all = protInfo_file.iloc[:,INSERT COLUMN NUMBER HERE DUMMY]

# Determine number of proteins
numprots = len(protnames)
# Initialize an ouput table
export = np.empty(0)

# This code computes the variation of the energy values
for iter in xrange(0, numprots, 2):
# in case the number of clusters is not constant
# numclusters = num_clusters_all[iter]
    prot_family = protnames[iter].split('_')[0]
    uint_array_1 = np.empty(0)
    uint_array_2 = np.empty(0)
    weights_array_1 = np.empty(0)
    weights_array_2 = np.empty(0)
    for pair_num in range(0, 2):
        numclusters = 50
        cluster_weights_file_name = protnames[iter + pair_num] + '_cluster_summary.dat'
        cluster_weights_file = pd.read_table(cluster_weights_file_name,
                                         delim_whitespace=True)
        weights = cluster_weights_file.iloc[:,2] 
        uint = np.empty(0)
        cut_weights = np.empty(0)
        stdDev = 0.0                    
        for z in range(0, numclusters):
            solvation_file_name = protnames[iter + pair_num] + '_c' + str(z) + '.txt'
            solvation_file = pd.read_table(solvation_file_name,
                                         delim_whitespace=True,
                                         header = None)
            solv = np.sum(solvation_file.iloc[:,1])
            allOnTerms_file_name = 'allOnTerms_' + protnames[iter + pair_num] +'_c' + str(z) + '.txt'
            allOnTerms_file = pd.read_table(allOnTerms_file_name,
                                         delim_whitespace=True,
                                         header = None)
            ucoul = allOnTerms_file.iloc[0,2]
            crfe = allOnTerms_file.iloc[0,0]
            upol = crfe - solv
            uint = np.append(uint, upol + ucoul)
            #if weights[z] > 0:
            #    uint = np.append(uint, upol + ucoul)
            #    cut_weights = np.append(cut_weights, weights[z])
        if pair_num == 0:
            uint_array_1 = uint
            weights_array_1 = weights
        elif pair_num == 1:
            uint_array_2 = uint
            weights_array_2 = weights
#        average = np.average(uint, weights=weights)
#        variance = np.average((uint-average)**2, weights=weights) 
#        stdDev = math.sqrt(variance)
    uint_array_1 = np.multiply(uint_array_1, weights_array_1)
    uint_array_2 = np.multiply(uint_array_2, weights_array_2)
    pval = stats.wilcoxon(uint_array_1, uint_array_2)
#    outputFileName  = prot_family + '_pval.txt'
#    f = open(outputFileName, 'w')
#    f.write(str(pval))
#    f.close       
    plt.scatter(uint_array_1, weights_array_1, color='blue')
    plt.scatter(uint_array_2, weights_array_2, color='red')
    plt.xlabel('Uint')
    plt.ylabel('Weights')
    plt.title(protnames[iter] + ' ' + protnames[iter + pair_num] + ' folded energy distributions')
    plt.show
    plt.savefig(protnames[iter] + '_' + protnames[iter + pair_num] + '_energy_distributions.pdf' )
    print protnames[iter]
    print protnames[iter + pair_num]
    print pval
#    print norm_test_1
#    print norm_test_2







