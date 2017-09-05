#!/usr/bin/env python
# Determine the fluctuations of energy 
# Trying to make this generalizable from the input of a "protein_info.txt" file
# that should contain the names of all the proteins. 
# Importing the libraries
import numpy as np
import pandas as pd
import sys 
import math

# Import the "protein_info.txt" file for name information 
protInfo_filename = sys.argv[1]
protInfo_file = pd.read_table(protInfo_filename, header = None , delim_whitespace=True)
protnames = protInfo_file.iloc[:,0]

# Define the subset for output file
subset = 'folded_state_fifty_clusters'

# Initialize number of clusters
# numclusters_all = 50

# Could also determine number of clusters from a column in "protein_info.txt"
numclusters_all = protInfo_file.iloc[:, 3]

# Determine number of proteins
numprots = len(protnames)

# Initialize an ouput table
export = np.empty(0)

# This code computes the variation of the energy values
for iter in range(0, numprots):
    stdDev = 0.0
# in case the number of clusters is not constant
    numclusters = numclusters_all[iter]   
# numclusters = 50
    cluster_weights_file_name = protnames[iter] + '_50_cluster_summary.dat'
    cluster_weights_file = pd.read_table(cluster_weights_file_name,
                                         delim_whitespace=True)
    weights = cluster_weights_file.iloc[0:numclusters, 2] 
    uint = np.empty(0)                     
    for z in range(0, numclusters):
        solvation_file_name = protnames[iter] + '_c' + str(z) + '.txt'
        solvation_file = pd.read_table(solvation_file_name,
                                         delim_whitespace=True,
                                         header = None)
        solv = np.sum(solvation_file.iloc[:,1])
        allOnTerms_file_name = 'allOnTerms_' + protnames[iter] +'_c' + str(z) + '.txt'
        allOnTerms_file = pd.read_table(allOnTerms_file_name,
                                         delim_whitespace=True,
                                         header = None)
        ucoul = allOnTerms_file.iloc[0,2]
        crfe = allOnTerms_file.iloc[0,0]
        upol = crfe - solv
        uint = np.append(uint, upol + ucoul)
    average = np.average(uint, weights=weights)
    variance = np.average((uint-average)**2, weights=weights) 
    stdDev = math.sqrt(variance)
#    outputFileName  = protnames[iter] + '_variation.txt'
    print uint 
#    f = open(outputFileName, 'w')
#    f.write(str(stdDev))
#    f.close
    





"""#This code reproduces the output values from the similar Mathematica code
for iter in range(0, numprots):
# in case the number of clusters is not constant
# numclusters = num_clusters_all[iter]   
    numclusters = 50
    cluster_weights_file_name = protnames[iter] + '_cluster_summary.dat'
    cluster_weights_file = pd.read_table(cluster_weights_file_name,
                                         delim_whitespace=True)
    cluster_weights = cluster_weights_file.iloc[:,2] 
    uint = np.empty(0)                     
    for z in range(0, numclusters):
        solvation_file_name = protnames[iter] + '_c' + str(z) + '.txt'
        solvation_file = pd.read_table(solvation_file_name,
                                         delim_whitespace=True,
                                         header = None)
        solv = np.sum(solvation_file.iloc[:,1])
        allOnTerms_file_name = 'allOnTerms_' + protnames[iter] +'_c' + str(z) + '.txt'
        allOnTerms_file = pd.read_table(allOnTerms_file_name,
                                         delim_whitespace=True,
                                         header = None)
        weight = cluster_weights[z]
        ucoul = allOnTerms_file.iloc[0,2]
        crfe = allOnTerms_file.iloc[0,0]
        upol = crfe - solv
        uint += (upol + ucoul)*weight
    export = np.append(export, uint)


outputfilename = subset + '_energy_results.txt'
np.savetxt(outputfilename, export, fmt="%.4f")
"""
