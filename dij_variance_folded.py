# Plot the variance in Rij between all amino acids |i-j| apart as a function of |i-j|
# Set the directory
import os
os.chdir('/Users/jonathanhuihui/Documents/results_and_spreadsheets/resurrected_trx/distance_variance_folded/')

# Import the libraries
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# Import protein names 
protein_names = pd.read_table('protein_names.txt', delim_whitespace =True, header=None).values
line_length = len(protein_names)
lineage = 0

# Initialize sequence separation fluctuation array
#seq_sep_fluctarray = np.zeros([4,2])

for lineage in range(0, 2): #modify later for whole branch (2 instead of 1)
    for iter in range(0,line_length ): # line_length instead of 1
        # Get the filenames
        filenames = glob.glob('./' + protein_names[iter, lineage] + '/*.dat')
        num_files = len(filenames)*10
        del filenames
        
        # Initialize data array
        data_array = pd.read_table('./' + protein_names[iter, lineage] + '/' + str(protein_names[iter, lineage]) + '_distance_' 
        + str(0) + '.dat', delim_whitespace =True, header=None).values 
        seq_max = len(data_array)        
        
        # Fill data array        
        for frame in range(10, num_files, 10):
            current_frame = pd.read_table('./' + protein_names[iter, lineage] + '/' + str(protein_names[iter, lineage]) + '_distance_' 
            + str(frame) + '.dat', delim_whitespace =True, header=None).values
            data_array = np.dstack((data_array, current_frame))

        # Initialize fluctuations array
        fluct_array = np.zeros([seq_max, seq_max])
        
        # Fill fluctuations array
        for x in range(1, seq_max):
            for y in range(0, x):
                pair_variance = np.var(data_array[x,y,:])
                pair_mean = np.mean(data_array[x,y,:])
                pair_mean = pair_mean**2
                fluct_array[x,y] = pair_variance/pair_mean
        
        seq_sep_fluct_array = np.zeros([seq_max])        
        index = np.array([0])
        
      
        # Get sequence separation statistics
        fluct_holder = 0.0
        pairs = 0.0
        position = 0
        for spacing in range(1, seq_max):
            while (position < seq_max) and (position + spacing < seq_max):
                fluct_holder += fluct_array[int(position + spacing), int(position)]
                pairs += 1
                position += 1
            position = 0
            seq_sep_fluct_array[spacing] = fluct_holder/pairs
            index = np.append(index, spacing)
            fluct_holder = 0
            pairs = 0
        y = seq_sep_fluct_array[2:]
        X = index[2:]
        slope, intercept, r_value, p_value, std_err = stats.linregress(X,y)

        np.savetxt(protein_names[iter, lineage] + '_fluct_array.txt', y)
  

        # Visualize the distribution
        plt.figure()
        plt.scatter(X, y, color = 'blue')
        plt.title(protein_names[iter, lineage] + '_fluctuation vs. sequence sep')
        plt.ylabel('Fluctuation')
        plt.xlabel('Sequence separation |i-j|')
        plt.show()
        
        
        
        
        
        
        
